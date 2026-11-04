// End to end: the built desktop app, its Rust backend and the bundled CLI on a real (synthetic) mzML.
//
// Playwright reaches the app's web view only where it speaks the Chrome DevTools protocol: WebView2 on Windows, given
// a debugging port through WEBVIEW2_ADDITIONAL_BROWSER_ARGUMENTS. macOS (WKWebView) and Linux (WebKitGTK) have no such
// port; there CI runs the bundled CLI from inside the built app instead (build.yml, "desktop" job).
//
//   DX_E2E_APP=path/to/diaspextract-gui.exe DX_E2E_INPUT=synth.mzML DX_E2E_VERSION=1.9.0 npx playwright test --project=e2e

import { chromium, expect, test, type Browser, type Page } from '@playwright/test'
import { spawn, type ChildProcess } from 'node:child_process'
import { existsSync, mkdtempSync, readFileSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'

const APP = process.env.DX_E2E_APP
const PORT = 9333
test.skip(!APP, 'set DX_E2E_APP to the built app (Windows)')

let app: ChildProcess
let browser: Browser
let page: Page

test.beforeAll(async () => {
  app = spawn(APP!, [], {
    env: { ...process.env, WEBVIEW2_ADDITIONAL_BROWSER_ARGUMENTS: `--remote-debugging-port=${PORT}` },
    stdio: 'inherit'
  })
  for (let i = 0; ; i++) {
    try {
      browser = await chromium.connectOverCDP(`http://127.0.0.1:${PORT}`)
      break
    } catch (e) {
      if (i > 120) throw e
      await new Promise((r) => setTimeout(r, 500))
    }
  }
  const ctx = browser.contexts()[0]
  page = ctx.pages()[0] ?? (await ctx.waitForEvent('page'))
  await page.waitForLoadState()
})

test.afterAll(async () => {
  await browser?.close().catch(() => {})
  app?.kill()
})

test('the app finds and starts its bundled CLI', async () => {
  await expect(page.getByTestId('cli-version')).toHaveText(`CLI ${process.env.DX_E2E_VERSION}`, { timeout: 60_000 })
})

test('a run writes spectra through the real CLI', async () => {
  const out = join(mkdtempSync(join(tmpdir(), 'dx-e2e-')), 'synth.pseudo.mzML')
  await page.getByLabel('Input').fill(process.env.DX_E2E_INPUT!)
  await page.getByLabel('Output').fill(out)
  await page.getByLabel('Threads').fill('2')
  await page.getByRole('button', { name: 'Run', exact: true }).click()
  await expect(page.getByTestId('status')).toContainText(' written to', { timeout: 120_000 })
  await expect(page.getByLabel('Log')).toContainText('window 1/1')
  expect(existsSync(out)).toBe(true)
  expect(readFileSync(out, 'utf8')).toMatch(/<spectrum /)
})

test('a run that fails is reported as failed', async () => {
  await page.getByLabel('Input').fill(join(tmpdir(), 'does-not-exist.mzML'))
  await page.getByRole('button', { name: 'Run', exact: true }).click()
  await expect(page.getByTestId('status')).toHaveText(/Failed: diaspextract exited with code/, { timeout: 60_000 })
})
