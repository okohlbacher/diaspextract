// The UI in Chromium against the mock backend (src/mockBridge.ts), steered through window.__dxMock.
import { expect, test, type Page } from '@playwright/test'

const run = (page: Page) => page.getByRole('button', { name: 'Run', exact: true })
const status = (page: Page) => page.getByTestId('status')
// Fields sit in collapsed sections; a user finds one with the filter, which opens every matching section.
const find = async (page: Page, name: string) => {
  await page.getByLabel('Filter parameters').fill(name)
  return page.getByLabel(name, { exact: true })
}
const mock = (page: Page, set: Record<string, unknown>) =>
  page.evaluate((s) => Object.assign(window.__dxMock, s), set)

test.beforeEach(async ({ page }) => {
  await page.goto('/')
  await page.evaluate(() => localStorage.clear())
  await page.reload()
  await expect(page.getByTestId('cli-version')).toHaveText('CLI 1.9.0')
})

test('Run waits for an input and an output', async ({ page }) => {
  await expect(run(page)).toBeDisabled()
  await expect(status(page)).toHaveText('Choose an input run.')
  await page.getByRole('button', { name: '.d folder…' }).click()
  await expect(page.getByLabel('Input')).toHaveValue('/data/HeLa_5min.d')
  await expect(page.getByLabel('Output')).toHaveValue('/data/HeLa_5min.pseudo.mzML')
  await expect(run(page)).toBeEnabled()
})

test('a new input replaces a suggested output, never a chosen one', async ({ page }) => {
  await page.getByLabel('Input').fill('/a/one.d')
  await expect(page.getByLabel('Output')).toHaveValue('/a/one.pseudo.mzML')
  await page.getByRole('button', { name: 'Run file…' }).click()
  await expect(page.getByLabel('Output')).toHaveValue('/data/orbitrap.pseudo.mzML')
  await page.getByLabel('Output').fill('/mine/keep.mzML')
  await page.getByLabel('Input').fill('/a/two.d')
  await expect(page.getByLabel('Output')).toHaveValue('/mine/keep.mzML')
})

test('a changed parameter is marked, counted and on the command line; reset undoes it', async ({ page }) => {
  await page.getByLabel('Input').fill('/a/r.d')
  const f = await find(page, 'gate:delta_rt')
  await expect(f).toHaveValue('4.0')
  await f.fill('6')
  await expect(page.getByTestId('changed-count')).toHaveText('1 changed')
  await expect(page.getByTestId('command-line')).toHaveText(
    'diaspextract -in /a/r.d -out /a/r.pseudo.mzML -threads 0 -gate:delta_rt 6')
  await page.locator('[data-param="gate:delta_rt"]').getByRole('button', { name: 'reset' }).click()
  await expect(page.getByTestId('changed-count')).toHaveText('all defaults')
  await expect(f).toHaveValue('4.0')
})

test('an out-of-range value blocks Run and says why', async ({ page }) => {
  await page.getByLabel('Input').fill('/a/r.d')
  await (await find(page, 'gate:min_correlation')).fill('2')
  await expect(page.locator('[data-param="gate:min_correlation"]')).toContainText('Must be at most 1.')
  await expect(run(page)).toBeDisabled()
  await expect(status(page)).toHaveText('-gate:min_correlation must be at most 1.')
})

test('choices and switches use their own controls', async ({ page }) => {
  await page.getByLabel('Input').fill('/a/r.d')
  await (await find(page, 'trace:detector')).selectOption('openms')
  await (await find(page, 'dnoise:ms1')).uncheck()
  await expect(page.getByTestId('command-line')).toHaveText(/-dnoise:ms1 false -trace:detector openms$/)
})

test('advanced options are hidden until asked for, and the filter searches names and help', async ({ page }) => {
  await expect(await find(page, 'gate:shrinkage')).toHaveCount(0)
  await page.getByLabel('Advanced').check()
  await expect(page.getByLabel('gate:shrinkage', { exact: true })).toBeVisible()
  await page.getByLabel('Advanced').uncheck()
  await page.getByLabel('Filter parameters').fill('valley')
  await expect(page.locator('[data-param]')).toHaveCount(2)
  await page.getByLabel('Filter parameters').fill('no such thing')
  await expect(page.getByText('No parameter matches')).toBeVisible()
})

test('a run shows progress, the log and the result, and can reveal the output', async ({ page }) => {
  await mock(page, { windows: 5, spectra: 4321, stepMs: 150 })
  await page.getByLabel('Input').fill('/a/r.d')
  await page.getByLabel('Threads').fill('8')
  await run(page).click()
  await expect(page.getByRole('button', { name: 'Cancel' })).toBeVisible()
  await expect(page.getByLabel('Input')).toBeDisabled()
  await expect(status(page)).toContainText(/window \d of 5/)
  await expect(status(page)).toContainText('4,321 spectra written to /a/r.pseudo.mzML', { timeout: 10_000 })
  await expect(page.getByLabel('Log')).toContainText('window 5/5')
  expect(await page.evaluate(() => window.__dxMock.runs)).toEqual([{ in: '/a/r.d', out: '/a/r.pseudo.mzML', threads: '8' }])
  await page.getByRole('button', { name: 'Show in folder' }).click()
  expect(await page.evaluate(() => window.__dxMock.revealed)).toEqual(['/a/r.pseudo.mzML'])
})

test('Cancel stops the run', async ({ page }) => {
  await mock(page, { windows: 50, stepMs: 200 })
  await page.getByLabel('Input').fill('/a/r.d')
  await run(page).click()
  await page.getByRole('button', { name: 'Cancel' }).click()
  await expect(status(page)).toHaveText('Cancelled.')
  await expect(run(page)).toBeEnabled()
})

test('a failing run says so', async ({ page }) => {
  await mock(page, { exitCode: 3, windows: 1 })
  await page.getByLabel('Input').fill('/a/r.d')
  await run(page).click()
  await expect(status(page)).toHaveText(/Failed: diaspextract exited with code 3/)
  await expect(page.getByLabel('Log')).toContainText('mock failure')
})

test('settings survive a restart', async ({ page }) => {
  await page.getByLabel('Input').fill('/a/r.d')
  await page.getByLabel('max_charge').fill('4')
  await page.reload()
  await expect(page.getByLabel('Input')).toHaveValue('/a/r.d')
  await expect(page.getByLabel('max_charge')).toHaveValue('4')
  await expect(page.getByTestId('changed-count')).toHaveText('1 changed')
})

test('a CLI that does not start is reported and blocks Run', async ({ page }) => {
  await page.addInitScript(() => {
    window.__dxMock = { version: null, windows: 1, spectra: 1, stepMs: 10, exitCode: 0, runs: [], revealed: [],
      picks: { dir: null, file: null, out: null } }
  })
  await page.reload()
  await expect(page.getByRole('alert')).toContainText('did not start')
  await page.getByLabel('Input').fill('/a/r.d')
  await expect(run(page)).toBeDisabled()
})
