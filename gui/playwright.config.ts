import { defineConfig } from '@playwright/test'

// unit: model.ts, no browser. ui: the React app in Chromium against the mock backend (vite dev server).
// e2e: the built desktop app with the real CLI, through WebView2's DevTools port (Windows; see tests/e2e.spec.ts).
export default defineConfig({
  reporter: process.env.CI ? [['list'], ['html', { open: 'never' }]] : 'list',
  retries: process.env.CI ? 1 : 0,
  projects: [
    { name: 'unit', testMatch: 'unit.spec.ts' },
    { name: 'ui', testMatch: 'ui.spec.ts', use: { baseURL: 'http://localhost:1420', browserName: 'chromium' } },
    { name: 'e2e', testMatch: 'e2e.spec.ts', timeout: 180_000 }
  ],
  webServer: process.env.DX_E2E_APP
    ? undefined
    : { command: 'npm run dev', url: 'http://localhost:1420', reuseExistingServer: !process.env.CI }
})
