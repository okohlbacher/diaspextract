# DIAspeXtract desktop app

A Tauri 2 window around the command-line tool. The app carries the packaged CLI (`bin/`, `lib/`, `share/OpenMS`) as a
resource and runs it as a child process; it never links OpenMS itself.

```
src/model.ts               parameters, validation, command line, output name, log parsing (pure, unit-tested)
src/params.generated.json  the option list, generated from the CLI's own -write_ini (npm run params)
src/App.tsx                the window
src/bridge.ts              the seam to the backend; src/mockBridge.ts stands in for it in a plain browser
src-tauri/src/lib.rs       finds the CLI, builds its argv from an allowlist, streams its output, cancels it
tests/                     Playwright: unit.spec.ts, ui.spec.ts (mock backend), e2e.spec.ts (the built app)
```

The form is generated, not written: names, types, defaults, limits and help texts come from the CLI, and CI fails when
`params.generated.json` differs from what the packaged CLI reports. The command line carries only the options that
differ from their defaults, so the CLI's defaults stay authoritative. The backend accepts only the options the
manifest declares and passes them as an argv array, never through a shell.

## Develop

```bash
npm ci
npm run dev                        # the UI in a browser, against the mock backend
cp -RL <unpacked DIAspeXtract-<platform>> src-tauri/resources/diaspextract
npm run params                     # after a CLI option changed
npx tauri dev                      # the real app with that CLI
npx tauri build --bundles app      # (dmg, nsis, deb on the other platforms)
```

`DIASPEXTRACT_BIN` points the app at another CLI binary.

## Test

```bash
npx playwright install chromium
npm test                           # unit and UI suites
(cd src-tauri && cargo test --lib) # argv building and the allowlist
```

The end-to-end suite drives the built app over the Chrome DevTools protocol, which only WebView2 offers, so it runs on
Windows (build.yml, desktop job, against the installed app):

```bash
DX_E2E_APP=...\diaspextract-gui.exe DX_E2E_INPUT=synth.mzML DX_E2E_VERSION=1.9.1 npm run test:e2e
```

On macOS and Linux CI starts the CLI from inside the finished app instead.
