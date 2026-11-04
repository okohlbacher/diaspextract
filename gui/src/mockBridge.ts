// A stand-in backend for a plain browser: `npm run dev` and the Playwright UI tests. It replays the CLI's log format
// (the lines model.ts reads) on a timer. Tests steer it through window.__dxMock before or during a run.

import type { Bridge, Done } from './bridge'
import type { Values } from './model'

export interface MockControl {
  version: string | null // null: the CLI "does not start"
  windows: number
  spectra: number
  stepMs: number // delay between log lines
  exitCode: number
  runs: Values[] // every run's values, for assertions
  picks: { dir: string | null; file: string | null; out: string | null }
  revealed: string[]
}

declare global {
  interface Window {
    __dxMock: MockControl
  }
}

export function mockBridge(): Bridge {
  const ctl: MockControl = (window.__dxMock ??= {
    version: '1.9.0',
    windows: 4,
    spectra: 1234,
    stepMs: 60,
    exitCode: 0,
    runs: [],
    picks: { dir: '/data/HeLa_5min.d', file: '/data/orbitrap.raw', out: null },
    revealed: []
  })
  const logs = new Set<(l: string) => void>()
  const dones = new Set<(d: Done) => void>()
  let timer: ReturnType<typeof setTimeout> | null = null
  let finish: ((d: Done) => void) | null = null

  const on = <T>(set: Set<(v: T) => void>) => (cb: (v: T) => void) => {
    set.add(cb)
    return () => void set.delete(cb)
  }

  return {
    probe: async () => ({
      bin: '/mock/bin/diaspextract',
      source: 'browser preview',
      version: ctl.version,
      detail: ctl.version ? '' : 'dyld: Library not loaded: @rpath/libOpenMS.dylib'
    }),
    run: async (values) => {
      if (finish) throw new Error('a run is already in progress')
      ctl.runs.push(values)
      const t0 = Date.now()
      const lines = [
        `Loaded 28 frames. [t=0s]`,
        ...Array.from({ length: ctl.windows }, (_, i) => `  window ${i + 1}/${ctl.windows} 590.0-620.0: spectra ${i + 1} [t=0s]`),
        ...(ctl.exitCode === 0
          ? [`Wrote ${ctl.spectra} pseudo-MS2 spectra to ${values.out} (mzML)`]
          : ['Error: Unexpected internal error (mock failure)'])
      ]
      const emit = (d: Done) => {
        if (timer) clearTimeout(timer)
        timer = null
        finish = null
        dones.forEach((cb) => cb(d))
      }
      finish = emit
      const step = (i: number) => {
        if (i < lines.length) {
          logs.forEach((cb) => cb(lines[i]))
          timer = setTimeout(() => step(i + 1), ctl.stepMs)
        } else emit({ code: ctl.exitCode, cancelled: false, elapsedMs: Date.now() - t0 })
      }
      timer = setTimeout(() => step(0), ctl.stepMs)
    },
    cancel: async () => {
      if (!finish) return false
      finish({ code: null, cancelled: true, elapsedMs: 0 })
      return true
    },
    pickRunDir: async () => ctl.picks.dir,
    pickRunFile: async () => ctl.picks.file,
    pickOutput: async () => ctl.picks.out,
    reveal: async (p) => void ctl.revealed.push(p),
    onLog: on(logs),
    onDone: on(dones)
  }
}
