// The one seam between the UI and the backend. main.tsx installs the Tauri implementation inside the app and the
// mock (mockBridge.ts) in a plain browser, which is how the UI tests and `npm run dev` run without the CLI.

import type { Values } from './model'

export interface BinaryInfo {
  bin: string
  source: string
  version: string | null
  detail: string
}
export interface Done {
  code: number | null
  cancelled: boolean
  elapsedMs: number
}

export interface Bridge {
  probe(): Promise<BinaryInfo>
  run(values: Values): Promise<void>
  cancel(): Promise<boolean>
  pickRunDir(): Promise<string | null>
  pickRunFile(): Promise<string | null>
  pickOutput(defaultPath?: string): Promise<string | null>
  reveal(path: string): Promise<void>
  onLog(cb: (line: string) => void): () => void
  onDone(cb: (d: Done) => void): () => void
}

export async function tauriBridge(): Promise<Bridge> {
  const { invoke } = await import('@tauri-apps/api/core')
  const { listen } = await import('@tauri-apps/api/event')
  const { open, save } = await import('@tauri-apps/plugin-dialog')
  const { revealItemInDir } = await import('@tauri-apps/plugin-opener')

  // listen() resolves later; hand back a synchronous unsubscribe that also covers "unsubscribed before it resolved"
  const sub = <T>(event: string, cb: (p: T) => void) => {
    let un: (() => void) | null = null
    let gone = false
    listen<T>(event, (e) => cb(e.payload)).then((u) => (gone ? u() : (un = u)))
    return () => {
      gone = true
      un?.()
    }
  }
  const one = (r: string | string[] | null) => (typeof r === 'string' ? r : null)

  return {
    probe: () => invoke('probe'),
    run: (values) => invoke('run', { values }),
    cancel: () => invoke('cancel'),
    pickRunDir: async () => one(await open({ directory: true, title: 'Bruker .d run' })),
    pickRunFile: async () =>
      one(await open({ filters: [{ name: 'Runs', extensions: ['mzML', 'mzml', 'mzpeak', 'raw', 'RAW'] }], title: 'mzML, mzPeak or Thermo .raw run' })),
    pickOutput: async (defaultPath) =>
      one(await save({ defaultPath, filters: [{ name: 'mzML', extensions: ['mzML'] }, { name: 'mzPeak', extensions: ['mzpeak'] }] })),
    reveal: (path) => revealItemInDir(path),
    onLog: (cb) => sub<string>('dx:log', cb),
    onDone: (cb) => sub<Done>('dx:done', cb)
  }
}
