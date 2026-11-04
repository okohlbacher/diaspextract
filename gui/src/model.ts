// Everything the UI decides without a backend: the parameter list, validation, the command line it shows, the output
// name it suggests and what it reads from the CLI's log. Pure functions, unit-tested in tests/unit.spec.ts.

import manifest from './params.generated.json' with { type: 'json' }

export interface Param {
  name: string
  section: string
  type: string
  default: string
  description: string
  advanced: boolean
  min?: number
  max?: number
  choices?: string[]
}
export interface Section {
  name: string
  description: string
}

export const PARAMS: Param[] = manifest.params as Param[]
export const SECTIONS: Section[] = [{ name: '', description: 'General' }, ...manifest.sections]
export const BY_NAME = new Map(PARAMS.map((p) => [p.name, p]))

// What the user set. in/out/threads are always present; a parameter is present only while it differs from the
// CLI's default, so the command line carries only the changes and the CLI's defaults stay the defaults.
export type Values = Record<string, string>

export const INITIAL: Values = { in: '', out: '', threads: '0' }

/** Why a value is not acceptable, or null. Mirrors the CLI's own checks, so a bad run is not even started. */
export function invalid(p: Param, v: string): string | null {
  if (p.choices) return p.choices.includes(v) ? null : `one of ${p.choices.map((c) => c || '(empty)').join(', ')}`
  if (p.type === 'int' || p.type === 'double') {
    const t = v.trim()
    if (t === '' || !Number.isFinite(Number(t))) return 'a number'
    const n = Number(t)
    if (p.type === 'int' && !Number.isInteger(n)) return 'a whole number'
    if (p.min !== undefined && n < p.min) return `at least ${p.min}`
    if (p.max !== undefined && n > p.max) return `at most ${p.max}`
  }
  return null
}

/** Set one value; going back to the default removes it. */
export function setValue(values: Values, name: string, v: string): Values {
  const next = { ...values }
  const p = BY_NAME.get(name)
  if (p && v === p.default) delete next[name]
  else next[name] = v
  return next
}

export function changed(values: Values): string[] {
  return Object.keys(values).filter((k) => BY_NAME.has(k))
}

export function problems(values: Values): string[] {
  const out: string[] = []
  if (!values.in?.trim()) out.push('Choose an input run.')
  if (!values.out?.trim()) out.push('Choose an output file.')
  if (!/^\d+$/.test(values.threads ?? '')) out.push('Threads must be a whole number (0 = all cores).')
  for (const k of changed(values)) {
    const why = invalid(BY_NAME.get(k)!, values[k])
    if (why) out.push(`-${k} must be ${why}.`)
  }
  return out
}

/** The argv the backend builds (src-tauri/src/lib.rs build_args), for display. Flags (bool) are bare. */
export function args(values: Values): string[] {
  const a = ['-in', values.in, '-out', values.out, '-threads', values.threads]
  for (const k of changed(values).sort()) {
    if (BY_NAME.get(k)!.type === 'bool') {
      if (values[k] === 'true') a.push(`-${k}`)
    } else a.push(`-${k}`, values[k])
  }
  return a
}

const quote = (s: string) => (s === '' ? "''" : /[^\w@%+=:,./\\-]/.test(s) ? `'${s.replace(/'/g, `'\\''`)}'` : s)
export const commandLine = (values: Values) => ['diaspextract', ...args(values)].map(quote).join(' ')

/** <input stem>.pseudo.mzML beside the input; a .d directory loses its .d. */
export function suggestOutput(input: string): string {
  const t = input.replace(/[\\/]+$/, '')
  if (!t) return ''
  const cut = Math.max(t.lastIndexOf('/'), t.lastIndexOf('\\'))
  return `${t.slice(0, cut + 1)}${t.slice(cut + 1).replace(/\.(d|mzml|mzpeak)$/i, '')}.pseudo.mzML`
}

export interface Progress {
  done: number
  total: number
}
/** "  window 3/32 590.0-620.0: spectra 41 [t=12s]" -> 3 of 32. */
export function progressOf(line: string): Progress | null {
  const m = /^\s*window (\d+)\/(\d+)\b/.exec(line)
  return m ? { done: Number(m[1]), total: Number(m[2]) } : null
}

export interface Written {
  spectra: number
  path: string
}
/** "Wrote 1234 pseudo-MS2 spectra to out.mzML (mzML)" */
export function writtenOf(line: string): Written | null {
  const m = /^Wrote (\d+) pseudo-MS2 spectra to (.+) \(\w+\)\s*$/.exec(line)
  return m ? { spectra: Number(m[1]), path: m[2] } : null
}

export function duration(ms: number): string {
  const s = Math.floor(ms / 1000)
  const hh = Math.floor(s / 3600)
  const mm = String(Math.floor((s % 3600) / 60)).padStart(2, '0')
  const ss = String(s % 60).padStart(2, '0')
  return hh ? `${hh}:${mm}:${ss}` : `${mm}:${ss}`
}
