import { expect, test } from '@playwright/test'
import {
  BY_NAME, INITIAL, PARAMS, args, changed, commandLine, invalid, problems, progressOf, setValue, suggestOutput, writtenOf
} from '../src/model.ts'

const base = { ...INITIAL, in: '/d/run.d', out: '/d/run.pseudo.mzML' }

test('the manifest is the CLI contract', () => {
  expect(PARAMS.length).toBeGreaterThan(40)
  for (const n of ['in', 'out', 'threads', 'ini', 'write_ini', 'diag:selftest_twins']) expect(BY_NAME.has(n)).toBe(false)
  expect(BY_NAME.get('gate:delta_rt')).toMatchObject({ type: 'double', min: 0 })
  expect(BY_NAME.get('trace:detector')?.choices).toEqual(['openms', 'integer'])
})

test('a value equal to the default is dropped, so only changes reach the command line', () => {
  let v = setValue(base, 'gate:delta_rt', '5')
  expect(changed(v)).toEqual(['gate:delta_rt'])
  v = setValue(v, 'gate:delta_rt', BY_NAME.get('gate:delta_rt')!.default)
  expect(changed(v)).toEqual([])
  expect(args(v)).toEqual(['-in', '/d/run.d', '-out', '/d/run.pseudo.mzML', '-threads', '0'])
})

test('args sorts options like the backend and passes flags bare', () => {
  let v = setValue(base, 'trace:detector', 'openms')
  v = setValue(v, 'assembly:min_fragments', '4')
  v = setValue(v, 'diag:spectrum_meta', 'true')
  expect(args(v).slice(6)).toEqual(['-assembly:min_fragments', '4', '-diag:spectrum_meta', '-trace:detector', 'openms'])
  expect(args(setValue(v, 'diag:spectrum_meta', 'false')).includes('-diag:spectrum_meta')).toBe(false)
})

test('validation follows the declared type and limits', () => {
  const p = (n: string) => BY_NAME.get(n)!
  expect(invalid(p('gate:min_correlation'), '0.5')).toBeNull()
  expect(invalid(p('gate:min_correlation'), '1.5')).toBe('at most 1')
  expect(invalid(p('assembly:min_fragments'), '2.5')).toBe('a whole number')
  expect(invalid(p('assembly:min_fragments'), 'x')).toBe('a number')
  expect(invalid(p('trace:detector'), 'magic')).toMatch(/^one of/)
  expect(problems({ ...INITIAL })).toEqual(['Choose an input run.', 'Choose an output file.'])
  expect(problems({ ...base, threads: '-1' })).toEqual(['Threads must be a whole number (0 = all cores).'])
  expect(problems(setValue(base, 'gate:delta_rt', '-3'))).toEqual(['-gate:delta_rt must be at least 0.'])
})

test('the shown command line quotes what a shell would split', () => {
  expect(commandLine({ ...base, in: "/My Data/it's.d" })).toContain(`-in '/My Data/it'\\''s.d'`)
  expect(commandLine({ ...base, in: 'C:\\runs\\a.d' })).toContain('-in C:\\runs\\a.d')
})

test('the suggested output sits beside the input', () => {
  expect(suggestOutput('/x/HeLa_5min.d')).toBe('/x/HeLa_5min.pseudo.mzML')
  expect(suggestOutput('/x/HeLa_5min.d/')).toBe('/x/HeLa_5min.pseudo.mzML')
  expect(suggestOutput('C:\\runs\\orbi.mzML')).toBe('C:\\runs\\orbi.pseudo.mzML')
  expect(suggestOutput('run.mzML')).toBe('run.pseudo.mzML')
  expect(suggestOutput('')).toBe('')
})

test('the log lines the UI reads', () => {
  expect(progressOf('  window 3/32 590.0-620.0: spectra 41 [t=12s]')).toEqual({ done: 3, total: 32 })
  expect(progressOf('Processing 32 isolation windows')).toBeNull()
  expect(writtenOf('Wrote 1234 pseudo-MS2 spectra to /o/a b.mzML (mzML)')).toEqual({ spectra: 1234, path: '/o/a b.mzML' })
})
