import { useEffect, useMemo, useRef, useState } from 'react'
import type { BinaryInfo, Bridge, Done } from './bridge'
import {
  BY_NAME, INITIAL, PARAMS, SECTIONS, changed, commandLine, duration, invalid, problems, progressOf, setValue,
  suggestOutput, writtenOf, type Param, type Progress, type Values, type Written
} from './model'
import logo from './logo.svg'

const STORE = 'dx.values.v1'
const MAX_LOG = 5000

// Saved values from an older CLI may name options that no longer exist, or hold ones that are now the default.
function restore(): Values {
  try {
    const saved = JSON.parse(localStorage.getItem(STORE) ?? '{}') as Values
    let v: Values = { ...INITIAL }
    for (const [k, x] of Object.entries(saved)) {
      if (k in INITIAL) v[k] = String(x)
      else if (BY_NAME.has(k)) v = setValue(v, k, String(x))
    }
    return v
  } catch {
    return { ...INITIAL }
  }
}

type Phase = { kind: 'idle' } | { kind: 'running'; since: number } | { kind: 'done'; done: Done; written: Written | null }

export default function App({ bridge }: { bridge: Bridge }) {
  const [info, setInfo] = useState<BinaryInfo | null>(null)
  const [values, setValues] = useState<Values>(restore)
  const [autoOut, setAutoOut] = useState(false) // the output is our suggestion, so a new input may replace it
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [query, setQuery] = useState('')
  const [phase, setPhase] = useState<Phase>({ kind: 'idle' })
  const [log, setLog] = useState<string[]>([])
  const [progress, setProgress] = useState<Progress | null>(null)
  const [now, setNow] = useState(Date.now())
  const [copied, setCopied] = useState(false)
  const written = useRef<Written | null>(null)
  const logBox = useRef<HTMLPreElement>(null)

  useEffect(() => {
    bridge.probe().then(setInfo)
  }, [bridge])
  useEffect(() => {
    localStorage.setItem(STORE, JSON.stringify(values))
  }, [values])
  useEffect(() => {
    const offLog = bridge.onLog((line) => {
      setLog((l) => (l.length >= MAX_LOG ? [...l.slice(1 - MAX_LOG), line] : [...l, line]))
      const p = progressOf(line)
      if (p) setProgress(p)
      written.current = writtenOf(line) ?? written.current
    })
    const offDone = bridge.onDone((done) => setPhase({ kind: 'done', done, written: written.current }))
    return () => {
      offLog()
      offDone()
    }
  }, [bridge])
  useEffect(() => {
    if (phase.kind !== 'running') return
    const t = setInterval(() => setNow(Date.now()), 500)
    return () => clearInterval(t)
  }, [phase.kind])
  useEffect(() => {
    // follow the log inside its box, not the page
    if (logBox.current) logBox.current.scrollTop = logBox.current.scrollHeight
  }, [log])

  const running = phase.kind === 'running'
  const issues = problems(values)
  const nChanged = changed(values).length

  const setInput = (path: string) => {
    setValues((v) => ({ ...v, in: path, ...(autoOut || !v.out ? { out: suggestOutput(path) } : {}) }))
    if (autoOut || !values.out) setAutoOut(true)
  }

  const start = async () => {
    setLog([])
    setProgress(null)
    written.current = null
    setPhase({ kind: 'running', since: Date.now() })
    setNow(Date.now())
    try {
      await bridge.run(values)
    } catch (e) {
      setLog([String(e)])
      setPhase({ kind: 'done', done: { code: -1, cancelled: false, elapsedMs: 0 }, written: null })
    }
  }

  const q = query.trim().toLowerCase()
  const visible = useMemo(
    () =>
      PARAMS.filter(
        (p) =>
          (showAdvanced || !p.advanced || p.name in values) &&
          (!q || p.name.toLowerCase().includes(q) || p.description.toLowerCase().includes(q))
      ),
    [showAdvanced, q, values]
  )

  return (
    <div className="app">
      <header>
        <img src={logo} alt="" width={40} height={40} />
        <div>
          <h1>DIAspeXtract</h1>
          <p className="sub">Pseudo-MS/MS spectra from DIA runs, for DDA search engines</p>
        </div>
        <div className="cli" title={info?.bin}>
          {info === null ? (
            <span className="badge">checking the CLI…</span>
          ) : info.version ? (
            <span className="badge ok" data-testid="cli-version">CLI {info.version}</span>
          ) : (
            <span className="badge bad">CLI not working</span>
          )}
        </div>
      </header>

      {info && !info.version && (
        <div className="banner" role="alert">
          The command-line tool at <code>{info.bin}</code> ({info.source}) did not start.
          {info.detail && <pre>{info.detail}</pre>}
        </div>
      )}

      <section className="card">
        <h2>Run</h2>
        <div className="row">
          <label htmlFor="in">Input</label>
          <input id="in" value={values.in} spellCheck={false} disabled={running}
            placeholder="a Bruker .d folder, an mzML, mzPeak or Thermo .raw file, or s3://bucket/run.mzpeak"
            onChange={(e) => setInput(e.target.value)} />
          <button disabled={running} onClick={async () => { const p = await bridge.pickRunDir(); if (p) setInput(p) }}>.d folder…</button>
          <button disabled={running} onClick={async () => { const p = await bridge.pickRunFile(); if (p) setInput(p) }}>Run file…</button>
        </div>
        <div className="row">
          <label htmlFor="out">Output</label>
          <input id="out" value={values.out} spellCheck={false} disabled={running} placeholder="spectra.mzML"
            onChange={(e) => { setAutoOut(false); setValues((v) => ({ ...v, out: e.target.value })) }} />
          <button disabled={running} onClick={async () => {
            const p = await bridge.pickOutput(values.out || undefined)
            if (p) { setAutoOut(false); setValues((v) => ({ ...v, out: p })) }
          }}>Save as…</button>
        </div>
        <div className="row">
          <label htmlFor="threads">Threads</label>
          <input id="threads" className="short" value={values.threads} disabled={running} inputMode="numeric"
            aria-invalid={!/^\d+$/.test(values.threads)}
            onChange={(e) => setValues((v) => ({ ...v, threads: e.target.value }))} />
          <span className="hint">0 = all cores</span>
        </div>
      </section>

      <section className="card">
        <div className="paramsbar">
          <h2>Parameters</h2>
          <span className="hint" data-testid="changed-count">
            {nChanged === 0 ? 'all defaults' : `${nChanged} changed`}
          </span>
          <input type="search" placeholder="Filter parameters" aria-label="Filter parameters" value={query}
            onChange={(e) => setQuery(e.target.value)} />
          <label className="check">
            <input type="checkbox" checked={showAdvanced} onChange={(e) => setShowAdvanced(e.target.checked)} />
            Advanced
          </label>
          <button disabled={running || nChanged === 0}
            onClick={() => setValues((v) => Object.fromEntries(Object.entries(v).filter(([k]) => k in INITIAL)))}>
            Reset all
          </button>
        </div>
        {SECTIONS.map((s) => {
          const ps = visible.filter((p) => p.section === s.name)
          if (ps.length === 0) return null
          const n = ps.filter((p) => p.name in values).length
          return (
            <details key={s.name} open={!!q || s.name === '' || n > 0 || undefined} data-section={s.name || 'general'}>
              <summary>
                {s.description} {s.name && <code>{s.name}</code>} {n > 0 && <span className="pill">{n} changed</span>}
              </summary>
              {ps.map((p) => (
                <Field key={p.name} p={p} value={values[p.name] ?? p.default} isSet={p.name in values}
                  disabled={running} onChange={(v) => setValues((x) => setValue(x, p.name, v))} />
              ))}
            </details>
          )
        })}
        {visible.length === 0 && <p className="hint">No parameter matches “{query}”.</p>}
      </section>

      <section className="card">
        <h2>Command line</h2>
        <div className="cmd">
          <code data-testid="command-line">{commandLine(values)}</code>
          <button onClick={() => navigator.clipboard?.writeText(commandLine(values)).then(() => {
            setCopied(true); setTimeout(() => setCopied(false), 1500)
          })}>{copied ? 'Copied' : 'Copy'}</button>
        </div>
      </section>

      <section className="card runbar">
        {running ? (
          <button className="danger" onClick={() => bridge.cancel()}>Cancel</button>
        ) : (
          <button className="primary" disabled={issues.length > 0 || !info?.version} onClick={start}>Run</button>
        )}
        <div className="status">
          {running && (
            <>
              <progress aria-label="Progress" max={progress?.total ?? undefined} value={progress?.done ?? undefined} />
              <span data-testid="status">
                {progress ? `window ${progress.done} of ${progress.total}` : 'starting…'} · {duration(now - phase.since)}
              </span>
            </>
          )}
          {phase.kind === 'idle' && issues.length > 0 && <span className="hint" data-testid="status">{issues[0]}</span>}
          {phase.kind === 'done' && <Outcome phase={phase} reveal={(p) => bridge.reveal(p)} />}
        </div>
      </section>

      {log.length > 0 && (
        <section className="card">
          <h2>Log</h2>
          <pre className="log" aria-label="Log" ref={logBox}>
            {log.join('\n')}
          </pre>
        </section>
      )}
    </div>
  )
}

function Outcome({ phase, reveal }: { phase: Extract<Phase, { kind: 'done' }>; reveal: (p: string) => void }) {
  const { done, written } = phase
  if (done.cancelled) return <span className="warn" data-testid="status">Cancelled.</span>
  if (done.code !== 0)
    return <span className="bad" data-testid="status">Failed: diaspextract exited with code {done.code ?? 'none'}. The log says why.</span>
  return (
    <span className="good" data-testid="status">
      Done in {duration(done.elapsedMs)}
      {written && <>: {written.spectra.toLocaleString('en')} {written.spectra === 1 ? 'spectrum' : 'spectra'} written to <code>{written.path}</code></>}
      {written && <button onClick={() => reveal(written.path)}>Show in folder</button>}
    </span>
  )
}

function Field({ p, value, isSet, disabled, onChange }: {
  p: Param; value: string; isSet: boolean; disabled: boolean; onChange: (v: string) => void
}) {
  const bad = invalid(p, value)
  const bool = p.type === 'bool' || (p.choices?.length === 2 && p.choices.includes('true') && p.choices.includes('false'))
  let control
  if (bool)
    control = <input type="checkbox" aria-label={p.name} checked={value === 'true'} disabled={disabled}
      onChange={(e) => onChange(String(e.target.checked))} />
  else if (p.choices)
    control = (
      <select aria-label={p.name} value={value} disabled={disabled} onChange={(e) => onChange(e.target.value)}>
        {p.choices.map((c) => <option key={c} value={c}>{c || '(none)'}</option>)}
      </select>
    )
  else
    control = <input aria-label={p.name} value={value} disabled={disabled} spellCheck={false} aria-invalid={!!bad}
      inputMode={p.type === 'int' || p.type === 'double' ? 'decimal' : undefined}
      onChange={(e) => onChange(e.target.value)} />
  return (
    <div className={`field${isSet ? ' set' : ''}${bad ? ' invalid' : ''}`} data-param={p.name}>
      <div className="name">
        <code>{p.name}</code>
        {p.advanced && <span className="pill dim">advanced</span>}
      </div>
      <div className="control">
        {control}
        {isSet && (
          <button className="link" disabled={disabled} title={`back to the default, ${p.default || '(none)'}`}
            onClick={() => onChange(p.default)}>reset</button>
        )}
        {!isSet && <span className="hint">default</span>}
      </div>
      <p className="desc">
        {p.description}
        {bad && <span className="bad"> Must be {bad}.</span>}
      </p>
    </div>
  )
}
