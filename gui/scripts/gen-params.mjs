// The GUI's parameter list, read from the tool's own -write_ini: the CLI is the one source of names, types, defaults,
// limits and help texts, so the form cannot drift from it. CI regenerates this from the bundled binary and fails on
// a difference.
//
//   node scripts/gen-params.mjs path/to/diaspextract > src/params.generated.json

import { execFileSync } from 'node:child_process'
import { mkdtempSync, readFileSync, rmSync } from 'node:fs'
import { tmpdir } from 'node:os'
import { join } from 'node:path'

const bin = process.argv[2] || process.env.DIASPEXTRACT_BIN || 'diaspextract'
const dir = mkdtempSync(join(tmpdir(), 'dx-ini-'))
const ini = join(dir, 'dx.ini')
try {
  execFileSync(bin, ['-write_ini', ini], { stdio: ['ignore', 'ignore', 'pipe'] })
} catch (e) {
  console.error(`failed to run ${bin} -write_ini: ${e.message}`)
  process.exit(1)
}
const xml = readFileSync(ini, 'utf8')
rmSync(dir, { recursive: true, force: true })

const unescape = (s) =>
  s.replace(/&quot;/g, '"').replace(/&apos;/g, "'").replace(/&lt;/g, '<').replace(/&gt;/g, '>').replace(/&amp;/g, '&')
const attrs = (tag) => Object.fromEntries([...tag.matchAll(/(\w+)="([^"]*)"/g)].map((m) => [m[1], unescape(m[2])]))

// TOPP's own plumbing, not DIAspeXtract options; the GUI sets -in, -out and -threads itself.
const SKIP = new Set(['in', 'out', 'threads', 'version', 'log', 'debug', 'no_progress', 'force', 'test', 'ini', 'write_ini'])

const params = []
const sections = []
const stack = [] // open NODE names; [0] is the tool, [1] the instance "1"
for (const m of xml.matchAll(/<NODE\s([^>]*?)>|<\/NODE>|<ITEM\s([^>]*?)\/>/g)) {
  if (m[0] === '</NODE>') { stack.pop(); continue }
  if (m[1] !== undefined) {
    const a = attrs(m[1])
    stack.push(a.name)
    if (stack.length === 3) sections.push({ name: a.name, description: a.description ?? '' })
    continue
  }
  const a = attrs(m[2])
  const name = [...stack.slice(2), a.name].join(':')
  if (stack.length < 2 || SKIP.has(name) || name.startsWith('diag:selftest')) continue
  const p = {
    name,
    section: stack[2] ?? '',
    type: a.type,
    // plain decimals: the INI writes 0.005 as 5.0e-03
    default: /^(int|double)$/.test(a.type) && /e/i.test(a.value) ? String(Number(a.value)) : a.value,
    description: a.description ?? '',
    advanced: a.advanced === 'true'
  }
  if (a.restrictions) {
    if (a.type === 'int' || a.type === 'double') {
      const [lo, hi] = a.restrictions.split(':')
      if (lo) p.min = Number(lo)
      if (hi) p.max = Number(hi)
    } else p.choices = a.restrictions.split(',')
  }
  params.push(p)
}
process.stdout.write(JSON.stringify({ sections, params }, null, 2) + '\n')
