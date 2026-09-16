#!/usr/bin/env python3
"""Summarise a DIASPEXTRACT_LEDGER trace: which structure holds the memory, when, and for how long.

    ledger_report.py <ledger.tsv> [more.tsv ...]

The tool samples every 250 ms: RSS, the glibc arena/free figures, live bytes per data-structure
category, the unattributed remainder (RSS minus the sum: allocator retention plus what is not
charged) and the phase marker. This prints, per file:

  * the peak of RSS and of every category, with the second it happened and the phase it was in,
  * the phase timeline: duration, RSS at entry/exit/peak, and the categories that dominate it,
  * the residency of each category: the fraction of the run it holds more than 5% of the peak RSS
    (a structure that is large AND long-lived is what a streaming design has to attack), and
  * the per-tile picture when the phases are tile markers.

Pure stdlib; the TSV is small (a 4-minute run is ~1,000 rows).
"""
import sys, collections

def load(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) != len(hdr): continue
            r = {}
            for k, v in zip(hdr, p):
                r[k] = v if k == "phase" else float(v)
            rows.append(r)
    return hdr, rows

def main(paths):
    for path in paths:
        hdr, rows = load(path)
        if not rows: print(f"{path}: empty"); continue
        # the tool wrote "_mb" until 2026-09-10 and "_mib" after it; the values were always MiB
        suf = "_mib" if any(h.endswith("_mib") for h in hdr) else "_mb"
        fixed = {"rss", "arena", "free", "sum", "unattributed"}
        cats = [h for h in hdr if h.endswith(suf) and h[: -len(suf)] not in fixed]
        R, S, U, F = "rss" + suf, "sum" + suf, "unattributed" + suf, "free" + suf
        t0, t1 = rows[0]["t_s"], rows[-1]["t_s"]
        print(f"\n=== {path}  ({t1 - t0:.0f} s sampled, {len(rows)} rows)")
        peak_rss = max(rows, key=lambda r: r[R])
        peak_sum = max(rows, key=lambda r: r[S])
        print(f"peak RSS      {peak_rss[R]:9.0f} MiB at t={peak_rss['t_s']:7.1f} s  phase={peak_rss['phase']}"
              f"   (charged {peak_rss[S]:.0f} MiB, unattributed {peak_rss[U]:.0f} MiB,"
              f" glibc free {peak_rss[F]:.0f} MiB)")
        print(f"peak charged  {peak_sum[S]:9.0f} MiB at t={peak_sum['t_s']:7.1f} s  phase={peak_sum['phase']}"
              f"   (RSS {peak_sum[R]:.0f} MiB)")
        print("\nper structure: peak, when, and how long it is above 5% of the peak RSS")
        thr = 0.05 * peak_rss[R]
        span = max(t1 - t0, 1e-9)
        for c in sorted(cats, key=lambda c: -max(r[c] for r in rows)):
            pk = max(rows, key=lambda r: r[c])
            if pk[c] <= 0: continue
            live = span * sum(1 for r in rows if r[c] > thr) / len(rows)
            at_peak = peak_rss[c]
            print(f"  {c[: -len(suf)]:<12} peak {pk[c]:8.0f} MiB at t={pk['t_s']:7.1f} ({pk['phase']:<10})"
                  f"  {at_peak:8.0f} MiB at the RSS peak   above 5% of peak RSS for {live:6.1f} s ({100.0 * live / span:4.1f}%)")
        print("\nphases")
        seq, cur = [], None
        for r in rows:
            if cur is None or r["phase"] != cur["phase"]:
                cur = {"phase": r["phase"], "t0": r["t_s"], "t1": r["t_s"], "rows": [r]}
                seq.append(cur)
            else:
                cur["t1"] = r["t_s"]; cur["rows"].append(r)
        for p in seq:
            rs = p["rows"]
            pk = max(rs, key=lambda r: r[R])
            dom = sorted(((max(r[c] for r in rs), c) for c in cats), reverse=True)[:3]
            dom_s = ", ".join(f"{c[: -len(suf)]} {v:.0f}" for v, c in dom if v > 0)
            print(f"  {p['phase']:<12} {p['t1'] - p['t0']:7.1f} s   RSS {rs[0][R]:7.0f} -> {rs[-1][R]:7.0f} MiB"
                  f"  (peak {pk[R]:7.0f})   top: {dom_s}")
        tiles = [p for p in seq if p["phase"].startswith(("tile", "loop", "write"))]
        if tiles:
            print("\nper tile (read | loop | write), RSS peak and the charged peak")
            byk = collections.defaultdict(dict)
            for p in tiles:
                kind = "read" if p["phase"].startswith("tile") else p["phase"][:4]
                k = "".join(ch for ch in p["phase"] if ch.isdigit())
                byk[k][kind] = p
            for k in sorted(byk, key=lambda x: int(x or 0)):
                parts = []
                for kind in ("read", "loop", "writ"):
                    p = byk[k].get(kind)
                    if not p: continue
                    pk = max(p["rows"], key=lambda r: r[R])
                    ps = max(p["rows"], key=lambda r: r[S])
                    parts.append(f"{kind} {p['t1'] - p['t0']:6.1f}s rss<={pk[R]:7.0f} charged<={ps[S]:7.0f}")
                print(f"  tile {k:>3}: " + " | ".join(parts))

if __name__ == "__main__":
    if len(sys.argv) < 2: raise SystemExit(__doc__)
    main(sys.argv[1:])
