#!/usr/bin/env python3
"""The mechanism, counted at the source: how many INFERRED PRECURSORS share a monoisotope, and what the two seed
traces behind such a pair look like. Uses the tool's own funnel dump (-diag:dump_ms1_tsv), which carries the
seed trace_idx for every precursor, so nothing has to be inferred from the emitted mzML.

Each duplicate pair is classified by the m/z RELATION of its two seed traces (traces.tsv):
  SAME-TRACE-MZ   the two seeds sit at the same m/z within --ppm: the MS1 detector produced two traces for one
                  feature, and the claim loop has no offset-0 case, so each one seeds a precursor
  ISO-STEP k      the seeds are k isotope steps apart: two envelope hypotheses that walked down to one mono
  OTHER           neither

Usage: funnel_dup.py <funnel.precursors.tsv> <funnel.traces.tsv> [--ppm=20] [--drt=1.5] [--im=0.01]
"""
import sys, bisect, collections
ISO = 1.0033548
prec_p, tr_p = sys.argv[1], sys.argv[2]
kw = dict(a.lstrip("-").split("=") for a in sys.argv[3:])
PPM, RT_TOL, IM_TOL = float(kw.get("ppm", 20)) * 1e-6, float(kw.get("drt", 1.5)), float(kw.get("im", 0.01))

tr = []
with open(tr_p) as f:
    next(f)
    for line in f:
        p = line.split("\t"); tr.append((float(p[0]), float(p[1]), float(p[2])))
P = []
with open(prec_p) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        P.append((float(p[0]), int(p[1]), float(p[2]), float(p[3]), int(p[4]), int(p[5])))
print(f"{len(tr):,} MS1 traces -> {len(P):,} precursors ({len(P)/max(1,len(tr)):.2f} per trace)", flush=True)
print(f"distinct seed traces used: {len({x[5] for x in P}):,}", flush=True)

byz = collections.defaultdict(list)
for i, p in enumerate(P): byz[p[1]].append(p)
for z in byz: byz[z].sort(key=lambda p: p[0])
keys = {z: [p[0] for p in v] for z, v in byz.items()}

cls = collections.Counter(); ndup = collections.Counter(); pairs = 0
for p in P:
    mono, z, rt, im, _, ti = p
    tol = mono * PPM
    arr, k = byz[z], keys[z]
    hits = [o for o in (arr[j] for j in range(bisect.bisect_left(k, mono - tol), bisect.bisect_right(k, mono + tol)))
            if o[5] != ti and abs(o[2] - rt) <= RT_TOL and abs(o[3] - im) <= IM_TOL]
    ndup[min(len(hits), 5)] += 1
    for o in hits:
        if o[5] <= ti: continue
        pairs += 1
        a, b = tr[ti][0], tr[o[5]][0]
        d = abs(a - b)
        if d <= max(a, b) * PPM: cls["SAME-TRACE-MZ"] += 1
        else:
            for kk in range(1, 6):
                if abs(d - kk * ISO / max(z, 1)) <= max(a, b) * PPM: cls["ISO-STEP %d" % kk] += 1; break
            else: cls["OTHER (d=%.3f)" % round(d, 1) if d < 5 else "OTHER (far)"] += 1

print(f"\n{pairs:,} same-monoisotope precursor PAIRS (same charge, {PPM*1e6:.0f} ppm, dRT<={RT_TOL}, dIM<={IM_TOL})")
for c, n in cls.most_common(10): print(f"  {c:<18}{n:>9,}  {100*n/max(1,pairs):>5.1f}%")
print("\nprecursors by number of same-mono competitors:")
for n, c in sorted(ndup.items()): print(f"  {n}{'+' if n == 5 else ''} competitors: {c:>9,}  {100*c/len(P):>5.1f}%")
