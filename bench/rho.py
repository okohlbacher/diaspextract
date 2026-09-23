#!/usr/bin/env python3
"""The one number that separates redundancy from coverage, and the safety property that was asserted but never
measured.

rho = of the IDENTIFIED spectra a rule removes, the fraction whose RETAINED twin carries the SAME peptide.
      Redundancy predicts rho near the both-identified chance level for this data, which the project's own
      isotope study puts at 96% (docs/ISOTOPE-DUPLICATION-2026-09-03.md); coverage predicts far lower.
orphan = removed spectra with NO retained twin at all -- the "removes only the weaker member" property, which is
      FALSE for the greedy pairwise rule on chains and TRUE by construction for the union-find clusters.

Usage: rho.py <iso.tsv> <results.sage.tsv> <drop.txt> [--im=0.01] [--drt=1.5]
"""
import sys, csv, bisect, collections
ISO, PPM = 1.0033548, 20.0
iso_p, sage_p, drop_p = sys.argv[1:4]
kw = dict((a.lstrip("-").split("=", 1) + ["1"])[:2] for a in sys.argv[4:])
IM_TOL, RT_TOL = float(kw.get("im", 0.01)), float(kw.get("drt", 1.5))

S = []
with open(iso_p) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        S.append((int(p[0]), float(p[1]), float(p[2]), int(p[3] or 0), float(p[4] or 0), int(p[5] or 0), int(p[7] or 0)))
drop = {int(x) for x in open(drop_p).read().split()}
pep = {}
for r in csv.DictReader(open(sage_p, newline=""), delimiter="\t"):
    if r["label"] == "1" and float(r["spectrum_q"]) <= 0.01: pep[int(r["scannr"].split("=")[-1])] = r["peptide"]

byz = collections.defaultdict(list)
for s in S: byz[s[3]].append(s)
for z in byz: byz[z].sort(key=lambda s: s[2])
mzs = {z: [s[2] for s in v] for z, v in byz.items()}

def twins(s):
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: return []
    tol = mz * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    return [o for o in (arr[j] for j in range(bisect.bisect_left(keys, mz - tol), bisect.bisect_right(keys, mz + tol)))
            if o[0] != idx and abs(o[1] - rt) <= RT_TOL and abs(o[4] - im) <= IM_TOL]

same = diff = none_id = orphan = n_id = 0
orphan_all = 0
peps_of = collections.defaultdict(set)
for i, p in pep.items(): peps_of[p].add(i)
for s in S:
    if s[0] not in drop: continue
    kept = [o for o in twins(s) if o[0] not in drop]
    if not kept: orphan_all += 1
    if s[0] not in pep: continue
    n_id += 1
    if not kept: orphan += 1; continue
    kp = {pep.get(o[0]) for o in kept} - {None}
    if pep[s[0]] in kp: same += 1
    elif kp: diff += 1
    else: none_id += 1

thin = sum(1 for p, sp in peps_of.items() if len(sp) == 1 and next(iter(sp)) in drop)
thin_tot = sum(1 for p, sp in peps_of.items() if len(sp) == 1)
print(f"removed {len(drop):,} spectra, {n_id:,} of them identified ({100*n_id/len(drop):.1f}%)")
print(f"  rho (retained twin carries the SAME peptide)  {same:>8,}  {100*same/max(1,n_id):>5.1f}%")
print(f"  retained twin identified, DIFFERENT peptide   {diff:>8,}  {100*diff/max(1,n_id):>5.1f}%")
print(f"  retained twin exists but is unidentified      {none_id:>8,}  {100*none_id/max(1,n_id):>5.1f}%")
print(f"  ORPHAN: no retained twin at all               {orphan:>8,}  {100*orphan/max(1,n_id):>5.1f}%   (all removed: {orphan_all:,} = {100*orphan_all/len(drop):.1f}%)")
print(f"singly-identified peptides whose one spectrum is removed: {thin:,} of {thin_tot:,} ({100*thin/max(1,thin_tot):.1f}%)"
      f"  vs identified-spectrum removal rate {100*n_id/max(1,len(pep)):.1f}%")
