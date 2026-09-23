#!/usr/bin/env python3
"""The control the matched-random arm cannot be: remove spectra that have a co-eluting, co-mobile, same-charge
NEIGHBOUR which is NOT on the isotope lattice. Such a spectrum is just as "well covered" as a twin -- it has a
similar friend at the same time and mobility -- but it is not a second claim on the same monoisotope. If the twin
rule only wins because twinned peptides have many spectra, this control wins just as much; if the same-monoisotope
relation is what matters, it does not.

Size-matched to k0_mid (314,857) by taking the candidates closest in m/z, and charge-matched by filling each charge
to that charge's k0_mid count where the candidates allow.
Usage: neigh_rule.py <iso.tsv> <drop_k0_mid.txt> <out.txt>
"""
import sys, bisect, collections
ISO, PPM = 1.0033548, 20.0
iso_p, ref_p, out_p = sys.argv[1:4]

S = []
with open(iso_p) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        S.append((int(p[0]), float(p[1]), float(p[2]), int(p[3] or 0), float(p[4] or 0), int(p[5] or 0), int(p[7] or 0)))
info = {s[0]: s for s in S}
ref = {int(x) for x in open(ref_p).read().split()}
target = collections.Counter(info[i][3] for i in ref)
print("k0_mid removed by charge:", dict(sorted(target.items())), flush=True)

byz = collections.defaultdict(list)
for s in S: byz[s[3]].append(s)
for z in byz: byz[z].sort(key=lambda s: s[2])
mzs = {z: [s[2] for s in v] for z, v in byz.items()}
weaker = lambda a, b: (a[6], a[5], -a[0]) < (b[6], b[5], -b[0])

cand = collections.defaultdict(list)   # charge -> (|dmz| to the neighbour, index)
for s in S:
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: continue
    arr, keys = byz[z], mzs[z]
    best = None
    for j in range(bisect.bisect_left(keys, mz - 3.0), bisect.bisect_right(keys, mz + 3.0)):
        o = arr[j]
        if o[0] == idx or abs(o[1] - rt) > 1.5 or abs(o[4] - im) > 0.01: continue
        d = o[2] - mz
        if any(abs(d - k * ISO / z) <= mz * PPM * 1e-6 for k in range(-3, 4)): continue   # on the lattice: not a control
        if not weaker(s, o): continue                                                     # remove the weaker one, as the rule does
        if best is None or abs(d) < best: best = abs(d)
    if best is not None: cand[z].append((best, idx))

R = set()
for z, c in cand.items():
    c.sort()
    take = min(target.get(z, 0), len(c))
    R.update(i for _, i in c[:take])
    if take < target.get(z, 0): print(f"  z={z}: only {len(c):,} candidates for {target[z]:,} wanted", flush=True)
short = len(ref) - len(R)
if short > 0:                                       # top up from the deepest pool so the SIZE matches exactly
    rest = sorted((d, i) for z in cand for d, i in cand[z] if i not in R)
    R.update(i for _, i in rest[:short])
with open(out_p, "w") as g: g.write("".join(f"{i}\n" for i in sorted(R)))
print(f"wrote {len(R):,} (k0_mid removed {len(ref):,}); overlap with k0_mid = {len(R & ref):,}", flush=True)
