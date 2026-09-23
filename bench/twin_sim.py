#!/usr/bin/env python3
"""Fragment-content similarity of k=0 twins -- the decisive merge-vs-drop question.

Two precursors at the SAME inferred monoisotope (src/diaspextract.cpp:3038 writes pc.mono_mz) are the same
precursor claimed twice. Their SPECTRA need not be: each is built from its own seed trace, so the fragment
sets are correlated against different XICs.
  HIGH cosine -> the two spectra say the same thing; dropping the weaker one costs nothing.
  LOW cosine  -> they carry complementary fragments; dropping loses content and MERGING is the fix.
Controls: an off-lattice co-eluting neighbour of the same charge, and an unrelated random pair.
Usage: twin_sim.py <iso.tsv> <pseudo.mzML> [--n=1500] [--im=0.01] [--drt=1.5] [--sage=results.sage.tsv]
"""
import sys, random, collections, math, base64, struct, re, bisect
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from iso_sim import cosine, fetch                      # same binning and decoder as the isotope study

ISO, PPM = 1.0033548, 20.0
iso_p, mz_p = sys.argv[1], sys.argv[2]
kw = dict(a.lstrip("-").split("=") for a in sys.argv[3:])
N, IM_TOL, RT_TOL = int(kw.get("n", 1500)), float(kw.get("im", 0.01)), float(kw.get("drt", 1.5))

S = []
with open(iso_p) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        S.append((int(p[0]), float(p[1]), float(p[2]), int(p[3] or 0), float(p[4] or 0), int(p[5] or 0), int(p[7] or 0)))
byz = collections.defaultdict(list)
for s in S: byz[s[3]].append(s)
for z in byz: byz[z].sort(key=lambda s: s[2])
mzs = {z: [s[2] for s in v] for z, v in byz.items()}

def near(s, k, im_tol, rt_tol):
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: return []
    t = mz - k * ISO / z; tol = t * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    return [o for o in (arr[j] for j in range(bisect.bisect_left(keys, t - tol), bisect.bisect_right(keys, t + tol)))
            if o[0] != idx and abs(o[1] - rt) <= rt_tol and abs(o[4] - im) <= im_tol]

rng = random.Random(7)
twins, neigh = [], []
pool = [s for s in S if s[3] > 0]
rng.shuffle(pool)
for s in pool:
    if len(twins) >= N: break
    t = near(s, 0, IM_TOL, RT_TOL)
    if t: twins.append((s, min(t, key=lambda o: abs(o[4] - s[4]))))
for s in pool:                                        # off-lattice control: co-eluting, same z, NOT on any lattice step
    if len(neigh) >= N: break
    if near(s, 0, IM_TOL, RT_TOL) or any(near(s, k, 0.05, 3.0) for k in (1, 2, 3)): continue
    z, mz, rt, im = s[3], s[2], s[1], s[4]
    arr, keys = byz[z], mzs[z]
    lo = bisect.bisect_left(keys, mz - 3.0)
    cand = [o for o in arr[lo:bisect.bisect_right(keys, mz + 3.0)]
            if o[0] != s[0] and abs(o[1] - rt) <= RT_TOL and abs(o[4] - im) <= IM_TOL
            and all(abs((o[2] - mz) - k * ISO / z) > mz * PPM * 1e-6 for k in (-3, -2, -1, 0, 1, 2, 3))]
    if cand: neigh.append((s, min(cand, key=lambda o: abs(o[2] - mz))))
rnd = [(pool[i], pool[-1 - i]) for i in range(min(N, len(pool) // 2))]

want = {x[0] for pr in twins + neigh + rnd for x in pr}
print(f"fetching {len(want):,} spectra", flush=True)
sp = fetch(mz_p, set(want))
pep = {}
if "sage" in kw:
    import csv
    for r in csv.DictReader(open(kw["sage"], newline=""), delimiter="\t"):
        if r["label"] == "1" and float(r["spectrum_q"]) <= 0.01: pep[int(r["scannr"].split("=")[-1])] = r["peptide"]

def report(name, pairs, extra=False):
    cs = [cosine(sp[a[0]], sp[b[0]]) for a, b in pairs if a[0] in sp and b[0] in sp]
    if not cs: print(f"{name}: none"); return
    cs.sort(); n = len(cs)
    print(f"{name:<12} n={n:<6} median {cs[n//2]:.3f}  mean {sum(cs)/n:.3f}  "
          f"q1 {cs[n//4]:.3f}  q3 {cs[3*n//4]:.3f}  >0.9 {100*sum(c>0.9 for c in cs)/n:.1f}%  >0.75 {100*sum(c>0.75 for c in cs)/n:.1f}%")
    if extra:
        both = [(a, b) for a, b in pairs if pep.get(a[0]) and pep.get(a[0]) == pep.get(b[0])]
        one = [(a, b) for a, b in pairs if bool(pep.get(a[0])) != bool(pep.get(b[0]))]
        for lbl, grp in (("both-same-peptide", both), ("only-one-identified", one)):
            c2 = [cosine(sp[a[0]], sp[b[0]]) for a, b in grp if a[0] in sp and b[0] in sp]
            if c2: c2.sort(); print(f"   {lbl:<22} n={len(c2):<5} median {c2[len(c2)//2]:.3f}  >0.9 {100*sum(c>0.9 for c in c2)/len(c2):.1f}%")
        dn = [abs(a[6] - b[6]) for a, b in pairs]; dp = [abs(a[5] - b[5]) for a, b in pairs]
        print(f"   envelope size differs by 0 in {100*sum(d==0 for d in dn)/len(dn):.0f}% of pairs; "
              f"peak count differs by 0 in {100*sum(d==0 for d in dp)/len(dp):.0f}%; median |dnpeaks| {sorted(dp)[len(dp)//2]}")

report("k0 TWIN", twins, extra=True)
report("neighbour", neigh)
report("random", rnd)
