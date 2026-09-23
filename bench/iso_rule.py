#!/usr/bin/env python3
"""Dump the exact spectrum indices a removal rule deletes, one per line, for bench/mzml_drop.py.
Same partner logic as iso_frontier.py (same lattice, same charge, co-eluting, co-mobile).
Usage: iso_rule.py <iso.tsv> <outdir>
"""
import sys, bisect, collections, random
ISO, PPM = 1.0033548, 20.0
iso_p, outdir = sys.argv[1], sys.argv[2].rstrip("/")

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

def partners(s, k, rt_tol, im_tol):
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: return []
    t = mz - k * ISO / z; tol = t * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    return [o for o in (arr[j] for j in range(bisect.bisect_left(keys, t - tol), bisect.bisect_right(keys, t + tol)))
            if o[0] != idx and abs(o[1] - rt) <= rt_tol and abs(o[4] - im) <= im_tol]

weaker = lambda a, b: (a[6], a[5], -a[0]) < (b[6], b[5], -b[0])

def k0(rt_tol, im_tol):
    R = set()
    for s in S:
        for o in partners(s, 0, rt_tol, im_tol):
            if weaker(s, o) and o[0] not in R: R.add(s[0]); break
    return R

def k1_heavy(im_tol):
    return {s[0] for s in S if partners(s, 1, 3.0, im_tol)}

out = {"k0_tight": k0(1.5, 0.005), "k0_mid": k0(1.5, 0.01), "k1h_tight": k1_heavy(0.005)}
# size- and charge-matched random control for k0_mid: the falsifier for "any 30% cut helps the q-values"
ids_by_z = collections.defaultdict(list); z_of = {s[0]: s[3] for s in S}
for s in S: ids_by_z[s[3]].append(s[0])
comp = collections.Counter(z_of[i] for i in out["k0_mid"])
rng = random.Random(1); R = set()
for z, c in comp.items(): R.update(rng.sample(ids_by_z[z], c))
out["rand_mid"] = R
for name, R in out.items():
    with open(f"{outdir}/drop_{name}.txt", "w") as f: f.write("".join(f"{i}\n" for i in sorted(R)))
    print(f"{name}\t{len(R)}", flush=True)
