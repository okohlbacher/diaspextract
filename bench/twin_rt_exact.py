#!/usr/bin/env python3
"""Do twins share their MS1 frame EXACTLY? The in-tool fold (1.4.0, mergeTwins_) links only spectra with bit-equal RT;
bench/twin_merge.py, which priced the fold, linked pairs up to 1.5 s apart. This counts, among the pairs the offline
relation links (same charge, 20 ppm of the earlier member, |dRT| <= 1.5 s, |d1/K0| <= tol), those whose RT is
text-identical in the mzML and those that differ, with the RT offsets of the latter. All identical => the two
relations are the same on that file.
Usage: twin_rt_exact.py <iso.tsv>        (iso.tsv from bench/iso_dup.py extract)
"""
import sys, bisect, collections
PPM = 20e-6
S = []
with open(sys.argv[1]) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        S.append((int(p[0]), p[1], float(p[1]), float(p[2]), int(p[3] or 0), float(p[4] or 0)))
byz = collections.defaultdict(list)
for s in S: byz[s[4]].append(s)
for z in byz: byz[z].sort(key=lambda s: s[3])
keys = {z: [s[3] for s in v] for z, v in byz.items()}
for im_tol in (0.005, 0.01):
    same = diff = 0; dts = collections.Counter()
    for s in S:
        if s[4] <= 0: continue
        tol = s[3] * PPM; arr = byz[s[4]]; k = keys[s[4]]
        for j in range(bisect.bisect_left(k, s[3] - tol), bisect.bisect_right(k, s[3] + tol)):
            o = arr[j]
            if o[0] <= s[0] or abs(o[2] - s[2]) > 1.5 or abs(o[5] - s[5]) > im_tol: continue
            if o[1] == s[1]: same += 1
            else: diff += 1; dts[round(abs(o[2] - s[2]), 3)] += 1
    print(f"im<={im_tol}: pairs with text-identical RT {same:,}; different RT {diff:,}  {dict(dts.most_common(5))}")
