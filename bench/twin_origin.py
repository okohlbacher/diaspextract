#!/usr/bin/env python3
"""Where do same-position twins come from? Classify each k=0 twin pair by the two members' DIA isolation window and
their distance to the nearest fixed RT-cell edge (tile:rt_sec, default 600 s, "a trace eluting across a cut is cut
there"). Three candidate origins, mutually exclusive on the window/RT test:
  CROSS-WINDOW  the two members were assembled in different isolation windows (a precursor in the overlap)
  CELL-CUT      same window, and both members sit within `--cell-tol` s of the same 600-s cell edge
  WITHIN        same window, away from any edge: two precursors from one feature in one pass -- ownership
Usage: twin_origin.py <iso.tsv> <pseudo.mzML> [--im=0.01] [--drt=1.5] [--edges=<csv from spx:tile_boundaries>] [--cell-tol=5] [--sage=results.tsv]
"""
import sys, re, bisect, collections
ISO, PPM = 1.0033548, 20.0
iso_p, mz_p = sys.argv[1], sys.argv[2]
kw = dict((a.lstrip("-").split("=", 1) + ["1"])[:2] for a in sys.argv[3:])
IM_TOL, RT_TOL = float(kw.get("im", 0.01)), float(kw.get("drt", 1.5))
EDGES = [float(x) for x in kw.get("edges", "").split(",") if x]   # the header's spx:tile_boundaries,
CELL_TOL = float(kw.get("cell-tol", 5))                           # not multiples of the nominal pitch

IDX = re.compile(rb'<spectrum [^>]*index="(\d+)"')
TGT = re.compile(rb'MS:1000827"[^>]*value="([^"]+)"')
LO = re.compile(rb'MS:1000828"[^>]*value="([^"]+)"')
HI = re.compile(rb'MS:1000829"[^>]*value="([^"]+)"')
win = {}
cur = None
with open(mz_p, "rb") as f:                       # one streaming pass; the file is 8 GB
    t = lo = hi = None
    for line in f:
        m = IDX.search(line)
        if m: cur, t, lo, hi = int(m.group(1)), None, None, None; continue
        if cur is None: continue
        m = TGT.search(line)
        if m: t = float(m.group(1)); continue
        m = LO.search(line)
        if m: lo = float(m.group(1)); continue
        m = HI.search(line)
        if m:
            hi = float(m.group(1))
            win[cur] = (round(t - lo, 3), round(t + hi, 3)); cur = None
print(f"windows read: {len(win):,} spectra, {len(set(win.values())):,} distinct isolation windows", flush=True)
edges = sorted({w for pair in set(win.values()) for w in pair})
print("window edges (first 20):", ["%.1f" % e for e in edges[:20]], flush=True)

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
pep = {}
if "sage" in kw:
    import csv
    for r in csv.DictReader(open(kw["sage"], newline=""), delimiter="\t"):
        if r["label"] == "1" and float(r["spectrum_q"]) <= 0.01: pep[int(r["scannr"].split("=")[-1])] = r["peptide"]

near_edge = lambda rt: any(abs(rt - e) <= CELL_TOL for e in EDGES)
cls = collections.Counter(); drt = collections.Counter(); ident = collections.Counter(); seen = set()
same_win_overlap = collections.Counter()
for s in S:
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: continue
    tol = mz * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    for j in range(bisect.bisect_left(keys, mz - tol), bisect.bisect_right(keys, mz + tol)):
        o = arr[j]
        if o[0] <= idx or abs(o[1] - rt) > RT_TOL or abs(o[4] - im) > IM_TOL: continue
        key = (idx, o[0])
        if key in seen: continue
        seen.add(key)
        wa, wb = win.get(idx), win.get(o[0])
        if wa is None or wb is None: c = "NO-WINDOW"
        elif wa != wb:
            c = "CROSS-WINDOW"
            same_win_overlap[(wa, wb)] += 1
        elif near_edge(rt) and near_edge(o[1]): c = "CELL-CUT"
        else: c = "WITHIN"
        cls[c] += 1
        drt[(c, round(abs(o[1] - rt), 1))] += 1
        pa, pb = pep.get(idx), pep.get(o[0])
        ident[(c, "both-same" if pa and pa == pb else "both-diff" if pa and pb else "one" if pa or pb else "neither")] += 1

tot = sum(cls.values())
print(f"\n{tot:,} k=0 twin PAIRS (same charge, {PPM:.0f} ppm, dRT<={RT_TOL}, dIM<={IM_TOL})\n")
for c, n in cls.most_common():
    print(f"  {c:<14}{n:>9,}  {100*n/tot:>5.1f}%")
if pep:
    print("\nidentification of the two members:")
    for c in cls:
        sub = {k[1]: v for k, v in ident.items() if k[0] == c}; t = sum(sub.values()) or 1
        print(f"  {c:<14}" + "  ".join(f"{k} {v:,} ({100*v/t:.0f}%)" for k, v in sorted(sub.items())))
print("\ndRT (s) within each class:")
for c in cls:
    sub = sorted((k[1], v) for k, v in drt.items() if k[0] == c)[:8]
    print(f"  {c:<14}" + "  ".join(f"{d}s:{v:,}" for d, v in sub))
if same_win_overlap:
    print("\ntop CROSS-WINDOW window pairs:")
    for (wa, wb), n in same_win_overlap.most_common(6):
        print(f"  [{wa[0]}, {wa[1]}] x [{wb[0]}, {wb[1]}]  {n:,}   overlap = {min(wa[1],wb[1]) - max(wa[0],wb[0]):.2f} m/z")
