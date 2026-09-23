#!/usr/bin/env python3
"""Merge same-position twins instead of dropping one: the fix that keeps the second copy's fragments.

Two precursors at the same inferred monoisotope (src/diaspextract.cpp:3038) are one precursor claimed twice, but
their spectra are NOT copies (peak counts differ in 69% of pairs, median 30 peaks). So the cluster's survivor is
rewritten with the UNION of the cluster's fragments -- peaks within --tol ppm of one another are one fragment,
intensity = max, m/z = intensity-weighted mean -- and the other members are deleted. Everything else is untouched;
`id="spectrum=N"` is kept so results join back to the unfiltered run.

    twin_merge.py <iso.tsv> <src.mzML> <dst.mzML> [--im=0.01] [--drt=1.5] [--tol=10] [--stats=path] [--drop-only]
"""
import sys, re, base64, struct, bisect, collections

ISO, PPM = 1.0033548, 20.0
iso_p, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
kw = dict((a.lstrip("-").split("=", 1) + ["1"])[:2] for a in sys.argv[4:])   # bare flags allowed
IM_TOL, RT_TOL, TOL = float(kw.get("im", 0.01)), float(kw.get("drt", 1.5)), float(kw.get("tol", 10)) * 1e-6

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

parent = {}
def find(a):
    while parent.get(a, a) != a: parent[a] = parent.get(parent[a], parent[a]); a = parent[a]
    return a
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

for s in S:                                            # cluster by the same relation the price list used
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: continue
    tol = mz * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    for j in range(bisect.bisect_left(keys, mz - tol), bisect.bisect_right(keys, mz + tol)):
        o = arr[j]
        if o[0] <= idx or abs(o[1] - rt) > RT_TOL or abs(o[4] - im) > IM_TOL: continue
        parent.setdefault(idx, idx); parent.setdefault(o[0], o[0]); union(idx, o[0])

info = {s[0]: s for s in S}
clusters = collections.defaultdict(list)
for i in list(parent): clusters[find(i)].append(i)
survivor, merge_into = {}, {}                          # strongest = most isotopes, then most peaks, then lowest index
for root, mem in clusters.items():
    best = max(mem, key=lambda i: (info[i][6], info[i][5], -i))
    survivor[best] = mem
    for i in mem:
        if i != best: merge_into[i] = best
print(f"{len(clusters):,} clusters over {len(parent):,} spectra; {len(merge_into):,} merged away, {len(survivor):,} survivors",
      flush=True)

IDX = re.compile(rb'<spectrum [^>]*index="(\d+)"')
DAL = re.compile(rb'(defaultArrayLength=")(\d+)(")')
ENC = re.compile(rb'(encodedLength=")(\d+)(")')
COUNT = re.compile(rb'(<spectrumList count=")(\d+)(")')

def decode(b64, width):
    raw = base64.b64decode(b64)
    return struct.unpack("<%d%s" % (len(raw) // width, "d" if width == 8 else "f"), raw)

def peaks_of(block):
    """(m/z, intensity) of one <spectrum> block; the writer emits uncompressed 64-bit m/z then 32-bit intensity."""
    kind, width, mz, it = None, 8, None, None
    for line in block:
        if b'name="m/z array"' in line: kind = "mz"
        elif b'name="intensity array"' in line: kind = "it"
        elif b'name="64-bit float"' in line: width = 8
        elif b'name="32-bit float"' in line: width = 4
        elif b"<binary>" in line:
            a, b = line.find(b"<binary>") + 8, line.find(b"</binary>")
            v = decode(line[a:b], width) if b > a else ()
            if kind == "mz": mz = v
            elif kind == "it": it = v
            kind = None
    return list(zip(mz or (), it or ()))

# pass 1: the fragments every cluster contributes, accumulated on the survivor
pool = collections.defaultdict(list)
need = set(merge_into) | set(survivor)
with open(src, "rb") as f:
    block, cur = None, -1
    for line in f:
        s = line.lstrip()
        if s.startswith(b"<spectrum "):
            cur = int(IDX.search(line).group(1)); block = [line] if cur in need else None
        elif block is not None:
            block.append(line)
            if s.startswith(b"</spectrum>"):
                pool[merge_into.get(cur, cur)].extend(peaks_of(block))
                block = None
print(f"pass 1 done: {len(pool):,} survivors carry {sum(len(v) for v in pool.values()):,} raw peaks", flush=True)

def fuse(ps):
    """Peaks within TOL of the running group's m/z are one fragment: intensity = max, m/z = intensity-weighted mean."""
    ps.sort()
    out, gm, gi, gw = [], 0.0, 0.0, 0.0
    for mz, it in ps:
        if gw and mz - gm / gw > (gm / gw) * TOL:
            out.append((gm / gw, gi)); gm = gi = gw = 0.0
        gm += mz * it; gw += it; gi = max(gi, it)
    if gw: out.append((gm / gw, gi))
    return out

merged = {} if "drop-only" in kw else {i: fuse(p) for i, p in pool.items() if i in survivor}
if "drop-only" in kw: print("drop-only: survivors keep their own peaks", flush=True)
total = 0
with open(src, "rb") as f:
    for line in f:
        if line.lstrip().startswith(b"<spectrum "):
            total += int(IDX.search(line).group(1)) not in merge_into

n_out = 0
with open(src, "rb") as f, open(dst, "wb") as out:
    block, cur, in_spec = None, -1, False
    for line in f:
        s = line.lstrip()
        if s.startswith(b"<spectrum "):
            cur, in_spec = int(IDX.search(line).group(1)), True
            block = [line]; continue
        if in_spec:
            block.append(line)
            if not s.startswith(b"</spectrum>"): continue
            in_spec = False
            if cur in merge_into: continue
            p = merged.get(cur)
            if p is not None:
                mzb = base64.b64encode(struct.pack("<%dd" % len(p), *[x[0] for x in p]))
                itb = base64.b64encode(struct.pack("<%df" % len(p), *[x[1] for x in p]))
                new, kind = [], None
                for L in block:
                    if b"<spectrum " in L:
                        L = DAL.sub(lambda m: m.group(1) + str(len(p)).encode() + m.group(3), L, count=1)
                    elif b'name="m/z array"' in L: kind = "mz"
                    elif b'name="intensity array"' in L: kind = "it"
                    elif b"<binary>" in L:
                        b64 = mzb if kind == "mz" else itb
                        a = L.find(b"<binary>")
                        L = L[:a] + b"<binary>" + b64 + b"</binary>\n"
                        # the encodedLength attribute sits on the enclosing binaryDataArray, already emitted
                        for k in range(len(new) - 1, -1, -1):
                            if b"encodedLength=" in new[k]:
                                new[k] = ENC.sub(lambda m: m.group(1) + str(len(b64)).encode() + m.group(3), new[k], count=1); break
                        kind = None
                    new.append(L)
                block = new
            block[0] = re.sub(rb'(index=")(\d+)(")', lambda m: m.group(1) + str(n_out).encode() + m.group(3), block[0], count=1)
            out.write(b"".join(block)); n_out += 1
            continue
        if s.startswith(b"<indexedmzML"): continue
        if s.startswith(b"<indexList"): break
        if b"<spectrumList" in line:
            line = COUNT.sub(lambda m: m.group(1) + str(total).encode() + m.group(3), line)
        out.write(line)
assert n_out == total, (n_out, total)
if "stats" in kw:
    with open(kw["stats"], "w") as g:
        for i, p in sorted(merged.items()): g.write("%d\t%d\t%d\n" % (i, info[i][5], len(p)))
print(f"wrote {n_out:,} spectra ({len(merge_into):,} merged away, {len(merged):,} rewritten) -> {dst}", flush=True)
