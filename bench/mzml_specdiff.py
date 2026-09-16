#!/usr/bin/env python3
"""Per-spectrum diff of two DIAspeXtract mzML outputs, aligned by (RT, precursor m/z, charge).

    mzml_specdiff.py A.mzML B.mzML [--tile SEC]

The semantic digest says "differs"; this says WHERE and HOW MUCH: spectra only in A / only in B,
identical spectra, spectra whose intensities differ (max relative difference binned), spectra whose
m/z lists differ, spectra with a different peak count. With --tile SEC every differing spectrum is
also binned by the distance of its RT to the nearest multiple of SEC (a tile boundary), which is the
direct test of "the tiling only changes spectra near the boundaries". Pure Python, streaming: the
base64 payloads are compared as bytes first and decoded only when they differ.
"""
import base64, re, struct, sys, zlib

RT = re.compile(rb'accession="MS:1000016"[^>]*value="([^"]+)"')
PMZ = re.compile(rb'accession="MS:1000744"[^>]*value="([^"]+)"')
CHG = re.compile(rb'accession="MS:1000041"[^>]*value="([^"]+)"')
NPK = re.compile(rb'defaultArrayLength="(\d+)"')
BIN = re.compile(rb'<binary>([^<]*)</binary>')
COMP = re.compile(rb'accession="MS:1000574"')   # zlib compression

def spectra(path):
    with open(path, "rb") as f:
        block = None
        for line in f:
            s = line.lstrip()
            if block is None:
                if s.startswith(b"<spectrum "): block = [line]
                elif s.startswith(b"</spectrumList>"): return
                continue
            block.append(line)
            if s.startswith(b"</spectrum>"):
                data = b"".join(block); block = None
                rt = float(RT.search(data).group(1))
                m = PMZ.search(data); pmz = float(m.group(1)) if m else 0.0
                m = CHG.search(data); z = int(m.group(1)) if m else 0
                n = int(NPK.search(data).group(1))
                bins = BIN.findall(data)
                yield (rt, pmz, z), n, bins, bool(COMP.search(data))

def decode(b64, comp, fmt):
    raw = base64.b64decode(b64)
    if comp: raw = zlib.decompress(raw)
    return struct.unpack("<%d%s" % (len(raw) // struct.calcsize(fmt), fmt), raw)

PPM = 10.0

def match_ppm(mz1, i1, mz2, i2, ppm=PPM):
    """Greedy merge of two m/z-sorted peak lists within `ppm`: the apex m/z of the SAME fragment
    moves with its apex frame's calibration factor, so an exact comparison counts an apex shift
    as a lost-plus-gained fragment. Returns (matched, a_only, b_only, matched with intensity
    within 0.1%)."""
    i = j = 0; m = ao = bo = same = 0
    while i < len(mz1) and j < len(mz2):
        a, b = mz1[i], mz2[j]
        if abs(a - b) <= a * ppm * 1e-6:
            m += 1
            if abs(i1[i] - i2[j]) <= 1e-3 * max(abs(i1[i]), 1e-30): same += 1
            i += 1; j += 1
        elif a < b: ao += 1; i += 1
        else: bo += 1; j += 1
    return m, ao + (len(mz1) - i), bo + (len(mz2) - j), same

def main(a, b, tile):
    ga, gb = spectra(a), spectra(b)
    na = nb = 0
    ka = next(ga, None); kb = next(gb, None)
    stats = {"identical": 0, "intensity_only": 0, "mz_differs": 0, "count_differs": 0, "only_A": 0, "only_B": 0}
    frag = {"common": 0, "a_only": 0, "b_only": 0}   # fragment-level, exact m/z, over spectra present in both
    tol = {"matched": 0, "a_only": 0, "b_only": 0, "int_same": 0}   # fragment-level within PPM
    relbins = {"<=1e-6": 0, "<=1e-4": 0, "<=1e-2": 0, ">1e-2": 0}
    dist = {"<=30s": 0, "<=63s": 0, "<=123s": 0, ">123s": 0}
    dist_all = dict(dist)
    def bin_dist(rt, d):
        if tile:
            r = rt % tile; x = min(r, tile - r)
            k = "<=30s" if x <= 30 else "<=63s" if x <= 63 else "<=123s" if x <= 123 else ">123s"
            d[k] += 1
    while ka or kb:
        if kb is None or (ka is not None and ka[0] < kb[0]):
            stats["only_A"] += 1; bin_dist(ka[0][0], dist); bin_dist(ka[0][0], dist_all); ka = next(ga, None); continue
        if ka is None or kb[0] < ka[0]:
            stats["only_B"] += 1; bin_dist(kb[0][0], dist); bin_dist(kb[0][0], dist_all); kb = next(gb, None); continue
        # same key
        (key, n1, b1, c1), (_, n2, b2, c2) = ka, kb
        bin_dist(key[0], dist_all)
        if b1 == b2 and c1 == c2 and n1 == n2:
            stats["identical"] += 1; frag["common"] += n1; tol["matched"] += n1; tol["int_same"] += n1
        elif n1 != n2 or decode(b1[0], c1, "d") != decode(b2[0], c2, "d"):
            mz1, mz2 = decode(b1[0], c1, "d"), decode(b2[0], c2, "d")
            i1, i2 = decode(b1[1], c1, "f"), decode(b2[1], c2, "f")
            s1, s2 = set(round(x, 4) for x in mz1), set(round(x, 4) for x in mz2)
            frag["common"] += len(s1 & s2); frag["a_only"] += len(s1 - s2); frag["b_only"] += len(s2 - s1)
            m, ao, bo, same = match_ppm(mz1, i1, mz2, i2)
            tol["matched"] += m; tol["a_only"] += ao; tol["b_only"] += bo; tol["int_same"] += same
            stats["count_differs" if n1 != n2 else "mz_differs"] += 1; bin_dist(key[0], dist)
        else:
            i1, i2 = decode(b1[1], c1, "f"), decode(b2[1], c2, "f")
            mx = max((abs(x - y) / max(abs(x), 1e-30) for x, y in zip(i1, i2)), default=0.0)
            stats["intensity_only"] += 1; bin_dist(key[0], dist); frag["common"] += len(i1); tol["matched"] += len(i1)
            tol["int_same"] += sum(1 for x, y in zip(i1, i2) if abs(x - y) <= 1e-3 * max(abs(x), 1e-30))
            relbins["<=1e-6" if mx <= 1e-6 else "<=1e-4" if mx <= 1e-4 else "<=1e-2" if mx <= 1e-2 else ">1e-2"] += 1
        ka = next(ga, None); kb = next(gb, None)
    tot = sum(stats.values())
    print(f"spectra compared: {tot}")
    for k, v in stats.items(): print(f"  {k:16s} {v:9d}  ({100.0 * v / max(tot, 1):.4f}%)")
    print("intensity-only differences, max relative:", relbins)
    tot_f = sum(frag.values())
    print(f"fragment-level over spectra present in both, exact m/z: common {frag['common']} ({100.0 * frag['common'] / max(tot_f, 1):.2f}% of the union), A-only {frag['a_only']}, B-only {frag['b_only']}")
    tot_t = tol["matched"] + tol["a_only"] + tol["b_only"]
    print(f"fragment-level within {PPM:g} ppm: matched {tol['matched']} ({100.0 * tol['matched'] / max(tot_t, 1):.2f}% of the union), of which intensity within 0.1%: {tol['int_same']} ({100.0 * tol['int_same'] / max(tol['matched'], 1):.2f}%); A-only {tol['a_only']}, B-only {tol['b_only']}")
    if tile:
        print(f"differing spectra by distance to the nearest {tile}-s boundary:", dist)
        print("all spectra by that distance:                      ", dist_all)

if __name__ == "__main__":
    args = sys.argv[1:]
    tile = 0.0
    if "--tile" in args:
        i = args.index("--tile"); tile = float(args[i + 1]); del args[i:i + 2]
    if len(args) != 2: raise SystemExit(__doc__)
    main(args[0], args[1], tile)
