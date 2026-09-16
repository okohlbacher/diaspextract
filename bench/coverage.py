#!/usr/bin/env python3
"""What do we LOSE vs the reference implementation? Properties of the coverage gap, on one raw acquisition.

This is the clean comparison the project has been missing. the reference implementation and DIAspeXtract both extract
pseudo-spectra from the SAME dataset D .d, and both were searched with the SAME sage_deiso.json against
the SAME fasta. So the peptide sets differ ONLY by extraction -- no second injection, no RT drift,
no DDA-conditioning. Precursor RT / m/z / ion mobility are directly comparable between the two.

  MISS   = the reference implementation peptides we do NOT identify   <- the coverage gap being characterised
  SHARED = identified by both                       <- the reference distribution

For every the reference implementation peptide we take its BEST PSM and read properties from two places:
  * Sage      : length, charge, m/z, RT, ion mobility, hyperscore, matched_peaks,
                matched_intensity_pct, ms2_intensity, missed cleavages, isotope error
  * the reference implementation : the emitting spectrum's PRECURSOR PEAK INTENSITY, TIC and base-peak intensity
                -- a direct abundance proxy in the original data, available for MISS and SHARED
                alike (our own mzML does not write precursor intensity, so the reference implementation's is the
                only intensity measure defined for BOTH classes; that is why it is used).

Then the class A/B split, which is finally well-posed here (same acquisition => same RT axis):
  A NOT EMITTED   no emitted spectrum of ours at that precursor coordinate  -> extraction gap
  B EMITTED, NOT IDENTIFIED                                                 -> assembly/ID gap

CHANCE FLOOR, mandatory: matching a coordinate inside a tolerance window succeeds sometimes by
coincidence. Every A/B number is reported beside a decoy-coordinate control (same RT/IM, m/z
shifted off by a non-physical offset), so "we emitted something there" cannot be read as real
when it is coincidence. This project retracted the MS1 funnel for exactly that omission.
"""
import csv, sys, bisect, statistics as st

OURS_MZML = sys.argv[1]     # our base pseudo.mzML
OURS_TSV  = sys.argv[2]     # our base results.sage.tsv
DT_MZML   = sys.argv[3]     # the reference implementation mzML
DT_TSV    = sys.argv[4]     # the reference implementation results.sage.tsv
PROTON = 1.00727646
PPM = 20.0


def best_psm_per_peptide(tsv):
    """Target PSMs at peptide_q<=0.01, best-scoring row per distinct peptide."""
    best = {}
    with open(tsv) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            try:
                if float(x["peptide_q"]) > 0.01:
                    continue
                prots = [p for p in x["proteins"].split(";") if p]
                if prots and all(p.startswith("rev_") for p in prots):
                    continue
                s = float(x["sage_discriminant_score"])
                pep = x["peptide"]
                if pep in best and best[pep]["s"] >= s:
                    continue
                z = int(x["charge"])
                best[pep] = {"s": s, "z": z, "len": int(x["peptide_len"]),
                             "mz": (float(x["expmass"]) + z * PROTON) / z,
                             "rt": float(x["rt"]), "im": float(x.get("ion_mobility") or 0),
                             "hyper": float(x["hyperscore"]),
                             "mp": float(x["matched_peaks"]),
                             "mip": float(x.get("matched_intensity_pct") or 0),
                             "ms2i": float(x.get("ms2_intensity") or 0),
                             "mc": int(x.get("missed_cleavages") or 0),
                             "iso": float(x.get("isotope_error") or 0),
                             "scan": x["scannr"]}
            except (ValueError, KeyError):
                continue
    return best


def spectrum_index(mzml, id_tag):
    """scannr-string -> (rt, mz, im, prec_intensity, tic). Streams; keeps only metadata."""
    out = {}
    cur = rt = mz = im = pint = tic = None
    with open(mzml) as f:
        for line in f:
            if "<spectrum " in line and id_tag in line:
                if cur is not None:
                    out[cur] = (rt, mz, im, pint, tic)
                i = line.index(id_tag) + len('id="')   # start of the native id itself
                j = line.index('"', i)
                cur = line[i:j]                        # e.g. 'frame=1234' == Sage scannr
                rt = mz = im = pint = tic = None
            elif cur is not None:
                if tic is None and 'name="total ion current"' in line:
                    k = line.index('value="') + 7; tic = float(line[k:line.index('"', k)])
                elif mz is None and 'name="selected ion m/z"' in line:
                    k = line.index('value="') + 7; mz = float(line[k:line.index('"', k)])
                elif rt is None and 'name="scan start time"' in line:
                    k = line.index('value="') + 7; rt = float(line[k:line.index('"', k)])
                elif im is None and 'name="inverse reduced ion mobility"' in line:
                    k = line.index('value="') + 7; im = float(line[k:line.index('"', k)])
                elif pint is None and 'name="peak intensity"' in line:
                    k = line.index('value="') + 7; pint = float(line[k:line.index('"', k)])
        if cur is not None:
            out[cur] = (rt, mz, im, pint, tic)
    return out


print("[load] sage results ...", flush=True)
ours = best_psm_per_peptide(OURS_TSV)
dt = best_psm_per_peptide(DT_TSV)
MISS = {p: v for p, v in dt.items() if p not in ours}
SHARED = {p: v for p, v in dt.items() if p in ours}
ONLY = {p: v for p, v in ours.items() if p not in dt}
print("[sets] the reference implementation %d, ours %d | SHARED %d, MISS %d, ours-only %d"
      % (len(dt), len(ours), len(SHARED), len(MISS), len(ONLY)), flush=True)

print("[load] the reference implementation spectra (precursor intensity) ...", flush=True)
dt_spec = spectrum_index(DT_MZML, 'id="frame=')
print("[load]   %d the reference implementation spectra" % len(dt_spec), flush=True)

# attach the reference implementation-side intensity to every the reference implementation peptide
for d in (MISS, SHARED):
    for p, v in d.items():
        s = dt_spec.get(v["scan"])
        v["pint"] = s[3] if s and s[3] is not None else None
        v["tic"] = s[4] if s and s[4] is not None else None


def cmp_prop(label, key, logscale=False):
    m = [v[key] for v in MISS.values() if v.get(key) is not None]
    h = [v[key] for v in SHARED.values() if v.get(key) is not None]
    if not m or not h:
        print("  %-24s (no data)" % label); return
    mm, hh = st.median(m), st.median(h)
    q = lambda v, p: sorted(v)[min(int(p * len(v)), len(v) - 1)]
    print("  %-24s MISS med %10.4g [%.4g-%.4g]   SHARED med %10.4g [%.4g-%.4g]   ratio %.2f"
          % (label, mm, q(m, .25), q(m, .75), hh, q(h, .25), q(h, .75),
             (mm / hh) if hh else float("nan")))


print("\n=== PROPERTIES: peptides the reference implementation finds and we MISS, vs those we SHARE ===")
cmp_prop("precursor intensity", "pint")
cmp_prop("spectrum TIC", "tic")
cmp_prop("ms2_intensity (Sage)", "ms2i")
cmp_prop("matched_intensity_pct", "mip")
cmp_prop("hyperscore (the reference implementation)", "hyper")
cmp_prop("matched_peaks", "mp")
cmp_prop("peptide length", "len")
cmp_prop("precursor m/z", "mz")
cmp_prop("charge", "z")
cmp_prop("ion mobility 1/K0", "im")
cmp_prop("RT (min)", "rt")
cmp_prop("missed cleavages", "mc")


def dist(label, key, bins):
    print("\n--- %s: MISS%% vs SHARED%% (enrichment) ---" % label)
    m = [v[key] for v in MISS.values() if v.get(key) is not None]
    h = [v[key] for v in SHARED.values() if v.get(key) is not None]
    print("  %-18s %8s %8s %9s" % ("bin", "MISS%", "SHARED%", "enrich"))
    for lo, hi in bins:
        fm = 100.0 * sum(1 for x in m if lo <= x < hi) / max(len(m), 1)
        fh = 100.0 * sum(1 for x in h if lo <= x < hi) / max(len(h), 1)
        e = (fm / fh) if fh > 0 else float("inf")
        print("  %-18s %7.1f%% %7.1f%% %8.2fx" % ("%g-%g" % (lo, hi), fm, fh, e))


dist("precursor intensity", "pint", [(0, 500), (500, 1000), (1000, 2500), (2500, 1e4), (1e4, 1e9)])
dist("peptide length", "len", [(7, 10), (10, 15), (15, 20), (20, 45)])
dist("precursor m/z", "mz", [(300, 500), (500, 700), (700, 900), (900, 1200), (1200, 2000)])
dist("charge", "z", [(1, 2), (2, 3), (3, 4), (4, 9)])
dist("hyperscore", "hyper", [(0, 15), (15, 25), (25, 40), (40, 200)])

# ---------------------------------------------------------------- class A/B split (well-posed)
print("\n[load] our emitted population ...", flush=True)
pop = []
cur = rt = mz = im = None
with open(OURS_MZML) as f:
    for line in f:
        if "<spectrum " in line and 'id="spectrum=' in line:
            if rt is not None and mz is not None:
                pop.append((mz, rt, im))
            rt = mz = im = None
            cur = 1
        elif cur is not None:
            if mz is None and 'name="selected ion m/z"' in line:
                k = line.index('value="') + 7; mz = float(line[k:line.index('"', k)])
            elif rt is None and 'name="scan start time"' in line:
                k = line.index('value="') + 7; rt = float(line[k:line.index('"', k)])
                # FAILURE 8: unit from the CV accessor, not the magnitude
                if 'unitName="second"' in line or 'UO:0000010' in line:
                    rt /= 60.0
            elif im is None and 'name="inverse reduced ion mobility"' in line:
                k = line.index('value="') + 7; im = float(line[k:line.index('"', k)])
    if rt is not None and mz is not None:
        pop.append((mz, rt, im))
pop.sort()
pmz = [p[0] for p in pop]
print("[load]   %d emitted spectra of ours" % len(pop), flush=True)
_prt = [p[1] for p in pop]
_drt = [v["rt"] for v in dt.values()]
print("[units] our emitted RT %.2f-%.2f min | the reference implementation peptide RT %.2f-%.2f min"
      % (min(_prt), max(_prt), min(_drt), max(_drt)), flush=True)
if not (min(_prt) <= max(_drt) and min(_drt) <= max(_prt) and max(_prt) < 5 * max(_drt)):
    sys.exit("ABORT: RT ranges do not overlap -- unit mismatch (FAILURE 8). Refusing to report.")


def emitted_near(mz, rt, im, rt_win, im_win=0.05):
    tol = mz * PPM * 1e-6
    lo = bisect.bisect_left(pmz, mz - tol)
    hi = bisect.bisect_right(pmz, mz + tol)
    for k in range(lo, hi):
        _, ert, eim = pop[k]
        if rt_win is not None and abs(ert - rt) > rt_win:
            continue
        if im_win is not None and im and eim and abs(eim - im) > im_win:
            continue
        return True
    return False


print("\n=== CLASS A/B: for peptides we MISS, did we emit a spectrum there at all? ===")
print("(same acquisition, so RT is directly comparable -- no drift confound)")
print("  %-12s %10s %10s %10s   %s" % ("RT window", "A not emit", "B emitted", "A %", "decoy-floor B%"))
for rw in (0.2, 0.5, 1.0, None):
    a = b = 0
    fb = 0
    for v in MISS.values():
        if emitted_near(v["mz"], v["rt"], v["im"], rw):
            b += 1
        else:
            a += 1
        # chance floor: same RT/IM, physically impossible m/z offset
        if emitted_near(v["mz"] + 7.3, v["rt"], v["im"], rw):
            fb += 1
    n = max(a + b, 1)
    lbl = "any RT" if rw is None else "%.1f min" % rw
    print("  %-12s %10d %10d %9.1f%%   %13.1f%%" % (lbl, a, b, 100.0 * a / n, 100.0 * fb / n))

# same split for SHARED as a positive control: these we DID identify, so B must dominate
print("\n  positive control (SHARED peptides, we identified them -> B must dominate):")
for rw in (0.5,):
    a = b = 0
    for v in SHARED.values():
        if emitted_near(v["mz"], v["rt"], v["im"], rw):
            b += 1
        else:
            a += 1
    print("  %-12s %10d %10d %9.1f%%" % ("%.1f min" % rw, a, b, 100.0 * a / max(a + b, 1)))
