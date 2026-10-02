#!/usr/bin/env python3
"""Entrapment-controlled OPEN-search benchmark.

Why this exists
---------------
Every number in README.md comes from a CLOSED Sage search scored by peptide_q against a DIA-NN
reference whose space is tryptic / z2-4 / <=2 mods -- structurally incapable of containing the
class of thing this tool exists to find. We emit 4.6 PSMs/peptide against the reference implementation's 1.46, so
redundancy is free score under that metric. That single fact predicts, in advance, every result
observed: each redundancy-reducing lever losing, and ms1split winning while simultaneously
raising co-isolation 52% and dropping charge agreement below baseline.

Six falsifications sharing one biased instrument are one falsification.

Three measurements, in increasing order of how much they depend on assumptions
-----------------------------------------------------------------------------
1. DELTA-MASS GEOMETRY -- needs no ground truth, no decoy, no entrapment.
   A charge misassignment fabricates

       dM = (z' - z) * (m/z - PROTON)

   which is LINEAR IN m/z. In the 2-D (precursor m/z, dM) plane those PSMs fall on RAYS through
   (PROTON, 0) with integer slopes +-1, +-2, ... A real modification is a constant mass shift and
   forms a HORIZONTAL BAND. The two are geometrically separable by construction, so the sloped-ray
   mass fraction measures charge error IN THE ACTUAL USE CASE -- strictly better than the
   DIA-NN-anchored proxy, whose decoy floor reached 91.9%.

2. ENTRAPMENT FDR -- needs a foreign proteome but no target-decoy assumption.
   Open search grows the candidate space ~10^3x and decoy-based FDR degrades accordingly
   (Kong 2017; Chick 2015). Entrapment is the accepted validation when the search space changes
   (Wen, Freestone, Noble, Kall, Nesvizhskii, Nat Methods 2025).

3. MODIFIED-PEPTIDE COUNT at <=1% entrapment FDR -- the primary metric.

Pre-registered falsification (written before the first run)
-----------------------------------------------------------
Closed-search ratio split_count/base on IH3 = 8411/5817 = 1.446.

  * modified-peptide ratio >= 1.446  -> "the closed benchmark rewards shotgunning" is FALSIFIED.
                                        The defaults are vindicated for BOTH purposes, the six
                                        falsified levers stand, redundancy is load-bearing, and
                                        we stop trying to reduce it.
  * ratio materially < 1.0           -> the shipped defaults are ACTIVELY HARMFUL for the stated
                                        purpose and must be reconsidered.
  * 1.0 <= ratio < 1.446             -> real but overstated by the closed metric.

Replicate spread is 0.17% (4 runs of an identical effective config: 8411/8408/8406/8397), so a
ratio difference below ~0.2% is noise.

POSITIVE CONTROL, non-negotiable. The envelope arm has 10,199 known 4->2 confusions. They MUST
appear as slope -2 ray mass in sage_open_wide.json, and must be markedly weaker in the count arm.
If the control does not fire, the instrument is not measuring charge error and NO other number
from this script may be quoted -- that is the same mistake as reporting a confident result from an
assay whose control was never run.
"""
import argparse, collections, csv, json, math, os, subprocess, sys
from pathlib import Path

PROTON = 1.00727646
HERE = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------- entrapment FASTA
def build_entrapment(target_fa, foreign_fa, out_fa, tag="ENTRAP_"):
    """Concatenate a foreign proteome, prefixing every foreign accession.

    Returns (n_target, n_entrap). Any PSM to a tagged protein is a KNOWN false positive: the
    sample is human, so an Arabidopsis peptide cannot be genuinely present.
    """
    n_t = n_e = 0
    with open(out_fa, "w") as out:
        with open(target_fa) as f:
            for line in f:
                if line.startswith(">"):
                    n_t += 1
                out.write(line)
        with open(foreign_fa) as f:
            for line in f:
                if line.startswith(">"):
                    n_e += 1
                    out.write(">" + tag + line[1:])
                else:
                    out.write(line)
    return n_t, n_e


def entrapment_fdr(rows, tag="ENTRAP_", ratio=1.0):
    """Empirical FDR from entrapment hits, independent of the decoy model.

    With target and entrapment databases of comparable size, an incorrect PSM is equally likely
    to land in either, so FP ~ 2 * n_entrap and FDR ~ 2 * n_entrap / n_total. `ratio` =
    n_entrap_proteins / n_target_proteins corrects for unequal sizes.

    Returns a score threshold achieving <=1% FDR, or None if unreachable -- None must propagate
    as "no threshold exists", never as a silently permissive one.
    """
    scored = sorted(rows, key=lambda r: -r["score"])
    n_e = n_t = 0
    thresh = None
    for r in scored:
        if r["entrap"]:
            n_e += 1
        else:
            n_t += 1
        if n_t == 0:
            continue
        fdr = (n_e / ratio) * 2.0 / (n_e + n_t)
        if fdr <= 0.01:
            thresh = r["score"]
    return thresh


# ------------------------------------------------------------------- delta-mass ray geometry
def ray_decomposition(psms, max_slope=4, band_tol=0.02, ray_tol_ppm=None, ray_tol_da=0.5):
    """Split delta-mass density into CHARGE-ERROR RAYS and MODIFICATION BANDS.

    A charge error z -> z' gives  dM = (z' - z)(mz - PROTON), so in (mz, dM) it is a straight
    line through (PROTON, 0) with INTEGER slope k = z' - z. A modification gives dM = const, a
    horizontal band. Assignment is by perpendicular-ish residual in dM:

        ray k   : |dM - k*(mz - PROTON)| <= ray_tol_da
        band    : |dM - round_to_nearest_observed_mode(dM)| <= band_tol

    k=0 is excluded -- it IS the horizontal axis and would swallow every unmodified PSM.

    Returns dict with intensity-free COUNTS (a PSM is one vote; weighting by score would let a
    few high-scoring PSMs dominate and is not what "fraction of identifications" means).
    """
    out = {"n": len(psms), "rays": collections.Counter(), "band": 0, "unassigned": 0}
    for p in psms:
        mz, dm = p["mz"], p["dm"]
        if abs(dm) < 0.1:                      # unmodified: not evidence either way
            continue
        hit = None
        for k in range(-max_slope, max_slope + 1):
            if k == 0:
                continue
            if abs(dm - k * (mz - PROTON)) <= ray_tol_da:
                # prefer the smallest |k| that explains it; larger k are near-degenerate at low mz
                if hit is None or abs(k) < abs(hit):
                    hit = k
        if hit is not None:
            out["rays"][hit] += 1
        elif abs(dm - round(dm)) <= band_tol or _near_known_mod(dm, band_tol):
            out["band"] += 1
        else:
            out["unassigned"] += 1
    shifted = sum(out["rays"].values()) + out["band"] + out["unassigned"]
    out["shifted"] = shifted
    out["ray_fraction"] = (sum(out["rays"].values()) / shifted) if shifted else None

    # CHANCE FLOOR for ray_fraction. A metric that ships without a null is unsupported: a 97.8%
    # ceiling against a 91.9% floor means nothing. A PSM
    # with an arbitrary dM lands on SOME ray by coincidence: 2*max_slope lines, each 2*ray_tol_da
    # wide, inside the searched dM range. Reported alongside, never subtracted silently.
    dms = [p["dm"] for p in psms if abs(p["dm"]) >= 0.1]
    if dms:
        span = max(dms) - min(dms)
        n_lines = 2 * max_slope                       # k = +-1..+-max_slope, k=0 excluded
        out["ray_fraction_chance"] = min(1.0, n_lines * 2.0 * ray_tol_da / span) if span > 0 else None
    else:
        out["ray_fraction_chance"] = None
    return out


KNOWN_MODS = {"oxidation": 15.9949, "acetyl": 42.0106, "deamidation": 0.98402,
              "pyroglu": -17.02655, "phospho": 79.96633, "carbamidomethyl": 57.02146}


def _near_known_mod(dm, tol):
    return any(abs(dm - v) <= tol for v in KNOWN_MODS.values())


# ------------------------------------------------------------------------------------ reading
def read_psms(tsv, tag="ENTRAP_"):
    """Read a Sage TSV into the minimal record this analysis needs.

    Aborts on a missing column rather than defaulting it. A .get(col, default) here is how
    'peptide_q missing -> every row fails -> a published bold zero' happens; that exact class of
    silent-wrong-number is why this file exists.
    """
    need = ["peptide", "proteins", "expmass", "calcmass", "charge", "sage_discriminant_score"]
    rows = []
    with open(tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        missing = [c for c in need if c not in (rd.fieldnames or [])]
        if missing:
            sys.exit("ABORT: %s lacks column(s) %s. Columns present: %s"
                     % (tsv, missing, rd.fieldnames))
        for x in rd:
            try:
                exp, calc, z = float(x["expmass"]), float(x["calcmass"]), int(x["charge"])
            except (TypeError, ValueError):
                continue
            if z <= 0:
                continue
            rows.append({
                "peptide": x["peptide"],
                "entrap": tag in x["proteins"],
                "score": float(x["sage_discriminant_score"]),
                "dm": exp - calc,                       # neutral-mass delta
                "mz": (exp + z * PROTON) / z,           # observed precursor m/z
                "z": z,
            })
    return rows


def analyse(tsv, label, ratio, tag="ENTRAP_"):
    rows = read_psms(tsv, tag)
    if not rows:
        return {"label": label, "error": "no PSMs parsed"}
    thr = entrapment_fdr(rows, tag, ratio)
    if thr is None:
        # Not a zero. Open search can fail to reach 1% FDR at all, and that IS the result.
        return {"label": label, "n_psms_raw": len(rows), "entrapment_1pct_threshold": None,
                "note": "no score threshold reaches 1% entrapment FDR"}
    keep = [r for r in rows if r["score"] >= thr and not r["entrap"]]
    shifted = [r for r in keep if abs(r["dm"]) > 0.1]
    geo = ray_decomposition(keep)
    return {
        "label": label,
        "n_psms_raw": len(rows),
        "entrapment_1pct_threshold": thr,
        "psms_at_1pct": len(keep),
        "modified_psms": len(shifted),
        "modified_peptides": len(set(r["peptide"] for r in shifted)),
        "peptides": len(set(r["peptide"] for r in keep)),
        "psms_per_peptide": len(keep) / max(len(set(r["peptide"] for r in keep)), 1),
        "geometry": {"ray_fraction": geo["ray_fraction"], "band": geo["band"],
                     "unassigned": geo["unassigned"], "shifted": geo["shifted"],
                     "rays": dict(geo["rays"])},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fasta", help="build the entrapment FASTA")
    f.add_argument("--target", required=True)
    f.add_argument("--foreign", required=True)
    f.add_argument("--out", required=True)

    a = sub.add_parser("analyse", help="score one Sage TSV")
    a.add_argument("tsv")
    a.add_argument("--label", required=True)
    a.add_argument("--ratio", type=float, default=1.0,
                   help="n_entrap_proteins / n_target_proteins")
    a.add_argument("--json", help="append the result here")

    args = ap.parse_args()
    if args.cmd == "fasta":
        nt, ne = build_entrapment(args.target, args.foreign, args.out)
        print("target %d, entrapment %d, ratio %.3f -> %s" % (nt, ne, ne / nt, args.out))
        return 0

    r = analyse(args.tsv, args.label, args.ratio)
    print(json.dumps(r, indent=2))
    if args.json:
        p = Path(args.json)
        acc = json.loads(p.read_text()) if p.exists() else []
        acc.append(r)
        p.write_text(json.dumps(acc, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
