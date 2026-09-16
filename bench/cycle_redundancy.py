#!/usr/bin/env python3
"""Where does the 9.51 PSMs/peptide come from -- within a cycle, or across cycles?

DIAspeXtract emits 9.51 PSMs per peptide against the reference implementation's 1.69. That ratio alone cannot
distinguish two very different situations:

  WITHIN-CYCLE duplication   several pseudo-spectra generated from ONE diaPASEF cycle all match
                             the same peptide. That is a deconvolution/dedup failure: one cycle
                             physically contains one elution point for one precursor.

  ACROSS-CYCLE spread        the peptide elutes over many cycles and each cycle contributes one
                             spectrum. That is NORMAL and the reference implementation must show it too.

The total factorises exactly:

    PSMs/peptide  =  (PSMs per (peptide,cycle))  x  (cycles per peptide)
                     ^ within-cycle multiplicity    ^ elution spread

Only the first term is a defect. Running it on both tools separates them: if the reference implementation's
redundancy is nearly all elution spread and ours is not, the excess is ours to fix.

Cycle assignment: diaPASEF cycle period is 1.385 s for these runs (1,342 cycles over
0.02-31.00 min, 16,105 MS2 frames / 12 window groups), read from the raw analysis.tdf rather
than assumed.
"""
import argparse, collections, csv, json, sys

DECOY_PREFIX = "rev_"


def load(tsv, cycle_s, fdr=0.01):
    """PSMs passing 1% PSM-level FDR under the same common procedure as joint_bench.py."""
    rows = []
    with open(tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        need = ["peptide", "proteins", "sage_discriminant_score", "rt", "scannr"]
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            sys.exit("ABORT: %s lacks %s" % (tsv, miss))
        for x in rd:
            try:
                rt, sc = float(x["rt"]), float(x["sage_discriminant_score"])
            except (TypeError, ValueError):
                continue
            rows.append({"peptide": x["peptide"], "rt": rt, "score": sc, "scannr": x["scannr"],
                         "decoy": all(p.startswith(DECOY_PREFIX)
                                      for p in x["proteins"].split(";") if p)})
    if not rows:
        sys.exit("ABORT: no rows parsed from %s" % tsv)

    # Sage reports rt in MINUTES here; detect rather than assume, since reading seconds as
    # minutes silently puts every PSM in cycle 0 and reports perfect within-cycle duplication.
    rt_max = max(r["rt"] for r in rows)
    unit = "min" if rt_max < 200 else "s"
    to_s = 60.0 if unit == "min" else 1.0

    ranked = sorted(rows, key=lambda r: -r["score"])
    n_t = n_d = 0
    keep = []
    acc = []
    for r in ranked:
        if r["decoy"]:
            n_d += 1
        else:
            n_t += 1
            acc.append(r)
        if n_t and (n_d / n_t) <= fdr:
            keep = list(acc)
    for r in keep:
        r["cycle"] = int(r["rt"] * to_s / cycle_s)
    return keep, unit, rt_max


def analyse(psms, label):
    peps = collections.Counter(r["peptide"] for r in psms)
    pc = collections.Counter((r["peptide"], r["cycle"]) for r in psms)
    n_psm, n_pep, n_pc = len(psms), len(peps), len(pc)

    # PSMs beyond the first in each (peptide, cycle) group are within-cycle duplicates
    dup_within = sum(v - 1 for v in pc.values())
    multi = collections.Counter(pc.values())

    within = n_psm / n_pc if n_pc else 0        # PSMs per (peptide, cycle)
    spread = n_pc / n_pep if n_pep else 0       # distinct cycles per peptide

    return {
        "label": label,
        "psms": n_psm, "peptides": n_pep, "peptide_cycle_pairs": n_pc,
        "psms_per_peptide": n_psm / n_pep if n_pep else 0,
        "within_cycle_multiplicity": within,
        "cycles_per_peptide": spread,
        "pct_psms_that_are_within_cycle_duplicates": 100.0 * dup_within / n_psm if n_psm else 0,
        "multiplicity_histogram": {str(k): multi[k] for k in sorted(multi)[:8]},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs", nargs="+", metavar="LABEL=results.sage.tsv")
    ap.add_argument("--cycle-seconds", type=float, default=1.385)
    a = ap.parse_args()

    out = []
    for p in a.pairs:
        label, _, path = p.partition("=")
        psms, unit, rtmax = load(path, a.cycle_seconds)
        r = analyse(psms, label)
        r["rt_unit_detected"] = unit
        r["rt_max"] = rtmax
        out.append(r)

    print("cycle period %.3f s\n" % a.cycle_seconds)
    print("%-14s %9s %9s %9s %9s %9s %9s" %
          ("run", "PSMs", "peptides", "PSM/pep", "within", "cycles", "%dup"))
    print("-" * 76)
    for r in out:
        print("%-14s %9d %9d %9.2f %9.2f %9.2f %8.1f%%" %
              (r["label"], r["psms"], r["peptides"], r["psms_per_peptide"],
               r["within_cycle_multiplicity"], r["cycles_per_peptide"],
               r["pct_psms_that_are_within_cycle_duplicates"]))
    print("-" * 76)
    print("within  = PSMs per (peptide, cycle)   -- >1 means the SAME cycle yielded the same")
    print("          peptide more than once, which one elution point cannot justify")
    print("cycles  = distinct cycles per peptide -- normal elution spread")
    print("%dup    = share of all PSMs that are within-cycle repeats\n")
    for r in out:
        print("%-14s multiplicity histogram (PSMs per peptide-cycle): %s"
              % (r["label"], r["multiplicity_histogram"]))
    print()
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    sys.exit(main())
