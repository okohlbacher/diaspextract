#!/usr/bin/env python3
"""Which collapses cost peptides, and which were free?

Pre-extraction isotope collapse cut IH2 emission 927,813 -> 691,144 (the reference implementation: 700,434) and 15% of
runtime, but lost 1,659 of 12,537 Sage peptides. That average is useless for deciding what to keep:
the question is whether the harmful removals are a distinguishable CATEGORY.

Method, entirely on existing output, no assumption about the rule that was applied:
  * a baseline precursor coordinate absent from the collapsed output was REMOVED (coordinates are
    bit-identical -- the collapse runs after inference, on the same list)
  * every removed precursor is assigned its SURVIVOR: the nearest kept precursor on the isotope
    lattice (same charge, co-eluting, co-mobile)
  * each removal is then classed by what it cost:
      dark      the removed spectrum was never identified            -> free
      covered   its peptide is still identified somewhere else       -> free
      LOST      its peptide is gone from the collapsed result        -> harmful
Then the three classes are compared on every feature available at collapse time, which is the only
way to get a rule that could be applied without knowing the answer.

Usage: iso_loss.py <base.tsv> <base.sage.tsv> <collapsed.tsv> <collapsed.sage.tsv> [--ppm 20] [--rt 3] [--im 0.02]
"""
import sys, csv, collections, bisect

ISO = 1.0033548


def load_tsv(p):
    rows = []
    with open(p) as f:
        next(f)
        for line in f:
            q = line.rstrip("\n").split("\t")
            rows.append((float(q[1]), float(q[2]), int(q[3]), float(q[4]), int(q[5]), int(q[7] or 0)))
    return rows


def ids(p):
    out = {}
    with open(p) as f:
        for x in csv.DictReader(f, delimiter="\t"):
            if x["label"] != "1" or float(x["spectrum_q"]) > 0.01: continue
            out[int(x["scannr"].split("=")[-1])] = x["peptide"]
    return out


def main(bt, bs, ct, cs, ppm=20.0, rt_tol=3.0, im_tol=0.02):
    base, coll = load_tsv(bt), load_tsv(ct)
    bid, cid = ids(bs), ids(cs)
    bpep, cpep = set(bid.values()), set(cid.values())
    lost = bpep - cpep
    print(f"baseline {len(base):,} spectra / {len(bpep):,} peptides;  collapsed {len(coll):,} / "
          f"{len(cpep):,};  lost {len(lost):,}, gained {len(cpep-bpep):,}\n")

    # a kept precursor keeps its exact coordinate, so identity is an exact key
    kept = set()
    for rt, mz, z, im, npk, ni in coll: kept.add((round(mz, 6), round(rt, 4), z, round(im, 6)))
    surv_mz = sorted((mz, i) for i, (rt, mz, z, im, npk, ni) in enumerate(coll))
    surv_only = [m for m, _ in surv_mz]

    removed = [i for i, (rt, mz, z, im, npk, ni) in enumerate(base)
               if (round(mz, 6), round(rt, 4), z, round(im, 6)) not in kept]
    print(f"{len(removed):,} baseline precursors removed ({100.0*len(removed)/len(base):.1f}%)\n")

    def survivor(i):
        """Nearest kept precursor on the isotope lattice: returns (k, collapsed index) or None."""
        rt, mz, z, im, npk, ni = base[i]
        if z <= 0: return None
        tol = mz * ppm * 1e-6
        for k in (1, 2, 3):
            for sgn in (-1, 1):
                t = mz + sgn * k * ISO / z
                lo = bisect.bisect_left(surv_only, t - tol); hi = bisect.bisect_right(surv_only, t + tol)
                for q in range(lo, hi):
                    j = surv_mz[q][1]
                    if coll[j][2] != z: continue
                    if abs(coll[j][0] - rt) <= rt_tol and abs(coll[j][3] - im) <= im_tol:
                        return (sgn * k, j)
        return None

    cls = collections.Counter()
    feat = collections.defaultdict(list)
    for i in removed:
        pep = bid.get(i)
        c = "dark" if pep is None else ("LOST" if pep in lost else "covered")
        cls[c] += 1
        sv = survivor(i)
        rt, mz, z, im, npk, ni = base[i]
        feat[c].append({"k": sv[0] if sv else 0, "z": z, "mz": mz, "niso": ni, "npk": npk,
                        "sv_niso": coll[sv[1]][5] if sv else -1,
                        "sv_npk": coll[sv[1]][4] if sv else -1,
                        "dim": abs(coll[sv[1]][3] - im) if sv else -1,
                        "drt": abs(coll[sv[1]][0] - rt) if sv else -1,
                        "sv_id": (sv is not None and sv[1] in cid)})
    tot = sum(cls.values())
    print(f"{'class':<10}{'removals':>10}{'share':>9}   what it cost")
    for c, lbl in (("dark", "never identified -- free"),
                   ("covered", "peptide still found elsewhere -- free"),
                   ("LOST", "peptide gone -- HARMFUL")):
        print(f"{c:<10}{cls[c]:>10,}{100.0*cls[c]/tot:>8.1f}%   {lbl}")
    print(f"\n(the {cls['LOST']:,} harmful removals cost {len(lost):,} distinct peptides)\n")

    keys = [("k", "isotope step (signed; - = survivor is HEAVIER)"), ("z", "charge"),
            ("niso", "isotope peaks of the removed"), ("sv_niso", "isotope peaks of the survivor"),
            ("npk", "fragments of the removed"), ("sv_npk", "fragments of the survivor"),
            ("mz", "precursor m/z"), ("drt", "|dRT| to survivor"), ("dim", "|dIM| to survivor")]
    print(f"{'feature':<40}{'dark':>12}{'covered':>12}{'LOST':>12}   (medians)")
    for k, lbl in keys:
        vals = []
        for c in ("dark", "covered", "LOST"):
            v = sorted(x[k] for x in feat[c])
            vals.append(v[len(v) // 2] if v else float("nan"))
        print(f"{lbl:<40}{vals[0]:>12.3f}{vals[1]:>12.3f}{vals[2]:>12.3f}")
    for c in ("dark", "covered", "LOST"):
        n = len(feat[c]) or 1
        print(f"{c:<10} survivor itself identified: {100.0*sum(1 for x in feat[c] if x['sv_id'])/n:5.1f}%"
              f"   no lattice survivor found: {100.0*sum(1 for x in feat[c] if x['k']==0)/n:5.1f}%")
    print("\nsigned step of the survivor, by class (share of removals):")
    for c in ("dark", "covered", "LOST"):
        h = collections.Counter(x["k"] for x in feat[c]); n = len(feat[c]) or 1
        print(f"  {c:<9}" + "  ".join(f"{k:+d}:{100.0*h[k]/n:4.1f}%" for k in sorted(h) if k != 0))


if __name__ == "__main__":
    if len(sys.argv) < 5: sys.exit(__doc__)
    kw = {}
    for a in sys.argv[5:]:
        k, v = a.lstrip("-").split("=")
        kw[{"ppm": "ppm", "rt": "rt_tol", "im": "im_tol"}[k]] = float(v)
    main(*sys.argv[1:5], **kw)
