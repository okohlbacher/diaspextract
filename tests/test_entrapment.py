#!/usr/bin/env python3
"""Self-checks for the corrected entrapment estimator (bench/entrapment.py).

Guards its two corrections: peptide-hypothesis ratio (not protein ratio) and foreign-ONLY entrapment counting.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bench"))
from entrapment import _bootstrap_fdr, peptide_ratio, digest_rule, _digest
import entrapment, random
import tempfile

def main():
    ratio = 0.6805
    lo, hi = _bootstrap_fdr([0] * 990 + [1] * 10, ratio)
    assert 0.0 < lo < hi < 10.0, (lo, hi)
    # corrected FDR = (n_entrap/ratio)/n_target -> ~1.5% here; the CI must bracket it
    point = 100.0 * (10 / ratio) / 990
    assert lo <= point <= hi, f"CI [{lo:.3f},{hi:.3f}] must bracket point estimate {point:.3f}"
    lo0, hi0 = _bootstrap_fdr([], ratio)
    assert lo0 != lo0 and hi0 != hi0, "empty input must yield NaN CI, not a number"
    lo_d, hi_d = _bootstrap_fdr([0] * 990 + [1] * 10, ratio)
    assert (lo, hi) == (lo_d, hi_d), "bootstrap must be deterministic (fixed seed)"
    # the peptide ratio digests as the search does: after E, a target without E is one peptide, a foreign one two
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(">t\nAAAAAAKAAAAAAK\n>ENTRAP_f\nAAAAAAEAAAAAAE\n")
    gluc = digest_rule({"database": {"enzyme": {"cleave_at": "E", "restrict": ""}}})
    r_try, r_glu = peptide_ratio(f.name), peptide_ratio(f.name, rule=gluc)
    os.unlink(f.name)
    assert (r_try, r_glu) == (0.5, 2.0), f"peptide ratio trypsin {r_try} (want 0.5), cleave after E {r_glu} (want 2.0)"
    # Sage's other digestion modes (crates/sage/src/enzyme.rs), on peptides of length 2..4 without a mass bound
    small = {"min_len": 2, "max_len": 4, "missed_cleavages": 0}
    def rule(**e):
        return digest_rule({"database": {"enzyme": dict(small, **e), "peptide_min_mass": 0, "peptide_max_mass": 1e9}})
    assert _digest("AAKAA", **rule()) == {"AAK", "AA"}, "trypsin, cut after K"
    assert _digest("AAKAA", **rule(c_terminal=False)) == {"AA", "KAA"}, "c_terminal false: cut before K"
    assert _digest("AAKPAA", **rule(c_terminal=False)) == {"AA", "KPAA"} and _digest("AAKPAA", **rule()) == set(), \
        "restrict P: blocks the cut after K; before K (Sage checks the residue at the cut) it does not"
    assert _digest("AAKAA", **rule(cleave_at="")) == {"AA", "AK", "KA", "AAK", "AKA", "KAA", "AAKA", "AKAA"}, "non-specific"
    assert _digest("AAKAA", **rule(semi_enzymatic=True)) == {"AAK", "AA", "AK"}, "semi-enzymatic"
    assert _digest("AAKAA", **rule(cleave_at="$", max_len=5)) == {"AAKAA"}, "no cleavage"
    # non-specific: every window of length 2..4 is a hypothesis; digested as whole proteins the ratio would be 1.0
    with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as f:
        f.write(">t\nAAA\n>ENTRAP_f\nCCCC\n")
    r_ns = peptide_ratio(f.name, rule=rule(cleave_at=""))
    assert r_ns == 3.0 / 2.0, f"non-specific ratio {r_ns}: CC, CCC, CCCC vs AA, AAA"
    old_max, entrapment.MAX_NONSPECIFIC = entrapment.MAX_NONSPECIFIC, 1
    try:
        peptide_ratio(f.name, rule=rule(cleave_at=""))
        raise AssertionError("a non-specific search beyond MAX_NONSPECIFIC must be refused")
    except SystemExit as x:
        assert "non-specific" in str(x)
    finally:
        entrapment.MAX_NONSPECIFIC = old_max
        os.unlink(f.name)
    # trypsin as before the rewrite: the cut list of 1.7.0's _digest on random sequences
    def old_digest(seq, mc=2, lo=7, hi=35, mlo=500.0, mhi=5000.0, cleave="KR", restrict="P"):
        cuts = [0] + [i+1 for i in range(len(seq)-1) if seq[i] in cleave and seq[i+1] not in (restrict or "")] + [len(seq)]
        peps, n = set(), len(cuts)-1
        for i in range(n):
            for j in range(i+1, min(i+2+mc, n+1)):
                p = seq[cuts[i]:cuts[j]]
                if lo <= len(p) <= hi and all(a in entrapment._AA for a in p):
                    m = 18.010565 + sum(entrapment._AA[a] for a in p)
                    if mlo <= m <= mhi: peps.add(p)
        return peps
    rng = random.Random(7)
    for _ in range(300):
        seq = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWYKRKRP") for _ in range(rng.randint(1, 120)))
        assert _digest(seq) == old_digest(seq), seq
    print(f"OK entrapment: point {point:.3f}%, CI [{lo:.3f}, {hi:.3f}], deterministic, NaN on empty; ratio follows the enzyme")
    return 0

if __name__ == "__main__":
    sys.exit(main())
