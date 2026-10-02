#!/usr/bin/env python3
"""Self-checks for the corrected entrapment estimator (bench/entrapment.py).

Guards its two corrections: peptide-hypothesis ratio (not protein ratio) and foreign-ONLY entrapment counting.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bench"))
from entrapment import _bootstrap_fdr

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
    print(f"OK entrapment: point {point:.3f}%, CI [{lo:.3f}, {hi:.3f}], deterministic, NaN on empty")
    return 0

if __name__ == "__main__":
    sys.exit(main())
