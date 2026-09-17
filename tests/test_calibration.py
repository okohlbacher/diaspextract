#!/usr/bin/env python3
"""Golden test for the TDF MzCalibration (ModelType 1) TOF->m/z model.

Guards the calibration that produced the 2026-09-01 +6-11% closed-search gain. The 270 golden cases
were produced by Bruker's own timsdata library and are stored in tests/calibration_golden.json, so this
test needs NO vendor library, NO cluster and NO raw data:
  * 60 on three diaPASEF files sharing one calibration vector (dataset D/dataset A/dataset B, 4 frames x 5
    TOF positions each);
  * 210 on six runs of the public PXD029836 (timsTOF Pro 2, acquisition software 2.0.53), whose C2 is
    NEGATIVE: 1418/1420/1422 share one vector, 1431, 1408 and 1347 carry one each (5 frames x 7 TOF
    positions each: first, middle and last frame, the coldest and the warmest).

Model (derived numerically against the vendor oracle). NOTE this is an INDEPENDENT re-implementation
of the model, not a copy of the shipped C++: it uses the textbook quadratic root where
src/TdfMzCalibration.h uses the cancellation-free form. That is deliberate -- two independent
implementations agreeing on the vendor numbers is a stronger check than one implementation tested
twice. tests/test_calibration_cpp.cpp is the one that exercises the code that actually ships.

    t_ns   = tof * DigitizerTimebase + DigitizerDelay
    C1_eff = C1 * (1 + dC1 * (T1_ref - T1_frame) / 1e6)          # temperature compensation
    t_ns   = C0 + (1e6 / sqrt(C1_eff)) * sqrt(m) + C2 * m        # solve the quadratic for sqrt(m)

Regressions this catches:
  * dropping the C2*m term -- the bug in the widely-copied open implementation (timsrust-calibration,
    adapted then disabled by mzdata): -11..-40 ppm on this data;
  * flipping the sign of a negative C2: 38..80 ppm on PXD029836 (checked per file). Taking the other root of the
    quadratic is not in that range: it lands beyond the parabola's vertex, near m/z 1e12;
  * dropping the dC1 temperature term: ~0.67 ppm. (m is proportional to C1_eff, NOT its square
    root, so dC1*dT/1e6 lands on m/z undiminished -- an earlier note here claimed ~0.33 ppm by
    applying a spurious halving; the measured guard value 0.669 confirms 0.67. The golden files
    span only ~0.03 K, so the guard threshold stays loose at 0.1 ppm over the whole set, and 0.05 ppm per
    file: the smallest per-file value is 0.074 ppm, PXD029836 1408, whose frames lie within 0.004 K of T1.)
Run: python3 tests/test_calibration.py
"""
import json, math, os, sys

TOL_PPM = 0.0001         # measured max 3.1e-5 ppm (this textbook root cancels on PXD029836; the shipped root: 2.6e-5);
                         # matches the C++ test's tolerance exactly
NO_C2_MIN_PPM = 5.0      # dropping C2 must be caught: it is worth >>5 ppm
NO_TEMP_MIN_PPM = 0.1    # dropping the temperature term must be visible
NO_TEMP_FILE_MIN_PPM = 0.05   # ... and visible in every file on its own

def tof_to_mz(tof, p, t1_frame, use_c2=True, use_temp=True):
    t_ns = tof * p["timebase"] + p["delay"]
    c1 = p["C1"] * (1.0 + p["dC1"] * (p["T1_ref"] - t1_frame) / 1e6) if use_temp else p["C1"]
    b = 1e6 / math.sqrt(c1)
    c2 = p["C2"] if use_c2 else 0.0
    if abs(c2) < 1e-12:                      # degenerate: pure sqrt law
        return ((t_ns - p["C0"]) / b) ** 2 - p.get("C4", 0.0)
    disc = b * b - 4.0 * c2 * (p["C0"] - t_ns)
    assert disc > 0.0, f"t = {t_ns} ns is beyond the calibration's range (negative discriminant)"
    u = (-b + math.sqrt(disc)) / (2.0 * c2)   # the rising-branch root for either sign of c2
    return u * u - p.get("C4", 0.0)       # C4: additive mass offset on the root

def ppm(a, b):
    return abs(a - b) / b * 1e6

def check_c4_offset():
    """C4 != 0 is an additive mass offset on the quadratic root (2019 timsTOF Pro, firmware 6.0.110,
    PXD017703 highsens files: C4 = -0.0905). Pinned against the vendor library on the first and the
    last frame of one file; dropping the term is -53..-938 ppm."""
    p = {"timebase": 0.2, "delay": 25224.4, "C0": 325.88517924839397, "C1": 152987.23707178448,
         "C2": 0.0030691291814964453, "C4": -0.09048112117134345, "T1_ref": 25.197341867575105, "dC1": 16.5}
    vendor = {25.201724358440018: [(2000.0, 98.00225185761244), (50000.0, 186.40848708590778), (150000.0, 461.14679357131746),
                                   (300000.0, 1102.6974507710581), (403000.0, 1702.6599405622635)],
              25.20983751539588:  [(2000.0, 98.00223875072383), (50000.0, 186.40846214421924), (150000.0, 461.14673185089845),
                                   (300000.0, 1102.6973031683212), (403000.0, 1702.6597126456059)]}
    no_c4 = dict(p, C4=0.0); worst, worst_no_c4 = 0.0, 0.0
    for t1, probes in vendor.items():
        for tof, mz in probes:
            worst = max(worst, ppm(tof_to_mz(tof, p, t1), mz))
            worst_no_c4 = max(worst_no_c4, ppm(tof_to_mz(tof, no_c4, t1), mz))
    assert worst < TOL_PPM, f"C4 case: {worst:.6f} ppm from the vendor library"
    assert worst_no_c4 > 50.0, f"C4 ablation only {worst_no_c4:.3f} ppm -- the case cannot catch a dropped C4"
    print(f"OK  C4 offset case: {worst:.6f} ppm vs vendor (10 probes, 2 frames); without C4 {worst_no_c4:.1f} ppm")

def check_linear_c2_zero():
    """C2 == 0 stored in the file (PXD017703, 2020 timsTOF Pro) is a real, purely linear-in-sqrt
    calibration and must be accepted; the quadratic collapses to sqrt(m) = (t - C0) / b."""
    p = {"timebase": 0.2, "delay": 25131.0, "C0": 315.70325869866065,
         "C1": 154272.1271422364, "C2": 0.0, "T1_ref": 25.63315397685876, "dC1": -0.2}
    import math
    for tof in (60000.0, 120000.0, 240000.0, 400000.0):
        got = tof_to_mz(tof, p, p["T1_ref"])
        t = tof * p["timebase"] + p["delay"]
        b = 1e6 / math.sqrt(p["C1"])
        expect = ((t - p["C0"]) / b) ** 2
        assert abs(got - expect) / expect < 1e-12, f"C2==0 closed form: tof {tof} got {got} want {expect}"
    print("OK  C2 == 0 (linear law, PXD017703 row) accepted and matches the closed form")

def main():
    check_linear_c2_zero()
    check_c4_offset()
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration_golden.json")
    files = json.load(open(path))
    n = 0
    worst = 0.0
    worst_noc2 = 0.0
    worst_notemp = 0.0
    negative = 0
    for f in files:
        assert f["model_type"] == 1, f"{f['file']}: golden set must be ModelType 1"
        f_worst = f_noc2 = f_flip = f_notemp = 0.0
        flipped = dict(f, C2=-f["C2"])
        for case in f["cases"]:
            got = tof_to_mz(case["tof"], f, case["t1"])
            e = ppm(got, case["mz"])
            assert e < TOL_PPM, (
                f"{f['file']} frame {case['frame']} tof {case['tof']:.0f}: "
                f"model {got:.6f} vs vendor {case['mz']:.6f} = {e:.5f} ppm > {TOL_PPM}")
            f_worst = max(f_worst, e)
            f_noc2 = max(f_noc2, ppm(tof_to_mz(case["tof"], f, case["t1"], use_c2=False), case["mz"]))
            f_flip = max(f_flip, ppm(tof_to_mz(case["tof"], flipped, case["t1"]), case["mz"]))
            f_notemp = max(f_notemp, ppm(tof_to_mz(case["tof"], f, case["t1"], use_temp=False), case["mz"]))
            n += 1
        # per file: every file must tell the model from the known-bad variants on its own
        assert (len(f["cases"]) >= 20 and f_noc2 > NO_C2_MIN_PPM and f_flip > NO_C2_MIN_PPM
                and f_notemp > NO_TEMP_FILE_MIN_PPM), (
            f"{f['file']}: {len(f['cases'])} cases, without C2 {f_noc2:.3f} ppm, C2 sign-flipped {f_flip:.3f} ppm, "
            f"without temperature {f_notemp:.4f} ppm")
        print(f"    {f['file']:<16} C2 {f['C2']:+.6e}  {len(f['cases']):2d} cases  max {f_worst:.7f} ppm  "
              f"(without C2 {f_noc2:.2f}, sign-flipped {f_flip:.2f}, without temperature {f_notemp:.3f} ppm)")
        worst, worst_noc2, worst_notemp = max(worst, f_worst), max(worst_noc2, f_noc2), max(worst_notemp, f_notemp)
        negative += f["C2"] < 0
    assert len(files) >= 9 and n >= 270 and negative >= 6, f"only {n} cases in {len(files)} files ({negative} with C2 < 0)"
    # the ablations must actually be detectable, or this test proves nothing
    assert worst_noc2 > NO_C2_MIN_PPM, (
        f"dropping C2 changed the result by only {worst_noc2:.3f} ppm -- the golden set cannot "
        f"distinguish the known-bad open implementation; regenerate it over a wider m/z range")
    assert worst_notemp > NO_TEMP_MIN_PPM, (
        f"dropping the temperature term changed the result by only {worst_notemp:.4f} ppm -- "
        f"the golden set has no temperature spread; regenerate across more frames")
    print(f"OK  {n} golden cases, {len(files)} files ({negative} with C2 < 0): max error {worst:.7f} ppm (tol {TOL_PPM})")
    print(f"    ablation guards: without C2 {worst_noc2:.2f} ppm, without temperature {worst_notemp:.4f} ppm")
    return 0

if __name__ == "__main__":
    sys.exit(main())
