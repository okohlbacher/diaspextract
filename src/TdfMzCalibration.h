#ifndef DIASPEXTRACT_TDF_MZ_CALIBRATION_H
#define DIASPEXTRACT_TDF_MZ_CALIBRATION_H
// Copyright (c) 2026, DIAspeXtract authors. BSD-3-Clause.
//
// Exact TOF -> m/z conversion for Bruker TDF (timsTOF) data, ModelType 1. Header-only and dependency-free:
// shared by the OpenMS loader patch (patches/openms-brukertims-mz-calibration.patch) and the C++ golden test
// (tests/test_calibration_cpp.cpp), so the code that ships is the code that is tested.
//
// Model, derived numerically against Bruker's timsdata library used as a local oracle (no vendor code was read,
// no vendor binary is redistributed; every constant comes from the user's own file):
//
//     t_ns   = tof * DigitizerTimebase + DigitizerDelay
//     C1_eff = C1 * (1 + dC1 * (T1_ref - T1_frame) / 1e6)        // digitizer temperature drift
//     t_ns   = C0 + (1e6 / sqrt(C1_eff)) * sqrt(m) + C2 * m      // solved for sqrt(m) =: u
//     m      = u^2 - C4                                          // C4: additive mass offset
//
// Verified to 2.5e-5 ppm against the vendor library (tests/calibration_golden.json: three files sharing ONE
// MzCalibration vector, 12 frames over 0.034 K); the C4 term was pinned the same way on a 2019 timsTOF Pro
// (tests/test_calibration_cpp.cpp). What the model does not cover (ModelType != 1, dC2 or C3 != 0, frames that
// reference more than one MzCalibration row) is REJECTED by unsupportedReason() and the loaders, never approximated.
// The open timsrust-calibration port drops the C2*m term (-11..-40 ppm on this data).

#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>

namespace diaspextract
{

/// Constants as stored in the TDF `MzCalibration` row plus the per-frame digitizer temperature.
struct TdfMzCalibration
{
  int model_type = 0;
  double digitizer_timebase = 0.0;
  double digitizer_delay = 0.0;
  double C0 = 0.0;
  double C1 = 0.0;      ///< must be > 0
  double C2 = 0.0;      ///< the quadratic term the open implementations drop
  double T1_ref = 0.0;  ///< reference digitizer temperature for the calibration
  double dC1 = 0.0;     ///< ppm/K drift of C1
  double dC2 = 0.0;     ///< drift of C2 -- NOT modelled; must be 0 (see isSupported)
  double C3 = 0.0;      ///< NOT modelled; must be 0 (never observed non-zero)
  double C4 = 0.0;      ///< additive mass offset on the quadratic root: m = u^2 - C4

  /// Why this calibration cannot be used, or empty if it can. Fail closed, never approximate.
  std::string unsupportedReason() const
  {
    // NaN defeats ordered comparisons (every < and > below is false for NaN), so screen first.
    for (double v : {digitizer_timebase, digitizer_delay, C0, C1, C2, T1_ref, dC1, dC2, C3, C4})
      if (std::isnan(v)) return "MzCalibration contains NaN";
    if (model_type != 1) return "MzCalibration ModelType " + std::to_string(model_type) + " != 1";
    // A merely-positive C1 is not enough: C1 < ~5.6e-297 overflows b*b to +inf and yields m/z 0.0
    // for every peak, silently, under an "exact" log line. Real values are ~1e5.
    if (!(C1 > 1.0) || !(C1 < 1e12)) return "MzCalibration C1 outside the plausible range (1, 1e12)";
    if (dC2 != 0.0) return "MzCalibration dC2 != 0 (temperature drift of C2 is not modelled)";
    if (C3 != 0.0) return "MzCalibration C3 != 0 (not modelled; never observed non-zero)";
    // C2 == 0.0 stored in the file is a real calibration (the 2020 timsTOF Pro firmware behind PXD017703 has no
    // quadratic term). A NULL C2, which sqlite3_column_double hands over as 0.0, is refused by the loaders
    // themselves. A negative C2 flips the root branch.
    if (C2 < 0.0) return "MzCalibration C2 < 0 (negative quadratic term is not a supported model)";
    if (!(digitizer_timebase > 0.0)) return "DigitizerTimebase <= 0";
    return std::string();
  }
  bool isSupported() const { return unsupportedReason().empty(); }

  /// Temperature-corrected 1e6/sqrt(C1_eff) for one frame; hoist out of per-peak loops.
  double frameFactor(double t1_frame) const
  {
    const double cf = 1.0 + (dC1 * (T1_ref - t1_frame)) / 1e6;
    return 1e6 / std::sqrt(C1 * cf);
  }

  /// TOF index -> m/z. @p b is frameFactor() for this frame.
  double tofToMz(double tof, double b) const
  {
    const double t = tof * digitizer_timebase + digitizer_delay;
    // Stable root of C2*u^2 + b*u + (C0 - t) = 0: u = 2*(t - C0) / (b + sqrt(disc)) avoids the textbook form's
    // cancellation (b ~ 2.5e3, C2 ~ 1e-3). Out-of-model inputs (t <= C0, a negative discriminant, m <= 0)
    // return NaN, never a plausible-looking 0.0 or the wrong mass a negative u would give.
    const double disc = b * b - 4.0 * C2 * (C0 - t);
    if (!(disc >= 0.0) || !(t > C0)) return std::numeric_limits<double>::quiet_NaN();
    const double denom = b + std::sqrt(disc);
    if (!(denom > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double u = 2.0 * (t - C0) / denom;
    const double m = u * u - C4;   // bit-identical to u*u when C4 == 0, i.e. on every earlier file
    return m > 0.0 ? m : std::numeric_limits<double>::quiet_NaN();
  }

  /// m/z -> TOF index (exact inverse of tofToMz; the model is closed-form in this direction).
  double mzToTof(double mz, double b) const
  {
    const double m = mz + C4;
    const double t = C0 + b * std::sqrt(m > 0.0 ? m : 0.0) + C2 * m;
    return (t - digitizer_delay) / digitizer_timebase;
  }
};

} // namespace diaspextract

#endif // DIASPEXTRACT_TDF_MZ_CALIBRATION_H
