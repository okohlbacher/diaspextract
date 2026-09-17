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
// (tests/test_calibration_cpp.cpp). A NEGATIVE C2 (timsTOF Pro 2, acquisition software 2.0.53: all 23 runs checked
// of the public PXD029836, four calibration vectors) is the same law and the same root: pinned on six of those runs
// (30 frames, 7.0e-6 ppm worst; dropping C2 is 19..40 ppm off, flipping its sign 38..80 ppm).
// What the model does not cover (ModelType != 1, dC2 or C3 != 0, frames that reference more than one MzCalibration
// row, a negative C2 outside its verified domain) is REJECTED by unsupportedReason() and the loaders, never
// approximated. The open timsrust-calibration port drops the C2*m term (-11..-40 ppm on this data).
//
// The negative-C2 domain. With C2 < 0 the law t(u) = C0 + b*u + C2*u^2 (u = sqrt(m + C4)) is a parabola that peaks at
// u = b / (2|C2|): no mass has a later flight time (the discriminant below is negative there and tofToMz returns NaN),
// and near the peak the root is ill-conditioned. At the root, sqrt(disc) is the slope dt/du, so unsupportedReason()
// requires disc >= b^2 / 4 -- the slope at least half its linear value -- over everything the tool can ask of the
// model: every TOF index up to kDomainMaxTof, every m/z up to kDomainMaxMz, and every frame temperature within
// kDomainT1Span of T1_ref. disc falls with the flight time and with the C1 factor, so the corner that decides is the
// largest TOF at the largest factor. On the six PXD029836 runs that corner has disc >= 0.984 b^2, the runs' own range
// (397,657 bins, frames within 0.011 K of T1_ref) >= 0.99992 b^2, and the root would only stop existing beyond TOF
// index 6.4e9, 16,000 times the digitizer's range. frameReason() lets a loader refuse a frame outside the checked
// temperature span: the tool's tdf reader (TdfLoad.h) does, the OpenMS loader patch predates it and does not -- a frame
// there converts as a positive-C2 frame at such a temperature would, and a peak beyond the parabola's peak comes back
// NaN, which the tool drops and counts (bad_mz). Positive and zero C2 take none of these branches.

#include <cmath>
#include <cstdio>
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
  double C2 = 0.0;      ///< the quadratic term the open implementations drop; negative only inside its domain
  double T1_ref = 0.0;  ///< reference digitizer temperature for the calibration
  double dC1 = 0.0;     ///< ppm/K drift of C1
  double dC2 = 0.0;     ///< drift of C2 -- NOT modelled; must be 0 (see isSupported)
  double C3 = 0.0;      ///< NOT modelled; must be 0 (never observed non-zero)
  double C4 = 0.0;      ///< additive mass offset on the quadratic root: m = u^2 - C4

  /// The negative-C2 domain (header comment): what the tool can hand the model. TdfLoad.h refuses a file whose
  /// DigitizerNumSamples (every stored TOF index is below it) or MzAcqRangeUpper exceeds these, and the compact
  /// store caps a peak's m/z at 42,000.
  static constexpr double kDomainMaxTof = 1e8;
  static constexpr double kDomainMaxMz = 1e5;
  /// Kelvin around T1_ref: wider than any temperature a working digitizer can have (T1_ref is ~25 C, so the span
  /// reaches below absolute zero); the golden files' frames span 0.004 to 0.034 K.
  static constexpr double kDomainT1Span = 1000.0;

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
    if (!(digitizer_timebase > 0.0)) return "DigitizerTimebase <= 0";
    // C2 == 0.0 stored in the file is a real calibration (the 2020 timsTOF Pro firmware behind PXD017703 has no
    // quadratic term). A NULL C2, which sqlite3_column_double hands over as 0.0, is refused by the loaders
    // themselves. A negative C2 is accepted inside its domain only.
    if (C2 < 0.0) return negativeC2Reason();
    return std::string();
  }
  bool isSupported() const { return unsupportedReason().empty(); }

  /// Why a frame at digitizer temperature @p t1_frame cannot be converted, or empty. Only a negative C2 has a
  /// temperature-dependent domain, and unsupportedReason() checked it for T1_ref +- kDomainT1Span: a frame outside
  /// that span (a corrupt T1) is refused rather than converted where nothing was checked.
  std::string frameReason(double t1_frame) const
  {
    if (!(C2 < 0.0) || std::fabs(t1_frame - T1_ref) <= kDomainT1Span) return std::string();
    return "Frames.T1 " + num(t1_frame) + " is not within " + num(kDomainT1Span) + " K of the MzCalibration T1 " + num(T1_ref)
           + ", the temperature span over which the negative C2 was checked";
  }

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
    // cancellation (b ~ 2.5e3, |C2| ~ 1e-3). It is the root on the rising branch for either sign of C2 and
    // continuous through C2 = 0, where it is (t - C0) / b. Out-of-model inputs (t <= C0, a negative discriminant --
    // with C2 < 0, a flight time beyond the parabola's peak -- m <= 0) return NaN, never a plausible-looking 0.0
    // or the wrong mass a negative u would give.
    const double disc = b * b - 4.0 * C2 * (C0 - t);
    if (!(disc >= 0.0) || !(t > C0)) return std::numeric_limits<double>::quiet_NaN();
    const double denom = b + std::sqrt(disc);
    if (!(denom > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double u = 2.0 * (t - C0) / denom;
    const double m = u * u - C4;   // bit-identical to u*u when C4 == 0, i.e. on every earlier file
    return m > 0.0 ? m : std::numeric_limits<double>::quiet_NaN();
  }

  /// m/z -> TOF index (exact inverse of tofToMz; the model is closed-form in this direction). With C2 < 0 it is the
  /// inverse on the rising branch, which unsupportedReason() verified up to kDomainMaxMz.
  double mzToTof(double mz, double b) const
  {
    const double m = mz + C4;
    const double t = C0 + b * std::sqrt(m > 0.0 ? m : 0.0) + C2 * m;
    return (t - digitizer_delay) / digitizer_timebase;
  }

private:
  static std::string num(double v)
  {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
  }

  /// The negative-C2 domain check (header comment). disc = b^2 - k*(t - C0) with k = 4|C2| is smallest at the
  /// largest flight time and the smallest frame factor b, i.e. the largest C1 factor within the temperature span.
  std::string negativeC2Reason() const
  {
    const double cf_lo = 1.0 - std::fabs(dC1) * kDomainT1Span / 1e6, cf_hi = 1.0 + std::fabs(dC1) * kDomainT1Span / 1e6;
    if (!(cf_lo > 0.0))
      return "MzCalibration C2 < 0 with dC1 " + num(dC1) + " ppm/K: the C1 factor leaves (0, inf) within " + num(kDomainT1Span)
             + " K of T1, so the negative C2 cannot be checked over that span";
    const double b = 1e6 / std::sqrt(C1 * cf_hi);
    const double k = -4.0 * C2;
    const double t_zero = C0 + b * b / k;          // disc = 0: beyond it no mass has this flight time
    const double t_half = C0 + 0.75 * b * b / k;   // disc = b^2 / 4: the slope dt/du is half its linear value
    const double u_half = b / k;                   // the root there, sqrt(m + C4)
    const double t_max = kDomainMaxTof * digitizer_timebase + digitizer_delay;
    if (t_max <= t_half && kDomainMaxMz + C4 <= u_half * u_half) return std::string();
    return "MzCalibration C2 " + num(C2) + " < 0 outside its verified domain: the TOF-to-m/z root stops existing at TOF index "
           + num((t_zero - digitizer_delay) / digitizer_timebase) + " and its slope halves at TOF index "
           + num((t_half - digitizer_delay) / digitizer_timebase) + " (m/z " + num(u_half * u_half - C4)
           + "); a negative C2 is accepted only while the slope stays above half up to TOF index " + num(kDomainMaxTof)
           + " and m/z " + num(kDomainMaxMz) + ", at every frame temperature within " + num(kDomainT1Span) + " K of T1";
  }
};

} // namespace diaspextract

#endif // DIASPEXTRACT_TDF_MZ_CALIBRATION_H
