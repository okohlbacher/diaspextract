// Copyright (c) 2026, DIAspeXtract authors. BSD-3-Clause.
//
// PORTED CODE: this file is a C++ port of the default MS1 path of dnoise v0.1.0
//   Garrett, Diedrich & Yates III, bioRxiv 2026.08.27.747603; github.com/pgarrett-scripts/dnoise,
//   tag v0.1.0 (commit fa750c9). MIT License, Copyright (c) 2026 Patrick Garrett.
// The MIT notice is reproduced in LICENSES/dnoise-MIT.txt and covers the ported logic: src/filter.rs,
// src/halo.rs, src/dia_ms1.rs, the gate construction in src/writer.rs and src/tdf/mod.rs, and the
// timsrust 0.4.2 metadata converters dnoise builds that gate with.
//
// Header-only and dependency-free for the same reason as TdfMzCalibration.h: the tool and the
// OpenMS-free unit test (tests/test_dnoise_ms1.cpp) compile the same code.
//
// The contract is BIT-IDENTITY with dnoise, not a better filter: for one MS1 frame, the keep mask is
// exactly the set of points dnoise v0.1.0 (default flags, diaPASEF) writes back -- for GlobalMetadata values
// that are finite decimal literals equal to zero or of magnitude in [DBL_MIN, DBL_MAX]. A value dnoise would
// parse otherwise (inf/infinity/nan in any case, overflow to inf, underflow to 0 or a subnormal) makes the port
// refuse the run where dnoise builds a degenerate gate. Every rule below was
// checked against the binary's output; each "obvious fix" changes that output measurably (points that
// differ on D frame 8724, 144,814 kept):
//   1. Stage order: the streak filter on the raw points, `iterations` passes, each summing only the
//      previous pass's survivors (1 pass: 24,906); then the halo ONCE, with the streak survivors as
//      its reference set (raw points as reference: 122); then the window gate ANDed (gate first: 199).
//   2. Raw u32 intensities from the frame blob, not the corrected ones the loaders return -- opentims computes
//      (uint32)(double(raw) * (100.0 / AccumulationTime) + 0.5), and the Bruker SDK corrects them too. At default
//      settings they only enter the halo ratio, yet on 24 TNBC 009 frames the corrected values flip 9 of 3.1
//      million kept points.
//   3. Streak: min_feature_length counts OCCUPIED scans of a gap-bridged run, not its span; points are
//      kept by span, which only matters for zero-intensity points.
//   4. Halo: the reference skips EVERY point in the point's own TOF column (that is what spares a
//      streak), ties are kept (strict >: 15), and an empty box keeps the point.
//   5. Gate: timsrust's two-point chord calibration (sqrt-linear over MzAcqRange/DigitizerNumSamples,
//      linear 1/K0 over max NumScans), NOT TdfMzCalibration.h; +-5 Da on the m/z range when the file
//      says "Bruker otofControl" (TOF edges move 317-1,953 bins without it); ScanNumEnd -- exclusive in
//      the schema -- converted as if inclusive (42); the union of every window row of every group.
//   6. Floating-point contraction. On D the 0.05 1/K0 pad is 59 scans to within an ulp, so gate scan
//      bounds are floor/ceil of values like 542.9999999999999; fusing `intercept + slope * scan` into
//      one FMA flips 3 of D's 24 boxes (13 points on 8724). The builder keeps that product apart
//      through a volatile and is compiled with contraction off; CMakeLists.txt passes -ffp-contract=off
//      as well. Value-changing flags the compiler announces in a macro are refused by the #error below.
//      What that catches depends on the compiler, because only an announced flag can be seen (measured in
//      CI): GCC 16 announces all four of -freciprocal-math, -funsafe-math-optimizations, -fno-signed-zeros
//      and -fassociative-math; GCC 11 announces none of them; clang announces only -funsafe-math-optimizations
//      and -ffp-model=fast, both on x86-64. Do not use them: -freciprocal-math,
//      -funsafe-math-optimizations and -ffp-model=fast each move D's boxes. tests/test_dnoise_ms1.cpp pins
//      D's boxes, the three FMA-sensitive ones included.
// Not ported, because dnoise's defaults never reach it on diaPASEF MS1: the mz_ppm width override,
// frame_half_width neighbourhoods, the selection-polygon gate (skipped whenever DiaFrameMsMsInfo and
// DiaFrameMsMsWindows map frames to windows), crop, smooth/watershed/box-centroid, and the MS2 paths.
// Also here, because the in-tool bit-identity rests on them: recover::, which inverts what the Bruker loader
// delivers (m/z, float32 1/K0, corrected intensity) back to the integers dnoise filters, and FrameCoverage,
// which holds the tool to denoising every MS1 frame exactly once.

#pragma once
#include <algorithm>
#include <array>
#include <cfloat>
#include <charconv>   // integer std::from_chars (parseU32)
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <locale>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Point 6: strict IEEE binary64 or no build, for every value-changing flag the compiler announces in a macro (GCC one per
// flag, clang on arm64 only __ARM_FP_FAST); none can be switched off per function. __NO_TRAPPING_MATH__ is left out: it is
// clang's default behaviour, under which every golden was measured.
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || defined(_M_FP_FAST) || \
    (defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD != 0) || defined(__RECIPROCAL_MATH__) || defined(__ASSOCIATIVE_MATH__) || \
    defined(__NO_SIGNED_ZEROS__) || defined(__ARM_FP_FAST)
#error "DnoiseMs1.h needs strict IEEE binary64 arithmetic to stay bit-identical with dnoise: build without -ffast-math/-Ofast, -ffinite-math-only, -funsafe-math-optimizations, -freciprocal-math, -fassociative-math, -fno-signed-zeros, -ffp-model=fast, /fp:fast and x87 excess precision (-mfpmath=387)"
#endif

namespace spx::dnoise
{

inline constexpr const char* kPortedFrom = "dnoise v0.1.0 (fa750c9)";

/// dnoise v0.1.0's MS1 knobs at its built-in defaults (params.rs:23-34, 52-60, 189-196; main.rs:367-437, 598-607).
/// Integer widths follow the Rust fields; unsigned arithmetic saturates or wraps where Rust's does.
struct Params
{
  // streak ("vertical IM feature") filter, filter.rs
  std::uint32_t mz_half_width = 3;          ///< window [c - w, c + w] around each TOF index c, saturating
  std::size_t min_feature_length = 5;       ///< occupied scans a run needs; bridged empty scans do not count
  std::size_t max_internal_gap = 2;         ///< empty scans a run may bridge (a larger step breaks it)
  std::uint64_t min_window_intensity = 0;   ///< window sum a scan needs to count as occupied (and > 0)
  std::uint64_t min_feature_intensity = 0;  ///< window sum over the run's whole span
  std::size_t iterations = 2;               ///< streak passes, each over the previous survivors; 0 = no streak pass (halo and gate still apply)
  // halo filter, halo.rs
  bool halo = true;
  double halo_peak_fraction = 0.15;         ///< keep iff I >= fraction * (max off-column I in the box)
  std::uint32_t halo_mz_idx_half_width = 80;
  std::size_t halo_scan_half_width = 2;
  // diaPASEF MS1 window gate, dia_ms1.rs + writer.rs:978-1021
  bool dia_ms1_window = true;
  double dia_ms1_mz_pad = 5.0;              ///< Da added on both sides of every isolation window
  double dia_ms1_im_pad = 0.05;             ///< 1/K0 added on both sides
};

/// What timsrust 0.4.2's MetadataReader reads for the gate: GlobalMetadata values exactly as stored (TEXT,
/// unparsed -- the builder parses them the way Rust does; std::nullopt for an absent key), and
/// SELECT MAX(NumScans) FROM Frames over ALL frames, MS1 and MS2 alike.
struct RunMeta
{
  std::optional<std::string> acquisition_software;         ///< AcquisitionSoftware
  std::optional<std::string> mz_acq_range_lower;           ///< MzAcqRangeLower
  std::optional<std::string> mz_acq_range_upper;           ///< MzAcqRangeUpper
  std::optional<std::string> digitizer_num_samples;        ///< DigitizerNumSamples
  std::optional<std::string> one_over_k0_acq_range_lower;  ///< OneOverK0AcqRangeLower
  std::optional<std::string> one_over_k0_acq_range_upper;  ///< OneOverK0AcqRangeUpper
  std::uint32_t max_num_scans = 0;
};

/// One row of tdf/mod.rs:246-249, i.e. of
///   SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth FROM DiaFrameMsMsWindows
///   WHERE IsolationMz IS NOT NULL AND IsolationWidth IS NOT NULL
/// No such table means no rows. DISTINCT does not change the union, so duplicates are harmless.
struct WindowRow
{
  std::int64_t scan_num_begin = 0;
  std::int64_t scan_num_end = 0;   ///< exclusive in the schema; dnoise converts it as if inclusive
  double isolation_mz = 0.0;
  double isolation_width = 0.0;
};

/// A padded isolation window as an inclusive integer (scan, TOF index) box (dia_ms1.rs TofScanBox).
struct GateBox
{
  std::uint32_t scan_lo = 0, scan_hi = 0, tof_lo = 0, tof_hi = 0;
};

/// The union of all padded windows as per-scan merged TOF intervals (dia_ms1.rs DiaMs1Gate). Immutable once
/// built, so one gate serves every thread.
class Ms1WindowGate
{
public:
  /// dia_ms1.rs:43-73. No gate without boxes or scans. A box with tof_hi < tof_lo is skipped; scan bounds are
  /// clamped onto the last scan; touching intervals (lo <= hi + 1) merge, which leaves membership unchanged.
  static std::optional<Ms1WindowGate> fromBoxes(const std::vector<GateBox>& boxes, std::uint32_t num_scans)
  {
    if (boxes.empty() || num_scans == 0) return std::nullopt;
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> per_scan(num_scans);
    for (const GateBox& b : boxes)
    {
      if (b.tof_hi < b.tof_lo) continue;
      const std::uint32_t lo = std::min(b.scan_lo, num_scans - 1), hi = std::min(b.scan_hi, num_scans - 1);
      for (std::uint32_t s = lo; s <= hi; ++s) per_scan[s].emplace_back(b.tof_lo, b.tof_hi);
    }
    Ms1WindowGate g;
    g.boxes_ = boxes;
    g.row_.assign(1, 0);
    for (auto& row : per_scan)
    {
      std::sort(row.begin(), row.end());
      for (const auto& iv : row)
      {
        const bool merge = g.iv_.size() > g.row_.back() &&
                           iv.first <= (g.iv_.back().second == std::numeric_limits<std::uint32_t>::max() ? g.iv_.back().second : g.iv_.back().second + 1);
        if (merge) g.iv_.back().second = std::max(g.iv_.back().second, iv.second);
        else g.iv_.push_back(iv);
      }
      g.row_.push_back(g.iv_.size());
    }
    return g;
  }

  /// dia_ms1.rs:76-84: inside some window; a scan past the run's last scan is outside.
  bool contains(std::uint32_t scan, std::uint32_t tof) const
  {
    if (scan >= row_.size() - 1) return false;
    const auto b = iv_.begin() + static_cast<std::ptrdiff_t>(row_[scan]);
    const auto e = iv_.begin() + static_cast<std::ptrdiff_t>(row_[scan + 1]);
    const auto it = std::partition_point(b, e, [tof](const std::pair<std::uint32_t, std::uint32_t>& iv) { return iv.second < tof; });
    return it != e && tof >= it->first;
  }

  /// The padded boxes in window-row order, before the per-scan clamp and merge (for logging and tests).
  const std::vector<GateBox>& boxes() const { return boxes_; }

private:
  Ms1WindowGate() = default;
  std::vector<GateBox> boxes_;
  std::vector<std::size_t> row_;                                ///< intervals of scan s: [row_[s], row_[s + 1])
  std::vector<std::pair<std::uint32_t, std::uint32_t>> iv_;     ///< inclusive [tof_lo, tof_hi], sorted per scan
};

namespace detail
{
/// Rust's `f as u32`: saturating, NaN -> 0.
inline std::uint32_t asU32(double v)
{
  if (!(v > 0.0)) return 0;
  if (v >= 4294967296.0) return std::numeric_limits<std::uint32_t>::max();
  return static_cast<std::uint32_t>(v);
}

/// `str::parse::<u32>()`: one optional '+', then decimal digits, no overflow, nothing else.
inline bool parseU32(const std::string& s, std::uint32_t& out)
{
  const char* b = s.data();
  const char* const e = b + s.size();
  if (b != e && *b == '+') ++b;
  if (b == e || *b < '0' || *b > '9') return false;
  const auto r = std::from_chars(b, e, out);
  return r.ec == std::errc() && r.ptr == e;
}

/// `str::parse::<f64>()` on the finite-decimal domain of the contract at the top: Rust's grammar
/// [+-]? (d+ | d+.d* | d*.d+) ([eE][+-]?d+)?, no whitespace, correctly rounded, locale-independent (a Qt application
/// may have set a decimal-comma C locale), and a value that is zero or of magnitude in [DBL_MIN, DBL_MAX]. inf, nan,
/// overflow, underflow to zero and subnormals -- all of which Rust accepts -- are refused, alike on every library.
inline bool parseF64(const std::string& s, double& out)
{
  std::size_t i = 0;
  if (i < s.size() && (s[i] == '+' || s[i] == '-')) ++i;
  const std::size_t start = i;
  bool nonzero = false;   // a non-zero MANTISSA digit: "0e999" is zero, "1e-999" underflows
  auto digits = [&s, &i, &nonzero](bool mantissa) {
    std::size_t k = 0;
    for (; i < s.size() && s[i] >= '0' && s[i] <= '9'; ++i, ++k) nonzero = nonzero || (mantissa && s[i] != '0');
    return k; };
  std::size_t mantissa = digits(true);
  if (i < s.size() && s[i] == '.') { ++i; mantissa += digits(true); }
  if (mantissa == 0) return false;
  if (i < s.size() && (s[i] == 'e' || s[i] == 'E'))
  {
    ++i;
    if (i < s.size() && (s[i] == '+' || s[i] == '-')) ++i;
    if (digits(false) == 0) return false;
  }
  if (i != s.size()) return false;
  const std::string body = (s[0] == '-' ? "-" : "") + s.substr(start);
  // ONE conversion path on every platform: the classic-locale stream, whose num_get converts through strtod in the C
  // locale (correctly rounded; the grammar is already checked, so it only converts). Libraries disagree on how they
  // report a range error -- failbit, or a silent 0 or subnormal -- so the domain is checked on the value instead.
  std::istringstream in(body);
  in.imbue(std::locale::classic());
  double v = 0.0;
  in >> v;
  const double a = std::fabs(v);
  if (in.fail() || (nonzero ? !(a >= DBL_MIN && a <= DBL_MAX) : v != 0.0)) return false;
  out = v;
  return true;
}
} // namespace detail

/// buildMs1WindowGate's result.
struct GateBuild
{
  std::optional<Ms1WindowGate> gate;  ///< no value: this run has no gate (not requested, no usable window rows, no scans)
  std::string error;                  ///< non-empty: the run is refused here -- a timsrust MetadataReader error, or a metadata value outside the port's finite-decimal domain
};

// Contraction off for everything the gate builder and the point recovery compute (point 6 at the top). Clang attaches its
// pragma to each operation, so it survives inlining, but an explicit -ffp-contract=fast overrides it -- hence also the
// volatile in the builder. GCC applies the optimize pragma to every function defined before pop_options.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("fp-contract=off")
#endif

/// The run's diaPASEF MS1 window gate, exactly as dnoise v0.1.0 builds it once per run (writer.rs:978-1021 with
/// timsrust 0.4.2 metadata_reader.rs:59-102, tof_to_mz.rs, scan_to_im.rs), in the Rust evaluation order.
inline GateBuild buildMs1WindowGate(const RunMeta& meta, const std::vector<WindowRow>& rows, const Params& p)
{
#if defined(__clang__)
  #pragma clang fp contract(off)
#endif
  GateBuild out;
  if (!p.dia_ms1_window) return out;
  // tdf/mod.rs:251-252 cast with `i64 as u32` (wrapping); :265 drops ScanNumEnd <= ScanNumBegin.
  std::vector<WindowRow> usable;
  for (const WindowRow& r : rows)
    if (static_cast<std::uint32_t>(r.scan_num_end) > static_cast<std::uint32_t>(r.scan_num_begin)) usable.push_back(r);
  if (usable.empty()) return out;   // writer.rs:984-986: no gate, and the metadata is never read

  auto number = [&out](const std::optional<std::string>& v, const char* key, double& x) {
    if (!v) out.error = std::string("GlobalMetadata has no ") + key;
    else if (!detail::parseF64(*v, x)) out.error = std::string("GlobalMetadata ") + key + " '" + *v + "' is not a finite decimal number of magnitude zero or in [DBL_MIN, DBL_MAX]";
    return out.error.empty();
  };
  double mz_min = 0.0, mz_max = 0.0, im_min = 0.0, im_max = 0.0;
  std::uint32_t tof_max_index = 0;
  if (!meta.acquisition_software) { out.error = "GlobalMetadata has no AcquisitionSoftware"; return out; }
  if (!number(meta.mz_acq_range_lower, "MzAcqRangeLower", mz_min) || !number(meta.mz_acq_range_upper, "MzAcqRangeUpper", mz_max) ||
      !number(meta.one_over_k0_acq_range_lower, "OneOverK0AcqRangeLower", im_min) ||
      !number(meta.one_over_k0_acq_range_upper, "OneOverK0AcqRangeUpper", im_max))
    return out;
  if (!meta.digitizer_num_samples || !detail::parseU32(*meta.digitizer_num_samples, tof_max_index))
  {
    out.error = "GlobalMetadata DigitizerNumSamples '" + meta.digitizer_num_samples.value_or("") + "' is missing or not a u32";
    return out;
  }
  if (*meta.acquisition_software == "Bruker otofControl") { mz_min -= 5.0; mz_max += 5.0; }   // metadata_reader.rs:67-70
  if (meta.max_num_scans == 0) return out;                                                     // writer.rs:988-991
  const std::uint32_t num_scans = meta.max_num_scans;

  const double tof_intercept = std::sqrt(mz_min);                                              // Tof2MzConverter::from_boundaries
  const double tof_slope = (std::sqrt(mz_max) - tof_intercept) / static_cast<double>(tof_max_index);
  const double scan_intercept = im_max;                                                        // Scan2ImConverter::from_boundaries
  const double scan_slope = (im_min - scan_intercept) / static_cast<double>(num_scans);
  std::vector<GateBox> boxes;
  for (const WindowRow& r : usable)
  {
    const double mz_lo = r.isolation_mz - r.isolation_width / 2.0;                            // tdf/mod.rs:258-259
    const double mz_hi = r.isolation_mz + r.isolation_width / 2.0;
    const double t0 = (std::sqrt(mz_lo - p.dia_ms1_mz_pad) - tof_intercept) / tof_slope;       // Tof2MzConverter::invert
    const double t1 = (std::sqrt(mz_hi + p.dia_ms1_mz_pad) - tof_intercept) / tof_slope;
    GateBox b;
    b.tof_lo = detail::asU32(std::fmax(std::floor(std::fmin(t0, t1)), 0.0));                  // f64::min/max ignore NaN, as fmin/fmax
    b.tof_hi = detail::asU32(std::fmax(std::ceil(std::fmax(t1, t0)), 0.0));
    // Scan2ImConverter::convert as two separate IEEE operations; the binary matches only the unfused result. The
    // product passes through a volatile, a barrier no optimiser can fuse across: clang does not honour its pragma
    // under an explicit -ffp-contract=fast.
    volatile double scaled0 = scan_slope * static_cast<double>(static_cast<std::uint32_t>(r.scan_num_begin));
    volatile double scaled1 = scan_slope * static_cast<double>(static_cast<std::uint32_t>(r.scan_num_end));
    const double im0 = scan_intercept + scaled0;
    const double im1 = scan_intercept + scaled1;
    const double s0 = (std::fmax(im0, im1) + p.dia_ms1_im_pad - scan_intercept) / scan_slope;  // Scan2ImConverter::invert
    const double s1 = (std::fmin(im0, im1) - p.dia_ms1_im_pad - scan_intercept) / scan_slope;
    b.scan_lo = detail::asU32(std::fmax(std::floor(std::fmin(s0, s1)), 0.0));
    b.scan_hi = std::min(detail::asU32(std::fmax(std::ceil(std::fmax(s0, s1)), 0.0)), num_scans - 1);
    boxes.push_back(b);
  }
  out.gate = Ms1WindowGate::fromBoxes(boxes, num_scans);
  return out;
}

/// [dnoise] Point recovery. The Bruker loader delivers m/z, a float32 1/K0 and a corrected intensity, not the (scan, TOF
/// bin, raw count) dnoise filters, so the tool (DnoiseRun in diaspextract.cpp) inverts all three per point and refuses the
/// run on anything that does not recover exactly. The arithmetic lives here -- dependency-free, contraction off, pinned by
/// tests/test_dnoise_ms1.cpp against the loader's forward maths on both benchmark files' metadata -- and the caller turns a
/// Refusal into its message.
namespace recover
{
inline constexpr double kTofTol = 1e-6;   ///< bins: the round trip is exact to 2.3e-10 on both benchmark files
inline constexpr float kImTol = 1e-6f;    ///< 1/K0: a float32 step is ~1.2e-7 there, adjacent scans >= 8.8e-4 apart

/// opentims' intensity correction of a frame, 100.0 / AccumulationTime (opentims.cpp:74), from the column's whole sqlite
/// text. opentims converts that text with atof inside its read_sql locale guard, so in the C locale; detail::parseF64 is
/// the same correctly rounded conversion on every text it accepts, whatever LC_NUMERIC the process has. false unless the
/// text is a finite decimal number > 0.
inline bool opentimsCorrection(const std::string& text, double& corr)
{
  double at = 0.0;
  if (!detail::parseF64(text, at) || !(at > 0.0)) return false;
  corr = 100.0 / at;
  return true;
}

/// Why the raw counts of MS1 frame @p fid cannot be recovered through its correction @p corr, or empty when they can. At
/// corr >= 1 (AccumulationTime <= 100 ms) raw counts below 2^24 do not collide (checked exhaustively at 19 values, 1 + 1 ulp
/// among them); below 1 they do, from raw 3 at 125 ms and from raw 1191 at 100.042 ms. 1e6 is a plausibility bound.
inline std::string correctionRefusal(std::size_t fid, double corr)
{
  if (corr >= 1.0 && corr < 1e6) return std::string();
  std::ostringstream m;   // the classic locale, as for every number this header reads
  m.imbue(std::locale::classic());
  m << "MS1 frame " << fid << ": ";
  if (corr < 1.0)
    m << "AccumulationTime " << 100.0 / corr << " ms exceeds 100 ms, so opentims' intensity correction 100/" << 100.0 / corr << " = "
      << corr << " < 1 maps several raw counts to one intensity and the raw counts dnoise filters cannot be recovered";
  else
    m << "opentims' intensity correction 100 / AccumulationTime = " << corr << " is not invertible";
  return m.str();
}

/// The float32 1/K0 of each of @p n_scans scans under one TimsCalibration row @p c (C0..C9) as the loader converts them --
/// RationalScan2ImConverter::applyFormula term for term, then frameToSpectrum's float cast -- into @p asc, ascending;
/// @p falls: 1/K0 falls with the scan, so @p asc is reversed. Empty, or why scanOf cannot be exact over these scans: the
/// values must rise strictly by more than 4 * kImTol, so no float32 lies within kImTol of two of them. 0 scans give an
/// empty table, which matches nothing. After a refusal @p asc and @p falls mean nothing.
inline std::string buildScanTable(const std::array<double, 10>& c, std::size_t n_scans, std::vector<float>& asc, bool& falls)
{
#if defined(__clang__)
  #pragma clang fp contract(off)
#endif
  asc.assign(n_scans, 0.0f);
  for (std::size_t s = 0; s < n_scans; ++s)
  {
    double V = c[2] + ((c[3] - c[2]) / c[1]) * ((double)s - c[4] - c[0]);
    if (V == 0.0) V = 1e-10;
    double denom = c[6] + c[7] / V;
    if (denom == 0.0) denom = 1e-10;
    asc[s] = static_cast<float>(1.0 / denom);
  }
  falls = n_scans > 1 && asc[1] < asc[0];
  if (falls) std::reverse(asc.begin(), asc.end());
  for (std::size_t s = 1; s < n_scans; ++s)
    if (!(asc[s] - asc[s - 1] > 4.0f * kImTol))   // strictly monotonic, finite, no two scans within reach of one value
      return "the TimsCalibration model does not separate every scan by more than " + std::to_string(4.0f * kImTol) + " 1/K0";
  return std::string();
}

/// The scan whose float32 1/K0 in @p asc (reversed when @p falls) lies within kImTol of @p v; -1 for none.
inline long scanOf(const std::vector<float>& asc, bool falls, float v)
{
  if (asc.empty()) return -1;
  const float* b = asc.data();
  for (std::size_t n = asc.size(); n > 1; ) { const std::size_t h = n / 2; if (b[h] <= v) b += h; n -= h; }   // the last entry <= v, else the first
  std::size_t i = (std::size_t)(b - asc.data());
  if (i + 1 < asc.size() && std::fabs(asc[i + 1] - v) < std::fabs(asc[i] - v)) ++i;
  if (!(std::fabs(asc[i] - v) <= kImTol)) return -1;
  return falls ? (long)(asc.size() - 1 - i) : (long)i;
}

/// Why a spectrum cannot be denoised as a raw frame, or one of its points not recovered.
enum class Refusal
{
  none,
  frame,           ///< its frame id is 0, past the Frames table, or not an MS1 frame (MsMsType 0)
  peaks,           ///< it holds other than the frame's NumPeaks points: not the whole raw frame
  im,              ///< it has points but no 1/K0 per point
  tof,             ///< the m/z is not the MzCalibration table model's value of a TOF bin (within kTofTol)
  scan,            ///< the 1/K0 is not the TimsCalibration model's value of one of the frame's scans (within kImTol)
  intensity,       ///< the intensity is not opentims' correction of exactly one raw count
  intensity_range  ///< the intensity is >= 2^24 = 16777216 (inf included): its float32 no longer identifies one raw count
};

/// The checks on a spectrum before any of its points, in this order. @p frames is the run's Frames table by Id (any row
/// type with msms_type and num_peaks: TdfDnoiseInputs::Frame in the tool), @p n the spectrum's points, @p has_im whether it
/// carries a 1/K0 array and @p im_n that array's size (ignored without one). A frame without points needs no 1/K0 array.
template <class FrameRow>
inline Refusal frameCheck(const std::vector<FrameRow>& frames, std::uint32_t fid, std::size_t n, bool has_im, std::size_t im_n)
{
  if (fid == 0 || fid >= frames.size() || frames[fid].msms_type != 0) return Refusal::frame;
  if (n != frames[fid].num_peaks) return Refusal::peaks;
  if (n != 0 && (!has_im || im_n != n)) return Refusal::im;
  return Refusal::none;
}

/// One point of a frame with factor @p b (Cal::mzToTof inverts the m/z model the loader calibrated with:
/// TdfMzCalibration in the tool), correction @p corr (correctionRefusal empty, so the search visits at most five raw
/// counts), scan table @p asc / @p falls (buildScanTable) and @p num_scans scans: the TOF bin, scan and raw count the
/// loader turned into @p mz, @p im and @p v. Checked in that order; on Refusal::none @p tof, @p scan and @p raw hold the
/// point. For a message: @p t is the fractional TOF, @p hits the raw counts whose correction is @p v (0 until the
/// intensity is checked).
template <class Cal>
inline Refusal recoverPoint(const Cal& cal, double b, double mz, float im, float v, double corr, const std::vector<float>& asc,
                            bool falls, std::uint32_t num_scans, std::uint32_t& tof, std::uint32_t& scan, std::uint32_t& raw,
                            double& t, int& hits)
{
#if defined(__clang__)
  #pragma clang fp contract(off)
#endif
  hits = 0;
  t = cal.mzToTof(mz, b);
  const double bin = std::nearbyint(t);
  if (!(std::fabs(t - bin) <= kTofTol && bin >= 0.0 && bin <= 4294967295.0)) return Refusal::tof;
  const long s = scanOf(asc, falls, im);
  if (s < 0 || (std::uint32_t)s >= num_scans) return Refusal::scan;
  if (v >= 16777216.0f) return Refusal::intensity_range;   // a float32 holds every count below 2^24 exactly, and no more
  std::uint32_t found = 0;
  int n = 0;
  if (v >= 0.0f && v == std::floor(v))
    for (double r = std::max(0.0, std::floor((v - 0.5) / corr) - 1.0), r_end = std::ceil((v + 0.5) / corr) + 1.0; r <= r_end; r += 1.0)
      if ((std::uint32_t)(r * corr + 0.5) == (std::uint32_t)v) { found = (std::uint32_t)r; ++n; }   // opentims.cpp:225: double(raw) * correction + 0.5, into a uint32
  hits = n;
  if (n != 1) return Refusal::intensity;
  tof = (std::uint32_t)bin; scan = (std::uint32_t)s; raw = found;
  return Refusal::none;
}
} // namespace recover

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

/// Points left after each stage of one frame, for logging.
struct FrameCounts
{
  std::uint64_t raw = 0;
  std::uint64_t after_streak = 0;
  std::uint64_t after_halo = 0;   ///< == after_streak with the halo off
  std::uint64_t kept = 0;         ///< after the window gate; == after_halo without one
};

/// Working memory of denoiseMs1Frame, one per thread, reused so a thread's frames do not reallocate.
/// Its contents between calls mean nothing.
struct Scratch
{
  std::vector<std::uint32_t> idx, idx2;             ///< radix sort by TOF
  std::vector<std::uint32_t> tof, scan, inten, id;  ///< live points in TOF order; id = input position
  std::vector<std::uint8_t> mark;
  std::vector<std::uint64_t> sum;                   ///< per-scan window sums, valid where stamp == the window
  std::vector<std::uint32_t> stamp;
  std::vector<std::uint32_t> touched;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> spans;
  std::vector<std::size_t> row, cursor;             ///< halo: survivors of scan s at [row[s], row[s + 1])
  std::vector<std::uint32_t> htof, hinten, hpos;
};

namespace detail
{
/// Drop the live points (TOF-ordered arrays of @p w) whose mark is 0; returns how many remain.
inline std::size_t compact(Scratch& w, std::size_t m)
{
  std::size_t k = 0;
  for (std::size_t j = 0; j < m; ++j)
    if (w.mark[j])
    {
      w.tof[k] = w.tof[j]; w.scan[k] = w.scan[j]; w.inten[k] = w.inten[j]; w.id[k] = w.id[j];
      ++k;
    }
  return k;
}

/// One streak pass (filter.rs:14-129) over the first @p m live points; returns the survivors left in place.
inline std::size_t streakPass(Scratch& w, std::size_t m, std::size_t num_scans, const Params& p)
{
  const std::uint32_t* const t = w.tof.data();
  const std::uint32_t* const s = w.scan.data();
  const std::uint32_t* const v = w.inten.data();
  const std::uint32_t hw = p.mz_half_width, u32max = std::numeric_limits<std::uint32_t>::max();
  w.mark.assign(m, 0);
  w.stamp.assign(num_scans, 0);
  w.sum.resize(num_scans);
  std::uint32_t window = 0;   // at most one window per live point, and m < 2^32
  std::size_t lo = 0, hi = 0;
  for (std::size_t k = 0; k < m;)
  {
    const std::uint32_t c = t[k];
    std::size_t kb = k + 1;
    while (kb < m && t[kb] == c) ++kb;   // the block of points at TOF c; only these are decided here
    // [c - w, c + w], saturating (filter.rs:48-51). Both ends only grow with c, so two pointers replace the
    // binary searches.
    const std::uint32_t lo_val = c >= hw ? c - hw : 0, hi_val = c <= u32max - hw ? c + hw : u32max;
    while (t[lo] < lo_val) ++lo;
    while (hi < m && t[hi] <= hi_val) ++hi;
    ++window;
    w.touched.clear();
    for (std::size_t j = lo; j < hi; ++j)
    {
      const std::uint32_t sc = s[j];
      if (w.stamp[sc] != window) { w.stamp[sc] = window; w.sum[sc] = 0; w.touched.push_back(sc); }
      w.sum[sc] += v[j];
    }
    auto occupied = [&w, &p](std::uint32_t sc) { return w.sum[sc] > 0 && w.sum[sc] >= p.min_window_intensity; };
    std::size_t n_occupied = 0;
    for (const std::uint32_t sc : w.touched) n_occupied += occupied(sc);
    // A run holds at most every occupied scan of the window, so fewer than min_feature_length accept nothing.
    if (n_occupied != 0 && n_occupied >= p.min_feature_length)
    {
      std::sort(w.touched.begin(), w.touched.end());
      w.spans.clear();
      std::size_t first = 0, last = 0, count = 0;   // positions in touched of the open run's first/last occupied scan
      auto close = [&]() {
        if (count < p.min_feature_length) return;   // occupied scans, not the span (filter.rs:86-92)
        if (p.min_feature_intensity > 0)            // the span's sum, sub-floor scans included (filter.rs:93-96)
        {
          std::uint64_t total = 0;
          for (std::size_t q = first; q <= last; ++q) total += w.sum[w.touched[q]];
          if (total < p.min_feature_intensity) return;
        }
        w.spans.emplace_back(w.touched[first], w.touched[last]);
      };
      bool open = false;
      for (std::size_t q = 0; q < w.touched.size(); ++q)
      {
        if (!occupied(w.touched[q])) continue;
        if (open && static_cast<std::size_t>(w.touched[q] - w.touched[last]) > p.max_internal_gap + 1) { close(); open = false; }
        if (!open) { open = true; first = q; count = 0; }
        last = q;
        ++count;
      }
      if (open) close();
      if (!w.spans.empty())
        for (std::size_t j = k; j < kb; ++j)   // kept iff the scan lies in an accepted run's span (filter.rs:108-118)
        {
          const std::uint32_t sc = s[j];
          const auto sp = std::partition_point(w.spans.begin(), w.spans.end(),
                                               [sc](const std::pair<std::uint32_t, std::uint32_t>& r) { return r.second < sc; });
          if (sp != w.spans.end() && sc >= sp->first) w.mark[j] = 1;
        }
    }
    k = kb;
  }
  return compact(w, m);
}
} // namespace detail

/// dnoise v0.1.0's default MS1 denoising of ONE frame (writer.rs:838-853). Inputs are the frame's points in any
/// order: 0-based scan index (< @p num_scans, the frame's NumScans), TOF index, and RAW blob intensity. @p gate is
/// the run's buildMs1WindowGate() result, or nullptr for none. On return keep[i] == 1 iff dnoise keeps point i.
/// A pure function of its inputs: @p w is working memory only, so any thread may process any frame.
/// Throws std::invalid_argument for a scan index >= num_scans (dnoise panics there), and -- before any allocation, and
/// before @p keep or an input point is touched -- for more than 2^32-1 points or scans: points are indexed as u32, and
/// num_scans <= 2^32-1 keeps `num_scans + 1` and `scan + 1` from wrapping (OpenMS is 64-bit only).
inline FrameCounts denoiseMs1Frame(const std::uint32_t* scan, const std::uint32_t* tof, const std::uint32_t* intensity,
                                   std::size_t n, std::size_t num_scans, const Params& p, const Ms1WindowGate* gate,
                                   Scratch& w, std::vector<std::uint8_t>& keep)
{
  FrameCounts c;
  c.raw = n;
  if (n > std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("dnoise: more than 2^32-1 points in one frame");
  if (num_scans > std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("dnoise: more than 2^32-1 scans in one frame");
  keep.assign(n, 0);
  if (n == 0) return c;
  std::uint32_t max_tof = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    if (scan[i] >= num_scans)
      throw std::invalid_argument("dnoise: scan index " + std::to_string(scan[i]) + " outside a frame of " + std::to_string(num_scans) + " scans");
    max_tof = std::max(max_tof, tof[i]);
  }

  // Points in TOF order (stable LSD radix sort, 11-bit digits); set semantics make the tie order irrelevant.
  w.idx.resize(n);
  w.idx2.resize(n);
  for (std::size_t i = 0; i < n; ++i) w.idx[i] = static_cast<std::uint32_t>(i);
  for (unsigned shift = 0; shift < 32 && (shift == 0 || (max_tof >> shift) != 0); shift += 11)
  {
    std::array<std::size_t, 2049> start{};
    for (std::size_t k = 0; k < n; ++k) ++start[((tof[w.idx[k]] >> shift) & 2047u) + 1];
    for (std::size_t d = 0; d < 2048; ++d) start[d + 1] += start[d];
    for (std::size_t k = 0; k < n; ++k) w.idx2[start[(tof[w.idx[k]] >> shift) & 2047u]++] = w.idx[k];
    w.idx.swap(w.idx2);
  }
  w.id.swap(w.idx);
  w.tof.resize(n);
  w.scan.resize(n);
  w.inten.resize(n);
  for (std::size_t k = 0; k < n; ++k)
  {
    const std::uint32_t i = w.id[k];
    w.tof[k] = tof[i]; w.scan[k] = scan[i]; w.inten[k] = intensity[i];
  }

  // 1. Streak filter (filter.rs:133-165). A pass that removes nothing is a fixed point: every further pass
  //    would see the same input, so stopping there is exact.
  std::size_t m = n;
  for (std::size_t pass = 0; pass < p.iterations && m > 0; ++pass)
  {
    const std::size_t before = m;
    m = detail::streakPass(w, m, num_scans, p);
    if (m == before) break;
  }
  c.after_streak = m;

  // 2. Halo (halo.rs:24-101), once, with the streak survivors as both the points judged and the reference;
  //    simultaneous, so a point it drops still counts as a neighbour's reference.
  if (p.halo && m > 0)
  {
    w.row.assign(num_scans + 1, 0);
    for (std::size_t j = 0; j < m; ++j) ++w.row[w.scan[j] + 1];
    for (std::size_t r = 0; r < num_scans; ++r) w.row[r + 1] += w.row[r];
    w.cursor.assign(w.row.begin(), w.row.end() - 1);
    w.htof.resize(m);
    w.hinten.resize(m);
    w.hpos.resize(m);
    for (std::size_t j = 0; j < m; ++j)   // stable by scan over TOF order: (scan, TOF) order
    {
      const std::size_t d = w.cursor[w.scan[j]]++;
      w.htof[d] = w.tof[j]; w.hinten[d] = w.inten[j]; w.hpos[d] = static_cast<std::uint32_t>(j);
    }
    w.mark.assign(m, 0);
    const std::size_t sw = p.halo_scan_half_width;
    const std::uint32_t mw = p.halo_mz_idx_half_width, u32max = std::numeric_limits<std::uint32_t>::max();
    const std::uint32_t* const ht = w.htof.data();
    for (std::size_t r = 0; r < num_scans; ++r)
      for (std::size_t i = w.row[r]; i < w.row[r + 1]; ++i)
      {
        const std::uint32_t ti = ht[i];
        const std::size_t v_lo = r >= sw ? r - sw : 0, v_hi = sw >= num_scans - 1 - r ? num_scans - 1 : r + sw;
        const std::uint32_t t_lo = ti >= mw ? ti - mw : 0, t_hi = ti <= u32max - mw ? ti + mw : u32max;
        std::uint32_t best = 0;
        for (std::size_t vs = v_lo; vs <= v_hi; ++vs)
        {
          const std::uint32_t* const e = ht + w.row[vs + 1];
          for (const std::uint32_t* q = std::lower_bound(ht + w.row[vs], e, t_lo); q != e && *q <= t_hi; ++q)
            if (*q != ti && w.hinten[static_cast<std::size_t>(q - ht)] > best) best = w.hinten[static_cast<std::size_t>(q - ht)];
        }
        // The Rust double compare, verbatim: no addition, so nothing to contract; best == 0 keeps.
        if (static_cast<double>(w.hinten[i]) >= p.halo_peak_fraction * static_cast<double>(best)) w.mark[w.hpos[i]] = 1;
      }
    m = detail::compact(w, m);
  }
  c.after_halo = m;

  // 3. The window gate, ANDed per point (writer.rs:846-853).
  for (std::size_t j = 0; j < m; ++j)
    if (gate == nullptr || gate->contains(w.scan[j], w.tof[j]))
    {
      keep[w.id[j]] = 1;
      ++c.kept;
    }
  return c;
}

/// [dnoise] That a run denoised every MS1 frame with points exactly once. Counts cannot tell: a frame denoised twice and a
/// missing one with equal NumPeaks cancel out. So one mark per frame id -- set serially as each batch's counts are merged,
/// checked once after the load, and reset with the counts when a failed stream load falls back to the resident one.
class FrameCoverage
{
public:
  /// Mark frame @p fid (an id frameCheck accepted) as denoised; false if it already was.
  bool mark(std::uint32_t fid)
  {
    if (fid >= seen_.size()) seen_.resize(static_cast<std::size_t>(fid) + 1, 0);
    if (seen_[fid]) return false;
    seen_[fid] = 1;
    return true;
  }
  /// The first MS1 frame with points in @p frames (by Id; any row type with msms_type and num_peaks) that was not marked,
  /// or 0 when every one was.
  template <class FrameRow>
  std::size_t firstMissing(const std::vector<FrameRow>& frames) const
  {
    for (std::size_t id = 1; id < frames.size(); ++id)
      if (frames[id].msms_type == 0 && frames[id].num_peaks > 0 && (id >= seen_.size() || !seen_[id])) return id;
    return 0;
  }
  void reset() { seen_.clear(); }

private:
  std::vector<std::uint8_t> seen_;   ///< by frame id
};

} // namespace spx::dnoise
