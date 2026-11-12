// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

// The extraction path for DIA input without ion mobility (Orbitrap, SWATH mzML). Separate from the IMS path by rule:
// it shares loading and writing with it, nothing else, and never changes what the IMS path emits.
//
// Precursors: Biosaur2 features of the whole run's MS1, plus isotope envelopes (z >= 2) of unfragmented precursors in
// each window's MS2 traces that MS1 does not already have.
// Per isolation window: (1) the window's precursors; (2) fragments = mass traces of the window's MS2 scans; (3) a
// fragment joins a precursor when its profile correlates with the precursor's over the precursor's elution span
// (Pearson r >= min_corr on the window's scan grid), or, with fewer than 3 grid points in that span, when its apex lies
// within one cycle of the precursor apex. A fragment's intensity is weighted by r^24 (Sage and MSFragger peptides and
// entrapment on three Orbitrap sets: stronger powers overfit Sage's scoring); the 200 heaviest are kept, one spectrum
// per precursor with >= min_frags (3) fragments.
// Parallelism: windows run in parallel, heaviest first; each window's precursors are OpenMP tasks, so threads that have
// no window left help with the precursors of the last ones. Results are index-addressed: the output does not depend on
// the thread count.
#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include "PortableMath.h"
#include "OmsString.h"
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

/// @brief The extraction path for DIA input without ion mobility.
namespace DIAspeXtract::NoIms
{

/// @brief Parameters of the non-IMS extraction (precursor detection, fragment traces, fragment-precursor grouping).
struct Params
{
  double ppm = 10.0;     ///< mass-trace and isotope m/z tolerance
  double min_corr = 0.3; ///< fragment-precursor Pearson r on the window grid
  // MS1 precursor detection (Biosaur2): the Orbitrap defaults; TOF input changes them (see noImsSpectra_)
  bool tof = false;              ///< Biosaur2's TOF noise filter (per-50-m/z-bin intensity floor)
  double min_intensity = 1.0;    ///< Biosaur2 mini: least peak intensity
  int min_hill_scans = 2;        ///< Biosaur2 minlh: fewest MS1 scans per hill
  int min_charge = 1;            ///< Biosaur2 cmin
  int min_frags = 3;             ///< fewest fragments a spectrum needs
  double min_frag_intensity = 0; ///< least apex intensity of a fragment trace
  int min_frag_points = 1;       ///< fewest scans of a fragment trace
  int max_frags = 200;           ///< most intense fragments kept per spectrum (by weighted intensity)
  double weight_pow = 24.0;      ///< fragment intensity x r^weight_pow
  double ms2_ppm = 0;            ///< MS2 mass-trace tolerance; 0 = ppm
  double split = 0;              ///< split MS2 traces at elution valleys below split x the lower flanking maximum; 0 = off
  bool env_all_charges = false;  ///< MS2 envelopes: a -1 isotope at one charge moves on to the next lower charge (else: no precursor)
  bool drop_precursor = false;   ///< drop the precursor's own isotope peaks (mono-1 .. mono+3) from its spectrum
  double comp_boost = 0;         ///< weight x this for fragments with a complementary partner (m1 + m2 = M + 2 H+, both z = 1); 0 = off
  bool env_pattern = false;      ///< MS2 envelopes: +1/mono intensity must fit averagine within x3; the best-fitting charge wins
  int hill_seeds = 0;            ///< fragment-driven tier: MS1 hills outside every feature become z 2 and 3 precursors, kept with at least this many
                                 ///< complementary fragment pairs; 0 = off
  bool tag_src = false;   ///< write userParam precursor_source (0 MS1 feature, 1 MS2 envelope, 3 MS1 hill with fragment pairs); hill seeds also get
                          ///< hill_pairs and hill_span (s)
  int pair_top = 50;      ///< hill seeds: complementary pairs are counted among this many most intense peaks
  bool seed_mono = false; ///< hill seeds: skip a charge when a hill one isotope step lighter overlaps in RT (not monoisotopic)
  double ms2_shift = 0;   ///< fragment m/z shift in ppm
  bool ms2_shift_auto = false; ///< estimate the shift from the y1 ions of tryptic peptides
  int span_pad = 0;            ///< widen the precursor's RT span by this many window scans on each side (correlation and apex gate); 0 = off
};

/// @brief Counts and timings of one run(), for the log.
struct Stats
{
  OpenMS::Size windows = 0, ms2_traces = 0, precursors = 0, spectra = 0,
               peaks = 0;        ///< isolation windows, MS2 mass traces, precursors, spectra written, peaks in them
  double busy_s = 0, wall_s = 0; ///< summed per-window busy time and wall time of run() (s)
  int threads = 1;               ///< OpenMP threads available (omp_get_max_threads)
};

/// @brief Building blocks of run(); not part of the tool's interface.
namespace Internal
{
  constexpr double C13 = 1.0033548, PROTON = 1.007276; ///< 13C-12C mass difference and proton mass (Da)
  inline double ms2_scale = 1.0;                       ///< fragment m/z recalibration factor, set once per run (noims:ms2_shift)
  /**
    @brief Keeps the first exception of a parallel region (call from a catch handler).
    @param[in,out] e Set to the current exception unless it already holds one.
  */
  inline void keepFirst(std::exception_ptr& e)
  {
#pragma omp critical(noims_err)
    if (! e) { e = std::current_exception(); }
  }
  constexpr double FEW_APEX_CYCLES = 1.0; ///< precursors with < 3 window scans: fragment apex within this many cycles of the precursor apex

  /// @brief A precursor candidate with its MS1 (or MS2 envelope) elution profile.
  struct Prec
  {
    double mz, rt, rt_lo, rt_hi, inten;  ///< m/z, apex RT, RT span (first and last point), apex intensity
    int z;                               ///< charge (0 = open or skipped, on hill seeds)
    std::vector<double> rt_pts, int_pts; ///< elution profile: RT and intensity of each point
    int src = 0;                         ///< source: 0 MS1 feature, 1 MS2 envelope, 3 MS1 hill seed
  };
  /**
    @brief A fragment trace.

    `pts` are (RT, intensity), `gpts` the same points as (grid index, intensity).
  */
  struct Frag
  {
    double mz, apex_rt, inten;                         ///< m/z, RT of the most intense point, its intensity
    std::vector<std::pair<double, double>> pts;        ///< (RT, intensity) points
    std::vector<std::pair<OpenMS::Size, double>> gpts; ///< (window grid index, intensity) points
  };

  /**
    @brief The RT of a trace's most intense point.
    @param[in] t The mass trace.
    @return The apex RT; the centroid RT for an empty trace.
  */
  inline double apexRT(const OpenMS::MassTrace& t)
  {
    double best = -1, rt = t.getCentroidRT();
    for (const auto& p : t)
    {
      if (p.getIntensity() > best)
      {
        best = p.getIntensity();
        rt = p.getRT();
      }
    }
    return rt;
  }

  /**
    @brief Mass traces of a map with OpenMS MassTraceDetection (no progress log, single-point traces kept).
    @param[in] m The spectra to trace (MS level 1 as MassTraceDetection reads it).
    @param[in] ppm Mass error tolerance (ppm).
    @return The detected mass traces.
  */
  inline std::vector<OpenMS::MassTrace> traces(const OpenMS::PeakMap& m, double ppm)
  {
    OpenMS::MassTraceDetection mtd;
    mtd.setLogType(OpenMS::ProgressLogger::NONE); // windows run in parallel: progress lines would interleave
    OpenMS::Param p = mtd.getDefaults();
    p.setValue("mass_error_ppm", ppm);
    p.setValue("min_trace_length", 0.0); // single points too: a fragment seen in one scan of a window
    mtd.setParameters(p);
    std::vector<OpenMS::MassTrace> mt;
    mtd.run(m, mt);
    return mt;
  }

  /**
    @brief Cut points (grid indices) of a trace's grid profile @p y (index g0 + k).

    Cuts at valleys of the 1-2-1 smoothed profile that fall below @p valley x both flanking maxima: one elution peak
    per segment.
    @param[in] y The profile on consecutive grid indices.
    @param[in] g0 Grid index of y[0].
    @param[in] valley Relative depth a valley must fall below.
    @return The grid indices of the valley minima, ascending.
  */
  inline std::vector<OpenMS::Size> valleys(const std::vector<double>& y, OpenMS::Size g0, double valley)
  {
    std::vector<OpenMS::Size> cut;
    double cur_max = 0, m = -1;
    OpenMS::Size mi = 0; // m < 0: still rising to the current maximum M
    for (OpenMS::Size k = 0; k < y.size(); ++k)
    {
      const double v = 0.25 * ((k ? y[k - 1] : 0.0) + 2 * y[k] + (k + 1 < y.size() ? y[k + 1] : 0.0));
      if (m < 0)
      {
        if (v >= cur_max) { cur_max = v; }
        else
        {
          m = v;
          mi = k;
        }
      }
      else if (v < m)
      {
        m = v;
        mi = k;
      }
      else if (m < valley * cur_max && m < valley * v)
      {
        cut.push_back(g0 + mi);
        cur_max = v;
        m = -1;
      }
      else if (v > cur_max)
      {
        cur_max = v;
        m = -1;
      }
    }
    return cut;
  }

  /**
    @brief Pearson r of two profiles sampled on the same grid.
    @param[in] a First profile.
    @param[in] b Second profile, same size as @p a.
    @return The correlation coefficient; 0 if either profile is constant.
  */
  inline double pearson(const std::vector<double>& a, const std::vector<double>& b)
  {
    const OpenMS::Size n = a.size();
    double ma = 0, mb = 0;
    for (OpenMS::Size i = 0; i < n; ++i)
    {
      ma += a[i];
      mb += b[i];
    }
    ma /= n;
    mb /= n;
    double sab = 0, saa = 0, sbb = 0;
    for (OpenMS::Size i = 0; i < n; ++i)
    {
      sab += (a[i] - ma) * (b[i] - mb);
      saa += (a[i] - ma) * (a[i] - ma);
      sbb += (b[i] - mb) * (b[i] - mb);
    }
    return (saa > 0 && sbb > 0) ? sab / std::sqrt(saa * sbb) : 0.0;
  }

  /**
    @brief Linear interpolation of (x, y) at t, 0 outside.
    @param[in] x Sample positions, ascending.
    @param[in] y Sample values, same size as @p x.
    @param[in] t Position to evaluate.
    @return The interpolated value, 0 outside [x.front(), x.back()].
  */
  inline double interp(const std::vector<double>& x, const std::vector<double>& y, double t)
  {
    if (x.empty() || t < x.front() || t > x.back()) { return 0.0; }
    const OpenMS::Size j = std::upper_bound(x.begin(), x.end(), t) - x.begin(); // >= 1, since t >= x.front()
    if (j >= x.size()) { return y.back(); }
    const double f = (t - x[j - 1]) / (x[j] - x[j - 1]);
    return y[j - 1] + f * (y[j] - y[j - 1]);
  }

  /**
    @brief Isotope envelopes (z = 2..4) from MS2 traces.

    A trace is a mono if it has a +1 isotope and no -1 isotope at the same charge.
    @param[in] mt The window's MS2 mass traces.
    @param[in] ppm Isotope m/z tolerance (ppm).
    @param[in] lo Lower bound of the mono m/z (inclusive).
    @param[in] hi Upper bound of the mono m/z (exclusive).
    @param[in] coelute Least apex RT tolerance for isotope partners (s); half the trace FWHM where larger.
    @param[in] all_charges A -1 isotope at one charge moves on to the next lower charge instead of rejecting the trace.
    @param[in] pattern Require the +1/mono intensity ratio to fit averagine within x3; the best-fitting charge wins.
    @return One precursor (src 1) per mono trace found.
  */
  // ponytail: greedy, highest charge first, no averagine check; add an intensity-pattern test if charge errors show up
  inline std::vector<Prec> envelopes(const std::vector<OpenMS::MassTrace>& mt,
                                     double ppm,
                                     double lo,
                                     double hi,
                                     double coelute,
                                     bool all_charges = false,
                                     bool pattern = false)
  {
    std::vector<OpenMS::Size> o(mt.size());
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](OpenMS::Size a, OpenMS::Size b) {
      return mt[a].getCentroidMZ() < mt[b].getCentroidMZ() || (mt[a].getCentroidMZ() == mt[b].getCentroidMZ() && a < b);
    });
    std::vector<double> mz(o.size()), art(o.size());
    for (OpenMS::Size k = 0; k < o.size(); ++k)
    {
      mz[k] = mt[o[k]].getCentroidMZ();
      art[k] = apexRT(mt[o[k]]);
    }
    auto partner = [&](OpenMS::Size k, double target) -> long { // co-eluting trace at target m/z
      const double tol = target * ppm * 1e-6;
      const OpenMS::MassTrace& a = mt[o[k]];
      const double half = std::max(a.getFWHM() * 0.5, coelute); // isotopes share scans: one scan of apex jitter
      long best = -1;
      double bi = -1;
      for (auto it = std::lower_bound(mz.begin(), mz.end(), target - tol); it != mz.end() && *it <= target + tol; ++it)
      {
        const OpenMS::Size j = it - mz.begin();
        if (std::fabs(art[j] - art[k]) > half) { continue; }
        const double in = mt[o[j]].getMaxIntensity(false);
        if (in > bi)
        {
          bi = in;
          best = (long)j;
        }
      }
      return best;
    };
    std::vector<Prec> out;
    auto emit = [&](OpenMS::Size k, int z) {
      const OpenMS::MassTrace& t = mt[o[k]];
      Prec p {mz[k], art[k], t.begin()->getRT(), t.rbegin()->getRT(), t.getMaxIntensity(false), z, {}, {}, 1};
      for (const auto& pk : t)
      {
        p.rt_pts.push_back(pk.getRT());
        p.int_pts.push_back(pk.getIntensity());
      }
      out.push_back(std::move(p));
    };
    auto in = [&](long j) { return (double)mt[o[j]].getMaxIntensity(false); };
    for (OpenMS::Size k = 0; pattern && k < o.size(); ++k)
    {
      if (! (mz[k] >= lo && mz[k] < hi)) { continue; }
      // ponytail: averagine as one number, I(+1)/I(mono) = 5.3e-4 x mass, within x3; fit the +2 peak too if charges stay wrong
      const double LN3 = 1.0986;
      int bz = 0;
      double bs = LN3;
      for (int z = 4; z >= 2; --z)
      {
        const long p1 = partner(k, mz[k] + C13 / z);
        if (p1 < 0) { continue; }
        const double ex = 5.3e-4 * z * (mz[k] - PROTON);
        const long m1 = partner(k, mz[k] - C13 / z);
        if (m1 >= 0 && std::fabs(std::log(in(k) / in(m1) / ex)) < LN3)
        {
          continue; // this trace is the +1 isotope at z
        }
        const double sc = std::fabs(std::log(in(p1) / in(k) / ex));
        if (sc < bs)
        {
          bs = sc;
          bz = z;
        }
      }
      if (bz) { emit(k, bz); }
    }
    for (OpenMS::Size k = 0; ! pattern && k < o.size(); ++k)
    {
      for (int z = 4; z >= 2; --z)
      {
        if (partner(k, mz[k] + C13 / z) < 0) { continue; }
        if (partner(k, mz[k] - C13 / z) >= 0)
        {
          if (all_charges) { continue; }
          break;
        } // not the lightest isotope at this charge
        if (mz[k] >= lo && mz[k] < hi) { emit(k, z); }
        break;
      }
    }
    return out;
  }

  /**
    @brief Precursors from OpenMS's Biosaur2 port, run once on the whole run's MS1 (averagine-checked isotope patterns).
    @param[in,out] ms1 The run's MS1 spectra; moved into the algorithm.
    @param[in] opts Parameters (ppm, min_charge, min_intensity, min_hill_scans, tof).
    @param[in] mzlo Lowest m/z considered.
    @param[in] mzhi Highest m/z considered.
    @param[out] seeds If not null, every hill is appended to it as a charge-open precursor (src 3), sorted by m/z.
    @return One precursor per feature with its monoisotopic hill's profile, in a total order by m/z.
  */
  inline std::vector<Prec> biosaur(OpenMS::PeakMap&& ms1, const Params& opts, double mzlo, double mzhi, std::vector<Prec>* seeds = nullptr)
  {
    std::vector<OpenMS::Biosaur2Algorithm::Hill> hills;
    std::vector<OpenMS::Biosaur2Algorithm::PeptideFeature> feats;
    {
      OpenMS::Biosaur2Algorithm b;
      OpenMS::Param p = b.getDefaults();
      p.setValue("htol", opts.ppm);
      p.setValue("itol", opts.ppm);
      p.setValue("minmz", mzlo);
      p.setValue("maxmz", mzhi);
      p.setValue("cmin", opts.min_charge);
      p.setValue("cmax", 4);
      p.setValue("mini", opts.min_intensity);
      p.setValue("minlh", opts.min_hill_scans);
      p.setValue("tof", opts.tof ? "true" : "false");
      // no per-peak ion mobility on this path: a FAIMS run's per-spectrum CV made Biosaur2 switch on its PASEF
      // centroiding and IM gating, both no-ops without an IM array but logged as if they ran
      p.setValue("paseftol", 0.0);
      b.setParameters(p);
      b.setMSData(std::move(ms1));
      OpenMS::FeatureMap fm; // unused; freed with the algorithm's copy of the MS1 data
      b.run(fm, hills, feats);
    }
    if (seeds) // every hill, charge open; run() drops those a feature explains
    {
      for (const auto& h : hills)
      {
        seeds->push_back(Prec {h.mz_weighted_mean, h.rt_apex, h.rt_start, h.rt_end, h.intensity_apex, 0, h.rt_values, h.intensities, 3});
      }
      std::sort(seeds->begin(), seeds->end(), [](const Prec& a, const Prec& b) {
        return std::tie(a.mz, a.rt, a.inten, a.rt_lo, a.rt_hi, a.rt_pts, a.int_pts)
               < std::tie(b.mz, b.rt, b.inten, b.rt_lo, b.rt_hi, b.rt_pts, b.int_pts);
      });
    }
    OpenMS::Size top = 0;
    for (const auto& h : hills)
    {
      top = std::max(top, h.hill_idx);
    }
    std::vector<OpenMS::Size> at(hills.empty() ? 0 : top + 1, SIZE_MAX); // hill_idx -> position in hills (last one wins, as a map would)
    for (OpenMS::Size i = 0; i < hills.size(); ++i)
    {
      at[hills[i].hill_idx] = i;
    }
    std::vector<Prec> out;
    out.reserve(feats.size());
    for (const auto& f : feats)
    {
      if (f.mono_hill_idx >= at.size() || at[f.mono_hill_idx] == SIZE_MAX) { continue; }
      const auto& h = hills[at[f.mono_hill_idx]];
      out.push_back(Prec {f.mz, f.rt_apex, h.rt_start, h.rt_end, f.intensity_apex, f.charge, h.rt_values, h.intensities});
    }
    std::sort(out.begin(), out.end(), [](const Prec& a, const Prec& b) { // a total order: Biosaur2's feature order is fixed, its ties are not
      return std::tie(a.mz, a.rt, a.z, a.inten, a.rt_lo, a.rt_hi, a.rt_pts, a.int_pts)
             < std::tie(b.mz, b.rt, b.z, b.inten, b.rt_lo, b.rt_hi, b.rt_pts, b.int_pts);
    });
    return out;
  }

  /**
    @brief Counts singly charged fragment pairs among the @p top most intense peaks of @p s that sum to the precursor (b + y = M + 2 H+).
    @param[in] s The spectrum.
    @param[in] neutral Neutral precursor mass M (Da).
    @param[in] ppm Tolerance on the pair sum (ppm).
    @param[in] top Number of most intense peaks considered.
    @return The number of disjoint complementary pairs.
  */
  inline int complementPairs(const OpenMS::MSSpectrum& s, double neutral, double ppm, OpenMS::Size top = 50)
  {
    std::vector<std::pair<double, double>> v; // (intensity, m/z)
    for (const auto& x : s)
    {
      v.emplace_back(x.getIntensity(), x.getMZ());
    }
    // ponytail: a plain count among the top peaks; score by intensity if the count alone is too coarse
    if (v.size() > top)
    {
      std::nth_element(v.begin(), v.begin() + top, v.end(), std::greater<>());
      v.resize(top);
    }
    std::vector<double> m;
    for (const auto& x : v)
    {
      m.push_back(x.second);
    }
    std::sort(m.begin(), m.end());
    const double sum = neutral + 2 * PROTON, tol = sum * ppm * 1e-6;
    int n = 0;
    for (long i = 0, j = (long)m.size() - 1; i < j;)
    {
      const double d = m[i] + m[j] - sum;
      // clang-format off: cpplint wants braces on both sides of else; the OpenMS style would put these blocks on one line
      if (d < -tol)
      {
        ++i;
      }
      else if (d > tol)
      {
        --j;
      }
      else
      {
        ++n;
        ++i;
        --j;
      }
      // clang-format on
    }
    return n;
  }

  /**
    @brief The spectrum of precursor @p p from the window's fragments (sorted by apex RT), or an empty one (no MS level).

    Fragments within the precursor's span are weighted by r^weight_pow of their correlation with the precursor
    profile (or by an apex-distance Gaussian for spans of fewer than 3 grid points); at most max_frags are kept.
    @param[in] p The precursor.
    @param[in] fr The window's fragments, sorted by apex RT.
    @param[in] grid RTs of the window's MS2 scans.
    @param[in] cycle MS2 cycle time (s).
    @param[in] opts Parameters.
    @param[in] lo Lower isolation window bound (m/z).
    @param[in] hi Upper isolation window bound (m/z).
    @param[in,out] raw If not null and non-empty, the weighted fragments to reuse (filtered by this precursor's mass
                       limit); if not null and empty, receives this precursor's weighted fragments.
    @return An MS2 centroid spectrum with one precursor, or an empty spectrum if fewer than min_frags fragments remain.
  */
  inline OpenMS::MSSpectrum assemble(const Prec& p,
                                     const std::vector<Frag>& fr,
                                     const std::vector<double>& grid,
                                     double cycle,
                                     const Params& opts,
                                     double lo,
                                     double hi,
                                     std::vector<std::pair<double, double>>* raw = nullptr)
  {
    OpenMS::MSSpectrum s;
    OpenMS::Size g0 = std::lower_bound(grid.begin(), grid.end(), p.rt_lo) - grid.begin();
    OpenMS::Size g1 = std::upper_bound(grid.begin(), grid.end(), p.rt_hi) - grid.begin();
    g0 -= std::min(g0, (OpenMS::Size)opts.span_pad);
    g1 = std::min(grid.size(), g1 + opts.span_pad); // the precursor is 0 outside its hill
    std::vector<double> pp;
    for (OpenMS::Size g = g0; g < g1; ++g)
    {
      pp.push_back(interp(p.rt_pts, p.int_pts, grid[g]));
    }
    const bool few = pp.size() < 3;
    const double neutral = p.z * (p.mz - PROTON);
    std::vector<std::pair<double, double>> pk; // (intensity, m/z)
    std::vector<double> fp(pp.size());
    // split traces (TOF) are short against the precursor span: their r comes from their own points (one-pass sums).
    // ponytail: unsplit traces keep the two-pass dense form, so Orbitrap output stays bit-identical
    const bool sparse = opts.split > 0;
    double sa = 0, saa = 0;
    for (double a : pp)
    {
      sa += a;
      saa += a * a;
    }
    const double nn = (double)pp.size(), va = nn * saa - sa * sa;
    const bool pad = opts.span_pad > 0 && ! few;
    const double a0 = few   ? p.rt - FEW_APEX_CYCLES * cycle
                      : pad ? grid[g0]
                            : p.rt_lo,
                 a1 = few   ? p.rt + FEW_APEX_CYCLES * cycle
                      : pad ? grid[g1 - 1]
                            : p.rt_hi;
    // raw: the weighted fragments depend on the charge only through the mass limit (drop_precursor off), so a
    // second charge of the same precursor filters the list of the higher charge instead of correlating again
    if (raw && ! raw->empty())
    {
      for (const auto& q : *raw)
      {
        if (q.second - PROTON <= neutral - 56.0) { pk.push_back(q); }
      }
    }
    else
    {
      const auto first = std::lower_bound(fr.begin(), fr.end(), a0, [](const Frag& f, double t) { return f.apex_rt < t; });
      for (OpenMS::Size k = first - fr.begin(); k < fr.size() && fr[k].apex_rt <= a1; ++k)
      {
        const Frag& f = fr[k];
        if (f.mz - PROTON > neutral - 56.0)
        {
          continue; // heavier than the precursor allows at z = 1
        }
        if (opts.drop_precursor)
        {
          const double iso = std::round((f.mz - p.mz) * p.z / C13);
          if (iso >= -1 && iso <= 3 && std::fabs(f.mz - p.mz - iso * C13 / p.z) <= f.mz * opts.ppm * 1e-6)
          {
            continue; // unfragmented precursor
          }
        }
        if (! few)
        {
          double r;
          if (sparse)
          {
            double sb = 0, sbb = 0, sab = 0;
            for (auto it = std::lower_bound(f.gpts.begin(), f.gpts.end(), std::make_pair(g0, -1.0)); it != f.gpts.end() && it->first < g1; ++it)
            {
              sb += it->second;
              sbb += it->second * it->second;
              sab += pp[it->first - g0] * it->second;
            }
            const double vb = nn * sbb - sb * sb;
            r = (va > 0 && vb > 0) ? (nn * sab - sa * sb) / std::sqrt(va * vb) : 0.0;
          }
          else
          {
            std::fill(fp.begin(), fp.end(), 0.0);
            for (const auto& q : f.gpts)
            {
              if (q.first >= g0 && q.first < g1) { fp[q.first - g0] = q.second; }
            }
            r = pearson(pp, fp);
          }
          if (r < opts.min_corr) { continue; }
          pk.emplace_back(f.inten * pmath::pow(r, opts.weight_pow), f.mz);
          continue;
        }
        const double d = (f.apex_rt - p.rt) / cycle;
        pk.emplace_back(f.inten * pmath::exp(-d * d * opts.weight_pow / 4.0), f.mz);
      }
      if (raw) { *raw = pk; }
    }
    if ((int)pk.size() < opts.min_frags) { return s; }
    if (opts.comp_boost > 0)
    {
      std::sort(pk.begin(), pk.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
      const double target = neutral + 2 * PROTON, tol = target * 2 * opts.ppm * 1e-6;
      std::vector<char> hit(pk.size(), 0);
      for (OpenMS::Size i = 0, j = pk.size(); i < pk.size(); ++i) // j: first peak too heavy to pair with i
      {
        while (j > 0 && pk[i].second + pk[j - 1].second > target + tol)
        {
          --j;
        }
        for (OpenMS::Size q = j; q-- > 0 && pk[i].second + pk[q].second >= target - tol;)
        {
          if (q != i)
          {
            hit[i] = hit[q] = 1; // a peak is not its own complement
          }
        }
      }
      for (OpenMS::Size i = 0; i < pk.size(); ++i)
      {
        if (hit[i]) { pk[i].first *= opts.comp_boost; }
      }
    }
    if ((int)pk.size() > opts.max_frags)
    {
      std::nth_element(pk.begin(), pk.begin() + opts.max_frags, pk.end(), std::greater<>());
      pk.resize(opts.max_frags);
    }
    // by m/z, equal m/z by intensity: nth_element leaves the kept peaks in an unspecified order
    std::sort(pk.begin(), pk.end(), [](const auto& a, const auto& b) { return std::tie(a.second, a.first) < std::tie(b.second, b.first); });
    s.setMSLevel(2);
    s.setType(OpenMS::SpectrumSettings::SpectrumType::CENTROID);
    s.setRT(p.rt);
    OpenMS::Precursor pr;
    pr.setMZ(p.mz);
    pr.setCharge(p.z);
    pr.setIntensity(p.inten);
    pr.setIsolationWindowLowerOffset(p.mz - lo);
    pr.setIsolationWindowUpperOffset(hi - p.mz);
    s.setPrecursors({pr});
    if (opts.tag_src) { s.setMetaValue("precursor_source", p.src); }
    s.reserve(pk.size());
    for (const auto& q : pk)
    {
      OpenMS::Peak1D x;
      x.setMZ(q.second * ms2_scale);
      x.setIntensity(q.first);
      s.push_back(x);
    }
    return s;
  }
} // namespace Internal

/**
  @brief Extracts pseudo-MS2 spectra from a DIA run without ion mobility.

  @p exp: the resident run without ion mobility (MS1 and MS2 centroided, one isolation window per MS2 scan, in RT
  order). Its spectra are moved out: @p exp is left empty. Fills @p out with one MS2 spectrum per precursor.
  @param[in,out] exp The run; left empty.
  @param[in] opts Parameters.
  @param[out] out Receives the spectra, sorted by RT, native IDs "index=<i>".
  @param[in] log Receives the progress and summary lines.
  @param[in] progress Reports one step per finished isolation window (its log type decides where to).
  @return Counts and timings of the run.
  @par Exceptions
  Rethrows the first exception raised while processing a window, of whatever type it was.
*/
inline Stats run(OpenMS::PeakMap& exp,
                 const Params& opts,
                 OpenMS::PeakMap& out,
                 const std::function<void(const DIAspeXtract::String&)>& log,
                 const OpenMS::ProgressLogger& progress)
{
  using clk = std::chrono::steady_clock;
  const auto t0 = clk::now();
  Stats st;
  st.threads = omp_get_max_threads();

  // windows from the MS2 isolation windows; MS1 scans indexed once
  std::map<std::pair<double, double>, std::vector<OpenMS::Size>> win; // (lo, hi) -> MS2 scan indices
  std::vector<OpenMS::Size> ms1;
  for (OpenMS::Size i = 0; i < exp.size(); ++i)
  {
    const OpenMS::MSSpectrum& s = exp[i];
    if (s.getMSLevel() == 1) { ms1.push_back(i); }
    else if (s.getMSLevel() == 2 && ! s.getPrecursors().empty())
    {
      const OpenMS::Precursor& pr = s.getPrecursors()[0];
      const double lo = std::round((pr.getMZ() - pr.getIsolationWindowLowerOffset()) * 100) / 100;
      const double hi = std::round((pr.getMZ() + pr.getIsolationWindowUpperOffset()) * 100) / 100;
      win[{lo, hi}].push_back(i);
    }
  }
  std::vector<std::pair<double, double>> keys;
  std::vector<std::vector<OpenMS::Size>> scans;
  for (auto& kv : win)
  {
    keys.push_back(kv.first);
    scans.push_back(std::move(kv.second));
  }
  win.clear();
  const OpenMS::Size n_windows = keys.size();
  st.windows = n_windows;
  if (n_windows == 0 || ms1.empty())
  {
    log("[noims] no MS1 scans or no MS2 isolation windows");
    exp.clear(true);
    return st;
  }
  // a precursor belongs to the window containing it whose center is nearest (windows may overlap); a window with
  // fewer than 3 MS2 scans gives no spectra and so owns nothing
  auto owner = [&](double mz) {
    long best = -1;
    double bd = 1e9;
    for (OpenMS::Size w = 0; w < n_windows; ++w)
    {
      if (mz >= keys[w].first && mz < keys[w].second && scans[w].size() >= 3)
      {
        const double d = std::fabs(mz - 0.5 * (keys[w].first + keys[w].second));
        if (d < bd)
        {
          bd = d;
          best = (long)w;
        }
      }
    }
    return best;
  };
  // median spacing of these scans (fallback when fewer than two)
  auto spacing = [&](const std::vector<OpenMS::Size>& v, double fallback) {
    std::vector<double> d;
    for (OpenMS::Size k = 1; k < v.size(); ++k)
    {
      d.push_back(exp[v[k]].getRT() - exp[v[k - 1]].getRT());
    }
    if (d.empty()) { return fallback; }
    std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
    return d[d.size() / 2];
  };
  // MS2 cycle: the scans of the first window with >= 3 of them (a short window would give a wrong cycle to all others)
  const auto cyc_win = std::find_if(scans.begin(), scans.end(), [](const auto& v) { return v.size() >= 3; });
  const double cycle = spacing(cyc_win != scans.end() ? *cyc_win : scans.front(), 1.0);
  log("[noims] " + DIAspeXtract::String(n_windows) + " isolation windows, " + DIAspeXtract::String(ms1.size()) + " MS1 scans every "
      + DIAspeXtract::String(spacing(ms1, cycle)) + " s, MS2 cycle " + DIAspeXtract::String(cycle) + " s");

  // the run's spectra move into Biosaur2's MS1 map and the windows' MS2 maps (MassTraceDetection reads MS1 only)
  std::vector<OpenMS::PeakMap> m2(n_windows);
  std::vector<std::vector<double>> grid(n_windows);
  std::vector<OpenMS::Size> wpeaks(n_windows, 0);
  for (OpenMS::Size w = 0; w < n_windows; ++w)
  {
    for (OpenMS::Size i : scans[w])
    {
      grid[w].push_back(exp[i].getRT());
      wpeaks[w] += exp[i].size();
      exp[i].setMSLevel(1);
      m2[w].addSpectrum(std::move(exp[i]));
    }
  }
  OpenMS::PeakMap m1;
  for (OpenMS::Size i : ms1)
  {
    m1.addSpectrum(std::move(exp[i]));
  }
  exp.clear(true);

  // fragment m/z shift (P.ms2_shift ppm); auto estimates the run's offset from the y1 ions of tryptic peptides
  // (K 147.11280, R 175.11895) in the raw MS2 scans.
  // ponytail: the 1.7 ppm margin was tuned on one TripleTOF run; fit it per instrument if other TOFs disagree
  Internal::ms2_scale = 1.0;
  if (opts.ms2_shift_auto || opts.ms2_shift != 0)
  {
    double shift = opts.ms2_shift;
    if (opts.ms2_shift_auto)
    {
      std::vector<double> d, core;
      for (const auto& m : m2)
      {
        for (const auto& sp : m)
        {
          for (double a : {147.11280, 175.11895})
          {
            for (auto it = sp.MZBegin(a * (1 - 15e-6)); it != sp.end() && it->getMZ() <= a * (1 + 15e-6); ++it)
            {
              d.push_back((it->getMZ() - a) / a * 1e6);
            }
          }
        }
      }
      const auto median = [](std::vector<double>& v) {
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
      };
      const double m0 = d.size() < 1000 ? 0 : median(d);
      for (double x : d)
      {
        if (std::fabs(x - m0) < 6)
        {
          core.push_back(x); // second pass: without the flat random background
        }
      }
      const double off = d.size() < 1000 ? 0 : median(core);
      shift = d.size() < 1000 ? 0 : -(off + 1.7);
      log("[noims] fragment m/z offset from y1 ions: " + DIAspeXtract::String(off) + " ppm (n=" + DIAspeXtract::String(d.size()) + "), shift "
          + DIAspeXtract::String(shift) + " ppm");
    }
    Internal::ms2_scale = 1.0 + shift * 1e-6;
  }

  double top = 0;
  for (const auto& k : keys)
  {
    top = std::max(top, k.second);
  }
  const auto tb = clk::now();
  std::vector<Internal::Prec> seeds;
  const std::vector<Internal::Prec> bio
    = Internal::biosaur(std::move(m1), opts, keys.front().first - 1.0, top + 1.0, opts.hill_seeds > 0 ? &seeds : nullptr);
  log("[noims] Biosaur2: " + DIAspeXtract::String(bio.size()) + " isotope-pattern features as precursors in "
      + DIAspeXtract::String(std::chrono::duration<double>(clk::now() - tb).count()) + " s");

  std::vector<OpenMS::Size> order(n_windows); // heaviest window first, so the largest ones do not start last
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](OpenMS::Size a, OpenMS::Size b) { return wpeaks[a] > wpeaks[b]; });
  std::vector<std::vector<OpenMS::MSSpectrum>> res(n_windows);
  std::vector<OpenMS::Size> nprec(n_windows, 0), ntr2(n_windows, 0);
  std::vector<double> busy(n_windows, 0.0), t_mtd(n_windows, 0.0), t_frag(n_windows, 0.0);
  std::vector<OpenMS::Size> nfrag(n_windows, 0);
  std::exception_ptr err;
  OpenMS::Size n_done = 0;
  progress.startProgress(0, (OpenMS::SignedSize)n_windows, "extracting pseudo-spectra (isolation windows)");
#pragma omp parallel for schedule(dynamic, 1)
  for (long o = 0; o < (long)n_windows; ++o)
  {
    try
    {
      const OpenMS::Size w = order[o];
      if (scans[w].size() < 3)
      {
        continue; // MassTraceDetection needs 3 scans
      }
      const auto tw = clk::now();
      const double lo = keys[w].first, hi = keys[w].second;
      std::vector<Internal::Prec> prec;
      for (auto it = std::lower_bound(bio.begin(), bio.end(), lo, [](const Internal::Prec& p, double m) { return p.mz < m; });
           it != bio.end() && it->mz < hi; ++it)
      {
        if (owner(it->mz) == (long)w) { prec.push_back(*it); }
      }

      // MS2 traces on the window's scan grid
      const std::vector<OpenMS::MassTrace> ft = Internal::traces(m2[w], opts.ms2_ppm > 0 ? opts.ms2_ppm : opts.ppm);
      ntr2[w] = ft.size();
      t_mtd[w] = std::chrono::duration<double>(clk::now() - tw).count();
      m2[w].clear(true);
      // unfragmented precursors seen in MS2, unless MS1 already has them (same z, 1 cycle, ppm)
      std::vector<Internal::Prec> q = Internal::envelopes(ft, opts.ppm, lo, hi, 1.01 * cycle, opts.env_all_charges, opts.env_pattern);
      const OpenMS::Size n1 = prec.size();
      for (auto& c : q)
      {
        if (owner(c.mz) != (long)w) { continue; }
        bool dup = false;
        for (OpenMS::Size i = 0; i < n1 && ! dup; ++i)
        {
          dup = prec[i].z == c.z && std::fabs(prec[i].mz - c.mz) <= c.mz * opts.ppm * 1e-6 && std::fabs(prec[i].rt - c.rt) <= cycle;
        }
        if (! dup) { prec.push_back(std::move(c)); }
      }
      nprec[w] = prec.size();
      const std::vector<double>& gr = grid[w];
      std::vector<Internal::Frag> fr;
      fr.reserve(ft.size());
      std::vector<double> prof;
      for (const auto& t : ft)
      {
        if (t.getMaxIntensity(false) < opts.min_frag_intensity || (int)t.getSize() < opts.min_frag_points) { continue; }
        Internal::Frag f {t.getCentroidMZ(), Internal::apexRT(t), t.getMaxIntensity(false), {}, {}};
        for (const auto& pk : t)
        {
          f.pts.emplace_back(pk.getRT(), pk.getIntensity());
          f.gpts.emplace_back(std::lower_bound(gr.begin(), gr.end(), pk.getRT() - 1e-6) - gr.begin(), pk.getIntensity());
        }
        std::vector<OpenMS::Size> cut;
        if (opts.split > 0 && f.gpts.size() >= 3)
        {
          const OpenMS::Size g0 = f.gpts.front().first;
          prof.assign(f.gpts.back().first - g0 + 1, 0.0);
          for (const auto& q : f.gpts)
          {
            prof[q.first - g0] = q.second;
          }
          cut = Internal::valleys(prof, g0, opts.split);
        }
        if (cut.empty())
        {
          fr.push_back(std::move(f));
          continue;
        }
        // one fragment per elution peak: its own apex, intensity and intensity-weighted m/z
        cut.push_back(SIZE_MAX);
        OpenMS::Size k = 0;
        for (OpenMS::Size c : cut)
        {
          Internal::Frag g {0, 0, -1, {}, {}};
          double sw = 0;
          for (auto it = t.begin() + k; k < f.gpts.size() && f.gpts[k].first < c; ++k, ++it)
          {
            g.pts.push_back(f.pts[k]);
            g.gpts.push_back(f.gpts[k]);
            g.mz += it->getMZ() * it->getIntensity();
            sw += it->getIntensity();
            if (it->getIntensity() > g.inten)
            {
              g.inten = it->getIntensity();
              g.apex_rt = it->getRT();
            }
          }
          if (g.gpts.empty()) { continue; }
          g.mz = sw > 0 ? g.mz / sw : t.getCentroidMZ();
          if (g.inten >= opts.min_frag_intensity && (int)g.gpts.size() >= opts.min_frag_points) { fr.push_back(std::move(g)); }
        }
      }
      std::sort(fr.begin(), fr.end(), [](const Internal::Frag& a, const Internal::Frag& b) {
        return std::tie(a.apex_rt, a.mz, a.inten, a.pts) < std::tie(b.apex_rt, b.mz, b.inten, b.pts);
      });
      for (auto& f : fr)
      {
        std::vector<std::pair<double, double>>().swap(f.pts); // pts only break sort ties; gpts serve the rest
      }
      nfrag[w] = fr.size();
      t_frag[w] = std::chrono::duration<double>(clk::now() - tw).count() - t_mtd[w];

      // one slot per precursor, filled by tasks that idle threads can take
      std::vector<OpenMS::MSSpectrum> slot(prec.size());
      std::exception_ptr terr;
      const long np = (long)prec.size();
#pragma omp taskloop grainsize(256) default(shared)
      for (long i = 0; i < np; ++i)
      {
        try
        {
          slot[i] = Internal::assemble(prec[i], fr, gr, cycle, opts, lo, hi);
        }
        catch (...)
        {
          Internal::keepFirst(terr);
        }
      }
      if (terr) { std::rethrow_exception(terr); }
      for (auto& s : slot)
      {
        if (s.getMSLevel() == 2) { res[w].push_back(std::move(s)); }
      }
      if (opts.hill_seeds > 0)
      {
        // fragment-driven tier, after and apart from the spectra above: hills that are no isotope peak (k = 0..5) of a
        // feature, tried as z 2 and 3; the charge with more complementary pairs is kept if it has enough of them
        std::vector<Internal::Prec> sd; // per seed: z 2, then z 3
        const bool seed_mono = opts.seed_mono, share = ! opts.drop_precursor;
        for (auto it = std::lower_bound(seeds.begin(), seeds.end(), lo, [](const Internal::Prec& p, double m) { return p.mz < m; });
             it != seeds.end() && it->mz < hi; ++it)
        {
          if (owner(it->mz) != (long)w) { continue; }
          bool used = false;
          for (int fz = 1; fz <= 4 && ! used; ++fz)
          {
            for (int k = 0; k <= 5 && ! used; ++k)
            {
              const double t = it->mz - k * Internal::C13 / fz, tol = t * opts.ppm * 1e-6;
              for (auto b = std::lower_bound(bio.begin(), bio.end(), t - tol, [](const Internal::Prec& p, double m) { return p.mz < m; });
                   b != bio.end() && b->mz <= t + tol && ! used; ++b)
              {
                used = b->z == fz && it->rt_lo <= b->rt_hi && b->rt_lo <= it->rt_hi;
              }
            }
          }
          if (used) { continue; }
          for (int z = 2; z <= 3; ++z)
          {
            sd.push_back(*it);
            sd.back().z = z;
            if (! seed_mono) { continue; }
            // noims:seed_mono: a hill one isotope step below means this one is not the lightest at z
            const double t = it->mz - Internal::C13 / z, tol = t * opts.ppm * 1e-6;
            for (auto b = std::lower_bound(seeds.begin(), seeds.end(), t - tol, [](const Internal::Prec& p, double m) { return p.mz < m; });
                 b != seeds.end() && b->mz <= t + tol && sd.back().z; ++b)
            {
              if (it->rt_lo <= b->rt_hi && b->rt_lo <= it->rt_hi)
              {
                sd.back().z = 0; // 0 = skipped
              }
            }
          }
        }
        std::vector<OpenMS::MSSpectrum> s2(sd.size());
        std::vector<int> np2(sd.size(), 0);
        const long ns = (long)sd.size();
#pragma omp taskloop grainsize(128) default(shared)
        for (long i = 0; i < ns; i += 2)
        {
          try
          {
            std::vector<std::pair<double, double>> raw;
            const auto one = [&](long j, std::vector<std::pair<double, double>>* r) {
              s2[j] = Internal::assemble(sd[j], fr, gr, cycle, opts, lo, hi, r);
              if (s2[j].getMSLevel() == 2)
              {
                np2[j] = Internal::complementPairs(s2[j], sd[j].z * (sd[j].mz - Internal::PROTON), opts.ms2_ppm > 0 ? opts.ms2_ppm : opts.ppm,
                                                   (OpenMS::Size)opts.pair_top);
              }
            };
            const bool z3 = sd[i + 1].z > 0, sh = z3 && share;
            if (z3) { one(i + 1, share ? &raw : nullptr); }
            if (sd[i].z > 0 && ! (sh && raw.empty()))
            {
              one(i, sh ? &raw : nullptr); // nothing at z 3: nothing at z 2
            }
          }
          catch (...)
          {
            Internal::keepFirst(terr);
          }
        }
        if (terr) { std::rethrow_exception(terr); }
        for (long i = 0; i < ns; i += 2)
        {
          const long j = np2[i + 1] > np2[i] ? i + 1 : i; // one charge per seed; z 2 on a tie
          if (np2[j] < opts.hill_seeds) { continue; }
          if (opts.tag_src)
          {
            s2[j].setMetaValue("hill_pairs", np2[j]);
            s2[j].setMetaValue("hill_span", sd[j].rt_hi - sd[j].rt_lo);
          } // quality of the seed, for offline ranking
          res[w].push_back(std::move(s2[j]));
          ++nprec[w];
        }
      }
      busy[w] = std::chrono::duration<double>(clk::now() - tw).count();
#pragma omp critical(noims_progress)
      progress.setProgress((OpenMS::SignedSize)++n_done);
    }
    catch (...)
    {
      Internal::keepFirst(err);
    }
  }
  progress.endProgress();
  if (err) { std::rethrow_exception(err); }

  for (OpenMS::Size w = 0; w < n_windows; ++w)
  {
    st.precursors += nprec[w];
    st.busy_s += busy[w];
    st.ms2_traces += ntr2[w];
    for (auto& s : res[w])
    {
      st.peaks += s.size();
      out.addSpectrum(std::move(s));
    }
  }
  // by RT; spectra of one RT keep the assembly order (window key, then precursor order), which is platform-independent
  std::stable_sort(out.getSpectra().begin(), out.getSpectra().end(),
                   [](const OpenMS::MSSpectrum& a, const OpenMS::MSSpectrum& b) { return a.getRT() < b.getRT(); });
  for (OpenMS::Size i = 0; i < out.size(); ++i)
  {
    out[i].setNativeID("index=" + DIAspeXtract::String(i));
  }
  st.spectra = out.size();
  if (const auto n = std::count_if(scans.begin(), scans.end(), [](const auto& v) { return v.size() < 3; }))
  {
    log("[noims] " + DIAspeXtract::String(n) + " of " + DIAspeXtract::String(n_windows)
        + " isolation windows have fewer than 3 MS2 scans: no spectra from them");
  }
  double sm = 0, sf = 0;
  OpenMS::Size nf = 0;
  for (OpenMS::Size w = 0; w < n_windows; ++w)
  {
    sm += t_mtd[w];
    sf += t_frag[w];
    nf += nfrag[w];
  }
  log("[noims] " + DIAspeXtract::String(nf) + " fragments; window busy: MS2 traces " + DIAspeXtract::String(sm) + " s, envelopes and fragments "
      + DIAspeXtract::String(sf) + " s, spectra " + DIAspeXtract::String(st.busy_s - sm - sf) + " s");
  st.wall_s = std::chrono::duration<double>(clk::now() - t0).count();
  log("[noims] " + DIAspeXtract::String(st.ms2_traces) + " MS2 traces, " + DIAspeXtract::String(st.precursors) + " precursors, "
      + DIAspeXtract::String(st.spectra) + " spectra, " + DIAspeXtract::String(st.peaks) + " peaks; " + DIAspeXtract::String(st.wall_s)
      + " s, window busy " + DIAspeXtract::String(st.busy_s) + " s on " + DIAspeXtract::String(st.threads) + " threads");
  return st;
}
} // namespace DIAspeXtract::NoIms
