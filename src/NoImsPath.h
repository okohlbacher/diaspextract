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
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iterator>
#include <limits>
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
  int env_ms1_support = 0;     ///< MS2 envelopes need their mono m/z in MS1 (15 ppm) in at least this many scans within 1.5 cycles of their
                               ///< apex; 0 = off. Without MS1 signal an envelope is mostly a fragment ion, not an unfragmented precursor.
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
  /// @brief A window's fragments in apex-RT order as flat arrays (one allocation per array, points contiguous in fragment order).
  struct Frags
  {
    std::vector<double> apex, mz, inten; ///< per fragment: apex RT, m/z, apex intensity
    std::vector<uint32_t> off;           ///< per fragment: first point; off[n] = number of points
    std::vector<uint32_t> g;             ///< per point: window grid index
    std::vector<double> in;              ///< per point: intensity
    std::vector<uint32_t> multi;         ///< fragments with more than one point, in apex order
    std::vector<uint32_t> s_off, s_idx;  ///< single-point fragments by grid index (s_off[g] .. s_off[g + 1]), each by m/z
    OpenMS::Size size() const { return apex.size(); }
    Frags() = default;
    void swap(Frags& o)
    {
      apex.swap(o.apex); mz.swap(o.mz); inten.swap(o.inten); off.swap(o.off); g.swap(o.g); in.swap(o.in);
      multi.swap(o.multi); s_off.swap(o.s_off); s_idx.swap(o.s_idx);
    }
    /// @brief Flattens @p fr (already in apex order) and frees it.
    explicit Frags(std::vector<Frag>& fr)
    {
      const OpenMS::Size n = fr.size();
      OpenMS::Size np = 0;
      for (const auto& f : fr) { np += f.gpts.size(); }
      if (np >= std::numeric_limits<uint32_t>::max()) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); } // 32-bit offsets
      apex.resize(n); mz.resize(n); inten.resize(n); off.resize(n + 1); g.resize(np); in.resize(np);
      OpenMS::Size j = 0;
      for (OpenMS::Size k = 0; k < n; ++k)
      {
        apex[k] = fr[k].apex_rt; mz[k] = fr[k].mz; inten[k] = fr[k].inten; off[k] = (uint32_t)j;
        for (const auto& q : fr[k].gpts) { g[j] = (uint32_t)q.first; in[j] = q.second; ++j; }
      }
      off[n] = (uint32_t)j;
      std::vector<Frag>().swap(fr);
      // single-point fragments per grid index, by m/z (assemble's fast path); the others in apex order
      uint32_t ng = 0;
      for (uint32_t x : g) { ng = std::max(ng, x + 1); }
      s_off.assign(ng + 1, 0);
      for (OpenMS::Size k = 0; k < n; ++k)
      {
        if (off[k + 1] - off[k] == 1) { ++s_off[g[off[k]] + 1]; }
        else { multi.push_back((uint32_t)k); }
      }
      for (uint32_t x = 0; x < ng; ++x) { s_off[x + 1] += s_off[x]; }
      s_idx.resize(s_off[ng]);
      std::vector<uint32_t> fill(s_off.begin(), s_off.end() - 1);
      for (OpenMS::Size k = 0; k < n; ++k)
      {
        if (off[k + 1] - off[k] == 1) { s_idx[fill[g[off[k]]]++] = (uint32_t)k; }
      }
      for (uint32_t x = 0; x < ng; ++x)
      {
        std::sort(s_idx.begin() + s_off[x], s_idx.begin() + s_off[x + 1], [&](uint32_t a, uint32_t b) { return mz[a] < mz[b] || (mz[a] == mz[b] && a < b); });
      }
    }
  };

  /**
    @brief A mass trace as the non-IMS path reads it: its points (RT order) in a shared array, centroid m/z, apex intensity.

    The subset of OpenMS::MassTrace this path uses, with the same values (FWHM is never estimated by
    MassTraceDetection: 0).
  */
  struct LTrace
  {
    const OpenMS::Peak2D* b; ///< first point
    OpenMS::Size n;          ///< number of points
    double cmz, maxi;        ///< centroid m/z (MassTrace::updateWeightedMeanMZ), largest raw intensity (getMaxIntensity(false))
    const OpenMS::Peak2D* begin() const { return b; }
    const OpenMS::Peak2D* end() const { return b + n; }
    std::reverse_iterator<const OpenMS::Peak2D*> rbegin() const { return std::reverse_iterator<const OpenMS::Peak2D*>(b + n); }
    OpenMS::Size getSize() const { return n; }
    double getCentroidMZ() const { return cmz; }
    double getFWHM() const { return 0.0; }
    double getMaxIntensity(bool) const { return maxi; }
  };
  /// @brief A window's mass traces: LTraces over one point array.
  struct Traces
  {
    std::vector<OpenMS::Peak2D> pts;
    std::vector<LTrace> t;
    OpenMS::Size size() const { return t.size(); }
    const LTrace& operator[](OpenMS::Size i) const { return t[i]; }
    std::vector<LTrace>::const_iterator begin() const { return t.begin(); }
    std::vector<LTrace>::const_iterator end() const { return t.end(); }
    /// @brief Appends trace bounds [p0, pts.size()) with its centroid m/z; pointers are fixed by finish().
    void add(OpenMS::Size p0, double cmz)
    {
      double mx = 0.0;
      for (OpenMS::Size i = p0; i < pts.size(); ++i)
      {
        if (pts[i].getIntensity() > mx) { mx = pts[i].getIntensity(); }
      }
      t.push_back(LTrace {reinterpret_cast<const OpenMS::Peak2D*>(p0), pts.size() - p0, cmz, mx});
    }
    void finish()
    {
      for (auto& x : t) { x.b = pts.data() + reinterpret_cast<OpenMS::Size>(x.b); }
    }
  };

  /**
    @brief The RT of a trace's most intense point.
    @param[in] t The mass trace.
    @return The apex RT (the first of equal maxima).
  */
  inline double apexRT(const LTrace& t)
  {
    double best = -1, rt = 0;
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
  inline Traces tracesOpenMS(const OpenMS::PeakMap& m, double ppm)
  {
    OpenMS::MassTraceDetection mtd;
    mtd.setLogType(OpenMS::ProgressLogger::NONE); // windows run in parallel: progress lines would interleave
    OpenMS::Param p = mtd.getDefaults();
    p.setValue("mass_error_ppm", ppm);
    p.setValue("min_trace_length", 0.0); // single points too: a fragment seen in one scan of a window
    mtd.setParameters(p);
    std::vector<OpenMS::MassTrace> mt;
    mtd.run(m, mt);
    Traces out;
    for (const auto& x : mt)
    {
      const OpenMS::Size p0 = out.pts.size();
      for (const auto& q : x) { out.pts.push_back(q); }
      out.add(p0, x.getCentroidMZ());
    }
    out.finish();
    return out;
  }

  /**
    @brief Mass traces of a map, the same as OpenMS MassTraceDetection with its defaults, mass_error_ppm = @p ppm and
    min_trace_length = 0 (tracesOpenMS), without its per-trace allocations; frees the map's spectra as it reads them.

    Falls back to tracesOpenMS when a spectrum has data arrays (ion mobility or FWHM arrays change what MassTraceDetection
    does; OpenMS also validates integer and string arrays against the peak count).
    @param[in,out] m The spectra to trace (all of MS level 1, sorted by m/z); left empty.
    @param[in] ppm Mass error tolerance (ppm).
    @return The traces in MassTraceDetection's order, with its centroid m/z.
  */
  inline Traces traces(OpenMS::PeakMap& m, double ppm)
  {
    for (const auto& s : m)
    {
      if (! s.getFloatDataArrays().empty() || ! s.getIntegerDataArrays().empty() || ! s.getStringDataArrays().empty() || s.getMSLevel() != 1)
      {
        Traces t = tracesOpenMS(m, ppm);
        m.clear(true);
        return t;
      }
    }
    constexpr double NOISE = 10.0, SNR = 3.0; // MassTraceDetection defaults noise_threshold_int, chrom_peak_snr
    constexpr OpenMS::Size OUTLIERS = 5;      // trace_termination_outliers (criterion "outlier")
    constexpr double MIN_SAMPLE_RATE = 0.5;
    const OpenMS::Size ns = m.size();
    if (ns < 3)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input map consists of too few MS1 spectra (less than 3!). Aborting...",
                                            DIAspeXtract::String(ns));
    }
    // step 1: peaks above the noise threshold (flat), apex candidates above SNR x noise
    // one 16-byte record per peak: the binary search's last probe, the intensity and the visited flag share a cache line
    struct PK
    {
      double mz;
      float in;
      uint32_t vis;
    };
    std::vector<double> rt(ns);
    std::vector<PK> pk;
    std::vector<OpenMS::Size> off(ns + 1, 0);
    struct Apex
    {
      double intensity;
      uint32_t scan, peak;
    };
    std::vector<Apex> ap;
    for (OpenMS::Size i = 0; i < ns; ++i)
    {
      OpenMS::MSSpectrum& s = m[i];
      rt[i] = s.getRT();
      uint32_t k = 0;
      for (const auto& q : s)
      {
        const double v = q.getIntensity();
        if (v > NOISE)
        {
          if (v > SNR * NOISE) { ap.push_back(Apex {v, (uint32_t)i, k}); }
          pk.push_back(PK {q.getMZ(), q.getIntensity(), 0});
          ++k;
        }
      }
      off[i + 1] = pk.size();
      s.clear(true); // frees the peaks (and the meta data)
      s.shrink_to_fit();
    }
    m.clear(true);
    // stable sort by intensity: an LSD radix sort on the float bits (positive floats order as their bit patterns; the
    // intensities are floats widened to double), stable like std::stable_sort
    {
      std::vector<Apex> tmp(ap.size());
      for (int sh = 0; sh < 32; sh += 8)
      {
        OpenMS::Size cnt[257] = {0};
        for (const auto& x : ap) { ++cnt[((std::bit_cast<uint32_t>((float)x.intensity) >> sh) & 255) + 1]; }
        for (int d = 0; d < 256; ++d) { cnt[d + 1] += cnt[d]; }
        for (const auto& x : ap) { tmp[cnt[(std::bit_cast<uint32_t>((float)x.intensity) >> sh) & 255]++] = x; }
        ap.swap(tmp);
      }
    }
    Traces out;
    std::vector<std::pair<OpenMS::Size, OpenMS::Size>> down, up; // (scan, flat peak index) gathered below and above the apex
    const PK* const PM = pk.data();
    // per scan, an m/z bin index (about 4 peaks per bin) that narrows the binary search of nearest(): bin(x) =
    // floor((x - lo) * inv), clamped, is monotonic in x, so the first peak >= mz lies in bin(mz) or is the first of a later bin
    std::vector<double> blo(ns), binv(ns);
    std::vector<uint32_t> bn(ns), boff(ns + 1, 0), bfirst;
    if (pk.size() / 4 + 2 * ns >= std::numeric_limits<uint32_t>::max()) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); } // 32-bit bin offsets
    for (OpenMS::Size i = 0; i < ns; ++i)
    {
      const OpenMS::Size np = off[i + 1] - off[i];
      const uint32_t nb = (uint32_t)std::max<OpenMS::Size>(1, np / 4);
      bn[i] = nb;
      boff[i + 1] = boff[i] + nb + 1;
      blo[i] = np ? PM[off[i]].mz : 0.0;
      const double span = np ? PM[off[i + 1] - 1].mz - blo[i] : 0.0;
      binv[i] = span > 0 ? nb / span : 0.0;
    }
    bfirst.resize(boff[ns]);
    auto binOf = [&](OpenMS::Size i, double mz) -> uint32_t {
      const double t = (mz - blo[i]) * binv[i];
      if (! (t > 0)) { return 0; }
      return t >= bn[i] ? bn[i] - 1 : (uint32_t)t;
    };
    for (OpenMS::Size i = 0; i < ns; ++i)
    {
      uint32_t* f = bfirst.data() + boff[i];
      uint32_t b = 0;
      for (OpenMS::Size j = off[i]; j < off[i + 1]; ++j)
      {
        const uint32_t pb = binOf(i, PM[j].mz);
        while (b <= pb) { f[b++] = (uint32_t)(j - off[i]); }
      }
      while (b <= bn[i]) { f[b++] = (uint32_t)(off[i + 1] - off[i]); }
    }
    // MSSpectrum::findNearest on scan i's peaks: the first peak at or above mz, or the one before it when that is not farther
    auto nearest = [&](OpenMS::Size i, double mz) -> OpenMS::Size {
      const PK* b = PM + off[i];
      const PK* e = PM + off[i + 1];
      const uint32_t* f = bfirst.data() + boff[i];
      const uint32_t q = binOf(i, mz);
      const PK* it = std::lower_bound(b + f[q], b + f[q + 1], mz, [](const PK& p, double x) { return p.mz < x; });
      if (it == b) { return off[i]; }
      if (it == e) { return off[i + 1] - 1; }
      return (std::fabs(it->mz - mz) < std::fabs((it - 1)->mz - mz)) ? (OpenMS::Size)(it - PM) : (OpenMS::Size)(it - 1 - PM);
    };
    for (auto a = ap.crbegin(); a != ap.crend(); ++a)
    {
      const OpenMS::Size apx = off[a->scan] + a->peak;
      if (pk[apx].vis) { continue; }
      const double amz = pk[apx].mz;
      const float ain = pk[apx].in;
      double centroid = amz, counter = (double)ain * amz, denom = ain;
      auto wmean = [&](double val, double w) { // MassTraceDetection::updateIterativeWeightedMean_
        const double ct = 1.0 + (w * val) / counter;
        const double dt = 1.0 + w / denom;
        centroid *= (ct / dt);
        counter *= ct;
        denom *= dt;
      };
      wmean(amz, ain);
      double sd = (centroid / 1e6) * ppm, wsum = ain;
      down.clear();
      up.clear();
      OpenMS::Size di = a->scan, ui = a->scan, d_scans = 0, u_scans = 0, d_miss = 0, u_miss = 0;
      bool d_on = true, u_on = true;
      auto step = [&](OpenMS::Size i, std::vector<std::pair<OpenMS::Size, OpenMS::Size>>& got, OpenMS::Size& miss) {
        if (off[i] == off[i + 1]) { return; } // empty spectra leave the counters alone
        const OpenMS::Size c = nearest(i, centroid);
        const double cmz = pk[c].mz;
        if (cmz <= centroid + 3 * sd && cmz >= centroid - 3 * sd && ! pk[c].vis)
        {
          const double cin = pk[c].in;
          got.emplace_back(i, c);
          wmean(cmz, cin);
          // updateWeightedSDEstimateRobust with the point as a Peak2D (float intensity)
          const float fin = (float)cin;
          const double d1 = pmath::log(wsum) + 2 * pmath::log(sd);
          const double d2 = pmath::log((double)fin) + 2 * pmath::log(std::abs(cmz - centroid));
          const double dd = std::sqrt(pmath::exp(d1) + pmath::exp(d2));
          const double ws = wsum + fin;
          const double tsd = dd / std::sqrt(ws);
          if (tsd > std::numeric_limits<double>::epsilon()) { sd = tsd; }
          wsum = ws;
          miss = 0;
        }
        else { ++miss; }
      };
      // the next scans' bin entries and peaks are cache misses: fetch them one and two steps ahead (no effect on results)
      auto ahead = [&](OpenMS::Size i1, OpenMS::Size i2, bool has2) {
        __builtin_prefetch(PM + off[i1] + bfirst[boff[i1] + binOf(i1, centroid)]);
        if (has2) { __builtin_prefetch(bfirst.data() + boff[i2] + binOf(i2, centroid)); }
      };
      while ((di > 0 && d_on) || (ui < ns - 1 && u_on))
      {
        if (di > 1 && d_on) { ahead(di - 2, di > 2 ? di - 3 : 0, di > 2); }
        if (ui + 2 < ns && u_on) { ahead(ui + 2, ui + 3 < ns ? ui + 3 : 0, ui + 3 < ns); }
        if (di > 0 && d_on)
        {
          step(di - 1, down, d_miss);
          --di;
          ++d_scans;
          if (d_miss > OUTLIERS) { d_on = false; }
        }
        if (ui < ns - 1 && u_on)
        {
          step(ui + 1, up, u_miss);
          ++ui;
          ++u_scans;
          if (u_miss > OUTLIERS) { u_on = false; }
        }
      }
      const OpenMS::Size npts = down.size() + 1 + up.size();
      const OpenMS::Size adjusted = d_scans + u_scans + 1 - d_miss - u_miss;
      if (! (static_cast<double>(npts) / static_cast<double>(adjusted) >= MIN_SAMPLE_RATE)) { continue; }
      // (min_trace_length 0 and no maximum: every RT range passes)
      pk[apx].vis = 1;
      for (const auto& c : down) { pk[c.second].vis = 1; }
      for (const auto& c : up) { pk[c.second].vis = 1; }
      const OpenMS::Size p0 = out.pts.size();
      auto put = [&](OpenMS::Size c, double r) {
        OpenMS::Peak2D q;
        q.setRT(r);
        q.setMZ(pk[c].mz);
        q.setIntensity(pk[c].in);
        out.pts.push_back(q);
      };
      for (auto it = down.rbegin(); it != down.rend(); ++it) { put(it->second, rt[it->first]); }
      put(apx, rt[a->scan]);
      for (const auto& c : up) { put(c.second, rt[c.first]); }
      double cmz;
      if (npts == 1) { cmz = out.pts[p0].getMZ(); }
      else
      { // MassTrace::updateWeightedMeanMZ
        double wsm = 0.0, tw = 0.0;
        for (OpenMS::Size i = p0; i < out.pts.size(); ++i)
        {
          const double w = out.pts[i].getIntensity();
          tw += w;
          wsm += w * out.pts[i].getMZ();
        }
        if (tw < std::numeric_limits<double>::epsilon())
        {
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "All weights were equal to zero! Empty trace? Aborting...",
                                                DIAspeXtract::String(tw));
        }
        cmz = wsm / tw;
      }
      out.add(p0, cmz);
    }
    out.finish();
    return out;
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
  inline std::vector<Prec> envelopes(const Traces& mt,
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
      const LTrace& a = mt[o[k]];
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
      const LTrace& t = mt[o[k]];
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
                                     const Frags& fr,
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
    thread_local std::vector<double> pp;
    pp.clear();
    for (OpenMS::Size g = g0; g < g1; ++g)
    {
      pp.push_back(interp(p.rt_pts, p.int_pts, grid[g]));
    }
    const bool few = pp.size() < 3;
    const double neutral = p.z * (p.mz - PROTON);
    thread_local std::vector<std::pair<double, double>> pk; // (intensity, m/z)
    pk.clear();
    thread_local std::vector<double> fp;
    fp.resize(pp.size());
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
      const double* const fmz = fr.mz.data();
      const double* const fap = fr.apex.data();
      const uint32_t* const foff = fr.off.data();
      const uint32_t* const fg = fr.g.data();
      const double* const fin = fr.in.data();
      const double mz_max = neutral - 56.0;
      const OpenMS::Size nfr = fr.size();
      // TOF fast path (split traces, padded span): single-point fragments by scan, where r depends on the scan only
      // (b cancels); a scan whose r is below min_corr by a margin far above rounding skips all its single-point fragments
      const bool fast = sparse && pad && ! few && opts.min_corr > 0 && ! opts.drop_precursor;
      if (fast && va > 0)
      {
        const double thr = opts.min_corr * std::sqrt(va * (nn - 1)) * (1 - 1e-9);
        const OpenMS::Size ge = std::min<OpenMS::Size>(g1, fr.s_off.empty() ? 0 : fr.s_off.size() - 1);
        for (OpenMS::Size g = g0; g < ge; ++g)
        {
          const double y = pp[g - g0], t = nn * y - sa;
          if (t + 1e-12 * (nn * std::fabs(y) + std::fabs(sa)) < thr) { continue; }
          for (uint32_t x = fr.s_off[g], xe = fr.s_off[g + 1]; x < xe; ++x)
          {
            const uint32_t k = fr.s_idx[x];
            const double f_mz = fmz[k];
            if (f_mz - PROTON > mz_max) { break; } // by m/z: all the rest are heavier
            const double b = fin[foff[k]];
            const double sb = b, sbb = b * b, sab = pp[g - g0] * b; // the general loop's sums for one point
            const double vb = nn * sbb - sb * sb;
            const double num = nn * sab - sa * sb;
            if (! (vb > 0) || ! (num > 0)) { continue; }
            if (num * num * (1 + 1e-9) < opts.min_corr * opts.min_corr * (va * vb)) { continue; }
            const double r = num / std::sqrt(va * vb);
            if (r < opts.min_corr) { continue; }
            pk.emplace_back(fr.inten[k] * pmath::pow(r, opts.weight_pow), f_mz);
          }
        }
      }
      const uint32_t* const fmulti = fr.multi.data();
      const OpenMS::Size nmulti = fr.multi.size();
      for (OpenMS::Size x = fast ? std::lower_bound(fr.multi.begin(), fr.multi.end(), a0, [&](uint32_t k, double t) { return fap[k] < t; }) - fr.multi.begin()
                                 : std::lower_bound(fr.apex.begin(), fr.apex.end(), a0) - fr.apex.begin();
           fast ? x < nmulti && fap[fmulti[x]] <= a1 : x < nfr && fap[x] <= a1; ++x)
      {
        const OpenMS::Size k = fast ? fmulti[x] : x;
        const double f_mz = fmz[k];
        if (f_mz - PROTON > mz_max)
        {
          continue; // heavier than the precursor allows at z = 1
        }
        if (opts.drop_precursor)
        {
          const double iso = std::round((f_mz - p.mz) * p.z / C13);
          if (iso >= -1 && iso <= 3 && std::fabs(f_mz - p.mz - iso * C13 / p.z) <= f_mz * opts.ppm * 1e-6)
          {
            continue; // unfragmented precursor
          }
        }
        if (! few)
        {
          double r;
          uint32_t j = foff[k];
          const uint32_t je = foff[k + 1];
          if (sparse)
          {
            double sb = 0, sbb = 0, sab = 0;
            while (j < je && fg[j] < g0) { ++j; }
            for (; j < je && fg[j] < g1; ++j)
            {
              const double b = fin[j];
              sb += b;
              sbb += b * b;
              sab += pp[fg[j] - g0] * b;
            }
            const double vb = nn * sbb - sb * sb;
            const double num = nn * sab - sa * sb;
            if (opts.min_corr > 0) // fragments that surely fail min_corr skip the square root and division
            {
              if (! (va > 0 && vb > 0) || ! (num > 0)) { continue; } // r <= 0 < min_corr
              // r^2 below min_corr^2 by a relative margin far above the rounding of r: r < min_corr as computed below too
              if (num * num * (1 + 1e-9) < opts.min_corr * opts.min_corr * (va * vb)) { continue; }
            }
            r = (va > 0 && vb > 0) ? num / std::sqrt(va * vb) : 0.0;
          }
          else
          {
            std::fill(fp.begin(), fp.end(), 0.0);
            for (; j < je; ++j)
            {
              if (fg[j] >= g0 && fg[j] < g1) { fp[fg[j] - g0] = fin[j]; }
            }
            r = pearson(pp, fp);
          }
          if (r < opts.min_corr) { continue; }
          pk.emplace_back(fr.inten[k] * pmath::pow(r, opts.weight_pow), f_mz);
          continue;
        }
        const double d = (fap[k] - p.rt) / cycle;
        pk.emplace_back(fr.inten[k] * pmath::exp(-d * d * opts.weight_pow / 4.0), f_mz);
      }
      if (raw) { *raw = pk; }
    }
    if ((int)pk.size() < opts.min_frags) { return s; }
    if (opts.comp_boost > 0)
    {
      // a peak with a complement (another peak q, m/z sum within target +- tol) gets the boost. The peaks go into 1-Th
      // buckets (lists in a reused head table), so the search needs no sort; what follows does not depend on pk's order
      const double target = neutral + 2 * PROTON, tol = target * 2 * opts.ppm * 1e-6;
      constexpr int NB = 1 << 14; // buckets of 1 Th: m/z 0 .. 16384
      thread_local std::vector<int> head(NB, -1);
      thread_local std::vector<int> nxt;
      thread_local std::vector<char> hit;
      const OpenMS::Size n = pk.size();
      nxt.resize(n);
      hit.assign(n, 0);
      bool inside = true;
      for (OpenMS::Size i = 0; i < n && inside; ++i) { inside = pk[i].second >= 0 && pk[i].second < NB; }
      if (inside)
      {
        for (OpenMS::Size i = 0; i < n; ++i)
        {
          const int b = (int)pk[i].second;
          nxt[i] = head[b];
          head[b] = (int)i;
        }
        for (OpenMS::Size i = 0; i < n; ++i)
        {
          const double a = pk[i].second;
          // clamp in double before the cast: an absurd tolerance must not overflow int
          const int b0 = (int)std::max(0.0, std::floor(target - tol - a - 1e-6)), b1 = (int)std::min((double)(NB - 1), std::floor(target + tol - a + 1e-6));
          for (int b = b0; b <= b1 && ! hit[i]; ++b)
          {
            for (int q = head[b]; q >= 0; q = nxt[q])
            {
              const double sum = a + pk[q].second;
              if ((OpenMS::Size)q != i && sum >= target - tol && ! (sum > target + tol))
              {
                hit[i] = 1; // a peak is not its own complement
                break;
              }
            }
          }
        }
        for (OpenMS::Size i = 0; i < n; ++i) { head[(int)pk[i].second] = -1; }
      }
      else // m/z outside the bucket range: the sorted two-pointer search
      {
        std::sort(pk.begin(), pk.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
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
      }
      for (OpenMS::Size i = 0; i < n; ++i)
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
  // MS1 scans as RT and sorted m/z, kept for the MS2 envelopes' MS1 support (Biosaur2 consumes m1)
  std::vector<double> ms1_rt;
  std::vector<std::vector<float>> ms1_mz;
  if (opts.env_ms1_support > 0)
  {
    for (const auto& s : m1)
    {
      ms1_rt.push_back(s.getRT());
      std::vector<float> v;
      v.reserve(s.size());
      for (const auto& p : s)
      {
        v.push_back((float)p.getMZ());
      }
      std::sort(v.begin(), v.end());
      ms1_mz.push_back(std::move(v));
    }
  }
  // number of MS1 scans within 1.5 cycles of rt with a peak within 15 ppm of mz
  const auto ms1Support = [&](double mz, double rt) {
    int n = 0;
    const double tol = mz * 15e-6;
    for (auto it = std::lower_bound(ms1_rt.begin(), ms1_rt.end(), rt - 1.5 * cycle); it != ms1_rt.end() && *it <= rt + 1.5 * cycle; ++it)
    {
      const auto& v = ms1_mz[it - ms1_rt.begin()];
      const auto p = std::lower_bound(v.begin(), v.end(), (float)(mz - tol));
      n += p != v.end() && *p <= mz + tol;
    }
    return n;
  };
  // Biosaur2 (MS1) and the windows' MS2 traces and fragments do not depend on each other: Biosaur2 runs as one task
  // (its parallel loops then run serially, as a nested region; its output does not depend on the thread count) while
  // the other threads trace the windows; the windows' precursors and spectra follow once both are done.
  std::vector<OpenMS::Size> order(n_windows); // heaviest window first, so the largest ones do not start last
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&](OpenMS::Size a, OpenMS::Size b) { return wpeaks[a] > wpeaks[b]; });
  std::vector<std::vector<OpenMS::MSSpectrum>> res(n_windows);
  std::vector<OpenMS::Size> nprec(n_windows, 0), ntr2(n_windows, 0);
  std::vector<double> busy(n_windows, 0.0), t_mtd(n_windows, 0.0), t_frag(n_windows, 0.0);
  std::vector<OpenMS::Size> nfrag(n_windows, 0), nsing(n_windows, 0);
  std::vector<Internal::Frags> wfr(n_windows);           // per window: fragments (phase A -> B)
  std::vector<std::vector<Internal::Prec>> wenv(n_windows); // per window: owned MS2 envelopes with MS1 support (phase A -> B)
  std::exception_ptr err;
  OpenMS::Size n_done = 0;
  std::vector<Internal::Prec> seeds, bio_v;
  double t_bio = 0;
  OpenMS::Size p1 = 0, p2 = 0; // MS1 and MS2 peaks
  for (const auto& sp : m1) { p1 += sp.size(); }
  for (OpenMS::Size w = 0; w < n_windows; ++w) { p2 += wpeaks[w]; }
  // overlap Biosaur2 with the MS2 traces only when there is enough MS2 work to hide it behind (serial, it takes about twice
  // its parallel time); NOIMS_OVERLAP=0/1 forces it (experiments)
  const char* ov_env = std::getenv("NOIMS_OVERLAP");
  const bool overlap = ov_env ? std::atoi(ov_env) != 0 : (double)p2 >= 15.0 * (double)p1 * std::max(1.0, st.threads / 16.0);
  log("[noims] " + DIAspeXtract::String(p1) + " MS1 peaks, " + DIAspeXtract::String(p2) + " MS2 peaks; Biosaur2 " + (overlap ? "alongside" : "before")
      + " the MS2 traces");
  if (! overlap)
  {
    const auto tb = clk::now();
    bio_v = Internal::biosaur(std::move(m1), opts, keys.front().first - 1.0, top + 1.0, opts.hill_seeds > 0 ? &seeds : nullptr);
    t_bio = std::chrono::duration<double>(clk::now() - tb).count();
  }
  const int max_levels = omp_get_max_active_levels();
  omp_set_max_active_levels(1); // Biosaur2's loops inside its task stay serial
  progress.startProgress(0, (OpenMS::SignedSize)n_windows, "extracting pseudo-spectra (isolation windows)");
  // phase B of one window (precursors, spectra, hill seeds): reads only that window's phase-A results and Biosaur2's output
  auto phase_b = [&](long o) {
        try
        {
          const OpenMS::Size w = order[o];
          if (scans[w].size() >= 3)
          {
            const auto tw = clk::now();
            const std::vector<Internal::Prec>& bio = bio_v;
            const double lo = keys[w].first, hi = keys[w].second;
            std::vector<Internal::Prec> prec;
            for (auto it = std::lower_bound(bio.begin(), bio.end(), lo, [](const Internal::Prec& p, double m) { return p.mz < m; });
                 it != bio.end() && it->mz < hi; ++it)
            {
              if (owner(it->mz) == (long)w) { prec.push_back(*it); }
            }
            // MS2 envelopes unless MS1 already has them (same z, 1 cycle, ppm)
            const OpenMS::Size n1 = prec.size();
            for (auto& c : wenv[w])
            {
              bool dup = false;
              for (OpenMS::Size i = 0; i < n1 && ! dup; ++i)
              {
                dup = prec[i].z == c.z && std::fabs(prec[i].mz - c.mz) <= c.mz * opts.ppm * 1e-6 && std::fabs(prec[i].rt - c.rt) <= cycle;
              }
              if (! dup && c.src == 1) { prec.push_back(std::move(c)); }
            }
            std::vector<Internal::Prec>().swap(wenv[w]);
            nprec[w] = prec.size();
            const std::vector<double>& gr = grid[w];
            const Internal::Frags& frs = wfr[w];
      // one slot per precursor, filled by tasks that idle threads can take
      std::vector<OpenMS::MSSpectrum> slot(prec.size());
      std::exception_ptr terr;
      const long np = (long)prec.size();
#pragma omp taskloop grainsize(256) default(shared)
      for (long i = 0; i < np; ++i)
      {
        try
        {
          slot[i] = Internal::assemble(prec[i], frs, gr, cycle, opts, lo, hi);
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
              s2[j] = Internal::assemble(sd[j], frs, gr, cycle, opts, lo, hi, r);
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
            Internal::Frags().swap(wfr[w]);
            busy[w] += std::chrono::duration<double>(clk::now() - tw).count();
          }
#pragma omp critical(noims_progress)
          progress.setProgress((OpenMS::SignedSize)++n_done);
        }
        catch (...)
        {
          Internal::keepFirst(err);
        }
  };
#pragma omp parallel
#pragma omp single
  {
#pragma omp task default(shared) if (overlap)
    if (overlap)
    {
      try
      {
        const auto tb = clk::now();
        bio_v = Internal::biosaur(std::move(m1), opts, keys.front().first - 1.0, top + 1.0, opts.hill_seeds > 0 ? &seeds : nullptr);
        t_bio = std::chrono::duration<double>(clk::now() - tb).count();
      }
      catch (...)
      {
        Internal::keepFirst(err);
      }
    }
    for (long o = 0; o < (long)n_windows; ++o)
    {
#pragma omp task default(shared) firstprivate(o)
      {
        try
        {
          const OpenMS::Size w = order[o];
          if (scans[w].size() >= 3) // MassTraceDetection needs 3 scans
          {
            const auto tw = clk::now();
            const double lo = keys[w].first, hi = keys[w].second;
      // MS2 traces on the window's scan grid
      const Internal::Traces ft = Internal::traces(m2[w], opts.ms2_ppm > 0 ? opts.ms2_ppm : opts.ppm);
      ntr2[w] = ft.size();
      t_mtd[w] = std::chrono::duration<double>(clk::now() - tw).count();
      m2[w].clear(true);
      // unfragmented precursors seen in MS2 (the MS1 duplicates are dropped once Biosaur2 is done)
      for (auto& c : Internal::envelopes(ft, opts.ppm, lo, hi, 1.01 * cycle, opts.env_all_charges, opts.env_pattern))
      {
        if (owner(c.mz) != (long)w) { continue; }
        c.src = opts.env_ms1_support > 0 && ms1Support(c.mz, c.rt) < opts.env_ms1_support ? -1 : 1; // -1: kept only as a duplicate
        wenv[w].push_back(std::move(c));
      }
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
      wfr[w] = Internal::Frags(fr);
      nsing[w] = wfr[w].s_idx.size();
      t_frag[w] = std::chrono::duration<double>(clk::now() - tw).count() - t_mtd[w];

            busy[w] = std::chrono::duration<double>(clk::now() - tw).count();
          }
          if (! overlap) // Biosaur2 is done: this window's phase B starts now, not after the slowest window's phase A
          {
#pragma omp task default(shared) firstprivate(o)
            phase_b(o);
          }
        }
        catch (...)
        {
          Internal::keepFirst(err);
        }
      }
    }
#pragma omp taskwait
    if (overlap) // Biosaur2 ran alongside: every window's phase B waits for it here; otherwise each window chained it (phase_b)
    {
      const bool phase_a_ok = ! err; // read once: phase-B tasks write err (under noims_err) while this loop submits
      for (long o = 0; o < (long)n_windows && phase_a_ok; ++o)
      {
#pragma omp task default(shared) firstprivate(o)
        phase_b(o);
      }
    }
#pragma omp taskwait
  } // the parallel region's barrier also waits for the phase-B tasks the windows chained
  omp_set_max_active_levels(max_levels);
  progress.endProgress();
  if (err) { std::rethrow_exception(err); }
  log("[noims] Biosaur2: " + DIAspeXtract::String(bio_v.size()) + " isotope-pattern features as precursors in " + DIAspeXtract::String(t_bio)
      + " s (" + (overlap ? "alongside" : "before") + " the MS2 traces)");

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
  OpenMS::Size n_sing = 0;
  for (OpenMS::Size w = 0; w < n_windows; ++w) { n_sing += nsing[w]; }
  log("[noims] " + DIAspeXtract::String(nf) + " fragments (" + DIAspeXtract::String(n_sing) + " single-point); window busy: MS2 traces " + DIAspeXtract::String(sm) + " s, envelopes and fragments "
      + DIAspeXtract::String(sf) + " s, spectra " + DIAspeXtract::String(st.busy_s - sm - sf) + " s");
  {
    OpenMS::Size wa = 0, wb = 0; // the windows with the longest phase A (traces, fragments) and phase B (spectra)
    for (OpenMS::Size w = 1; w < n_windows; ++w)
    {
      if (t_mtd[w] + t_frag[w] > t_mtd[wa] + t_frag[wa]) { wa = w; }
      if (busy[w] - t_mtd[w] - t_frag[w] > busy[wb] - t_mtd[wb] - t_frag[wb]) { wb = w; }
    }
    if (n_windows > 0)
    {
      log("[noims] slowest window: phase A " + DIAspeXtract::String(t_mtd[wa] + t_frag[wa]) + " s (MS2 traces " + DIAspeXtract::String(t_mtd[wa])
          + " s), phase B " + DIAspeXtract::String(busy[wb] - t_mtd[wb] - t_frag[wb]) + " s");
    }
  }
  st.wall_s = std::chrono::duration<double>(clk::now() - t0).count();
  log("[noims] " + DIAspeXtract::String(st.ms2_traces) + " MS2 traces, " + DIAspeXtract::String(st.precursors) + " precursors, "
      + DIAspeXtract::String(st.spectra) + " spectra, " + DIAspeXtract::String(st.peaks) + " peaks; " + DIAspeXtract::String(st.wall_s)
      + " s, window busy " + DIAspeXtract::String(st.busy_s) + " s on " + DIAspeXtract::String(st.threads) + " threads");
  return st;
}
} // namespace DIAspeXtract::NoIms
