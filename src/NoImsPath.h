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
#include "PortableMath.h"
#include "OmsString.h"
#include <omp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <functional>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

namespace spx_noims
{
using diaspextract::String;
  using namespace OpenMS;

  struct Params
  {
    double ppm = 10.0;          ///< mass-trace and isotope m/z tolerance
    double min_corr = 0.3;      ///< fragment-precursor Pearson r on the window grid
    // MS1 precursor detection (Biosaur2): the Orbitrap defaults; TOF input changes them (see noImsSpectra_)
    bool tof = false;           ///< Biosaur2's TOF noise filter (per-50-m/z-bin intensity floor)
    double min_intensity = 1.0; ///< Biosaur2 mini: least peak intensity
    int min_hill_scans = 2;     ///< Biosaur2 minlh: fewest MS1 scans per hill
    int min_charge = 1;         ///< Biosaur2 cmin
    int min_frags = 3;          ///< fewest fragments a spectrum needs
    double min_frag_intensity = 0; ///< least apex intensity of a fragment trace
    int min_frag_points = 1;    ///< fewest scans of a fragment trace
    int max_frags = 200;        ///< most intense fragments kept per spectrum (by weighted intensity)
    double weight_pow = 24.0;   ///< fragment intensity x r^weight_pow
    double ms2_ppm = 0;         ///< MS2 mass-trace tolerance; 0 = ppm
    double split = 0;           ///< split MS2 traces at elution valleys below split x the lower flanking maximum; 0 = off
    bool env_all_charges = false; ///< MS2 envelopes: a -1 isotope at one charge moves on to the next lower charge (else: no precursor)
    bool drop_precursor = false;  ///< drop the precursor's own isotope peaks (mono-1 .. mono+3) from its spectrum
  };

  struct Stats { size_t windows = 0, ms2_traces = 0, precursors = 0, spectra = 0, peaks = 0; double busy_s = 0, wall_s = 0; int threads = 1; };

  namespace detail
  {
    constexpr double C13 = 1.0033548, PROTON = 1.007276;
    constexpr double FEW_APEX_CYCLES = 1.0; ///< precursors with < 3 window scans: fragment apex within this many cycles of the precursor apex

    struct Prec { double mz, rt, rt_lo, rt_hi, inten; int z; std::vector<double> rt_pts, int_pts; };
    /// a fragment trace; `pts` are (RT, intensity), `gpts` the same points as (grid index, intensity)
    struct Frag { double mz, apex_rt, inten; std::vector<std::pair<double, double>> pts; std::vector<std::pair<size_t, double>> gpts; };

    inline double apexRT(const MassTrace& t)
    {
      double best = -1, rt = t.getCentroidRT();
      for (const auto& p : t) if (p.getIntensity() > best) { best = p.getIntensity(); rt = p.getRT(); }
      return rt;
    }

    inline std::vector<MassTrace> traces(const PeakMap& m, double ppm)
    {
      MassTraceDetection mtd;
      mtd.setLogType(ProgressLogger::NONE);   // windows run in parallel: progress lines would interleave
      Param p = mtd.getDefaults();
      p.setValue("mass_error_ppm", ppm);
      p.setValue("min_trace_length", 0.0);   // single points too: a fragment seen in one scan of a window
      mtd.setParameters(p);
      std::vector<MassTrace> mt;
      mtd.run(m, mt);
      return mt;
    }

    /// cut points (grid indices) of a trace's grid profile @p y (index g0 + k), at valleys of the 1-2-1 smoothed profile
    /// that fall below @p valley x both flanking maxima: one elution peak per segment
    inline std::vector<size_t> valleys(const std::vector<double>& y, size_t g0, double valley)
    {
      std::vector<size_t> cut;
      double M = 0, m = -1; size_t mi = 0;   // m < 0: still rising to the current maximum M
      for (size_t k = 0; k < y.size(); ++k)
      {
        const double v = 0.25 * ((k ? y[k - 1] : 0.0) + 2 * y[k] + (k + 1 < y.size() ? y[k + 1] : 0.0));
        if (m < 0) { if (v >= M) M = v; else { m = v; mi = k; } }
        else if (v < m) { m = v; mi = k; }
        else if (m < valley * M && m < valley * v) { cut.push_back(g0 + mi); M = v; m = -1; }
        else if (v > M) { M = v; m = -1; }
      }
      return cut;
    }

    /// Pearson r of two profiles sampled on the same grid
    inline double pearson(const std::vector<double>& a, const std::vector<double>& b)
    {
      const size_t n = a.size();
      double ma = 0, mb = 0;
      for (size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
      ma /= n; mb /= n;
      double sab = 0, saa = 0, sbb = 0;
      for (size_t i = 0; i < n; ++i) { sab += (a[i] - ma) * (b[i] - mb); saa += (a[i] - ma) * (a[i] - ma); sbb += (b[i] - mb) * (b[i] - mb); }
      return (saa > 0 && sbb > 0) ? sab / std::sqrt(saa * sbb) : 0.0;
    }

    /// linear interpolation of (x, y) at t, 0 outside
    inline double interp(const std::vector<double>& x, const std::vector<double>& y, double t)
    {
      if (x.empty() || t < x.front() || t > x.back()) return 0.0;
      const size_t j = std::upper_bound(x.begin(), x.end(), t) - x.begin();   // >= 1, since t >= x.front()
      if (j >= x.size()) return y.back();
      const double f = (t - x[j - 1]) / (x[j] - x[j - 1]);
      return y[j - 1] + f * (y[j] - y[j - 1]);
    }

    /// isotope envelopes (z = 2..4) from MS2 traces; a trace is a mono if it has a +1 isotope and no -1 isotope at the same charge
    // ponytail: greedy, highest charge first, no averagine check; add an intensity-pattern test if charge errors show up
    inline std::vector<Prec> envelopes(const std::vector<MassTrace>& mt, double ppm, double lo, double hi, double coelute, bool all_charges = false)
    {
      std::vector<size_t> o(mt.size());
      std::iota(o.begin(), o.end(), 0);
      std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return mt[a].getCentroidMZ() < mt[b].getCentroidMZ() || (mt[a].getCentroidMZ() == mt[b].getCentroidMZ() && a < b); });
      std::vector<double> mz(o.size()), art(o.size());
      for (size_t k = 0; k < o.size(); ++k) { mz[k] = mt[o[k]].getCentroidMZ(); art[k] = apexRT(mt[o[k]]); }
      auto partner = [&](size_t k, double target) -> long {   // co-eluting trace at target m/z
        const double tol = target * ppm * 1e-6;
        const MassTrace& a = mt[o[k]];
        const double half = std::max(a.getFWHM() * 0.5, coelute);   // isotopes share scans: one scan of apex jitter
        long best = -1; double bi = -1;
        for (auto it = std::lower_bound(mz.begin(), mz.end(), target - tol); it != mz.end() && *it <= target + tol; ++it)
        {
          const size_t j = it - mz.begin();
          if (std::fabs(art[j] - art[k]) > half) continue;
          const double in = mt[o[j]].getMaxIntensity(false);
          if (in > bi) { bi = in; best = (long)j; }
        }
        return best;
      };
      std::vector<Prec> out;
      for (size_t k = 0; k < o.size(); ++k)
      {
        for (int z = 4; z >= 2; --z)
        {
          if (partner(k, mz[k] + C13 / z) < 0) continue;
          if (partner(k, mz[k] - C13 / z) >= 0) { if (all_charges) continue; break; }   // not the lightest isotope at this charge
          if (mz[k] >= lo && mz[k] < hi)
          {
            const MassTrace& t = mt[o[k]];
            Prec p{mz[k], art[k], t.begin()->getRT(), t.rbegin()->getRT(), t.getMaxIntensity(false), z, {}, {}};
            for (const auto& pk : t) { p.rt_pts.push_back(pk.getRT()); p.int_pts.push_back(pk.getIntensity()); }
            out.push_back(std::move(p));
          }
          break;
        }
      }
      return out;
    }

    /// precursors from OpenMS's Biosaur2 port, run once on the whole run's MS1 (averagine-checked isotope patterns)
    inline std::vector<Prec> biosaur(PeakMap&& ms1, const Params& P, double mzlo, double mzhi)
    {
      std::vector<Biosaur2Algorithm::Hill> hills;
      std::vector<Biosaur2Algorithm::PeptideFeature> feats;
      {
        Biosaur2Algorithm b;
        Param p = b.getDefaults();
        p.setValue("htol", P.ppm); p.setValue("itol", P.ppm);
        p.setValue("minmz", mzlo); p.setValue("maxmz", mzhi);
        p.setValue("cmin", P.min_charge); p.setValue("cmax", 4);
        p.setValue("mini", P.min_intensity); p.setValue("minlh", P.min_hill_scans);
        p.setValue("tof", P.tof ? "true" : "false");
        // no per-peak ion mobility on this path: a FAIMS run's per-spectrum CV made Biosaur2 switch on its PASEF
        // centroiding and IM gating, both no-ops without an IM array but logged as if they ran
        p.setValue("paseftol", 0.0);
        b.setParameters(p);
        b.setMSData(std::move(ms1));
        FeatureMap fm;   // unused; freed with the algorithm's copy of the MS1 data
        b.run(fm, hills, feats);
      }
      Size top = 0;
      for (const auto& h : hills) top = std::max(top, h.hill_idx);
      std::vector<size_t> at(hills.empty() ? 0 : top + 1, SIZE_MAX);   // hill_idx -> position in hills (last one wins, as a map would)
      for (size_t i = 0; i < hills.size(); ++i) at[hills[i].hill_idx] = i;
      std::vector<Prec> out;
      out.reserve(feats.size());
      for (const auto& f : feats)
      {
        if (f.mono_hill_idx >= at.size() || at[f.mono_hill_idx] == SIZE_MAX) continue;
        const auto& h = hills[at[f.mono_hill_idx]];
        out.push_back(Prec{f.mz, f.rt_apex, h.rt_start, h.rt_end, f.intensity_apex, f.charge, h.rt_values, h.intensities});
      }
      std::sort(out.begin(), out.end(), [](const Prec& a, const Prec& b) {   // a total order: Biosaur2's feature order is fixed, its ties are not
        return std::tie(a.mz, a.rt, a.z, a.inten, a.rt_lo, a.rt_hi, a.rt_pts, a.int_pts) < std::tie(b.mz, b.rt, b.z, b.inten, b.rt_lo, b.rt_hi, b.rt_pts, b.int_pts); });
      return out;
    }

    /// the spectrum of precursor @p p from the window's fragments (sorted by apex RT), or an empty one (no MS level)
    inline MSSpectrum assemble(const Prec& p, const std::vector<Frag>& fr, const std::vector<double>& grid, double cycle,
                               const Params& P, double lo, double hi)
    {
      MSSpectrum s;
      const size_t g0 = std::lower_bound(grid.begin(), grid.end(), p.rt_lo) - grid.begin();
      const size_t g1 = std::upper_bound(grid.begin(), grid.end(), p.rt_hi) - grid.begin();
      std::vector<double> pp;
      for (size_t g = g0; g < g1; ++g) pp.push_back(interp(p.rt_pts, p.int_pts, grid[g]));
      const bool few = pp.size() < 3;
      const double neutral = p.z * (p.mz - PROTON);
      std::vector<std::pair<double, double>> pk;   // (intensity, m/z)
      std::vector<double> fp(pp.size());
      // split traces (TOF) are short against the precursor span: their r comes from their own points (one-pass sums).
      // ponytail: unsplit traces keep the two-pass dense form, so Orbitrap output stays bit-identical
      const bool sparse = P.split > 0;
      double sa = 0, saa = 0;
      for (double a : pp) { sa += a; saa += a * a; }
      const double nn = (double)pp.size(), va = nn * saa - sa * sa;
      const double a0 = few ? p.rt - FEW_APEX_CYCLES * cycle : p.rt_lo, a1 = few ? p.rt + FEW_APEX_CYCLES * cycle : p.rt_hi;
      const auto first = std::lower_bound(fr.begin(), fr.end(), a0, [](const Frag& f, double t) { return f.apex_rt < t; });
      for (size_t k = first - fr.begin(); k < fr.size() && fr[k].apex_rt <= a1; ++k)
      {
        const Frag& f = fr[k];
        if (f.mz - PROTON > neutral - 56.0) continue;   // heavier than the precursor allows at z = 1
        if (P.drop_precursor)
        {
          const double iso = std::round((f.mz - p.mz) * p.z / C13);
          if (iso >= -1 && iso <= 3 && std::fabs(f.mz - p.mz - iso * C13 / p.z) <= f.mz * P.ppm * 1e-6) continue;   // unfragmented precursor
        }
        if (!few)
        {
          double r;
          if (sparse)
          {
            double sb = 0, sbb = 0, sab = 0;
            for (auto it = std::lower_bound(f.gpts.begin(), f.gpts.end(), std::make_pair(g0, -1.0)); it != f.gpts.end() && it->first < g1; ++it)
            { sb += it->second; sbb += it->second * it->second; sab += pp[it->first - g0] * it->second; }
            const double vb = nn * sbb - sb * sb;
            r = (va > 0 && vb > 0) ? (nn * sab - sa * sb) / std::sqrt(va * vb) : 0.0;
          }
          else
          {
            std::fill(fp.begin(), fp.end(), 0.0);
            for (const auto& q : f.gpts) if (q.first >= g0 && q.first < g1) fp[q.first - g0] = q.second;
            r = pearson(pp, fp);
          }
          if (r < P.min_corr) continue;
          pk.emplace_back(f.inten * pmath::pow(r, P.weight_pow), f.mz);
          continue;
        }
        const double d = (f.apex_rt - p.rt) / cycle;
        pk.emplace_back(f.inten * pmath::exp(-d * d * P.weight_pow / 4.0), f.mz);
      }
      if ((int)pk.size() < P.min_frags) return s;
      if ((int)pk.size() > P.max_frags)
      { std::nth_element(pk.begin(), pk.begin() + P.max_frags, pk.end(), std::greater<>()); pk.resize(P.max_frags); }
      // by m/z, equal m/z by intensity: nth_element leaves the kept peaks in an unspecified order
      std::sort(pk.begin(), pk.end(), [](const auto& a, const auto& b) { return std::tie(a.second, a.first) < std::tie(b.second, b.first); });
      s.setMSLevel(2);
      s.setType(SpectrumSettings::SpectrumType::CENTROID);
      s.setRT(p.rt);
      Precursor pr;
      pr.setMZ(p.mz); pr.setCharge(p.z); pr.setIntensity(p.inten);
      pr.setIsolationWindowLowerOffset(p.mz - lo); pr.setIsolationWindowUpperOffset(hi - p.mz);
      s.setPrecursors({pr});
      s.reserve(pk.size());
      for (const auto& q : pk) { Peak1D x; x.setMZ(q.second); x.setIntensity(q.first); s.push_back(x); }
      return s;
    }
  }

  /// @p exp: the resident run without ion mobility (MS1 and MS2 centroided, one isolation window per MS2 scan, in RT
  /// order). Its spectra are moved out: @p exp is left empty. Fills @p out with one MS2 spectrum per precursor.
  inline Stats run(PeakMap& exp, const Params& P, PeakMap& out, const std::function<void(const String&)>& log)
  {
    using namespace detail;
    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    Stats st;
    st.threads = omp_get_max_threads();

    // windows from the MS2 isolation windows; MS1 scans indexed once
    std::map<std::pair<double, double>, std::vector<size_t>> win;   // (lo, hi) -> MS2 scan indices
    std::vector<size_t> ms1;
    for (size_t i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum& s = exp[i];
      if (s.getMSLevel() == 1) ms1.push_back(i);
      else if (s.getMSLevel() == 2 && !s.getPrecursors().empty())
      {
        const Precursor& pr = s.getPrecursors()[0];
        const double lo = std::round((pr.getMZ() - pr.getIsolationWindowLowerOffset()) * 100) / 100;
        const double hi = std::round((pr.getMZ() + pr.getIsolationWindowUpperOffset()) * 100) / 100;
        win[{lo, hi}].push_back(i);
      }
    }
    std::vector<std::pair<double, double>> keys;
    std::vector<std::vector<size_t>> scans;
    for (auto& kv : win) { keys.push_back(kv.first); scans.push_back(std::move(kv.second)); }
    win.clear();
    const size_t W = keys.size();
    st.windows = W;
    if (W == 0 || ms1.empty()) { log("[noims] no MS1 scans or no MS2 isolation windows"); exp.clear(true); return st; }
    // a precursor belongs to the window containing it whose center is nearest (windows may overlap); a window with
    // fewer than 3 MS2 scans gives no spectra and so owns nothing
    auto owner = [&](double mz) {
      long best = -1; double bd = 1e9;
      for (size_t w = 0; w < W; ++w)
        if (mz >= keys[w].first && mz < keys[w].second && scans[w].size() >= 3)
        { const double d = std::fabs(mz - 0.5 * (keys[w].first + keys[w].second)); if (d < bd) { bd = d; best = (long)w; } }
      return best;
    };
    // median spacing of these scans (fallback when fewer than two)
    auto spacing = [&](const std::vector<size_t>& v, double fallback) {
      std::vector<double> d;
      for (size_t k = 1; k < v.size(); ++k) d.push_back(exp[v[k]].getRT() - exp[v[k - 1]].getRT());
      if (d.empty()) return fallback;
      std::nth_element(d.begin(), d.begin() + d.size() / 2, d.end());
      return d[d.size() / 2];
    };
    // MS2 cycle: the scans of the first window with >= 3 of them (a short window would give a wrong cycle to all others)
    const auto cyc_win = std::find_if(scans.begin(), scans.end(), [](const auto& v) { return v.size() >= 3; });
    const double cycle = spacing(cyc_win != scans.end() ? *cyc_win : scans.front(), 1.0);
    log("[noims] " + String(W) + " isolation windows, " + String(ms1.size()) + " MS1 scans every " + String(spacing(ms1, cycle))
        + " s, MS2 cycle " + String(cycle) + " s");

    // the run's spectra move into Biosaur2's MS1 map and the windows' MS2 maps (MassTraceDetection reads MS1 only)
    std::vector<PeakMap> m2(W);
    std::vector<std::vector<double>> grid(W);
    std::vector<size_t> wpeaks(W, 0);
    for (size_t w = 0; w < W; ++w)
      for (size_t i : scans[w])
      {
        grid[w].push_back(exp[i].getRT());
        wpeaks[w] += exp[i].size();
        exp[i].setMSLevel(1);
        m2[w].addSpectrum(std::move(exp[i]));
      }
    PeakMap m1;
    for (size_t i : ms1) m1.addSpectrum(std::move(exp[i]));
    exp.clear(true);

    double top = 0;
    for (const auto& k : keys) top = std::max(top, k.second);
    const auto tb = clk::now();
    const std::vector<Prec> bio = biosaur(std::move(m1), P, keys.front().first - 1.0, top + 1.0);
    log("[noims] Biosaur2: " + String(bio.size()) + " isotope-pattern features as precursors in "
        + String(std::chrono::duration<double>(clk::now() - tb).count()) + " s");

    std::vector<size_t> order(W);   // heaviest window first, so the largest ones do not start last
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return wpeaks[a] > wpeaks[b]; });
    std::vector<std::vector<MSSpectrum>> res(W);
    std::vector<size_t> nprec(W, 0), ntr2(W, 0);
    std::vector<double> busy(W, 0.0), t_mtd(W, 0.0), t_frag(W, 0.0);
    std::vector<size_t> nfrag(W, 0);
    std::exception_ptr err;
    #pragma omp parallel for schedule(dynamic, 1)
    for (long o = 0; o < (long)W; ++o)
    try
    {
      const size_t w = order[o];
      if (scans[w].size() < 3) continue;   // MassTraceDetection needs 3 scans
      const auto tw = clk::now();
      const double lo = keys[w].first, hi = keys[w].second;
      std::vector<Prec> prec;
      for (auto it = std::lower_bound(bio.begin(), bio.end(), lo, [](const Prec& p, double m) { return p.mz < m; });
           it != bio.end() && it->mz < hi; ++it)
        if (owner(it->mz) == (long)w) prec.push_back(*it);

      // MS2 traces on the window's scan grid
      const std::vector<MassTrace> ft = traces(m2[w], P.ms2_ppm > 0 ? P.ms2_ppm : P.ppm);
      ntr2[w] = ft.size();
      t_mtd[w] = std::chrono::duration<double>(clk::now() - tw).count();
      m2[w].clear(true);
      // unfragmented precursors seen in MS2, unless MS1 already has them (same z, 1 cycle, ppm)
      std::vector<Prec> q = envelopes(ft, P.ppm, lo, hi, 1.01 * cycle, P.env_all_charges);
      const size_t n1 = prec.size();
      for (auto& c : q)
      {
        if (owner(c.mz) != (long)w) continue;
        bool dup = false;
        for (size_t i = 0; i < n1 && !dup; ++i)
          dup = prec[i].z == c.z && std::fabs(prec[i].mz - c.mz) <= c.mz * P.ppm * 1e-6 && std::fabs(prec[i].rt - c.rt) <= cycle;
        if (!dup) prec.push_back(std::move(c));
      }
      nprec[w] = prec.size();
      const std::vector<double>& gr = grid[w];
      std::vector<Frag> fr;
      fr.reserve(ft.size());
      std::vector<double> prof;
      for (const auto& t : ft)
      {
        if (t.getMaxIntensity(false) < P.min_frag_intensity || (int)t.getSize() < P.min_frag_points) continue;
        Frag f{t.getCentroidMZ(), apexRT(t), t.getMaxIntensity(false), {}, {}};
        for (const auto& pk : t)
        {
          f.pts.emplace_back(pk.getRT(), pk.getIntensity());
          f.gpts.emplace_back(std::lower_bound(gr.begin(), gr.end(), pk.getRT() - 1e-6) - gr.begin(), pk.getIntensity());
        }
        std::vector<size_t> cut;
        if (P.split > 0 && f.gpts.size() >= 3)
        {
          const size_t g0 = f.gpts.front().first;
          prof.assign(f.gpts.back().first - g0 + 1, 0.0);
          for (const auto& q : f.gpts) prof[q.first - g0] = q.second;
          cut = valleys(prof, g0, P.split);
        }
        if (cut.empty()) { fr.push_back(std::move(f)); continue; }
        // one fragment per elution peak: its own apex, intensity and intensity-weighted m/z
        cut.push_back(SIZE_MAX);
        size_t k = 0;
        for (size_t c : cut)
        {
          Frag g{0, 0, -1, {}, {}};
          double sw = 0;
          for (auto it = t.begin() + k; k < f.gpts.size() && f.gpts[k].first < c; ++k, ++it)
          {
            g.pts.push_back(f.pts[k]); g.gpts.push_back(f.gpts[k]);
            g.mz += it->getMZ() * it->getIntensity(); sw += it->getIntensity();
            if (it->getIntensity() > g.inten) { g.inten = it->getIntensity(); g.apex_rt = it->getRT(); }
          }
          if (g.gpts.empty()) continue;
          g.mz = sw > 0 ? g.mz / sw : t.getCentroidMZ();
          if (g.inten >= P.min_frag_intensity && (int)g.gpts.size() >= P.min_frag_points) fr.push_back(std::move(g));
        }
      }
      std::sort(fr.begin(), fr.end(), [](const Frag& a, const Frag& b) { return std::tie(a.apex_rt, a.mz, a.inten, a.pts) < std::tie(b.apex_rt, b.mz, b.inten, b.pts); });
      for (auto& f : fr) std::vector<std::pair<double, double>>().swap(f.pts);   // pts only break sort ties; gpts serve the rest
      nfrag[w] = fr.size();
      t_frag[w] = std::chrono::duration<double>(clk::now() - tw).count() - t_mtd[w];

      // one slot per precursor, filled by tasks that idle threads can take
      std::vector<MSSpectrum> slot(prec.size());
      std::exception_ptr terr;
      const long np = (long)prec.size();
      #pragma omp taskloop grainsize(256) default(shared)
      for (long i = 0; i < np; ++i)
      {
        try { slot[i] = assemble(prec[i], fr, gr, cycle, P, lo, hi); }
        catch (...)
        {
          #pragma omp critical(noims_err)
          if (!terr) terr = std::current_exception();
        }
      }
      if (terr) std::rethrow_exception(terr);
      for (auto& s : slot) if (s.getMSLevel() == 2) res[w].push_back(std::move(s));
      busy[w] = std::chrono::duration<double>(clk::now() - tw).count();
    }
    catch (...)
    {
      #pragma omp critical(noims_err)
      if (!err) err = std::current_exception();
    }
    if (err) std::rethrow_exception(err);

    for (size_t w = 0; w < W; ++w)
    {
      st.precursors += nprec[w]; st.busy_s += busy[w]; st.ms2_traces += ntr2[w];
      for (auto& s : res[w]) { st.peaks += s.size(); out.addSpectrum(std::move(s)); }
    }
    // by RT; spectra of one RT keep the assembly order (window key, then precursor order), which is platform-independent
    std::stable_sort(out.getSpectra().begin(), out.getSpectra().end(), [](const MSSpectrum& a, const MSSpectrum& b) { return a.getRT() < b.getRT(); });
    for (size_t i = 0; i < out.size(); ++i) out[i].setNativeID("index=" + String(i));
    st.spectra = out.size();
    if (const auto n = std::count_if(scans.begin(), scans.end(), [](const auto& v) { return v.size() < 3; }))
      log("[noims] " + String(n) + " of " + String(W) + " isolation windows have fewer than 3 MS2 scans: no spectra from them");
    double sm = 0, sf = 0; size_t nf = 0;
    for (size_t w = 0; w < W; ++w) { sm += t_mtd[w]; sf += t_frag[w]; nf += nfrag[w]; }
    log("[noims] " + String(nf) + " fragments; window busy: MS2 traces " + String(sm) + " s, envelopes and fragments " + String(sf)
        + " s, spectra " + String(st.busy_s - sm - sf) + " s");
    st.wall_s = std::chrono::duration<double>(clk::now() - t0).count();
    log("[noims] " + String(st.ms2_traces) + " MS2 traces, " + String(st.precursors) + " precursors, " + String(st.spectra) + " spectra, " + String(st.peaks) + " peaks; "
        + String(st.wall_s) + " s, window busy " + String(st.busy_s) + " s on " + String(st.threads) + " threads");
    return st;
  }
}
