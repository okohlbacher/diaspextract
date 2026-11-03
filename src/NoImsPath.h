// The extraction path for DIA input without ion mobility (Orbitrap, SWATH mzML). Separate from the IMS path by rule:
// it shares loading and writing with it, nothing else, and never changes what the IMS path emits.
//
// Precursors: Biosaur2 features of the whole run's MS1 (noims mono 2, the default; mono 0 = MS1 mass traces of each
// window's m/z slice grouped into isotope envelopes, charge 1-4, the mono the lightest trace of the group).
// Per isolation window, in parallel: (1) the window's precursors; (2) fragments = mass traces of the window's MS2
// scans; (3) a fragment joins a precursor when its profile correlates with the precursor's over the precursor's
// elution span (Pearson r >= noims:min_corr on the window's scan grid), or, with fewer than 3 grid points in that
// span, when its apex lies within one cycle of the precursor apex. A fragment's intensity is weighted by r^24 (Sage and MSFragger
// peptides and entrapment on three Orbitrap sets: stronger powers overfit Sage's scoring); the 200 heaviest are kept,
// one spectrum per precursor with >= 3 fragments.
#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include "PortableMath.h"
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
  using namespace OpenMS;

  struct Params
  {
    double ppm = 10.0;          ///< mass-trace and isotope m/z tolerance
    double min_corr = 0.3;      ///< fragment-precursor Pearson r on the window grid
    int max_peaks = 200;        ///< most intense fragments kept per spectrum (by weighted intensity)
    double ms1_noise = 0.0;     ///< MS1 peaks below this intensity are ignored (0 = MassTraceDetection's default)
    int ms1_min_scans = 2;      ///< MS1 traces need at least this many MS1 scans
    bool ms1_split = true;      ///< split MS1 traces at valleys (ElutionPeakDetection)
    double few_apex_cycles = 1.0; ///< precursors with < 3 window scans: fragment apex within this many cycles of the precursor apex
    int ms2_min_scans = 1;      ///< fragment traces need at least this many of the window's scans (1 = single points too; 0 = auto: 1 if sparse, else 2)
    int singles = -1;           ///< MS1 traces in no isotope envelope become precursors at z = 2 and z = 3 (-1 = auto: only if sparse)
    double weight_pow = 24.0;   ///< fragment intensity x r^weight_pow (few-point profiles: x exp(-(apex distance / cycle)^2)^(weight_pow/4))
    int min_frags = 3;          ///< fewest fragments a spectrum needs
    bool ms2_prec = true;       ///< isotope envelopes (z >= 2) of unfragmented precursors in the window's MS2 traces (DIA-Umpire Q3-like)
    int mono = 2;               ///< MS1 precursors: 0 greedy envelopes, 1 averagine shift test (avgEnvelopes), 2 Biosaur2 features
  };

  struct Stats { size_t windows = 0, ms1_traces = 0, ms2_traces = 0, precursors = 0, spectra = 0, peaks = 0; double busy_s = 0, wall_s = 0; int threads = 1; };

  namespace detail
  {
    constexpr double C13 = 1.0033548, PROTON = 1.007276;

    struct Prec { double mz, rt, rt_lo, rt_hi, inten; int z; std::vector<double> rt_pts, int_pts; };
    struct Frag { double mz, apex_rt, inten; std::vector<std::pair<double, double>> pts; };

    inline double apexRT(const MassTrace& t)
    {
      double best = -1, rt = t.getCentroidRT();
      for (const auto& p : t) if (p.getIntensity() > best) { best = p.getIntensity(); rt = p.getRT(); }
      return rt;
    }

    inline std::vector<MassTrace> traces(const PeakMap& m, double ppm, double min_len, double noise, bool split)
    {
      MassTraceDetection mtd;
      Param p = mtd.getDefaults();
      p.setValue("mass_error_ppm", ppm);
      p.setValue("min_trace_length", min_len);
      if (noise > 0) p.setValue("noise_threshold_int", noise);
      mtd.setParameters(p);
      std::vector<MassTrace> mt;
      mtd.run(m, mt);
      if (!split || mt.empty()) return mt;
      ElutionPeakDetection epd;
      Param q = epd.getDefaults();
      q.setValue("width_filtering", "off");
      epd.setParameters(q);
      std::vector<MassTrace> out;
      epd.detectPeaks(mt, out);
      return out;
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

    /// isotope envelopes from MS1 traces; a trace is a mono if it has a +1 isotope and no -1 isotope at the same charge
    // ponytail: greedy, highest charge first, no averagine check; add an intensity-pattern test if charge errors show up
    inline std::vector<Prec> envelopes(const std::vector<MassTrace>& mt, double ppm, double lo, double hi, double coelute,
                                       int zmin = 1, bool singles = false)
    {
      std::vector<size_t> o(mt.size());
      std::iota(o.begin(), o.end(), 0);
      std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return mt[a].getCentroidMZ() < mt[b].getCentroidMZ() || (mt[a].getCentroidMZ() == mt[b].getCentroidMZ() && a < b); });
      std::vector<double> mz(o.size()), art(o.size());
      for (size_t k = 0; k < o.size(); ++k) { mz[k] = mt[o[k]].getCentroidMZ(); art[k] = apexRT(mt[o[k]]); }
      auto partner = [&](size_t k, double target) -> long {   // co-eluting trace at target m/z
        const double tol = target * ppm * 1e-6;
        const MassTrace& a = mt[o[k]];
        const double half = std::max(a.getFWHM() * 0.5, coelute);   // isotopes share MS1 scans: one scan of apex jitter
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
      std::vector<char> used(o.size(), 0);   // part of an envelope (mono or isotope)
      auto make = [&](size_t k, int z) {
        const MassTrace& t = mt[o[k]];
        Prec p{mz[k], art[k], t.begin()->getRT(), t.rbegin()->getRT(), t.getMaxIntensity(false), z, {}, {}};
        for (const auto& pk : t) { p.rt_pts.push_back(pk.getRT()); p.int_pts.push_back(pk.getIntensity()); }
        out.push_back(std::move(p));
      };
      for (size_t k = 0; k < o.size(); ++k)
      {
        for (int z = 4; z >= zmin; --z)
        {
          const long nx = partner(k, mz[k] + C13 / z);
          if (nx < 0) continue;
          if (partner(k, mz[k] - C13 / z) >= 0) { used[k] = 1; break; }   // not the lightest isotope at this charge: not a mono
          used[k] = 1;
          for (long j = nx, i = 1; j >= 0 && i <= 3; j = partner((size_t)j, mz[j] + C13 / z), ++i) used[(size_t)j] = 1;
          if (mz[k] >= lo && mz[k] < hi) make(k, z);
          break;
        }
      }
      if (singles)
        for (size_t k = 0; k < o.size(); ++k)
          if (!used[k] && mz[k] >= lo && mz[k] < hi) { make(k, 2); make(k, 3); }
      return out;
    }

    /// [mono] averagine isotope k of a peptide of this neutral mass (Poisson, lambda = 0.000594 m: the IMS path's rule)
    inline double avg(double neutral, int k)
    {
      const double lam = std::max(0.10, 0.000594 * neutral);
      double p = pmath::exp(-lam);
      for (int i = 1; i <= k; ++i) p *= lam / i;
      return p;
    }

    /// [mono] isotope envelopes chosen by an averagine shift test (Dinosaur / Monocle): seeds brightest first; for every
    /// charge z at which the seed has a +-1 isotope partner and every monoisotope position c in [-3, 0] (the seed is then
    /// isotope -c), the summed co-eluting intensity at positions c-1 .. c+5 is compared with [0, A0 .. A5] (averagine at the
    /// mono's mass): score = cosine x the averagine fraction observed. The best (z, c) with an observed mono wins and claims
    /// the traces of its positions c .. c+3. The precursor's elution profile is the seed's (the brightest trace of the envelope).
    inline std::vector<Prec> avgEnvelopes(const std::vector<MassTrace>& mt, double ppm, double lo, double hi, double coelute,
                                          int zmin, bool singles)
    {
      std::vector<size_t> o(mt.size());
      std::iota(o.begin(), o.end(), 0);
      std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return mt[a].getCentroidMZ() < mt[b].getCentroidMZ() || (mt[a].getCentroidMZ() == mt[b].getCentroidMZ() && a < b); });
      std::vector<double> mz(o.size()), art(o.size()), in(o.size());
      for (size_t k = 0; k < o.size(); ++k) { mz[k] = mt[o[k]].getCentroidMZ(); art[k] = apexRT(mt[o[k]]); in[k] = mt[o[k]].getMaxIntensity(false); }
      // summed intensity of the traces co-eluting with seed k at its isotope position `pos` (charge z); `best` = the brightest
      auto PS = [&](size_t k, int z, int pos, long* best) -> double {
        const double target = mz[k] + pos * C13 / z, tol = target * ppm * 1e-6;
        const double half = std::max(mt[o[k]].getFWHM() * 0.5, coelute);
        double sum = 0, bi = -1;
        for (auto it = std::lower_bound(mz.begin(), mz.end(), target - tol); it != mz.end() && *it <= target + tol; ++it)
        {
          const size_t j = it - mz.begin();
          if (std::fabs(art[j] - art[k]) > half) continue;
          sum += in[j];
          if (in[j] > bi) { bi = in[j]; if (best) *best = (long)j; }
        }
        return sum;
      };
      std::vector<size_t> seeds(o.size());
      std::iota(seeds.begin(), seeds.end(), 0);
      std::sort(seeds.begin(), seeds.end(), [&](size_t a, size_t b) { return in[a] > in[b] || (in[a] == in[b] && a < b); });
      std::vector<Prec> out;
      std::vector<char> used(o.size(), 0);
      struct Hyp { double score; int z, c; long mono; };
      for (size_t k : seeds)
      {
        if (used[k]) continue;
        Hyp h1{0, 0, 0, -1};
        for (int z = 4; z >= zmin; --z)
        {
          double ps[11];   // positions -4 .. +6 relative to the seed
          long bt[11];
          for (int p = -4; p <= 6; ++p) { bt[p + 4] = -1; ps[p + 4] = p == 0 ? in[k] : PS(k, z, p, &bt[p + 4]); }
          if (ps[4 + 1] <= 0 && ps[4 - 1] <= 0) continue;   // no isotope partner at this charge
          for (int c = -3; c <= 0; ++c)
          {
            if (ps[c + 4] <= 0) continue;   // the mono must be observed
            const double neutral = z * (mz[k] + c * C13 / z - PROTON);
            double sab = 0, saa = 0, sbb = 0, seen = 0, all = 0;
            for (int i = -1; i <= 5; ++i)
            {
              const double a = i < 0 ? 0.0 : avg(neutral, i), b = ps[c + i + 4];
              sab += a * b; saa += a * a; sbb += b * b;
              if (i >= 0) { all += a; if (b > 0) seen += a; }
            }
            const double s = sab / std::sqrt(saa * sbb) * seen / all;
            const Hyp h{s, z, c, c == 0 ? (long)k : bt[c + 4]};
            if (s > h1.score) h1 = h;
          }
        }
        if (h1.mono < 0) { if (singles && mz[k] >= lo && mz[k] < hi) for (int z : {2, 3}) { /* unpaired: as the greedy path */
            const MassTrace& t = mt[o[k]];
            Prec p{mz[k], art[k], t.begin()->getRT(), t.rbegin()->getRT(), in[k], z, {}, {}};
            for (const auto& pk : t) { p.rt_pts.push_back(pk.getRT()); p.int_pts.push_back(pk.getIntensity()); }
            out.push_back(std::move(p)); }
          used[k] = 1; continue; }
        if (const double m = mz[(size_t)h1.mono]; m >= lo && m < hi)
        {
          const MassTrace& t = mt[o[k]];
          Prec p{m, art[k], t.begin()->getRT(), t.rbegin()->getRT(), in[k], h1.z, {}, {}};
          for (const auto& pk : t) { p.rt_pts.push_back(pk.getRT()); p.int_pts.push_back(pk.getIntensity()); }
          out.push_back(std::move(p));
        }
        used[k] = 1;
        for (int i = 0; i <= 3; ++i)   // claim the winner's positions
        {
          long b = -1;
          if (h1.c + i == 0) continue;
          if (PS(k, h1.z, h1.c + i, &b) > 0 && b >= 0) used[(size_t)b] = 1;
        }
      }
      return out;
    }

    /// [mono C] precursors from OpenMS's Biosaur2 port, run once on the whole run's MS1 (averagine-checked isotope patterns)
    inline std::vector<Prec> biosaur(const PeakMap& exp, const std::vector<size_t>& ms1, double ppm, double mzlo, double mzhi)
    {
      PeakMap m;
      for (size_t i : ms1) m.addSpectrum(exp[i]);
      Biosaur2Algorithm b;
      Param p = b.getDefaults();
      p.setValue("htol", ppm); p.setValue("itol", ppm);
      p.setValue("minmz", mzlo); p.setValue("maxmz", mzhi);
      p.setValue("cmin", 1); p.setValue("cmax", 4);
      b.setParameters(p);
      b.setMSData(std::move(m));
      FeatureMap fm;
      std::vector<Biosaur2Algorithm::Hill> hills;
      std::vector<Biosaur2Algorithm::PeptideFeature> feats;
      b.run(fm, hills, feats);
      std::map<Size, size_t> at;
      for (size_t i = 0; i < hills.size(); ++i) at[hills[i].hill_idx] = i;
      std::vector<Prec> out;
      for (const auto& f : feats)
      {
        auto it = at.find(f.mono_hill_idx);
        if (it == at.end()) continue;
        const auto& h = hills[it->second];
        Prec p{f.mz, f.rt_apex, h.rt_start, h.rt_end, f.intensity_apex, f.charge, h.rt_values, h.intensities};
        out.push_back(std::move(p));
      }
      std::sort(out.begin(), out.end(), [](const Prec& a, const Prec& b) {   // a total order: Biosaur2's feature order is fixed, its ties are not
        return std::tie(a.mz, a.rt, a.z, a.inten, a.rt_lo, a.rt_hi, a.rt_pts, a.int_pts) < std::tie(b.mz, b.rt, b.z, b.inten, b.rt_lo, b.rt_hi, b.rt_pts, b.int_pts); });
      return out;
    }
  }

  /// @p exp: the resident run without ion mobility (MS1 and MS2 centroided). Fills @p out with one MS2 spectrum per precursor.
  inline Stats run(const PeakMap& exp, const Params& P, PeakMap& out, const std::function<void(const String&)>& log)
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
    for (const auto& kv : win) keys.push_back(kv.first);
    st.windows = keys.size();
    if (keys.empty() || ms1.empty()) { log("[noims] no MS1 scans or no MS2 isolation windows"); return st; }
    // a precursor belongs to the window containing it whose center is nearest (windows may overlap); a window with
    // fewer than 3 MS2 scans gives no spectra and so owns nothing
    auto owner = [&](double mz) {
      long best = -1; double bd = 1e9;
      for (size_t w = 0; w < keys.size(); ++w)
        if (mz >= keys[w].first && mz < keys[w].second && win.at(keys[w]).size() >= 3)
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
    const auto cyc_win = std::find_if(win.begin(), win.end(), [](const auto& kv) { return kv.second.size() >= 3; });
    const double cycle = spacing(cyc_win != win.end() ? cyc_win->second : win.begin()->second, 1.0);
    const double ms1_dt = spacing(ms1, cycle);
    log("[noims] " + String(keys.size()) + " isolation windows, " + String(ms1.size()) + " MS1 scans every " + String(ms1_dt)
        + " s, MS2 cycle " + String(cycle) + " s");

    // sparse: each window is sampled much less often than MS1 (narrow-window schemes): fragments seen once and precursors
    // without an isotope partner carry information there; on dense data they only add noise
    const bool sparse = cycle > 2.0 * ms1_dt;
    const bool singles = P.singles < 0 ? sparse : P.singles != 0;
    const int ms2_min = P.ms2_min_scans > 0 ? P.ms2_min_scans : (sparse ? 1 : 2);
    log(String("[noims] ") + (sparse ? "sparse" : "dense") + " sampling: single-scan fragments " + (ms2_min == 1 ? "on" : "off")
        + ", unpaired MS1 traces as z = 2/3 " + (singles ? "on" : "off"));
    std::vector<Prec> bio;   // [mono C] Biosaur2 features of the whole run
    if (P.mono == 2)
    {
      double top = 0;
      for (const auto& k : keys) top = std::max(top, k.second);
      const auto tb = clk::now();
      bio = biosaur(exp, ms1, P.ppm, keys.front().first - 1.0, top + 1.0);
      log("[noims] Biosaur2: " + String(bio.size()) + " isotope-pattern features as precursors in "
          + String(std::chrono::duration<double>(clk::now() - tb).count()) + " s");
    }
    std::vector<std::vector<MSSpectrum>> res(keys.size());
    std::vector<size_t> nprec(keys.size(), 0), ntr1(keys.size(), 0), ntr2(keys.size(), 0);
    std::vector<double> busy(keys.size(), 0.0);
    std::vector<char> short_win(keys.size(), 0);
    std::exception_ptr err;
    #pragma omp parallel for schedule(dynamic, 1)
    for (long w = 0; w < (long)keys.size(); ++w)
    try
    {
      const auto tw = clk::now();
      if (win.at(keys[w]).size() < 3) { short_win[w] = 1; continue; }   // MassTraceDetection needs 3 scans
      const double lo = keys[w].first, hi = keys[w].second;
      std::vector<Prec> prec;
      if (P.mono == 2)
      { for (const Prec& p : bio) if (p.mz >= lo && p.mz < hi) prec.push_back(p); }
      else
      {
        PeakMap m1;   // MS1 slice: the window's m/z range plus room for a mono's isotopes
        for (size_t i : ms1)
        {
          const MSSpectrum& s = exp[i];
          MSSpectrum c;
          c.setRT(s.getRT()); c.setMSLevel(1);
          for (auto it = s.MZBegin(lo - 1.5); it != s.end() && it->getMZ() < hi + 1.5; ++it) c.push_back(*it);
          m1.addSpectrum(std::move(c));
        }
        std::vector<MassTrace> mt1 = traces(m1, P.ppm, (P.ms1_min_scans - 1) * ms1_dt * 0.99, P.ms1_noise, P.ms1_split);
        ntr1[w] = mt1.size();
        prec = P.mono == 1 ? avgEnvelopes(mt1, P.ppm, lo, hi, 1.01 * ms1_dt, 1, singles)
                           : envelopes(mt1, P.ppm, lo, hi, 1.01 * ms1_dt, 1, singles);
      }
      prec.erase(std::remove_if(prec.begin(), prec.end(), [&](const Prec& p) { return owner(p.mz) != w; }), prec.end());
      nprec[w] = prec.size();

      // MS2 traces and the window's scan grid
      PeakMap m2;
      std::vector<double> grid;
      for (size_t i : win.at(keys[w]))
      { m2.addSpectrum(exp[i]); m2.getSpectra().back().setMSLevel(1); grid.push_back(exp[i].getRT()); }   // MassTraceDetection reads MS1 only
      std::vector<MassTrace> ft = traces(m2, P.ppm, (ms2_min - 1) * cycle * 0.99, 0.0, false);
      ntr2[w] = ft.size();
      m2.clear(true);
      if (P.ms2_prec)   // unfragmented precursors seen in MS2, unless MS1 already has them (same z, 1 cycle, ppm)
      {
        std::vector<Prec> q = envelopes(ft, P.ppm, lo, hi, 1.01 * cycle, 2, false);
        const size_t n1 = prec.size();
        for (auto& c : q)
        {
          if (owner(c.mz) != w) continue;
          bool dup = false;
          for (size_t i = 0; i < n1 && !dup; ++i)
            dup = prec[i].z == c.z && std::fabs(prec[i].mz - c.mz) <= c.mz * P.ppm * 1e-6 && std::fabs(prec[i].rt - c.rt) <= cycle;
          if (!dup) prec.push_back(std::move(c));
        }
        nprec[w] = prec.size();
      }
      std::vector<Frag> fr;
      fr.reserve(ft.size());
      for (const auto& t : ft)
      {
        Frag f{t.getCentroidMZ(), apexRT(t), t.getMaxIntensity(false), {}};
        for (const auto& pk : t) f.pts.emplace_back(pk.getRT(), pk.getIntensity());
        fr.push_back(std::move(f));
      }
      std::sort(fr.begin(), fr.end(), [](const Frag& a, const Frag& b) { return std::tie(a.apex_rt, a.mz, a.inten, a.pts) < std::tie(b.apex_rt, b.mz, b.inten, b.pts); });

      for (const Prec& p : prec)
      {
        const size_t g0 = std::lower_bound(grid.begin(), grid.end(), p.rt_lo) - grid.begin();
        const size_t g1 = std::upper_bound(grid.begin(), grid.end(), p.rt_hi) - grid.begin();
        std::vector<double> pp;
        for (size_t g = g0; g < g1; ++g) pp.push_back(interp(p.rt_pts, p.int_pts, grid[g]));
        const bool few = pp.size() < 3;
        const double neutral = p.z * (p.mz - PROTON);
        std::vector<std::pair<double, double>> pk;   // (intensity, m/z)
        const double a0 = few ? p.rt - P.few_apex_cycles * cycle : p.rt_lo, a1 = few ? p.rt + P.few_apex_cycles * cycle : p.rt_hi;
        const auto first = std::lower_bound(fr.begin(), fr.end(), a0, [](const Frag& f, double t) { return f.apex_rt < t; });
        for (size_t k = first - fr.begin(); k < fr.size() && fr[k].apex_rt <= a1; ++k)
        {
          const Frag& f = fr[k];
          if (f.mz - PROTON > neutral - 56.0) continue;   // heavier than the precursor allows at z = 1
          if (!few)
          {
            std::vector<double> fp(pp.size(), 0.0);
            for (const auto& q : f.pts)
            {
              const size_t g = std::lower_bound(grid.begin(), grid.end(), q.first - 1e-6) - grid.begin();
              if (g >= g0 && g < g1) fp[g - g0] = q.second;
            }
            const double r = pearson(pp, fp);
            if (r < P.min_corr) continue;
            pk.emplace_back(f.inten * pmath::pow(r, P.weight_pow), f.mz);
            continue;
          }
          const double d = (f.apex_rt - p.rt) / cycle;
          pk.emplace_back(f.inten * pmath::exp(-d * d * P.weight_pow / 4.0), f.mz);
        }
        if ((int)pk.size() < P.min_frags) continue;
        if ((int)pk.size() > P.max_peaks)
        { std::nth_element(pk.begin(), pk.begin() + P.max_peaks, pk.end(), std::greater<>()); pk.resize(P.max_peaks); }
        // by m/z, equal m/z by intensity: nth_element leaves the kept peaks in an unspecified order
        std::sort(pk.begin(), pk.end(), [](const auto& a, const auto& b) { return std::tie(a.second, a.first) < std::tie(b.second, b.first); });
        MSSpectrum s;
        s.setMSLevel(2);
        s.setType(SpectrumSettings::SpectrumType::CENTROID);
        s.setRT(p.rt);
        Precursor pr;
        pr.setMZ(p.mz); pr.setCharge(p.z); pr.setIntensity(p.inten);
        pr.setIsolationWindowLowerOffset(p.mz - lo); pr.setIsolationWindowUpperOffset(hi - p.mz);
        s.setPrecursors({pr});
        for (const auto& q : pk) { Peak1D x; x.setMZ(q.second); x.setIntensity(q.first); s.push_back(x); }
        res[w].push_back(std::move(s));
      }
      busy[w] = std::chrono::duration<double>(clk::now() - tw).count();
    }
    catch (...)
    {
      #pragma omp critical(noims_err)
      if (!err) err = std::current_exception();
    }
    if (err) std::rethrow_exception(err);

    for (size_t w = 0; w < keys.size(); ++w)
    {
      st.precursors += nprec[w]; st.busy_s += busy[w]; st.ms1_traces += ntr1[w]; st.ms2_traces += ntr2[w];
      for (auto& s : res[w]) { st.peaks += s.size(); out.addSpectrum(std::move(s)); }
    }
    // by RT; spectra of one RT keep the assembly order (window key, then precursor order), which is platform-independent
    std::stable_sort(out.getSpectra().begin(), out.getSpectra().end(), [](const MSSpectrum& a, const MSSpectrum& b) { return a.getRT() < b.getRT(); });
    for (size_t i = 0; i < out.size(); ++i) out[i].setNativeID("index=" + String(i));
    st.spectra = out.size();
    if (const auto n = std::count(short_win.begin(), short_win.end(), 1))
      log("[noims] " + String(n) + " of " + String(keys.size()) + " isolation windows have fewer than 3 MS2 scans: no spectra from them");
    st.wall_s = std::chrono::duration<double>(clk::now() - t0).count();
    if (P.mono) log("[noims] monoisotope: " + String(P.mono == 1 ? "averagine shift test" : "Biosaur2 features"));
    log("[noims] " + String(st.ms1_traces) + " MS1 traces (window slices), " + String(st.ms2_traces) + " MS2 traces, " + String(st.precursors) + " precursors, " + String(st.spectra) + " spectra, " + String(st.peaks) + " peaks; "
        + String(st.wall_s) + " s, window busy " + String(st.busy_s) + " s on " + String(st.threads) + " threads");
    return st;
  }
}
