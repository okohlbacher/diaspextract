// Golden test for the SHIPPED C++ calibration (src/TdfMzCalibration.h) against values produced by
// Bruker's own timsdata library. Runs anywhere: no vendor library, no cluster, no raw data, no
// OpenMS -- only a C++17 compiler. Build: c++ -std=c++17 -I src tests/test_calibration_cpp.cpp
#include "TdfMzCalibration.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using Cal = diaspextract::TdfMzCalibration;

static bool has(const std::string& s, const char* part) { return s.find(part) != std::string::npos; }

/// Dataset D's row (the first golden file): valid, so each negative path below changes exactly one thing.
static Cal validD()
{
  Cal c;
  c.model_type = 1; c.digitizer_timebase = 0.125; c.digitizer_delay = 25655.375;
  c.C0 = 279.3262846272992; c.C1 = 155279.13067653627; c.C2 = 0.001260061434461731;
  c.T1_ref = 25.693668980735552; c.dC1 = 20.0;
  return c;
}

/// PXD029836 run 1418's row (timsTOF Pro 2, acquisition software 2.0.53): the negative C2 the golden set pins.
static Cal pxd1418()
{
  Cal c;
  c.model_type = 1; c.digitizer_timebase = 0.19999999999999998; c.digitizer_delay = 25779.8;
  c.C0 = 314.1341896383674; c.C1 = 154199.217568937; c.C2 = -0.0010422663259118895;
  c.T1_ref = 25.618823611312585; c.dC1 = 20.0;
  return c;
}

/// The reason text's number format (TdfMzCalibration::num), so a test can find the value it computed itself.
static std::string g6(double v)
{
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

// Minimal extraction from the golden JSON: it is machine-generated with a fixed shape, so a tiny
// scanner beats adding a JSON dependency to a test that must build everywhere.
static bool nextNumber(const std::string& s, const std::string& key, size_t& pos, double& out)
{
  const size_t k = s.find("\"" + key + "\":", pos);
  if (k == std::string::npos) return false;
  size_t v = s.find(':', k) + 1;
  while (v < s.size() && (s[v] == ' ' || s[v] == '\n')) ++v;
  out = std::strtod(s.c_str() + v, nullptr);
  pos = v;
  return true;
}

int main(int argc, char** argv)
{
  const std::string path = argc > 1 ? argv[1] : "tests/calibration_golden.json";
  std::ifstream f(path);
  if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); return 2; }
  std::stringstream ss; ss << f.rdbuf();
  const std::string j = ss.str();
  // The scanner below reads eight calibration keys and leaves C4 at 0 (true of every golden file so far). A golden
  // file that carries C4 must not be tested silently at C4 = 0: teach the scanner the key first (review 2026-09-17).
  if (j.find("\"C4\":") != std::string::npos)
  { std::fprintf(stderr, "%s carries a C4 key, which this test does not read yet\n", path.c_str()); return 1; }

  size_t pos = 0;
  int n_files = 0, n_cases = 0, n_negative = 0;
  double worst = 0.0, worst_noc2 = 0.0, worst_notemp = 0.0, worst_back = 0.0;
  while (true)
  {
    const size_t fstart = j.find("\"model_type\":", pos);
    if (fstart == std::string::npos) break;
    const size_t name_at = j.rfind("\"file\": \"", fstart);
    const std::string name = name_at == std::string::npos ? "?" : j.substr(name_at + 9, j.find('"', name_at + 9) - name_at - 9);
    diaspextract::TdfMzCalibration cal;
    size_t p = fstart;
    double v = 0;
    nextNumber(j, "model_type", p, v);        cal.model_type = (int)v;
    nextNumber(j, "timebase", p, v);          cal.digitizer_timebase = v;
    nextNumber(j, "delay", p, v);             cal.digitizer_delay = v;
    nextNumber(j, "C0", p, v);                cal.C0 = v;
    nextNumber(j, "C1", p, v);                cal.C1 = v;
    nextNumber(j, "C2", p, v);                cal.C2 = v;
    nextNumber(j, "T1_ref", p, v);            cal.T1_ref = v;
    nextNumber(j, "dC1", p, v);               cal.dC1 = v;
    if (!cal.isSupported())
    { std::fprintf(stderr, "golden file %d unsupported: %s\n", n_files, cal.unsupportedReason().c_str()); return 1; }

    const size_t fend = j.find("\"model_type\":", fstart + 10);
    size_t cp = j.find("\"cases\":", fstart);
    int f_cases = 0;
    double f_worst = 0.0, f_noc2 = 0.0, f_flip = 0.0, f_notemp = 0.0;
    while (true)
    {
      const size_t c = j.find("\"frame\":", cp);
      if (c == std::string::npos || (fend != std::string::npos && c > fend)) break;
      size_t q = c;
      double frame = 0, t1 = 0, tof = 0, mz = 0;
      nextNumber(j, "frame", q, frame);
      nextNumber(j, "t1", q, t1);
      nextNumber(j, "tof", q, tof);
      nextNumber(j, "mz", q, mz);
      if (!cal.frameReason(t1).empty())
      { std::fprintf(stderr, "FAIL %s frame %.0f refused: %s\n", name.c_str(), frame, cal.frameReason(t1).c_str()); return 1; }
      const double b = cal.frameFactor(t1);
      const double got = cal.tofToMz(tof, b);
      const double ppm = std::fabs(got - mz) / mz * 1e6;
      if (!(ppm <= 1e-4))   // same tolerance as the python twin; measured max is 2.6e-5 (NaN fails too)
      {
        std::fprintf(stderr, "FAIL %s frame %.0f tof %.0f: got %.6f want %.6f (%.6f ppm)\n",
                     name.c_str(), frame, tof, got, mz, ppm);
        return 1;
      }
      f_worst = std::fmax(f_worst, ppm);
      // round trip must return the same TOF, to the dnoise port's tolerance (recover::kTofTol)
      const double back = cal.mzToTof(got, b);
      worst_back = std::fmax(worst_back, std::fabs(back - tof));
      if (!(std::fabs(back - tof) <= 1e-6))
      { std::fprintf(stderr, "FAIL %s round trip: tof %.4f -> mz -> %.10f\n", name.c_str(), tof, back); return 1; }
      // ablations: dropping C2, flipping its sign or dropping the temperature term must be visible, else the golden
      // set proves nothing
      diaspextract::TdfMzCalibration no_c2 = cal; no_c2.C2 = 0.0;
      diaspextract::TdfMzCalibration flip = cal; flip.C2 = -cal.C2;
      f_noc2 = std::fmax(f_noc2, std::fabs(no_c2.tofToMz(tof, b) - mz) / mz * 1e6);
      f_flip = std::fmax(f_flip, std::fabs(flip.tofToMz(tof, b) - mz) / mz * 1e6);
      f_notemp = std::fmax(f_notemp, std::fabs(cal.tofToMz(tof, cal.frameFactor(cal.T1_ref)) - mz) / mz * 1e6);
      ++n_cases; ++f_cases;
      cp = q;
    }
    // per file, so a file whose cases cannot tell the model from the known-bad ports fails on its own; the temperature
    // floor is per file too (the smallest measured is 0.074 ppm, PXD029836 1408: its frames lie within 0.004 K of T1)
    if (f_cases < 20 || !(f_noc2 >= 5.0) || !(f_flip >= 5.0) || !(f_notemp >= 0.05))
    { std::fprintf(stderr, "FAIL %s: %d cases, C2 ablation %.3f ppm, sign flip %.3f ppm, temperature ablation %.4f ppm\n",
                   name.c_str(), f_cases, f_noc2, f_flip, f_notemp); return 1; }
    worst_notemp = std::fmax(worst_notemp, f_notemp);
    std::printf("    %-16s C2 %+.6e  %2d cases  max %.7f ppm  (C2 dropped %.2f, flipped %.2f, temperature dropped %.3f ppm)\n",
                name.c_str(), cal.C2, f_cases, f_worst, f_noc2, f_flip, f_notemp);
    worst = std::fmax(worst, f_worst);
    worst_noc2 = std::fmax(worst_noc2, f_noc2);
    n_negative += cal.C2 < 0.0;
    ++n_files;
    pos = fstart + 10;
  }
  // 60 cases on dataset D's vector, 210 on six PXD029836 runs (four negative-C2 vectors)
  if (n_files < 9 || n_cases < 270 || n_negative < 6)
  { std::fprintf(stderr, "only %d golden cases in %d files (%d with C2 < 0) parsed\n", n_cases, n_files, n_negative); return 1; }
  if (worst_noc2 < 5.0)
  { std::fprintf(stderr, "C2 ablation only %.3f ppm -- golden set cannot catch the known-bad port\n", worst_noc2); return 1; }
  if (worst_notemp < 0.1)
  { std::fprintf(stderr, "temperature ablation only %.4f ppm -- the golden set has no temperature spread\n", worst_notemp); return 1; }

  // negative paths must be REJECTED, not approximated -- each one a valid row with exactly one defect, and the
  // reason must name it (a row refused for something else would pass a bare isSupported() check vacuously)
  if (!validD().isSupported()) { std::fprintf(stderr, "dataset D's row must be supported: %s\n", validD().unsupportedReason().c_str()); return 1; }
  {
    struct Neg { const char* what; Cal cal; const char* reason; };
    std::vector<Neg> negs;
    Cal c = validD(); c.model_type = 2;          negs.push_back({"ModelType 2", c, "ModelType 2"});
    c = validD(); c.dC2 = 1e-9;                   negs.push_back({"dC2 != 0", c, "dC2 != 0"});
    c = validD(); c.C3 = 1e-9;                    negs.push_back({"C3 != 0", c, "C3 != 0"});
    c = validD(); c.C1 = 1e-300;                  negs.push_back({"implausible C1", c, "C1 outside"});
    c = validD(); c.C1 = 0.0;                     negs.push_back({"C1 <= 0", c, "C1 outside"});
    c = validD(); c.C0 = std::nan("");            negs.push_back({"NaN C0", c, "NaN"});
    c = validD(); c.C2 = std::nan("");            negs.push_back({"NaN C2", c, "NaN"});
    c = validD(); c.digitizer_timebase = 0.0;     negs.push_back({"timebase 0", c, "DigitizerTimebase"});
    for (const Neg& n : negs)
    {
      const std::string why = n.cal.unsupportedReason();
      if (!has(why, n.reason)) { std::fprintf(stderr, "%s must be rejected for it, got '%s'\n", n.what, why.c_str()); return 1; }
    }
  }
  // C4 != 0 is an additive mass offset on the quadratic root (2019 timsTOF Pro, firmware 6.0.110;
  // the PXD017703 highsens files carry C4 = -0.0905 and every frame references that row). Pinned
  // against the vendor library on the first and the last frame of one such file; dropping the
  // term is -53..-938 ppm, so the ablation guard below is what makes this case worth having.
  {
    diaspextract::TdfMzCalibration c4;
    c4.model_type = 1; c4.digitizer_timebase = 0.2; c4.digitizer_delay = 25224.4;
    c4.C0 = 325.88517924839397; c4.C1 = 152987.23707178448; c4.C2 = 0.0030691291814964453;
    c4.C4 = -0.09048112117134345; c4.T1_ref = 25.197341867575105; c4.dC1 = 16.5;
    if (!c4.isSupported()) { std::fprintf(stderr, "C4 != 0 must be ACCEPTED: %s\n", c4.unsupportedReason().c_str()); return 1; }
    const double t1[2] = {25.201724358440018, 25.20983751539588};          // frame 1, frame 68428
    const double tof[5] = {2000.0, 50000.0, 150000.0, 300000.0, 403000.0};
    const double vendor[2][5] = {{98.00225185761244, 186.40848708590778, 461.14679357131746, 1102.6974507710581, 1702.6599405622635},
                                 {98.00223875072383, 186.40846214421924, 461.14673185089845, 1102.6973031683212, 1702.6597126456059}};
    diaspextract::TdfMzCalibration no_c4 = c4; no_c4.C4 = 0.0;
    double worst_c4 = 0.0, worst_no_c4 = 0.0;
    for (int f = 0; f < 2; ++f)
      for (int i = 0; i < 5; ++i)
      {
        const double b4 = c4.frameFactor(t1[f]), got = c4.tofToMz(tof[i], b4);
        worst_c4 = std::fmax(worst_c4, std::fabs(got - vendor[f][i]) / vendor[f][i] * 1e6);
        worst_no_c4 = std::fmax(worst_no_c4, std::fabs(no_c4.tofToMz(tof[i], b4) - vendor[f][i]) / vendor[f][i] * 1e6);
        const double back = c4.mzToTof(got, b4);
        if (std::fabs(back - tof[i]) > 1e-6) { std::fprintf(stderr, "C4 round trip: %.4f -> %.4f\n", tof[i], back); return 1; }
      }
    if (worst_c4 > 1e-4) { std::fprintf(stderr, "C4 case: %.6f ppm from the vendor library\n", worst_c4); return 1; }
    if (worst_no_c4 < 50.0) { std::fprintf(stderr, "C4 ablation only %.3f ppm -- the case cannot catch a dropped C4\n", worst_no_c4); return 1; }
    std::printf("OK  C4 offset case: %.6f ppm vs vendor (10 probes, 2 frames); without C4 %.1f ppm\n", worst_c4, worst_no_c4);
  }
  // A NEGATIVE C2 is the same law (the golden cases above pin six PXD029836 runs against the vendor library), accepted
  // inside its domain: the rising branch, with the slope dt/du at least half its linear value, for every TOF index up to
  // 1e8, every m/z up to 1e5 and every frame temperature within 1000 K of T1_ref.
  {
    const Cal px = pxd1418();
    if (!px.isSupported()) { std::fprintf(stderr, "PXD029836's negative C2 must be ACCEPTED: %s\n", px.unsupportedReason().c_str()); return 1; }
    // the corner the check decides on: dC1 > 0, so the smallest frame factor is at T1_ref - 1000 K
    const double b_lo = px.frameFactor(px.T1_ref - Cal::kDomainT1Span), b_hi = px.frameFactor(px.T1_ref + Cal::kDomainT1Span);
    if (!(b_lo < px.frameFactor(px.T1_ref) && px.frameFactor(px.T1_ref) < b_hi)) { std::fprintf(stderr, "corner: frame factor order\n"); return 1; }
    const double t_max = Cal::kDomainMaxTof * px.digitizer_timebase + px.digitizer_delay;
    const double ratio = (b_lo * b_lo + 4.0 * px.C2 * (t_max - px.C0)) / (b_lo * b_lo);
    if (!(ratio > 0.98 && ratio < 0.99)) { std::fprintf(stderr, "PXD029836 corner disc / b^2 = %.6f, expected 0.987\n", ratio); return 1; }
    // over the whole domain the conversion is finite, rising and round-trips, at both temperature extremes
    double worst_rt = 0.0;
    for (double b : {b_lo, b_hi})
    {
      double prev = 0.0;
      for (double tof = 0.0; tof <= Cal::kDomainMaxTof; tof += Cal::kDomainMaxTof / 4096.0)
      {
        const double mz = px.tofToMz(tof, b);
        if (!(std::isfinite(mz) && mz > prev)) { std::fprintf(stderr, "negative C2: tof %.0f -> %.9g after %.9g\n", tof, mz, prev); return 1; }
        // bins: 3e-8 at TOF index 8.4e7, as for C2 >= 0. Checked per point: fmax drops a NaN operand, so a NaN inverse
        // would never reach worst_rt (review 2026-09-17)
        const double rt = std::fabs(px.mzToTof(mz, b) - tof);
        if (!(rt <= 1e-6)) { std::fprintf(stderr, "negative C2 round trip: tof %.0f -> m/z %.9g -> %.3g bins off\n", tof, mz, rt); return 1; }
        worst_rt = std::fmax(worst_rt, rt);
        prev = mz;
      }
      for (double mz : {50.0, 1700.0, Cal::kDomainMaxMz})   // the inverse stays on the rising branch up to 1e5
      {
        const double tof = px.mzToTof(mz, b);
        if (!(std::fabs(px.tofToMz(tof, b) - mz) / mz < 1e-12)) { std::fprintf(stderr, "negative C2: m/z %.1f does not round trip\n", mz); return 1; }
      }
    }
    if (!(worst_rt <= 1e-6)) { std::fprintf(stderr, "negative C2 round trip: %.3g bins\n", worst_rt); return 1; }
    // a stored -0.0 is zero: accepted, and the same conversion bit for bit
    Cal nz = validD(), z = validD(); nz.C2 = -0.0; z.C2 = 0.0;
    if (!nz.isSupported()) { std::fprintf(stderr, "C2 = -0.0 must be ACCEPTED: %s\n", nz.unsupportedReason().c_str()); return 1; }
    for (double tof : {0.0, 1000.0, 300000.0, 634072.0})
      if (nz.tofToMz(tof, nz.frameFactor(25.7)) != z.tofToMz(tof, z.frameFactor(25.7))) { std::fprintf(stderr, "C2 = -0.0 differs from 0.0\n"); return 1; }

    // frames: a negative C2 refuses a T1 outside the checked span (and a NaN); positive and zero C2 never refuse one
    for (double dt : {0.0, 999.9, -999.9, Cal::kDomainT1Span, -Cal::kDomainT1Span})
      if (!px.frameReason(px.T1_ref + dt).empty()) { std::fprintf(stderr, "T1_ref %+g K must be converted: %s\n", dt, px.frameReason(px.T1_ref + dt).c_str()); return 1; }
    for (double t1 : {px.T1_ref + 1000.001, px.T1_ref - 1000.001, 1e9, -1e9, std::nan("")})
      if (!has(px.frameReason(t1), "Frames.T1")) { std::fprintf(stderr, "T1 %g must be refused for a negative C2\n", t1); return 1; }
    for (double t1 : {1e9, -1e9, std::nan("")})
    {
      Cal zero = validD(); zero.C2 = 0.0;
      if (!validD().frameReason(t1).empty() || !zero.frameReason(t1).empty()) { std::fprintf(stderr, "C2 >= 0 must not refuse a frame\n"); return 1; }
    }

    // a negative C2 whose root stops existing INSIDE the domain is refused, and the reason names the TOF index where the
    // discriminant turns negative; the ablation shows the refusal is not vacuous: without it that TOF converts to NaN
    Cal deep = px; deep.C2 = -0.1;
    const std::string why = deep.unsupportedReason();
    const double k = -4.0 * deep.C2, bw = deep.frameFactor(deep.T1_ref - Cal::kDomainT1Span);
    const double tof_zero = (deep.C0 + bw * bw / k - deep.digitizer_delay) / deep.digitizer_timebase;
    if (!(tof_zero < Cal::kDomainMaxTof) || !has(why, "outside its verified domain") || !has(why, ("stops existing at TOF index " + g6(tof_zero)).c_str()))
    { std::fprintf(stderr, "C2 = -0.1 must be refused naming TOF index %s, got '%s'\n", g6(tof_zero).c_str(), why.c_str()); return 1; }
    if (!std::isnan(deep.tofToMz(Cal::kDomainMaxTof, bw)) || !std::isfinite(deep.tofToMz(0.99 * tof_zero, bw)))
    { std::fprintf(stderr, "C2 = -0.1: the discriminant must turn negative at TOF index %.6g, no earlier\n", tof_zero); return 1; }
    // the margin: C2 = -0.07 still has a root at every TOF index up to 1e8, but its slope falls below half there
    Cal shallow = px; shallow.C2 = -0.07;
    if (shallow.isSupported() || !std::isfinite(shallow.tofToMz(Cal::kDomainMaxTof, bw)) || !has(shallow.unsupportedReason(), "slope halves"))
    { std::fprintf(stderr, "C2 = -0.07 must be refused by the slope margin: '%s'\n", shallow.unsupportedReason().c_str()); return 1; }
    Cal inside = px; inside.C2 = -0.05;   // 48 times PXD029836's: still inside
    if (!inside.isSupported()) { std::fprintf(stderr, "C2 = -0.05 must be ACCEPTED: %s\n", inside.unsupportedReason().c_str()); return 1; }
    // the boundary itself, stated independently: the slope is half its linear value at TOF index 1e8 and T1_ref - 1000 K
    // (dC1 > 0) when |C2| = 3 b^2 / (16 (t_max - C0)); a wrong corner, span or margin moves it by far more than 1e-6
    const double c2_edge = 3.0 * b_lo * b_lo / (16.0 * (t_max - px.C0));
    Cal edge_in = px, edge_out = px; edge_in.C2 = -c2_edge * (1.0 - 1e-6); edge_out.C2 = -c2_edge * (1.0 + 1e-6);
    if (!edge_in.isSupported() || edge_out.isSupported())
    { std::fprintf(stderr, "the domain edge |C2| = %.9g: just inside '%s', just outside '%s'\n", c2_edge, edge_in.unsupportedReason().c_str(), edge_out.unsupportedReason().c_str()); return 1; }
    // the m/z side: a fine timebase keeps every TOF index early, yet the branch peaks below m/z 1e5
    Cal fine = px; fine.C2 = -2.5; fine.digitizer_timebase = 0.001;
    const double kf = -4.0 * fine.C2, bf = fine.frameFactor(fine.T1_ref - Cal::kDomainT1Span);
    const bool tof_side_ok = Cal::kDomainMaxTof * fine.digitizer_timebase + fine.digitizer_delay <= fine.C0 + 0.75 * bf * bf / kf;
    if (!tof_side_ok || fine.isSupported() || !has(fine.unsupportedReason(), "m/z 100000"))
    { std::fprintf(stderr, "C2 = -2.5 at 1 ps must be refused on the m/z side (TOF side ok %d): '%s'\n", int(tof_side_ok), fine.unsupportedReason().c_str()); return 1; }
    Cal minus_inf = px; minus_inf.C2 = -std::numeric_limits<double>::infinity();
    if (minus_inf.isSupported()) { std::fprintf(stderr, "C2 = -inf must be rejected\n"); return 1; }
    // a temperature coefficient that takes the C1 factor through zero within the span cannot be checked
    Cal hot = px; hot.dC1 = 1000.0;
    Cal hot_pos = validD(); hot_pos.dC1 = 1000.0;
    if (!has(hot.unsupportedReason(), "dC1") || !hot_pos.isSupported())
    { std::fprintf(stderr, "dC1 1000: refused for C2 < 0 ('%s'), accepted for C2 > 0 ('%s')\n", hot.unsupportedReason().c_str(), hot_pos.unsupportedReason().c_str()); return 1; }
    std::printf("OK  negative C2: PXD029836 accepted (corner disc/b^2 %.4f, round trip %.1e bins up to TOF index 1e8); C2 -0.1 refused at TOF index %s, "
                "-0.07 by the margin, -2.5 at 1 ps on m/z; edge |C2| %.6g; T1 outside +-1000 K refused\n", ratio, worst_rt, g6(tof_zero).c_str(), c2_edge);
  }
  // C2 == 0.0 STORED in the file is a real calibration -- the 2020 timsTOF Pro firmware behind
  // PXD017703 ships t = C0 + C1_eff*sqrt(m) with no quadratic term, every frame referencing it --
  // and must be ACCEPTED. A NULL C2 is what must be refused; the loader converts NULL to NaN, which
  // the "NaN must be rejected" case below covers. Pin the linear law in closed form and round trip.
  {
    diaspextract::TdfMzCalibration lin;
    lin.model_type = 1; lin.digitizer_timebase = 0.2; lin.digitizer_delay = 25131.0;
    lin.C0 = 315.70325869866065; lin.C1 = 154272.1271422364; lin.C2 = 0.0;
    lin.T1_ref = 25.63315397685876; lin.dC1 = -0.2;               // the PXD017703 row 1 values
    if (!lin.isSupported()) { std::fprintf(stderr, "C2 == 0 stored must be ACCEPTED: %s\n", lin.unsupportedReason().c_str()); return 1; }
    const double b = lin.frameFactor(lin.T1_ref);
    for (double tof : {60000.0, 120000.0, 240000.0, 400000.0})
    {
      const double t = tof * lin.digitizer_timebase + lin.digitizer_delay;
      const double u = (t - lin.C0) / b;                            // the quadratic collapses to this
      const double expect = u * u, got = lin.tofToMz(tof, b);
      if (!(std::fabs(got - expect) / expect < 1e-12)) { std::fprintf(stderr, "C2==0 closed form: tof %.0f got %.9f want %.9f\n", tof, got, expect); return 1; }
      const double back = lin.mzToTof(got, b);
      if (std::fabs(back - tof) > 1e-6) { std::fprintf(stderr, "C2==0 round trip: %.4f -> %.4f\n", tof, back); return 1; }
    }
  }

  // unphysical TOF must NOT return a plausible mass (silent-wrongness guard)
  diaspextract::TdfMzCalibration ok;
  ok.model_type = 1; ok.digitizer_timebase = 0.125; ok.digitizer_delay = 25655.375;
  ok.C0 = 279.3262846272992; ok.C1 = 155279.13067653627; ok.C2 = 0.001260061434461731;
  ok.T1_ref = 25.693668980735552; ok.dC1 = 20.0;
  if (!ok.isSupported()) { std::fprintf(stderr, "reference calibration must be supported\n"); return 1; }
  const double bb = ok.frameFactor(ok.T1_ref);
  // tof = 0 is still physical here: t = delay (25655) is far above C0 (279), and the model correctly
  // returns the acquisition floor, ~m/z 100 (MzAcqRangeLower). Assert that rather than assuming 0.
  const double mz_at_zero = ok.tofToMz(0.0, bb);
  if (!(mz_at_zero > 99.0 && mz_at_zero < 101.0))
  { std::fprintf(stderr, "tof=0 should give the ~100 m/z acquisition floor, got %.6f\n", mz_at_zero); return 1; }
  if (!(ok.tofToMz(1e5, bb) > mz_at_zero))
  { std::fprintf(stderr, "m/z must increase with tof\n"); return 1; }
  // the t <= C0 guard: force it with a calibration whose zero sits above the arrival time
  diaspextract::TdfMzCalibration late = ok; late.C0 = 1e9;
  if (!std::isnan(late.tofToMz(1e5, late.frameFactor(late.T1_ref))))
  { std::fprintf(stderr, "t <= C0 must yield NaN, not a plausible-looking mass\n"); return 1; }

  std::printf("OK  %d golden cases, %d files (%d with C2 < 0): max %.7f ppm; C2 ablation %.2f ppm, temperature %.3f ppm; "
              "round trip %.1e bins; negative paths pass\n", n_cases, n_files, n_negative, worst, worst_noc2, worst_notemp, worst_back);
  return 0;
}
