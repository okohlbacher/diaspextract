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

  size_t pos = 0;
  int n_files = 0, n_cases = 0;
  double worst = 0.0, worst_noc2 = 0.0;
  while (true)
  {
    const size_t fstart = j.find("\"model_type\":", pos);
    if (fstart == std::string::npos) break;
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
      const double b = cal.frameFactor(t1);
      const double got = cal.tofToMz(tof, b);
      const double ppm = std::fabs(got - mz) / mz * 1e6;
      if (ppm > 1e-4)   // same tolerance as the python twin; measured max is 2.6e-5
      {
        std::fprintf(stderr, "FAIL frame %.0f tof %.0f: got %.6f want %.6f (%.6f ppm)\n",
                     frame, tof, got, mz, ppm);
        return 1;
      }
      worst = std::fmax(worst, ppm);
      // round trip must return the same TOF
      const double back = cal.mzToTof(got, b);
      if (std::fabs(back - tof) > 1e-3)
      { std::fprintf(stderr, "FAIL round trip: tof %.4f -> mz -> %.4f\n", tof, back); return 1; }
      // ablation: dropping C2 must be clearly visible, else the golden set proves nothing
      diaspextract::TdfMzCalibration no_c2 = cal; no_c2.C2 = 0.0;
      worst_noc2 = std::fmax(worst_noc2, std::fabs(no_c2.tofToMz(tof, b) - mz) / mz * 1e6);
      ++n_cases;
      cp = q;
    }
    ++n_files;
    pos = fstart + 10;
  }
  if (n_cases < 30) { std::fprintf(stderr, "only %d golden cases parsed\n", n_cases); return 1; }
  if (worst_noc2 < 5.0)
  { std::fprintf(stderr, "C2 ablation only %.3f ppm -- golden set cannot catch the known-bad port\n", worst_noc2); return 1; }

  // negative paths must be REJECTED, not approximated
  diaspextract::TdfMzCalibration bad;
  bad.model_type = 2; bad.C1 = 1.0; bad.digitizer_timebase = 0.125;
  if (bad.isSupported()) { std::fprintf(stderr, "ModelType 2 must be rejected\n"); return 1; }
  bad.model_type = 1; bad.dC2 = 1e-9;
  if (bad.isSupported()) { std::fprintf(stderr, "dC2 != 0 must be rejected\n"); return 1; }
  bad.dC2 = 0.0; bad.C3 = 1e-9;
  if (bad.isSupported()) { std::fprintf(stderr, "C3 != 0 must be rejected\n"); return 1; }
  bad.C3 = 0.0; bad.C4 = 0.0;
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
  bad.C2 = -1e-6;
  if (bad.isSupported()) { std::fprintf(stderr, "C2 < 0 must be rejected\n"); return 1; }
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
  bad.C2 = 1e-3; bad.C1 = 1e-300;   // positive but implausible: b*b overflows -> silent m/z 0
  if (bad.isSupported()) { std::fprintf(stderr, "implausible C1 must be rejected\n"); return 1; }
  bad.C1 = 155279.0; bad.C0 = std::nan("");
  if (bad.isSupported()) { std::fprintf(stderr, "NaN must be rejected\n"); return 1; }
  bad.C0 = 279.3; bad.C1 = 0.0;
  if (bad.isSupported()) { std::fprintf(stderr, "C1 <= 0 must be rejected\n"); return 1; }

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

  std::printf("OK  %d golden cases, %d files: max %.6f ppm; C2 ablation %.2f ppm; round trip + negative paths pass\n",
              n_cases, n_files, worst, worst_noc2);
  return 0;
}
