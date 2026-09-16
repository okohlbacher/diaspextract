// Golden and rule tests for the dnoise v0.1.0 MS1 port (src/DnoiseMs1.h). Runs anywhere: no OpenMS, no vendor data,
// no cluster -- only a C++20 compiler. Build: c++ -std=c++20 -ffp-contract=off -I src tests/test_dnoise_ms1.cpp
//
// Bit-identity with the dnoise binary itself was established on real frames of both files below (SDK-read and
// blob-decoded, every point). What is pinned here is what that identity rests on: both files' gate boxes -- D's three
// FMA-sensitive edges and TNBC 009's "Bruker otofControl" widening included --, each streak, halo and gate rule on a
// synthetic frame, and the optimised implementation against a direct transliteration of the Rust on random frames.
// Also the tool's point recovery (recover::) against the loader's forward maths, and its frame bookkeeping.
#include "DnoiseMs1.h"
#include "TdfMzCalibration.h"   // the m/z model the patched loader calibrates with, which recover::recoverPoint inverts
#include <algorithm>
#include <array>
#include <cfloat>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace dn = spx::dnoise;
using u32 = std::uint32_t;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::fprintf(stderr, "FAIL line %d: ", __LINE__); std::fprintf(stderr, __VA_ARGS__); std::fputc('\n', stderr); ++g_fail; } } while (0)

// ---- the gate inputs of two real files, from read-only queries of their analysis.tdf -------------------------------
// Rows are SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth FROM DiaFrameMsMsWindows, as exact
// doubles. Boxes are {scan_lo, scan_hi, tof_lo, tof_hi}: D's reproduced the dnoise binary on real frames; both sets
// come from an independent CPython implementation of the spec.

static dn::RunMeta metaD()   // dataset D
{
  dn::RunMeta m;
  m.acquisition_software = "timsTOF";
  m.mz_acq_range_lower = "99.990834";
  m.mz_acq_range_upper = "1700.000000";
  m.digitizer_num_samples = "634073";
  m.one_over_k0_acq_range_lower = "0.600000";
  m.one_over_k0_acq_range_upper = "1.400000";
  m.max_num_scans = 944;
  return m;
}
static const std::vector<dn::WindowRow> kRowsD = {
  {34, 602, 0x1.676c28f5c28f6p+9, 0x1.70cccccccccc0p+4},  {602, 944, 0x1.8e28f5c28f5c2p+8, 0x1.1a8f5c28f5c28p+7},
  {34, 579, 0x1.734f5c28f5c29p+9, 0x1.a800000000000p+4},  {579, 944, 0x1.e431eb851eb85p+8, 0x1.0651eb851eb80p+5},
  {34, 568, 0x1.7ff1eb851eb85p+9, 0x1.a0a3d70a3d700p+4},  {568, 944, 0x1.fff1eb851eb84p+8, 0x1.8b5c28f5c2900p+4},
  {34, 557, 0x1.8d151eb851eb8p+9, 0x1.c828f5c28f5c0p+4},  {557, 944, 0x1.0b45c28f5c28fp+9, 0x1.67d70a3d70a40p+4},
  {34, 545, 0x1.9ab6666666666p+9, 0x1.c028f5c28f5c0p+4},  {545, 944, 0x1.155ae147ae148p+9, 0x1.3d70a3d70a400p+4},
  {34, 523, 0x1.a975c28f5c28fp+9, 0x1.07d70a3d70a40p+5},  {523, 944, 0x1.1e9eb851eb852p+9, 0x1.33851eb851ec0p+4},
  {34, 511, 0x1.b9d7ae147ae14p+9, 0x1.1466666666680p+5},  {511, 944, 0x1.282b851eb851fp+9, 0x1.4fae147ae1480p+4},
  {34, 488, 0x1.cc38f5c28f5c2p+9, 0x1.47c28f5c28f40p+5},  {488, 944, 0x1.320c28f5c28f6p+9, 0x1.487ae147ae140p+4},
  {34, 477, 0x1.e1bccccccccccp+9, 0x1.78b851eb85200p+5},  {477, 944, 0x1.3bfb851eb851fp+9, 0x1.535c28f5c2900p+4},
  {34, 443, 0x1.fc228f5c28f5cp+9, 0x1.e400000000000p+5},  {443, 944, 0x1.463999999999ap+9, 0x1.5c28f5c28f5c0p+4},
  {34, 409, 0x1.10928f5c28f5cp+10, 0x1.6628f5c28f5c0p+6}, {409, 944, 0x1.50e7ae147ae14p+9, 0x1.6f5c28f5c2900p+4},
  {34, 295, 0x1.3cd5c28f5c28fp+10, 0x1.0a8f5c28f5c28p+8}, {295, 944, 0x1.5c270a3d70a3ep+9, 0x1.807ae147ae180p+4},
};
static const std::vector<dn::GateBox> kBoxesD = {
  {0, 661, 335024, 347539}, {542, 943, 161592, 238906}, {0, 638, 343405, 357004}, {519, 943, 233745, 253491},
  {0, 627, 352940, 366144}, {508, 943, 248497, 264073}, {0, 616, 362146, 376020}, {497, 943, 259193, 273460},
  {0, 605, 372090, 385555}, {485, 943, 268677, 281541}, {0, 583, 381689, 396648}, {464, 943, 276838, 289229},
  {0, 571, 392855, 408070}, {451, 943, 284601, 297525}, {0, 547, 404348, 421405}, {428, 943, 292974, 305502},
  {0, 537, 417763, 436436}, {418, 943, 301024, 313628}, {0, 503, 432880, 455334}, {384, 943, 309220, 321844},
  {0, 468, 451880, 482490}, {349, 943, 317506, 330397}, {0, 355, 479174, 558154}, {235, 943, 326129, 339222},
};

static dn::RunMeta metaTnbc()   // PXD047793, 4223_TIMS2_009_2_Slot2-9_1_2883 (TNBC 009)
{
  dn::RunMeta m;
  m.acquisition_software = "Bruker otofControl";
  m.mz_acq_range_lower = "99.994278";
  m.mz_acq_range_upper = "1700.000000";
  m.digitizer_num_samples = "399967";
  m.one_over_k0_acq_range_lower = "0.602016";
  m.one_over_k0_acq_range_upper = "1.601104";
  m.max_num_scans = 918;
  return m;
}
static const std::vector<dn::WindowRow> kRowsTnbc = {
  {0, 665, 763.0, 26.0},  {666, 918, 413.0, 26.0}, {0, 647, 788.0, 26.0},  {648, 918, 438.0, 26.0},
  {0, 628, 813.0, 26.0},  {629, 918, 463.0, 26.0}, {0, 615, 838.0, 26.0},  {616, 918, 488.0, 26.0},
  {0, 596, 863.0, 26.0},  {597, 918, 513.0, 26.0}, {0, 578, 888.0, 26.0},  {579, 918, 538.0, 26.0},
  {0, 559, 913.0, 26.0},  {560, 918, 563.0, 26.0}, {0, 541, 938.0, 26.0},  {542, 918, 588.0, 26.0},
  {0, 523, 963.0, 26.0},  {524, 918, 613.0, 26.0}, {0, 504, 988.0, 26.0},  {505, 918, 638.0, 26.0},
  {0, 486, 1013.0, 26.0}, {487, 918, 663.0, 26.0}, {0, 472, 1038.0, 26.0}, {473, 918, 688.0, 26.0},
  {0, 454, 1063.0, 26.0}, {455, 918, 713.0, 26.0}, {0, 436, 1088.0, 26.0}, {437, 918, 738.0, 26.0},
};
static const std::vector<dn::GateBox> kBoxesTnbc = {
  {0, 711, 222496, 230760}, {620, 917, 128416, 139650}, {0, 693, 228255, 236387}, {602, 917, 136268, 147176},
  {0, 674, 233921, 241927}, {583, 917, 143890, 154499}, {0, 661, 239499, 247384}, {570, 917, 151300, 161634},
  {0, 642, 244992, 252762}, {551, 917, 158516, 168595}, {0, 624, 250404, 258064}, {533, 917, 165552, 175394},
  {0, 605, 255739, 263294}, {514, 917, 172420, 182041}, {0, 587, 261001, 268454}, {496, 917, 179133, 188547},
  {0, 569, 266191, 273547}, {478, 917, 185701, 194920}, {0, 550, 271313, 278575}, {459, 917, 192131, 201169},
  {0, 532, 276369, 283541}, {441, 917, 198433, 207299}, {0, 518, 281363, 288448}, {427, 917, 204615, 213317},
  {0, 500, 286295, 293296}, {409, 917, 210681, 219230}, {0, 482, 291169, 298089}, {391, 917, 216640, 225043},
};

static bool same(const dn::GateBox& a, const dn::GateBox& b)
{
  return a.scan_lo == b.scan_lo && a.scan_hi == b.scan_hi && a.tof_lo == b.tof_lo && a.tof_hi == b.tof_hi;
}

static void checkBoxes(const char* name, const dn::RunMeta& meta, const std::vector<dn::WindowRow>& rows,
                       const std::vector<dn::GateBox>& want)
{
  const dn::GateBuild b = dn::buildMs1WindowGate(meta, rows, dn::Params{});
  CHECK(b.error.empty() && b.gate.has_value(), "%s: no gate (%s)", name, b.error.c_str());
  if (!b.gate) return;
  const std::vector<dn::GateBox>& got = b.gate->boxes();
  CHECK(got.size() == want.size(), "%s: %zu boxes, want %zu", name, got.size(), want.size());
  for (std::size_t i = 0; i < std::min(got.size(), want.size()); ++i)
    CHECK(same(got[i], want[i]), "%s box %zu: scans [%u,%u] tof [%u,%u], want [%u,%u] [%u,%u]", name, i + 1, got[i].scan_lo,
          got[i].scan_hi, got[i].tof_lo, got[i].tof_hi, want[i].scan_lo, want[i].scan_hi, want[i].tof_lo, want[i].tof_hi);
}

static void testGate()
{
  checkBoxes("D", metaD(), kRowsD, kBoxesD);
  checkBoxes("TNBC 009", metaTnbc(), kRowsTnbc, kBoxesTnbc);

  // D boxes 11, 15 and 16 sit on razor edges (s1 = 582.0000000000001 and 547.0, s0 = 428.9999999999999): with
  // `intercept + slope * scan` fused into one FMA they become [0,582], [0,548] and [429,943]. Named for the report.
  const dn::GateBuild d = dn::buildMs1WindowGate(metaD(), kRowsD, dn::Params{});
  if (d.gate && d.gate->boxes().size() == 24)
  {
    const std::vector<dn::GateBox>& b = d.gate->boxes();
    CHECK(b[10].scan_hi == 583 && b[14].scan_hi == 547 && b[15].scan_lo == 428,
          "D boxes 11/15/16: scan_hi %u (583), scan_hi %u (547), scan_lo %u (428) -- floating-point contraction?",
          b[10].scan_hi, b[14].scan_hi, b[15].scan_lo);
  }
  // The goldens can tell the quirks apart: honouring the exclusive ScanNumEnd, or dropping the otofControl widening,
  // changes the boxes.
  std::vector<dn::WindowRow> exclusive = kRowsD;
  for (dn::WindowRow& r : exclusive) --r.scan_num_end;
  const dn::GateBuild e = dn::buildMs1WindowGate(metaD(), exclusive, dn::Params{});
  std::size_t changed = 0;
  for (std::size_t i = 0; e.gate && i < e.gate->boxes().size(); ++i) changed += !same(e.gate->boxes()[i], kBoxesD[i]);
  CHECK(changed > 0, "an exclusive ScanNumEnd must change D's boxes");
  dn::RunMeta plain = metaTnbc();
  plain.acquisition_software = "timsTOF";
  const dn::GateBuild t = dn::buildMs1WindowGate(plain, kRowsTnbc, dn::Params{});
  std::size_t unchanged = 0;
  for (std::size_t i = 0; t.gate && i < t.gate->boxes().size(); ++i)
    unchanged += t.gate->boxes()[i].tof_lo == kBoxesTnbc[i].tof_lo || t.gate->boxes()[i].tof_hi == kBoxesTnbc[i].tof_hi;
  CHECK(t.gate && unchanged == 0, "without the otofControl rule every TNBC TOF edge must move (%zu did not)", unchanged);

  // dia_ms1.rs semantics of the lookup.
  const auto g = dn::Ms1WindowGate::fromBoxes({{10, 20, 1000, 2000}}, 100);
  CHECK(g && g->contains(15, 1500) && g->contains(10, 1000) && g->contains(20, 2000), "box corners are inclusive");
  CHECK(g && !g->contains(15, 999) && !g->contains(15, 2001) && !g->contains(9, 1500) && !g->contains(21, 1500),
        "one index past a corner is outside");
  CHECK(!dn::Ms1WindowGate::fromBoxes({}, 100) && !dn::Ms1WindowGate::fromBoxes({{0, 1, 2, 3}}, 0), "no boxes or no scans: no gate");
  const auto m = dn::Ms1WindowGate::fromBoxes({{10, 20, 1000, 1500}, {12, 18, 1501, 2000}, {0, 0, 100, 200}, {0, 0, 202, 300}}, 100);
  CHECK(m && m->contains(15, 1500) && m->contains(15, 1501) && !m->contains(11, 1501), "touching intervals merge, per scan");
  CHECK(m && m->contains(0, 200) && !m->contains(0, 201) && m->contains(0, 202), "a one-index gap stays a gap");
  const auto c = dn::Ms1WindowGate::fromBoxes({{150, 200, 10, 20}, {0, 5, 30, 25}}, 100);
  CHECK(c && c->contains(99, 15) && !c->contains(98, 15) && !c->contains(100, 15),
        "a box past the last scan is clamped onto it; a scan past the last is outside");
  CHECK(c && !c->contains(0, 27) && !c->contains(0, 30), "a box with tof_hi < tof_lo is skipped");

  // When the builder gives no gate, and when it refuses the run.
  dn::Params off;
  off.dia_ms1_window = false;
  const dn::GateBuild a = dn::buildMs1WindowGate(metaD(), kRowsD, off);
  CHECK(!a.gate && a.error.empty(), "dia_ms1_window false: no gate, no error");
  const dn::GateBuild n = dn::buildMs1WindowGate(dn::RunMeta{}, {{10, 10, 500.0, 20.0}, {20, 5, 500.0, 20.0}}, dn::Params{});
  CHECK(!n.gate && n.error.empty(), "no row with ScanNumEnd > ScanNumBegin: no gate, and the metadata is never read");
  dn::RunMeta bad = metaD();
  bad.digitizer_num_samples.reset();
  const dn::GateBuild x = dn::buildMs1WindowGate(bad, kRowsD, dn::Params{});
  CHECK(!x.gate && !x.error.empty(), "a missing DigitizerNumSamples refuses the run");
  bad = metaD();
  bad.one_over_k0_acq_range_lower = "0,600000";
  const dn::GateBuild y = dn::buildMs1WindowGate(bad, kRowsD, dn::Params{});
  CHECK(!y.gate && !y.error.empty(), "a decimal comma refuses the run");
  bad = metaD();
  bad.one_over_k0_acq_range_lower = "1e-999";   // Rust parses 0 here and dnoise builds a degenerate gate; the port refuses
  const dn::GateBuild u = dn::buildMs1WindowGate(bad, kRowsD, dn::Params{});
  CHECK(!u.gate && u.error.find("OneOverK0AcqRangeLower '1e-999'") != std::string::npos, "an underflowing value refuses the run (%s)", u.error.c_str());
  bad = metaD();
  bad.max_num_scans = 0;
  const dn::GateBuild z = dn::buildMs1WindowGate(bad, kRowsD, dn::Params{});
  CHECK(!z.gate && z.error.empty(), "no scans: no gate");
}

static void testParse()
{
  const std::pair<const char*, double> accept[] = {
    {"99.990834", 99.990834}, {"1700.000000", 1700.0}, {"0.602016", 0.602016}, {"+1.5", 1.5}, {"-0.25", -0.25},
    {"1.", 1.0}, {".5", 0.5}, {"1e3", 1000.0}, {"2.5E-1", 0.25}, {"1e+2", 100.0}, {"0.1", 0.1},
    {"9007199254740993", 9007199254740992.0},                                             // halfway: ties to even
    {"1.00000000000000011102230246251565404236316680908203125", 1.0},                     // exactly halfway to the next double
    {"1.00000000000000011102230246251565404236316680908203126", 0x1.0000000000001p+0},    // just above: rounds up
    {"0e999", 0.0}, {"-0", -0.0},                                                         // zero, in any spelling
    {"2.2250738585072014e-308", DBL_MIN}, {"1.7976931348623157e308", DBL_MAX},            // the domain's two ends
  };
  for (const auto& [s, want] : accept)
  {
    double v = -7.0;
    CHECK(dn::detail::parseF64(s, v) && v == want, "parseF64(\"%s\") = %.17g, want %.17g", s, v, want);
  }
  double neg0 = 1.0;
  CHECK(dn::detail::parseF64("-0", neg0) && std::signbit(neg0), "parseF64(\"-0\") keeps the sign");
  // The last seven are outside the finite-decimal domain, which Rust accepts: refused, and alike on every library.
  for (const char* s : {"", " 1.5", "1.5 ", "1,5", ".", "e5", "1e", "1e+", "+-1", "--1", "0x1p3", "inf", "nan",
                        "infinity", "1e999", "1.5f", "1.2.3",
                        "INF", "-nan", "Infinity", "1e-999", "-1e-999", "1e-320", "2.2250738585072011e-308"})
  {
    double v = 0.0;
    CHECK(!dn::detail::parseF64(s, v), "parseF64(\"%s\") must fail", s);
  }
  u32 w = 0;
  CHECK(dn::detail::parseU32("634073", w) && w == 634073u, "parseU32 634073");
  CHECK(dn::detail::parseU32("+7", w) && w == 7u && dn::detail::parseU32("4294967295", w) && w == 4294967295u, "parseU32 + and max");
  for (const char* s : {"", "+", "-1", "4294967296", "1.0", " 1", "12a", "++1"})
    CHECK(!dn::detail::parseU32(s, w), "parseU32(\"%s\") must fail", s);
  CHECK(dn::detail::asU32(std::nan("")) == 0 && dn::detail::asU32(-5.0) == 0 && dn::detail::asU32(0.9) == 0 &&
          dn::detail::asU32(3.0) == 3 && dn::detail::asU32(4294967295.0) == 4294967295u &&
          dn::detail::asU32(4294967296.0) == 4294967295u && dn::detail::asU32(1e300) == 4294967295u &&
          dn::detail::asU32(std::numeric_limits<double>::infinity()) == 4294967295u,
        "asU32 is Rust's saturating `as u32`");
}

// ---- frames -----------------------------------------------------------------------------------------------------------
struct Pt
{
  u32 scan, tof, inten;
};
struct Run
{
  std::vector<std::uint8_t> keep;
  dn::FrameCounts counts;
};

static dn::Scratch g_scratch;   // shared by every call on purpose: working memory must never leak from one frame to the next

static Run run(const std::vector<Pt>& f, std::size_t num_scans, const dn::Params& p, const dn::Ms1WindowGate* gate = nullptr)
{
  std::vector<u32> s(f.size()), t(f.size()), v(f.size());
  for (std::size_t i = 0; i < f.size(); ++i)
  {
    s[i] = f[i].scan; t[i] = f[i].tof; v[i] = f[i].inten;
  }
  Run r;
  r.counts = dn::denoiseMs1Frame(s.data(), t.data(), v.data(), f.size(), num_scans, p, gate, g_scratch, r.keep);
  return r;
}

static void column(std::vector<Pt>& f, u32 tof, std::initializer_list<u32> scans, u32 inten = 100)
{
  for (u32 s : scans) f.push_back({s, tof, inten});
}

static std::size_t kept(const std::vector<std::uint8_t>& k, std::size_t from, std::size_t to)
{
  return static_cast<std::size_t>(std::count(k.begin() + static_cast<std::ptrdiff_t>(from), k.begin() + static_cast<std::ptrdiff_t>(to), 1));
}

static void testStreak()
{
  dn::Params p;
  p.halo = false;
  std::vector<Pt> f;
  column(f, 1000, {10, 11, 12, 13, 14});   // [0,5)   5 occupied scans: kept
  column(f, 2000, {10, 11, 12, 13});       // [5,9)   4: dropped
  column(f, 3000, {10, 12, 14, 16});       // [9,13)  4 occupied over a span of 7: dropped, the span does not count
  column(f, 4000, {10, 12, 14, 16, 18});   // [13,18) 5 occupied, one empty scan between each: kept
  column(f, 5000, {10, 11, 12, 15, 16});   // [18,23) 2 empty scans (13, 14) are bridged: one run of 5, kept
  column(f, 6000, {10, 11, 12, 16, 17});   // [23,28) 3 empty scans split it into 3 + 2: dropped
  column(f, 7000, {10, 11, 12});           // [28,31) 7003 is inside [7000-3, 7000+3] and 7000 inside 7003's window:
  column(f, 7003, {13, 14});               // [31,33) both see 5 occupied scans
  column(f, 8000, {10, 11, 12});           // [33,36) 8004 is one index too far
  column(f, 8004, {13, 14});               // [36,38)
  Run r = run(f, 100, p);
  CHECK(kept(r.keep, 0, 5) == 5 && kept(r.keep, 5, 9) == 0, "min_feature_length 5");
  CHECK(kept(r.keep, 9, 13) == 0 && kept(r.keep, 13, 18) == 5, "a run's length is its occupied scans, not its span");
  CHECK(kept(r.keep, 18, 23) == 5 && kept(r.keep, 23, 28) == 0, "max_internal_gap 2 bridges 2 empty scans, not 3");
  CHECK(kept(r.keep, 28, 33) == 5 && kept(r.keep, 33, 38) == 0, "the window is [c-3, c+3], inclusive");

  // Pass 2 sums only pass 1's survivors: 9003 survives pass 1 on the strength of 9000 and 9006, which do not.
  std::vector<Pt> g;
  column(g, 9000, {10, 11});
  column(g, 9003, {12});
  column(g, 9006, {13, 14});
  p.iterations = 1;
  CHECK(run(g, 100, p).keep == std::vector<std::uint8_t>({0, 0, 1, 0, 0}), "one pass keeps the bridge only");
  p.iterations = 2;
  CHECK(kept(run(g, 100, p).keep, 0, 5) == 0, "the second pass drops it");
  p.iterations = 1000;
  CHECK(kept(run(g, 100, p).keep, 0, 5) == 0, "iterations 1000");
  p.iterations = 0;
  CHECK(kept(run(g, 100, p).keep, 0, 5) == 5, "iterations 0: the streak filter keeps every point");

  p = dn::Params{};
  p.halo = false;
  std::vector<Pt> z;
  column(z, 1000, {10, 11, 12, 13});
  z.push_back({14, 1000, 0});              // [0,5)  a zero-intensity point leaves its scan empty: 4 occupied, dropped
  column(z, 2000, {10, 11, 13, 14, 15});
  z.push_back({12, 2000, 0});              // [5,11) inside the accepted span [10,15]: kept
  z.push_back({20, 2000, 0});              // [11]   outside it: dropped
  r = run(z, 100, p);
  CHECK(kept(r.keep, 0, 5) == 0 && kept(r.keep, 5, 11) == 6 && r.keep[11] == 0,
        "zero-intensity points are empty scans, yet kept inside a span");

  std::vector<Pt> q;
  column(q, 3000, {10, 11, 13, 14});
  q.push_back({12, 3000, 10});             // below min_window_intensity 50: not occupied, but inside the span
  p.min_window_intensity = 50;
  CHECK(kept(run(q, 100, p).keep, 0, 5) == 0, "a scan below min_window_intensity is not occupied (4 < 5)");
  p.min_feature_length = 4;
  CHECK(kept(run(q, 100, p).keep, 0, 5) == 5, "the sub-floor point lies inside the accepted span");
  p.min_feature_intensity = 410;
  CHECK(kept(run(q, 100, p).keep, 0, 5) == 5, "the span total includes the sub-floor scan (4 x 100 + 10)");
  p.min_feature_intensity = 411;
  CHECK(kept(run(q, 100, p).keep, 0, 5) == 0, "min_feature_intensity 411 > 410");

  p = dn::Params{};
  p.halo = false;
  const u32 top = std::numeric_limits<u32>::max();
  std::vector<Pt> e;
  column(e, 0, {10, 11});
  column(e, 2, {12, 13, 14});              // the window of TOF 0 is [0, 3], not a wrapped-around one
  column(e, top, {10, 11});
  column(e, top - 2, {12, 13, 14});        // and the last TOF index's is [top - 3, top]
  CHECK(kept(run(e, 100, p).keep, 0, 10) == 10, "windows saturate at both ends of the TOF axis");

  bool threw = false;
  try { run({{100, 5, 1}}, 100, p); } catch (const std::invalid_argument&) { threw = true; }
  CHECK(threw, "a scan index >= num_scans throws");
  // Oversized frames are refused before any allocation or read: the arrays below hold ONE point, and with the streak off
  // the halo would index num_scans + 1 rows (a wrap at SIZE_MAX, scan + 1 at 2^32).
  const u32 one_scan = 0, one_tof = 5, one_inten = 1;
  dn::Params halo_only;
  halo_only.iterations = 0;
  auto refused = [&](std::size_t n, std::size_t num_scans) {
    std::vector<std::uint8_t> k = {7};
    try { dn::denoiseMs1Frame(&one_scan, &one_tof, &one_inten, n, num_scans, halo_only, nullptr, g_scratch, k); }
    catch (const std::invalid_argument&) { return k == std::vector<std::uint8_t>{7}; }
    return false; };
  CHECK(refused(1, std::numeric_limits<std::size_t>::max()), "num_scans = SIZE_MAX is refused before any allocation");
  CHECK(refused(1, std::size_t(1) << 32), "num_scans = 2^32 is refused before any allocation");
  CHECK(refused(std::size_t(std::numeric_limits<u32>::max()) + 1, 100), "2^32 points are refused before keep or an input point is touched");
  const Run empty = run({}, 100, p);
  CHECK(empty.keep.empty() && empty.counts.raw == 0 && empty.counts.kept == 0, "an empty frame");
}

static void testHalo()
{
  dn::Params p;
  p.iterations = 0;   // the streak keeps everything, so these frames exercise the halo alone
  const u32 top = std::numeric_limits<u32>::max();
  const std::vector<Pt> f = {
    {10, 1000, 10000}, {11, 1000, 5}, {11, 1001, 5},                  // [0,3)   the own column is skipped, the next one is not
    {20, 2000, 100}, {20, 2050, 15}, {21, 2060, 14},                   // [3,6)   15 >= 0.15 * 100 is kept (a tie), 14 is not
    {30, 3000, 1}, {40, 4000, 1}, {41, 4000, 1000},                    // [6,9)   an empty box, or only the own column: kept
    {50, 5000, 100000}, {53, 5010, 1}, {50, 5081, 1}, {52, 4920, 1},   // [9,13)  the box is scans +-2 x TOF +-80, inclusive
    {60, 6000, 10}, {60, 6010, 100}, {61, 5920, 1},                    // [13,16) simultaneous: 6000 is dropped, yet drops 5920
    {0, 8000, 1}, {1, 8001, 1000}, {99, 8100, 1}, {97, 8101, 1000},    // [16,20) the scan box is clamped at both ends
    {20, top, 100000}, {20, top - 80, 1}, {30, 0, 100000}, {30, 80, 1} // [20,24) the TOF box saturates
  };
  const std::vector<std::uint8_t> want = {1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0};
  const Run r = run(f, 100, p);
  for (std::size_t i = 0; i < want.size() && i < r.keep.size(); ++i)
    CHECK(r.keep[i] == want[i], "halo point %zu (scan %u, tof %u, I %u): kept %d, want %d", i, f[i].scan, f[i].tof,
          f[i].inten, r.keep[i], want[i]);

  // The reference set is the streak survivors: the bright point is alone in its window, so the streak removes it,
  // and with it its halo. Measured against raw points, the column would go (10 < 0.15 * 100000).
  std::vector<Pt> g;
  column(g, 7000, {10, 11, 12, 13, 14}, 10);
  g.push_back({12, 7010, 100000});
  const Run s = run(g, 100, dn::Params{});
  CHECK(kept(s.keep, 0, 5) == 5 && s.keep[5] == 0 && s.counts.after_streak == 5 && s.counts.after_halo == 5,
        "the halo's reference set is the streak survivors, not the raw points");

  // For fraction 0.15 the double compare equals the integer test 100 * I >= 15 * best, so either reading is exact.
  std::size_t disagree = 0;
  auto around = [&disagree](std::uint64_t best) {
    const std::uint64_t k = 15 * best / 100;
    for (std::uint64_t i = k ? k - 1 : 0; i <= k + 1 && i <= 0xFFFFFFFFu; ++i)
      disagree += (static_cast<double>(i) >= 0.15 * static_cast<double>(best)) != (100 * i >= 15 * best);
  };
  for (std::uint64_t best = 0; best <= 200000; ++best) around(best);
  for (std::uint64_t best = 0xFFFFFFFFu; best > 200000; best = best * 7 / 9) around(best);
  CHECK(disagree == 0, "0.15 double compare vs 100*I >= 15*best: %zu disagreements", disagree);
}

static void testGateAnd()
{
  // ANDed last: scans 10 and 11 are outside the gate yet still carry the streak of scans 12-14, which on its own
  // (3 occupied scans) would not survive -- gating first would drop all five.
  const auto gate = dn::Ms1WindowGate::fromBoxes({{12, 50, 8990, 9000}}, 100);
  std::vector<Pt> f;
  column(f, 9000, {10, 11, 12, 13, 14});
  column(f, 9500, {10, 11, 12, 13, 14});
  const Run r = run(f, 100, dn::Params{}, gate ? &*gate : nullptr);
  CHECK(r.keep == std::vector<std::uint8_t>({0, 0, 1, 1, 1, 0, 0, 0, 0, 0}), "the gate is ANDed after streak and halo");
  CHECK(r.counts.raw == 10 && r.counts.after_streak == 10 && r.counts.after_halo == 10 && r.counts.kept == 3, "stage counts");
  dn::Params p0;
  p0.iterations = 0;
  const Run r0 = run(f, 100, p0, gate ? &*gate : nullptr);
  CHECK(r0.counts.after_halo == 10 && r0.counts.kept == 3, "iterations 0 still runs the halo and the gate");
}

// ---- a direct transliteration of dnoise (filter.rs, halo.rs, writer.rs order) -----------------------------------------
// Sort plus binary search per window, and a touched list that may hold duplicates, exactly as the Rust does -- unlike
// the header's two-pointer sweep, early outs and radix sorts, so the random frames compare two implementations.
static std::vector<std::uint8_t> refFilterOnce(const std::vector<Pt>& f, std::size_t num_scans, const dn::Params& p)
{
  const std::size_t n = f.size();
  std::vector<std::uint8_t> keep(n, 0);
  if (n == 0) return keep;
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&f](std::size_t a, std::size_t b) { return f[a].tof < f[b].tof; });
  std::vector<u32> st(n);
  for (std::size_t k = 0; k < n; ++k) st[k] = f[order[k]].tof;
  std::vector<std::uint64_t> profile(num_scans, 0);
  std::vector<std::size_t> touched, occupied;
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  const u32 w = p.mz_half_width, top = std::numeric_limits<u32>::max();
  for (std::size_t k = 0; k < n;)
  {
    const u32 c = st[k];
    std::size_t block_hi = k + 1;
    while (block_hi < n && st[block_hi] == c) ++block_hi;
    const u32 lo_val = c >= w ? c - w : 0, hi_val = c <= top - w ? c + w : top;
    const std::size_t w_lo = static_cast<std::size_t>(std::lower_bound(st.begin(), st.end(), lo_val) - st.begin());
    const std::size_t w_hi = static_cast<std::size_t>(std::upper_bound(st.begin(), st.end(), hi_val) - st.begin());
    for (std::size_t j = w_lo; j < w_hi; ++j)
    {
      const std::size_t s = f[order[j]].scan;
      if (profile[s] == 0) touched.push_back(s);
      profile[s] += f[order[j]].inten;
    }
    occupied.clear();
    for (std::size_t s : touched)
      if (profile[s] > 0 && profile[s] >= p.min_window_intensity) occupied.push_back(s);
    std::sort(occupied.begin(), occupied.end());
    spans.clear();
    if (!occupied.empty())
    {
      auto close = [&](std::size_t rs, std::size_t re) {
        std::size_t count = 0;
        std::uint64_t total = 0;
        for (std::size_t s = rs; s <= re; ++s)
        {
          count += profile[s] > 0 && profile[s] >= p.min_window_intensity;
          total += profile[s];
        }
        if (count >= p.min_feature_length && total >= p.min_feature_intensity) spans.emplace_back(rs, re);
      };
      std::size_t run_start = occupied[0], prev = occupied[0];
      for (std::size_t q = 1; q < occupied.size(); ++q)
      {
        if (occupied[q] - prev > p.max_internal_gap + 1)
        {
          close(run_start, prev);
          run_start = occupied[q];
        }
        prev = occupied[q];
      }
      close(run_start, prev);
    }
    for (std::size_t j = k; j < block_hi; ++j)
    {
      const std::size_t s = f[order[j]].scan;
      std::size_t i = 0;
      while (i < spans.size() && spans[i].second < s) ++i;
      if (i < spans.size() && s >= spans[i].first) keep[order[j]] = 1;
    }
    for (std::size_t s : touched) profile[s] = 0;
    touched.clear();
    k = block_hi;
  }
  return keep;
}

static std::vector<std::uint8_t> refHalo(const std::vector<Pt>& f, std::size_t num_scans, const dn::Params& p)
{
  const std::size_t n = f.size();
  std::vector<std::uint8_t> keep(n, 1);
  if (n == 0 || num_scans == 0) return keep;
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(),
            [&f](std::size_t a, std::size_t b) { return std::tie(f[a].scan, f[a].tof) < std::tie(f[b].scan, f[b].tof); });
  std::vector<std::size_t> block_start(num_scans, 0), block_len(num_scans, 0);
  for (std::size_t k = 0; k < n;)
  {
    const std::size_t s = f[order[k]].scan, start = k;
    while (k < n && f[order[k]].scan == s) ++k;
    block_start[s] = start;
    block_len[s] = k - start;
  }
  const std::int64_t sw = static_cast<std::int64_t>(p.halo_scan_half_width);
  const u32 mw = p.halo_mz_idx_half_width, top = std::numeric_limits<u32>::max();
  for (std::size_t i = 0; i < n; ++i)
  {
    const std::int64_t si = f[order[i]].scan;
    const u32 mi = f[order[i]].tof;
    const std::size_t lo_scan = static_cast<std::size_t>(std::max<std::int64_t>(si - sw, 0));
    const std::size_t hi_scan = static_cast<std::size_t>(std::min<std::int64_t>(si + sw, static_cast<std::int64_t>(num_scans) - 1));
    const u32 lo_mz = mi >= mw ? mi - mw : 0, hi_mz = mi <= top - mw ? mi + mw : top;
    u32 best = 0;
    for (std::size_t v = lo_scan; v <= hi_scan; ++v)
    {
      if (block_len[v] == 0) continue;
      std::size_t left = block_start[v], right = block_start[v] + block_len[v];
      const std::size_t be = right;
      while (left < right)
      {
        const std::size_t mid = (left + right) / 2;
        if (f[order[mid]].tof < lo_mz) left = mid + 1;
        else right = mid;
      }
      for (std::size_t j = left; j < be && f[order[j]].tof <= hi_mz; ++j)
        if (f[order[j]].tof != mi && f[order[j]].inten > best) best = f[order[j]].inten;
    }
    keep[order[i]] = static_cast<double>(f[order[i]].inten) >= p.halo_peak_fraction * static_cast<double>(best);
  }
  return keep;
}

struct RefResult
{
  std::vector<std::uint8_t> keep;
  std::size_t streak = 0, halo = 0;
};

static RefResult refFrame(const std::vector<Pt>& f, std::size_t num_scans, const dn::Params& p, const dn::Ms1WindowGate* gate)
{
  const std::size_t n = f.size();
  RefResult r;
  r.keep.assign(n, 1);
  if (n != 0 && p.iterations != 0)
    for (std::size_t it = 0; it < p.iterations; ++it)
    {
      std::vector<std::size_t> active;
      for (std::size_t i = 0; i < n; ++i)
        if (r.keep[i]) active.push_back(i);
      if (active.empty()) break;
      std::vector<Pt> sub;
      for (std::size_t a : active) sub.push_back(f[a]);
      const std::vector<std::uint8_t> mask = refFilterOnce(sub, num_scans, p);
      std::vector<std::uint8_t> next(n, 0);
      for (std::size_t j = 0; j < active.size(); ++j)
        if (mask[j]) next[active[j]] = 1;
      const bool any = std::find(next.begin(), next.end(), 1) != next.end();
      r.keep = next;
      if (!any) break;
    }
  r.streak = kept(r.keep, 0, n);
  if (p.halo)
  {
    std::vector<std::size_t> idx;
    for (std::size_t i = 0; i < n; ++i)
      if (r.keep[i]) idx.push_back(i);
    std::vector<Pt> sub;
    for (std::size_t i : idx) sub.push_back(f[i]);
    const std::vector<std::uint8_t> h = refHalo(sub, num_scans, p);
    for (std::size_t k = 0; k < idx.size(); ++k)
      if (!h[k]) r.keep[idx[k]] = 0;
  }
  r.halo = kept(r.keep, 0, n);
  if (gate)
    for (std::size_t i = 0; i < n; ++i) r.keep[i] = r.keep[i] && gate->contains(f[i].scan, f[i].tof);
  return r;
}

struct Rng   // xorshift64: the same frames on every platform
{
  std::uint64_t s;
  u32 below(std::size_t n)
  {
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return static_cast<u32>(s % n);
  }
};

static void testRandomAgainstReference(std::size_t& frames, std::size_t& points, std::size_t& kept_points)
{
  Rng rng{0x9E3779B97F4A7C15ull};
  const u32 widths[] = {0, 1, 3, 7}, halo_mz[] = {0, 5, 80};
  const std::size_t lengths[] = {0, 1, 3, 5}, gaps[] = {0, 1, 2, 5}, iterations[] = {0, 1, 2, 3, 7}, halo_scan[] = {0, 1, 2, 50};
  const std::uint64_t wmins[] = {0, 0, 20, 60}, fmins[] = {0, 0, 50, 300};
  const double fractions[] = {0.15, 0.15, 0.5, 0.0, 1.0};
  const u32 top = std::numeric_limits<u32>::max();
  int bad = 0;
  const int trials = 3000;
  for (int trial = 0; trial < trials; ++trial)
  {
    const std::size_t num_scans = 1 + rng.below(40);
    const u32 base = trial % 3 == 0 ? 0 : trial % 3 == 1 ? 100000 : top - 150;
    const u32 span = 5 + rng.below(400);
    auto tofAt = [&](u32 off) { return off > top - base ? top : base + off; };
    std::vector<Pt> f(rng.below(600));
    for (Pt& q : f)
    {
      q.scan = rng.below(num_scans);
      q.tof = tofAt(rng.below(span));
      q.inten = rng.below(5) == 0 ? 0 : rng.below(100);
    }
    for (std::size_t k = 0; k < 20 && !f.empty(); ++k)   // duplicate (scan, TOF) pairs
    {
      Pt d = f[rng.below(f.size())];
      d.inten = rng.below(100);
      f[rng.below(f.size())] = d;
    }
    dn::Params p;
    p.mz_half_width = widths[rng.below(std::size(widths))];
    p.min_feature_length = lengths[rng.below(std::size(lengths))];
    p.max_internal_gap = gaps[rng.below(std::size(gaps))];
    p.min_window_intensity = wmins[rng.below(std::size(wmins))];
    p.min_feature_intensity = fmins[rng.below(std::size(fmins))];
    p.iterations = iterations[rng.below(std::size(iterations))];
    p.halo = rng.below(4) != 0;
    p.halo_peak_fraction = fractions[rng.below(std::size(fractions))];
    p.halo_mz_idx_half_width = halo_mz[rng.below(std::size(halo_mz))];
    p.halo_scan_half_width = halo_scan[rng.below(std::size(halo_scan))];
    std::optional<dn::Ms1WindowGate> gate;
    if (rng.below(2))
    {
      std::vector<dn::GateBox> boxes(1 + rng.below(4));
      for (dn::GateBox& b : boxes)
      {
        b.scan_lo = rng.below(50); b.scan_hi = rng.below(50);
        b.tof_lo = tofAt(rng.below(span)); b.tof_hi = tofAt(rng.below(span));
      }
      gate = dn::Ms1WindowGate::fromBoxes(boxes, static_cast<u32>(num_scans + rng.below(3)));
    }
    const dn::Ms1WindowGate* g = gate ? &*gate : nullptr;
    const RefResult want = refFrame(f, num_scans, p, g);
    const Run got = run(f, num_scans, p, g);
    bool ok = got.keep == want.keep && got.counts.after_streak == want.streak && got.counts.after_halo == want.halo &&
              got.counts.kept == kept(want.keep, 0, f.size());
    const Run back = run(std::vector<Pt>(f.rbegin(), f.rend()), num_scans, p, g);   // input order must not matter
    for (std::size_t i = 0; ok && i < f.size(); ++i) ok = back.keep[f.size() - 1 - i] == got.keep[i];
    if (!ok && ++bad <= 3)
      std::fprintf(stderr, "random frame %d (%zu points, %zu scans; w %u L %zu G %zu wmin %llu fmin %llu it %zu halo %d %.2f/%u/%zu gate %d) "
                           "differs from the transliteration\n", trial, f.size(), num_scans, p.mz_half_width, p.min_feature_length,
                   p.max_internal_gap, static_cast<unsigned long long>(p.min_window_intensity),
                   static_cast<unsigned long long>(p.min_feature_intensity), p.iterations, int(p.halo), p.halo_peak_fraction,
                   p.halo_mz_idx_half_width, p.halo_scan_half_width, int(g != nullptr));
    ++frames;
    points += f.size();
    kept_points += kept(want.keep, 0, f.size());
  }
  CHECK(bad == 0, "%d of %d random frames differ from the transliteration", bad, trials);
  CHECK(kept_points > points / 20 && kept_points < points / 2, "random frames too trivial: %zu of %zu kept", kept_points, points);
}

// ---- point recovery (recover::) ---------------------------------------------------------------------------------------
// The loader's forward maths, kept apart from the header's inverse: TdfMzCalibration::tofToMz (the patched loader's model),
// RationalScan2ImConverter::applyFormula followed by frameToSpectrum's float cast, and opentims' intensity correction.
// Inputs are metadata only -- MzCalibration, frame 1's T1, TimsCalibration, AccumulationTime -- of the two files above.
struct RecoverFile
{
  const char* name;
  diaspextract::TdfMzCalibration cal;
  double t1;                        ///< frame 1's T1
  long n_bins;                      ///< DigitizerNumSamples
  std::array<double, 10> tims;      ///< TimsCalibration C0..C9
  u32 n_scans;
  const char* accumulation_time;    ///< as sqlite renders it
  std::array<float, 4> im_gold;     ///< the loader's float32 1/K0 of scans 0, 1, n/2 and n-1
  u32 no_preimage;                  ///< the first corrected value no raw count gives
};

static RecoverFile recoverD()
{
  RecoverFile f{"D", {}, 25.727120813485598, 634073,
                {1, 943, 234.09826168388614, 95.59372590131042, 33.9622641509434, 1, -0.0031464178402676644, 167.9496150068565, 16.646316600032645, 2241.865411900982},
                944, "99.958", {0x1.6e4c14p+0f, 0x1.6e124p+0f, 0x1.03cfc8p+0f, 0x1.33a7e6p-1f}, 1190};
  f.cal.model_type = 1; f.cal.digitizer_timebase = 0.125; f.cal.digitizer_delay = 25655.375; f.cal.C0 = 279.3262846272992;
  f.cal.C1 = 155279.13067653627; f.cal.C2 = 0.001260061434461731; f.cal.T1_ref = 25.693668980735552; f.cal.dC1 = 20.0;
  return f;
}
static RecoverFile recoverTnbc()
{
  RecoverFile f{"TNBC 009", {}, 25.393315017452863, 399967,
                {1, 917, 217.9811466010841, 74.83820441307512, 33.0, 1, 0.02256006096271849, 131.32054866963384, 12.961442137235451, 2637.0446006790567},
                918, "99.953", {0x1.a33458p+0f, 0x1.a2ec16p+0f, 0x1.201aa4p+0f, 0x1.343b72p-1f}, 1064};
  f.cal.model_type = 1; f.cal.digitizer_timebase = 0.2; f.cal.digitizer_delay = 25030.2; f.cal.C0 = 314.59993079240667;
  f.cal.C1 = 155508.97802581752; f.cal.C2 = 0.0; f.cal.T1_ref = 25.41837029204028; f.cal.dC1 = 20.0;
  return f;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("fp-contract=off")
#endif
static float loaderIm(const std::array<double, 10>& c, double s)   // RationalScan2ImConverter.cpp:47-55, then a float
{
#if defined(__clang__)
  #pragma clang fp contract(off)
#endif
  double V = c[2] + ((c[3] - c[2]) / c[1]) * (s - c[4] - c[0]);
  if (V == 0.0) V = 1e-10;
  double denom = c[6] + c[7] / V;
  if (denom == 0.0) denom = 1e-10;
  return static_cast<float>(1.0 / denom);
}
static u32 opentimsCorrect(u32 raw, double corr)   // opentims.cpp:225
{
#if defined(__clang__)
  #pragma clang fp contract(off)
#endif
  return static_cast<u32>(static_cast<double>(raw) * corr + 0.5);
}
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

struct Row { int msms_type; u32 num_peaks; };
static void testRecover()
{
  namespace rc = dn::recover;
  using R = rc::Refusal;
  u32 tof = 0, scan = 0, raw = 0;
  double t = 0.0;
  int hits = 0;

  // intensity correction: 100 / AccumulationTime, the text parsed whole and in the C locale
  double corr = 0.0;
  CHECK(rc::opentimsCorrection("99.958", corr) && corr == 0x1.001b896437031p+0, "corr('99.958') = %a", corr);
  CHECK(rc::opentimsCorrection("99.953", corr) && corr == 0x1.001ed0ffd561bp+0, "corr('99.953') = %a", corr);
  CHECK(rc::opentimsCorrection("100", corr) && corr == 1.0, "corr('100') = %a", corr);
  for (const char* s : {"99,953", "abc", "0", "-0", "-99.958", "", " 99.958", "99.958 ", "99.958ms", "inf", "nan", "0x1p3", "1e-999"})
    CHECK(!rc::opentimsCorrection(s, corr), "opentimsCorrection(\"%s\") must refuse", s);
  bool comma_locale = false;
  if (const char* prev = std::setlocale(LC_NUMERIC, nullptr))
  {
    const std::string saved = prev;
    if (std::setlocale(LC_NUMERIC, "de_DE.UTF-8") && std::strtod("99.953", nullptr) == 99.0)   // atof would read 99 here
    {
      comma_locale = true;
      CHECK(rc::opentimsCorrection("99.953", corr) && corr == 100.0 / 99.953, "under de_DE, corr('99.953') = %a", corr);
    }
    std::setlocale(LC_NUMERIC, saved.c_str());
  }
  CHECK(rc::correctionRefusal(7, 0x1.001b896437031p+0).empty() && rc::correctionRefusal(7, 1.0).empty() &&
          rc::correctionRefusal(7, std::nextafter(1.0, 2.0)).empty(), "corr >= 1 is invertible");
  const std::string r125 = rc::correctionRefusal(7, 100.0 / 125.0);
  CHECK(r125.find("MS1 frame 7: AccumulationTime 125 ms exceeds 100 ms") == 0 && r125.find("= 0.8 < 1") != std::string::npos,
        "AccumulationTime 125 ms: '%s'", r125.c_str());
  CHECK(rc::opentimsCorrection("100.042", corr) && corr < 1.0 && !rc::correctionRefusal(1, corr).empty(), "AccumulationTime 100.042 ms is refused");
  CHECK(!rc::correctionRefusal(7, 0.5).empty() && !rc::correctionRefusal(7, std::nextafter(1.0, 0.0)).empty() &&
          !rc::correctionRefusal(7, 1e6).empty() && !rc::correctionRefusal(7, std::numeric_limits<double>::infinity()).empty(),
        "corr 0.5, 1 - 1 ulp, 1e6 and inf are refused");

  for (const RecoverFile& f : {recoverD(), recoverTnbc()})
  {
    const double b = f.cal.frameFactor(f.t1);
    double c = 0.0;
    rc::opentimsCorrection(f.accumulation_time, c);
    std::vector<float> asc;
    bool falls = false;
    const std::string why = rc::buildScanTable(f.tims, f.n_scans, asc, falls);
    CHECK(why.empty() && falls && asc.size() == f.n_scans, "%s: scan table (%s), falls %d", f.name, why.c_str(), int(falls));
    const float im0 = loaderIm(f.tims, 0.0);

    // TOF: every digitizer bin comes back; half a bin and +1 ppm off are refused
    long tof_bad = 0, half_bad = 0, ppm_bad = 0;
    for (long bin = 0; bin < f.n_bins; ++bin)
    {
      const double mz = f.cal.tofToMz(static_cast<double>(bin), b);
      tof_bad += rc::recoverPoint(f.cal, b, mz, im0, 0.0f, c, asc, falls, f.n_scans, tof, scan, raw, t, hits) != R::none || tof != static_cast<u32>(bin);
      half_bad += rc::recoverPoint(f.cal, b, f.cal.tofToMz(static_cast<double>(bin) + 0.5, b), im0, 0.0f, c, asc, falls, f.n_scans, tof, scan, raw, t, hits) != R::tof;
      ppm_bad += rc::recoverPoint(f.cal, b, mz * (1.0 + 1e-6), im0, 0.0f, c, asc, falls, f.n_scans, tof, scan, raw, t, hits) != R::tof;
    }
    CHECK(tof_bad == 0 && half_bad == 0 && ppm_bad == 0, "%s TOF: %ld of %ld bins lost, %ld half-bin and %ld +1 ppm m/z not refused",
          f.name, tof_bad, f.n_bins, half_bad, ppm_bad);

    // scan: the loader's float32 goldens; every scan comes back, also 1 ulp off; half a scan off matches none
    const u32 at[4] = {0, 1, f.n_scans / 2, f.n_scans - 1};
    for (int k = 0; k < 4; ++k)
      CHECK(loaderIm(f.tims, at[k]) == f.im_gold[static_cast<std::size_t>(k)], "%s scan %u: 1/K0 %a, want %a", f.name, at[k],
            loaderIm(f.tims, at[k]), f.im_gold[static_cast<std::size_t>(k)]);
    long scan_bad = 0, ulp_bad = 0, half_hit = 0;
    for (u32 s = 0; s < f.n_scans; ++s)
    {
      const float v = loaderIm(f.tims, s);
      scan_bad += rc::scanOf(asc, falls, v) != static_cast<long>(s);
      ulp_bad += rc::scanOf(asc, falls, std::nextafter(v, 10.0f)) != static_cast<long>(s) || rc::scanOf(asc, falls, std::nextafter(v, 0.0f)) != static_cast<long>(s);
      half_hit += rc::scanOf(asc, falls, loaderIm(f.tims, s + 0.5)) >= 0;
    }
    CHECK(scan_bad == 0 && ulp_bad == 0 && half_hit == 0, "%s scans: %ld lost, %ld lost 1 ulp off, %ld half-scan values matched", f.name, scan_bad, ulp_bad, half_hit);
    const double mz = f.cal.tofToMz(1000.0, b);
    const float last = loaderIm(f.tims, f.n_scans - 1);
    CHECK(rc::recoverPoint(f.cal, b, mz, last, 0.0f, c, asc, falls, f.n_scans - 1, tof, scan, raw, t, hits) == R::scan,
          "%s: scan %u in a frame of %u scans is refused", f.name, f.n_scans - 1, f.n_scans - 1);
    CHECK(rc::recoverPoint(f.cal, b, mz, last, 0.0f, c, asc, falls, f.n_scans, tof, scan, raw, t, hits) == R::none && scan == f.n_scans - 1 && tof == 1000,
          "%s: the last scan of a whole frame", f.name);

    // intensity: every raw count whose correction stays below 2^24 comes back uniquely -- all below 2^20, every 7th above
    long raw_bad = 0, raws = 0;
    for (u32 r = 0;; r += r < (1u << 20) ? 1 : 7)
    {
      const u32 v = opentimsCorrect(r, c);
      if (v >= 16777216u) break;
      raw_bad += rc::recoverPoint(f.cal, b, mz, im0, static_cast<float>(v), c, asc, falls, f.n_scans, tof, scan, raw, t, hits) != R::none || raw != r || hits != 1;
      ++raws;
    }
    CHECK(raw_bad == 0 && raws > 3000000, "%s intensity: %ld of %ld raw counts do not come back uniquely", f.name, raw_bad, raws);
    CHECK(rc::recoverPoint(f.cal, b, mz, im0, static_cast<float>(f.no_preimage), c, asc, falls, f.n_scans, tof, scan, raw, t, hits) == R::intensity && hits == 0,
          "%s: corrected value %u has no raw count (%d found)", f.name, f.no_preimage, hits);
  }

  // intensity edge cases, on D's frame geometry
  const RecoverFile d = recoverD();
  const double b = d.cal.frameFactor(d.t1), mz = d.cal.tofToMz(1000.0, b);
  std::vector<float> asc;
  bool falls = false;
  rc::buildScanTable(d.tims, d.n_scans, asc, falls);
  const float im0 = loaderIm(d.tims, 0.0);
  auto point = [&](float v, double c) { return rc::recoverPoint(d.cal, b, mz, im0, v, c, asc, falls, d.n_scans, tof, scan, raw, t, hits); };
  rc::opentimsCorrection("99.958", corr);
  CHECK(point(0.0f, corr) == R::none && raw == 0 && hits == 1, "v 0 is raw 0");
  CHECK(point(12.5f, corr) == R::intensity && hits == 0 && point(-1.0f, corr) == R::intensity && point(std::nanf(""), corr) == R::intensity,
        "a fraction, a negative and NaN are refused");
  CHECK(point(16777216.0f, corr) == R::intensity_range && point(std::numeric_limits<float>::infinity(), corr) == R::intensity_range &&
          point(16777215.0f, 1.0) == R::none && raw == 16777215, "2^24 and inf are out of range, 2^24 - 1 is not");
  // opentims' operation order decides ties: at '99.968', raw 1562 corrects to 1562 ((raw * 100) / at would give 1563)
  rc::opentimsCorrection("99.968", corr);
  CHECK(opentimsCorrect(1562, corr) == 1562 && point(1562.0f, corr) == R::none && raw == 1562 && hits == 1 && point(1563.0f, corr) == R::intensity && hits == 0,
        "the '99.968' tie: 1562 recovers, 1563 has no raw count");
  rc::opentimsCorrection("125", corr);
  CHECK(point(2.0f, corr) == R::intensity && hits == 2, "corr 0.8: raw counts 2 and 3 both give 2 (%d found)", hits);
  rc::opentimsCorrection("50", corr);
  CHECK(opentimsCorrect(1u << 23, corr) == 16777216u && point(16777216.0f, corr) == R::intensity_range, "corr 2: raw 2^23 gives 2^24, refused");
  // the corr >= 1 threshold: no two raw counts collide at corr 1, 1 + 1 ulp or 100/96 (tie-prone), below 2^20
  long premise_bad = 0;
  for (const double c : {1.0, std::nextafter(1.0, 2.0), 100.0 / 96.0})
    for (u32 r = 0; r < (1u << 20); ++r)
      premise_bad += point(static_cast<float>(opentimsCorrect(r, c)), c) != R::none || raw != r;
  CHECK(premise_bad == 0, "corr >= 1: %ld raw counts below 2^20 do not come back", premise_bad);

  // scan tables span the scans their MS1 frames have (TdfDnoiseInputs::ms1_max_scans): TNBC's row separates its 918 scans
  // but not 40,000 (its pole is near scan 38,722); a row with its pole inside 918 scans, or a flat one, separates none,
  // yet over 0 scans -- no MS1 frame references it -- it is an empty table that matches nothing
  const RecoverFile tn = recoverTnbc();
  CHECK(rc::buildScanTable(tn.tims, 918, asc, falls).empty() && !rc::buildScanTable(tn.tims, 40000, asc, falls).empty(),
        "TNBC's row: 918 scans separate, 40,000 do not");
  std::array<double, 10> pole = tn.tims, flat = d.tims;
  pole[4] -= 38300.0;
  flat[3] = flat[2];
  CHECK(!rc::buildScanTable(pole, 918, asc, falls).empty() && !rc::buildScanTable(flat, 944, asc, falls).empty(), "a pole or a flat row is refused");
  CHECK(rc::buildScanTable(pole, 0, asc, falls).empty() && asc.empty() && rc::scanOf(asc, falls, im0) == -1, "over 0 scans: an empty table");
  CHECK(rc::buildScanTable(d.tims, 1, asc, falls).empty() && asc.size() == 1 && !falls && rc::scanOf(asc, falls, im0) == 0, "one scan");

  // frame checks, in order
  const std::vector<Row> rows = {{-1, 0}, {0, 5}, {9, 5}, {0, 0}};
  CHECK(rc::frameCheck(rows, 0, 5, true, 5) == R::frame && rc::frameCheck(rows, 4, 5, true, 5) == R::frame && rc::frameCheck(rows, 2, 5, true, 5) == R::frame,
        "frame 0, a frame past the table and an MS2 frame are refused");
  CHECK(rc::frameCheck(rows, 1, 4, true, 4) == R::peaks && rc::frameCheck(rows, 1, 6, true, 6) == R::peaks, "a spectrum must hold NumPeaks points");
  CHECK(rc::frameCheck(rows, 1, 5, false, 0) == R::im && rc::frameCheck(rows, 1, 5, true, 4) == R::im, "points need a 1/K0 each");
  CHECK(rc::frameCheck(rows, 1, 5, true, 5) == R::none && rc::frameCheck(rows, 3, 0, false, 0) == R::none, "a whole frame; an empty frame needs no 1/K0 array");

  if (g_fail == 0) std::printf("    point recovery%s\n", comma_locale ? " (also under a de_DE LC_NUMERIC)" : " (no de_DE locale here: the LC_NUMERIC case is skipped)");
}

static void testCoverage()
{
  const std::vector<Row> frames = {{-1, 0}, {0, 5}, {9, 7}, {0, 0}, {0, 2}};   // 1 MS1, 2 MS2, 3 an MS1 frame without points, 4 MS1
  dn::FrameCoverage cov;
  CHECK(cov.firstMissing(frames) == 1, "nothing marked: frame 1 is missing");
  CHECK(cov.mark(1) && !cov.mark(1), "a frame denoised twice is refused");
  CHECK(cov.firstMissing(frames) == 4, "frame 4 is missing, not the MS2 frame 2 or the empty MS1 frame 3");
  CHECK(cov.mark(4) && cov.firstMissing(frames) == 0, "complete");
  CHECK(cov.mark(3) && !cov.mark(3), "an empty MS1 frame that arrives is marked, and refused the second time");
  cov.reset();
  CHECK(cov.firstMissing(frames) == 1 && cov.mark(1) && cov.mark(4) && cov.firstMissing(frames) == 0, "reset clears every mark");
}

int main()
{
  int before = g_fail;
  testGate();
  testParse();
  if (g_fail == before)
    std::printf("OK  gate: D 24 boxes (FMA-sensitive 11/15/16 intact), TNBC 009 28 boxes (otofControl); lookup, refusal and parse rules\n");
  before = g_fail;
  testStreak();
  testHalo();
  testGateAnd();
  if (g_fail == before)
    std::printf("OK  streak, halo and gate rules on synthetic frames; 0.15 double compare == integer rule\n");
  before = g_fail;
  std::size_t frames = 0, points = 0, kept_points = 0;
  testRandomAgainstReference(frames, points, kept_points);
  if (g_fail == before)
    std::printf("OK  %zu random frames (%zu points, %zu kept) identical to a transliteration of the Rust, in either input order\n",
                frames, points, kept_points);
  before = g_fail;
  testRecover();
  testCoverage();
  if (g_fail == before)
    std::printf("OK  point recovery on D and TNBC 009 metadata: every TOF bin and scan, raw counts to 2^24; correction, frame and coverage rules\n");
  if (g_fail) std::fprintf(stderr, "%d check(s) failed\n", g_fail);
  return g_fail ? 1 : 0;
}
