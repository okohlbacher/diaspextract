// Tests of src/TdfLoad.h: every analysis.tdf read the tool makes is whole or refused, and each refusal names what it
// refuses. Needs SQLite, not OpenMS: CMake builds it only where it finds SQLite3, with a shim include directory that maps
// <OpenMS/FORMAT/TdfMzCalibration.h> onto src/TdfMzCalibration.h. The fixtures are synthetic tdfs shaped like dataset D's,
// written through the sqlite3 C API into a temporary directory; corrupted ones get 12 bytes of 0xff over a page header.
#include "TdfLoad.h"
#include <chrono>
#include <clocale>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace dn = spx::dnoise;
using diaspextract::TdfDnoiseInputs;

static int g_fail = 0;
#define CHECK(cond, ...) do { if (!(cond)) { std::fprintf(stderr, "FAIL line %d: ", __LINE__); std::fprintf(stderr, __VA_ARGS__); std::fputc('\n', stderr); ++g_fail; } } while (0)

static std::filesystem::path g_dir;
static constexpr long kPageSize = 4096;

[[noreturn]] static void fixtureFailed(const std::string& what)
{
  std::fprintf(stderr, "fixture failed: %s\n", what.c_str());
  std::exit(2);
}

/// D's metadata, MzCalibration and TimsCalibration rows; @p n_frames frames (every fifth MS1, 944 scans, 1000 points,
/// AccumulationTime 99.958); @p n_windows window rows in a table named @p windows (nullptr: no table). Frames columns carry
/// no type, so a test can store TEXT or REAL where Bruker stores INTEGER.
static std::string dShaped(int n_frames = 10, int n_windows = 24, const char* windows = "DiaFrameMsMsWindows")
{
  std::string s =
    "PRAGMA page_size = " + std::to_string(kPageSize) + ";"
    "CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);"
    "INSERT INTO GlobalMetadata VALUES ('AcquisitionSoftware', 'timsTOF'), ('DigitizerNumSamples', '634073'), ('MzAcqRangeLower', '99.990834'),"
    "  ('MzAcqRangeUpper', '1700.000000'), ('OneOverK0AcqRangeLower', '0.600000'), ('OneOverK0AcqRangeUpper', '1.400000');"
    "CREATE TABLE MzCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, DigitizerTimebase, DigitizerDelay, T1, dC1, C0, C1, C2, dC2, C3, C4);"
    "INSERT INTO MzCalibration VALUES (1, 1, 0.125, 25655.375, 25.693668980735552, 20, 279.3262846272992, 155279.13067653627, 0.001260061434461731, 0, 0, 0);"
    "CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9);"
    "INSERT INTO TimsCalibration VALUES (1, 2, 1, 943, 234.09826168388614, 95.59372590131042, 33.9622641509434, 1, -0.0031464178402676644,"
    "  167.9496150068565, 16.646316600032645, 2241.865411900982);"
    "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration, MzCalibration, T1);"
    "WITH RECURSIVE k(i) AS (SELECT 1 UNION ALL SELECT i + 1 FROM k WHERE i < " + std::to_string(n_frames) + ")"
    "  INSERT INTO Frames SELECT i, CASE WHEN i % 5 = 1 THEN 0 ELSE 9 END, 944, 1000, 99.958, 1, 1, 25.727120813485598 + i * 1e-6 FROM k;";
  if (windows)
    s += std::string("CREATE TABLE ") + windows + " (WindowGroup INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL);"
         "WITH RECURSIVE k(i) AS (SELECT 0 UNION ALL SELECT i + 1 FROM k WHERE i < " + std::to_string(n_windows - 1) + ")"
         "  INSERT INTO " + windows + " SELECT 1 + i % 12, 34 + (i * 7) % 400, 500 + (i * 13) % 400, 400.0 + 0.5 * i, 25.0, 30.0 FROM k;";
  return s;
}

/// @p table rebuilt with only the columns @p keep: a column dropped without ALTER TABLE DROP COLUMN, which needs SQLite 3.35.
static std::string keepColumns(const std::string& table, const char* keep)
{
  return "CREATE TABLE " + table + "_kept AS SELECT " + keep + " FROM " + table + "; DROP TABLE " + table + "; ALTER TABLE " + table +
         "_kept RENAME TO " + table + ";";
}

static std::string makeTdf(const std::string& name, const std::string& sql)
{
  const std::string path = (g_dir / name).string();
  std::filesystem::remove(path);
  sqlite3* db = nullptr;
  if (sqlite3_open(path.c_str(), &db) != SQLITE_OK) fixtureFailed("cannot create " + path);
  char* err = nullptr;
  if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &err) != SQLITE_OK)
    fixtureFailed(name + ": " + (err ? err : "?") + "\n" + sql);
  sqlite3_close(db);
  return path;
}

static long long rootPage(const std::string& path, const char* table)
{
  sqlite3* db = nullptr;
  sqlite3_stmt* st = nullptr;
  long long page = 0;
  if (sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) == SQLITE_OK &&
      sqlite3_prepare_v2(db, "SELECT rootpage FROM sqlite_master WHERE name = ?", -1, &st, nullptr) == SQLITE_OK &&
      sqlite3_bind_text(st, 1, table, -1, SQLITE_STATIC) == SQLITE_OK && sqlite3_step(st) == SQLITE_ROW)
    page = sqlite3_column_int64(st, 0);
  sqlite3_finalize(st);
  sqlite3_close(db);
  if (page < 1) fixtureFailed(std::string("no root page for ") + table + " in " + path);
  return page;
}

static std::vector<unsigned char> readPage(const std::string& path, long long page)
{
  std::ifstream in(path, std::ios::binary);
  in.seekg((page - 1) * kPageSize);
  std::vector<unsigned char> b(static_cast<std::size_t>(kPageSize));
  if (!in.read(reinterpret_cast<char*>(b.data()), kPageSize)) fixtureFailed("cannot read page " + std::to_string(page) + " of " + path);
  return b;
}

static void corruptPage(const std::string& path, long long page)
{
  std::fstream io(path, std::ios::in | std::ios::out | std::ios::binary);
  io.seekp((page - 1) * kPageSize);
  const char ff[12] = {'\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff', '\xff'};
  if (!io.write(ff, sizeof ff)) fixtureFailed("cannot corrupt page " + std::to_string(page) + " of " + path);
}

/// Corrupt the middle leaf of @p table, whose root must be an interior table page (0x05) over leaf pages (0x0d).
static void corruptMiddleLeaf(const std::string& path, const char* table)
{
  const std::vector<unsigned char> root = readPage(path, rootPage(path, table));
  if (root[0] != 0x05) fixtureFailed(std::string(table) + "'s root is not an interior page");
  const auto be = [](const unsigned char* p, int n) { std::uint32_t v = 0; for (int i = 0; i < n; ++i) v = (v << 8) | p[i]; return v; };
  const std::uint32_t cells = be(&root[3], 2);
  std::vector<std::uint32_t> kids;
  for (std::uint32_t i = 0; i < cells; ++i) kids.push_back(be(&root[be(&root[12 + 2 * i], 2)], 4));   // each cell starts with its left child
  kids.push_back(be(&root[8], 4));                                                                  // and the right-most child
  const std::uint32_t leaf = kids[kids.size() / 2];
  if (kids.size() < 3 || readPage(path, leaf)[0] != 0x0d) fixtureFailed(std::string(table) + ": the middle child is not a leaf");
  corruptPage(path, leaf);
}

static bool has(const std::string& s, const char* part) { return s.find(part) != std::string::npos; }

static bool loads(const std::string& path, TdfDnoiseInputs& in, std::string& why) { return diaspextract::loadTdfDnoise(path, in, why); }

/// loadTdfDnoise must refuse @p sql's tdf with a reason that contains @p part.
static void refused(const char* name, const std::string& sql, const char* part)
{
  TdfDnoiseInputs in;
  std::string why;
  const bool ok = loads(makeTdf(std::string(name) + ".tdf", sql), in, why);
  CHECK(!ok && has(why, part), "%s: loadTdfDnoise %s (%s), want a refusal naming '%s'", name, ok ? "accepted" : "refused", why.c_str(), part);
}

static void testWindows()
{
  TdfDnoiseInputs in;
  std::string why;
  const std::string ok = makeTdf("ok24.tdf", dShaped());
  CHECK(loads(ok, in, why) && in.windows.size() == 24 && in.meta.max_num_scans == 944 && in.frames.size() == 11 && in.tims_cal.size() == 1,
        "ok24: %s, %zu windows, %u scans, %zu frame slots", why.c_str(), in.windows.size(), in.meta.max_num_scans, in.frames.size());
  const dn::GateBuild g = dn::buildMs1WindowGate(in.meta, in.windows, dn::Params{});
  CHECK(g.error.empty() && g.gate && g.gate->boxes().size() == 24, "ok24: gate (%s)", g.error.c_str());
  CHECK(in.meta.acquisition_software == std::optional<std::string>("timsTOF") && in.meta.digitizer_num_samples == std::optional<std::string>("634073"),
        "ok24: GlobalMetadata as stored text");

  const std::string absent = makeTdf("absent.tdf", dShaped(10, 0, nullptr));
  CHECK(loads(absent, in, why) && in.windows.empty() && !dn::buildMs1WindowGate(in.meta, in.windows, dn::Params{}).gate,
        "absent: no DiaFrameMsMsWindows table is no gate, not a refusal (%s)", why.c_str());
  CHECK(loads(makeTdf("lower.tdf", dShaped(10, 24, "diaframemsmswindows")), in, why) && in.windows.size() == 24,
        "the table name resolves case-insensitively: %s, %zu windows", why.c_str(), in.windows.size());
  refused("nocol", dShaped() + keepColumns("DiaFrameMsMsWindows", "WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, CollisionEnergy"), "IsolationWidth");
  refused("null_begin", dShaped() + "UPDATE DiaFrameMsMsWindows SET ScanNumBegin = NULL WHERE rowid = 3;", "NULL ScanNumBegin");
  refused("null_end", dShaped() + "UPDATE DiaFrameMsMsWindows SET ScanNumEnd = NULL WHERE rowid = 3;", "ScanNumEnd");

  // corruption: a whole one-page table and one leaf of eighteen used to read as no gate and as half a gate
  const std::string c24 = makeTdf("corrupt24.tdf", dShaped());
  corruptPage(c24, rootPage(c24, "DiaFrameMsMsWindows"));
  CHECK(!loads(c24, in, why) && has(why, "malformed") && has(why, "DiaFrameMsMsWindows"), "corrupt24: %s", why.c_str());
  CHECK(loads(makeTdf("ok3000.tdf", dShaped(10, 3000)), in, why) && in.windows.size() == 3000, "ok3000: %s, %zu windows", why.c_str(), in.windows.size());
  const std::string c3000 = makeTdf("corruptmid3000.tdf", dShaped(10, 3000));
  corruptMiddleLeaf(c3000, "DiaFrameMsMsWindows");
  CHECK(!loads(c3000, in, why) && has(why, "malformed"), "corruptmid3000: %s (%zu windows)", why.c_str(), in.windows.size());
}

static void testFrames()
{
  TdfDnoiseInputs in;
  std::string why;
  CHECK(loads(makeTdf("frames_ok.tdf", dShaped()), in, why), "D-shaped: %s", why.c_str());
  const TdfDnoiseInputs::Frame f1 = in.frames.size() > 2 ? in.frames[1] : TdfDnoiseInputs::Frame{}, f2 = in.frames.size() > 2 ? in.frames[2] : TdfDnoiseInputs::Frame{};
  double corr = 0.0;
  CHECK(dn::recover::opentimsCorrection("99.958", corr) && f1.msms_type == 0 && f1.num_scans == 944 && f1.num_peaks == 1000 && f1.corr == corr && f1.cal == 0 &&
          f2.msms_type == 9 && in.frames[0].msms_type == -1, "D-shaped frames: MS1 %d/%u/%u corr %a, MS2 %d", f1.msms_type, f1.num_scans, f1.num_peaks, f1.corr, f2.msms_type);
  CHECK(in.ms1_max_scans == std::vector<std::uint32_t>{944}, "ms1_max_scans");
  // the tool instantiates recover:: and FrameCoverage on TdfDnoiseInputs::Frame
  CHECK(dn::recover::frameCheck(in.frames, 1, 1000, true, 1000) == dn::recover::Refusal::none &&
          dn::recover::frameCheck(in.frames, 2, 1000, true, 1000) == dn::recover::Refusal::frame, "frameCheck on the loaded table");
  dn::FrameCoverage cov;
  CHECK(cov.firstMissing(in.frames) == 1 && cov.mark(1) && cov.mark(6) && cov.firstMissing(in.frames) == 0, "FrameCoverage on the loaded table");

  // every integer is typed and ranged before it is narrowed
  refused("scans_neg", dShaped() + "UPDATE Frames SET NumScans = -1 WHERE Id = 3;", "NumScans");
  refused("scans_2p32", dShaped() + "UPDATE Frames SET NumScans = 4294967297 WHERE Id = 3;", "NumScans");
  refused("scans_text", dShaped() + "UPDATE Frames SET NumScans = '944' WHERE Id = 3;", "NumScans");
  refused("scans_real", dShaped() + "UPDATE Frames SET NumScans = 100.5 WHERE Id = 3;", "NumScans");
  refused("peaks_2p32", dShaped() + "UPDATE Frames SET NumPeaks = 4294967301 WHERE Id = 3;", "NumPeaks");
  refused("peaks_neg", dShaped() + "UPDATE Frames SET NumPeaks = -1 WHERE Id = 3;", "NumPeaks");
  refused("type_2p31", dShaped() + "UPDATE Frames SET MsMsType = 2147483648 WHERE Id = 3;", "MsMsType");
  refused("type_text", dShaped() + "UPDATE Frames SET MsMsType = '0' WHERE Id = 3;", "MsMsType");
  refused("cal_2p32", dShaped() + "UPDATE Frames SET TimsCalibration = 4294967297 WHERE Id = 3;", "TimsCalibration");
  refused("cal_null", dShaped() + "UPDATE Frames SET TimsCalibration = NULL WHERE Id = 3;", "TimsCalibration");
  refused("cal_unknown", dShaped() + "UPDATE Frames SET TimsCalibration = 7 WHERE Id = 3;", "unknown TimsCalibration");
  refused("calid_2p32", dShaped() + "INSERT INTO TimsCalibration SELECT 4294967297, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration;",
          "TimsCalibration.Id");
  refused("ms1_no_scans", dShaped() + "UPDATE Frames SET NumScans = 0 WHERE Id = 6;", "no scans");
  CHECK(loads(makeTdf("ms1_empty.tdf", dShaped() + "UPDATE Frames SET NumScans = 0, NumPeaks = 0 WHERE Id = 6;"), in, why) && in.frames[6].num_peaks == 0,
        "an MS1 frame with neither scans nor points is accepted: %s", why.c_str());
  CHECK(loads(makeTdf("ms2_no_scans.tdf", dShaped() + "UPDATE Frames SET NumScans = 0 WHERE Id = 7;"), in, why), "an MS2 frame without scans is not the port's: %s", why.c_str());
  refused("no_frames", dShaped() + "DELETE FROM Frames;", "no Frames rows");
  refused("no_metadata", dShaped() + "DROP TABLE GlobalMetadata;", "GlobalMetadata");
  refused("no_timscal", dShaped() + "DROP TABLE TimsCalibration;", "linear");
  refused("no_timscal_col", dShaped() + keepColumns("Frames", "Id, MsMsType, NumScans, NumPeaks, AccumulationTime, MzCalibration, T1"), "TimsCalibration");
  refused("modeltype", dShaped() + "UPDATE TimsCalibration SET ModelType = 1;", "ModelType 1");
  {
    std::string w;
    CHECK(!diaspextract::loadTdfDnoise((g_dir / "does_not_exist.tdf").string(), in, w) && has(w, "cannot open"), "a missing file: %s", w.c_str());
  }

  // a scan table per TimsCalibration row over its own MS1 frames' scans: MS1 frames on a TNBC-shaped row with 918 scans,
  // MS2 frames on a flat row with 40,000 -- the file's MAX(NumScans), which the gate still uses
  const std::string two = dShaped() +
    "UPDATE TimsCalibration SET C3 = C2;"
    "INSERT INTO TimsCalibration VALUES (2, 2, 1, 917, 217.9811466010841, 74.83820441307512, 33.0, 1, 0.02256006096271849, 131.32054866963384,"
    "  12.961442137235451, 2637.0446006790567);"
    "UPDATE Frames SET TimsCalibration = 2, NumScans = 918 WHERE MsMsType = 0;"
    "UPDATE Frames SET NumScans = 40000 WHERE MsMsType <> 0;";
  CHECK(loads(makeTdf("two_rows.tdf", two), in, why) && in.meta.max_num_scans == 40000 && in.ms1_max_scans == std::vector<std::uint32_t>({0, 918}) &&
          in.frames[1].cal == 1 && in.frames[2].cal == 0, "two rows: %s, max %u", why.c_str(), in.meta.max_num_scans);
  std::vector<float> asc;
  bool falls = false;
  CHECK(in.tims_cal.size() == 2 && in.ms1_max_scans.size() == 2 && dn::recover::buildScanTable(in.tims_cal[0], in.ms1_max_scans[0], asc, falls).empty() &&
          dn::recover::buildScanTable(in.tims_cal[1], in.ms1_max_scans[1], asc, falls).empty() &&
          !dn::recover::buildScanTable(in.tims_cal[1], in.meta.max_num_scans, asc, falls).empty() &&
          !dn::recover::buildScanTable(in.tims_cal[0], in.meta.max_num_scans, asc, falls).empty(),
        "each row's table over its MS1 scans builds; over the file's MAX neither would");

  // corruption in the middle of Frames: refused at setup, and by the calibration read (whose T1 loop used to stop early)
  const std::string ok3000 = makeTdf("frames3000.tdf", dShaped(3000));
  {
    diaspextract::TdfMzCalibration cal;
    std::vector<double> t1;
    std::string w;
    CHECK(loads(ok3000, in, why) && in.frames.size() == 3001, "frames3000: %s", why.c_str());
    CHECK(diaspextract::loadTdfCalibration(ok3000, cal, t1, w) && t1.size() == 3001 && t1[3000] == 25.727120813485598 + 3000 * 1e-6 && cal.C1 == 155279.13067653627,
          "loadTdfCalibration frames3000: %s (%zu T1)", w.c_str(), t1.size());
    const std::string bad = makeTdf("framescorrupt3000.tdf", dShaped(3000));
    corruptMiddleLeaf(bad, "Frames");
    CHECK(!loads(bad, in, why) && has(why, "malformed"), "a corrupt Frames leaf: %s", why.c_str());
    CHECK(!diaspextract::loadTdfCalibration(bad, cal, t1, w) && has(w, "malformed"), "loadTdfCalibration, a corrupt Frames leaf: %s", w.c_str());
  }
}

static void testAccumulationTime()
{
  TdfDnoiseInputs in;
  std::string why;
  int k = 0;
  for (const char* text : {"99,953", "abc", "0", "-99.958", " 99.958", "99.958ms", ""})
    refused(("acc_bad" + std::to_string(k++)).c_str(), dShaped() + "UPDATE Frames SET AccumulationTime = '" + text + "' WHERE Id = 6;", "AccumulationTime");
  refused("acc_null", dShaped() + "UPDATE Frames SET AccumulationTime = NULL WHERE Id = 6;", "AccumulationTime");
  CHECK(loads(makeTdf("acc_text.tdf", dShaped() + "UPDATE Frames SET AccumulationTime = '99.953' WHERE Id = 6;"), in, why) &&
          in.frames[6].corr == 100.0 / 99.953 && in.frames[1].corr == 100.0 / 99.958, "AccumulationTime as TEXT and as REAL: %s", why.c_str());
  // above 100 ms the correction loads, and the setup refusal names it
  CHECK(loads(makeTdf("acc_125.tdf", dShaped() + "UPDATE Frames SET AccumulationTime = '125' WHERE Id = 6;"), in, why) && in.frames[6].corr == 0.8 &&
          has(dn::recover::correctionRefusal(6, in.frames[6].corr), "AccumulationTime 125 ms exceeds 100 ms") &&
          dn::recover::correctionRefusal(1, in.frames[1].corr).empty(), "AccumulationTime 125: %s", why.c_str());
  CHECK(loads(makeTdf("acc_100042.tdf", dShaped() + "UPDATE Frames SET AccumulationTime = 100.042 WHERE Id = 6;"), in, why) && in.frames[6].corr < 1.0 &&
          !dn::recover::correctionRefusal(6, in.frames[6].corr).empty(), "AccumulationTime 100.042 (REAL): %s", why.c_str());
  // opentims parses under a C-locale guard; so must the port, whatever LC_NUMERIC the process has
  bool comma = false;
  if (const char* prev = std::setlocale(LC_NUMERIC, nullptr))
  {
    const std::string saved = prev;
    if (std::setlocale(LC_NUMERIC, "de_DE.UTF-8") && std::strtod("99.953", nullptr) == 99.0)
    {
      comma = true;
      CHECK(loads(makeTdf("acc_de.tdf", dShaped() + "UPDATE Frames SET AccumulationTime = '99.953' WHERE Id = 6;"), in, why) &&
              in.frames[6].corr == 100.0 / 99.953 && in.frames[1].corr == 100.0 / 99.958, "under de_DE: %s", why.c_str());
      refused("acc_de_comma", dShaped() + "UPDATE Frames SET AccumulationTime = '99,953' WHERE Id = 6;", "AccumulationTime");
      diaspextract::TdfAxisBounds bd;   // the preflight's GlobalMetadata bounds parse in the C locale too
      CHECK(diaspextract::loadTdfAxisBounds(makeTdf("bounds_de.tdf", dShaped()), bd, why) && bd.mz_lo == 99.990834 && bd.mz_hi == 1700.0 && bd.n_bins == 634073,
            "loadTdfAxisBounds under de_DE: %s", why.c_str());
    }
    std::setlocale(LC_NUMERIC, saved.c_str());
  }
  if (!comma) std::printf("    (no de_DE locale here: the LC_NUMERIC case is skipped)\n");
}

static void testCalibrationAndBounds()
{
  diaspextract::TdfMzCalibration cal;
  std::vector<double> t1;
  std::string why;
  const std::string ok = makeTdf("cal_ok.tdf", dShaped());
  CHECK(diaspextract::loadTdfCalibration(ok, cal, t1, why) && t1.size() == 11 && cal.digitizer_timebase == 0.125 && cal.T1_ref == 25.693668980735552,
        "loadTdfCalibration: %s", why.c_str());
  CHECK(diaspextract::loadTdfCalibration(makeTdf("cal_nocol.tdf", dShaped() + keepColumns("Frames", "Id, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration, T1")), cal, t1, why),
        "no Frames.MzCalibration column and one MzCalibration row: %s", why.c_str());
  CHECK(!diaspextract::loadTdfCalibration(makeTdf("cal_nocol2.tdf", dShaped() + keepColumns("Frames", "Id, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration, T1") +
                                                 "INSERT INTO MzCalibration SELECT 2, ModelType, DigitizerTimebase, DigitizerDelay, T1, dC1, C0, C1, C2, dC2, C3, C4 FROM MzCalibration;"),
                                         cal, t1, why) && has(why, "need exactly 1"), "no column and two rows: %s", why.c_str());
  // a generated Frames.MzCalibration column (SQLite 3.31+) is the frames' reference like a stored one: PRAGMA table_info skips it
  if (sqlite3_libversion_number() >= 3031000)
    CHECK(diaspextract::loadTdfCalibration(makeTdf("cal_generated.tdf", dShaped() +
                                                   "INSERT INTO MzCalibration SELECT 2, ModelType, DigitizerTimebase, DigitizerDelay, T1, dC1, C0, 2 * C1, C2, dC2, C3, C4 FROM MzCalibration;"
                                                   "ALTER TABLE Frames RENAME TO F0;"
                                                   "CREATE TABLE Frames (Id INTEGER PRIMARY KEY, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration, T1,"
                                                   "  MzCalibration GENERATED ALWAYS AS (1));"
                                                   "INSERT INTO Frames SELECT Id, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration, T1 FROM F0; DROP TABLE F0;"),
                                           cal, t1, why) && cal.C1 == 155279.13067653627 && t1.size() == 11,
          "a generated MzCalibration column and two rows: %s", why.c_str());
  CHECK(!diaspextract::loadTdfCalibration(makeTdf("cal_two.tdf", dShaped() + "INSERT INTO MzCalibration SELECT 2, ModelType, DigitizerTimebase, DigitizerDelay, T1, dC1, C0, C1, C2, dC2, C3, C4 FROM MzCalibration;"
                                                 "UPDATE Frames SET MzCalibration = 2 WHERE Id = 4;"), cal, t1, why) && has(why, "2 distinct"),
        "frames on two MzCalibration rows: %s", why.c_str());
  CHECK(!diaspextract::loadTdfCalibration(makeTdf("cal_notable.tdf", dShaped() + "DROP TABLE MzCalibration;"), cal, t1, why) && has(why, "MzCalibration"),
        "no MzCalibration table: %s", why.c_str());

  diaspextract::TdfAxisBounds bd;
  CHECK(diaspextract::loadTdfAxisBounds(ok, bd, why) && bd.n_bins == 634073 && bd.mz_lo == 99.990834 && bd.mz_hi == 1700.0, "loadTdfAxisBounds: %s", why.c_str());
  CHECK(!diaspextract::loadTdfAxisBounds(makeTdf("bounds_notable.tdf", dShaped() + "DROP TABLE GlobalMetadata;"), bd, why) && has(why, "GlobalMetadata"),
        "no GlobalMetadata: %s", why.c_str());
  // the two bounds are the negative-C2 domain's: 1e8 bins and m/z 1e5 are accepted, one more is not
  CHECK(diaspextract::loadTdfAxisBounds(makeTdf("bounds_max.tdf", dShaped() + "UPDATE GlobalMetadata SET Value = '100000000' WHERE Key = 'DigitizerNumSamples';"
                                                "UPDATE GlobalMetadata SET Value = '100000' WHERE Key = 'MzAcqRangeUpper';"), bd, why) &&
          bd.n_bins == 100000000 && bd.mz_hi == 1e5, "DigitizerNumSamples 1e8, MzAcqRangeUpper 1e5: %s", why.c_str());
  CHECK(!diaspextract::loadTdfAxisBounds(makeTdf("bounds_bins.tdf", dShaped() + "UPDATE GlobalMetadata SET Value = '100000001' WHERE Key = 'DigitizerNumSamples';"), bd, why) &&
          has(why, "DigitizerNumSamples"), "DigitizerNumSamples 1e8 + 1: %s", why.c_str());
  CHECK(!diaspextract::loadTdfAxisBounds(makeTdf("bounds_mz.tdf", dShaped() + "UPDATE GlobalMetadata SET Value = '100000.001' WHERE Key = 'MzAcqRangeUpper';"), bd, why) &&
          has(why, "not usable"), "MzAcqRangeUpper above 1e5: %s", why.c_str());

  // a negative C2 (PXD029836 run 1418's row, timsTOF Pro 2) loads; one whose root stops existing inside the domain is
  // refused, naming the TOF index; a frame whose T1 lies outside the checked span is refused for a negative C2 only
  const std::string pxd = "UPDATE MzCalibration SET DigitizerTimebase = 0.19999999999999998, DigitizerDelay = 25779.8, T1 = 25.618823611312585,"
                          "  dC1 = 20, C0 = 314.1341896383674, C1 = 154199.217568937, C2 = -0.0010422663259118895;"
                          "UPDATE Frames SET T1 = 25.614796222167577 - Id * 1e-6;";
  CHECK(diaspextract::loadTdfCalibration(makeTdf("cal_c2neg.tdf", dShaped() + pxd), cal, t1, why) && cal.C2 == -0.0010422663259118895 &&
          cal.digitizer_timebase == 0.19999999999999998 && t1.size() == 11 && t1[10] == 25.614796222167577 - 10 * 1e-6,
        "a negative C2 inside its domain: %s", why.c_str());
  CHECK(!diaspextract::loadTdfCalibration(makeTdf("cal_c2deep.tdf", dShaped() + pxd + "UPDATE MzCalibration SET C2 = -0.1;"), cal, t1, why) &&
          has(why, "C2 -0.1 < 0 outside its verified domain") && has(why, "stops existing at TOF index 7.93471e+07"),
        "C2 -0.1: %s", why.c_str());
  CHECK(!diaspextract::loadTdfCalibration(makeTdf("cal_c2neg_hot.tdf", dShaped() + pxd + "UPDATE Frames SET T1 = 1100 WHERE Id = 7;"), cal, t1, why) &&
          has(why, "frame 7: Frames.T1 1100 is not within 1000 K of the MzCalibration T1 25.6188"), "a negative C2, frame 7 at 1100 C: %s", why.c_str());
  CHECK(diaspextract::loadTdfCalibration(makeTdf("cal_c2neg_null.tdf", dShaped() + pxd + "UPDATE Frames SET T1 = NULL WHERE Id = 7;"), cal, t1, why) &&
          t1[7] == cal.T1_ref, "a negative C2, a NULL T1 is the reference temperature: %s", why.c_str());
  CHECK(diaspextract::loadTdfCalibration(makeTdf("cal_c2pos_hot.tdf", dShaped() + "UPDATE Frames SET T1 = 1100 WHERE Id = 7;"), cal, t1, why) && t1[7] == 1100.0,
        "a positive C2 reads frame 7 at 1100 C as before: %s", why.c_str());
}

int main()
{
  g_dir = std::filesystem::temp_directory_path() /
          ("spx_test_tdf_load_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directories(g_dir);
  std::printf("SQLite %s, fixtures in %s\n", sqlite3_libversion(), g_dir.string().c_str());
  int before = g_fail;
  testWindows();
  if (g_fail == before) std::printf("OK  DiaFrameMsMsWindows: whole (24, 3000), absent = no gate, any case; a missing column, a NULL scan, a corrupt page refused\n");
  before = g_fail;
  testFrames();
  if (g_fail == before) std::printf("OK  Frames and TimsCalibration: typed and ranged before narrowing; per-row scan need; corrupt leaf refused\n");
  before = g_fail;
  testAccumulationTime();
  if (g_fail == before) std::printf("OK  AccumulationTime: whole C-locale decimal > 0; above 100 ms the setup refusal\n");
  before = g_fail;
  testCalibrationAndBounds();
  if (g_fail == before) std::printf("OK  loadTdfCalibration and loadTdfAxisBounds: whole reads, schema-decided fallback; a negative C2 inside its domain, frame T1 span\n");
  std::error_code ec;
  std::filesystem::remove_all(g_dir, ec);
  if (g_fail) std::fprintf(stderr, "%d check(s) failed\n", g_fail);
  return g_fail ? 1 : 0;
}
