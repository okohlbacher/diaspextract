// analysis.tdf readers (sqlite3 C API): the m/z calibration, its axis bounds, and the dnoise port's inputs.
// Kept out of TdfMzCalibration.h, which is sqlite-free and compiled by the golden test on every CI platform.
#pragma once
#include <OpenMS/FORMAT/TdfMzCalibration.h>   // the copy the OpenMS patch installs; same header as src/TdfMzCalibration.h
#include <sqlite3.h>
#include "DnoiseMs1.h"   // [dnoise] RunMeta, WindowRow: the MS1 port's gate inputs; recover::opentimsCorrection
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace diaspextract
{
// ---------------------------------------------------------------------------------------------
// sqlite3 directly, not SQLiteCpp: SQLiteCpp is an OpenMS in-tree extern absent from an installed OpenMS,
// and a missing SQLite must be a compile error, never a silent detector switch (docs/BASELINE.md).
// Every read is whole or refused with sqlite3_errmsg: each prepare is checked and each row loop must end in
// SQLITE_DONE (an error merely ends a loop, so a partial read would look complete). Optional tables and
// columns are read from the schema, never inferred from a failed query.
// ---------------------------------------------------------------------------------------------
/// Run-level flight-time axis bounds from GlobalMetadata: DigitizerNumSamples (every bin index is below it)
/// and the acquisition m/z range. A NULL, trailing garbage, non-finite or unordered value refuses.
struct TdfAxisBounds { long n_bins = 0; double mz_lo = 0.0, mz_hi = 0.0; };
inline bool loadTdfAxisBounds(const std::string& tdf, TdfAxisBounds& out, std::string& why)
{
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(tdf.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
  { why = std::string("cannot open ") + tdf + ": " + (db ? sqlite3_errmsg(db) : "out of memory"); sqlite3_close(db); return false; }
  sqlite3_stmt* st = nullptr; out = TdfAxisBounds{}; why.clear();
  auto num = [&](const unsigned char* txt, int nbytes, const char* key, double& v) {
    if (!txt) { why = std::string(key) + " is NULL"; return false; }
    const std::string s(reinterpret_cast<const char*>(txt), (size_t)nbytes);   // the WHOLE value (an embedded NUL is not an end), in the C locale
    if (!spx::dnoise::detail::parseF64(s, v)) { why = std::string(key) + " is not a finite number: '" + s + "'"; return false; }
    return true; };
  if (sqlite3_prepare_v2(db, "SELECT Key, Value FROM GlobalMetadata WHERE Key IN ('MzAcqRangeLower','MzAcqRangeUpper','DigitizerNumSamples')", -1, &st, nullptr) != SQLITE_OK)
  { why = std::string("cannot read GlobalMetadata: ") + sqlite3_errmsg(db); sqlite3_finalize(st); sqlite3_close(db); return false; }
  bool ok = true; bool have_lo = false, have_hi = false, have_n = false;
  int rc = SQLITE_DONE;
  while (ok && (rc = sqlite3_step(st)) == SQLITE_ROW)
  {
    const unsigned char* kt = sqlite3_column_text(st, 0); if (!kt) continue;
    const std::string k = reinterpret_cast<const char*>(kt); double v = 0.0;
    { bool& have = (k == "MzAcqRangeLower") ? have_lo : (k == "MzAcqRangeUpper") ? have_hi : have_n;
      if (have) { why = k + " appears twice in GlobalMetadata"; ok = false; break; }   // a repeated key would let the later row win silently
      have = true; }
    if (!num(sqlite3_column_text(st, 1), sqlite3_column_bytes(st, 1), k.c_str(), v)) { ok = false; break; }
    if (k == "MzAcqRangeLower") out.mz_lo = v;
    else if (k == "MzAcqRangeUpper") out.mz_hi = v;
    else
    {
      // an integral count in [2, 1e8], checked BEFORE the narrowing cast (1e100 -> long is UB)
      if (!(v >= 2.0 && v <= 100000000.0 && std::trunc(v) == v)) { why = "DigitizerNumSamples is not an integer in [2, 1e8]: " + std::to_string(v); ok = false; break; }
      out.n_bins = (long)v;
    }
  }
  if (ok && rc != SQLITE_DONE) { why = std::string("reading GlobalMetadata: ") + sqlite3_errmsg(db); ok = false; }
  sqlite3_finalize(st); sqlite3_close(db);
  if (!ok) return false;
  if (!(have_lo && have_hi && have_n)) { why = "GlobalMetadata lacks MzAcqRangeLower/MzAcqRangeUpper/DigitizerNumSamples (" + std::to_string(have_lo + have_hi + have_n) + " of 3)"; return false; }
  // plausibility: a TIMS acquisition ends far below 1e5 m/z; a corrupt but finite value beyond
  // that would drive the flight-time conversion past the axis (review 3a, round 2)
  if (!(out.mz_lo > 0.0) || !(out.mz_hi > out.mz_lo) || !(out.mz_hi <= 1.0e5))
  { why = "GlobalMetadata bounds are not usable: MzAcqRange " + std::to_string(out.mz_lo) + "-" + std::to_string(out.mz_hi) + ", DigitizerNumSamples " + std::to_string(out.n_bins); return false; }
  return true;
}

/// The MzCalibration row the FRAMES reference (Frames.MzCalibration; PXD017703 files carry two rows), refused
/// when frames reference more than one. t1_by_frame[id] = Frames.T1 (T1_ref when NULL); index 0 is never a frame.
inline bool loadTdfCalibration(const std::string& tdf, TdfMzCalibration& cal,
                               std::vector<double>& t1_by_frame, std::string& why)
{
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(tdf.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
  { why = std::string("cannot open ") + tdf + ": " + (db ? sqlite3_errmsg(db) : "out of memory"); sqlite3_close(db); return false; }
  sqlite3_stmt* st = nullptr;
  auto fail = [&](const std::string& w) { why = w; sqlite3_finalize(st); sqlite3_close(db); return false; };
  auto sqlfail = [&](const std::string& what) { return fail(what + ": " + sqlite3_errmsg(db)); };   // errmsg is read before fail() finalizes
  int rc = SQLITE_OK;
  bool frames_have_cal = false;   // PRAGMA table_xinfo: one row per column, generated ones too (table_info skips them), the name in column 1, matched ASCII case-insensitively as sqlite resolves it
  if (sqlite3_prepare_v2(db, "PRAGMA table_xinfo(Frames)", -1, &st, nullptr) != SQLITE_OK) return sqlfail("cannot read the Frames schema");
  while ((rc = sqlite3_step(st)) == SQLITE_ROW)
  { const unsigned char* c = sqlite3_column_text(st, 1); frames_have_cal = frames_have_cal || (c && sqlite3_stricmp(reinterpret_cast<const char*>(c), "MzCalibration") == 0); }
  if (rc != SQLITE_DONE) return sqlfail("reading the Frames schema");
  sqlite3_finalize(st); st = nullptr;
  long long cal_id = -1; int ncal = 0;
  if (frames_have_cal)
  {
    if (sqlite3_prepare_v2(db, "SELECT DISTINCT MzCalibration FROM Frames", -1, &st, nullptr) != SQLITE_OK) return sqlfail("cannot read Frames.MzCalibration");
    while ((rc = sqlite3_step(st)) == SQLITE_ROW) { cal_id = sqlite3_column_int64(st, 0); ++ncal; }
    if (rc != SQLITE_DONE) return sqlfail("reading Frames.MzCalibration");
    sqlite3_finalize(st); st = nullptr;
    if (ncal != 1) return fail("frames reference " + std::to_string(ncal) + " distinct MzCalibration rows");
  }
  else
  {
    // Frames has no MzCalibration column (a minimal or very old tdf): fall back to the single-row
    // rule rather than refusing -- the reference is then unambiguous only if there IS one row.
    if (sqlite3_prepare_v2(db, "SELECT Id FROM MzCalibration", -1, &st, nullptr) != SQLITE_OK) return sqlfail("Frames has no MzCalibration column, and MzCalibration cannot be read");
    while ((rc = sqlite3_step(st)) == SQLITE_ROW) { cal_id = sqlite3_column_int64(st, 0); ++ncal; }
    if (rc != SQLITE_DONE) return sqlfail("reading MzCalibration");
    sqlite3_finalize(st); st = nullptr;
    if (ncal != 1) return fail("Frames has no MzCalibration column and MzCalibration has " + std::to_string(ncal) + " rows (need exactly 1)");
  }
  if (sqlite3_prepare_v2(db, "SELECT ModelType, DigitizerTimebase, DigitizerDelay, C0, C1, C2, T1, dC1, dC2, C3, C4 "
                             "FROM MzCalibration WHERE Id = ?", -1, &st, nullptr) != SQLITE_OK)
    return sqlfail("cannot read MzCalibration");
  rc = sqlite3_bind_int64(st, 1, cal_id);
  if (rc == SQLITE_OK) rc = sqlite3_step(st);
  if (rc == SQLITE_DONE) return fail("no MzCalibration row with Id " + std::to_string(cal_id));
  if (rc != SQLITE_ROW) return sqlfail("reading MzCalibration row " + std::to_string(cal_id));
  if (sqlite3_column_type(st, 5) == SQLITE_NULL) return fail("MzCalibration C2 is NULL (missing, not zero)");
  cal.model_type = sqlite3_column_int(st, 0);
  cal.digitizer_timebase = sqlite3_column_double(st, 1); cal.digitizer_delay = sqlite3_column_double(st, 2);
  cal.C0 = sqlite3_column_double(st, 3); cal.C1 = sqlite3_column_double(st, 4);
  cal.C2 = sqlite3_column_double(st, 5);
  cal.T1_ref = sqlite3_column_double(st, 6); cal.dC1 = sqlite3_column_double(st, 7); cal.dC2 = sqlite3_column_double(st, 8);
  cal.C3 = sqlite3_column_double(st, 9); cal.C4 = sqlite3_column_double(st, 10);
  sqlite3_finalize(st); st = nullptr;
  if (!cal.isSupported()) return fail(cal.unsupportedReason());
  t1_by_frame.clear();
  if (sqlite3_prepare_v2(db, "SELECT Id, T1 FROM Frames", -1, &st, nullptr) != SQLITE_OK) return sqlfail("cannot read Frames.T1");
  while ((rc = sqlite3_step(st)) == SQLITE_ROW)
  {
    const long long id = sqlite3_column_int64(st, 0);
    if (id < 0 || id > 10000000) return fail("implausible Frames.Id");
    if ((size_t)id >= t1_by_frame.size()) t1_by_frame.resize((size_t)id + 1, cal.T1_ref);
    t1_by_frame[(size_t)id] = (sqlite3_column_type(st, 1) == SQLITE_NULL) ? cal.T1_ref : sqlite3_column_double(st, 1);
  }
  if (rc != SQLITE_DONE) return sqlfail("reading Frames.T1");
  sqlite3_finalize(st); sqlite3_close(db);
  if (t1_by_frame.size() < 2) { why = "no Frames rows"; return false; }
  return true;
}

// ---------------------------------------------------------------------------------------------
// [dnoise] What the MS1 port (DnoiseMs1.h) reads besides the m/z calibration:
//  * the gate's inputs exactly as dnoise v0.1.0 reads them: six GlobalMetadata values as stored TEXT,
//    MAX(NumScans) over ALL frames, and the DiaFrameMsMsWindows rows (RunMeta / WindowRow);
//  * what the Bruker loader made of each raw point, so the tool can invert it: per frame MsMsType, NumScans, NumPeaks,
//    opentims' intensity correction (recover::opentimsCorrection of sqlite's text of AccumulationTime) and the
//    TimsCalibration row of the rational 1/K0 model. That text is the loader's only where both share one SQLite, as
//    in-tree; a standalone build links its own, and SQLite versions render a REAL that 15 digits do not round-trip
//    differently. In a brute force over every raw count below 2^24 the few-ulp corr change refused points but never
//    recovered a wrong count, and every AccumulationTime seen (49 .d files) round-trips at 15 digits.
// Integers are type- and range-checked before narrowing: sqlite3_column_int64 converts TEXT/REAL silently, a cast wraps.
// Refuses wherever tryCreateRationalConverter would NOT install that model: the loader's 1/K0 is then linear.
// ---------------------------------------------------------------------------------------------
struct TdfDnoiseInputs
{
  spx::dnoise::RunMeta meta;
  std::vector<spx::dnoise::WindowRow> windows;
  /// corr: recover::opentimsCorrection of the frame's AccumulationTime (finite, > 0; whether it is invertible is the caller's check)
  struct Frame { int msms_type = -1; std::uint32_t num_scans = 0, num_peaks = 0; double corr = 0.0; std::uint32_t cal = 0; };
  std::vector<Frame> frames;                     ///< by Frames.Id (msms_type -1: no such frame); cal indexes tims_cal
  std::vector<std::array<double, 10>> tims_cal;  ///< C0..C9 of each TimsCalibration row
  std::vector<std::uint32_t> ms1_max_scans;      ///< per tims_cal row: the largest NumScans of an MS1 frame referencing it (0: none does)
};
inline bool loadTdfDnoise(const std::string& tdf, TdfDnoiseInputs& out, std::string& why)
{
  out = TdfDnoiseInputs{}; why.clear();
  sqlite3* db = nullptr;
  if (sqlite3_open_v2(tdf.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK)
  { why = std::string("cannot open ") + tdf + ": " + (db ? sqlite3_errmsg(db) : "out of memory"); sqlite3_close(db); return false; }
  sqlite3_stmt* st = nullptr;
  auto fail = [&](const std::string& w) { why = w; sqlite3_finalize(st); sqlite3_close(db); return false; };
  auto sqlfail = [&](const std::string& what, const char* tail = "") { return fail(what + ": " + sqlite3_errmsg(db) + tail); };   // errmsg is read before fail() finalizes
  auto int_in = [&](int col, long long lo, long long hi, long long& v) {   // a stored INTEGER in [lo, hi]
    if (sqlite3_column_type(st, col) != SQLITE_INTEGER) return false;
    v = sqlite3_column_int64(st, col);
    return v >= lo && v <= hi; };
  auto cell = [&](int col) {   // a stored value as a refusal names it: NULL spelled out, TEXT quoted
    const int type = sqlite3_column_type(st, col);
    if (type == SQLITE_NULL) return std::string("NULL");
    const unsigned char* x = sqlite3_column_text(st, col);
    const std::string s = x ? std::string(reinterpret_cast<const char*>(x), (size_t)sqlite3_column_bytes(st, col)) : std::string();
    return type == SQLITE_INTEGER || type == SQLITE_FLOAT ? s : "'" + s + "'"; };
  const char* const linear = " (the loader converts 1/K0 with the linear GlobalMetadata model)";
  int rc = SQLITE_OK;
  if (sqlite3_prepare_v2(db, "SELECT Key, Value FROM GlobalMetadata", -1, &st, nullptr) != SQLITE_OK) return sqlfail("cannot read GlobalMetadata");
  while ((rc = sqlite3_step(st)) == SQLITE_ROW)
  {
    const unsigned char* k = sqlite3_column_text(st, 0);
    if (!k) continue;
    const std::string key = reinterpret_cast<const char*>(k);
    std::optional<std::string>* slot = key == "AcquisitionSoftware" ? &out.meta.acquisition_software
                                     : key == "MzAcqRangeLower" ? &out.meta.mz_acq_range_lower : key == "MzAcqRangeUpper" ? &out.meta.mz_acq_range_upper
                                     : key == "DigitizerNumSamples" ? &out.meta.digitizer_num_samples
                                     : key == "OneOverK0AcqRangeLower" ? &out.meta.one_over_k0_acq_range_lower
                                     : key == "OneOverK0AcqRangeUpper" ? &out.meta.one_over_k0_acq_range_upper : nullptr;
    const unsigned char* v = sqlite3_column_text(st, 1);
    if (slot && v) *slot = std::string(reinterpret_cast<const char*>(v), (size_t)sqlite3_column_bytes(st, 1));
  }
  if (rc != SQLITE_DONE) return sqlfail("reading GlobalMetadata");
  sqlite3_finalize(st); st = nullptr;
  if (sqlite3_prepare_v2(db, "SELECT MAX(NumScans) FROM Frames", -1, &st, nullptr) != SQLITE_OK) return sqlfail("cannot read Frames.NumScans");
  if ((rc = sqlite3_step(st)) != SQLITE_ROW) return sqlfail("reading MAX(Frames.NumScans)");
  { long long v = 0;   // NULL: an empty Frames table, refused below. A TEXT cell anywhere is the MAX: it sorts above INTEGER.
    if (sqlite3_column_type(st, 0) != SQLITE_NULL && !int_in(0, 0, 1000000, v))
      return fail("MAX(Frames.NumScans) " + cell(0) + " is not an INTEGER in [0, 1000000]");
    out.meta.max_num_scans = (std::uint32_t)v; }
  sqlite3_finalize(st); st = nullptr;
  // DiaFrameMsMsWindows is optional (no table: dnoise has no gate). Table names resolve case-insensitively.
  if (sqlite3_prepare_v2(db, "SELECT 1 FROM sqlite_master WHERE type IN ('table', 'view') AND name = 'DiaFrameMsMsWindows' COLLATE NOCASE", -1, &st, nullptr) != SQLITE_OK)
    return sqlfail("cannot read the schema");
  rc = sqlite3_step(st);
  if (rc != SQLITE_ROW && rc != SQLITE_DONE) return sqlfail("reading the schema");
  const bool have_windows = rc == SQLITE_ROW;
  sqlite3_finalize(st); st = nullptr;
  if (have_windows)
  {
    if (sqlite3_prepare_v2(db, "SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth FROM DiaFrameMsMsWindows "
                               "WHERE IsolationMz IS NOT NULL AND IsolationWidth IS NOT NULL", -1, &st, nullptr) != SQLITE_OK)
      return sqlfail("cannot read DiaFrameMsMsWindows");
    while ((rc = sqlite3_step(st)) == SQLITE_ROW)
    {
      if (sqlite3_column_type(st, 0) == SQLITE_NULL || sqlite3_column_type(st, 1) == SQLITE_NULL)   // column_int64 would read scan 0
        return fail("DiaFrameMsMsWindows has a row with a NULL ScanNumBegin or ScanNumEnd");
      out.windows.push_back({sqlite3_column_int64(st, 0), sqlite3_column_int64(st, 1), sqlite3_column_double(st, 2), sqlite3_column_double(st, 3)});
    }
    if (rc != SQLITE_DONE) return sqlfail("reading DiaFrameMsMsWindows");
    sqlite3_finalize(st); st = nullptr;
  }
  std::vector<std::uint32_t> cal_id;
  if (sqlite3_prepare_v2(db, "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration", -1, &st, nullptr) != SQLITE_OK)
    return sqlfail("cannot read TimsCalibration", linear);
  while ((rc = sqlite3_step(st)) == SQLITE_ROW)
  {
    if (sqlite3_column_int(st, 1) != 2) return fail("TimsCalibration ModelType " + std::to_string(sqlite3_column_int(st, 1)) + linear);
    long long id = 0;
    if (!int_in(0, 0, 4294967295LL, id)) return fail("TimsCalibration.Id " + cell(0) + " is not an INTEGER in [0, 4294967295]");
    std::array<double, 10> c{};
    for (int k = 0; k < 10; ++k) c[(size_t)k] = sqlite3_column_double(st, 2 + k);
    cal_id.push_back((std::uint32_t)id); out.tims_cal.push_back(c);
  }
  if (rc != SQLITE_DONE) return sqlfail("reading TimsCalibration");
  sqlite3_finalize(st); st = nullptr;
  if (out.tims_cal.empty()) return fail(std::string("TimsCalibration is empty") + linear);
  out.ms1_max_scans.assign(out.tims_cal.size(), 0);
  if (sqlite3_prepare_v2(db, "SELECT Id, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration FROM Frames", -1, &st, nullptr) != SQLITE_OK)
    return sqlfail("cannot read Frames (Id, MsMsType, NumScans, NumPeaks, AccumulationTime, TimsCalibration)",
                   " (without a TimsCalibration column the loader converts 1/K0 with the linear GlobalMetadata model)");
  while ((rc = sqlite3_step(st)) == SQLITE_ROW)
  {
    long long id = 0, type = 0, scans = 0, peaks = 0, cid = 0;
    if (!int_in(0, 1, 10000000, id)) return fail("implausible Frames.Id " + cell(0));
    auto bad = [&](int col, const char* name, const std::string& want) {
      return fail("frame " + std::to_string(id) + ": Frames." + name + " " + cell(col) + " is not " + want); };
    if (!int_in(1, INT_MIN, INT_MAX, type)) return bad(1, "MsMsType", "an INTEGER in the int range");
    if (!int_in(2, 0, out.meta.max_num_scans, scans)) return bad(2, "NumScans", "an INTEGER in [0, " + std::to_string(out.meta.max_num_scans) + "]");
    if (!int_in(3, 0, 4294967295LL, peaks)) return bad(3, "NumPeaks", "an INTEGER in [0, 4294967295]");
    if (type == 0 && peaks > 0 && scans < 1)
      return fail("frame " + std::to_string(id) + ": an MS1 frame with " + std::to_string(peaks) + " points (Frames.NumPeaks) and no scans (Frames.NumScans 0)");
    const unsigned char* acc = sqlite3_column_text(st, 4);
    if (!acc) return fail("frame " + std::to_string(id) + " has no AccumulationTime");
    double corr = 0.0;
    if (!spx::dnoise::recover::opentimsCorrection(std::string(reinterpret_cast<const char*>(acc), (size_t)sqlite3_column_bytes(st, 4)), corr))
      return fail("frame " + std::to_string(id) + ": Frames.AccumulationTime '" + std::string(reinterpret_cast<const char*>(acc), (size_t)sqlite3_column_bytes(st, 4))
                  + "' is not a positive decimal number");
    if (sqlite3_column_type(st, 5) == SQLITE_NULL) return fail("frame " + std::to_string(id) + " has no TimsCalibration row" + linear);
    if (!int_in(5, 0, 4294967295LL, cid)) return bad(5, "TimsCalibration", "an INTEGER in [0, 4294967295]");
    size_t ci = cal_id.size();
    for (size_t i = 0; i < cal_id.size(); ++i) if (cal_id[i] == (std::uint32_t)cid) ci = i;   // a repeated Id: the last row, as the loader's map keeps it
    if (ci == cal_id.size()) return fail("frame " + std::to_string(id) + " references an unknown TimsCalibration row" + linear);
    if ((size_t)id >= out.frames.size()) out.frames.resize((size_t)id + 1);
    TdfDnoiseInputs::Frame& f = out.frames[(size_t)id];
    f.msms_type = (int)type;
    f.num_scans = (std::uint32_t)scans;
    f.num_peaks = (std::uint32_t)peaks;
    f.corr = corr;
    f.cal = (std::uint32_t)ci;
    if (f.msms_type == 0) out.ms1_max_scans[ci] = std::max(out.ms1_max_scans[ci], f.num_scans);
  }
  if (rc != SQLITE_DONE) return sqlfail("reading Frames");
  sqlite3_finalize(st); sqlite3_close(db);
  if (out.frames.size() < 2) { why = "no Frames rows"; return false; }
  return true;
}
} // namespace diaspextract
