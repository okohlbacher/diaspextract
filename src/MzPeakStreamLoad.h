// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

// DIAspeXtract: streaming diaPASEF load from an .mzpeak archive through the mzPeak C++ library
// (github.com/OpenMS/mzpeak, fork okohlbacher/mzpeak-openms). Mirrors BrukerTimsFile::
// loadDIAStreaming's contract so PickCompactConsumer sees the same stream: every MS1 spectrum
// first (frame order), then one MSSpectrum per (MS2 frame, isolation window) holding only the
// peaks inside that window's 1/K0 band, handed over serially in frame order.
//
// mzPeak stores a diaPASEF frame as ONE spectrum carrying N precursors (one per window, each
// with an ion-mobility band) and a per-peak ion-mobility array; ims-compact archives carry raw
// TOF and the library reconstructs m/z with the archive's own two-point transform. Peaks arrive
// mobility-major, so each per-window spectrum is m/z-sorted here.
//
// ponytail: decode is parallel over contiguous frame ranges with one Index per thread (decodes
// serialise on a reader's mutex); the hand-off is serial. No caching beyond the library's own.
#pragma once
#ifdef DIASPEXTRACT_WITH_MZPEAK

  #include "TdfLoad.h"
  #include <mzpeak/open.h>
  #include <mzpeak/index.h>
  #include <mzpeak/spectra.h>
  #include <mzpeak/spectrum.h>
  #include <mzpeak/spectrum_metadata.h>

  #include <zip.h>
  #include <zlib.h>
  #include <arrow/api.h>
  #include <parquet/file_reader.h>
  #include <parquet/arrow/schema.h>
  #include <cstdio>
  #include <cstdlib>
  #ifdef _WIN32
    #include <fcntl.h> // [win] mkstemp's stand-in: _mktemp_s + an exclusive _sopen_s in %TEMP%
    #include <io.h>
    #include <share.h>
    #include <sys/stat.h>
    #include <filesystem>
  #else
    #include <unistd.h>
  #endif

  #include <OpenMS/KERNEL/MSSpectrum.h>
  #include <OpenMS/METADATA/Precursor.h>
  #include <OpenMS/CONCEPT/Exception.h>
  #include <OpenMS/IONMOBILITY/IMDataConverter.h>
  #include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>

  #include <algorithm>
  #include <cmath>
  #include <exception>
  #include <limits>
  #include <memory>
  #include <string>
  #include <vector>
  #include <cstring>
  #include <stdexcept>
  #include <omp.h>

namespace DIAspeXtract
{
/**
  @brief One isolation window of an mzPeak spectrum with its ion-mobility band.

  @p lo, @p hi are the isolation window's m/z bounds, @p im_lo, @p im_hi its 1/K0 band. An MS1 spectrum carries one
  window with lo == hi == 0 and an unbounded mobility band.
*/
struct MzPeakWin
{ double lo, hi, im_lo, im_hi; };

/**
  @brief A missing, corrupt or unsupported part of an .mzpeak archive.

  The streaming load rethrows it instead of falling back to a whole-run read: the archive itself is unusable.
*/
struct MzPeakError : OpenMS::Exception::InvalidValue
{ using OpenMS::Exception::InvalidValue::InvalidValue; };

/**
  @brief Calibration provenance of the last mzPeak load (the exact TDF model, or the archive's two-point transform), stamped as spx:mz_calibration.

  @return A reference to the process-wide provenance string ("unset" before any load).
*/
inline std::string& lastMzPeakCalibration()
{
  static std::string s = "unset";
  return s;
}


/**
  @brief Exact TOF->m/z for an ims-compact archive.

  Exact TOF->m/z for an ims-compact archive: (c0, c1) from the peaks table's column metadata, the
  MzCalibration row and Frames.T1 from the embedded vendor tdf. Throws on anything missing (fail closed).
*/
struct MzPeakExactMz
{
  double c0 = 0, c1 = 0;              ///< two-point transform parameters: sqrt(m/z) = c0 + c1 * tof
  diaspextract::TdfMzCalibration cal; ///< MzCalibration model from the vendor tdf
  std::vector<double> t1_by_frame;    ///< index = Frames.Id
  double acq_lo = std::numeric_limits<double>::quiet_NaN(), acq_hi = std::numeric_limits<double>::quiet_NaN(); ///< GlobalMetadata MzAcqRange
  long n_bins = 0;        ///< number of TOF bins of the digitizer (0 when unknown)
  std::string bounds_why; ///< DigitizerNumSamples, or why not

  /**
    @brief Read a byte range of one zip member.

    @p from < 0 counts from the member's end; @p len < 0 reads to it.

    @param[in] z The open archive
    @param[in] name Member name
    @param[in] from First byte to read (negative: offset from the member's end)
    @param[in] len Number of bytes to read (negative: up to the member's end)
    @return The bytes read
    @throws MzPeakError if the member is missing, cannot be opened or seeked, the range exceeds it, or the read is short
  */
  static std::vector<char> readMember(zip_t* z, const char* name, zip_int64_t from = 0, zip_int64_t len = -1)
  {
    zip_int64_t idx = zip_name_locate(z, name, 0);
    if (idx < 0) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: member missing: ") + name), ""); }
    zip_stat_t st;
    if (zip_stat_index(z, (zip_uint64_t)idx, 0, &st) != 0 || ! (st.valid & ZIP_STAT_SIZE))
    {
      throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: cannot stat member ") + name), "");
    }
    const zip_int64_t size = (zip_int64_t)st.size;
    if (from < 0) { from += size; }
    if (len < 0) { len = size - from; }
    if (from < 0 || len < 0 || from > size || len > size - from) // sizes from the archive: never allocate past the member
    {
      throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                        std::string(std::string("mzpeak: a read beyond the member ") + name + " (corrupt archive?)"), "");
    }
    zip_file_t* f = zip_fopen_index(z, (zip_uint64_t)idx, 0);
    if (! f) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: cannot open member ") + name), ""); }
    if (from > 0 && zip_fseek(f, from, SEEK_SET) != 0)
    {
      zip_fclose(f);
      throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                        std::string(std::string("mzpeak: cannot seek in ") + name + " (member not stored?)"), "");
    }
    // the advertised size is not trusted for the allocation: the buffer grows only as data actually arrive (64 MiB at a
    // time), so a member claiming more than the archive holds ends in a short read, not in a huge allocation
    constexpr zip_int64_t kChunk = zip_int64_t(1) << 26;
    std::vector<char> buf;
    zip_int64_t got = 0;
    while (got < len)
    {
      buf.resize((OpenMS::Size)std::min(len, got + kChunk));
      zip_int64_t n = zip_fread(f, buf.data() + got, (zip_uint64_t)((zip_int64_t)buf.size() - got));
      if (n <= 0) { break; }
      got += n;
    }
    zip_fclose(f);
    if (got != len) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: short read of ") + name), ""); }
    return buf;
  }

  /**
    @brief Load the transform parameters and the vendor calibration of an archive.

    Reads (c0, c1) from the peaks parquet footer, then the MzCalibration row, Frames.T1 and the axis bounds from the
    embedded vendor/analysis.tdf.gz (or the sidecar named by DIASPEXTRACT_MZPEAK_TDF), gunzipped to a temporary file.

    @param[in] path Path of the .mzpeak archive
    @throws MzPeakError if the archive, the transform parameters or the vendor calibration cannot be read
  */
  explicit MzPeakExactMz(const std::string& path)
  {
    int err = 0;
    zip_t* z = zip_open(path.c_str(), ZIP_RDONLY, &err);
    if (! z) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: zip_open failed on " + path), ""); }
    try
    {
      // (1) two-point transform params from the parquet footer of spectra_peaks.parquet (stored member)
      std::vector<char> tail = readMember(z, "spectra_peaks.parquet", -8, 8);
      uint32_t flen;
      std::memcpy(&flen, tail.data(), 4);
      if (std::string(tail.data() + 4, 4) != "PAR1")
      {
        throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: peaks member is not parquet"), "");
      }
      // readMember refuses a footer length beyond the member, before allocating it
      std::vector<char> footer = readMember(z, "spectra_peaks.parquet", -8 - (zip_int64_t)flen, flen);
      uint32_t flen2 = flen;
      std::shared_ptr<parquet::FileMetaData> md = parquet::FileMetaData::Make(footer.data(), &flen2);
      std::shared_ptr<arrow::Schema> schema;
      parquet::ArrowReaderProperties props;
      PARQUET_THROW_NOT_OK(parquet::arrow::FromParquetSchema(md->schema(), props, md->key_value_metadata(), &schema));
      auto point = schema->GetFieldByName("point");
      if (! point) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: no 'point' column"), ""); }
      const auto point_t = std::dynamic_pointer_cast<arrow::StructType>(point->type());
      if (! point_t) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: the 'point' column is not a struct"), ""); }
      auto tof = point_t->GetFieldByName("tof");
      if (! tof || ! tof->metadata())
      {
        throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: no tof column / metadata (not an ims-compact archive?)"),
                          "");
      }
      auto tp = tof->metadata()->Get("mzpeak:transform_params");
      if (! tp.ok())
      {
        throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: tof column has no mzpeak:transform_params"), "");
      }
      if (std::sscanf(tp->c_str(), "%lf,%lf", &c0, &c1) != 2 || ! (c1 > 0))
      {
        throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: bad transform_params " + *tp), "");
      }

      // (2) the vendor tdf: gunzip vendor/analysis.tdf.gz to a temp file, read MzCalibration + Frames.T1
      // --no-vendor archives have no embedded tdf: accept a sidecar via DIASPEXTRACT_MZPEAK_TDF, an analysis.tdf.gz or
      // the plain analysis.tdf (an SQLite file is read in place)
      std::vector<char> gz;
      std::string tdf_path;
      const char* side = std::getenv("DIASPEXTRACT_MZPEAK_TDF");
      if (side)
      {
        FILE* sf = std::fopen(side, "rb");
        if (! sf)
        {
          throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                            std::string(std::string("mzpeak: cannot read DIASPEXTRACT_MZPEAK_TDF ") + side), "");
        }
        char b[1 << 16];
        OpenMS::Size n;
        while ((n = std::fread(b, 1, sizeof b, sf)) > 0)
        {
          gz.insert(gz.end(), b, b + n);
        }
        std::fclose(sf);
        if (gz.size() >= 16 && std::memcmp(gz.data(), "SQLite format 3", 16) == 0)
        {
          tdf_path = side;
          std::vector<char>().swap(gz);
        }
      }
      else
      {
        gz = readMember(z, "vendor/analysis.tdf.gz");
      }
  #ifdef _WIN32
      const std::string tmp_s = (std::filesystem::temp_directory_path() / "spx_mzpeak_tdf_XXXXXX").string();
      std::vector<char> tmp_b(tmp_s.begin(), tmp_s.end());
      tmp_b.push_back('\0');
      char* tmpl = tmp_b.data();
      int fd = -1;
      if (_mktemp_s(tmpl, tmp_b.size()) == 0) { _sopen_s(&fd, tmpl, _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY, _SH_DENYNO, _S_IREAD | _S_IWRITE); }
      if (fd < 0) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: mkstemp"), ""); }
      _close(fd);
  #else
      char tmpl[] = "/tmp/spx_mzpeak_tdf_XXXXXX";
      int fd = mkstemp(tmpl);
      if (fd < 0) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: mkstemp"), ""); }
      close(fd);
  #endif
      struct Rm
      {
        const char* p;
        ~Rm()
        { std::remove(p); }
      } rm {tmpl}; // the uncompressed tdf goes on every way out
      if (tdf_path.empty())
      {
        z_stream zs {};
        FILE* out = std::fopen(tmpl, "wb");
        if (! out) { throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: cannot write ") + tmpl), ""); }
        if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK)
        {
          std::fclose(out);
          throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: gunzip of analysis.tdf.gz failed"), "");
        }
        // a real analysis.tdf is well below 1 GB; a hostile gzip must not fill the temp directory
        const unsigned long long kMaxTdf = 1ull << 32;
        // the output is counted here, not read from zs.total_out: that is a uLong, 32 bits on Windows, wrapping below the cap
        zs.next_in = (Bytef*)gz.data();
        zs.avail_in = (uInt)gz.size();
        std::vector<char> ob(1 << 20);
        int rc;
        unsigned long long total = 0;
        do
        {
          zs.next_out = (Bytef*)ob.data();
          zs.avail_out = (uInt)ob.size();
          rc = inflate(&zs, Z_NO_FLUSH);
          const OpenMS::Size got = ob.size() - zs.avail_out;
          total += got;
          if ((rc != Z_OK && rc != Z_STREAM_END) || total > kMaxTdf || std::fwrite(ob.data(), 1, got, out) != got)
          {
            std::fclose(out);
            inflateEnd(&zs);
            throw MzPeakError(
              __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              std::string(total > kMaxTdf ? "mzpeak: analysis.tdf.gz expands beyond 4 GiB" : "mzpeak: gunzip of analysis.tdf.gz failed"), "");
          }
        } while (rc != Z_STREAM_END);
        inflateEnd(&zs);
        if (std::fclose(out) != 0)
        {
          throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string(std::string("mzpeak: cannot write ") + tmpl), "");
        }
        tdf_path = tmpl;
      }
      {
        std::string why;
        if (! DIAspeXtract::loadTdfCalibration(tdf_path, cal, t1_by_frame, why))
        {
          throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzpeak: " + why), "");
        }
        std::string w2;
        DIAspeXtract::TdfAxisBounds bd;
        if (DIAspeXtract::loadTdfAxisBounds(tdf_path, bd, w2))
        {
          acq_lo = bd.mz_lo;
          acq_hi = bd.mz_hi;
          n_bins = bd.n_bins;
        }
        else
        {
          bounds_why = w2;
        }
      }
    }
    catch (...)
    {
      zip_close(z);
      throw;
    }
    zip_close(z);
  }

  /**
    @brief Recover the integer tof from the archive's two-point m/z and re-apply the exact model.

    @param[in,out] mz m/z values of one frame, replaced by the exact values
    @param[in] frame_id Vendor frame id (Frames.Id) selecting T1
    @throws OpenMS::Exception::InvalidValue if the frame has no T1, an m/z is not positive, or the result is not finite
  */
  void exact(std::vector<double>& mz, long frame_id) const
  {
    // Fail rather than fall back to the reference T1: substituting it for an unknown frame stamps
    // masses calibrated at the wrong temperature as "exact", which is silent wrongness the
    // provenance record would then vouch for.
    if (frame_id < 0 || (OpenMS::Size)frame_id >= t1_by_frame.size() || std::isnan(t1_by_frame[(OpenMS::Size)frame_id])) // NaN: a gap in the table
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak exact calibration: no tdf T1 for frame",
                                            std::to_string(frame_id));
    }
    const double b = cal.frameFactor(t1_by_frame[(OpenMS::Size)frame_id]);
    for (double& v : mz)
    {
      if (! std::isfinite(v) || v <= 0.0)
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak exact calibration: non-positive archive m/z",
                                              std::to_string(v));
      }
      const double tof = std::llround((std::sqrt(v) - c0) / c1);
      const double e = cal.tofToMz(tof, b);
      if (! std::isfinite(e))
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak exact calibration produced a non-finite m/z at tof",
                                              std::to_string(tof));
      }
      v = e;
    }
  }
};

/**
  @brief One decoded frame split per window (index-aligned with `wins` of that frame).

  Keeps the peaks inside each window's 1/K0 band with positive intensity, sorted by m/z with the ion-mobility array
  carried along. Windows without peaks stay empty.

  @param[in] sp The decoded frame
  @param[in] rt Retention time of the frame
  @param[in] frame_idx Spectrum index in the archive (written to the native ID)
  @param[in] wins Isolation windows of the frame
  @param[out] out One spectrum per window
  @param[in] exact Exact calibration, or nullptr to keep the archive's m/z
  @param[in] frame_id Vendor frame id, or -1 when unknown
  @throws OpenMS::Exception::InvalidValue if the array lengths disagree or the exact calibration fails
*/
inline void frameToSpectra(const MzPeak::Spectrum& sp,
                           double rt,
                           OpenMS::Size frame_idx,
                           const std::vector<MzPeakWin>& wins,
                           std::vector<OpenMS::MSSpectrum>& out,
                           const MzPeakExactMz* exact,
                           long frame_id)
{
  std::vector<double> mz_exact;
  if (exact)
  {
    mz_exact = sp.mz();
    exact->exact(mz_exact, frame_id);
  }
  const std::vector<double>& mz = exact ? mz_exact : sp.mz();
  const std::vector<float>& in = sp.intensity();
  const std::vector<double>& im = sp.ion_mobility_array();
  // Array lengths must agree. Truncating to the shorter one turns a malformed archive into a
  // smaller frame that looks valid, and the peaks that vanish are never reported anywhere.
  if (mz.size() != in.size())
  {
    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak frame has mismatched m/z and intensity array lengths",
                                          std::to_string(mz.size()) + " vs " + std::to_string(in.size()));
  }
  const OpenMS::Size n = mz.size();
  if (im.size() < n)
  {
    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak frame has a short ion-mobility array",
                                          std::to_string(im.size()));
  }
  out.assign(wins.size(), OpenMS::MSSpectrum());
  for (OpenMS::Size w = 0; w < wins.size(); ++w)
  {
    const MzPeakWin& win = wins[w];
    OpenMS::MSSpectrum& spec = out[w];
    spec.setRT(rt);
    spec.setMSLevel(wins.size() == 1 && ! (win.lo < win.hi) ? 1 : 2);
    spec.setDriftTimeUnit(OpenMS::DriftTimeUnit::VSSC);
    // "frame=<Id>" keys the per-frame calibration (without it the integer detector falls back); "mzpeak=" is the archive index.
    spec.setNativeID("mzpeak=" + std::to_string(frame_idx) + (frame_id >= 0 ? " frame=" + std::to_string(frame_id) : std::string())
                     + " window=" + std::to_string(w));
    if (spec.getMSLevel() == 2)
    {
      OpenMS::Precursor prec;
      prec.setMZ(0.5 * (win.lo + win.hi));
      prec.setIsolationWindowLowerOffset(0.5 * (win.hi - win.lo));
      prec.setIsolationWindowUpperOffset(0.5 * (win.hi - win.lo));
      spec.setPrecursors({prec});
      spec.setMetaValue("ion mobility lower limit", win.im_lo);
      spec.setMetaValue("ion mobility upper limit", win.im_hi);
    }
    OpenMS::DataArrays::FloatDataArray im_array;
    OpenMS::IMDataConverter::setIMUnit(im_array, OpenMS::DriftTimeUnit::VSSC);
    spec.reserve(n / wins.size() + 64);
    im_array.reserve(n / wins.size() + 64);
    for (OpenMS::Size k = 0; k < n; ++k)
    {
      const double k0 = im[k];
      if (k0 < win.im_lo || k0 > win.im_hi) { continue; }
      if (! (in[k] > 0)) { continue; }
      spec.emplace_back(mz[k], in[k]);
      im_array.push_back(static_cast<float>(k0));
    }
    if (spec.empty()) { continue; }
    // mobility-major on disk -> sort by m/z, carrying the IM array along
    std::vector<OpenMS::Size> ord(spec.size());
    for (OpenMS::Size k = 0; k < ord.size(); ++k)
    {
      ord[k] = k;
    }
    std::stable_sort(ord.begin(), ord.end(), [&](OpenMS::Size a, OpenMS::Size b) { return spec[a].getMZ() < spec[b].getMZ(); });
    std::vector<OpenMS::Peak1D> sorted;
    sorted.reserve(spec.size());
    OpenMS::DataArrays::FloatDataArray im_sorted;
    OpenMS::IMDataConverter::setIMUnit(im_sorted, OpenMS::DriftTimeUnit::VSSC);
    im_sorted.reserve(spec.size());
    for (OpenMS::Size k : ord)
    {
      sorted.push_back(spec[k]);
      im_sorted.push_back(im_array[k]);
    }
    spec.clear(false); // peaks only; RT/precursor/meta stay
    for (const OpenMS::Peak1D& pk : sorted)
    {
      spec.push_back(pk);
    }
    spec.getFloatDataArrays().push_back(std::move(im_sorted));
    spec.setIMPeakType(OpenMS::IMPeakType::IM_PROFILE);
  }
}

/**
  @brief [frames] The archive's run metadata from ONE sweep (no peak decode).

  [frames] The archive's run metadata from ONE sweep (no peak decode): per spectrum its RT,
  vendor frame id, isolation windows with their 1/K0 bands, and number_of_peaks (MS:1003059,
  -1 when the archive does not record it); the MS1/MS2 index lists; the exact calibration.
  Cached by path so the frozen frame tables and every tile of a tiled read share it: 72k
  metadata() calls once per run, never per tile.
*/
struct MzPeakRunMeta
{
  std::string path;                         ///< archive path (the cache key)
  std::vector<double> rt;                   ///< retention time per spectrum
  std::vector<long> frame_id;               ///< vendor frame id per spectrum (-1 when the native ID has none)
  std::vector<std::vector<MzPeakWin>> wins; ///< isolation windows per spectrum
  std::vector<OpenMS::Size> ms1, ms2;       ///< indices of the MS1 spectra and of the MS2 spectra with windows
  std::vector<long> n_peaks;                ///< number_of_peaks per spectrum (-1 when not recorded)
  OpenMS::Size n_peaks_missing_ms2 = 0;     ///< MS2 spectra without number_of_peaks (they count as present)
  std::shared_ptr<MzPeakExactMz> exact;     ///< exact calibration, or null when DIASPEXTRACT_MZPEAK_EXACT=0
  std::string calibration;                  ///< lastMzPeakCalibration() as the sweep left it (restored on a cache hit)
};
/**
  @brief The cached run metadata of the last swept archive.

  @return A reference to the process-wide cache slot (null before the first sweep)
*/
inline std::shared_ptr<const MzPeakRunMeta>& lastMzPeakMeta()
{
  static std::shared_ptr<const MzPeakRunMeta> p;
  return p;
}

/**
  @brief Sweep an archive's spectrum metadata once, or return the cached sweep for the same path.

  Also sets lastMzPeakCalibration() and, unless DIASPEXTRACT_MZPEAK_EXACT=0, loads the exact TDF calibration.
  Each MS2 window is matched to the selected ion nearest its isolation target to obtain its 1/K0 band.

  @param[in] path Path of the .mzpeak archive
  @return The run metadata (shared with the cache)
  @throws OpenMS::Exception::InvalidValue if the exact calibration cannot be recovered or an MS2 window has no matching ion-mobility band;
          a standard exception from the mzPeak, Arrow or zip libraries passes through (mzPeakRunMeta() translates it)
*/
inline std::shared_ptr<const MzPeakRunMeta> mzPeakRunMetaUnchecked(const std::string& path)
{
  if (lastMzPeakMeta() && lastMzPeakMeta()->path == path)
  {
    // a failed sweep of another archive in between may have reset it
    lastMzPeakCalibration() = lastMzPeakMeta()->calibration;
    return lastMzPeakMeta();
  }
  auto meta = std::make_shared<MzPeakRunMeta>();
  meta->path = path;
  MzPeak::Index index = MzPeak::open(path);
  MzPeak::Spectra spectra = index.spectra();
  const OpenMS::Size n = spectra.size();
  lastMzPeakCalibration() = "mzpeak_two_point_transform (library-applied MS:1003825 from the archive; not the TDF MzCalibration model) archive="
                            + path.substr(path.find_last_of('/') + 1);
  std::vector<double>& rt = meta->rt;
  rt.assign(n, 0.0);
  std::vector<long>& frame_id = meta->frame_id;
  frame_id.assign(n, -1);
  std::vector<std::vector<MzPeakWin>>& wins = meta->wins;
  wins.assign(n, {});
  std::vector<OpenMS::Size>& ms1 = meta->ms1;
  std::vector<OpenMS::Size>& ms2 = meta->ms2;
  meta->n_peaks.assign(n, -1);
  std::shared_ptr<MzPeakExactMz>& exact = meta->exact;
  // Exact TDF model by default, fail closed like the .d path; the refusal below names the cost and both escapes.
  const char* ex = std::getenv("DIASPEXTRACT_MZPEAK_EXACT");
  if (! (ex && std::string(ex) == "0"))
  {
    try
    {
      exact = std::make_shared<MzPeakExactMz>(path);
    }
    catch (const std::exception& e)
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            std::string("mzPeak input: cannot recover the exact TDF calibration (") + e.what()
                                              + "). The archive's two-point transform is less exact and identifies fewer peptides; set "
                                                "DIASPEXTRACT_MZPEAK_TDF=<analysis.tdf.gz> or DIASPEXTRACT_MZPEAK_EXACT=0 to accept it.",
                                            path);
    }
    lastMzPeakCalibration() = "tdf_table_modeltype1 (recovered from the archive's two-point transform + embedded vendor/analysis.tdf.gz) archive="
                              + path.substr(path.find_last_of('/') + 1);
  }
  for (OpenMS::Size i = 0; i < n; ++i)
  {
    MzPeak::Spectrum sp = spectra[i];
    const auto& md = sp.metadata();
    rt[i] = md.retention_time.value_or(0.0);
    {
      const std::string& id = md.id;
      const auto pos = id.find("frame=");
      if (pos != std::string::npos) { frame_id[i] = std::atol(id.c_str() + pos + 6); }
    }
    if (md.number_of_peaks) { meta->n_peaks[i] = (long)*md.number_of_peaks; }
    if (md.ms_level.value_or(1) == 1)
    {
      wins[i] = {MzPeakWin {0, 0, -1e9, 1e9}};
      ms1.push_back(i);
      continue;
    }
    if (! md.number_of_peaks) { ++meta->n_peaks_missing_ms2; }
    // Selected ions carry the per-window 1/K0 band. mzpeak-convert 0.9.x writes precursor_index
    // as NULL and the band as NAME-ONLY CV params; the reader then attaches every ion of the
    // frame to the FIRST precursor. So: flatten all ions of the spectrum and match each window
    // to the ion whose selected_ion_mz is its isolation target.
    struct Ion
    {
      double mz, lo, hi;
      bool used = false;
    };
    std::vector<Ion> ions;
    for (const auto& p : md.precursors)
    {
      for (const auto& si : p.selected_ions)
      {
        Ion ion {si.selected_ion_mz.value_or(-1.0), -1e9, 1e9};
        if (si.ion_mobility_lower_limit && si.ion_mobility_upper_limit)
        {
          ion.lo = *si.ion_mobility_lower_limit;
          ion.hi = *si.ion_mobility_upper_limit;
        }
        for (const auto& cv : si.parameters)
        {
          if (! cv.name || ! cv.value) { continue; }
          // clang-format off: cpplint wants braces on both sides of else; the OpenMS style would put these blocks on one line
          if (cv.name->find("ion mobility lower limit") != std::string::npos)
          {
            ion.lo = std::stod(*cv.value);
          }
          else if (cv.name->find("ion mobility upper limit") != std::string::npos)
          {
            ion.hi = std::stod(*cv.value);
          }
          // clang-format on
        }
        ions.push_back(ion);
      }
    }
    for (const auto& p : md.precursors)
    {
      if (! p.isolation_window.target_mz) { continue; }
      const double c = *p.isolation_window.target_mz;
      const double lo = c - p.isolation_window.lower_offset.value_or(0.f), hi = c + p.isolation_window.upper_offset.value_or(0.f);
      // the nearest ion not taken by an earlier window: two windows with one centre (and different widths or
      // mobility bands) get their own ions, in the order the archive lists them
      Ion* best = nullptr;
      for (Ion& ion : ions)
      {
        if (! ion.used && (! best || std::fabs(ion.mz - c) < std::fabs(best->mz - c))) { best = &ion; }
      }
      if (! best || std::fabs(best->mz - c) > 0.05 || best->lo < -1e8 || best->hi > 1e8)
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "mzPeak MS2 window without a matching ion-mobility band (spectrum " + std::to_string(i) + ")",
                                              std::to_string(c));
      }
      best->used = true;
      wins[i].push_back(MzPeakWin {lo, hi, best->lo, best->hi});
    }
    if (! wins[i].empty()) { ms2.push_back(i); }
  }
  meta->calibration = lastMzPeakCalibration();
  lastMzPeakMeta() = meta;
  return meta;
}

/**
  @brief [noims] Whether the archive's spectra carry a per-peak ion-mobility array.

  The streaming reader is for diaPASEF
  archives; one without ion mobility (Orbitrap, FAIMS) is read whole and takes the non-IMS path. Without this its
  stream preflight refused it (no exact TDF calibration) under the default integer detector, before any fallback.
  The first non-empty one of the first 16 spectra decides. An archive this cannot read counts as one with ion
  mobility, so its stream attempt and errors stay what they were.

  @param[in] path Path of the .mzpeak archive
  @return false only if the first non-empty spectrum has no ion-mobility array
*/
inline bool mzPeakHasIonMobility(const std::string& path)
{
  try
  {
    MzPeak::Index index = MzPeak::open(path);
    MzPeak::Spectra spectra = index.spectra();
    std::vector<OpenMS::Size> ids(std::min<OpenMS::Size>(spectra.size(), 16));
    for (OpenMS::Size i = 0; i < ids.size(); ++i)
    {
      ids[i] = i;
    }
    for (const MzPeak::Spectrum& sp : spectra.get_spectra_batch(ids))
    {
      if (! sp.mz().empty()) { return ! sp.ion_mobility_array().empty(); }
    }
  }
  catch (...)
  {
  }
  return true;
}

/**
  @brief Stream an .mzpeak run into `consumer` (MS1 first, then MS2 per (frame, window), frame order).

  Frames are decoded in parallel in batches of contiguous ranges (one reader per thread); non-empty spectra are
  handed to the consumer serially in frame order.

  @param[in] path Path of the .mzpeak archive
  @param[in,out] consumer Receives every non-empty spectrum
  @param[in] threads Requested number of decode threads (capped by omp_get_max_threads())
  @throws OpenMS::Exception::InvalidValue if the metadata sweep or the decode of a frame fails (rethrown from the worker);
          a standard exception from the mzPeak, Arrow or zip libraries passes through (loadMzPeakStreaming() translates it)
*/
inline void loadMzPeakStreamingUnchecked(const std::string& path, OpenMS::FullSwathFileConsumer& consumer, int threads)
{
  // 1) the run metadata (cached sweep)
  const std::shared_ptr<const MzPeakRunMeta> meta = mzPeakRunMetaUnchecked(path);
  const MzPeakRunMeta& m = *meta;

  // 2) decode in parallel over contiguous ranges (row-group locality), hand off serially in order
  const int nthr = std::min(std::max(1, threads), std::max(1, omp_get_max_threads()));
  auto run = [&](const std::vector<OpenMS::Size>& ids) {
    const OpenMS::Size per = 12; // frames per thread per batch (~1 row group)
    const OpenMS::Size batch = std::max<OpenMS::Size>(per, per * (OpenMS::Size)nthr);
    for (OpenMS::Size b0 = 0; b0 < ids.size(); b0 += batch)
    {
      const OpenMS::Size b1 = std::min(ids.size(), b0 + batch);
      const int nchunk = (int)((b1 - b0 + per - 1) / per);
      std::vector<std::vector<std::vector<OpenMS::MSSpectrum>>> out((OpenMS::Size)nchunk);
      std::exception_ptr err;
  #pragma omp parallel for schedule(dynamic, 1) num_threads(nthr)
      for (int c = 0; c < nchunk; ++c)
      {
        try
        {
          const OpenMS::Size c0 = b0 + (OpenMS::Size)c * per, c1 = std::min(b1, c0 + per);
          // one reader per THREAD: MzPeak::open re-reads the zip directory and every parquet footer
          static thread_local std::unique_ptr<MzPeak::Index> tidx; // Index is neither movable nor copyable:
          static thread_local std::string tidx_path;               // construct it from the prvalue (elided)
          if (! tidx || tidx_path != path)
          {
            tidx.reset(new MzPeak::Index(MzPeak::open(path)));
            tidx_path = path;
          }
          MzPeak::Spectra tsp = tidx->spectra();
          std::vector<OpenMS::Size> want(ids.begin() + (std::ptrdiff_t)c0, ids.begin() + (std::ptrdiff_t)c1);
          std::vector<MzPeak::Spectrum> got = tsp.get_spectra_batch(want);
          // A short batch is an error, never silently fewer frames.
          if (got.size() != want.size())
          {
            throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "mzPeak reader returned fewer frames than requested",
                                                  std::to_string(got.size()) + " of " + std::to_string(want.size()));
          }
          auto& slot = out[(OpenMS::Size)c];
          slot.resize(want.size());
          for (OpenMS::Size k = 0; k < want.size(); ++k)
            frameToSpectra(got[k], m.rt[want[k]], want[k], m.wins[want[k]], slot[k], m.exact.get(), m.frame_id[want[k]]);
        }
        catch (...)
        {
  #pragma omp critical(mzpeak_err)
          if (! err) err = std::current_exception();
        }
      }
      if (err) std::rethrow_exception(err);
      for (auto& chunk : out)
      {
        for (auto& frame : chunk)
        {
          for (auto& spec : frame)
            if (! spec.empty()) consumer.consumeSpectrum(spec);
        }
      }
    }
  };
  run(m.ms1);
  run(m.ms2);
}

/**
  @brief Run @p f; a standard exception from the mzPeak, Arrow or zip libraries becomes an MzPeakError.

  OpenMS exceptions pass through unchanged, so their type and source location stay as thrown.

  @param[in] f The call to run
  @return What @p f returns
  @throws MzPeakError for a standard exception thrown by @p f
*/
template<class F>
decltype(auto) withMzPeakErrors(F&& f)
{
  try
  {
    return f();
  }
  catch (const OpenMS::Exception::BaseException&)
  {
    throw;
  }
  catch (const std::exception& e)
  {
    throw MzPeakError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, std::string("mzPeak reader: ") + e.what(), "");
  }
}

/**
  @brief The run metadata of an .mzpeak archive (see mzPeakRunMetaUnchecked()), with OpenMS exceptions only.

  @param[in] path Path of the .mzpeak archive
  @return The cached or newly swept metadata
  @throws OpenMS::Exception::InvalidValue (or MzPeakError, a subclass) if the archive cannot be read as an ims-compact run
*/
inline std::shared_ptr<const MzPeakRunMeta> mzPeakRunMeta(const std::string& path)
{
  return withMzPeakErrors([&] { return mzPeakRunMetaUnchecked(path); });
}

/**
  @brief Stream an .mzpeak run into @p consumer (see loadMzPeakStreamingUnchecked()), with OpenMS exceptions only.

  @param[in] path Path of the .mzpeak archive
  @param[in,out] consumer Receives every non-empty spectrum
  @param[in] threads Requested number of decode threads
  @throws OpenMS::Exception::InvalidValue (or MzPeakError, a subclass) if the metadata sweep or a frame decode fails
*/
inline void loadMzPeakStreaming(const std::string& path, OpenMS::FullSwathFileConsumer& consumer, int threads)
{
  withMzPeakErrors([&] { loadMzPeakStreamingUnchecked(path, consumer, threads); });
}
} // namespace DIAspeXtract
#endif // DIASPEXTRACT_WITH_MZPEAK
