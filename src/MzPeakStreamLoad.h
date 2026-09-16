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
#include <unistd.h>

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

namespace spx
{
  struct MzPeakWin { double lo, hi, im_lo, im_hi; };

  /// Calibration provenance of the last mzPeak load (the exact TDF model, or the archive's two-point transform), stamped as spx:mz_calibration.
  inline std::string& lastMzPeakCalibration() { static std::string s = "unset"; return s; }


  /// Exact TOF->m/z for an ims-compact archive: (c0, c1) from the peaks table's column metadata, the
  /// MzCalibration row and Frames.T1 from the embedded vendor tdf. Throws on anything missing (fail closed).
  struct MzPeakExactMz
  {
    double c0 = 0, c1 = 0;
    diaspextract::TdfMzCalibration cal;
    std::vector<double> t1_by_frame;          // index = Frames.Id
    double acq_lo = std::numeric_limits<double>::quiet_NaN(), acq_hi = std::numeric_limits<double>::quiet_NaN();   // GlobalMetadata MzAcqRange
    long n_bins = 0; std::string bounds_why;                                                                        // DigitizerNumSamples, or why not

    static std::vector<char> readMember_(zip_t* z, const char* name, zip_int64_t from = 0, zip_int64_t len = -1)
    {
      zip_int64_t idx = zip_name_locate(z, name, 0);
      if (idx < 0) throw std::runtime_error(std::string("mzpeak: member missing: ") + name);
      zip_stat_t st; zip_stat_index(z, (zip_uint64_t)idx, 0, &st);
      if (len < 0) len = (zip_int64_t)st.size - from;
      zip_file_t* f = zip_fopen_index(z, (zip_uint64_t)idx, 0);
      if (!f) throw std::runtime_error(std::string("mzpeak: cannot open member ") + name);
      if (from > 0 && zip_fseek(f, from, SEEK_SET) != 0) { zip_fclose(f); throw std::runtime_error(std::string("mzpeak: cannot seek in ") + name + " (member not stored?)"); }
      std::vector<char> buf((size_t)len); zip_int64_t got = 0;
      while (got < len) { zip_int64_t n = zip_fread(f, buf.data() + got, (zip_uint64_t)(len - got)); if (n <= 0) break; got += n; }
      zip_fclose(f);
      if (got != len) throw std::runtime_error(std::string("mzpeak: short read of ") + name);
      return buf;
    }

    explicit MzPeakExactMz(const std::string& path)
    {
      int err = 0; zip_t* z = zip_open(path.c_str(), ZIP_RDONLY, &err);
      if (!z) throw std::runtime_error("mzpeak: zip_open failed on " + path);
      try
      {
        // (1) two-point transform params from the parquet footer of spectra_peaks.parquet (stored member)
        zip_int64_t idx = zip_name_locate(z, "spectra_peaks.parquet", 0);
        if (idx < 0) throw std::runtime_error("mzpeak: spectra_peaks.parquet missing");
        zip_stat_t st; zip_stat_index(z, (zip_uint64_t)idx, 0, &st);
        std::vector<char> tail = readMember_(z, "spectra_peaks.parquet", (zip_int64_t)st.size - 8, 8);
        uint32_t flen; std::memcpy(&flen, tail.data(), 4);
        if (std::string(tail.data() + 4, 4) != "PAR1") throw std::runtime_error("mzpeak: peaks member is not parquet");
        std::vector<char> footer = readMember_(z, "spectra_peaks.parquet", (zip_int64_t)st.size - 8 - flen, flen);
        uint32_t flen2 = flen;
        std::shared_ptr<parquet::FileMetaData> md = parquet::FileMetaData::Make(footer.data(), &flen2);
        std::shared_ptr<arrow::Schema> schema;
        parquet::ArrowReaderProperties props;
        PARQUET_THROW_NOT_OK(parquet::arrow::FromParquetSchema(md->schema(), props, md->key_value_metadata(), &schema));
        auto point = schema->GetFieldByName("point");
        if (!point) throw std::runtime_error("mzpeak: no 'point' column");
        auto tof = std::dynamic_pointer_cast<arrow::StructType>(point->type())->GetFieldByName("tof");
        if (!tof || !tof->metadata()) throw std::runtime_error("mzpeak: no tof column / metadata (not an ims-compact archive?)");
        auto tp = tof->metadata()->Get("mzpeak:transform_params");
        if (!tp.ok()) throw std::runtime_error("mzpeak: tof column has no mzpeak:transform_params");
        if (std::sscanf(tp->c_str(), "%lf,%lf", &c0, &c1) != 2 || !(c1 > 0)) throw std::runtime_error("mzpeak: bad transform_params " + *tp);

        // (2) the vendor tdf: gunzip vendor/analysis.tdf.gz to a temp file, read MzCalibration + Frames.T1
        // --no-vendor archives have no embedded tdf: accept a sidecar via DIASPEXTRACT_MZPEAK_TDF=<analysis.tdf.gz>
        std::vector<char> gz;
        if (const char* side = std::getenv("DIASPEXTRACT_MZPEAK_TDF"))
        {
          FILE* sf = std::fopen(side, "rb"); if (!sf) throw std::runtime_error(std::string("mzpeak: cannot read DIASPEXTRACT_MZPEAK_TDF ") + side);
          char b[1 << 16]; size_t n; while ((n = std::fread(b, 1, sizeof b, sf)) > 0) gz.insert(gz.end(), b, b + n); std::fclose(sf);
        }
        else gz = readMember_(z, "vendor/analysis.tdf.gz");
        char tmpl[] = "/tmp/spx_mzpeak_tdf_XXXXXX"; int fd = mkstemp(tmpl); if (fd < 0) throw std::runtime_error("mzpeak: mkstemp");
        close(fd);
        { z_stream zs{}; inflateInit2(&zs, 16 + MAX_WBITS); FILE* out = std::fopen(tmpl, "wb");
          zs.next_in = (Bytef*)gz.data(); zs.avail_in = (uInt)gz.size(); std::vector<char> ob(1 << 20); int rc;
          do { zs.next_out = (Bytef*)ob.data(); zs.avail_out = (uInt)ob.size(); rc = inflate(&zs, Z_NO_FLUSH);
               if (rc != Z_OK && rc != Z_STREAM_END) { std::fclose(out); inflateEnd(&zs); throw std::runtime_error("mzpeak: gunzip of analysis.tdf.gz failed"); }
               std::fwrite(ob.data(), 1, ob.size() - zs.avail_out, out); } while (rc != Z_STREAM_END);
          inflateEnd(&zs); std::fclose(out); }
        { std::string why;
          if (!diaspextract::loadTdfCalibration(std::string(tmpl), cal, t1_by_frame, why))
            throw std::runtime_error("mzpeak: " + why);
          std::string w2; diaspextract::TdfAxisBounds bd;
          if (diaspextract::loadTdfAxisBounds(std::string(tmpl), bd, w2)) { acq_lo = bd.mz_lo; acq_hi = bd.mz_hi; n_bins = bd.n_bins; } else bounds_why = w2; }
        std::remove(tmpl);
      }
      catch (...) { zip_close(z); throw; }
      zip_close(z);
    }

    /// Recover the integer tof from the archive's two-point m/z and re-apply the exact model.
    void exact(std::vector<double>& mz, long frame_id) const
    {
      // Fail rather than fall back to the reference T1: substituting it for an unknown frame stamps
      // masses calibrated at the wrong temperature as "exact", which is silent wrongness the
      // provenance record would then vouch for.
      if (frame_id < 0 || (size_t)frame_id >= t1_by_frame.size())
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "mzPeak exact calibration: no tdf T1 for frame", std::to_string(frame_id));
      const double b = cal.frameFactor(t1_by_frame[(size_t)frame_id]);
      for (double& v : mz)
      {
        if (!std::isfinite(v) || v <= 0.0)
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "mzPeak exact calibration: non-positive archive m/z", std::to_string(v));
        const double tof = std::llround((std::sqrt(v) - c0) / c1);
        const double e = cal.tofToMz(tof, b);
        if (!std::isfinite(e))
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "mzPeak exact calibration produced a non-finite m/z at tof", std::to_string(tof));
        v = e;
      }
    }
  };

  /// One decoded frame split per window (index-aligned with `wins` of that frame).
  inline void frameToSpectra_(const MzPeak::Spectrum& sp, double rt, std::size_t frame_idx,
                              const std::vector<MzPeakWin>& wins, std::vector<OpenMS::MSSpectrum>& out,
                              const MzPeakExactMz* exact, long frame_id)
  {
    using namespace OpenMS;
    std::vector<double> mz_exact;
    if (exact) { mz_exact = sp.mz(); exact->exact(mz_exact, frame_id); }
    const std::vector<double>& mz = exact ? mz_exact : sp.mz();
    const std::vector<float>& in = sp.intensity();
    const std::vector<double>& im = sp.ion_mobility_array();
    // Array lengths must agree. Truncating to the shorter one turns a malformed archive into a
    // smaller frame that looks valid, and the peaks that vanish are never reported anywhere.
    if (mz.size() != in.size())
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "mzPeak frame has mismatched m/z and intensity array lengths",
            std::to_string(mz.size()) + " vs " + std::to_string(in.size()));
    const std::size_t n = mz.size();
    if (im.size() < n)
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "mzPeak frame has a short ion-mobility array", std::to_string(im.size()));
    out.assign(wins.size(), MSSpectrum());
    for (std::size_t w = 0; w < wins.size(); ++w)
    {
      const MzPeakWin& win = wins[w];
      MSSpectrum& spec = out[w];
      spec.setRT(rt);
      spec.setMSLevel(wins.size() == 1 && !(win.lo < win.hi) ? 1 : 2);
      spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      // "frame=<Id>" keys the per-frame calibration (without it the integer detector falls back); "mzpeak=" is the archive index.
      spec.setNativeID("mzpeak=" + std::to_string(frame_idx) + (frame_id >= 0 ? " frame=" + std::to_string(frame_id) : std::string())
                       + " window=" + std::to_string(w));
      if (spec.getMSLevel() == 2)
      {
        Precursor prec;
        prec.setMZ(0.5 * (win.lo + win.hi));
        prec.setIsolationWindowLowerOffset(0.5 * (win.hi - win.lo));
        prec.setIsolationWindowUpperOffset(0.5 * (win.hi - win.lo));
        spec.setPrecursors({prec});
        spec.setMetaValue("ion mobility lower limit", win.im_lo);
        spec.setMetaValue("ion mobility upper limit", win.im_hi);
      }
      DataArrays::FloatDataArray im_array;
      IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
      spec.reserve(n / wins.size() + 64); im_array.reserve(n / wins.size() + 64);
      for (std::size_t k = 0; k < n; ++k)
      {
        const double k0 = im[k];
        if (k0 < win.im_lo || k0 > win.im_hi) continue;
        if (!(in[k] > 0)) continue;
        spec.emplace_back(mz[k], in[k]);
        im_array.push_back(static_cast<float>(k0));
      }
      if (spec.empty()) continue;
      // mobility-major on disk -> sort by m/z, carrying the IM array along
      std::vector<std::size_t> ord(spec.size());
      for (std::size_t k = 0; k < ord.size(); ++k) ord[k] = k;
      std::stable_sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return spec[a].getMZ() < spec[b].getMZ(); });
      std::vector<Peak1D> sorted; sorted.reserve(spec.size());
      DataArrays::FloatDataArray im_sorted; IMDataConverter::setIMUnit(im_sorted, DriftTimeUnit::VSSC); im_sorted.reserve(spec.size());
      for (std::size_t k : ord) { sorted.push_back(spec[k]); im_sorted.push_back(im_array[k]); }
      spec.clear(false);                                          // peaks only; RT/precursor/meta stay
      for (const Peak1D& pk : sorted) spec.push_back(pk);
      spec.getFloatDataArrays().push_back(std::move(im_sorted));
      spec.setIMPeakType(IMPeakType::IM_PROFILE);
    }
  }

  /// [frames] The archive's run metadata from ONE sweep (no peak decode): per spectrum its RT,
  /// vendor frame id, isolation windows with their 1/K0 bands, and number_of_peaks (MS:1003059,
  /// -1 when the archive does not record it); the MS1/MS2 index lists; the exact calibration.
  /// Cached by path so the frozen frame tables and every tile of a tiled read share it: 72k
  /// metadata() calls once per run, never per tile.
  struct MzPeakRunMeta
  {
    std::string path;
    std::vector<double> rt;
    std::vector<long> frame_id;
    std::vector<std::vector<MzPeakWin>> wins;
    std::vector<std::size_t> ms1, ms2;
    std::vector<long> n_peaks;
    std::size_t n_peaks_missing_ms2 = 0;   ///< MS2 spectra without number_of_peaks (they count as present)
    std::shared_ptr<MzPeakExactMz> exact;
    std::string calibration;   ///< lastMzPeakCalibration() as the sweep left it (restored on a cache hit)
  };
  inline std::shared_ptr<const MzPeakRunMeta>& lastMzPeakMeta() { static std::shared_ptr<const MzPeakRunMeta> p; return p; }

  inline std::shared_ptr<const MzPeakRunMeta> mzPeakRunMeta(const std::string& path)
  {
    if (lastMzPeakMeta() && lastMzPeakMeta()->path == path)
    {
      // a failed sweep of another archive in between may have reset it (frame-tables review)
      lastMzPeakCalibration() = lastMzPeakMeta()->calibration;
      return lastMzPeakMeta();
    }
    auto meta = std::make_shared<MzPeakRunMeta>();
    meta->path = path;
    MzPeak::Index index = MzPeak::open(path);
    MzPeak::Spectra spectra = index.spectra();
    const std::size_t n = spectra.size();
    lastMzPeakCalibration() = "mzpeak_two_point_transform (library-applied MS:1003825 from the archive; not the TDF MzCalibration model) archive=" + path.substr(path.find_last_of('/') + 1);
    std::vector<double>& rt = meta->rt; rt.assign(n, 0.0);
    std::vector<long>& frame_id = meta->frame_id; frame_id.assign(n, -1);
    std::vector<std::vector<MzPeakWin>>& wins = meta->wins; wins.assign(n, {});
    std::vector<std::size_t>& ms1 = meta->ms1; std::vector<std::size_t>& ms2 = meta->ms2;
    meta->n_peaks.assign(n, -1);
    std::shared_ptr<MzPeakExactMz>& exact = meta->exact;
    // Exact TDF model by default, fail closed like the .d path; the refusal below names the cost and both escapes.
    const char* ex = std::getenv("DIASPEXTRACT_MZPEAK_EXACT");
    if (!(ex && std::string(ex) == "0"))
    {
      try { exact = std::make_shared<MzPeakExactMz>(path); }
      catch (const std::exception& e)
      {
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          std::string("mzPeak input: cannot recover the exact TDF calibration (") + e.what() + "). The archive's two-point transform is ~7 ppm off and costs ~12% peptides; set DIASPEXTRACT_MZPEAK_TDF=<analysis.tdf.gz> or DIASPEXTRACT_MZPEAK_EXACT=0 to accept it.", path);
      }
      lastMzPeakCalibration() = "tdf_table_modeltype1 (recovered from the archive's two-point transform + embedded vendor/analysis.tdf.gz) archive=" + path.substr(path.find_last_of('/') + 1);
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      MzPeak::Spectrum sp = spectra[i];
      const auto& md = sp.metadata();
      rt[i] = md.retention_time.value_or(0.0);
      { const std::string& id = md.id; const auto pos = id.find("frame="); if (pos != std::string::npos) frame_id[i] = std::atol(id.c_str() + pos + 6); }
      if (md.number_of_peaks) meta->n_peaks[i] = (long)*md.number_of_peaks;
      if (md.ms_level.value_or(1) == 1) { wins[i] = {MzPeakWin{0, 0, -1e9, 1e9}}; ms1.push_back(i); continue; }
      if (!md.number_of_peaks) ++meta->n_peaks_missing_ms2;
      // Selected ions carry the per-window 1/K0 band. mzpeak-convert 0.9.x writes precursor_index
      // as NULL and the band as NAME-ONLY CV params; the reader then attaches every ion of the
      // frame to the FIRST precursor. So: flatten all ions of the spectrum and match each window
      // to the ion whose selected_ion_mz is its isolation target.
      struct Ion { double mz, lo, hi; };
      std::vector<Ion> ions;
      for (const auto& p : md.precursors)
        for (const auto& si : p.selected_ions)
        {
          Ion ion{si.selected_ion_mz.value_or(-1.0), -1e9, 1e9};
          if (si.ion_mobility_lower_limit && si.ion_mobility_upper_limit) { ion.lo = *si.ion_mobility_lower_limit; ion.hi = *si.ion_mobility_upper_limit; }
          for (const auto& cv : si.parameters)
          {
            if (!cv.name || !cv.value) continue;
            if (cv.name->find("ion mobility lower limit") != std::string::npos) ion.lo = std::stod(*cv.value);
            else if (cv.name->find("ion mobility upper limit") != std::string::npos) ion.hi = std::stod(*cv.value);
          }
          ions.push_back(ion);
        }
      for (const auto& p : md.precursors)
      {
        if (!p.isolation_window.target_mz) continue;
        const double c = *p.isolation_window.target_mz;
        const double lo = c - p.isolation_window.lower_offset.value_or(0.f), hi = c + p.isolation_window.upper_offset.value_or(0.f);
        const Ion* best = nullptr;
        for (const Ion& ion : ions) if (!best || std::fabs(ion.mz - c) < std::fabs(best->mz - c)) best = &ion;
        if (!best || std::fabs(best->mz - c) > 0.05 || best->lo < -1e8 || best->hi > 1e8)
          throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "mzPeak MS2 window without a matching ion-mobility band (spectrum " + std::to_string(i) + ")", std::to_string(c));
        wins[i].push_back(MzPeakWin{lo, hi, best->lo, best->hi});
      }
      if (!wins[i].empty()) ms2.push_back(i);
    }
    meta->calibration = lastMzPeakCalibration();
    lastMzPeakMeta() = meta;
    return meta;
  }

  /// Stream an .mzpeak run into `consumer` (MS1 first, then MS2 per (frame, window), frame order).
  inline void loadMzPeakStreaming(const std::string& path, OpenMS::FullSwathFileConsumer& consumer, int threads)
  {
    // 1) the run metadata (cached sweep)
    const std::shared_ptr<const MzPeakRunMeta> meta = mzPeakRunMeta(path);
    const MzPeakRunMeta& m = *meta;

    // 2) decode in parallel over contiguous ranges (row-group locality), hand off serially in order
    const int nthr = std::min(std::max(1, threads), std::max(1, omp_get_max_threads()));
    auto run = [&](const std::vector<std::size_t>& ids)
    {
      const std::size_t per = 12;                                  // frames per thread per batch (~1 row group)
      const std::size_t batch = std::max<std::size_t>(per, per * (std::size_t)nthr);
      for (std::size_t b0 = 0; b0 < ids.size(); b0 += batch)
      {
        const std::size_t b1 = std::min(ids.size(), b0 + batch);
        const int nchunk = (int)((b1 - b0 + per - 1) / per);
        std::vector<std::vector<std::vector<OpenMS::MSSpectrum>>> out((std::size_t)nchunk);
        std::exception_ptr err;
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nthr)
        for (int c = 0; c < nchunk; ++c)
        {
          try
          {
            const std::size_t c0 = b0 + (std::size_t)c * per, c1 = std::min(b1, c0 + per);
            // one reader per THREAD: MzPeak::open re-reads the zip directory and every parquet footer
            static thread_local std::unique_ptr<MzPeak::Index> tidx;   // Index is neither movable nor copyable:
            static thread_local std::string tidx_path;                  // construct it from the prvalue (elided)
            if (!tidx || tidx_path != path) { tidx.reset(new MzPeak::Index(MzPeak::open(path))); tidx_path = path; }
            MzPeak::Spectra tsp = tidx->spectra();
            std::vector<std::size_t> want(ids.begin() + (std::ptrdiff_t)c0, ids.begin() + (std::ptrdiff_t)c1);
            std::vector<MzPeak::Spectrum> got = tsp.get_spectra_batch(want);
            // A short batch is an error, never silently fewer frames.
            if (got.size() != want.size())
              throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "mzPeak reader returned fewer frames than requested",
                    std::to_string(got.size()) + " of " + std::to_string(want.size()));
            auto& slot = out[(std::size_t)c]; slot.resize(want.size());
            for (std::size_t k = 0; k < want.size(); ++k)
              frameToSpectra_(got[k], m.rt[want[k]], want[k], m.wins[want[k]], slot[k], m.exact.get(), m.frame_id[want[k]]);
          }
          catch (...)
          {
            #pragma omp critical(mzpeak_err)
            if (!err) err = std::current_exception();
          }
        }
        if (err) std::rethrow_exception(err);
        for (auto& chunk : out) for (auto& frame : chunk) for (auto& spec : frame)
          if (!spec.empty()) consumer.consumeSpectrum(spec);
      }
    };
    run(m.ms1);
    run(m.ms2);
  }
} // namespace spx
#endif // DIASPEXTRACT_WITH_MZPEAK
