// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: OpenMS-the reference implementation project $
// --------------------------------------------------------------------------

#include "MzPeakStreamLoad.h"   // [mzpeak] streaming .mzpeak input (DIASPEXTRACT_WITH_MZPEAK)
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <sys/resource.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include "TdfLoad.h"   // loadTdfCalibration: sqlite3, tool-only
#include "DnoiseMs1.h"   // [dnoise] the dnoise v0.1.0 MS1 port (header-only; tests/test_dnoise_ms1.cpp compiles the same code)
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/VersionInfo.h>   // verboseVersion_ carries the OpenMS build alongside ours

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cassert>
#include <type_traits>
#include <array>
#include <functional>
#include <cstdio>
#include <fstream>
#include <atomic>
#include <mutex>
#include <condition_variable>   // [outer] the MS1 pick worker
#include <thread>
#include <map>
#include <string>
#include <vector>
#include <cstdlib>     // getenv
#include <cstring>     // memcpy
#include <optional>
#include <memory>      // [outer] unique_ptr
#include <time.h>      // clock_gettime(CLOCK_THREAD_CPUTIME_ID): [dnoise] filter CPU
#ifdef __GLIBC__
#include <malloc.h>    // mallinfo2: retained-vs-live at the milestones [mem]
#endif

// The four OpenMS patches in patches/ are build prerequisites. This is the compile-time check for the
// third: without defaulted moves an rvalue MassTrace binds to the copy constructor (not noexcept), and
// every "move" of a trace in the band/split gathering is a deep copy of its points.
static_assert(std::is_nothrow_move_constructible_v<OpenMS::MassTrace> && std::is_nothrow_move_assignable_v<OpenMS::MassTrace>,
              "OpenMS MassTrace has no noexcept move operations: apply patches/openms-masstrace-move.patch to the OpenMS tree");
#include <exception>   // exception_ptr
#include <omp.h>
#include <iomanip>
#include <sstream>

// The release version. The standalone CMake build passes it from project(VERSION); the in-tree
// OpenMS build does not, so the fallback here is the same number.
#ifndef DIASPEXTRACT_VERSION
#define DIASPEXTRACT_VERSION "1.6.0"
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_DIAspeXtract DIAspeXtract

  @brief Extracts pseudo-MS/MS ("pseudo-DDA") spectra from diaPASEF (ion-mobility DIA) data.

  Reconstructs DDA-like MS2 spectra from ion-mobility DIA (timsTOF / diaPASEF) data by
  exploiting the fact that, in PASEF, fragments are recorded at the ion-mobility (1/K0)
  elution coordinate of their precursor. Precursor and fragment mass traces are detected
  IM-aware (@ref MassTraceDetection with an ion-mobility tolerance), then fragments are
  assigned to a precursor when they co-localize in ion mobility, retention time, and
  elution-profile (Pearson) correlation. The result is written as searchable mzML.

  This is an OpenMS-native, BSD-3 analogue of the published diaPASEF pseudo-MS/MS extractors. It is
  a discovery front-end (enables open / semi-tryptic / PTM searches); it does not perform the database
  search, FDR, library building or quantification.

  See docs/BASELINE.md for the measured decisions behind every default.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_DIAspeXtract.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_DIAspeXtract.html
*/

/// @cond TOPPCLASSES

namespace
{
  /// Seconds since tool start. [perf]
  static double phase_clock_()
  {
    static const auto t0 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
  static std::string clk_() { char b[32]; snprintf(b, sizeof b, " [t=%.0fs]", phase_clock_()); return b; }

  /// [perf-instr] Per-phase wall/CPU/RSS accounting, printed unconditionally at the end.
  struct PhaseStat { double wall = 0; double cpu = 0; long rss_end_mb = 0; int n = 0; bool nested = false; };
  /// [perf-instr] A phase opened inside another (the tile loop's read/sort/write inside WINDOW_LOOP)
  /// is marked enclosed and left out of the TOTAL. Phases are opened on the master thread, and with pipelined
  /// tiles also by the writer task while the master reads (atomic depth, locked phaseAdd_).
  inline std::atomic<int>& phase_depth_() { static std::atomic<int> d{0}; return d; }
  inline double& flush_wall_() { static double v = 0; return v; }   // [perf-load] pick inside LOAD
  /// [task-err] An exception escaping an OpenMP task terminates the process with no message. Task bodies catch
  /// everything and hand the first exception here; the code after the loop rethrows it once the tasks have joined.
  inline void taskErr_(std::exception_ptr& first)
  {
    #pragma omp critical(task_err)
    if (!first) first = std::current_exception();
  }
  /// [task-err] e2e falsifier only: DIASPEXTRACT_TASK_THROW=<site> throws inside that task, so the suite can show the run
  /// ends with an error message and a non-zero exit instead of an abort.
  inline void taskFault_(const char* site)
  {
    static const char* want = std::getenv("DIASPEXTRACT_TASK_THROW");
    if (want && String(want) == site)
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("[task-err] injected failure in task ") + site, site);
  }
  /// [pipe] Helper tasks of teamFor_ created and not yet started, over every call: capped at the team size, so
  /// stale helpers (started after their loop ran out) cannot pile up towards the libgomp cliff (64 x nthreads).
  inline std::atomic<int>& teamHelpersQueued_() { static std::atomic<int> q{0}; return q; }
  /// [pipe] A loop that the pipelined tile loop runs on the window team's driver (the tile reader and writer).
  /// Outside a parallel region the call sites keep their `omp parallel for`. Inside the team a `parallel for`
  /// would open a NESTED team -- fresh threads on the cores the windows use -- so the loop runs here instead:
  /// the caller pulls `chunk` indices at a time from a shared cursor, and up to one helper task per team
  /// thread pulls with it once a thread is free. The caller never reaches a task scheduling point, so it is
  /// never handed a window task (libgomp lets a thread waiting at a taskgroup end run ANY child of its task),
  /// and it never waits for a helper to START: a helper that starts after the cursor ran out touches nothing
  /// but its own reference to the cursor, so only helpers that claimed a chunk are waited for (seq_cst: a
  /// helper raises `busy` before it claims). f(i, slot): `slot` = the puller's omp_get_thread_num(), unique
  /// among the concurrent pullers (a thread pulls one chunk at a time and never nests a pull), for per-slot
  /// scratch where a thread number was used. Every caller writes index-addressed results only, so the output
  /// is the parallel for's. The first exception is rethrown here, after every claimed chunk has finished.
  /// Helper tasks queue behind the band and score tasks, and while windows hold the threads they barely run
  /// (a first pipelined build: a heavy tile of a 60-min run read in 31-43 s at 64 threads, 5-7 s with the barrier). So the
  /// loop is also published for teamAssist_: window code pulls a chunk at its natural break points, which
  /// gives the I/O a share of the window threads without a nested team (libgomp starts fresh threads for every
  /// nested region, measured in a second build).
  struct PullJob_
  {
    std::atomic<long> cur{0}; std::atomic<int> busy{0}; std::exception_ptr err;
    long n = 0, chunk = 1;
    void (*call)(void*, long, int) = nullptr; void* fp = nullptr;
    /// up to max_chunks chunks (all with max_chunks < 0), returns how many it claimed; `busy` is raised BEFORE
    /// the first claim (see above)
    long pull(const long max_chunks)
    {
      const int slot = omp_get_thread_num();
      busy.fetch_add(1);
      long a = 0, k = 0;
      for (; (max_chunks < 0 || k < max_chunks) && (a = cur.fetch_add(chunk)) < n; ++k)
      { try { for (long i = a, e = std::min(a + chunk, n); i < e; ++i) call(fp, i, slot); }
        catch (...) { taskErr_(err); cur.store(n); } }
      busy.fetch_sub(1);
      return k;
    }
  };
  inline std::mutex& pullMu_() { static std::mutex m; return m; }
  inline std::vector<std::shared_ptr<PullJob_>>& pullJobs_() { static std::vector<std::shared_ptr<PullJob_>> v; return v; }
  inline std::atomic<int>& pullLive_() { static std::atomic<int> n{0}; return n; }
  /// [pipe] Help a live teamFor_ loop (the pipelined tile reader or writer) with up to max_chunks chunks; true if
  /// it did any. One relaxed load when none is live. Called inside the window team only (slot = thread number).
  inline bool teamAssist_(const long max_chunks)
  {
    if (pullLive_().load(std::memory_order_relaxed) == 0) return false;
    std::shared_ptr<PullJob_> j;   // the copy keeps the job alive; its loop body is the caller's until it returns
    { std::lock_guard<std::mutex> g(pullMu_());
      for (const auto& p : pullJobs_()) if (p->cur.load(std::memory_order_relaxed) < p->n) { j = p; break; } }
    return j && j->pull(max_chunks) > 0;
  }
  template <class F> static void teamFor_(long n, long chunk, F&& f)
  {
    if (n <= 0) return;
    using FT = std::remove_reference_t<F>;
    auto st = std::make_shared<PullJob_>();
    st->n = n; st->chunk = chunk; st->fp = (void*)&f;
    st->call = [](void* p, long i, int slot) { (*static_cast<FT*>(p))(i, slot); };
    { std::lock_guard<std::mutex> g(pullMu_()); pullJobs_().push_back(st); }
    ++pullLive_();
    const int team = std::max(1, omp_get_num_threads());
    const int nt = (int)std::min<long>((n + chunk - 1) / chunk, (long)team);
    for (int t = 1; t < nt && teamHelpersQueued_().load() < team; ++t)
    {
      ++teamHelpersQueued_();
      #pragma omp task default(none) firstprivate(st)
      { --teamHelpersQueued_(); st->pull(-1); }
    }
    st->pull(-1);
    { std::lock_guard<std::mutex> g(pullMu_());
      auto& v = pullJobs_(); v.erase(std::remove(v.begin(), v.end(), st), v.end()); }
    --pullLive_();
    while (st->busy.load() != 0) std::this_thread::yield();
    if (st->err) std::rethrow_exception(st->err);
  }
  /// [review] True on the pipelined tile loop's driver while it drives. Every window task is a child of the driver's
  /// implicit task, so a taskloop's taskgroup end on the driver (libgomp then runs ANY queued child of the waiting task)
  /// could hand it a whole window, stalling the reads and admissions: the [det] toTof() in the tile read reached one.
  /// parTasks_ keeps every pass the driver reaches inline instead (nt = 1 is forTasks_'s serial branch: the same result).
  inline bool& onPipeDriver_() { static thread_local bool v = false; return v; }
  /// [serial] Task count for one of a window's data-parallel passes over n items: at most 16 and none smaller than
  /// min_per; 1 (run inline) outside a parallel region and on the pipe driver. A bounded COUNT, never a grainsize over
  /// data: libgomp runs a whole taskloop inline once the team's task count passes 64 x nthreads ([gomp-cliff]).
  inline long parTasks_(size_t n, size_t min_per)
  {
    const long cap = omp_in_parallel() && !onPipeDriver_() ? std::min(16, omp_get_num_threads()) : 1;
    return std::max(1L, std::min(cap, (long)(n / std::max<size_t>(min_per, 1))));
  }
  /// [serial] body(c) for c in [0, nt): nt tasks of the current team (the creating thread waits in the taskloop's
  /// taskgroup and runs its own children, so a busy pool degrades to the serial loop), or inline for nt <= 1.
  /// The first exception is rethrown once every task has joined.
  template <class F> void forTasks_(long nt, F&& body)
  {
    if (nt <= 1) { body(0L); return; }
    std::exception_ptr err;
    #pragma omp taskloop num_tasks(nt) default(shared)
    for (long c = 0; c < nt; ++c)
    { try { body(c); } catch (...) { taskErr_(err); } }
    if (err) std::rethrow_exception(err);
  }
  /// [serial] Merge-path co-rank: how many of the first r elements of the STABLE merge of a[0, m) and b[0, k) come
  /// from a (a first on ties, as std::merge takes them).
  template <class T, class Less> size_t coRank_(size_t r, const T* a, size_t m, const T* b, size_t k, const Less& less)
  {
    size_t lo = r > k ? r - k : 0, hi = std::min(r, m);
    while (lo < hi)
    {
      const size_t i = lo + (hi - lo) / 2, j = r - i;   // i < m and 1 <= j <= k here
      if (!less(b[j - 1], a[i])) lo = i + 1; else hi = i;
    }
    return lo;
  }
  /// [serial] Sort a[0, n) as nt tasks: nt chunks std::sort-ed, then pairwise std::merge levels, each pair's output
  /// cut at merge-path co-ranks into pieces of about n / nt, ping-ponging with one scratch buffer of n elements.
  /// The result equals std::sort's whenever `less` is a TOTAL order (no two elements equivalent: the sorted sequence
  /// is then unique). Under a weaker order equivalent elements may come out in another order.
  template <class T, class Less> void parallelSort_(T* a, size_t n, const Less& less, long nt)
  {
    static_assert(std::is_trivially_copyable<T>::value, "parallelSort_ moves elements through raw scratch memory");
    if (nt <= 1 || n < 2 * (size_t)nt) { std::sort(a, a + n, less); return; }
    vector<size_t> run(nt + 1);
    for (long c = 0; c <= nt; ++c) run[c] = (size_t)c * n / (size_t)nt;
    forTasks_(nt, [&](long c) { std::sort(a + run[c], a + run[c + 1], less); });
    // uninitialised scratch, first touched by the merge tasks (an implicit-lifetime type: malloc creates the objects)
    std::unique_ptr<T, void (*)(T*)> buf(static_cast<T*>(std::malloc(n * sizeof(T))), [](T* p) { std::free(p); });
    if (!buf) throw std::bad_alloc();
    T* src = a; T* dst = buf.get();
    struct Piece { size_t s, mid, e, r0, r1; };   // merge src[s, mid) with src[mid, e); output ranks [r0, r1)
    while (run.size() > 2)
    {
      const size_t nr = run.size() - 1;
      vector<size_t> nxt{0};
      vector<Piece> pcs;
      for (size_t p = 0; p < nr; p += 2)
      {
        const size_t s = run[p], mid = run[p + 1], e = p + 1 < nr ? run[p + 2] : mid;   // an odd last run is copied
        const size_t len = e - s, np = std::max<size_t>(1, (len * (size_t)nt + n - 1) / n);
        for (size_t q = 0; q < np; ++q) pcs.push_back({s, mid, e, q * len / np, (q + 1) * len / np});
        nxt.push_back(e);
      }
      forTasks_((long)pcs.size(), [&](long c) {
        const Piece& p = pcs[(size_t)c];
        const T* A = src + p.s; const T* B = src + p.mid;
        const size_t m = p.mid - p.s, k = p.e - p.mid;
        const size_t i0 = coRank_(p.r0, A, m, B, k, less), i1 = coRank_(p.r1, A, m, B, k, less);
        std::merge(A + i0, A + i1, B + (p.r0 - i0), B + (p.r1 - i1), dst + p.s + p.r0, less);
      });
      std::swap(src, dst); run.swap(nxt);
    }
    if (src != a)
      forTasks_(nt, [&](long c) { const size_t lo = (size_t)c * n / (size_t)nt, hi = (size_t)(c + 1) * n / (size_t)nt;
                                  std::copy(src + lo, src + hi, a + lo); });
  }

  /// [outer] std::sort's result, computed in parallel, for a comparator that is a strict TOTAL order on the elements
  /// present (no two equivalent): the sorted sequence is then unique, so any correct sort returns exactly what std::sort
  /// returns. Chunks are sorted concurrently and merged pairwise; each merge is cut into equal output pieces by merge path
  /// (exactly std::merge's sequence). Serial std::sort when small, on one thread, or inside a parallel region. Callers
  /// guarantee the total order by construction (a unique last key) or check the result (sortExact_).
  template <typename T, typename Less>
  void parallelSortTotal_(std::vector<T>& v, Less less)
  {
    const size_t n = v.size();
    const int nt = omp_get_max_threads();
    if (n < ((size_t)1 << 15) || nt < 2 || omp_in_parallel()) { std::sort(v.begin(), v.end(), less); return; }
    size_t nc = 1;                                              // a power of two, about one chunk per thread
    while (nc < (size_t)nt && n / (2 * nc) >= 4096) nc *= 2;
    auto bound = [n, nc](size_t c) { return c * n / nc; };     // chunk c is [bound(c), bound(c + 1))
    std::vector<T> buf(n);
    T* src = v.data(); T* dst = buf.data();
    #pragma omp parallel for schedule(dynamic, 1)
    for (long c = 0; c < (long)nc; ++c) std::sort(src + bound((size_t)c), src + bound((size_t)c + 1), less);
    for (size_t w = 1; w < nc; w *= 2)                          // runs of w chunks, merged pairwise
    {
      const size_t nm = nc / (2 * w), np = std::max<size_t>(1, (2 * (size_t)nt + nm - 1) / nm);   // merges, pieces per merge
      #pragma omp parallel for schedule(dynamic, 1)
      for (long q = 0; q < (long)(nm * np); ++q)
      {
        const size_t m = (size_t)q / np, p = (size_t)q % np;
        const size_t a0 = bound(2 * m * w), b0 = bound((2 * m + 1) * w), b1 = bound((2 * m + 2) * w);
        const T* a = src + a0; const T* b = src + b0;
        const size_t na = b0 - a0, nb = b1 - b0, k0 = p * (na + nb) / np, k1 = (p + 1) * (na + nb) / np;
        // merge path: how many of std::merge(a, b)'s first k outputs come from a (a[i] precedes b[j] iff !less(b[j], a[i]))
        auto corank = [&](size_t k) {
          size_t lo = k > nb ? k - nb : 0, hi = std::min(k, na);
          while (lo < hi) { const size_t i = lo + (hi - lo) / 2; if (!less(b[k - i - 1], a[i])) lo = i + 1; else hi = i; }
          return lo;
        };
        const size_t i0 = corank(k0), i1 = corank(k1);
        std::merge(a + i0, a + i1, b + (k0 - i0), b + (k1 - i1), dst + a0 + k0, less);
      }
      std::swap(src, dst);
    }
    if (src != v.data()) v.swap(buf);
  }

  /// [outer] Sorts @p v exactly as std::sort(v.begin(), v.end(), less) does, in parallel whenever that is provably the
  /// same result: a permutation is sorted in parallel (index as the last key, so a total order) and, when no two
  /// neighbours are equivalent under @p less (every element then has one place, and std::sort's sequence is that one),
  /// gathered. With an equivalent pair, @p v -- still untouched -- goes through std::sort, so the order of equal keys is
  /// std::sort's own. Returns the number of equal-key neighbours of the result (the same count either way).
  template <typename T, typename Less>
  size_t sortExact_(std::vector<T>& v, Less less)
  {
    const size_t n = v.size();
    auto serial = [&]() { std::sort(v.begin(), v.end(), less); size_t t = 0;
                          for (size_t i = 1; i < n; ++i) if (!less(v[i - 1], v[i])) ++t;
                          return t; };
    if (n < ((size_t)1 << 15) || n > (size_t)std::numeric_limits<uint32_t>::max() || omp_get_max_threads() < 2 || omp_in_parallel())
      return serial();
    std::vector<uint32_t> perm(n);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n; ++i) perm[(size_t)i] = (uint32_t)i;
    parallelSortTotal_(perm, [&v, &less](uint32_t a, uint32_t b) {
      if (less(v[a], v[b])) return true;
      if (less(v[b], v[a])) return false;
      return a < b; });
    size_t ties = 0;
    #pragma omp parallel for schedule(static) reduction(+ : ties)
    for (long i = 1; i < (long)n; ++i) if (!less(v[perm[(size_t)i - 1]], v[perm[(size_t)i]])) ++ties;
    if (ties) return serial();
    std::vector<T> out(n);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n; ++i) out[(size_t)i] = v[perm[(size_t)i]];
    v.swap(out);
    return 0;
  }

#if defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE == 13
  /// [outer] libstdc++ 13's std::__introsort_loop with the right part of each cut handed to a task: the two parts are
  /// disjoint and each is processed by the same calls in the same order as in the serial loop, so the array ends
  /// exactly as std::sort leaves it, equal keys included. Small parts run inline (bounded task count).
  template <typename It, typename Size, typename Cmp>
  void introsortLoopPar_(It first, It last, Size depth, Cmp cmp)
  {
    while (last - first > 16)   // libstdc++'s _S_threshold
    {
      if (depth == 0) { std::__partial_sort(first, last, last, cmp); return; }
      --depth;
      It cut = std::__unguarded_partition_pivot(first, last, cmp);
      if (last - cut >= 16384)
      {
        #pragma omp task default(shared) firstprivate(cut, last, depth, cmp)
        introsortLoopPar_(cut, last, depth, cmp);
      }
      else introsortLoopPar_(cut, last, depth, cmp);
      last = cut;
    }
  }
#endif

  /// [outer] std::sort(v.begin(), v.end(), less), equal keys in std::sort's own order, in parallel where the standard
  /// library is libstdc++ 13 (the release compiler): its std::__sort is introsort loop + final insertion sort, replayed
  /// here with tasks (introsortLoopPar_) and the same final pass. Any other library, small inputs or a caller already
  /// inside a parallel region: std::sort itself.
  template <typename T, typename Less>
  void stdSortPar_(std::vector<T>& v, Less less)
  {
#if defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE == 13
    if (v.size() >= ((size_t)1 << 15) && omp_get_max_threads() > 1 && !omp_in_parallel())
    {
      // DIASPEXTRACT_OUTER_VERIFY: also sort a copy with std::sort and report whether the two sequences are identical
      static const bool verify = std::getenv("DIASPEXTRACT_OUTER_VERIFY") != nullptr;
      std::vector<T> ref;
      if (verify) { ref = v; std::sort(ref.begin(), ref.end(), less); }
      auto cmp = __gnu_cxx::__ops::__iter_comp_iter(less);
      const auto first = v.begin(), last = v.end();
      const auto depth = std::__lg(last - first) * 2;   // as std::__sort
      #pragma omp parallel
      #pragma omp single
      {
        #pragma omp taskgroup
        introsortLoopPar_(first, last, depth, cmp);
      }
      std::__final_insertion_sort(first, last, cmp);
      if (verify) std::fprintf(stderr, "[outer] stdSortPar_ verify: %zu elements, %s\n", v.size(),
                               ref == v ? "identical to std::sort" : "DIFFERENT from std::sort");
      return;
    }
#endif
    std::sort(v.begin(), v.end(), less);
  }
  /// [gate-dump] Diagnostics, read-only: with DIASPEXTRACT_GATE_DUMP=<out.tsv> and DIASPEXTRACT_GATE_TARGETS=<tsv: id mz
  /// rt_sec im z>, every precursor within 20 ppm (mono or +-1 isotope at the target's charge), 20 s and 0.05 1/K0 of a
  /// target writes one PRECURSOR row (fragments passing the gate) and one row per fragment trace within 0.05 1/K0 and
  /// 15 s, with the first gate it fails (plain Pearson gate rules only, RT window as applied: gate:short_traces,
  /// gate:colocation and gate:shrinkage are not reflected in the fate). The spectra are unchanged; off unless both variables are set.
  struct GateDump_
  {
    struct T { std::string id; double mz, rt, im; int z; };
    std::vector<T> t;   ///< sorted by mz
    std::ofstream out;
    std::mutex mu;
    bool on = false;
    GateDump_()
    {
      const char* f = std::getenv("DIASPEXTRACT_GATE_DUMP");
      const char* g = std::getenv("DIASPEXTRACT_GATE_TARGETS");
      if (!f || !g) return;
      std::ifstream in(g);
      T x;
      while (in >> x.id >> x.mz >> x.rt >> x.im >> x.z) t.push_back(x);
      std::sort(t.begin(), t.end(), [](const T& a, const T& b) { return a.mz < b.mz; });
      out.open(f);
      out << "target\tiso_k\tpc_mz\tpc_z\tpc_rt\tpc_im\tpc_iso\tpc_guessed\tfrag_mz\tdim\tdrt\tnp\toverlap\tcorr\tinten\tfate\n";
      on = !t.empty() && out.good();
    }
    static GateDump_& get() { static GateDump_ d; return d; }
  };
  inline std::map<std::string, PhaseStat>& phase_stats_()
  { static std::map<std::string, PhaseStat> m; return m; }
  inline std::vector<std::string>& phase_order_()
  { static std::vector<std::string> v; return v; }
  static double cpu_seconds_()
  {
    struct rusage ru;                       // RUSAGE_SELF sums ALL threads -> parallel efficiency
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
    return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6
         + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
  }
  static long rss_mb_()
  {
    std::ifstream f("/proc/self/status"); std::string k;   // deliberate: Linux-only, returns 0 elsewhere
    while (f >> k) { if (k == "VmRSS:") { long v; f >> v; return v / 1024; } f.ignore(1 << 20, '\n'); }
    return 0;
  }
  static void phaseAdd_(const std::string& name, double w0, double c0, bool nested)
  {
    // [pipe] the pipelined tile loop times the tile read (driver) and the tile write (writer task) concurrently
    static std::mutex m; std::lock_guard<std::mutex> g(m);
    auto& st = phase_stats_()[name];
    if (st.n++ == 0) phase_order_().push_back(name);
    st.wall += phase_clock_() - w0; st.cpu += cpu_seconds_() - c0;
    if (nested) st.nested = true;
    st.rss_end_mb = rss_mb_();
  }
  struct Phase
  {
    std::string name; double w0, c0; int d;
    explicit Phase(std::string n) : name(std::move(n)), w0(phase_clock_()), c0(cpu_seconds_()), d(phase_depth_()++) {}
    ~Phase() { --phase_depth_(); phaseAdd_(name, w0, c0, d > 0); }
  };
  static void report_phases_(double total_wall)
  {
    if (phase_order_().empty()) return;
    OPENMS_LOG_INFO << "\n[perf] phase breakdown (wall s | % of total | CPU s | par x | RSS MB at end)\n";
    double acc = 0;
    for (const auto& n : phase_order_())
    {
      const auto& p = phase_stats_()[n];
      if (!p.nested) acc += p.wall;   // enclosed rows are already inside their parent's wall
      char b[220];
      const char* tag = p.nested ? (p.n > 1 ? "  (summed, enclosed)" : "  (enclosed)")
                                 : (p.n > 1 ? "  (summed)" : "");
      snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  %9.1f  %5.1fx  %8ld%s\n", n.c_str(), p.wall,
               total_wall > 0 ? 100.0 * p.wall / total_wall : 0.0, p.cpu,
               p.wall > 0 ? p.cpu / p.wall : 0.0, p.rss_end_mb, tag);
      OPENMS_LOG_INFO << b;
    }
    char b[160];
    snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  (unattributed %.1f s)\n", "TOTAL(measured)", acc,
             total_wall > 0 ? 100.0 * acc / total_wall : 0.0, total_wall - acc);
    OPENMS_LOG_INFO << b;
  }

  bool log_overlap_ = false;   ///< [B0] gate:coelution == "logoverlap"
  double im_weight_sigma_ = 0.0;     ///< [recall] Gaussian IM-proximity weight on emitted fragment intensity (0=off)
  double mono_guard_ = 0.0;          ///< [mono-guard] averagine slack on leftward isotope steps (0=off)
  bool mono_select_ = false;         ///< [mono-select] pick the mono by averagine fit over the isotope run
  int mono_check_ = 0;               ///< [monocheck] charge:mono_profile_check: 0 false, 1 true (v1 profile + 1/K0 steps), 2 v2
  double mono_check_min_corr_ = 0.6, mono_check_max_dim_ = 0.015;   ///< [monocheck] charge:mono_profile_min_corr, _max_dim
  double mono_pos_min_ratio_ = 0.2, mono_pos_im_box_ = 0.0, mono_pos_extend_ = 0.6;   ///< [monocheck-v2] charge:mono_position_* (im_box 0 = auto)
  bool mono_pos_iterate_ = true;     ///< [monocheck-v2] charge:mono_position_iterate
  String mono_check_stamp_;          ///< [monocheck] the spx:monocheck run stamp, filled by inferPrecursors_
  bool iso_closest_ = false;         ///< [isotol] charge:iso_partner == "closest": the nearest candidate in RT and 1/K0, not the first in m/z
  bool var_support_ = false;         ///< [Q1] Pearson variance over union support, not full grid G (0-pad fix)
  bool short_conc_ = false;          ///< [short] gate:short_traces == "concentric": span test for traces too short for Pearson
  double short_center_tol_ = 0.5;    ///< [short] gate:short_center_tol, fraction of the precursor half-span
  int short_min_pts_ = 1;            ///< [short] gate:short_min_points
  bool coloc_ = false;               ///< [coloc] gate:colocation == "ranges": each apex inside the other span, span overlap
  double coloc_overlap_ = 0.3;       ///< [coloc] gate:colocation_overlap, fraction of the shorter span
  bool shrink_ = false;              ///< [shrink] gate:shrinkage == "fisher": Fisher-z posterior mean instead of a hard point cut
  double shrink_r0_ = 0.0, shrink_sd_ = 0.5;   ///< [shrink] prior N(atanh r0, sd^2) on the Fisher z of the true correlation
  int shrink_min_pts_ = 2;           ///< [shrink] gate:shrink_min_points, overlapping frames needed to compute r at all
  bool adaptive_ = false;            ///< [adaptive] gate:adaptive == "faint": the relaxed modes only for faint precursors
  bool adaptive_auto_ = false;       ///< [adaptive] gate:adaptive == "auto": faint mode only when the run's faint fraction >= auto_faint_frac_
  double auto_faint_frac_ = 0.45;    ///< [adaptive] gate:auto_faint_fraction
  int faint_frames_ = 4;             ///< [adaptive] gate:faint_max_frames
  double faint_drt_ = 6.0, faint_cp_ = 1.0;   ///< [adaptive] gate:faint_delta_rt, gate:faint_corr_power
  double corr_power_ = 0.0;          ///< [Q2] emitted fragment intensity *= corr^corr_power (0=off; engine-agnostic)
  int cap_rank_ = 0;                 ///< [cap] assembly:cap_rank: 0 corr_inten (score x intensity), 1 umpire (r^2 ln I), 2 evidence (r ln I (1+overlap))
  double im_slice_merge_ = 0.015;   ///< [merge] trace:im_slice_merge, 1/K0 link (group span <= 2x); negative = off
  bool merge_same_bin_ = true;       ///< [merge] assembly:merge_same_bin
  int frag_int_ = 0;                 ///< [fragint] assembly:frag_intensity: 0 apex (legacy), 1 core, 2 sum
  int frag_core_half_ = 1;           ///< [fragint] (assembly:frag_core_frames - 1) / 2
  bool anchor_rescue_ = false;       ///< [anchor] gate:anchor_rescue: weak fragments of non-faint precursors judged against an anchor consensus
  int anchor_min_count_ = 3;         ///< [anchor] gate:anchor_min_count
  double anchor_min_corr_ = 0.8;     ///< [anchor] gate:anchor_min_corr
  double anchor_rescue_corr_ = 0.7;  ///< [anchor] gate:anchor_rescue_corr
  int anchor_max_rescued_ = 50;      ///< [anchor] gate:anchor_max_rescued
  bool post_fusion_cap_ = false;     ///< [recap] assembly:post_fusion_cap: carry the cap key per peak, re-cap after the twin fold
  bool twin_ext_mark_ = false;       ///< [twin-rt] assembly:twin_rt_match ext: mark spectra whose mono v2's extension placed (spx_mono_ext, never written)
  double foreign_w_ = 1.0;           ///< [foreign] assembly:foreign_precursor_weight (1 = off, 0 removes)
  double foreign_ppm_ = 10.0;        ///< [foreign] assembly:foreign_precursor_ppm
  double foreign_dim_ = 0.015;       ///< [foreign] assembly:foreign_precursor_dim
  int foreign_min_iso_ = 2;          ///< [foreign] assembly:foreign_precursor_min_iso
  double fuse_ppm_ = 0.0;            ///< [fuse] assembly:fragment_fuse_ppm (0 = off)
  bool fuse_linkage_ = false;        ///< [fuse] assembly:fragment_fuse_mode: false anchored, true linkage

  /// glibc arena accounting at a milestone: what the allocator holds in arenas vs mmapped blocks
  /// and how much of it is free-but-retained. RSS alone cannot tell live data from retained pages.
  static std::string mem_()
  {
#ifdef __GLIBC__
    const struct mallinfo2 mi = mallinfo2();
    return " [malloc arena=" + std::to_string(mi.arena >> 20) + " MB mmap=" + std::to_string(mi.hblkhd >> 20)
         + " MB free=" + std::to_string(mi.fordblks >> 20) + " MB]";
#else
    return "";
#endif
  }

  // [ledger] DIASPEXTRACT_LEDGER=<file>: a sampler thread appends one TSV row every 250 ms -- wall
  // clock, VmRSS, glibc arena/free, live bytes per category, the unattributed remainder (RSS minus
  // the sum: mostly allocator retention) and the phase marker. Charges are atomic adds at allocation
  // boundaries, never per peak; with the variable unset every call is one branch.
  enum LedCat { LC_SLAB, LC_MS1MAP, LC_MS1TRACE, LC_PRECURSOR, LC_FROZEN, LC_PREP, LC_BAND,
                LC_MS1BAND, LC_PARENT, LC_WSTORE, LC_FRAG, LC_SPECTRA, LC_CARRY, LC_WRITER,
                LC_SCRATCH, LC_N };
  inline const char* ledName_(int c)
  {
    // band_arena: the MS2 integer path's per-window arenas; ms1_band: detectTraces_'s copy of the MS1 map.
    static const char* const n[LC_N] = {"slab", "ms1_map", "ms1_trace", "precursor", "frozen", "prep",
                                        "band_arena", "ms1_band", "parent", "wstore", "frag_aux",
                                        "spectra", "carry", "writer", "scratch"};
    return n[c];
  }
  inline std::array<std::atomic<long long>, LC_N>& ledger_() { static std::array<std::atomic<long long>, LC_N> a{}; return a; }
  inline const char* ledFile_() { static const char* f = std::getenv("DIASPEXTRACT_LEDGER"); return f; }
  inline std::atomic<const char*>& ledPhaseName_() { static std::atomic<const char*> p{"start"}; return p; }   // stored pointers need static lifetime: the sampler reads them from another thread
  inline void ledAdd_(int c, long long b) { if (ledFile_()) ledger_()[(size_t)c].fetch_add(b, std::memory_order_relaxed); }
  inline void ledSet_(int c, long long b) { if (ledFile_()) ledger_()[(size_t)c].store(b, std::memory_order_relaxed); }
  inline void ledPhaseTile_(const char* what, size_t k)
  {
    if (!ledFile_()) return;
    static std::mutex m; static std::string buf;
    std::lock_guard<std::mutex> g(m);
    buf = std::string(what) + std::to_string(k);
    ledPhaseName_().store(buf.c_str(), std::memory_order_relaxed);
  }
  /// A scoped charge: the subtract cannot be skipped by an early return or an exception.
  struct LedCharge
  {
    int c; long long b;
    LedCharge(int cat, long long bytes) : c(cat), b(bytes) { ledAdd_(c, b); }
    void change(long long nb) { ledAdd_(c, nb - b); b = nb; }
    ~LedCharge() { ledAdd_(c, -b); }
    LedCharge(const LedCharge&) = delete;
    LedCharge& operator=(const LedCharge&) = delete;
  };
  /// [ledger] live bytes of a PeakMap: peaks at the padded 20 B (Peak1D 16 + the IM float) plus
  /// the per-spectrum header and its settings tail. O(spectra), never O(peaks).
  inline long long ledMapBytes_(const MSExperiment& m)
  {
    long long b = (long long)m.size() * 700;
    for (const auto& sp : m) b += (long long)sp.size() * 20;
    return b;
  }
  /// [ms1-hist] DIASPEXTRACT_MS1_HIST=1: picked MS1 peaks counted by intensity band (<100, 100-300, 300-1000,
  /// 1000-10000, >=10000; 100 and 300 are MS1 tracing's default noise and seed thresholds). Counting only.
  inline bool ms1HistOn_() { static const bool on = []{ const char* e = std::getenv("DIASPEXTRACT_MS1_HIST"); return e && *e && *e != '0'; }(); return on; }
  inline std::array<std::atomic<long long>, 5>& ms1Hist_() { static std::array<std::atomic<long long>, 5> h{}; return h; }
  /// [pspan] precursors reaching assembly by MS1 span in window frames: 1, 2, 3, 4, 5-6, 7-10, 11-20, >20; [8] = faint (adaptive)
  inline std::array<std::atomic<long long>, 9>& pspanHist_() { static std::array<std::atomic<long long>, 9> h{}; return h; }
  /// [merge] [0] picked MS2 peaks merged away in toTof (one frame, one bin, 1/K0 within trace:im_slice_merge);
  /// [1] frames that had such a merge; [2] fragments merged away at assembly (one bin in one spectrum)
  inline std::array<std::atomic<long long>, 3>& mergeCount_() { static std::array<std::atomic<long long>, 3> c{}; return c; }
  /// [fragint] [0] gate-passing (and gate:anchor_rescue rescued) fragments given the core/sum amplitude; [1] those of
  /// them with no point there (fallback)
  inline std::array<std::atomic<long long>, 2>& fragIntCount_() { static std::array<std::atomic<long long>, 2> c{}; return c; }
  /// [anchor] gate:anchor_rescue: [0] non-faint precursors scored, [1] of them anchored, [2] candidates of the anchored,
  /// [3] rescued, [4] passed beyond gate:anchor_max_rescued, [5] precursors with a rescue, [6] rescued peaks in written
  /// spectra (after the same-bin merge and the cap, before the twin fold), [7] written spectra holding one, [8..10] rescued
  /// fragments the gate had rejected at the point floor, the overlap guard, the correlation cut
  inline std::array<std::atomic<long long>, 11>& anchorCount_() { static std::array<std::atomic<long long>, 11> c{}; return c; }
  /// [foreign] per assembly (one precursor in one window): [8] assemblies checked (with fragments), [0] of them with a
  /// foreign candidate, [1] foreign envelopes matched, [2] assemblies with a fragment attenuated, [3] fragments
  /// attenuated, [4] candidates skipped as this precursor's own ion, [5] assemblies of [2] with at least
  /// assembly:min_fragments before the removal (weight 0) and fewer after it and the same-bin merge, [6] matched
  /// fragments' signed m/z error summed in 0.001 ppm, [7] matches behind [6]
  inline std::array<std::atomic<long long>, 9>& foreignCount_() { static std::array<std::atomic<long long>, 9> c{}; return c; }
  /// [fuse] assembly:fragment_fuse_ppm, per assembly (one precursor in one window) with at least one fragment: [0]
  /// assemblies, [1] fragments in, [2] peaks out, [3] groups of >= 2 fused, [4] fragments in those groups, [5] assemblies
  /// with at least assembly:min_fragments before the fusion and fewer after it
  inline std::array<std::atomic<long long>, 6>& fuseCount_() { static std::array<std::atomic<long long>, 6> c{}; return c; }
  /// RAII in main_: one sampler for the run, joined before any return path leaves the function.
  struct LedRun
  {
    std::atomic<bool> stop{false};
    std::thread t;
    LedRun()
    {
      if (!ledFile_()) return;
      t = std::thread([this]() {
        std::ofstream f(ledFile_());
        // MiB: every value below is a byte count >> 20.
        f << "t_s\trss_mib\tarena_mib\tfree_mib\tsum_mib\tunattributed_mib";
        for (int c = 0; c < LC_N; ++c) f << '\t' << ledName_(c) << "_mib";
        f << "\tphase\n";
        while (true)
        {
          const bool last = stop.load(std::memory_order_relaxed);
          long long arena = 0, freeb = 0;
#ifdef __GLIBC__
          { const struct mallinfo2 mi = mallinfo2(); arena = (long long)(mi.arena >> 20); freeb = (long long)(mi.fordblks >> 20); }
#endif
          const long rss = rss_mb_();
          long long sum = 0, v[LC_N];
          for (int c = 0; c < LC_N; ++c) { v[c] = ledger_()[(size_t)c].load(std::memory_order_relaxed); sum += v[c]; }
          const long long sum_mb = sum >> 20;
          char buf[64];
          std::snprintf(buf, sizeof buf, "%.2f", phase_clock_());
          f << buf << '\t' << rss << '\t' << arena << '\t' << freeb << '\t' << sum_mb << '\t' << (rss - sum_mb);
          for (int c = 0; c < LC_N; ++c) f << '\t' << (v[c] >> 20);
          f << '\t' << ledPhaseName_().load(std::memory_order_relaxed) << '\n';
          f.flush();
          if (last) break;
          std::this_thread::sleep_for(std::chrono::milliseconds(250));
        }
      });
    }
    ~LedRun()
    {
      if (!t.joinable()) return;
      ledPhaseName_().store("end", std::memory_order_relaxed); stop.store(true, std::memory_order_relaxed); t.join();
    }
  };
  /// " [RSS cur, peak hwm]" from /proc/self/status; "" where unsupported. [mem]
  static std::string rss_()
  {
    std::ifstream st("/proc/self/status");
    if (!st) return "";
    std::string line, cur = "?", peak = "?";
    while (std::getline(st, line))
    {
      auto val = [&line]() {
        std::string v = line.substr(line.find(':') + 1);
        size_t a = v.find_first_not_of(" \t");
        return a == std::string::npos ? std::string("?") : v.substr(a);
      };
      if (line.rfind("VmRSS:", 0) == 0) cur = val();
      else if (line.rfind("VmHWM:", 0) == 0) peak = val();
    }
    return " [RSS " + cur + ", peak " + peak + "]";
  }

  /// [dyn-mem] MemAvailable in bytes: the kernel's estimate of what can be allocated without swapping on a
  /// shared node. Our own RSS is NOT added back: after perf:malloc_trim it is the windows in flight, memory
  /// already spoken for. 0 if unreadable (the caller then admits one window at a time).
  static size_t availableBytes_()
  {
    std::ifstream mi("/proc/meminfo");
    std::string line;
    size_t avail_kb = 0;
    while (std::getline(mi, line))
      if (line.rfind("MemAvailable:", 0) == 0)
      { avail_kb = strtoull(line.c_str() + 13, nullptr, 10); break; }
    return avail_kb * 1024ull;
  }

  /// The global RT axis: every distinct frame retention time (MS1 and MS2), sorted. No RT is stored per
  /// point: each TraceStore maps its frames to indices on this axis (rt_index). Built once after loading,
  /// read-only thereafter.
  inline vector<double>& rtAxis() { static vector<double> v; return v; }

  /// Index of the element of sorted @p a nearest to @p x (a tie keeps the later one).
  inline uint32_t nearestIndex(const vector<double>& a, double x)
  {
    auto it = lower_bound(a.begin(), a.end(), x);
    if (it == a.end()) return (uint32_t)(a.size() - 1);
    if (it != a.begin() && (*it - x) > (x - *(it - 1))) --it;
    return (uint32_t)(it - a.begin());
  }

  /// Frame index of an RT. Used only while BUILDING traces: the value comes from a frame, so the
  /// match is exact, and the nearest-neighbour fallback only guards a last-ulp difference.
  inline uint32_t rtIndex(double rt) { return rtAxis().empty() ? 0 : nearestIndex(rtAxis(), rt); }

  struct Trace
  {
    double mz = 0.0;
    double im = 0.0;   ///< centroid ion mobility (1/K0)
    /// The flight-time bin this trace was measured on, and the calibration factor of its apex
    /// frame. The pair is the AUTHORITATIVE m/z: `mz` below is a cached convenience for the gates
    /// that still compare in m/z space, and the value actually written out is recomputed from these
    /// two at export. 0 means the trace came from the OpenMS detector, which has no bin.
    uint32_t tof = 0;
    uint32_t bframe = 0;   ///< store-local frame whose factor calibrates `tof` (see bOf())
    /// The profile is a span of the owning TraceStore's arena, [off, off+len), one slot per frame
    /// from `frame0` (store-local: a window's frames are not contiguous on the global RT axis).
    /// ZERO MEANS MISSING, valid because every observed intensity is > noise_threshold >= 0;
    /// consumers walk the real points in frame order.
    uint32_t frame0 = 0;   ///< first frame of the span, store-local
    uint32_t off = 0;      ///< offset of the span in the store's arena
    uint16_t len = 0;      ///< frames spanned (a gradient is ~1,400; asserted < 65536)
    uint16_t npts = 0;     ///< REAL points in the span
    uint16_t apex = 0;     ///< position of the apex within the span
    /// [soa-rt] The RT is a frame index. At the shipped defaults every producer reports one of the store's
    /// frame times (the integer detector its apex frame, ElutionPeakDetection its smoothed-maximum point);
    /// an off-grid RT (EPD skipped) is snapped and counted, see rtOffGrid(). Read it with rtOf(); never
    /// derive it from `apex` (the RAW maximum, a different frame). SIGNED: trimToSpan moves `frame0` by
    /// +lo and this by -lo, so the reported RT survives any trim.
    int16_t rt_at = 0;     ///< frame whose time is this trace's RT, relative to frame0 (see rtOf())
    size_t np() const { return npts; }
    size_t span() const { return len; }
    void freeProfile() { len = 0; npts = 0; }        ///< the arena is the store's; nothing to free
  };

  /// The frame table and the arenas that a set of traces point into: one per window, one for the
  /// MS1 map. `rt_index`/`b` are per FRAME (store-local index); `inten` holds every trace's span
  /// back to back (zero = missing); `bins` holds the per-frame flight-time bin for the integer
  /// path only and is freed once valley splitting has assigned each child its apex bin.
  struct TraceStore
  {
    vector<uint32_t> rt_index;   ///< frame -> global RT axis
    vector<double>   b;          ///< frame -> calibration factor (empty on the OpenMS path)
    vector<double>   frame_rt;   ///< frame -> RT value, for toTrace()'s frame lookup
    vector<float>    inten;      ///< arena of spans
    /// [bin-delta] A point's flight-time bin as the SIGNED offset from its trace's own bin
    /// (`Trace::tof`): a trace spans a few bins, int16 leaves ~30,000 of headroom, and the writer
    /// throws rather than truncate. Reconstruct with `t.tof + bins[t.off + k]`; a binless trace
    /// (tof 0, zero-filled slots) reconstructs to 0.
    vector<int16_t> bins;        ///< arena of per-frame bin OFFSETS (integer path, until EPD)
    size_t frames() const { return rt_index.size(); }
    long long bytes() const { return (long long)(inten.capacity() * 4 + bins.capacity() * 2 + rt_index.capacity() * 4 + b.capacity() * 8 + frame_rt.capacity() * 8); }
    void setFrames(const vector<uint32_t>& ri, const vector<double>* bb)
    {
      rt_index = ri; if (bb) b = *bb; else b.clear();
      frame_rt.resize(ri.size());
      for (size_t f = 0; f < ri.size(); ++f) frame_rt[f] = rtAxis()[ri[f]];
    }
    /// Store-local frame of an RT value (exact by construction; nearest guards a last-ulp slip).
    uint32_t frameOf(double rt) const { return nearestIndex(frame_rt, rt); }
    /// Append `other`'s arenas and return the offset they now start at (for rebasing).
    uint32_t absorb(TraceStore& other)
    {
      const size_t base = inten.size();
      if (base + other.inten.size() >= (size_t)std::numeric_limits<uint32_t>::max())
        throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      inten.insert(inten.end(), other.inten.begin(), other.inten.end());
      if (!other.bins.empty()) { bins.resize(base, 0); bins.insert(bins.end(), other.bins.begin(), other.bins.end()); }
      // Once ANY source has contributed bins, keep bins index-parallel to inten: a later binless
      // source would otherwise leave bins short, and a bins[t.off + k] read would run off the end.
      // No current caller mixes the two, which is exactly why this is worth pinning down here.
      else if (!bins.empty()) bins.resize(inten.size(), 0);
      vector<float>().swap(other.inten); vector<int16_t>().swap(other.bins);
      return (uint32_t)base;
    }
  };
  static_assert(sizeof(Trace) == 40, "Trace must stay 40 bytes: it is the largest structure in the tool "
                                    "(2.1e9 records on a 2-h file), and a silent 8-byte regrowth costs ~16 GiB");

  /// [serial] `dst.absorb()` of every source in order, each source's records rebased (off += its base; frame0 and
  /// bframe += its frame shift) and appended to `out`: the serial loop's result byte for byte (sizes and capacities
  /// included), with the copies run as bounded tasks at prefix-summed offsets. The parallel form needs an empty
  /// `dst` and `out` and every source's bins index-parallel to its inten (or absent) -- what every caller has; any
  /// other input runs the serial loop. The sources are emptied either way.
  inline void absorbAll_(TraceStore& dst, vector<Trace>& out, const vector<std::pair<TraceStore*, vector<Trace>*>>& src,
                         const vector<uint32_t>& fshift)
  {
    const size_t ns = src.size();
    vector<size_t> ib(ns + 1, 0), tb(ns + 1, 0);   // arena / record prefix sums
    bool par = dst.inten.empty() && dst.bins.empty() && out.empty(), any_bins = false;
    for (size_t s = 0; s < ns; ++s)
    {
      const TraceStore& o = *src[s].first;
      if (!o.bins.empty()) { any_bins = true; par &= (o.bins.size() == o.inten.size()); }
      ib[s + 1] = ib[s] + o.inten.size(); tb[s + 1] = tb[s] + src[s].second->size();
    }
    if (!par)
    {
      out.reserve(out.size() + tb[ns]);
      for (size_t s = 0; s < ns; ++s)
      {
        const uint32_t base = dst.absorb(*src[s].first), fs = fshift[s];
        for (auto& t : *src[s].second) { t.off += base; t.frame0 += fs; t.bframe += fs; out.push_back(std::move(t)); }
        vector<Trace>().swap(*src[s].second);
      }
      return;
    }
    // absorb() throws at the first source that would reach 2^32 - 1 floats; the prefix sums are monotone
    if (ib[ns] >= (size_t)std::numeric_limits<uint32_t>::max())
      throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // absorb()'s bins: once any source carries bins they are index-parallel to inten, zero where a source has none
    dst.inten.resize(ib[ns]); if (any_bins) dst.bins.resize(ib[ns], 0);
    out.resize(tb[ns]);
    // call f(s, a, e) for the part [a, e) of source s that falls in slice c of nt of the concatenation `pre`
    auto slices = [ns](const vector<size_t>& pre, long c, long nt, auto&& f) {
      const size_t N = pre[ns], lo = (size_t)c * N / (size_t)nt, hi = (size_t)(c + 1) * N / (size_t)nt;
      size_t s = (size_t)(std::upper_bound(pre.begin(), pre.end(), lo) - pre.begin()) - 1;
      for (; s < ns && pre[s] < hi; ++s)
      { const size_t a = std::max(lo, pre[s]) - pre[s], e = std::min(hi, pre[s + 1]) - pre[s]; if (a < e) f(s, a, e); }
    };
    const long nt = parTasks_(ib[ns] + 8 * tb[ns], 1 << 20);   // a record is ~8 floats of copying
    forTasks_(nt, [&](long c) {
      slices(ib, c, nt, [&](size_t s, size_t a, size_t e) {
        const TraceStore& o = *src[s].first;
        std::copy(o.inten.begin() + a, o.inten.begin() + e, dst.inten.begin() + ib[s] + a);
        if (!o.bins.empty()) std::copy(o.bins.begin() + a, o.bins.begin() + e, dst.bins.begin() + ib[s] + a);
      });
      slices(tb, c, nt, [&](size_t s, size_t a, size_t e) {
        const uint32_t base = (uint32_t)ib[s], fs = fshift[s];
        const vector<Trace>& r = *src[s].second;
        for (size_t q = a; q < e; ++q) { Trace t = r[q]; t.off += base; t.frame0 += fs; t.bframe += fs; out[tb[s] + q] = t; }
      });
    });
    for (size_t s = 0; s < ns; ++s)
    { vector<float>().swap(src[s].first->inten); vector<int16_t>().swap(src[s].first->bins); vector<Trace>().swap(*src[s].second); }
  }

  /// Span accessors take the owning store explicitly (no per-record back-pointer); bounds asserted in debug builds.
  inline float xv(const TraceStore& s, const Trace& t, size_t k)
  { assert(t.off + k < s.inten.size()); return s.inten[t.off + k]; }
  inline bool real(const TraceStore& s, const Trace& t, size_t k) { return xv(s, t, k) > 0.0f; }
  inline uint32_t gframe(const TraceStore& s, const Trace& t, size_t k)
  { assert(t.frame0 + k < s.rt_index.size()); return s.rt_index[t.frame0 + k]; }
  inline double rtAtSpan(const TraceStore& s, const Trace& t, size_t k) { return rtAxis()[gframe(s, t, k)]; }
  /// [soa-rt] The trace's reported RT. `frame_rt` is one double per frame of the owning store
  /// (~38 KB for a window), so this is a single load from a table that stays hot, against the
  /// 8 bytes per record it replaces. Identical by construction: frame_rt[f] IS rtAxis()[ri[f]],
  /// the same array element the producers read.
  inline double rtOf(const TraceStore& s, const Trace& t)
  { const int64_t f = (int64_t)t.frame0 + t.rt_at;
    assert(f >= 0 && (size_t)f < s.frame_rt.size()); return s.frame_rt[(size_t)f]; }
  /// Traces whose RT is NOT one of the store's frame times: only reachable when
  /// ElutionPeakDetection is skipped (trace:ms{1,2}_split_valleys = 0, or a single trace), where
  /// MassTraceDetection reports an intensity-weighted mean instead. Counted, not silent.
  inline std::atomic<size_t>& rtOffGrid() { static std::atomic<size_t> n{0}; return n; }
  /// [soa-rt diag] is the stray RT strictly BETWEEN two frame times (a weighted
  /// mean, which no frame index can represent) or merely a hair off one (a representable value)?
  inline std::atomic<size_t>& rtStrictlyBetween() { static std::atomic<size_t> n{0}; return n; }
  /// Record `rt` as a frame offset. Snaps to the nearest frame if the value is off-grid (counted).
  inline void setRtFrame(const TraceStore& s, Trace& t, double rt)
  {
    if (s.frame_rt.empty()) { t.rt_at = 0; return; }
    const uint32_t f = s.frameOf(rt);
    if (s.frame_rt[f] != rt)
    {
      rtOffGrid().fetch_add(1, std::memory_order_relaxed);
      // strictly interior to the frame grid => an interpolated value, unrepresentable as an index
      if (f > 0 && f + 1 < s.frame_rt.size() && rt > s.frame_rt[f - 1] && rt < s.frame_rt[f + 1])
        rtStrictlyBetween().fetch_add(1, std::memory_order_relaxed);
    }
    const int64_t d = (int64_t)f - (int64_t)t.frame0;
    if (d < -32768 || d > 32767) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    t.rt_at = (int16_t)d;
  }
  /// The calibration factor of the trace's own apex/measurement frame. Empty table (the OpenMS
  /// detector has no per-frame factors) yields 0.0, which is what the field held there before, and
  /// exportMz_ returns the cached m/z for those traces because their `tof` is 0.
  inline double bOf(const TraceStore& s, const Trace& t) { return s.b.empty() ? 0.0 : s.b[t.bframe]; }
  /// The trace's apex intensity: the arena value at `apex` (makeSpan's maximum; trimming moves `apex`
  /// with the data). 0 for a released profile (len == 0), which is never scored.
  inline double intensityOf(const TraceStore& s, const Trace& t)
  { return t.len ? (double)s.inten[t.off + t.apex] : 0.0; }

  /// The real points of a trace, packed in frame order, for the consumers that were written over
  /// point lists (smoother, FWHM): identical sequences, so identical arithmetic.
  inline void packReal(const TraceStore& s, const Trace& t, vector<double>& rt, vector<double>& v)
  {
    rt.clear(); v.clear();
    for (size_t k = 0; k < t.span(); ++k)
      if (real(s, t, k)) { rt.push_back(rtAtSpan(s, t, k)); v.push_back((double)xv(s, t, k)); }
  }

  /// [fwhm] FWHM (seconds) from the trace's own XIC: apex, then outwards to half maximum with linear interpolation;
  /// a side that never reaches half maximum ends at the first/last point. -1 below 3 points.
  inline double xicFwhm(const TraceStore& s, const Trace& tr)
  {
    if (tr.np() < 3) return -1.0;
    vector<double> prt, pv; packReal(s, tr, prt, pv);
    size_t ap = 0;
    for (size_t i = 1; i < pv.size(); ++i) if (pv[i] > pv[ap]) ap = i;
    const double half = pv[ap] * 0.5;
    double lo = prt[0], hi = prt[pv.size() - 1];
    for (size_t i = ap; i > 0; --i)
      if (pv[i - 1] <= half)
      {
        const double d = pv[i] - pv[i - 1];
        const double f = d > 0 ? (pv[i] - half) / d : 0.0;
        lo = prt[i] - f * (prt[i] - prt[i - 1]);
        break;
      }
    for (size_t i = ap; i + 1 < pv.size(); ++i)
      if (pv[i + 1] <= half)
      {
        const double d = pv[i] - pv[i + 1];
        const double f = d > 0 ? (pv[i] - half) / d : 0.0;
        hi = prt[i] + f * (prt[i + 1] - prt[i]);
        break;
      }
    return hi - lo;
  }

  /// Trim a trace to at most `cap` seconds around its apex, AFTER detection and valley splitting, so
  /// peak ownership is unchanged (rationale: trace:max_span_sec). Only the span fields change and no
  /// bytes are copied; the arena keeps the trimmed bytes until a compaction pass.
  inline void trimToSpan(const TraceStore& s, Trace& t, double cap)
  {
    if (cap <= 0.0 || t.span() == 0) return;
    const size_t a = t.apex;
    size_t lo = a, hi = a;
    // grow symmetrically outward from the apex while the span fits, preferring the earlier frame on
    // a tie so the result does not depend on which side is tested first
    for (;;)
    {
      const bool can_lo = lo > 0, can_hi = hi + 1 < t.span();
      if (!can_lo && !can_hi) break;
      const double d_lo = can_lo ? rtAtSpan(s, t, hi) - rtAtSpan(s, t, lo - 1) : 1e30;
      const double d_hi = can_hi ? rtAtSpan(s, t, hi + 1) - rtAtSpan(s, t, lo) : 1e30;
      if (d_lo <= d_hi) { if (!can_lo || d_lo > cap) { if (!can_hi || d_hi > cap) break; ++hi; } else --lo; }
      else              { if (!can_hi || d_hi > cap) { if (!can_lo || d_lo > cap) break; --lo; } else ++hi; }
    }
    // Zero-trim the ends. makeSpan starts and ends a span on a real point by construction, but
    // growing outward from the apex can stop on a gap, and a leading or trailing zero is pure
    // padding: it costs 4 B, it is skipped by every consumer, and it makes `span()` overstate the
    // trace. Interior zeros are the gaps and MUST stay -- they are the presence encoding.
    while (lo < hi && xv(s, t, lo) == 0.0f) ++lo;
    while (hi > lo && xv(s, t, hi) == 0.0f) --hi;
    if (lo == 0 && hi + 1 == t.span()) return;                  // already within the cap, ends real
    uint16_t np = 0;
    for (size_t k = lo; k <= hi; ++k) if (xv(s, t, k) > 0.0f) ++np;
    // [soa-rt] rt_at names EPD's smoothed apex, not `apex`, and may lie outside the kept span. It is
    // not clamped (a clamp moved the digest, docs/BASELINE.md): frame0 and rt_at move by opposite
    // amounts, so the reported RT is unchanged.
    const long rr = (long)t.rt_at - (long)lo;
    if (rr < -32768L || rr > 32767L) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    t.frame0 += (uint32_t)lo; t.off += (uint32_t)lo;
    t.len = (uint16_t)(hi - lo + 1); t.apex = (uint16_t)(a - lo); t.npts = np;
    t.rt_at = (int16_t)rr;
  }

  /// Drop the spans nobody references and compact the arena IN PLACE in OFFSET order (ms1_traces is
  /// m/z-sorted, so container order would overwrite spans not yet copied). Spans must be disjoint --
  /// validated before a byte moves. A freed trace keeps its scalars; returns the profiles released.
  inline Size compactUnreferenced(vector<Trace>& tr, TraceStore& store, const vector<bool>& needed)
  {
    vector<uint32_t> ids;
    Size freed = 0;
    for (size_t i = 0; i < tr.size(); ++i)
    {
      if (!needed[i]) { if (tr[i].np() > 0) ++freed; tr[i].freeProfile(); continue; }
      if (tr[i].len > 0) ids.push_back((uint32_t)i);     // a zero-length span owns no bytes
    }
    // [outer] valid spans have distinct offsets (equal ones overlap and throw below either way), so with the index as a
    // last key the order is total and the parallel sort returns std::sort's
    parallelSortTotal_(ids, [&](uint32_t a, uint32_t b) { return tr[a].off != tr[b].off ? tr[a].off < tr[b].off : a < b; });
    uint64_t prev_end = 0;
    for (uint32_t i : ids)
    {
      const uint64_t beg = tr[i].off, end = beg + (uint64_t)tr[i].len;
      if (end > store.inten.size() || beg < prev_end)
        throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "MS1 spans overlap or exceed the arena");
      prev_end = end;
    }
    size_t w = 0;
    for (uint32_t i : ids)
    {
      Trace& t = tr[i];
      if (w != t.off) std::memmove(store.inten.data() + w, store.inten.data() + t.off, (size_t)t.len * sizeof(float));
      t.off = (uint32_t)w; w += t.len;
    }
    store.inten.resize(w); store.inten.shrink_to_fit();
    return freed;
  }

  /// Build a span from (store-local frame, intensity) points into `store`, returning the trace
  /// with `st` UNSET: the caller sets it once the store is immutable. Points may arrive unsorted.
  inline Trace makeSpan(vector<pair<uint32_t, float>>& pts, TraceStore& store)
  {
    Trace t;
    if (pts.empty()) return t;
    sort(pts.begin(), pts.end());
    const uint32_t f0 = pts.front().first, f1 = pts.back().first;
    const size_t L = (size_t)f1 - f0 + 1;
    if (L >= 65536) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    t.frame0 = f0; t.len = (uint16_t)L; t.off = (uint32_t)store.inten.size();
    store.inten.resize(store.inten.size() + L, 0.0f);
    float best = -1.0f; uint16_t np = 0;
    for (const auto& q : pts)
    {
      const size_t k = q.first - f0;
      if (store.inten[t.off + k] == 0.0f) ++np;      // one point per frame; a repeat overwrites
      store.inten[t.off + k] = q.second;
      if (q.second > best) { best = q.second; t.apex = (uint16_t)k; }
    }
    t.npts = np;
    return t;
  }

  /// A precursor hypothesis after isotope/charge inference.
  struct Precursor_
  {
    double mono_mz = 0.0;
    int charge = 0;          ///< 0 == unknown
    double rt = 0.0;
    double im = 0.0;
    size_t trace_idx = 0;    ///< index into ms1 traces
    bool guessed = false;    ///< [way-4] charge was DEFAULTED, not isotope-supported
    bool mono_ext = false;   ///< [twin-rt] the monoisotope was placed by monocheck v2's extension (B)
    int n_isotopes = 0;      ///< [Q4] isotope peaks behind the mono call (0/1 = weak, mono ±1.00335 uncertain)
  };

  /// [monocheck-v2] What charge:mono_profile_check v2 needs of a count walk, one per precursor of inferPrecursors_
  /// (z = 0: nothing to test): the winning charge, its leftward partners' m/z (nearest first), the monoisotope before
  /// v2, and whether the leftward walk stopped at charge:mono_averagine_guard instead of on a miss.
  struct MonoV2In_ { int z = 0, nl = 0; bool guard_stop = false; double mono = 0.0; std::array<double, 5> left_mz{}; };

  /// MassTraceDetection locates the per-peak IM array by the exact name
  /// Constants::UserParam::ION_MOBILITY ("Ion Mobility"). Converters name it via CV term and
  /// PeakPickerIM names it "Ion Mobility Centroid"; rename whichever exists so tracing is IM-aware.
  void ensureIMArrayName(MSSpectrum& s)
  {
    auto& fdas = s.getFloatDataArrays();
    if (s.containsIMData())
    {
      fdas[s.getIMData().first].setName(Constants::UserParam::ION_MOBILITY);
      return;
    }
    for (auto& fda : fdas)
    {
      if (fda.getName() == Constants::UserParam::ION_MOBILITY_CENTROID)
      {
        fda.setName(Constants::UserParam::ION_MOBILITY);
        return;
      }
    }
  }

  /// Reported-m/z estimator, set once from `trace:mz_estimator` before any trace is converted.
  /// A file-scope value because toTrace() is a free helper with no access to the tool's parameters.
  inline std::string& mzEstimator() { static std::string v = "apex"; return v; }

  Trace toTrace(const MassTrace& mt, TraceStore& store)
  {
    // Reported m/z of a mass trace. "apex" (the default) takes the m/z of the most intense peak;
    // "mean" is the OpenMS intensity-weighted centroid; "median" the median m/z. Set by
    // `trace:mz_estimator` -- see docs/BASELINE.md for the comparison that chose apex.
    const std::string& mz_est = mzEstimator();
    double mz;
    if (mz_est == "apex" && mt.getSize() > 0) mz = mt[mt.findMaxByIntPeak(false)].getMZ();
    else if (mz_est == "median" && mt.getSize() > 0) { MassTrace m2 = mt; m2.updateMedianMZ(); mz = m2.getCentroidMZ(); }
    else mz = mt.getCentroidMZ();
    // One peak per frame (MassTraceDetection takes at most one candidate per spectrum), so the
    // profile is a span over the store's frame table with zeros where frames were missed.
    const Size np = mt.getSize();
    vector<pair<uint32_t, float>> pts;
    pts.reserve(np);
    for (Size i = 0; i < np; ++i) pts.emplace_back(store.frameOf(mt[i].getRT()), (float)mt[i].getIntensity());
    Trace t = makeSpan(pts, store);
    t.mz = mz; t.im = mt.getCentroidIM();
    setRtFrame(store, t, mt.getCentroidRT());   // [soa-rt] a frame time, stored as an offset
    return t;
  }

  /// [outer] The span length makeSpan gives @p mt: its points' frame range, 0 without points. Throws as makeSpan does.
  inline size_t spanLength(const MassTrace& mt, const TraceStore& store)
  {
    const Size np = mt.getSize();
    if (np == 0) return 0;
    uint32_t f0 = std::numeric_limits<uint32_t>::max(), f1 = 0;
    for (Size i = 0; i < np; ++i) { const uint32_t f = store.frameOf(mt[i].getRT()); f0 = std::min(f0, f); f1 = std::max(f1, f); }
    const size_t L = (size_t)f1 - f0 + 1;
    if (L >= 65536) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    return L;
  }

  /// [outer] toTrace with the span placed at arena offset @p off (already sized and zero) instead of appended: the same
  /// statements as toTrace + makeSpan, so the same Trace and the same arena bytes when @p off is where the append would
  /// have put it. @p pts is the caller's scratch.
  inline Trace toTraceAt(const MassTrace& mt, TraceStore& store, size_t off, vector<pair<uint32_t, float>>& pts)
  {
    const std::string& mz_est = mzEstimator();
    double mz;
    if (mz_est == "apex" && mt.getSize() > 0) mz = mt[mt.findMaxByIntPeak(false)].getMZ();
    else if (mz_est == "median" && mt.getSize() > 0) { MassTrace m2 = mt; m2.updateMedianMZ(); mz = m2.getCentroidMZ(); }
    else mz = mt.getCentroidMZ();
    const Size np = mt.getSize();
    pts.clear();
    for (Size i = 0; i < np; ++i) pts.emplace_back(store.frameOf(mt[i].getRT()), (float)mt[i].getIntensity());
    Trace t;
    if (!pts.empty())
    {
      sort(pts.begin(), pts.end());
      const uint32_t f0 = pts.front().first, f1 = pts.back().first;
      const size_t L = (size_t)f1 - f0 + 1;
      if (L >= 65536) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      t.frame0 = f0; t.len = (uint16_t)L; t.off = (uint32_t)off;
      float best = -1.0f; uint16_t npt = 0;
      for (const auto& q : pts)
      {
        const size_t k = q.first - f0;
        if (store.inten[t.off + k] == 0.0f) ++npt;      // one point per frame; a repeat overwrites
        store.inten[t.off + k] = q.second;
        if (q.second > best) { best = q.second; t.apex = (uint16_t)k; }
      }
      t.npts = npt;
    }
    t.mz = mz; t.im = mt.getCentroidIM();
    setRtFrame(store, t, mt.getCentroidRT());   // [soa-rt] a frame time, stored as an offset
    return t;
  }


  /// [compact] Quantised per-peak store, 10 B/peak, far finer than every downstream tolerance:
  ///   m/z : uint32 at 1e-5 Da     (0.007 ppm at m/z 1400)
  ///   1/K0: uint16 over [0.4,1.8] (2.1e-5)
  ///   intensity: float32, unchanged
  constexpr double MZ_Q = 1e5;      ///< m/z quantum = 1e-5 Da

  /// The m/z axis is the instrument's flight-time bin: lossless, the same bin for the same ion in every
  /// frame, and a ppm tolerance is a handful of integer bins. The per-frame factor `b` (digitizer
  /// temperature) turns a bin into m/z at report time. See docs/MZ-AXIS-DESIGN.md.
  using TofIdx = uint32_t;

  /// The flight-time axis for one run: the calibration model plus each frame's temperature factor.
  /// Fail-closed: without a calibration there is no TOF axis and the caller must not pretend there
  /// is one.
  struct TofAxis
  {
    diaspextract::TdfMzCalibration cal;
    vector<double> b_by_frame;        ///< frame factor per FRAME ID (index 0 unused, as in the tdf)
    bool ok = false;
    double acq_lo = std::numeric_limits<double>::quiet_NaN();   ///< the acquisition's m/z range (GlobalMetadata)
    double acq_hi = std::numeric_limits<double>::quiet_NaN();
    long n_bins = 0;                                            ///< DigitizerNumSamples: every bin index is below it
    String bounds_why;                                          ///< why the bounds are missing, if they are

    double factor(size_t frame_id) const
    { return frame_id < b_by_frame.size() ? b_by_frame[frame_id] : cal.frameFactor(cal.T1_ref); }
    double mzOf(TofIdx tof, double b) const { return cal.tofToMz((double)tof, b); }
    TofIdx tofOf(double mz, double b) const
    {
      const double t = cal.mzToTof(mz, b);
      if (!(t > 0.0)) return 0u;
      return t >= 4294967294.0 ? 4294967294u : (TofIdx)llround(t);   // saturate: a wrap here would collapse a band partition silently
    }
  };

  inline TofAxis& tofAxis() { static TofAxis a; return a; }

  /// Build the flight-time axis from a calibration and the per-frame T1 (index = Frames.Id).
  inline void setTofAxis(const diaspextract::TdfMzCalibration& cal, const std::vector<double>& t1)
  {
    TofAxis& ax = tofAxis();
    ax.cal = cal;
    // A frame without a positive factor gets the reference one. Index 0 is the "unknown frame"
    // sentinel (Frames.Id starts at 1; frameIdOf() returns 0 when it cannot parse) and stays 0.0.
    const double bref = ax.cal.frameFactor(ax.cal.T1_ref);
    ax.b_by_frame.assign(t1.size(), 0.0);
    for (size_t i = 1; i < t1.size(); ++i) { const double f = ax.cal.frameFactor(t1[i]); ax.b_by_frame[i] = f > 0.0 ? f : bref; }
    ax.ok = true;
  }
  /// Load the axis from an analysis.tdf; false (axis unusable) when anything is missing.
  inline bool loadTofAxis(const std::string& tdf_path, String& why)
  {
    diaspextract::TdfMzCalibration cal;
    std::vector<double> t1;
    if (!diaspextract::loadTdfCalibration(tdf_path, cal, t1, why)) return false;
    setTofAxis(cal, t1);
    { std::string w2; diaspextract::TdfAxisBounds bd;
      if (diaspextract::loadTdfAxisBounds(tdf_path, bd, w2)) { tofAxis().acq_lo = bd.mz_lo; tofAxis().acq_hi = bd.mz_hi; tofAxis().n_bins = bd.n_bins; }
      else tofAxis().bounds_why = w2; }
    return true;
  }
  constexpr double IM_LO = 0.4, IM_HI = 1.8;
  constexpr double IM_Q = 65535.0 / (IM_HI - IM_LO);

  /// Peaks that cannot be represented exactly are DROPPED and COUNTED: clamping 1/K0 would stack
  /// outliers on two artificial mobilities, and NaN into llround is UB.
  struct CompactStats
  {
    size_t no_im_array = 0, size_mismatch = 0, bad_mz = 0, bad_im = 0, kept = 0;
    size_t unmapped_frame = 0;   ///< nativeID carried no parseable "frame=": cannot be calibrated
    /// merge per-frame stats after a parallel batch [perf-load]
    CompactStats& operator+=(const CompactStats& o)
    {
      no_im_array += o.no_im_array; size_mismatch += o.size_mismatch;
      bad_mz += o.bad_mz; bad_im += o.bad_im; kept += o.kept;
      unmapped_frame += o.unmapped_frame;
      return *this;
    }
  };

  /// The unsigned integer after `key` in a native ID; 0 when absent.
  inline uint32_t nidField(const String& nid, const char* key)
  {
    const size_t i = nid.find(key);
    if (i == std::string::npos) return 0;
    uint32_t v = 0; size_t k = i + strlen(key);
    for (; k < nid.size() && isdigit((unsigned char)nid[k]); ++k) v = v * 10 + (uint32_t)(nid[k] - '0');
    return v;
  }
  /// The vendor Frames.Id ("frame=", written by the .d loader and the mzPeak reader; "mzpeak=" is an archive index). 0 = unknown.
  inline uint32_t frameIdOf(const String& nid) { return nidField(nid, "frame="); }
  /// "windowGroup=<WindowGroup>": which acquisition group this frame belongs to. A group holds each
  /// isolation m/z at most once, so (m/z, group) names one ion-mobility slice of the window; the
  /// scan range alone does not (two slices can share ScanNumBegin = 0 and differ in ScanNumEnd).
  /// 0 when absent (mzPeak input), which keeps that path keyed by m/z alone, as before.
  inline uint32_t windowGroupOf(const String& nid) { return nidField(nid, "windowGroup="); }

  /// [tile-2d] One frame's picked peaks written into a caller-provided region sized by the raw count
  /// (a pick never emits more points). No per-frame allocation: vectors allocated on picking threads
  /// and freed on the main thread strand free-list in their arenas.
  size_t compactifyInto(const MSSpectrum& s, CompactStats& st, uint32_t* mzq, float* inten, uint16_t* imq, size_t room)
  {
    if (!s.containsIMData()) { st.no_im_array += s.size(); return 0; }
    const auto& im = s.getFloatDataArrays()[s.getIMData().first];
    const Size n = s.size();
    if (im.size() != n) { st.size_mismatch += n; return 0; }   // no out-of-bounds read
    if (n > room) throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "picked frame larger than its raw region", std::to_string(n));
    size_t k = 0;
    for (Size i = 0; i < n; ++i)
    {
      const double mz = s[i].getMZ();
      const double imv = (double)im[i];
      if (!std::isfinite(mz) || mz <= 0.0 || mz * MZ_Q > 4.2e9) { ++st.bad_mz; continue; }
      if (!std::isfinite(imv) || imv < IM_LO || imv > IM_HI)    { ++st.bad_im; continue; }
      mzq[k] = (uint32_t)llround(mz * MZ_Q);
      inten[k] = s[i].getIntensity();
      imq[k] = (uint16_t)llround((imv - IM_LO) * IM_Q);
      ++k; ++st.kept;
    }
    return k;
  }

  /// [det] Order-insensitive 64-bit digest of a stage's numeric output, for localising
  /// run-to-run nondeterminism: sort the bit patterns, then FNV-1a. Logged under DIASPEXTRACT_DET=1.
  static uint64_t detDigest_(vector<uint64_t> bits)
  {
    std::sort(bits.begin(), bits.end());
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t b : bits) { h ^= b; h *= 1099511628211ULL; }
    return h;
  }
  static void pushBits_(vector<uint64_t>& b, double v) { uint64_t u; std::memcpy(&u, &v, 8); b.push_back(u); }
  static bool detOn_() { static const bool on = std::getenv("DIASPEXTRACT_DET") != nullptr; return on; }
  static uint64_t traceDigest_(const TraceStore& st, const vector<Trace>& ts)
  {
    // im (the fragment sort key and scoring's search axis) and tof (the exported m/z) are hashed too.
    vector<uint64_t> b; b.reserve(ts.size() * 5);
    for (const auto& t : ts)
    {
      pushBits_(b, t.mz); pushBits_(b, rtOf(st, t)); pushBits_(b, t.im); pushBits_(b, intensityOf(st, t));
      b.push_back((uint64_t)t.tof);
    }
    return detDigest_(std::move(b));
  }

  /// A window's peaks in index space, frame-major, appended by the loader one frame at a time and
  /// moved into the detector. Until toTof() runs, `tof` holds the 1e-5 Da m/z quantum (`mzq`), which
  /// the OpenMS-detector path reads as is. See docs/MZ-AXIS-DESIGN.md.
  struct PeakSlab
  {
    vector<uint32_t> frame_off;    ///< size F+1: frame f owns [frame_off[f], frame_off[f+1])
    vector<uint32_t> rt_index;     ///< size F: the frame's index on the GLOBAL RT axis (indexRt)
    vector<TofIdx>   tof;          ///< flight-time bin, ascending within a frame (m/z quanta while mzq)
    vector<double>   b;            ///< size F: the frame's calibration factor (bin -> m/z) (toTof)
    vector<float>    inten;
    vector<uint16_t> imq;          ///< the store's quantised 1/K0, kept for export
    vector<double>   rt;           ///< size F until indexRt: the frame's retention time
    vector<uint32_t> frame_id;     ///< size F until toTof: vendor frame Id, the calibration key
    bool mzq = true;               ///< tof[] holds m/z quanta, not bins, until toTof()
    size_t frames() const { return frame_off.empty() ? 0 : frame_off.size() - 1; }
    size_t peaks()  const { return tof.size(); }
    size_t bytes()  const { return tof.size() * 4 + inten.size() * 4 + imq.size() * 2; }
    /// Append one picked frame. Capacity is projected from the running peaks-per-frame x
    /// `expected_frames` (the MS1 frame count: one MS2 frame per window per cycle, and MS1 arrives
    /// first on both readers), over-reserved 2.5x: untouched reserve costs no RSS, a reallocation
    /// copies the window. Never shrunk.
    void append(double frame_rt, uint32_t fid, const uint32_t* mzq_, const float* inten_, const uint16_t* imq_, size_t n,
                size_t expected_frames)
    {
      if (frame_off.empty()) frame_off.push_back(0);
      const size_t F = frames();
      auto grow = [&](auto& v, size_t add) {
        if (v.size() + add <= v.capacity()) return;
        const double per_frame = F ? (double)v.size() / (double)F : (double)add;
        const size_t want = (size_t)(per_frame * (double)std::max(expected_frames, F + 1) * 2.5) + add;
        v.reserve(std::max(v.size() + add, want));
      };
      if (n)   // an EMPTY frame (a padded table frame) carries no peak pointers
      {
        grow(tof, n); grow(inten, n); grow(imq, n);
        tof.insert(tof.end(), mzq_, mzq_ + n);
        inten.insert(inten.end(), inten_, inten_ + n);
        imq.insert(imq.end(), imq_, imq_ + n);
      }
      rt.push_back(frame_rt); frame_id.push_back(fid);
      frame_off.push_back((uint32_t)tof.size());
    }
    size_t reserved() const { return (tof.capacity() - tof.size()) * 4 + (inten.capacity() - inten.size()) * 4 + (imq.capacity() - imq.size()) * 2; }
  };

  /// [cell] A window's FROZEN frame table, in RT order: frame times, vendor ids, RT-axis indices,
  /// calibration factors, and the MS1-frame -> window-frame map the Pearson uses. Store and scoring
  /// support come from here in every mode, so a tile holding only some frames scores exactly as the
  /// whole run does. Filled from the reader's tables (.d streaming source) or from the resident slabs.
  struct FrozenWin
  {
    vector<double>   frame_rt;
    vector<uint32_t> frame_id;
    vector<uint32_t> rt_index;
    vector<double>   b;
    vector<int>      nearest_local;   ///< MS1-store frame -> nearest window frame within delta_rt, or -1
    size_t frames() const { return frame_rt.size(); }
  };

  /// The MS1-frame -> window-frame map over ALL window frames (every frame carries fragment
  /// support: G == F, measured on 28/28 TNBC and 24/24 D windows), ties to the LATER frame.
  inline vector<int> nearestLocal(const vector<double>& win_rt, const vector<double>& ms1_rt, double delta_rt)
  {
    vector<int> nl(ms1_rt.size(), -1);
    size_t hi = 0;
    for (size_t m = 0; m < ms1_rt.size(); ++m)
    {
      const double t = ms1_rt[m];
      while (hi < win_rt.size() && win_rt[hi] < t) ++hi;
      double bd = delta_rt; int gi = -1;
      if (hi > 0)             { const double d = fabs(win_rt[hi - 1] - t); if (d <= bd) { bd = d; gi = (int)(hi - 1); } }
      if (hi < win_rt.size()) { const double d = fabs(win_rt[hi] - t);     if (d <= bd) { bd = d; gi = (int)hi; } }
      nl[m] = gi;
    }
    return nl;
  }

  /// Ordered, incremental digest of a window's slab: frame table (rt index, factor bits, count) then
  /// every (tof, intensity bits, imq) tuple in frame order. Unlike detDigest_ it is NOT sorted, so a
  /// reordered equal-intensity peak or a tuple swapped between frames changes it -- which is what a
  /// representation change of the store must not do. [det]
  static uint64_t slabDigest_(const PeakSlab& sl, uint64_t seed = 1469598103934665603ULL)
  {
    uint64_t h = seed;
    auto mix = [&](uint64_t v) { h ^= v; h *= 1099511628211ULL; };
    for (size_t f = 0; f < sl.frames(); ++f)
    {
      uint64_t bb; std::memcpy(&bb, &sl.b[f], 8);
      mix(sl.rt_index[f]); mix(bb); mix(sl.frame_off[f + 1] - sl.frame_off[f]);
      for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
      { uint32_t ib; std::memcpy(&ib, &sl.inten[k], 4); mix(sl.tof[k]); mix(ib); mix(sl.imq[k]); }
    }
    return h;
  }

  /// Frame times -> indices on the GLOBAL RT axis, once the axis exists; the frames must be in RT
  /// order (the extension loop steps frames by index and takes the RT span as last minus first --
  /// the OpenMS path sorts its spectra for the same reason). The loader delivers them ordered; this
  /// checks rather than assumes.
  void indexRt(PeakSlab& sl)
  {
    const size_t F = sl.frames();
    sl.rt_index.resize(F);
    for (size_t i = 0; i < F; ++i) sl.rt_index[i] = rtIndex(sl.rt[i]);
    for (size_t i = 1; i < F; ++i)
      if (sl.rt_index[i] < sl.rt_index[i - 1])
        throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "window frames are not in retention-time order at frame", std::to_string(i));
    vector<double>().swap(sl.rt);
  }

  /// [cell] One RT cell of a window's slab -- frames [f0, f1) -- as a slab of its own, so the
  /// detector's seeds, `visited`, frame_live and tof index are the cell's alone. Traces found on
  /// it carry cell-local frame indices; the caller rebases them by f0 into the window store.
  PeakSlab sliceFrames(const PeakSlab& sl, size_t f0, size_t f1)
  {
    PeakSlab s;
    const uint32_t p0 = sl.frame_off[f0], p1 = sl.frame_off[f1];
    s.frame_off.assign(sl.frame_off.begin() + f0, sl.frame_off.begin() + f1 + 1);
    for (auto& o : s.frame_off) o -= p0;
    s.rt_index.assign(sl.rt_index.begin() + f0, sl.rt_index.begin() + f1);
    s.b.assign(sl.b.begin() + f0, sl.b.begin() + f1);
    s.tof.assign(sl.tof.begin() + p0, sl.tof.begin() + p1);
    s.inten.assign(sl.inten.begin() + p0, sl.inten.begin() + p1);
    s.imq.assign(sl.imq.begin() + p0, sl.imq.begin() + p1);
    s.mzq = sl.mzq;
    return s;
  }

  /// m/z quanta -> flight-time bins, in place. The store keeps m/z at a 1e-5 Da quantum, ~400x finer
  /// than one TOF bin at m/z 600, so inverting the calibration recovers the ORIGINAL bin exactly --
  /// no approximation is introduced here, and the loader does not have to know the axis.
  ///
  /// [merge] The IM picker's 1/K0 box (gate:delta_im) cuts an ion whose mobility peak is wider into several centroids per
  /// frame that land in ONE bin; the tracer takes one peak per frame, so each slice became its own trace and a spectrum
  /// could carry the same m/z twice. Here, where the bin first exists, a
  /// frame's peaks that share a bin and chain within trace:im_slice_merge in 1/K0 become one peak: intensity summed,
  /// 1/K0 intensity-weighted. Frame-local, so thread and tile invariant; peaks stay ascending by bin within a frame.
  void toTof(PeakSlab& sl)
  {
    if (!sl.mzq) return;
    const long long led_b0 = (long long)sl.bytes();         // [ledger] the load charged the pre-merge bytes
    const TofAxis& ax = tofAxis();
    const size_t F = sl.frames();
    sl.b.resize(F);
    // single linkage along 1/K0 with a link of im_slice_merge_, the whole group within 2x that (one ion's 2-3 slices;
    // measured slice spacing median 0.008, q90 0.014, max 0.019 at the 0.01 picker box)
    const long tol_q = im_slice_merge_ >= 0.0 ? (long)llround(im_slice_merge_ * IM_Q) : -1, span_q = 2 * tol_q;
    // [serial] Frames convert and merge independently, so they run as tasks over frame chunks, each frame compacted
    // within its own range (its cursor starts at its first peak and never passes the read position, as the global
    // cursor never did); the frames are then closed up in order -- the one-cursor loop's layout, byte for byte.
    const long nt = parTasks_(sl.tof.size(), 1 << 17);
    vector<uint32_t> nkept(F);                               // the frame's peaks after the merge
    vector<long long> mc(nt, 0), fmc(nt, 0);                 // per task: peaks merged away, frames with a merge
    forTasks_(nt, [&](long ch) {
    vector<pair<uint16_t, uint32_t>> run;                    // (imq, slab index) of one same-bin run
    vector<std::tuple<uint32_t, float, uint16_t>> cl;        // (first slab index, intensity, imq) of the run's groups
    long long merged = 0, frames_merged = 0;
    for (size_t i = (size_t)ch * F / (size_t)nt; i < (size_t)(ch + 1) * F / (size_t)nt; ++i)
    {
      const double b = ax.factor(sl.frame_id[i]);
      sl.b[i] = b;
      const uint32_t lo = sl.frame_off[i], hi = sl.frame_off[i + 1];
      for (uint32_t k = lo; k < hi; ++k)
        sl.tof[k] = ax.tofOf((double)sl.tof[k] / MZ_Q, b);
      uint32_t w = lo;
      bool any = false;
      for (uint32_t k = lo; k < hi;)
      {
        uint32_t e = k + 1;
        while (e < hi && sl.tof[e] == sl.tof[k]) ++e;
        if (e - k == 1 || tol_q < 0)
        {
          for (uint32_t j = k; j < e; ++j) { sl.tof[w] = sl.tof[j]; sl.inten[w] = sl.inten[j]; sl.imq[w] = sl.imq[j]; ++w; }
          k = e; continue;
        }
        run.clear();
        for (uint32_t j = k; j < e; ++j) run.emplace_back(sl.imq[j], j);
        std::sort(run.begin(), run.end());                   // by 1/K0, then slab order: deterministic
        cl.clear();                                          // every read of the run happens before any write
        for (size_t r = 0; r < run.size();)
        {
          size_t q = r + 1;
          while (q < run.size() && (long)run[q].first - (long)run[q - 1].first <= tol_q
                 && (long)run[q].first - (long)run[r].first <= span_q) ++q;
          if (q - r == 1) { cl.emplace_back(run[r].second, sl.inten[run[r].second], run[r].first); r = q; continue; }   // untouched
          double isum = 0.0, wim = 0.0;
          uint32_t first = run[r].second;
          for (size_t u = r; u < q; ++u)
          {
            first = std::min(first, run[u].second);
            const float v = sl.inten[run[u].second];
            if (v > 0.0f && std::isfinite(v)) { isum += v; wim += (double)v * run[u].first; }
          }
          cl.emplace_back(first, (float)isum, isum > 0.0 ? (uint16_t)std::min(65535.0, std::max(0.0, std::round(wim / isum))) : run[r].first);
          merged += (long long)(q - r - 1); any = true;
          r = q;
        }
        // groups written at their first member's slab position: with nothing merged the slab is byte-identical,
        // so the tracer's same-bin tie-break (slab order) is unchanged wherever no merge happened
        std::sort(cl.begin(), cl.end(), [](const auto& x, const auto& y) { return std::get<0>(x) < std::get<0>(y); });
        const TofIdx t = sl.tof[k];
        for (const auto& c : cl) { sl.tof[w] = t; sl.inten[w] = std::get<1>(c); sl.imq[w] = std::get<2>(c); ++w; }
        k = e;
      }
      nkept[i] = w - lo;
      frames_merged += any;
    }
    mc[ch] = merged; fmc[ch] = frames_merged;
    });
    long long merged = 0, frames_merged = 0;
    for (long c = 0; c < nt; ++c) { merged += mc[c]; frames_merged += fmc[c]; }
    uint32_t w = 0;
    for (size_t i = 0; i < F; ++i)                           // close up: every destination lies at or before its source
    {
      const uint32_t lo = sl.frame_off[i], n = nkept[i];
      if (w != lo)
      {
        std::copy(sl.tof.begin() + lo, sl.tof.begin() + lo + n, sl.tof.begin() + w);
        std::copy(sl.inten.begin() + lo, sl.inten.begin() + lo + n, sl.inten.begin() + w);
        std::copy(sl.imq.begin() + lo, sl.imq.begin() + lo + n, sl.imq.begin() + w);
      }
      sl.frame_off[i] = w; w += n;
    }
    if (F) sl.frame_off[F] = w;
    sl.tof.resize(w); sl.inten.resize(w); sl.imq.resize(w);
    ledAdd_(LC_SLAB, (long long)sl.bytes() - led_b0);         // [ledger] give back what the merge freed
    mergeCount_()[0] += merged; mergeCount_()[1] += frames_merged;
    vector<uint32_t>().swap(sl.frame_id);
    sl.mzq = false;
  }


  /// Per-cell, read-only preparation shared by every band task of detectTracesInteger_: per-frame
  /// "has an above-noise peak" (OpenMS deletes sub-noise peaks before asking whether a frame is
  /// empty), the intensity-ordered seed list segmented by band, and (banded) the tof index.
  struct TracePrep
  {
    vector<char>     frame_live;  ///< F: has at least one member
    vector<uint32_t> order;       ///< seeds, intensity desc, later-index-first among equals
    vector<size_t>   seg;         ///< nb+1 offsets into `order`: band b owns [seg[b], seg[b+1])
    /// [tofidx] Per-frame direct-address index over the flight-time bins: bucket `(tof - tlo) >> tsh`
    /// of frame f starts at peak `fidx[f * (nbk+1) + bucket]`, so best() narrows its search by arithmetic.
    vector<uint32_t> fidx;
    TofIdx tlo = 0; unsigned tsh = 0; uint32_t nbk = 0;
    inline void range(size_t f, TofIdx lo_t, uint32_t a, uint32_t b, uint32_t& s0, uint32_t& s1) const
    {
      if (fidx.empty()) { s0 = a; s1 = b; return; }
      const uint64_t d = lo_t <= tlo ? 0 : ((uint64_t)(lo_t - tlo) >> tsh);
      const uint32_t k = d >= nbk ? nbk : (uint32_t)d;
      const uint32_t* row = &fidx[f * ((size_t)nbk + 1)];
      s0 = row[k]; s1 = (k + 1 <= nbk) ? row[k + 1] : b;
    }
  };
  /// [tofidx] Bucket every frame's flight-time bins so `best()` can address its scan start instead
  /// of binary-searching for it. ~32 peaks per bucket; the shift is global so bucket boundaries
  /// line up across frames. Exact: the bucket start is the first peak with `tof >= bucket_lo`, so
  /// the sub-range handed to lower_bound provably contains the same position it would have found.
  void buildTofIndex(const PeakSlab& sl, TracePrep& tp, TofIdx tlo, TofIdx thi)
  {
    const size_t F = sl.frames(), P = sl.tof.size();
    if (F == 0 || P == 0 || thi < tlo) return;
    const uint64_t span = (uint64_t)thi - tlo + 1;
    const size_t ppf = std::max<size_t>(1, P / F);
    uint32_t want = 64; while (want < 4096 && (size_t)want * 32 < ppf) want <<= 1;
    unsigned sh = 0; while ((span >> sh) > want) ++sh;
    const uint32_t nbk = (uint32_t)std::min<uint64_t>(want, (span >> sh) + 1);
    // one row of nbk+1 offsets per frame; the extra entry is the row's end
    if ((double)F * (nbk + 1) > 3e8) return;                 // pathological geometry: skip the index
    tp.tlo = tlo; tp.tsh = sh; tp.nbk = nbk;
    tp.fidx.assign(F * ((size_t)nbk + 1), 0);
    #pragma omp taskloop default(shared)
    for (long fi = 0; fi < (long)F; ++fi)
    {
      const size_t f = (size_t)fi;
      const uint32_t a = sl.frame_off[f], e = sl.frame_off[f + 1];
      uint32_t* row = &tp.fidx[f * ((size_t)nbk + 1)];
      uint32_t k = a;
      for (uint32_t bk = 0; bk <= nbk; ++bk)
      {
        const uint64_t blo = (uint64_t)bk << sh;
        while (k < e && (uint64_t)(sl.tof[k] - tlo) < blo) ++k;
        row[bk] = k;
      }
    }
  }

  /// `bnd` are the nb+1 flight-time band edges the detector will be given (bnd[0] = tlo,
  /// bnd[nb] = thi+1); pass {0, 0} for the unbanded path, which wants one globally sorted list.
  /// [seg] The seed list is segmented by band and each segment sorted in parallel. Exact because the
  /// comparator is a total order: filter_b(sort(all)) == sort(filter_b(all)).
  TracePrep prepareTracing(const PeakSlab& sl, double noise, double snr, const vector<TofIdx>& bnd,
                           TofIdx tlo, TofIdx thi)
  {
    TracePrep tp;
    const size_t F = sl.frames(), P = sl.tof.size();
    // Membership IS inten > noise (no per-peak map); seeds are counted first so `order` is sized exactly.
    tp.frame_live.assign(F, 0);
    const double min_apex = snr * noise;
    // OpenMS: stable_sort ascending by intensity, then iterate in reverse -> among equal
    // intensities the LATER peak seeds first.
    auto cmp = [&](uint32_t a, uint32_t b) {
      if (sl.inten[a] != sl.inten[b]) return sl.inten[a] > sl.inten[b];
      return a > b;
    };
    const int nb = (int)bnd.size() - 1;
    const bool banded = nb >= 2;                         // unbanded: one globally sorted list
    const size_t nseg = banded ? (size_t)nb : 1;
    // Band of a peak: the same half-open [bnd[b], bnd[b+1]) partition the detector applies, found
    // by upper_bound so that empty bands (bnd[b] == bnd[b+1], reachable when nb exceeds the tof
    // range) are skipped exactly as the detector's own range test skips them.
    // Clamped at both ends: a bin outside [bnd[0], bnd[nb]) would otherwise index cnt/cur out of
    // bounds. The window's edges are asserted to enclose every cell's bins before this runs (a
    // loud throw); the clamp is the memory-safety floor under that assert.
    auto bandOf = [&](TofIdx t) {
      if (!banded) return (size_t)0;
      const auto it = std::upper_bound(bnd.begin(), bnd.end(), t);
      if (it == bnd.begin()) return (size_t)0;
      return std::min((size_t)(it - bnd.begin() - 1), (size_t)(nb - 1));
    };
    // [serial] Count, then distribute, as tasks over peak-balanced FRAME chunks (frame_live stays chunk-private):
    // chunk c counts its seeds per band, the prefix over (band, chunk) is its write cursor, so every band receives
    // its seeds in ascending peak index -- what the one-pass loops wrote. Chunk c owns peaks [kb(c), kb(c + 1)).
    const long nt = parTasks_(P, 1 << 17);
    vector<size_t> fc(nt + 1, F);
    fc[0] = 0;
    for (long c = 1; c < nt; ++c)
      fc[c] = (size_t)(std::lower_bound(sl.frame_off.begin(), sl.frame_off.begin() + F, (uint32_t)((size_t)c * P / (size_t)nt)) - sl.frame_off.begin());
    auto kb = [&](long c) { return fc[c] >= F ? (uint32_t)P : sl.frame_off[fc[c]]; };
    vector<size_t> cnt((size_t)nt * nseg, 0);
    forTasks_(nt, [&](long c) {
      size_t* cc = &cnt[(size_t)c * nseg];
      for (size_t f = fc[c]; f < fc[c + 1]; ++f)
        for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
          if ((double)sl.inten[k] > noise) { tp.frame_live[f] = 1; if ((double)sl.inten[k] > min_apex) ++cc[bandOf(sl.tof[k])]; }
    });
    tp.seg.assign(nseg + 1, 0);
    vector<size_t> cur((size_t)nt * nseg);
    size_t n_seed = 0;
    for (size_t b = 0; b < nseg; ++b)
    {
      tp.seg[b] = n_seed;
      for (long c = 0; c < nt; ++c) { cur[(size_t)c * nseg + b] = n_seed; n_seed += cnt[(size_t)c * nseg + b]; }
    }
    tp.seg[nseg] = n_seed;
    tp.order.resize(n_seed);
    forTasks_(nt, [&](long c) {
      size_t* cc = &cur[(size_t)c * nseg];
      for (uint32_t k = kb(c), ke = kb(c + 1); k < ke; ++k)
        if ((double)sl.inten[k] > noise && (double)sl.inten[k] > min_apex) tp.order[cc[bandOf(sl.tof[k])]++] = k;
    });
    // [serial] A segment holding many seeds (the dense band) is itself sorted as tasks. The comparator is a TOTAL
    // order (the index breaks every tie), so the sorted segment is unique: parallelSort_ returns std::sort's result.
    auto sortSeg = [&](size_t b) {
      const size_t n = tp.seg[b + 1] - tp.seg[b];
      parallelSort_(tp.order.data() + tp.seg[b], n, cmp, parTasks_(n, 1 << 18));
    };
    if (!banded) { sortSeg(0); return tp; }
    std::exception_ptr sort_err;   // [task-err]
    #pragma omp taskloop grainsize(1) default(shared)
    for (int b = 0; b < nb; ++b)
    { try { sortSeg((size_t)b); } catch (...) { taskErr_(sort_err); } }
    if (sort_err) std::rethrow_exception(sort_err);
    buildTofIndex(sl, tp, tlo, thi);
    return tp;
  }

  /// Mass-trace detection on integer arrays with the semantics of OpenMS MassTraceDetection::run_
  /// (written from that source; a version written from a description lost most peptides):
  ///   * seeds in intensity order, ties later-index-first;
  ///   * the two directions extend interleaved, one frame down then one up, sharing one centroid;
  ///   * accept the CLOSEST above-noise peak (in m/z) within +-3 sd of the running centroid and inside
  ///     the mobility tolerance; if it is taken the frame is a MISS (the next-closest is not tried);
  ///   * sd starts at centroid*ppm and is re-estimated from accepted peaks;
  ///   * a direction stops after more than 5 consecutive misses; an empty frame is not a miss;
  ///   * valid iff RT span >= min length and size / (visited - trailing misses) >= min_sample_rate;
  ///     peaks are marked taken only for a valid trace; the seed counts twice in the centroid.
  /// It does not reproduce the OpenMS detector peak for peak: judge it on identified peptides with both
  /// engines, never on a digest (docs/MZ-AXIS-DESIGN.md). Valley splitting runs afterwards.
  vector<Trace> detectTracesInteger_(const PeakSlab& sl, const TracePrep& tp, TraceStore& store, double mass_ppm,
                                     double im_tol, double noise, double min_len_sec,
                                     TofIdx band_lo = 0, TofIdx band_hi = 0,
                                     int band = -1)
  {
    vector<Trace> out;
    const size_t F = sl.frames(), P = sl.tof.size();
    if (F == 0 || P == 0) return out;
    const bool banded = band_hi > band_lo;
    const TofAxis& ax = tofAxis();
    const double msr = 0.5;   // MassTraceDetection's default min_sample_rate
    const double imq_tol = im_tol * IM_Q;
    const Size max_missed = 5;                                            // trace_termination_outliers

    vector<uint64_t> visited((P + 63) / 64, 0);   // one BIT per peak: 8x less than the char array
    auto seen = [&](size_t k) { return (visited[k >> 6] >> (k & 63)) & 1u; };
    auto mark = [&](size_t k) { visited[k >> 6] |= (uint64_t)1 << (k & 63); };
    auto mzAt = [&](uint32_t k, size_t f) { return ax.mzOf(sl.tof[k], sl.b[f]); };
    auto iwm = [](double v, double w, double& c, double& cnt, double& den) {   // updateIterativeWeightedMean_
      const double ct = 1.0 + (w * v) / cnt, dt = 1.0 + w / den;
      c *= (ct / dt); cnt += w * v; den += w;
    };
    auto sdRobust = [](double mz, double w, double mean, double& sd, double& wsum) {   // updateWeightedSDEstimateRobust
      const double d1 = std::log(wsum) + 2.0 * std::log(sd);
      const double d2 = std::log(w) + 2.0 * std::log(std::abs(mz - mean));
      const double denom = std::sqrt(std::exp(d1) + std::exp(d2));
      const double ws = wsum + w;
      const double t = denom / std::sqrt(ws);
      if (t > std::numeric_limits<double>::epsilon()) sd = t;
      wsum = ws;
    };
    /// findBestPeak_ + isPeakAcceptable_: the closest peak in m/z within +-3 sd of the running
    /// centroid whose mobility is inside the gate; -1 for "no candidate", and `taken` set when the
    /// closest one is already used, which OpenMS counts as a miss rather than trying the next.
    /// The +-3 sd window is converted to a BIN range once per frame step, so the walk is bounded
    /// by two calibration calls and a binary search; only candidates inside that window are
    /// converted back to m/z for the comparison.
    auto best = [&](size_t f, double cmz, double sd, double cim, bool& taken) -> long {
      taken = false;
      const uint32_t a = sl.frame_off[f], b = sl.frame_off[f + 1];
      if (b <= a) return -1;
      const double bf = sl.b[f];
      const double lo_mz = cmz - 3.0 * sd, hi_mz = cmz + 3.0 * sd;
      const TofIdx lo_t = ax.tofOf(lo_mz > 0.0 ? lo_mz : 1e-9, bf);
      const TofIdx hi_t = ax.tofOf(hi_mz, bf);
      uint32_t r0 = a, r1 = b; tp.range(f, lo_t, a, b, r0, r1);   // [tofidx] computed, not searched
      const uint32_t s0 = (uint32_t)(lower_bound(sl.tof.begin() + r0, sl.tof.begin() + r1, lo_t) - sl.tof.begin());
      long bk = -1; double bd = std::numeric_limits<double>::infinity();
      for (uint32_t k = s0; k < b && sl.tof[k] <= hi_t; ++k)
      {
        if (!((double)sl.inten[k] > noise)) continue;
        if (std::abs((double)sl.imq[k] - cim) > imq_tol) continue;
        // "Closest" in m/z, as OpenMS does, not in bins: across the centroid the calibration's
        // curvature can flip which is nearer. Only candidates inside the bin window pay the conversion.
        const double d = std::abs(ax.mzOf(sl.tof[k], bf) - cmz);
        if (d < bd) { bd = d; bk = (long)k; }          // strict <: first of equals wins
      }
      if (bk >= 0 && seen((size_t)bk)) { taken = true; return -1; }
      return bk;
    };

    struct Dir { Size missed = 0, scans = 0; bool active = true; };
    vector<pair<uint32_t, uint32_t>> got;   // (frame, peak); marked only if the trace is valid
    // [seg] Only this band's own segment of the seed list; the range test below stays as the
    // invariant it always was (and is what the unbanded/whole-list call still relies on).
    const bool segd = band >= 0 && (size_t)band + 1 < tp.seg.size();
    const size_t s0 = segd ? tp.seg[band] : 0, s1 = segd ? tp.seg[band + 1] : tp.order.size();
    for (size_t si = s0; si < s1; ++si)
    {
      const uint32_t seed = tp.order[si];
      if (banded && (sl.tof[seed] < band_lo || sl.tof[seed] >= band_hi)) continue;   // not this band's core
      if (seen(seed)) continue;
      const size_t f0 = (uint32_t)(std::upper_bound(sl.frame_off.begin(), sl.frame_off.end(), seed) - sl.frame_off.begin() - 1);
      const double w0 = sl.inten[seed], m0 = mzAt(seed, f0);
      double cmz = m0, cnt = w0 * m0, den = w0;
      iwm(m0, w0, cmz, cnt, den);                             // seed counted twice, as run_ does
      double cim = sl.imq[seed], cnt_im = w0 * cim, den_im = w0;
      iwm((double)sl.imq[seed], w0, cim, cnt_im, den_im);
      double sd = (cmz / 1e6) * mass_ppm, wsum = w0;
      got.clear(); got.emplace_back((uint32_t)f0, seed);
      Dir dn, up;
      long fd = (long)f0, fu = (long)f0;
      while ((fd > 0 && dn.active) || (fu + 1 < (long)F && up.active))
      {
        for (int side = 0; side < 2; ++side)
        {
          Dir& d = side == 0 ? dn : up;
          long& fi = side == 0 ? fd : fu;
          if (!d.active) continue;
          if (side == 0 ? !(fi > 0) : !(fi + 1 < (long)F)) continue;
          const size_t f = (size_t)(side == 0 ? fi - 1 : fi + 1);
          // "empty" AFTER sub-noise deletion, as OpenMS tests it; an empty frame is not a miss.
          if (tp.frame_live[f])
          {
            bool taken = false;
            const long k = best(f, cmz, sd, cim, taken);
            if (k >= 0)
            {
              const double mz = mzAt((uint32_t)k, f), w = sl.inten[k];
              iwm(mz, w, cmz, cnt, den);
              iwm((double)sl.imq[k], w, cim, cnt_im, den_im);
              sdRobust(mz, w, cmz, sd, wsum);   // reestimate_mt_sd: the tool never turns it off
              got.emplace_back((uint32_t)f, (uint32_t)k);
              d.missed = 0;
            }
            else ++d.missed;
          }
          fi += (side == 0 ? -1 : 1);
          ++d.scans;
          if (d.missed > max_missed) d.active = false;
        }
      }
      // isTraceValid_
      uint32_t fmin = got[0].first, fmax = fmin;
      for (const auto& g : got) { fmin = std::min(fmin, g.first); fmax = std::max(fmax, g.first); }
      const double span_s = rtAxis()[sl.rt_index[fmax]] - rtAxis()[sl.rt_index[fmin]];
      if (span_s < min_len_sec) continue;
      const Size total = dn.scans + up.scans + 1;
      const Size adjusted = total > dn.missed + up.missed ? total - dn.missed - up.missed : 1;
      if ((double)got.size() / (double)adjusted < msr) continue;
      for (const auto& g : got) mark(g.second);               // valid: now they are taken

      // the span: (store-local frame, intensity) points, plus the per-frame bin arena the valley
      // splitter needs to give each child ITS measured apex bin
      vector<pair<uint32_t, float>> pts;
      pts.reserve(got.size());
      uint32_t kap = got[0].second, fap = got[0].first;
      for (const auto& g : got) { pts.emplace_back(g.first, sl.inten[g.second]);
                                  if (sl.inten[g.second] > sl.inten[kap]) { kap = g.second; fap = g.first; } }
      Trace t = makeSpan(pts, store);
      const uint32_t t_tof = sl.tof[kap];             // this trace's own bin: the offsets' origin
      store.bins.resize(store.inten.size(), 0);
      for (const auto& g : got)
      {
        const int64_t d = (int64_t)sl.tof[g.second] - (int64_t)t_tof;   // [bin-delta]
        if (d < -32768 || d > 32767) throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        store.bins[t.off + (g.first - t.frame0)] = (int16_t)d;
      }
      t.tof = t_tof; t.bframe = fap;   // got-order first max: the frame the bin was measured on
      t.mz = cmz;                                              // the weighted centroid, as OpenMS reports
      t.rt_at = (int16_t)((int64_t)fap - (int64_t)t.frame0);   // [soa-rt] was rtAxis()[sl.rt_index[fap]]
      t.im = IM_LO + cim / IM_Q;
      out.push_back(std::move(t));
    }
    return out;
  }

  /// The acquisition's cycle time: the smallest POSITIVE gap between consecutive frame times. Every gap
  /// is a multiple of the cycle, so skipped frames cannot inflate it (a median can). It is the valley
  /// splitter's scan time (trace:split_scan_time=frame) and the tile:rt_sec cell unit; 0 when no two
  /// frames differ, and the splitter then uses each trace's own average.
  inline double frameSpacing(const vector<double>& frame_rt)
  {
    double m = 0.0;
    for (size_t i = 1; i < frame_rt.size(); ++i)
    { const double d = frame_rt[i] - frame_rt[i - 1]; if (d > 0.0 && (m == 0.0 || d < m)) m = d; }
    return m;
  }

  /// ElutionPeakDetection as a valley splitter: it splits, never deletes (width_filtering off, no SNR filter).
  inline void setSplitter(ElutionPeakDetection& epd, double chrom_fwhm, double scan_time)
  {
    epd.setLogType(ProgressLogger::NONE);
    Param ep = epd.getParameters();
    ep.setValue("chrom_fwhm", chrom_fwhm);
    ep.setValue("width_filtering", "off");
    ep.setValue("masstrace_snr_filtering", "false");
    ep.setValue("scan_time", scan_time);   // 0 = the trace's own average cycle time
    epd.setParameters(ep);
  }

  /// Integer traces -> MassTrace -> ElutionPeakDetection -> Trace, each child's flight-time bin
  /// re-attached from its apex frame's calibration. Chunked and parallel like the OpenMS path's split;
  /// only the traces become OpenMS objects, never the peak arrays.
  vector<Trace> splitIntegerTraces(vector<Trace>& in, double split_valleys, TraceStore& wst, int nchunk, double scan_time)
  {
    // No `in.size() < 2` guard: a lone parent must split as it would among others (tile invariance).
    if (split_valleys <= 0.0 || in.empty()) return std::move(in);
    const TofAxis& ax = tofAxis();
    // [serial] Every parent splits on its own (below), so neither the chunk count nor the cut points can change the
    // output. A large window (>= 8192 parents per chunk) gets up to one chunk per team thread instead of `nchunk`,
    // and the cuts sit at equal cumulative cost (real points + 8 per trace) rather than equal counts: a window's
    // parents come band by band, most intense (longest) first, so equal counts put the longest traces together.
    const size_t n_in = in.size();
    const int team = omp_in_parallel() ? omp_get_num_threads() : 1;
    nchunk = std::max(nchunk, std::min<int>(team, (int)std::min<size_t>(n_in / 8192, 1 << 20)));
    nchunk = std::max(1, std::min<int>(nchunk, (int)n_in));
    vector<size_t> cut(nchunk + 1, n_in);
    cut[0] = 0;
    { uint64_t tot = 0, acc = 0; int c = 1;
      for (const Trace& t : in) tot += (uint64_t)t.npts + 8;
      for (size_t q = 0; q < n_in && c < nchunk; ++q)
      { acc += (uint64_t)in[q].npts + 8;
        while (c < nchunk && acc * (uint64_t)nchunk >= tot * (uint64_t)c) cut[c++] = q + 1; } }
    vector<vector<Trace>> parts(nchunk);
    vector<TraceStore> stores(nchunk);                         // chunk-private arenas, merged after
    std::exception_ptr task_err;   // [task-err]
    #pragma omp taskloop grainsize(1) default(shared)
    for (int ci = 0; ci < nchunk; ++ci)
    { try {
      const size_t lo = cut[ci], hi = cut[ci + 1];
      TraceStore& cst = stores[ci];
      cst.setFrames(wst.rt_index, &wst.b);
      // One MassTrace in flight at a time. Per-trace detectPeaks equals the batch call: with
      // width_filtering=off and masstrace_snr_filtering=false no cross-trace routine runs, and the
      // patched batch path, nested as it is here, is itself a per-trace loop appending in the same
      // order. epd is configured once per chunk: per trace that would be ~4e4 Param round-trips.
      ElutionPeakDetection epd; setSplitter(epd, split_valleys, scan_time);
      vector<Trace>& out = parts[ci];
      out.reserve((hi - lo) + (hi - lo) / 3);   // ~1.32 children per parent, measured
      vector<Peak2D> pk;                                      // reused: a vector, not a per-point list
      vector<MassTrace> split;                                // reused: detectPeaks clears it itself
      for (size_t q = lo; q < hi; ++q)
      {
        const Trace& t = in[q];
        pk.clear();
        for (size_t k = 0; k < t.span(); ++k)                 // REAL points only: EPD smooths what it is given
        {
          if (!real(wst, t, k)) continue;
          Peak2D pt; pt.setRT(rtAtSpan(wst, t, k)); pt.setIntensity(xv(wst, t, k));
          // [bin-delta] reconstruct the absolute bin; identical to the value this used to store
          const uint32_t bin = wst.bins.empty() ? 0u : (uint32_t)((int64_t)t.tof + wst.bins[t.off + k]);
          pt.setMZ(bin ? ax.mzOf(bin, wst.b[t.frame0 + k]) : t.mz);
          pk.push_back(pt);
        }
        MassTrace mt(pk);
        mt.setCentroidIM(t.im);
        epd.detectPeaks(mt, split);
        // Each child's bin comes from its APEX PEAK m/z (the parent's per-frame bin was written into
        // every Peak2D above) by inverting the calibration -- never from the parent list, whose frame
        // ranges overlap.
        for (const MassTrace& cm : split)
        {
          Trace c = toTrace(cm, cst);
          if (c.np() == 0) continue;
          const uint32_t cf = c.frame0 + c.apex;               // apex frame, window-local
          c.bframe = cf;
          const double cb = wst.b.empty() ? ax.factor(0) : wst.b[cf];
          const double amz = cm.getSize() ? cm[cm.findMaxByIntPeak(false)].getMZ() : c.mz;
          c.tof = ax.tofOf(amz, cb);
          out.push_back(std::move(c));
        }
      }
    } catch (...) { taskErr_(task_err); } }
    if (task_err) std::rethrow_exception(task_err);
    // the parents' spans are dead: rebuild the window arena from the chunk arenas
    vector<Trace>().swap(in);
    vector<float>().swap(wst.inten); vector<int16_t>().swap(wst.bins);
    size_t tot = 0, arena = 0; bool any_bins = false;
    for (int ci = 0; ci < nchunk; ++ci)
    { tot += parts[ci].size(); arena += stores[ci].inten.size(); any_bins |= !stores[ci].bins.empty(); }
    // Reserve the merged arena once (absorb() does not); the chunk arenas stay live until absorbed.
    wst.inten.reserve(arena); if (any_bins) wst.bins.reserve(arena);
    vector<Trace> out; out.reserve(tot);
    vector<std::pair<TraceStore*, vector<Trace>*>> src(nchunk);   // [serial] absorbed in chunk order, copies as tasks
    for (int ci = 0; ci < nchunk; ++ci) src[ci] = {&stores[ci], &parts[ci]};
    absorbAll_(wst, out, src, vector<uint32_t>(nchunk, 0));
    return out;
  }

  /// Rebuild a PeakMap for ONE window from its slab (still in m/z quanta), for MassTraceDetection.
  /// The slab is consumed.
  PeakMap materializeWindow(PeakSlab& sl)
  {
    if (!sl.mzq) throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "materializeWindow needs a slab in m/z quanta", "tof");
    PeakMap m;
    const size_t F = sl.frames();
    m.reserve(F);
    for (size_t f = 0; f < F; ++f)
    {
      MSSpectrum s;
      s.setRT(rtAxis()[sl.rt_index[f]]);
      s.setMSLevel(1); // MassTraceDetection traces MS1-level only
      const size_t n = sl.frame_off[f + 1] - sl.frame_off[f];
      s.reserve(n);
      OpenMS::DataArrays::FloatDataArray ima;
      ima.setName(Constants::UserParam::ION_MOBILITY);
      ima.reserve(n);
      for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
      {
        Peak1D p;
        p.setMZ((double)sl.tof[k] / MZ_Q);
        p.setIntensity(sl.inten[k]);
        s.push_back(p);
        ima.push_back((float)(IM_LO + (double)sl.imq[k] / IM_Q));
      }
      s.getFloatDataArrays().push_back(std::move(ima));
      m.addSpectrum(std::move(s));
    }
    sl = PeakSlab();   // consumed
    return m;
  }

  /// Isolation window identity: (lo*100, hi*100, WindowGroup). The group is part of it because one
  /// diaPASEF scheme (PXD017703 "py3") acquires the SAME m/z window in two groups with shifted,
  /// overlapping mobility slices; keyed by m/z alone the halves would merge into one non-monotone
  /// window. With one slice per m/z the key orders and partitions exactly as (lo, hi) did.
  using WinKey = std::array<int, 3>;
  inline WinKey winKey(double lo, double hi, uint32_t group) { return {(int)llround(lo * 100.0), (int)llround(hi * 100.0), (int)group}; }

  /// A refusal of the input: structural frame tables [frames], a non-finite MS1 m/z [ms1-finite], the dnoise MS1 port
  /// [dnoise]. Rethrown past the stream load's fallback, which would otherwise load the whole run again onto the partial
  /// slabs and MS1 map it leaves behind, only to refuse it there.
  struct InputRefusal : Exception::InvalidValue { using Exception::InvalidValue::InvalidValue; };

  /// The MS1 band partition's halo in units of trace:mass_error_ppm: detectTraces_ widens every band by
  /// this many tolerances on each side, and the load-time MS1 prune is exact only because its hop is
  /// derived from the same number. There is deliberately only this one. [ms1-prune]
  constexpr double kBandHaloPpmMul = 20.0;

  /// MSSpectrum::select right-sizes the peaks but reserves every data array at its OLD size (OpenMS MSSpectrum.cpp:49,
  /// 71, 93): without this a selected spectrum keeps its Ion Mobility array's capacity for every peak it lost, 4 B each
  /// (+2.3 GB resident on dataset D's pruned MS1 map). [ms1-prune] [dnoise]
  inline void shrinkDataArrays(MSSpectrum& s)
  {
    for (auto& a : s.getFloatDataArrays()) a.shrink_to_fit();
    for (auto& a : s.getIntegerDataArrays()) a.shrink_to_fit();
    for (auto& a : s.getStringDataArrays()) a.shrink_to_fit();
  }

  /// [ms1-finite] Refuses an input with a non-finite MS1 m/z. Nothing after the pick is defined on one, prune on or off:
  /// sortByPosition leaves finite m/z out of order around a NaN while isSorted() stays true, a NaN band edge leaves two
  /// band cores owning nothing, and MassTraceDetection's m/z search matches anything.
  inline void requireFiniteMs1Mz(const MSSpectrum& s)
  {
    for (Size i = 0; i < s.size(); ++i)
      if (!std::isfinite(s[i].getMZ()))
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                         "[ms1] MS1 spectrum has a non-finite m/z (" + std::to_string(s[i].getMZ()) + ") at peak " + std::to_string(i)
                         + " of " + std::to_string(s.size()) + " after peak picking; the input is malformed and the run is refused", s.getNativeID());
  }

  /// [ms1-prune] Picked MS1 centroids at or below trace:noise_threshold_int, dropped at load: MassTraceDetection never
  /// traces them (its work map keeps (double)intensity > noise only). Two things still read them, and both stay exact:
  ///  (i)  the band EDGES are quantiles of a stride-97 m/z sample of the whole map, so the sample is taken here, from the
  ///       sorted spectrum BEFORE a peak is removed, and detectTraces_ uses it;
  ///  (ii) a band spectrum EXISTS iff any peak lies in the band's range, and one whose peaks are all sub-threshold is still
  ///       a SCAN to MassTraceDetection (it counts against min_sample_rate). So a spectrum keeps its first and last peak plus
  ///       the fewest sub-threshold witnesses that leave no two consecutive kept peaks more than hop (relative m/z) apart.
  ///       A band range with its halo is at least 2 hop / (1 - hop) of its lower bound wide, so a range that holds any peak
  ///       holds a kept one, with a factor-2 margin over floating-point rounding.
  /// Never removes a spectrum or strips the IM array. Measured: docs/BASELINE.md, "MS1 prune, output-identical".
  struct Ms1Prune
  {
    struct Counts { long long picked = 0, survivors = 0, kept = 0, at_noise = 0, nan = 0, bad_mz = 0; };
    bool on = false;              ///< perf:ms1_prune
    bool no_witness = false;      ///< DIASPEXTRACT_MS1_PRUNE_NO_WITNESS=1: survivors only -- a negative control, NOT exact
    bool no_chain = false;        ///< DIASPEXTRACT_MS1_PRUNE_NO_CHAIN=1: interior chain links dropped, anchors kept -- a negative control, NOT exact
    double noise = 0.0;           ///< trace:noise_threshold_int, the very double MassTraceDetection is given
    double hop = 0.0;             ///< trace:mass_error_ppm x 1e-6 x kBandHaloPpmMul
    std::vector<double> sample;   ///< the stride-97 m/z sample of every spectrum pruned so far, taken unpruned
    size_t frames = 0, sample_n = 0;   ///< spectra pruned; the sample size their picked counts imply
    Counts total;

    void reset() { std::vector<double>().swap(sample); frames = sample_n = 0; total = Counts(); }

    /// One spectrum: sort by m/z (a no-op on a picked one), write the m/z values detectTraces_ would
    /// sample to @p samp, then keep the survivors and the witnesses. The chain runs once to count, and a second
    /// time to collect only when it removes a peak: a spectrum it would not shrink costs no index vector and no select.
    void frame(MSSpectrum& s, double* samp, Counts& c) const
    {
      s.sortByPosition();
      const Size n = s.size();
      for (Size i = 0; i < n; i += 97) *samp++ = s[i].getMZ();
      const Size none = (Size)-1;
      auto reach = [&](Size k, double x) { const double m = s[k].getMZ(); return x <= m + hop * m; };
      bool finite = true;
      auto chain = [&](bool count, auto&& take) {   // take(i) for every kept peak, in ascending order
        Size last = none, cand = none;   // the last kept peak; the farthest sub-threshold peak within its reach
        for (Size i = 0; i < n; ++i)
        {
          const double v = s[i].getIntensity();   // float -> double, as MassTraceDetection reads it
          const bool surv = v > noise;            // its own test; the complement also drops NaN and == noise
          const double x = s[i].getMZ();
          if (count) { c.survivors += surv; c.at_noise += (v == noise); c.nan += (v != v); finite = finite && std::isfinite(x); }
          if (no_witness) { if (surv) take(i); continue; }
          if (last != none && reach(last, x)) { if (surv) { take(i); last = i; cand = none; } else cand = i; continue; }
          if (cand != none)                       // out of reach: the candidate becomes a witness first
          {
            if (!no_chain) take(cand);            // (the NO_CHAIN control walks on from it without keeping it)
            last = cand; cand = none;
            if (reach(last, x)) { if (surv) { take(i); last = i; } else cand = i; continue; }
          }
          take(i); last = i;                      // the first peak, or the far side of a gap nothing bridges
        }
        if (cand != none) take(cand);             // the last peak
      };
      Size kept = 0;
      chain(true, [&kept](Size) { ++kept; });
      c.picked += (long long)n;
      // A non-finite m/z hides from is_sorted, so the band builder's lower_bound over it lands where the spectrum's
      // LENGTH decides -- which the prune changes. Both load paths refuse such a spectrum right after its pick, so it
      // never reaches this line and bad_mz stays 0.
      if (!finite) requireFiniteMs1Mz(s);   // [ms1-finite] unreachable, both load paths refuse first; throws if a future entry point skips them
      c.kept += (long long)kept;
      if (kept == n) return;
      std::vector<Size> keep;
      keep.reserve(kept);
      chain(false, [&keep](Size i) { keep.push_back(i); });
      s.select(keep);        // peaks and every data array in lockstep
      shrinkDataArrays(s);   // select right-sizes only the peaks
    }

    /// Prune every spectrum of @p v that @p take selects, in parallel. A spectrum's sample size does not
    /// depend on the prune, so the offsets are fixed first; the order is irrelevant (detectTraces_ sorts).
    void batch(std::vector<MSSpectrum>& v, bool (*take)(const MSSpectrum&))
    {
      std::vector<char> sel(v.size());
      std::vector<size_t> off(v.size());
      size_t at = sample.size();
      for (size_t i = 0; i < v.size(); ++i) { sel[i] = take(v[i]); off[i] = at; if (sel[i]) at += (v[i].size() + 96) / 97; }
      sample.resize(at);
      std::vector<Counts> cs(v.size());
      std::exception_ptr err;
      #pragma omp parallel for schedule(dynamic, 1)
      for (long i = 0; i < (long)v.size(); ++i)
      {
        if (!sel[(size_t)i]) continue;
        try { frame(v[(size_t)i], sample.data() + off[(size_t)i], cs[(size_t)i]); }
        catch (...)
        {
          #pragma omp critical(ms1_prune_err)
          if (!err) err = std::current_exception();
        }
      }
      if (err) std::rethrow_exception(err);
      for (size_t i = 0; i < v.size(); ++i)
      {
        if (!sel[i]) continue;
        const Counts& c = cs[i];
        ++frames; sample_n += (size_t)(c.picked + 96) / 97;
        total.picked += c.picked; total.survivors += c.survivors; total.kept += c.kept; total.at_noise += c.at_noise; total.nan += c.nan; total.bad_mz += c.bad_mz;
      }
    }
  };

  /// [dnoise] dnoise v0.1.0's default MS1 path (DnoiseMs1.h) on every MS1 frame of a Bruker .d, on the raw
  /// points, before the pick. The loader delivers m/z, float32 1/K0 and a corrected intensity, so the
  /// (TOF bin, scan, raw count) dnoise filters are recovered per point by spx::dnoise::recover and the run
  /// is REFUSED on anything that does not recover exactly -- as it is when the loader did not use the
  /// MzCalibration table model, a spectrum is not a whole raw frame, or an MS1 frame with points does
  /// not reach the filter exactly once. Kept points are selected in place; an emptied frame stays.
  struct DnoiseRun
  {
    struct Counts { long long frames = 0, emptied = 0, raw = 0, streak = 0, halo = 0, kept = 0; double cpu = 0.0; uint32_t fid = 0; };
    /// one per thread for one batch: dnoise's Scratch plus the recovered integers (~45 B per raw point in flight)
    struct Work { spx::dnoise::Scratch w; std::vector<uint32_t> scan, tof, raw; std::vector<uint8_t> keep; std::vector<Size> idx; };
    bool on = false;                          ///< dnoise:ms1 true on a Bruker .d, set up
    spx::dnoise::Params p;
    std::optional<spx::dnoise::Ms1WindowGate> gate;
    diaspextract::TdfDnoiseInputs tdf;
    std::vector<std::vector<float>> im_asc;   ///< per TimsCalibration row: the float32 1/K0 of each scan its MS1 frames have, ascending
    std::vector<char> im_falls;               ///< the row's 1/K0 falls with the scan, so its im_asc is reversed
    Counts total;
    spx::dnoise::FrameCoverage seen;          ///< the frame ids add() merged, each once

    static double threadCpu_() { timespec ts; return clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0 ? (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9 : 0.0; }

    /// Read the tdf, build the gate and the scan tables. Empty = ready; anything else is why the run is refused.
    std::string setup(const std::string& tdf_path)
    {
      std::string why;
      if (!diaspextract::loadTdfDnoise(tdf_path, tdf, why)) return why;
      spx::dnoise::GateBuild g = spx::dnoise::buildMs1WindowGate(tdf.meta, tdf.windows, p);
      if (!g.error.empty()) return "the dnoise MS1 port cannot build its window gate from this file: " + g.error;
      gate = std::move(g.gate);
      for (size_t id = 1; id < tdf.frames.size(); ++id)   // above 100 ms raw counts are approximate (recoverPoint), not refused
        if (tdf.frames[id].msms_type == 0)
          if (std::string r = spx::dnoise::recover::correctionRefusal(id, tdf.frames[id].corr); !r.empty()) return r;
      // one table per TimsCalibration row, over the most scans an MS1 frame using it has (tdf.ms1_max_scans): a row no MS1
      // frame uses stays empty and matches nothing. The gate keeps the run's MAX(NumScans), as dnoise does.
      im_asc.assign(tdf.tims_cal.size(), {});
      im_falls.assign(tdf.tims_cal.size(), 0);
      for (size_t ci = 0; ci < tdf.tims_cal.size(); ++ci)
      {
        bool falls = false;
        const std::string r = spx::dnoise::recover::buildScanTable(tdf.tims_cal[ci], tdf.ms1_max_scans[ci], im_asc[ci], falls);
        if (!r.empty()) return r;
        im_falls[ci] = falls;
      }
      reset();
      on = true;
      return std::string();
    }

    /// The scans each row's table covers, for the setup log: "944", or "944,0" with a row no MS1 frame uses.
    String scanTables() const
    {
      String out;
      for (size_t ci = 0; ci < tdf.ms1_max_scans.size(); ++ci) out += (ci ? "," : "") + String(tdf.ms1_max_scans[ci]);
      return out;
    }

    /// The loader must have calibrated m/z with the model the TOF recovery inverts. Serially, once per batch.
    void checkLoader() const
    {
      if (BrukerTimsFile::lastMzCalibration() != "tdf_table_modeltype1")
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[dnoise] MS1 denoising recovers each point's TOF bin through the MzCalibration "
                          "table model, and the loader calibrated m/z otherwise; -dnoise:ms1 false runs without MS1 denoising", BrukerTimsFile::lastMzCalibration());
    }

    /// One MS1 spectrum as the loader delivered it, filtered in place. @p k is the calling thread's, @p c the frame's.
    void frame(MSSpectrum& s, Work& k, Counts& c) const
    {
      using spx::dnoise::recover::Refusal;
      const double cpu0 = threadCpu_();
      const uint32_t fid = frameIdOf(s.getNativeID());
      const Size n = s.size();
      auto refuse = [&](const String& what) {
        return InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[dnoise] " + what + "; -dnoise:ms1 false runs without MS1 denoising", s.getNativeID()); };
      const bool has_im = s.containsIMData();
      switch (spx::dnoise::recover::frameCheck(tdf.frames, fid, n, has_im, has_im ? s.getFloatDataArrays()[s.getIMData().first].size() : 0))
      {
        case Refusal::none: break;
        case Refusal::frame: throw refuse("an MS1 spectrum is not an MS1 frame of the input's analysis.tdf");
        case Refusal::peaks:
          throw refuse("frame " + String(fid) + " arrived with " + String(n) + " of its " + String(tdf.frames[fid].num_peaks) + " raw points: only a whole raw frame can be denoised (no frame aggregation, no loader centroiding)");
        default: throw refuse("frame " + String(fid) + " carries no 1/K0 per point");   // Refusal::im
      }
      const diaspextract::TdfDnoiseInputs::Frame& fr = tdf.frames[fid];
      ++c.frames; c.raw += (long long)n; c.fid = fid;
      if (n == 0) { c.cpu += threadCpu_() - cpu0; return; }
      const auto& im = s.getFloatDataArrays()[s.getIMData().first];
      const TofAxis& ax = tofAxis();
      const double b = ax.factor(fid), corr = fr.corr;
      const std::vector<float>& asc = im_asc[fr.cal];
      const bool falls = im_falls[fr.cal] != 0;
      k.scan.resize(n); k.tof.resize(n); k.raw.resize(n);
      for (Size i = 0; i < n; ++i)
      {
        const float v = s[i].getIntensity();
        double t = 0.0; int hits = 0;
        switch (spx::dnoise::recover::recoverPoint(ax.cal, b, s[i].getMZ(), im[i], v, corr, asc, falls, fr.num_scans, k.tof[i], k.scan[i], k.raw[i], t, hits))
        {
          case Refusal::none: break;
          case Refusal::tof:
            throw refuse("m/z " + String(s[i].getMZ()) + " in frame " + String(fid) + " is not the MzCalibration table model's value of a TOF bin (" + String(t) + "): the loader calibrated m/z otherwise");
          case Refusal::scan:
            throw refuse("1/K0 " + String(im[i]) + " in frame " + String(fid) + " is not the TimsCalibration model's value of one of its " + String(fr.num_scans) + " scans: the loader converted 1/K0 otherwise");
          case Refusal::intensity_range:
            throw refuse("intensity " + String(v) + " in frame " + String(fid) + " is >= 2^24 = 16777216: the loader's float32 intensity no longer identifies one raw count");
          default:   // Refusal::intensity
            throw refuse("intensity " + String(v) + " in frame " + String(fid) + " is not opentims' correction of any raw count");
        }
      }
      const spx::dnoise::FrameCounts fc = spx::dnoise::denoiseMs1Frame(k.scan.data(), k.tof.data(), k.raw.data(), n, fr.num_scans, p, gate ? &*gate : nullptr, k.w, k.keep);
      c.streak += (long long)fc.after_streak; c.halo += (long long)fc.after_halo; c.kept += (long long)fc.kept; c.emptied += (fc.kept == 0);
      if (fc.kept != n)
      {
        k.idx.clear();
        for (Size i = 0; i < n; ++i) if (k.keep[i]) k.idx.push_back(i);
        s.select(k.idx);        // the peaks and every data array in lockstep
        shrinkDataArrays(s);    // select leaves every data array at its raw capacity
      }
      c.cpu += threadCpu_() - cpu0;
    }

    /// A batch's per-spectrum counts, serially: each denoised frame (a slot with frames 0 is an MS2 spectrum's) marked
    /// once -- a frame delivered twice is refused, it would be filtered and traced twice.
    void add(const std::vector<Counts>& v)
    {
      for (const Counts& c : v)
      {
        if (c.frames == 0) continue;
        if (!seen.mark(c.fid))
          throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[dnoise] MS1 frame " + String(c.fid) + " reached the filter twice: the loader delivered it more than once; -dnoise:ms1 false runs without MS1 denoising", "frame=" + String(c.fid));
        total.frames += c.frames; total.emptied += c.emptied; total.raw += c.raw; total.streak += c.streak; total.halo += c.halo; total.kept += c.kept; total.cpu += c.cpu;
      }
    }

    /// After the load, serially: every MS1 frame with points (Frames.NumPeaks > 0) reached the filter. Counts cannot
    /// tell -- a frame delivered twice and a missing one with equal NumPeaks cancel out -- so this reads the marks.
    void requireComplete() const
    {
      if (const size_t miss = seen.firstMissing(tdf.frames))
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[dnoise] MS1 frame " + String(miss) + " (" + String(tdf.frames[miss].num_peaks)
                          + " raw points in Frames.NumPeaks) never reached the filter: the loader did not deliver every MS1 frame; -dnoise:ms1 false runs without MS1 denoising", "frame=" + String(miss));
    }

    /// The counts and marks back to none: at setup, and when a failed stream load falls back to the resident one, which
    /// delivers every frame again.
    void reset() { total = Counts(); seen.reset(); }

    /// The applied parameters and the run's gate size, for the log and spx:dnoise_ms1_params.
    String params() const
    {
      return "mz_half_width=" + String(p.mz_half_width) + " min_feature_length=" + String(p.min_feature_length) + " max_internal_gap=" + String(p.max_internal_gap)
             + " iterations=" + String(p.iterations) + " min_window_intensity=" + String(p.min_window_intensity) + " min_feature_intensity=" + String(p.min_feature_intensity)
             + " halo=" + (p.halo ? "true" : "false") + " halo_peak_fraction=" + String(p.halo_peak_fraction) + " halo_mz_idx_half_width=" + String(p.halo_mz_idx_half_width)
             + " halo_scan_half_width=" + String(p.halo_scan_half_width) + " dia_ms1_window=" + (p.dia_ms1_window ? "true" : "false")
             + " dia_ms1_mz_pad=" + String(p.dia_ms1_mz_pad) + " dia_ms1_im_pad=" + String(p.dia_ms1_im_pad) + " gate_boxes=" + String(gate ? gate->boxes().size() : (size_t)0);
    }

    String summary() const
    {
      char b[320];
      snprintf(b, sizeof b, "[dnoise] %s MS1: %lld frames (%lld emptied), %lld raw points -> streak %lld -> halo %lld -> window gate: kept %lld MS1 points (%.2f%%); filter CPU %.1f s",
               spx::dnoise::kPortedFrom, total.frames, total.emptied, total.raw, total.streak, total.halo, total.kept,
               total.raw > 0 ? 100.0 * (double)total.kept / (double)total.raw : 0.0, total.cpu);
      return String(b);
    }
  };

  /// [stream] Pick-and-compact consumer: each frame is picked and compacted on arrival, so the raw run
  /// is never resident. Windows are keyed from the spectrum's own precursor, not swath_nr.
  class PickCompactConsumer : public FullSwathFileConsumer
  {
  public:
    PickCompactConsumer(PeakPickerIM& picker, map<WinKey, PeakSlab>& windows,
                        PeakMap& ms1, CompactStats& stats)
      : picker_(picker), windows_(windows), ms1_(ms1), stats_(stats) {}
    /// [outer] the MS1 worker ends with the consumer, also when an exception leaves the load early
    ~PickCompactConsumer() { stopMs1Worker_(); }

    size_t frames_seen = 0;
    long long slab_charged_ = 0;   ///< [ledger] LC_SLAB bytes this consumer's MS2 batches charged (the pipelined tile read replaces them)
    String ms1_pick_mode_;   ///< [outer] how the MS1 batches were picked, for the log (set by finish())

    /// [frames] A window's FROZEN frame table (rt, vendor id) with the lockstep cursor: every table
    /// frame the reader skipped (no raw peaks, an empty mobility slice) is appended EMPTY with the
    /// table's rt and id, so the slab's frame sequence IS the table's in every mode (an empty frame
    /// is a scan, not a miss).
    struct FrameTab
    {
      std::vector<double> rt; std::vector<uint32_t> id;
      size_t cur = 0, delivered = 0, padded = 0;
      /// pad up to the delivered frame (its rt AND vendor id equal the table row's -- an rt alone is
      /// ambiguous when two rows share one, e.g. an mzPeak archive without retention times); false
      /// = not in the table
      bool step(PeakSlab& sl, double t, uint32_t fid)
      {
        while (cur < rt.size() && (rt[cur] != t || id[cur] != fid)) { sl.append(rt[cur], id[cur], nullptr, nullptr, nullptr, 0, rt.size()); ++cur; ++padded; }
        if (cur == rt.size()) return false;
        ++cur; ++delivered; return true;
      }
      void tail(PeakSlab& sl)
      { for (; cur < rt.size(); ++cur, ++padded) sl.append(rt[cur], id[cur], nullptr, nullptr, nullptr, 0, rt.size()); }
    };
    using FrameTables = std::map<WinKey, FrameTab>;
    /// The reader's frozen tables, fetched at the first MS2 flush: both streaming readers publish
    /// their metadata before their first MS2 frame (the .d loader exports its frame table and
    /// windows before its MS2 loop, the mzPeak reader sweeps the archive before decoding). Unset
    /// (an mzML) = no tables, no padding: the slab is its own table.
    std::function<FrameTables()> tablefn;
    const FrameTables& tables() const { return tables_; }
    /// [ms1-prune] pass 1's run-level MS1 prune and its edge sample; null = every picked peak is kept
    Ms1Prune* ms1_prune = nullptr;
    /// [dnoise] the run's MS1 denoising of the raw frames before their pick; null or off = none
    DnoiseRun* dnoise = nullptr;

    /// [perf-load] The readers hand frames over serially: pick+compact a buffered BATCH in parallel, then
    /// append strictly by buffer index (arrival order).
    void flush_()
    {
      if (!tables_init_) { tables_init_ = true; if (tablefn) tables_ = tablefn(); }
      const int n = (int)buf_.size();
      if (n == 0) return;
      waitMs1Worker_();   // [review] DIASPEXTRACT_MS1_ASYNC: the worker appends to ms1_ (read below) and adds to flush_wall_()
      const double _tf0 = phase_clock_();
      // [ledger] the batch's RAW frames, 20 B per raw peak (Peak1D 16 + IM float 4), freed one by one in the pick loop
      long long raw_b = 0;
      for (const auto& sp : buf_) raw_b += (long long)sp.size() * 20;
      LedCharge _raw(LC_SCRATCH, raw_b);
      // per-frame regions of one reusable scratch, sized by the RAW peak counts
      scratch_off_.resize((size_t)n + 1);
      { size_t total = 0; for (int i = 0; i < n; ++i) { scratch_off_[(size_t)i] = total; total += buf_[(size_t)i].size(); } scratch_off_[(size_t)n] = total;
        if (scratch_mzq_.size() < total) { scratch_mzq_.resize(total); scratch_inten_.resize(total); scratch_imq_.resize(total); } }
      static std::atomic<long long> scratch_charged{0};   // [ledger] the reusable compaction scratch, once
      { const long long want = (long long)(scratch_mzq_.capacity() * 4 + scratch_inten_.capacity() * 4 + scratch_imq_.capacity() * 2);
        ledAdd_(LC_SCRATCH, want - scratch_charged.exchange(want)); }
      std::vector<size_t> kept((size_t)n);
      std::vector<double> frt((size_t)n);
      std::vector<uint32_t> fid((size_t)n);
      std::vector<CompactStats> st((size_t)n);          // per-frame stats, merged serially below
      PeakPickerIM picker = picker_;                    // Param is a deep-copy value type
      std::exception_ptr pick_err;                      // an exception escaping an OMP region aborts
      auto pickOne = [&](const int i, PeakPickerIM& pk)
      {
        try
        {
          auto& s = buf_[(size_t)i];
          ensureIMArrayName(s);
          if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) pk.pickIMCluster(s);
          ensureIMArrayName(s);
          frt[(size_t)i] = s.getRT();
          fid[(size_t)i] = frameIdOf(s.getNativeID());
          if (fid[(size_t)i] == 0) ++st[(size_t)i].unmapped_frame;
          const size_t o = scratch_off_[(size_t)i];
          kept[(size_t)i] = compactifyInto(s, st[(size_t)i], &scratch_mzq_[o], &scratch_inten_[o], &scratch_imq_[o], scratch_off_[(size_t)i + 1] - o);
          s.clear(true);                                // raw frame dies here, inside the loop
        }
        catch (...)
        {
          #pragma omp critical(pick_err)
          if (!pick_err) pick_err = std::current_exception();
        }
      };
      if (!omp_in_parallel())
      {
        // chunk 1: with chunk 8 a 64-frame batch made only 8 chunks = 8 workers
        #pragma omp parallel for firstprivate(picker) schedule(dynamic, 1)
        for (int i = 0; i < n; ++i) pickOne(i, picker);
      }
      else
      {
        // [pipe] the pipelined tile reader, on the window team's driver: one picker copy per puller slot, made
        // by the slot's first frame (a slot has one puller at a time)
        std::vector<std::unique_ptr<PeakPickerIM>> pk((size_t)std::max(1, omp_get_num_threads()));
        teamFor_(n, 1, [&](long i, int slot) {
          std::unique_ptr<PeakPickerIM>& p = pk[(size_t)slot];
          if (!p) p.reset(new PeakPickerIM(picker));
          pickOne((int)i, *p); });
      }
      if (pick_err) std::rethrow_exception(pick_err);
      // Append in arrival order per window, the windows of this batch in parallel: a window's
      // frames stay in order (one thread walks them), distinct windows never touch each other's
      // slab, and the copy of ~1 GB per batch into cold pages no longer runs on one thread.
      struct Batch { PeakSlab* slab; FrameTab* tab; std::vector<int> idx; };
      std::vector<Batch> per_win;
      { std::map<WinKey, size_t> at;
        for (int i = 0; i < n; ++i)
        {
          stats_ += st[(size_t)i];
          const WinKey& k = key_[(size_t)i];
          auto it = at.find(k);
          if (it == at.end())
          { it = at.emplace(k, per_win.size()).first;
            auto tt = tables_.find(k);
            if (!tables_.empty() && tt == tables_.end())
              throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[frames] a delivered MS2 window is not in the reader's frame table",
                                    String(k[0] / 100.0) + "-" + String(k[1] / 100.0) + " group " + String(k[2]));
            per_win.push_back(Batch{&windows_[k], tt == tables_.end() ? nullptr : &tt->second, {}}); }
          per_win[it->second].idx.push_back(i);
        } }
      const size_t n_ms1 = ms1_.size();
      std::atomic<long long> slab_b{0};   // [ledger] what this batch actually touched
      std::exception_ptr pad_err;
      auto appendWin = [&](const long w)
      {
        Batch& bw = per_win[(size_t)w];
        for (int i : bw.idx)
        {
          const size_t exp_f = bw.tab ? bw.tab->rt.size() : n_ms1;   // the window's frame count: exact when a table exists
          if (bw.tab && !bw.tab->step(*bw.slab, frt[(size_t)i], fid[(size_t)i]))   // [frames] lockstep with the frozen table
          {
            #pragma omp critical(pad_err)
            if (!pad_err) pad_err = std::make_exception_ptr(InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "[frames] a delivered MS2 frame (rt, frame id) is not in its window's frozen frame table (window "
                  + String(key_[(size_t)i][0] / 100.0) + "-" + String(key_[(size_t)i][1] / 100.0) + ", frame " + String(fid[(size_t)i]) + ")", String(frt[(size_t)i], 6)));
            break;
          }
          const size_t o = scratch_off_[(size_t)i];
          bw.slab->append(frt[(size_t)i], fid[(size_t)i], &scratch_mzq_[o], &scratch_inten_[o], &scratch_imq_[o], kept[(size_t)i], exp_f);
          slab_b += (long long)kept[(size_t)i] * 10 + 24;   // [ledger] touched bytes: 10 B/peak + the frame's table row
        }
      };
      if (!omp_in_parallel())
      {
        #pragma omp parallel for schedule(dynamic, 1)
        for (long w = 0; w < (long)per_win.size(); ++w) appendWin(w);
      }
      else teamFor_((long)per_win.size(), 1, [&](long w, int) { appendWin(w); });   // [pipe] as the pick above
      ledAdd_(LC_SLAB, slab_b); slab_charged_ += slab_b;
      if (pad_err) std::rethrow_exception(pad_err);
      flush_wall_() += phase_clock_() - _tf0;
      buf_.clear(); key_.clear();
      loadTrim_();
    }

  protected:
    void consumeMS1Spectrum_(MapType::SpectrumType& s) override
    {
      // MS1 stays a PeakMap (MassTraceDetection needs it whole); buffered and picked in parallel like MS2, order kept by index.
      // [dnoise] the loader chose its m/z model before its first frame: refuse at that frame, not after a whole batch
      // (flushMS1_ checks every batch again and stays the authority)
      if (dnoise && dnoise->on && !dn_loader_checked_) { dnoise->checkLoader(); dn_loader_checked_ = true; }
      ms1_buf_.push_back(std::move(s));
      ++frames_seen;
      if (ms1_buf_.size() >= ms1BatchNow_()) flushMS1Async_();
    }

    /// Picks the buffered MS1 frames now, after the worker's batch (if any): the order of the MS1 map is the arrival order.
    void flushMS1_()
    {
      waitMs1Worker_();
      pickMS1Batch_(ms1_buf_);
    }

    /// [outer] Hands the buffered MS1 frames to the worker thread and returns to the reader, which decodes the next batch
    /// meanwhile. The previous batch is picked and appended first, so batches still reach the MS1 map in arrival order,
    /// and every frame is processed exactly as by flushMS1_: batch composition does not change the output (kBatch_).
    void flushMS1Async_()
    {
      if (!ms1Async_()) { flushMS1_(); return; }
      waitMs1Worker_();
      ms1_inflight_.swap(ms1_buf_);   // ms1_buf_ takes the picked (empty) buffer back
      ms1_buf_.clear();
      if (!w_th_.joinable())
      {
        const int nthr = std::max(1, omp_get_max_threads());   // a new thread starts from the global default team size
        w_th_ = std::thread([this, nthr] {
          omp_set_num_threads(nthr);
          std::unique_lock<std::mutex> lk(w_mu_);
          for (;;)
          {
            w_cv_.wait(lk, [this] { return w_job_ || w_quit_; });
            if (!w_job_) return;   // quit, nothing pending
            lk.unlock();
            std::exception_ptr err;
            try { pickMS1Batch_(ms1_inflight_); } catch (...) { err = std::current_exception(); }
            lk.lock();
            w_err_ = err; w_job_ = false;
            w_cv_.notify_all();
          }
        });
      }
      { std::lock_guard<std::mutex> g(w_mu_); w_job_ = true; }
      w_cv_.notify_all();
    }

    /// [outer] Waits for the worker's batch and rethrows its error, if any.
    void waitMs1Worker_()
    {
      if (!w_th_.joinable()) return;
      std::unique_lock<std::mutex> lk(w_mu_);
      w_cv_.wait(lk, [this] { return !w_job_; });
      if (w_err_) { std::exception_ptr e = w_err_; w_err_ = nullptr; lk.unlock(); std::rethrow_exception(e); }
    }

    /// [outer] Lets the worker finish its batch and ends it (finish() and the destructor; an unreported error is dropped
    /// here only when an exception is already unwinding the load).
    void stopMs1Worker_()
    {
      if (!w_th_.joinable()) return;
      { std::lock_guard<std::mutex> g(w_mu_); w_quit_ = true; }
      w_cv_.notify_all();
      w_th_.join();
      w_quit_ = false;
    }

    void pickMS1Batch_(std::vector<MapType::SpectrumType>& buf)   // one MS1 batch, picked and appended to the MS1 map, then emptied
    {
      const int n = (int)buf.size();
      if (n == 0) return;
      const double _tf0 = phase_clock_();
      // [ledger] the buffered RAW MS1 frames, 20 B per raw peak, held until picked: pass 1's largest transient.
      long long raw_b = 0;
      for (const auto& sp : buf) raw_b += (long long)sp.size() * 20;
      LedCharge _raw(LC_SCRATCH, raw_b);
      PeakPickerIM picker = picker_;
      std::exception_ptr pick_err;
      // [dnoise] per thread for this batch only (freed with it), per frame counts merged below
      const bool dn = dnoise && dnoise->on;
      if (dn) dnoise->checkLoader();
      std::vector<DnoiseRun::Work> dn_work(dn ? (size_t)std::max(1, omp_get_max_threads()) : 0);
      std::vector<DnoiseRun::Counts> dn_cnt(dn ? (size_t)n : 0);
      #pragma omp parallel for firstprivate(picker) schedule(dynamic, 1)
      for (int i = 0; i < n; ++i)
      {
        try
        {
          auto& s = buf[(size_t)i];
          ensureIMArrayName(s);
          if (dn) dnoise->frame(s, dn_work.at((size_t)omp_get_thread_num()), dn_cnt[(size_t)i]);   // [dnoise] the raw points, before the pick
          if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
          ensureIMArrayName(s);
          requireFiniteMs1Mz(s);   // [ms1-finite] before the prune, prune on or off
        }
        catch (...)
        {
          #pragma omp critical(pick_err_ms1)
          if (!pick_err) pick_err = std::current_exception();
        }
      }
      if (pick_err) std::rethrow_exception(pick_err);
      if (dn) dnoise->add(dn_cnt);
      if (ms1HistOn_())
        for (const auto& sp : buf)
        {
          long long c[5] = {0, 0, 0, 0, 0};
          for (const auto& pk : sp) { const float v = pk.getIntensity(); ++c[v < 100 ? 0 : v < 300 ? 1 : v < 1000 ? 2 : v < 10000 ? 3 : 4]; }
          for (int k = 0; k < 5; ++k) ms1Hist_()[(size_t)k].fetch_add(c[k], std::memory_order_relaxed);
        }
      // [ms1-prune] after [ms1-hist] (which counts what was PICKED), before the charge below (which is then
      // the kept bytes)
      if (ms1_prune && ms1_prune->on) ms1_prune->batch(buf, [](const MSSpectrum&) { return true; });
      { long long mb = 0; for (const auto& s : buf) mb += (long long)s.size() * 20; ledAdd_(LC_MS1MAP, mb); }
      for (auto& s : buf) ms1_.addSpectrum(std::move(s));
      buf.clear();
      flush_wall_() += phase_clock_() - _tf0;
      loadTrim_();
    }

    void consumeSwathSpectrum_(MapType::SpectrumType& s, size_t /*swath_nr*/) override
    {
      const auto& prec = s.getPrecursors();
      if (prec.empty()) return;                      // no isolation window -> not routable
      const double c = prec[0].getMZ();
      const double lo = c - prec[0].getIsolationWindowLowerOffset();
      const double hi = c + prec[0].getIsolationWindowUpperOffset();
      flushMS1_();                                   // MS1 arrives first; pick it before any MS2 work
      key_.push_back(winKey(lo, hi, windowGroupOf(s.getNativeID())));
      buf_.push_back(std::move(s));                  // pick+compact deferred to flush_()
      ++frames_seen;
      if (buf_.size() >= kBatch_()) flush_();
    }

    void ensureMapsAreFilled_() override { finish(); }

  public:
    /// Flush the buffered tails. Callers MUST call this after the load: FullSwathFileConsumer reaches
    /// ensureMapsAreFilled_ only via retrieveSwathMaps, which this consumer never uses, so without it the
    /// last `n mod kBatch` spectra are dropped (docs/BASELINE.md, step 0b).
    void finish()
    {
      flushMS1_(); flush_();
      stopMs1Worker_();
      ms1_pick_mode_ = ms1Async_() ? "beside the reader on a worker thread, batches of " + String(ms1BatchNow_())
                                   : "inside the reader's hand-over, batches of " + String(kBatchMS1_());
      // [frames] the tail of every delivered window's table; a table window that delivered nothing
      // gets no slab (it is counted as absent by the caller)
      for (auto& kv : tables_)
      { auto w = windows_.find(kv.first); if (w != windows_.end()) kv.second.tail(w->second); }
    }

  private:
    /// [ledger] DIASPEXTRACT_LOAD_TRIM=N: malloc_trim(0) every Nth flush (0 = never). The load's retention
    /// sits in per-thread arenas, which only an explicit malloc_trim reaches (MALLOC_TRIM_THRESHOLD_ does
    /// not). Priced in docs/BASELINE.md, "The composition on D".
    static int loadTrimEvery_() { static const int v = []{ const char* e = std::getenv("DIASPEXTRACT_LOAD_TRIM");
                                                           return e ? std::atoi(e) : 0; }(); return v; }
    static void loadTrim_()
    {
      const int n = loadTrimEvery_();
      if (n <= 0) return;
      static std::atomic<long long> k{0};
      const bool due = (k.fetch_add(1, std::memory_order_relaxed) % n == n - 1);
#ifdef __GLIBC__
      if (due) malloc_trim(0);
#else
      (void)due;   // every other trim site in this file is glibc-only for the same reason
#endif
    }
    /// Frames per pick batch, 256 by default (~1 GB raw buffer); DIASPEXTRACT_PICK_BATCH overrides.
    /// Batch composition does not change the output.
    static size_t kBatch_() { static const size_t v = []{ const char* e = std::getenv("DIASPEXTRACT_PICK_BATCH");
      const long x = e ? std::atol(e) : 0; return x > 0 ? (size_t)x : (size_t)256; }(); return v; }
    /// [batch] DIASPEXTRACT_PICK_BATCH_MS1: the MS1-only flush batch (unset = kBatch_), so shrinking it for memory does not
    /// re-batch the MS2 tile reads (a smaller MS2 batch only costs CPU). Output-invariant; priced in docs/BASELINE.md "The composition on D".
    static size_t kBatchMS1_() { static const size_t v = []{ const char* e = std::getenv("DIASPEXTRACT_PICK_BATCH_MS1");
      const long x = e ? std::atol(e) : 0; return x > 0 ? (size_t)x : (size_t)0; }();
      return v ? v : kBatch_(); }
    std::vector<MapType::SpectrumType> buf_;
    std::vector<WinKey> key_;
    FrameTables tables_; bool tables_init_ = false;   ///< [frames] the reader's frozen tables + cursors
    std::vector<size_t> scratch_off_;                 ///< per-frame regions of the batch scratch (raw sizes)
    std::vector<uint32_t> scratch_mzq_; std::vector<float> scratch_inten_; std::vector<uint16_t> scratch_imq_;
    PeakPickerIM& picker_;
    std::vector<MapType::SpectrumType> ms1_buf_;   // MS1 frames awaiting the parallel pick
    bool dn_loader_checked_ = false;               ///< [dnoise] checkLoader ran on the first MS1 frame
    /// [outer] The MS1 worker: picks ms1_inflight_ while the reader fills ms1_buf_. w_job_ is set by the reader's thread
    /// and cleared by the worker, both under w_mu_; the worker touches only its batch, the MS1 map, the prune, the dnoise
    /// counts and the ledger (atomic), none of which the reader's thread uses before waitMs1Worker_ has returned.
    std::vector<MapType::SpectrumType> ms1_inflight_;
    std::thread w_th_;
    std::mutex w_mu_;
    std::condition_variable w_cv_;
    bool w_job_ = false, w_quit_ = false;
    std::exception_ptr w_err_;
    /// [outer] DIASPEXTRACT_MS1_ASYNC=1: the MS1 batches are picked on the worker, beside the reader's decode of the next
    /// batch. OFF by default: the worker's team is a second OpenMP team on the same places as the reader's (libgomp binds
    /// a foreign thread's team from place 0), so at 64 threads a 60-min run's LOAD(stream) got 0.3 s slower, and its extra thread
    /// arenas raised a 2-h HeLa run's peak RSS by 6 % for 0.5-0.9 s (timed 2026-10-01). Default: the MS1 batches
    /// are picked inside the reader's hand-over, one full batch at a time (the behaviour before 2026-10-01).
    static bool ms1Async_() { static const bool v = []{ const char* e = std::getenv("DIASPEXTRACT_MS1_ASYNC"); return e && *e && *e != '0'; }(); return v; }
    /// [outer] Beside the reader two batches are in memory (one being picked, one filling), so each is half the
    /// synchronous batch and the raw-frame residency stays what it was.
    static size_t ms1BatchNow_() { return ms1Async_() ? std::max<size_t>(1, kBatchMS1_() / 2) : kBatchMS1_(); }
    map<WinKey, PeakSlab>& windows_;
    PeakMap& ms1_;
    CompactStats& stats_;
  };

  /// What the scorer needs per window: the Pearson denominator G (the frozen table's frame count),
  /// per-fragment mean and 1/norm over G, and the MS1 -> window nearest-frame map.
  struct FragStats
  {
    vector<double> mean, invnorm;   ///< per fragment, over the full grid of G frames
    const vector<int>* nearest_local = nullptr;   ///< MS1-store frame -> nearest window frame within delta_rt, or -1 (the FrozenWin's)
    size_t G = 0;
  };

  /// The Pearson support is the window's FROZEN frame table: G = F (every frame), and the
  /// MS1 -> window map is the table's, so neither depends on which frames are resident.
  FragStats buildFragStats(const vector<Trace>& frags, const TraceStore& wst, const FrozenWin& fw)
  {
    FragStats g;
    g.G = fw.frames();
    const double G = (double)g.G;
    g.nearest_local = &fw.nearest_local;
    g.mean.resize(frags.size()); g.invnorm.resize(frags.size());
    #pragma omp taskloop default(shared)
    for (long ii = 0; ii < (long)frags.size(); ++ii)
    {
      const Trace& f = frags[(size_t)ii];
      double sum = 0, sumsq = 0;
      for (size_t k = 0; k < f.span(); ++k) { const double v = xv(wst, f, k); if (v > 0.0) { sum += v; sumsq += v * v; } }
      g.mean[ii] = sum / G;
      const double var = sumsq - G * g.mean[ii] * g.mean[ii];
      g.invnorm[ii] = var > 0 ? 1.0 / sqrt(var) : 0.0;  // 0 => constant/degenerate -> never correlates
    }
    return g;
  }

  /// Canonical total-order sort (RT, precursor m/z, charge, then the fragment sequence):
  /// thread-count-independent and free of RT-tie ambiguity (RT-only sortSpectra is NOT canonical).
  /// Sorts a compact key array and permutes once: std::sort moves what it is handed, and an
  /// MSSpectrum is ~750 bytes against 32 for a key. The keys start in the spectra's order, so
  /// introsort makes the same decisions as sorting the spectra would, ties included. [par-Crit-2]
  /// Tiles are RT-disjoint and RT is the first key, so per-tile sorted blocks concatenate to the
  /// order one sort of everything gives.
  static void canonicalSort_(vector<MSSpectrum>& v)
  {
    struct SortKey { double rt, mz; int charge; uint32_t n, idx; };
    vector<SortKey> keys(v.size());
    auto key = [&](const long i)
    {
      const MSSpectrum& a = v[i];
      keys[i].rt = a.getRT();
      keys[i].mz = a.getPrecursors().empty() ? 0.0 : a.getPrecursors()[0].getMZ();
      keys[i].charge = a.getPrecursors().empty() ? 0 : a.getPrecursors()[0].getCharge();
      keys[i].n = (uint32_t)a.size();
      keys[i].idx = (uint32_t)i;
    };
    if (!omp_in_parallel())
    {
#pragma omp parallel for schedule(static)
      for (long i = 0; i < (long)v.size(); ++i) key(i);
    }
    else teamFor_((long)v.size(), 4096, [&](long i, int) { key(i); });   // [pipe] the pipelined tile writer
    sort(keys.begin(), keys.end(), [&v](const SortKey& ka, const SortKey& kb) {
      if (ka.rt != kb.rt) return ka.rt < kb.rt;
      if (ka.mz != kb.mz) return ka.mz < kb.mz;
      if (ka.charge != kb.charge) return ka.charge < kb.charge;
      // fragment-sequence tiebreaks so equal-precursor spectra have a canonical order
      if (ka.n != kb.n) return ka.n < kb.n;
      const MSSpectrum& a = v[ka.idx];
      const MSSpectrum& b = v[kb.idx];
      for (Size i = 0; i < a.size(); ++i)
      {
        if (a[i].getMZ() != b[i].getMZ()) return a[i].getMZ() < b[i].getMZ();
        if (a[i].getIntensity() != b[i].getIntensity()) return a[i].getIntensity() < b[i].getIntensity();
      }
      // equal peaks: the fields that still serialise differently -- two overlapping isolation
      // windows can emit one precursor with the same peaks and different offsets (review, step 2)
      const auto& pa = a.getPrecursors(); const auto& pb = b.getPrecursors();
      if (pa.empty() != pb.empty()) return pa.empty();
      if (!pa.empty())
      {
        if (pa[0].getIsolationWindowLowerOffset() != pb[0].getIsolationWindowLowerOffset()) return pa[0].getIsolationWindowLowerOffset() < pb[0].getIsolationWindowLowerOffset();
        if (pa[0].getIsolationWindowUpperOffset() != pb[0].getIsolationWindowUpperOffset()) return pa[0].getIsolationWindowUpperOffset() < pb[0].getIsolationWindowUpperOffset();
        if (pa[0].getDriftTime() != pb[0].getDriftTime()) return pa[0].getDriftTime() < pb[0].getDriftTime();
      }
      for (const char* key : {"spx_guessed", "spx_n_isotopes"})
      {
        const int ia = a.metaValueExists(key) ? (int)a.getMetaValue(key) : -1, ib = b.metaValueExists(key) ? (int)b.getMetaValue(key) : -1;
        if (ia != ib) return ia < ib;
      }
      return false;
    });
    vector<MSSpectrum> ord; ord.reserve(v.size());
    for (const SortKey& k : keys) ord.push_back(std::move(v[k.idx]));
    v.swap(ord);
  }

  /// [twins] One precursor claimed twice is emitted as two spectra with complementary fragments: two seeds that
  /// resolved to the same monoisotope (the claim loop's `used[]` owns traces, not mono positions), or one
  /// precursor in the 1-m/z overlap of two isolation windows. Folds each such cluster into one spectrum.
  /// Relation: same MS1 frame (RT bit-equal; every precursor RT is a frame time), same charge > 0, precursor
  /// m/z within 20 ppm of the member earlier in `v`, |d1/K0| <= im_tol; clusters are its connected
  /// components. Survivor: most isotopes, then most peaks, then the earliest. Its peaks become the union of
  /// the cluster's, fused within 10 ppm (m/z intensity-weighted, intensity = max); every other field is its
  /// own. Same frame => same RT cell => same tile under every grouping, so the output stays grouping- and
  /// thread-invariant. This is bench/twin_merge.py, which priced it: docs/SAME-POSITION-DUPLICATION-2026-09-20.md.
  /// `v` in canonical order; survivors keep their places. Returns the number of spectra removed.
  /// [recap] `keyed`: every spectrum carries its peaks' cap keys as float data array 0; a fused peak keeps the max key.
  /// [twin-rt] `rt_tol` > 0 widens "same MS1 frame" to |dRT| <= rt_tol (assembly:twin_rt_frames): two calls of one
  /// monoisotope at nearby frames (seeds of one ion whose apexes differ, or v2's extension onto another precursor's
  /// monoisotope) fold like same-frame twins. Such a pair must also lie in one RT cell: the cell that owns a precursor
  /// (`cells`: the cell edges, `cell_drt`: tile_drt; the tile loop's own test, (T - rt) > cell_drt), so a pair never
  /// straddles a tile under any grouping and the output stays grouping- and thread-invariant. rt_tol 0 is the
  /// same-frame relation (|dRT| <= 0 and one frame is one cell), bit for bit. `st` (optional) adds {pairs linked at one
  /// frame, pairs linked across frames, near pairs refused at a cell edge, spectra the same-frame relation alone would
  /// have removed, clusters spanning more than rt_tol}. `near_mode`, how a pair from two frames must match: 0 precursor
  /// m/z within TWIN_PPM, as at one frame; 1 within BIN_PPM (one flight-time bin: one MS1 trace, or one bin at two
  /// frames, is the monoisotope of both); 2 as 1, and one of the two carries the userParam spx_mono_ext (its mono was
  /// placed by monocheck v2's extension), which mode 2 removes from every spectrum.
  static size_t mergeTwins_(vector<MSSpectrum>& v, double im_tol, bool keyed = false, double rt_tol = 0.0,
                            const vector<double>* cells = nullptr, double cell_drt = 0.0, std::array<size_t, 5>* st = nullptr,
                            int near_mode = 0)
  {
    const double TWIN_PPM = 20.0, FUSE_PPM = 10.0;   // ponytail: the measured values, no option until another is measured
    // [twin-rt] one flight-time bin across frames: the per-frame calibration moves one bin's m/z by < 0.01 ppm, distinct
    // monoisotope traces of one charge and 1/K0 at nearby frames sit >= 0.3 ppm apart (HeLa 10 ng: no pair in between)
    const double BIN_PPM = 0.05;
    struct K { double rt, mz, im; int z; uint32_t i, cell; bool ext; };
    auto cellOf = [&](double rt) -> uint32_t {   // the first edge T with (T - rt) > cell_drt; edges ascend, so the test is monotone
      if (!cells || cells->empty()) return 0;
      return (uint32_t)(std::partition_point(cells->begin(), cells->end(), [&](double T) { return !(T - rt > cell_drt); }) - cells->begin());
    };
    vector<K> ks; ks.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
      const auto& p = v[i].getPrecursors();
      if (!p.empty() && p[0].getCharge() > 0) ks.push_back({v[i].getRT(), p[0].getMZ(), p[0].getDriftTime(), p[0].getCharge(), (uint32_t)i, cellOf(v[i].getRT()),
                                                               near_mode == 2 && v[i].metaValueExists("spx_mono_ext")});
    }
    // one charge, then m/z: a pair is found from its lower-m/z member's 20-ppm window whatever its frames
    sort(ks.begin(), ks.end(), [](const K& a, const K& b) {
      if (a.z != b.z) return a.z < b.z;
      if (a.mz != b.mz) return a.mz < b.mz;
      if (a.rt != b.rt) return a.rt < b.rt;
      return a.i < b.i;
    });
    vector<uint32_t> up(v.size()), up1;   // up1: the same-frame relation alone, for the counter only
    for (size_t i = 0; i < up.size(); ++i) up[i] = (uint32_t)i;
    if (st) up1 = up;
    auto rootIn = [](vector<uint32_t>& u, uint32_t a) { while (u[a] != a) { u[a] = u[u[a]]; a = u[a]; } return a; };
    auto root = [&](uint32_t a) { return rootIn(up, a); };
    size_t n_same = 0, n_near = 0, n_edge = 0;
    for (size_t g0 = 0, g1; g0 < ks.size(); g0 = g1)
    {
      for (g1 = g0 + 1; g1 < ks.size() && ks[g1].z == ks[g0].z; ++g1) {}   // one charge
      for (size_t a = g0; a < g1; ++a)
      {
        // the probe is the EARLIER member of a pair, with its own tolerance, as in bench/twin_merge.py
        const double tol = ks[a].mz * TWIN_PPM * 1e-6, lo = ks[a].mz - tol, hi = ks[a].mz + tol;
        auto link = [&](size_t b) {
          if (ks[b].i <= ks[a].i || !(fabs(ks[b].im - ks[a].im) <= im_tol) || !(fabs(ks[b].rt - ks[a].rt) <= rt_tol)) return;
          const bool same = ks[b].rt == ks[a].rt;
          if (!same && near_mode >= 1 && fabs(ks[b].mz - ks[a].mz) > ks[a].mz * BIN_PPM * 1e-6) return;
          if (!same && near_mode == 2 && !ks[a].ext && !ks[b].ext) return;
          if (ks[b].cell != ks[a].cell) { ++n_edge; return; }
          ++(same ? n_same : n_near);
          { const uint32_t ra = root(ks[a].i), rb = root(ks[b].i); if (ra != rb) up[ra] = rb; }
          if (st && same) { const uint32_t ra = rootIn(up1, ks[a].i), rb = rootIn(up1, ks[b].i); if (ra != rb) up1[ra] = rb; }
        };
        for (size_t b = a + 1; b < g1 && ks[b].mz <= hi; ++b) link(b);
        for (size_t b = a; b-- > g0 && ks[b].mz >= lo; ) link(b);
      }
    }
    if (st)
    {
      vector<uint8_t> seen(v.size(), 0);
      size_t n_roots = 0;
      for (uint32_t i = 0; i < (uint32_t)v.size(); ++i) { const uint32_t r = rootIn(up1, i); if (!seen[r]) { seen[r] = 1; ++n_roots; } }
      (*st)[0] += n_same; (*st)[1] += n_near; (*st)[2] += n_edge; (*st)[3] += v.size() - n_roots;
    }
    // survivor per cluster: (isotopes, peaks) strictly greater replaces, so a tie keeps the earliest
    const uint32_t NONE = UINT32_MAX;
    vector<uint32_t> n_in(v.size(), 0), best(v.size(), NONE);
    for (uint32_t i = 0; i < (uint32_t)v.size(); ++i) ++n_in[root(i)];
    auto isotopes = [&v](uint32_t i) { return v[i].metaValueExists("spx_n_isotopes") ? (int)v[i].getMetaValue("spx_n_isotopes") : 0; };
    for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
    {
      const uint32_t r = root(i);
      if (n_in[r] < 2) continue;
      const uint32_t b = best[r];
      if (b == NONE || std::make_pair(isotopes(i), v[i].size()) > std::make_pair(isotopes(b), v[b].size())) best[r] = i;
    }
    vector<uint32_t> surv;
    vector<int> slot(v.size(), -1);
    for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
    {
      const uint32_t r = root(i);
      if (n_in[r] >= 2 && best[r] == i) { slot[i] = (int)surv.size(); surv.push_back(i); }
    }
    vector<vector<uint32_t>> members(surv.size());
    for (uint32_t i = 0; i < (uint32_t)v.size(); ++i) { const uint32_t r = root(i); if (n_in[r] >= 2) members[(size_t)slot[best[r]]].push_back(i); }
    if (st && rt_tol > 0.0)   // [twin-rt] clusters whose members span more than rt_tol: chained through intermediate frames
      for (const auto& m : members)
      {
        double r0 = v[m[0]].getRT(), r1 = r0;
        for (uint32_t i : m) { r0 = std::min(r0, v[i].getRT()); r1 = std::max(r1, v[i].getRT()); }
        if (r1 - r0 > rt_tol) ++(*st)[4];
      }
    // fuse: each thread reads and writes one cluster's spectra only
    auto fuse = [&](const long c)
    {
      if (keyed)
      {
        // [recap] the same fusion over (m/z, intensity, key) triples. Peaks equal in m/z and intensity add equal terms
        // and always land in one fused peak, so the peaks are the unkeyed fold's; each fused peak keeps its members' max key.
        vector<std::array<double, 3>> ps;
        for (uint32_t i : members[(size_t)c])
        {
          const auto& fd = v[i].getFloatDataArrays();
          for (Size j = 0; j < v[i].size(); ++j) ps.push_back({v[i][j].getMZ(), (double)v[i][j].getIntensity(), fd.empty() ? 0.0 : (double)fd[0][j]});
        }
        sort(ps.begin(), ps.end());
        vector<pair<double, double>> fused;
        OpenMS::DataArrays::FloatDataArray kk;
        kk.setName("spx_cap_key");
        const double NO_KEY = -std::numeric_limits<double>::infinity();
        double gm = 0.0, gw = 0.0, gi = 0.0, gk = NO_KEY;
        for (const auto& pk : ps)
        {
          if (gw != 0.0 && pk[0] - gm / gw > (gm / gw) * (FUSE_PPM * 1e-6)) { fused.emplace_back(gm / gw, gi); kk.push_back((float)gk); gm = gw = gi = 0.0; gk = NO_KEY; }
          gm += pk[0] * pk[1]; gw += pk[1]; gi = std::max(gi, pk[1]); gk = std::max(gk, pk[2]);
        }
        if (gw != 0.0) { fused.emplace_back(gm / gw, gi); kk.push_back((float)gk); }
        MSSpectrum& s = v[surv[(size_t)c]];
        s.clear(false);
        s.reserve(fused.size());
        for (const auto& f : fused) s.emplace_back(f.first, (float)f.second);
        s.getFloatDataArrays().clear();
        s.getFloatDataArrays().push_back(std::move(kk));   // moved, not copied (assign(1, x) copies)
        return;
      }
      vector<pair<double, double>> ps;
      for (uint32_t i : members[(size_t)c]) for (const auto& pk : v[i]) ps.emplace_back(pk.getMZ(), (double)pk.getIntensity());
      sort(ps.begin(), ps.end());
      vector<pair<double, double>> fused;
      double gm = 0.0, gw = 0.0, gi = 0.0;
      for (const auto& pk : ps)
      {
        if (gw != 0.0 && pk.first - gm / gw > (gm / gw) * (FUSE_PPM * 1e-6)) { fused.emplace_back(gm / gw, gi); gm = gw = gi = 0.0; }
        gm += pk.first * pk.second; gw += pk.second; gi = std::max(gi, pk.second);
      }
      if (gw != 0.0) fused.emplace_back(gm / gw, gi);
      MSSpectrum& s = v[surv[(size_t)c]];
      s.clear(false);   // peaks only: RT, precursor and userParams are the survivor's
      s.reserve(fused.size());
      for (const auto& f : fused) s.emplace_back(f.first, (float)f.second);
    };
    if (!omp_in_parallel())
    {
      std::exception_ptr fold_err;   // [task-err] an exception escaping the region (bad_alloc) would abort without a message
#pragma omp parallel for schedule(dynamic, 64)
      for (long c = 0; c < (long)surv.size(); ++c) { try { fuse(c); } catch (...) { taskErr_(fold_err); } }
      if (fold_err) std::rethrow_exception(fold_err);
    }
    else teamFor_((long)surv.size(), 64, [&](long c, int) { fuse(c); });   // [pipe] the pipelined tile writer; teamFor_ rethrows the first exception
    size_t w = 0;
    for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
    {
      const uint32_t r = root(i);
      if (n_in[r] >= 2 && best[r] != i) continue;
      if (w != i) v[w] = std::move(v[i]);
      ++w;
    }
    const size_t removed = v.size() - w;
    v.erase(v.begin() + (long)w, v.end());
    if (near_mode == 2) for (auto& sp : v) if (sp.metaValueExists("spx_mono_ext")) sp.removeMetaValue("spx_mono_ext");   // [twin-rt] never written
    return removed;
  }

  /// [recap] assembly:post_fusion_cap: the twin fold unions capped spectra, so a folded one can exceed max_frags. Keeps
  /// the max_frags peaks with the highest carried cap key (ties: higher intensity, then lower m/z), in m/z order, and
  /// drops the key array of every spectrum (it is never written). Per spectrum, so thread-invariant.
  /// Returns {spectra re-capped, peaks removed, largest spectrum before}.
  static std::array<size_t, 3> capAfterFusion_(vector<MSSpectrum>& v, Size max_frags)
  {
    size_t n_sp = 0, n_pk = 0, n_max = 0;
    // the counters are parameters, not captures: inside the reduction loop they must be the thread's copies
    auto capOne = [&v, max_frags](const long i, size_t& n_sp, size_t& n_pk, size_t& n_max)
    {
      MSSpectrum& s = v[(size_t)i];
      n_max = std::max(n_max, (size_t)s.size());
      if (s.size() > max_frags && !s.getFloatDataArrays().empty())
      {
        const auto& k = s.getFloatDataArrays()[0];
        vector<Size> idx(s.size());
        for (Size j = 0; j < idx.size(); ++j) idx[j] = j;
        partial_sort(idx.begin(), idx.begin() + (long)max_frags, idx.end(), [&](Size a, Size b) {
          if (k[a] != k[b]) return k[a] > k[b];
          if (s[a].getIntensity() != s[b].getIntensity()) return s[a].getIntensity() > s[b].getIntensity();
          if (s[a].getMZ() != s[b].getMZ()) return s[a].getMZ() < s[b].getMZ();
          return a < b;
        });
        idx.resize(max_frags);
        sort(idx.begin(), idx.end());   // the peaks are m/z-sorted: ascending positions keep that order
        n_pk += s.size() - max_frags; ++n_sp;
        s.getFloatDataArrays().clear();
        s.select(idx);
      }
      s.getFloatDataArrays().clear();
    };
    if (!omp_in_parallel())
    {
      std::exception_ptr cap_err;   // [task-err] an exception escaping the region (bad_alloc) would abort without a message
#pragma omp parallel for schedule(dynamic, 256) reduction(+ : n_sp, n_pk) reduction(max : n_max)
      for (long i = 0; i < (long)v.size(); ++i) { try { capOne(i, n_sp, n_pk, n_max); } catch (...) { taskErr_(cap_err); } }
      if (cap_err) std::rethrow_exception(cap_err);
    }
    else
    {
      // [pipe] the pipelined tile writer: per-slot partial counters, summed here
      vector<std::array<size_t, 3>> part((size_t)std::max(1, omp_get_num_threads()), std::array<size_t, 3>{0, 0, 0});
      teamFor_((long)v.size(), 256, [&](long i, int slot) { auto& p = part[(size_t)slot]; capOne(i, p[0], p[1], p[2]); });
      for (const auto& p : part) { n_sp += p[0]; n_pk += p[1]; n_max = std::max(n_max, p[2]); }
    }
    return {n_sp, n_pk, n_max};
  }

  /// [tile-2e] Streams spectra to an mzML in BLOCKS: the header once, then every block encoded in
  /// parallel chunks (the tool-side twin of patches/openms-mzml-parallel-write.patch) and appended in
  /// order, the index offsets rebased as each chunk lands. Blocks arrive in the canonical order (a
  /// tile's sorted spectra; tiles are RT-disjoint), so the file is in the order the bulk writer
  /// produced. With the total known up front (setExpectedSize) the count is written exactly and the
  /// bytes from <spectrumList> on are the bulk writer's; without it a zero-padded placeholder is
  /// overwritten at the end -- same width, so the indexedmzML offsets stay valid.
  class TileWriter : public PlainMSDataWritingConsumer
  {
  public:
    explicit TileWriter(const String& filename) : PlainMSDataWritingConsumer(filename) {}
    ~TileWriter() override
    {
      if (count_off_ >= 0 && ofs_.is_open())
      {
        const std::streampos end = ofs_.tellp();
        ofs_.seekp(count_off_);
        const char fill = ofs_.fill('0');
        ofs_ << std::setw(10) << spectra_written_;
        ofs_.fill(fill);
        ofs_.seekp(end);
      }
      // the base destructor closes <spectrumList> and writes the footer from spectra_offsets_
    }

    void writeBlock(const vector<MSSpectrum>& block)
    {
      if (block.empty()) return;
      if (!started_writing_)
      {
        MapType dummy; dummy = settings_; dummy.addSpectrum(block[0]);
        Internal::MzMLHandler::writeHeader_(ofs_, dummy, dps_, *validator_);
        started_writing_ = true;
      }
      if (!writing_spectra_)
      {
        ofs_ << "\t\t<spectrumList count=\"";
        if (spectra_expected_ > 0) ofs_ << spectra_expected_;
        else { count_off_ = ofs_.tellp(); ofs_ << "0000000000"; }
        ofs_ << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
        writing_spectra_ = true;
      }
      // The emitted spectra carry no native IDs, and the bulk writer therefore renumbered them as
      // "spectrum=<index>" (its uniqueness check trips on the duplicates); the running index here is
      // that same number, so renewing reproduces its bytes exactly.
      // A handler is not reentrant (writeSpectrum_ memoises CV-term validation in cached_terms_),
      // so each thread encodes with its own; dps_ is shared read-only -- every spectrum carries the
      // one data-processing entry the header registered, so writeSpectrum_ never appends to it.
      // [pipe] On the pipelined tile writer (inside the window team) the encoders are per puller slot of
      // teamFor_, at most the team size; the chunks and their order are the same either way.
      const bool in_team = omp_in_parallel();
      const Size CHUNK = 256;
      const Size n_thr = (Size)std::max(1, in_team ? omp_get_num_threads() : omp_get_max_threads());
      while (enc_.size() < n_thr) { enc_.emplace_back(new TileWriter("/dev/null")); enc_.back()->options_ = options_; }
      const Size wave = CHUNK * n_thr;
      for (Size w0 = 0; w0 < block.size(); w0 += wave)
      {
        const Size w1 = std::min(w0 + wave, block.size());
        const Size n_chunk = (w1 - w0 + CHUNK - 1) / CHUNK;
        vector<std::string> text(n_chunk);
        vector<vector<pair<std::string, Int64>>> offs(n_chunk);
        std::exception_ptr err;
        auto encode = [&](const long c, TileWriter& h)
        {
          try
          {
            std::ostringstream buf; buf.copyfmt(ofs_);          // identical formatting state, or the numbers differ
            const Size lo = w0 + (Size)c * CHUNK, hi = std::min(lo + CHUNK, w1);
            const Size o0 = h.spectra_offsets_.size();
            for (Size i = lo; i < hi; ++i) h.writeSpectrum_(buf, block[i], spectra_written_ + i, *h.validator_, /*renew_native_ids=*/true, dps_);
            text[(Size)c] = buf.str();
            offs[(Size)c].assign(h.spectra_offsets_.begin() + o0, h.spectra_offsets_.end());   // buffer-relative
            h.spectra_offsets_.resize(o0);
          }
          catch (...)
          {
            #pragma omp critical(tile_writer_err)
            if (!err) err = std::current_exception();
          }
        };
        if (!in_team)
        {
          #pragma omp parallel for schedule(dynamic, 1)
          for (long c = 0; c < (long)n_chunk; ++c) encode(c, *enc_[(Size)omp_get_thread_num()]);
        }
        else teamFor_((long)n_chunk, 1, [&](long c, int slot) { encode(c, *enc_[(Size)slot]); });
        if (err) std::rethrow_exception(err);
        for (Size c = 0; c < n_chunk; ++c)
        {
          const Int64 base = (Int64)ofs_.tellp();                 // where this chunk lands in the FILE
          ofs_ << text[c];
          for (const auto& p : offs[c]) spectra_offsets_.emplace_back(p.first, base + p.second);
          std::string().swap(text[c]);
        }
      }
      spectra_written_ += block.size();
    }

  private:
    std::streamoff count_off_ = -1;
    vector<std::unique_ptr<TileWriter>> enc_;   ///< per-thread encoders (never write their /dev/null)
  };
}

class TOPPDIAspeXtract : public TOPPBase
{
public:
  TOPPDIAspeXtract() :
    // official=false: a standalone tool, absent from ToolHandler's official list (TOPPBase rejects an unregistered official name).
    TOPPBase("DIAspeXtract", "Extracts pseudo-MS/MS spectra from diaPASEF (ion-mobility DIA) data.", false)
  {
    // Otherwise TOPPBase reports OpenMS's version as ours: version_ feeds the ini and the mzML software entry,
    // verboseVersion_ the --help Version line, which carries both numbers because a bug report needs the OpenMS build too.
    version_ = DIASPEXTRACT_VERSION;
    verboseVersion_ = String(DIASPEXTRACT_VERSION) + " (OpenMS " + VersionInfo::getVersion() + ")";
  }

  // [docurl] patches/openms-topp-external.patch makes this virtual. TOPPBase's default is an OpenMS doxygen page that
  // does not exist for a tool outside the release: the first line of --help and the docurl of every CTD 404'd.
  String getDocumentationURL() const override { return "https://github.com/okohlbacher/diaspextract"; }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input diaPASEF data (ion-mobility DIA; 1/K0 / VSSC): mzML, mzPeak, or a Bruker .d directory.");
#ifdef DIASPEXTRACT_WITH_MZPEAK
    setValidFormats_("in", {"mzML", "mzpeak", "d"});
#else
    setValidFormats_("in", {"mzML", "d"});      // stock OpenMS does not know the mzpeak format
#endif
    registerOutputFile_("out", "<file>", "", "Output pseudo-MS2 spectra. The FORMAT FOLLOWS THE "
                        "EXTENSION: .mzpeak (the default in a build with mzPeak support) or .mzML (the only format otherwise). mzPeak is columnar and much smaller; "
                        "mzML is what DDA search engines read today, so pass an .mzML name (or "
                        "-out_type mzML) if the next step is a search. mzPeak is written in one piece, "
                        "so it runs as one tile (the whole run's memory).");
#ifdef DIASPEXTRACT_WITH_MZPEAK
    setValidFormats_("out", {"mzpeak", "mzML"});
#else
    setValidFormats_("out", {"mzML"});
#endif
    registerStringOption_("out_type", "<type>", "", "Force the output format instead of taking it "
                          "from the extension of -out. Empty = follow the extension; mzPeak if the "
                          "extension says nothing (mzML in a build without mzPeak support). mzPeak output runs as one tile (see -out).", false);
    setValidStrings_("out_type", {"", "mzpeak", "mzML"});

    registerTOPPSubsection_("gate", "Precursor-to-fragment co-localization gate");
    registerDoubleOption_("gate:delta_im", "<1/K0>", 0.01, "Max |precursor - fragment| ion mobility difference (1/K0). (the reference implementation deltaApexIM default 0.01)", false);
    registerDoubleOption_("gate:delta_rt", "<sec>", 3.0, "Max |precursor - fragment| apex RT difference (seconds).", false);
    registerDoubleOption_("gate:min_correlation", "<0..1>", 0.3, "Min Pearson correlation of elution profiles. (the reference implementation ms1MS2Corr default 0.3; raise for cleaner but sparser spectra)", false);
    registerIntOption_("gate:min_correlation_points", "<n>", 3, "Min overlapping XIC points required to accept a correlation.", false, true);

    registerTOPPSubsection_("assembly", "Pseudo-spectrum assembly");
    registerIntOption_("assembly:min_fragments", "<n>", 3, "Emit a pseudo-spectrum only if it has at least this many fragments.", false);
    registerIntOption_("assembly:max_fragments", "<n>", 500, "Keep at most this many (top-ranked) fragments per pseudo-spectrum (applied before the twin fold, see assembly:twin_im_tolerance; assembly:post_fusion_cap re-applies it after the fold).", false);
    registerDoubleOption_("assembly:im_weight_sigma", "<1/K0>", 0.005, "Multiply each fragment's intensity by exp(-dIM^2 / 2 sigma^2), dIM its 1/K0 distance from the precursor; this also changes which fragments survive assembly:max_fragments. 0 = off, as in 1.5.0. Default 0.005 since 1.6.0: together with assembly:foreign_precursor_weight 0, +5.8 % Sage and +4.2 % MSFragger peptides over 1.5.0 on 43 public runs at unchanged entrapment (CHANGELOG, v1.6.0); earlier measurements in docs/BASELINE.md 'Confirmed WINS post-freeze'.", false);
    registerDoubleOption_("assembly:corr_power", "<k>", 2.0, "Multiply each emitted fragment's intensity by corr^k (k >= 0), corr its co-elution score with the precursor. It changes the emitted intensity, not the cap key; but assembly:fragment_fuse_ppm (on by default, mode anchored) groups the fragments strongest first by emitted intensity, so it can change which pieces fuse and with them which peaks survive assembly:max_fragments. 0 = off (the behaviour before 2026-07-28). Default 2: +8-10% Sage on three files (docs/BASELINE.md 'Confirmed WINS post-freeze') and +4% MSFragger on the one file measured with it ('CONTENT CANDIDATE 1 FALSIFIED').", false);
    registerDoubleOption_("assembly:twin_im_tolerance", "<1/K0>", 0.005, "Fold spectra that claim one precursor twice into one: same MS1 frame (or nearby frames, see assembly:twin_rt_frames), same charge, precursor m/z within 20 ppm, 1/K0 within this tolerance. The member with the most isotopes, then the most peaks, keeps its precursor and takes the union of the cluster's fragments (peaks within 10 ppm fused), so a folded spectrum can exceed assembly:max_fragments. Twins come from two seeds that resolve to one monoisotope, or from one precursor in the overlap of two isolation windows. -1 = off (the spectra before 1.4.0); 0 still folds twins of identical 1/K0. CHANGES OUTPUT. Default 0.005: -23.8% spectra for -0.21% Sage / -0.18% MSFragger peptides on PXD029836, entrapment FDR 0.98% with the fold against 0.98% without (docs/SAME-POSITION-DUPLICATION-2026-09-20.md).", false);
    registerIntOption_("assembly:twin_rt_frames", "<n>", 0, "Let the twin fold (assembly:twin_im_tolerance) also join two spectra of one monoisotope from nearby MS1 frames: charge:mono_profile_check v2's extension can move a monoisotope onto another precursor's monoisotope trace at another frame, which the same-frame fold cannot join. A pair folds across frames with the same charge, precursor m/z in one flight-time bin (within 0.05 ppm) and, by default, one of the two monoisotopes placed by the extension (assembly:twin_rt_match), 1/K0 within assembly:twin_im_tolerance, apexes at most n MS1 cycles apart (|dRT| <= (n + 0.5) x the cycle) and in one retention-time cell (tile:rt_sec; a pair across a cell edge stays two spectra). The survivor rule and the fusion are the fold's. 0 = auto: the run's median MS1 peak width (FWHM of the precursors' seed traces with at least 3 points) in MS1 cycles, rounded, at least 1. -1 = off: the same MS1 frame only (the fold before 1.6.0). Logged as [twins], stamped in spx:twin_merge. CHANGES OUTPUT.", false);
    registerStringOption_("assembly:twin_rt_match", "<mode>", "ext", "assembly:twin_rt_frames: which pairs from two different MS1 frames fold. 'ext': precursor m/z in one flight-time bin (within 0.05 ppm) and the monoisotope of at least one of the two placed by charge:mono_profile_check v2's extension, i.e. the extension moved it onto another precursor's monoisotope (nothing folds across frames without v2). 'bin': any pair in one flight-time bin (one MS1 trace, or one bin at two frames, is the monoisotope of both; distinct monoisotope traces of one charge and 1/K0 at nearby frames lie 0.3 ppm or more apart). 'ppm': within 20 ppm, as at one frame; it also folds the spectra of one ion seeded from different flight-time bins. 'bin' and 'ppm' cost peptides in the screen.", false, true);
    registerIntOption_("assembly:default_charge", "<z>", 2, "Charge given to a precursor left without a charge call (no isotope partner, or a z=1 call on the 3+ mobility band): > 0 assigns that charge, 0 leaves it unset for the search engine. These are the guessed precursors, so it has no effect unless -assembly:require_isotope_support false keeps them.", false);

    registerTOPPSubsection_("trace", "Mass-trace detection (see MassTraceDetection)");
    registerDoubleOption_("trace:mass_error_ppm", "<ppm>", 15.0, "m/z tolerance for trace detection.", false);
    registerDoubleOption_("trace:noise_threshold_int", "<int>", 100.0, "MS1 (precursor) noise intensity threshold.", false);
    registerDoubleOption_("trace:ms2_noise_threshold_int", "<int>", 10.0, "MS2 (fragment) noise intensity threshold.", false);
    registerDoubleOption_("trace:min_length_sec", "<sec>", -1.0, "Minimum MS1 mass-trace length in seconds (last minus first point). "
                          "Negative = auto (the default since 1.6.0): trace:min_length_fwhm x the run's median MS1 peak FWHM, at least half an MS1 cycle, so a trace needs 2 MS1 points. "
                          "The FWHM is that of the brightest tenth of a fixed sample of MS1 traces (1-Th m/z slices every 25 Th, traced like the run but with no minimum length). "
                          "Logged as [minlen], stamped spx:ms1_min_length. 0 = no length filter; 3 = the value before 1.6.0. CHANGES OUTPUT.", false);
    registerDoubleOption_("trace:min_length_fwhm", "<x>", 0.4, "trace:min_length_sec auto: the minimum MS1 trace length as this fraction of the run's median MS1 peak FWHM.", false, true);
    registerDoubleOption_("trace:ms1_chrom_peak_snr", "<x>", 3.0, "Apex signal-to-noise multiplier for MS1 traces (effective apex threshold = this x noise_threshold_int). OpenMS default 3.0.", false, true);
    registerDoubleOption_("trace:ms2_chrom_peak_snr", "<x>", 1.0, "Apex signal-to-noise multiplier for MS2 (fragment) traces: the apex threshold is this x trace:ms2_noise_threshold_int (OpenMS default 3.0). Lower recovers weaker fragments.", false);
    registerDoubleOption_("trace:ms2_min_length_sec", "<sec>", 0.0, "Minimum MS2 (fragment) mass-trace length (seconds). 0 = no length filter (the reference implementation-style relaxed MS2); MS1 keeps trace:min_length_sec.", false);
    registerTOPPSubsection_("tile", "Retention-time cells and tiles (integer detector)");   // without it --helphelp threw ElementNotFound
    registerDoubleOption_("tile:rt_sec", "<sec>", 600.0, "Pitch (seconds) of the fixed retention-time cells the integer detector traces in, cut at the run's MS1 frame times and stamped in the header (spx:tile_boundaries). The spectra are digest-identical for any grouping into tiles and any thread count at a fixed perf:trace_bands (the header stamps the grouping); a trace eluting across a cut is cut there. -1 = one cell (the whole run). CHANGES OUTPUT; the price of 600: docs/BASELINE.md '2026-09-09: the cell-atomic grid, step 2'.", false, true);
    registerStringOption_("trace:band_edges", "<mode>", "acquisition", "Where the integer detector's per-window flight-time band edges come from. 'acquisition': tlo = 0 and thi = the larger of the digitizer's last bin (GlobalMetadata DigitizerNumSamples) and MzAcqRangeUpper through the window's frame calibrations, run-level metadata a reader holding one tile can compute. A calibration tdf without those bounds is refused; an input with no tdf falls back to -trace:detector openms, which has no bands. 'slab': the window's resident peaks' extremes (the default before 2026-09-10); on a .d it needs -perf:stream_load false. Measured neutral: docs/BASELINE.md 'Cell grid, step 3a'.", false);
    setValidStrings_("trace:band_edges", {"slab", "acquisition"});
    registerIntOption_("tile:cells_per_tile", "<n>", 1, "Group the retention-time cells into TILES of n cells each and run the window loop tile by tile: "
                       "detection on the tile's frames only, scoring of the precursors the tile owns (an MS1 frame time in [T_k - gate:delta_rt, T_{k+1} - gate:delta_rt)) "
                       "with the previous tile's boundary fragments carried over, the tile's spectra sorted and written before the next tile starts. "
                       "Output is digest-identical for any tile count (the cell grid); memory is one tile's. At the defaults TNBC 009 peaks at 99.70 GB in one "
                       "tile and 21.07 GB in 13 tiles with a trim between tiles (wall 9:42 -> 9:57, same node), dataset D at 37.44 -> 21.99 GB in 3 tiles "
                       "(3:15 -> 3:47). 0 = all cells in one tile; mzPeak output always runs as one tile (it is written in one piece).", false, true);
    registerStringOption_("trace:split_scan_time", "<mode>", "frame", "Scan time behind the valley splitter's window (chrom_fwhm / scan time). 'frame': the run's MS1 cycle time, one value for every trace. 'trace': each trace's own average cycle time (OpenMS's behaviour; the default before 2026-09-09, Sage -5.7% on dataset D): docs/BASELINE.md '2026-09-09: the valley splitter's scan time'.", false, true);
    setValidStrings_("trace:split_scan_time", {"frame", "trace"});
    registerStringOption_("trace:detector", "<mode>", "integer", "Mass-trace detector. 'integer' works on the instrument's flight-time bins "
                          "and needs the vendor calibration (without it: 'openms', logged). 'openms' runs OpenMS MassTraceDetection. Different "
                          "algorithms: they share ~85% of identified peptides, neither wins on every file and engine, integer uses ~40% less memory. docs/MZ-AXIS-DESIGN.md.", false);
    setValidStrings_("trace:detector", {"openms", "integer"});
    registerStringOption_("trace:mz_estimator", "<mode>", "apex", "Reported m/z of a mass trace: 'apex' the most intense peak's, 'mean' the intensity-weighted centroid (OpenMS), 'median' the median. CHANGES OUTPUT. Fragment traces of the integer detector always report the calibrated m/z of their apex bin, whatever this says; with trace:detector openms the estimator applies to them too. apex measured best on both engines: docs/BASELINE.md 'APEX estimator: 3-FILE RESULT'.", false);
    setValidStrings_("trace:mz_estimator", {"apex", "mean", "median"});
    registerDoubleOption_("trace:ms2_split_valleys", "<chrom_fwhm>", 7.0, "Split MS2 mass traces at chromatographic valleys (ElutionPeakDetection, width filtering off). Value = chrom_fwhm (seconds); 0 = off.", false);
    registerFlag_("diag:selftest_twins", "Run assertions on the twin fold (mergeTwins_) and exit: the relation (also across nearby frames), the cluster closure, the survivor order, the 10-ppm fusion and the keyed fold with the post-fusion cap (assembly:post_fusion_cap)", true);
    registerFlag_("diag:selftest_arena", "Run assertions on the MS1 arena compaction (compactUnreferenced) and exit: kept spans must survive a container order that differs from their offset order", true);
    registerDoubleOption_("trace:ms1_split_valleys", "<chrom_fwhm>", 7.0, "Split MS1 mass traces at chromatographic valleys (ElutionPeakDetection), so two peptides eluting apart within the m/z tolerance do not merge into one precursor with the wrong monoisotope and charge. Value = chrom_fwhm (seconds); 0 = off.", false);
    registerDoubleOption_("trace:max_span_sec", "<sec>", 120.0, "Trim each mass trace (precursor and fragment) to at most this many seconds "
                          "around its apex, after detection and valley splitting. 0 = no cap. CHANGES OUTPUT (peptides within 0.08% from 0 "
                          "to 240 on dataset D: docs/BASELINE.md '0c. trace:max_span_sec sweep').", false);
    setMinFloat_("trace:max_span_sec", 0.0);

    // NB: "threads" is reserved by TOPPBase (the standard -threads option), so this lives under "perf".

    registerStringOption_("perf:malloc_trim", "<true/false>", "true",
                          "Return the allocator's free pages to the OS (glibc only) at the two phase boundaries and between tiles. "
                          "The window loop's allocations otherwise stack on top of everything the loading and MS1 "
                          "phases freed but the allocator kept: measured on a 2-h acquisition in one tile that floor is 138 GB "
                          "at the loop start, and trimming it first takes the run's PEAK from 265 to 186 GB. "
                          "IT IS NOT FREE: the loop then faults those pages back in cold (minor faults +26%), which "
                          "costs about 8.5% of wall time -- so this trades ~30% of peak memory for ~8% of runtime, "
                          "measured interleaved on one node. Between tiles (after tile k is written, before tile k+1 is read) it "
                          "takes about 4 GB off TNBC 009 and 0.5 GB off dataset D at the defaults. On by default because memory is "
                          "what bounds this tool (it decides how many isolation windows fit at once, and whether the run fits the "
                          "machine at all). Turn it off when time is the binding constraint. The spectra are digest-identical "
                          "either way.", false);
    setValidStrings_("perf:malloc_trim", {"true", "false"});
    registerStringOption_("perf:stream_load", "<true/false>", "true",
                          "Read the input frame by frame, picking each frame on arrival: a Bruker .d, or an .mzpeak in a build with mzPeak "
                          "support (mzML is read whole either way). false reads the whole run first: ~1.75x the memory and ~1.7x the wall on a "
                          "large .d, and CHANGES OUTPUT (+2.0% Sage peptides at unchanged entrapment FDR): docs/BASELINE.md 'THE READER QUESTION, SETTLED'.", false);
    setValidStrings_("perf:stream_load", {"true", "false"});
    registerIntOption_("trace:native_ms1_neighbors", "<n>", 0, "Sum n neighbouring MS1 frames on each side in the .d reader, before peak picking: 0 = off, 1 = 3-frame sum, 2 = 5-frame sum. Needs perf:stream_load true and -dnoise:ms1 false. CHANGES OUTPUT: at 1 (2026-07-27, older defaults) Sage +7.8% for 2.7x the spectra.", false);

    registerTOPPSubsection_("perf", "Concurrency / memory tradeoff");
    registerStringOption_("gate:coelution", "<mode>", "pearson", "Fragment-precursor co-elution statistic: 'pearson' over the full RT grid, or 'logoverlap' (AlphaDIA-style profile overlap in log space where the precursor has signal; falsified: docs/BASELINE.md 'What is already falsified'). The scales differ: retune gate:min_correlation when switching. 'logoverlap' needs gate:shrinkage off (the default fisher refuses it).", false);
    setValidStrings_("gate:coelution", {"pearson", "logoverlap"});
    registerStringOption_("gate:short_traces", "<mode>", "concentric", "Fragment traces with fewer than gate:min_correlation_points points, or fewer overlapping frames: 'off' drops them; 'concentric' keeps one whose detected RT span lies inside the precursor's (one cycle of slack) and whose span centre is within max(1 cycle, gate:short_center_tol x the precursor half-span) of the precursor's. Faint ions are the tips of their elution peaks: their spans shrink with intensity but stay centred on the apex. Experimental.", false, true);
    setValidStrings_("gate:short_traces", {"off", "concentric"});
    registerDoubleOption_("gate:short_center_tol", "<frac>", 0.5, "gate:short_traces concentric: centre tolerance as a fraction of the precursor half-span (at least one cycle).", false, true);
    setMinFloat_("gate:short_center_tol", 0.0);
    registerIntOption_("gate:short_min_points", "<n>", 1, "gate:short_traces concentric: fewest points a short fragment trace needs.", false, true);
    setMinInt_("gate:short_min_points", 1);
    registerDoubleOption_("gate:fragment_delta_im", "<1/K0>", -1.0, "Max |precursor - fragment| 1/K0 difference in the co-elution gate alone; negative = gate:delta_im, which also sets the tracing tolerance (widening that merges ions of different mobility). Experimental.", false, true);
    registerStringOption_("gate:anchor_rescue", "<true/false>", "false", "Rescue weak fragments of well-sampled precursors (signal in more than gate:faint_max_frames window frames, the faint test of gate:adaptive in every mode). When at least gate:anchor_min_count fragments pass the gate with a correlation of gate:anchor_min_corr or more, their profiles, each scaled to unit maximum, are summed into a consensus over the window frames. A fragment trace inside the RT and 1/K0 windows that the gate rejected only for too few shared points or by gate:min_correlation, with at least 2 points, then joins when its apex is within one cycle of the consensus apex and its points correlate with the consensus at those frames (Pearson over its own points, not the full grid) by gate:anchor_rescue_corr or more; a 2-point trace instead needs both points in the consensus core (at least half its maximum) rising or falling with it, and scores gate:anchor_rescue_corr. The emitted intensity and the cap key use that score, capped at gate:anchor_min_corr; at most gate:anchor_max_rescued per spectrum, best score first. Logged as [anchor], stamped spx:anchor. Experimental.", false, true);
    setValidStrings_("gate:anchor_rescue", {"true", "false"});
    registerIntOption_("gate:anchor_min_count", "<n>", 3, "gate:anchor_rescue: fewest anchors (accepted fragments at gate:anchor_min_corr or more) a precursor needs for the rescue.", false, true);
    setMinInt_("gate:anchor_min_count", 1);
    registerDoubleOption_("gate:anchor_min_corr", "<r>", 0.8, "gate:anchor_rescue: least gate correlation of an anchor; also caps the correlation score of a rescued fragment (its cap key can still exceed an anchor's through intensity and shared frames).", false, true);
    setMinFloat_("gate:anchor_min_corr", 0.0);
    setMaxFloat_("gate:anchor_min_corr", 1.0);
    registerDoubleOption_("gate:anchor_rescue_corr", "<r>", 0.7, "gate:anchor_rescue: least Pearson correlation of a rejected fragment's points with the anchor consensus.", false, true);
    setMinFloat_("gate:anchor_rescue_corr", 0.0);
    setMaxFloat_("gate:anchor_rescue_corr", 1.0);
    registerIntOption_("gate:anchor_max_rescued", "<n>", 50, "gate:anchor_rescue: most rescued fragments per spectrum (best score first, ties by m/z).", false, true);
    setMinInt_("gate:anchor_max_rescued", 0);
    registerDoubleOption_("trace:im_slice_merge", "<1/K0>", 0.015, "Integer detector: before tracing, merge a frame's picked MS2 peaks that fall in one flight-time bin into one peak (intensity summed, 1/K0 intensity-weighted). Peaks join when their 1/K0 is within this distance of the previous one and the group spans at most twice this distance: the IM picker (box gate:delta_im) cuts an ion with a wider mobility peak into several same-bin centroids, which would otherwise become separate traces and duplicate m/z in a spectrum. Negative = off.", false, true);
    setMaxFloat_("trace:im_slice_merge", 1.0);   // 1/K0 units: anything wider than the mobility axis is a typo
    registerStringOption_("assembly:merge_same_bin", "<true/false>", "true", "Merge fragments of one flight-time bin (one m/z; exact m/z on the OpenMS detector) in a spectrum into one peak before assembly:max_fragments: members overlapping the strongest one in retention time are summed, RT-disjoint pieces dropped; the peak keeps the strongest member's m/z and the best cap key of the summed members. Catches what trace:im_slice_merge cannot (valley splits, band edges, tile carries, slices beyond its span). Two peaks at one m/z break the fragment matching of OpenMS-based search engines.", false, true);
    setValidStrings_("assembly:merge_same_bin", {"true", "false"});
    registerStringOption_("assembly:post_fusion_cap", "<true/false>", "false", "Re-apply assembly:max_fragments after the twin fold (assembly:twin_im_tolerance), whose union of the members' capped fragments can exceed it. Every peak carries its assembly:cap_rank key (as a 32-bit float) through the same-bin merge and the fold (a fused peak keeps its members' best key); a spectrum over the cap keeps its top-ranked peaks (ties: higher intensity, then lower m/z). Logged as [recap], stamped spx:post_fusion_cap. CHANGES OUTPUT. Experimental.", false, true);
    setValidStrings_("assembly:post_fusion_cap", {"true", "false"});
    registerStringOption_("assembly:cap_rank", "<mode>", "evidence", "Which fragments survive assembly:max_fragments: 'corr_inten' ranks by gate score x intensity (the score is the shrunk or concentric one when those modes are on); 'umpire' by max(r,0)^2 x ln(1+I) with r the raw co-elution score before shrinkage (Pearson, or the logoverlap score under gate:coelution logoverlap; DIA-Umpire's key); 'evidence' by max(r,0) x ln(1+I) x (1 + overlapping frames), so long, correlated support outranks an intense one-frame coincidence. Experimental.", false, true);
    setValidStrings_("assembly:cap_rank", {"corr_inten", "umpire", "evidence"});
    registerStringOption_("assembly:frag_intensity", "<mode>", "apex", "The base amplitude of each emitted fragment, before the 1/K0 weight (assembly:im_weight_sigma) and corr^assembly:corr_power: 'apex' the fragment trace's maximum in one frame; 'core' its intensity summed over the assembly:frag_core_frames window frames centred on the precursor's apex frame (the same frames for every fragment of the precursor; a missing point counts 0); 'sum' its intensity summed over the precursor's span, the window frames from its first to its last point (gaps included). A fragment with no point in those frames takes its point nearest the precursor's apex frame, so none is dropped or written at 0. The gate and the assembly:cap_rank key (still the apex) are unchanged; the later steps that weigh by intensity follow the new values (the same-bin merge's strongest member, the anchored grouping and the sums of assembly:fragment_fuse_ppm, the twin fold's intensity-weighted m/z), so peaks can move, fuse differently and change which survive the cap. Experimental.", false, true);
    setValidStrings_("assembly:frag_intensity", {"apex", "core", "sum"});
    registerIntOption_("assembly:frag_core_frames", "<n>", 3, "assembly:frag_intensity core: the number of window frames (odd) centred on the precursor's apex frame.", false, true);
    setMinInt_("assembly:frag_core_frames", 1);
    setMaxInt_("assembly:frag_core_frames", 15);
    registerDoubleOption_("assembly:foreign_precursor_weight", "<w>", 0.0, "Multiply the fragments that lie on the residual isotope envelope of ANOTHER precursor isolated with this one by w, before assembly:max_fragments (0 removes them; 1 = off, as in 1.5.0); the cap key follows the attenuated intensity. Default 0 since 1.6.0. A foreign precursor is one of the run's isotope-supported precursors whose mono lies in the same isolation window, with its apex within the gate's RT window of this precursor's, 1/K0 within assembly:foreign_precursor_dim, and neither a twin (mono within 0.01) nor this ion called again (either mono on the other's isotope grid, up to 3 steps, within the larger of assembly:foreign_precursor_ppm and 30 ppm). Its envelope counts when at least assembly:foreign_precursor_min_iso of its positions mono + k x 1.00335/z, k = 0..3, carry a fragment within assembly:foreign_precursor_ppm. Fragments within that larger tolerance of this precursor's own envelope (mono to M+3) are never touched. Logged as [foreign], stamped spx:foreign. CHANGES OUTPUT.", false, false);
    setMinFloat_("assembly:foreign_precursor_weight", 0.0);
    setMaxFloat_("assembly:foreign_precursor_weight", 1.0);
    registerDoubleOption_("assembly:foreign_precursor_ppm", "<ppm>", 10.0, "assembly:foreign_precursor_weight: m/z tolerance between a foreign envelope position and a fragment (ppm).", false, true);
    setMinFloat_("assembly:foreign_precursor_ppm", 0.0);
    registerDoubleOption_("assembly:foreign_precursor_dim", "<1/K0>", 0.015, "assembly:foreign_precursor_weight: largest 1/K0 distance between this precursor and a foreign one.", false, true);
    setMinFloat_("assembly:foreign_precursor_dim", 0.0);
    registerIntOption_("assembly:foreign_precursor_min_iso", "<n>", 2, "assembly:foreign_precursor_weight: how many of a foreign envelope's four positions (mono to M+3) must carry a fragment before its fragments are attenuated.", false, true);
    setMinInt_("assembly:foreign_precursor_min_iso", 1);
    setMaxInt_("assembly:foreign_precursor_min_iso", 4);
    registerDoubleOption_("assembly:fragment_fuse_ppm", "<ppm>", 25.0, "Fuse a spectrum's fragments that lie within this many ppm into one peak, after the same-bin merge and before assembly:min_fragments and assembly:max_fragments (0 = off, as in 1.5.0; default 25 with assembly:fragment_fuse_mode anchored since 1.6.0): the IM picker keeps 1-2 flight-time bins of an ion that spans 4-7, so one fragment ion becomes several traces and several peaks a few ppm apart. A fused peak's intensity is the sum of its members' emitted intensities (whatever their retention-time overlap, unlike the same-bin merge), its m/z their intensity-weighted mean and its cap key (assembly:cap_rank) the key of one fragment applied to the whole ion: the best-keyed member's score (and, under evidence, its shared frames) with the sum of the members' IM-weighted apex intensities, never below that member's own key (non-negative keys, up to rounding); a fused ion takes one cap slot and the slots it frees refill. A spectrum that fuses below assembly:min_fragments is not written. Grouping by assembly:fragment_fuse_mode. Logged as [fuse], stamped spx:fragment_fuse. CHANGES OUTPUT.", false, false);
    setMinFloat_("assembly:fragment_fuse_ppm", 0.0);
    setMaxFloat_("assembly:fragment_fuse_ppm", 50.0);
    registerStringOption_("assembly:fragment_fuse_mode", "<mode>", "anchored", "assembly:fragment_fuse_ppm: 'anchored' takes the fragments by emitted intensity, strongest first (ties: lower m/z); each one not yet fused becomes an anchor and takes every unfused fragment within the tolerance of its own m/z. 'linkage' walks the fragments in m/z order and joins a fragment to the group of the previous one when it lies within the tolerance of it and within twice the tolerance of the group's first member (single linkage, span capped).", false, true);
    setValidStrings_("assembly:fragment_fuse_mode", {"anchored", "linkage"});
    registerStringOption_("gate:adaptive", "<mode>", "faint", "'faint': apply gate:shrinkage, gate:colocation, gate:short_traces, gate:faint_delta_rt and gate:faint_corr_power only to faint precursors, those with signal in at most gate:faint_max_frames window frames (a faint ion is the short tip of its elution peak); every other precursor keeps the plain gate (gate:delta_rt, gate:min_correlation_points, assembly:corr_power). 'auto': 'faint' when at least gate:auto_faint_fraction of the run's precursors have at most gate:faint_max_frames MS1 trace points, else the plain gate for every precursor (gate:shrinkage, gate:colocation and gate:short_traces off; the choice is stamped as spx:gate_adaptive auto:on|auto:off). 'off': gate:shrinkage, gate:colocation and gate:short_traces apply to every precursor at gate:delta_rt and assembly:corr_power (gate:faint_delta_rt and gate:faint_corr_power unused). Experimental.", false, true);
    setValidStrings_("gate:adaptive", {"off", "faint", "auto"});
    registerDoubleOption_("gate:auto_faint_fraction", "<frac>", 0.45, "gate:adaptive auto: the run gets the faint gate when at least this fraction of its precursors span at most gate:faint_max_frames MS1 points (a sparse run: HeLa 10-100 ng 0.57-0.69, dense tissue 0.21); otherwise the plain gate. Logged and stamped spx:gate_adaptive.", false, true);
    setMinFloat_("gate:auto_faint_fraction", 0.0);
    setMaxFloat_("gate:auto_faint_fraction", 1.0);
    registerIntOption_("gate:faint_max_frames", "<n>", 4, "gate:adaptive faint: most window frames with precursor signal for a faint precursor.", false, true);
    setMinInt_("gate:faint_max_frames", 1);
    registerDoubleOption_("gate:faint_delta_rt", "<sec>", 6.0, "gate:adaptive faint: the apex RT window for faint precursors.", false, true);
    setMinFloat_("gate:faint_delta_rt", 0.0);
    registerDoubleOption_("gate:faint_corr_power", "<k>", 1.0, "gate:adaptive faint: assembly:corr_power for faint precursors.", false, true);
    setMinFloat_("gate:faint_corr_power", 0.0);
    registerStringOption_("gate:colocation", "<mode>", "ranges", "'ranges': keep a fragment only if its apex lies inside the precursor's detected RT span and the precursor apex inside the fragment's (one cycle of slack), and the spans overlap by at least gate:colocation_overlap of the shorter one (DIA-Umpire's apex-in-range and RT-overlap tests). Experimental.", false, true);
    setValidStrings_("gate:colocation", {"off", "ranges"});
    registerDoubleOption_("gate:colocation_overlap", "<frac>", 0.3, "gate:colocation ranges: least span overlap as a fraction of the shorter span.", false, true);
    setMinFloat_("gate:colocation_overlap", 0.0);
    registerStringOption_("gate:shrinkage", "<mode>", "fisher", "'fisher': gate on the posterior mean of the co-elution score instead of a hard gate:min_correlation_points cut: z = atanh(r), r the gate:coelution pearson statistic (full grid by default), with the overlapping frames n as a heuristic effective sample size (variance 1/max(n-3,1)) and prior N(atanh(gate:shrink_prior_r), gate:shrink_prior_sd^2). Scores on few frames are pulled toward the prior, so a prior below gate:min_correlation demands more from short traces; a prior equal to it reduces the mode to a plain cut at gate:shrink_min_points. Not with gate:coelution logoverlap. Experimental. Replaces gate:min_correlation_points with gate:shrink_min_points.", false, true);
    setValidStrings_("gate:shrinkage", {"off", "fisher"});
    registerDoubleOption_("gate:shrink_prior_r", "<r>", 0.6, "gate:shrinkage fisher: prior mean correlation; 0.6 shrinks toward moderate co-elution, 0 toward none (and loses most faint fragments).", false, true);
    registerDoubleOption_("gate:shrink_prior_sd", "<sd>", 0.3, "gate:shrinkage fisher: prior standard deviation in Fisher-z units (smaller = stronger shrinkage).", false, true);
    setMinFloat_("gate:shrink_prior_sd", 0.01);
    registerIntOption_("gate:shrink_min_points", "<n>", 1, "gate:shrinkage fisher: fewest overlapping frames to compute a correlation; at 1 a one-frame overlap is scored almost at the prior (the data get weight 1/(1 + 1/gate:shrink_prior_sd^2)), so pair it with gate:adaptive faint.", false, true);
    setMinInt_("gate:shrink_min_points", 1);
    registerFlag_("gate:variance_support", "Compute the Pearson co-elution score over the union support of the two profiles instead of the full RT grid (gate:coelution pearson only; the scale changes, retune gate:min_correlation). Falsified 2026-07-27, kept for reproducibility: docs/BASELINE.md 'What is already falsified'.");
    registerIntOption_("perf:ms1_trace_bands", "<n>", 48, "Split MS1 mass-trace detection into N m/z bands traced concurrently (the halo partition of perf:trace_bands). Not exact: the band count changes the spectrum list slightly. Default 48 since 2026-09-07, gated on datasets D and A: docs/BASELINE.md 'ms1_trace_bands=48 adopted'. 1 = off; 0 is read as 1.", false);
    registerStringOption_("perf:ms1_prune", "<true/false>", "true",
                          "Drop picked MS1 centroids at or below trace:noise_threshold_int as they load, keeping the witnesses that leave every MS1 "
                          "band spectrum as it was. Digest-identical spectra, smaller MS1 map (peak at one cell per tile, dnoise off: dataset D 31.99 -> "
                          "25.90 GB, TNBC 009 69.11 -> 31.97 GB). false = keep every centroid. docs/BASELINE.md 'MS1 prune, output-identical'.", false);
    setValidStrings_("perf:ms1_prune", {"true", "false"});
    registerIntOption_("perf:trace_bands", "<n>", 12, "Split each window's fragment m/z range into N bands traced concurrently, lifting parallelism past the window count. The band count changes the output: docs/BASELINE.md 'perf:trace_bands 48: FAILS the gate'. 1 = off. 0 = auto (windows x bands >= threads), which the integer detector refuses at a positive tile:rt_sec (the default) and on the .d streaming source.", false);
    registerDoubleOption_("perf:mem_fraction", "<f>", 0.75, "Fraction of currently-FREE RAM (/proc/meminfo MemAvailable) the window loop may commit. Concurrency is re-decided at every window admission, so the run adapts to other jobs on a shared node. One window is always admitted regardless, so a single oversized window cannot deadlock.", false);
    registerIntOption_("perf:max_concurrent_windows", "<n>", 0, "UPPER bound on isolation windows processed concurrently (memory admission may run fewer). Peak RAM scales with windows IN FLIGHT (each holds its frames + traces + grid), not with the total window count, so lowering this trades wall time for RAM. 0 = unlimited (use all threads).", false);

    registerTOPPSubsection_("diag", "Diagnostics ('debug' is reserved by TOPP)");
    registerStringOption_("diag:dump_ms1_tsv", "<prefix>", "", "Write MS1 traces and inferred precursors to <prefix>.traces.tsv / <prefix>.precursors.tsv so precursor loss can be attributed to a stage (no trace / no hypothesis / no spectrum) rather than inferred from the final count. Empty = off.", false, true);

    registerIntOption_("max_charge", "<n>", 5, "Maximum precursor charge considered during isotope inference.", false);

    // NOTE: the subsection MUST be registered or printUsage_() throws ElementNotFound — which turns
    // any user parameter typo into a FATAL crash instead of a usage message. [bugfix]
    registerTOPPSubsection_("charge", "Isotope / charge-state inference");
    registerIntOption_("charge:min_charge", "<n>", 1, "Minimum precursor charge to emit. Default 1 since 2026-09-08: identified z=1 peptides on tryptic data are genuine 1+ ions (mobility on the 1+ trend line; 14-20% of Sage peptides on a 2-hour acquisition), but their stratum carries a higher entrapment FDR than a pooled 1% cut implies, so control FDR per charge downstream. "
                       "-charge:min_charge 2 -charge:im_charge_veto false restores the old charge handling (the veto runs first). docs/BASELINE.md 'charge:min_charge 1 + the ion-mobility charge veto'.", false);
    registerStringOption_("charge:scoring", "<mode>", "count", "Charge/monoisotope inference: 'count', the partner-count walk (ties favour the lower charge); 'envelope', averagine cosine x isotope co-elution (the default until 2026-07-21; 49.7% vs 71.8% charge agreement with a DIA-NN reference on dataset B). docs/charge-inference.md.", false);
    setValidStrings_("charge:scoring", {"envelope", "count"});
    registerStringOption_("charge:im_charge_veto", "<true/false>", "true", "Correct charge halving by ion mobility: fit the run's 2+ and 3+ 1/K0 trend lines from its confident calls (>= 3 isotopes); a z=1 call within charge:im_veto_band sigma of the 2+ line (and nearer it than the 3+ line) is re-called 2+, "
                          "one on the 3+ line loses its charge and becomes a guessed precursor (dropped by assembly:require_isotope_support). Genuine 1+ ions sit on their own line and stay. docs/BASELINE.md 'The veto's 3+ arm dropped'.", false);
    setValidStrings_("charge:im_charge_veto", {"true", "false"});
    registerDoubleOption_("charge:im_veto_band", "<sigma>", 2.5, "Half-width of the 2+/3+ mobility band, in MAD-sigma of "
                          "the confident calls' residuals (measured residual sd ~0.037 1/K0 against a 1+/2+ separation of ~0.29).", false);
    registerStringOption_("assembly:require_isotope_support", "<true/false>", "true",
                          "Drop precursor hypotheses with no isotope partner (the guessed singletons). false keeps them: on dataset D (2026-09-04) "
                          "that doubled the spectra, ran 6.9x slower and identified fewer peptides on both engines. docs/BASELINE.md (defaults audit, top of file).", false);
    setValidStrings_("assembly:require_isotope_support", {"true", "false"});
    registerDoubleOption_("charge:mono_averagine_guard", "<slack>", 0.0, "Reject a leftward isotope step of the charge:scoring=count walk when the lighter peak is below slack x heavier peak / lambda, averagine's loosest ratio (lambda = 0.000594 x neutral mass; 1.0 = exactly that bound, smaller values add slack). 0 = off. Falsified: it trades the -1.003 Da open-search artefact for a +1.003 Da one, docs/BASELINE.md 'Monoisotope -1.003 Da open-search artefact'.", false);
    registerFlag_("charge:mono_averagine_select", "Choose the monoisotope by averagine cosine over the whole isotope run of the charge:scoring=count walk instead of its furthest-left peak; only the reported precursor mass moves. Missed its pre-registered criterion: docs/BASELINE.md 'Monoisotope -1.003 Da open-search artefact'.");
    registerStringOption_("charge:mono_profile_check", "<false/true/v2>", "v2", "Check the monoisotope of the charge:scoring=count walk, for the winning charge only (default v2 since 1.6.0; false = off, as in 1.5.0). "
                          "'true' (v1): a lighter trace becomes the monoisotope only if its MS1 elution profile correlates with the seed trace's (Pearson >= charge:mono_profile_min_corr over the frames both have; below 3 shared frames the correlation test passes unchecked) and its 1/K0 is within charge:mono_profile_max_dim of the seed's. "
                          "The walk stops at the first failed step: the monoisotope is the last accepted peak, and the rejected trace and the ones beyond it stay unclaimed for other precursors. "
                          "'v2': isotope-position intensity sums (charge:mono_position_*). PS(k) = summed apex intensity of every MS1 trace, claimed or not, at the seed's isotope position k (the partner tolerance in m/z, gate:delta_rt, charge:mono_position_im_box 1/K0 of the seed, by default from the run's own isotope scatter). "
                          "The walk's leftmost step j is rejected when PS(-j) / PS(-j+1) x lambda < charge:mono_position_min_ratio (lambda = max(0.1, 0.000594 x neutral mass), averagine's M-1/M0 bound; an empty heavier position counts as 0); "
                          "the monoisotope moves one isotope right, the rejected trace stays claimed, and (charge:mono_position_iterate true) the new leftmost step is tested in turn. If nothing was rejected and the walk ended on a miss (not a charge:mono_averagine_guard stop), the next lighter position is accepted without claiming when its ratio is >= charge:mono_position_extend (0 = off). "
                          "v2 judges each walk at its precursor's final charge, after charge:im_charge_veto: a z=1 walk re-called 2+ stepped 1.003 m/z, two isotopes at 2+, and is read at 2+ spacing. "
                          "Both 'true' and 'v2' keep the charge, the seed and the envelope size (spx_n_isotopes). Logged as [monocheck], stamped spx:monocheck. CHANGES OUTPUT.", false, false);
    setValidStrings_("charge:mono_profile_check", {"false", "true", "v2"});
    registerDoubleOption_("charge:mono_profile_min_corr", "<r>", 0.6, "charge:mono_profile_check true: minimum Pearson correlation of a leftward partner's MS1 profile with the seed's.", false, true);
    registerDoubleOption_("charge:mono_profile_max_dim", "<1/K0>", 0.015, "charge:mono_profile_check true: maximum 1/K0 difference of a leftward partner from the seed.", false, true);
    registerDoubleOption_("charge:mono_position_min_ratio", "<q>", 0.2, "charge:mono_profile_check v2: the walk's leftmost step is rejected when its averagine-scaled position ratio PS(-j) / PS(-j+1) x lambda is below this.", false, true);
    registerDoubleOption_("charge:mono_position_im_box", "<1/K0>", 0.0, "charge:mono_profile_check v2: 1/K0 half-width around the seed of the isotope-position intensity sums. 0 = auto: 3.5 x the 90th percentile of the 1/K0 gap between the brightest tenth of the count-walk seeds and their nearest MS1 trace one isotope heavier, clamped to [0.01, 0.03] (0.018 and 0.022 on the 2019 HeLa py3 10 and 100 ng runs of PXD017703, whose isotope traces scatter widely in 1/K0). A value > 0 fixes the box. The count walk accepts partners up to charge:iso_im_tolerance away (auto: 0.01-0.05 per run), so a walk's own partner can lie outside the box and then counts as an empty position.", false, true);
    registerDoubleOption_("charge:mono_position_extend", "<q>", 0.6, "charge:mono_profile_check v2: when no step was rejected and the walk ended on a miss, move the monoisotope one isotope lighter (the most intense trace there, not claimed) if that position's ratio is >= this. 0 = off.", false, true);
    registerStringOption_("charge:mono_position_iterate", "<true/false>", "true", "charge:mono_profile_check v2: after a rejected leftmost step, test the new leftmost step too, and repeat while steps fail. false tests the leftmost step only.", false, true);
    setValidStrings_("charge:mono_position_iterate", {"true", "false"});
    registerDoubleOption_("charge:iso_im_tolerance", "<1/K0>", 0.0, "1/K0 tolerance for matching isotope partners (default 0 = auto since 1.6.0; 0.05 = the fixed value of 1.5.0). "
                          "auto measures it per run before the isotope walk: a stride sample of the MS1 traces is probed at each charge's +1 isotope position (the partner m/z tolerance, gate:delta_rt) "
                          "and at two control positions 50 ppm either side, where no isotope of the probe at that charge can sit; the nearest trace's 1/K0 gap is recorded at each. "
                          "The tolerance is the one that admits the most isotope-position traces net of twice the control traces (true partners admitted minus chance partners admitted), "
                          "on a 0.001 grid, clamped to [0.01, 0.05]: about 0.01 on dense short-gradient runs, up to 0.05 on sparse low-load runs. A value > 0 fixes it. Logged as [isotol], stamped spx:iso_im_tolerance. CHANGES OUTPUT.", false);
    registerStringOption_("charge:iso_partner", "<first/closest>", "first", "Which candidate inside the isotope-partner box (m/z tolerance, gate:delta_rt, charge:iso_im_tolerance) the walk takes: "
                          "'first' in m/z order (as in 1.5.0), or 'closest' to the seed in apex RT and 1/K0, each scaled by its tolerance (ties to the first in m/z). Screened 2026-10-01: closest lost 1.6-5.0% MSFragger peptides on two 10 ng HeLa runs and wrote 14-15% more spectra on Evosep for about the same peptides; not recommended. CHANGES OUTPUT.", false, true);
    setValidStrings_("charge:iso_partner", {"first", "closest"});

    registerTOPPSubsection_("dnoise", "MS1 denoising of the raw frames before peak picking (dnoise v0.1.0's MS1 path; Bruker .d only)");
    registerStringOption_("dnoise:ms1", "<true/false>", "true",
                          "Denoise every MS1 frame on its raw points before picking: a port of dnoise v0.1.0's default MS1 path, bit-identical up to AccumulationTime 100 ms "
                          "(Garrett, Diedrich & Yates III, bioRxiv 2026.08.27.747603; MIT, LICENSES/dnoise-MIT.txt). Bruker .d only: other inputs "
                          "pass through, logged and stamped spx:dnoise_ms1. Above 100 ms several raw counts share one intensity and the nearest is taken. A run whose raw points do not recover is refused, as is "
                          "trace:native_ms1_neighbors > 0. CHANGES OUTPUT (dataset D -4.3% Sage peptides, TNBC 009 +1.5%) and lowers peak memory "
                          "13-21% at one cell per tile: CHANGELOG, 2026-09-10 'MS1 denoising is on by default'. false = the output before 2026-09-10.", false);
    setValidStrings_("dnoise:ms1", {"true", "false"});
    registerIntOption_("dnoise:mz_half_width", "<bins>", 3, "Streak filter: TOF bins on each side of a point's own bin summed per scan (dnoise --mz-half-width).", false, true);
    registerIntOption_("dnoise:min_feature_length", "<scans>", 5, "Streak filter: OCCUPIED scans a run needs; bridged empty scans do not count (dnoise --min-feature-length).", false, true);
    registerIntOption_("dnoise:max_internal_gap", "<scans>", 2, "Streak filter: empty scans a run may bridge (dnoise --max-internal-gap).", false, true);
    registerIntOption_("dnoise:iterations", "<n>", 2, "Streak filter passes, each over the previous pass's survivors; 0 skips the streak filter, and the halo filter and the window gate still run unless dnoise:halo / dnoise:dia_ms1_window are false (dnoise --iterations).", false, true);
    registerIntOption_("dnoise:min_window_intensity", "<counts>", 0, "Streak filter: raw window sum a scan needs to count as occupied (dnoise --min-window-intensity).", false, true);
    registerIntOption_("dnoise:min_feature_intensity", "<counts>", 0, "Streak filter: raw sum a run needs over its span (dnoise --min-feature-intensity).", false, true);
    registerStringOption_("dnoise:halo", "<true/false>", "true", "Halo filter, once, on the streak survivors: keep a point at >= halo_peak_fraction of the most intense survivor in its box outside its own TOF bin (false = dnoise --no-halo).", false, true);
    setValidStrings_("dnoise:halo", {"true", "false"});
    registerDoubleOption_("dnoise:halo_peak_fraction", "<f>", 0.15, "Halo filter: the fraction (dnoise --halo-peak-fraction).", false, true);
    registerIntOption_("dnoise:halo_mz_idx_half_width", "<bins>", 80, "Halo filter: the box's half-width in TOF bins (dnoise --halo-mz-idx-half-width).", false, true);
    registerIntOption_("dnoise:halo_scan_half_width", "<scans>", 2, "Halo filter: the box's half-width in scans (dnoise --halo-scan-half-width).", false, true);
    registerStringOption_("dnoise:dia_ms1_window", "<true/false>", "true", "Keep only MS1 points inside one of the run's diaPASEF isolation windows (every DiaFrameMsMsWindows row), padded, in dnoise's own two-point calibration (false = dnoise --no-dia-ms1-window).", false, true);
    setValidStrings_("dnoise:dia_ms1_window", {"true", "false"});
    registerDoubleOption_("dnoise:dia_ms1_mz_pad", "<Th>", 5.0, "Isolation-window gate: m/z added on both sides of every window (dnoise --dia-ms1-mz-pad).", false, true);
    registerDoubleOption_("dnoise:dia_ms1_im_pad", "<1/K0>", 0.05, "Isolation-window gate: 1/K0 added on both sides of every window (dnoise --dia-ms1-im-pad).", false, true);

    // Bounds, not documentation: without these a negative count wraps through `Size` (min_fragments
    // -1 becomes ~1.8e19 and nothing is ever emitted) and a negative tolerance silently gates
    // everything out. Both fail as an empty-but-valid result, the worst failure mode here.
    setMinInt_("assembly:min_fragments", 1);
    setMinInt_("assembly:max_fragments", 1);
    setMinInt_("assembly:default_charge", 0);
    setMinFloat_("assembly:twin_im_tolerance", -1.0);   // any negative value = off
    setMinInt_("assembly:twin_rt_frames", -1);          // -1 = off (same frame), 0 = auto
    setValidStrings_("assembly:twin_rt_match", {"ext", "bin", "ppm"});
    setMinInt_("max_charge", 1);
    setMinInt_("perf:trace_bands", 0);      // 0 = auto (documented)
    setMinFloat_("gate:delta_im", 0.0);
    setMinFloat_("gate:delta_rt", 0.0);
    setMinFloat_("gate:min_correlation", -1.0);
    setMaxFloat_("gate:min_correlation", 1.0);
    setMinFloat_("charge:iso_im_tolerance", 0.0);
    setMinFloat_("charge:mono_profile_min_corr", -1.0);
    setMaxFloat_("charge:mono_profile_min_corr", 1.0);
    setMinFloat_("charge:mono_profile_max_dim", 0.0);
    setMinFloat_("charge:mono_position_min_ratio", 0.0);
    setMinFloat_("charge:mono_position_im_box", 0.0);
    setMinFloat_("trace:min_length_fwhm", 0.0);   // [minlen] 0 = the 2-point floor
    setMinFloat_("charge:mono_position_extend", 0.0);
    setMinFloat_("trace:ms1_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_noise_threshold_int", 0.0);
    setMinFloat_("trace:noise_threshold_int", 0.0);
    setMinInt_("perf:ms1_trace_bands", 0);      // 0 is read as 1 (off): MS1 has no auto
    setMinFloat_("perf:mem_fraction", 0.05);
    setMaxFloat_("perf:mem_fraction", 0.95);
    // [dnoise] dnoise's fields are unsigned, and a negative fraction or pad is not one
    for (const char* o : {"dnoise:mz_half_width", "dnoise:min_feature_length", "dnoise:max_internal_gap", "dnoise:iterations",
                          "dnoise:min_window_intensity", "dnoise:min_feature_intensity", "dnoise:halo_mz_idx_half_width", "dnoise:halo_scan_half_width"})
      setMinInt_(o, 0);
    for (const char* o : {"dnoise:halo_peak_fraction", "dnoise:dia_ms1_mz_pad", "dnoise:dia_ms1_im_pad"}) setMinFloat_(o, 0.0);
  }

  /// IM-aware MassTraceDetection (banded) + valley splitting on @p map. The apex gate is snr * noise
  /// while membership needs only > noise, so MS1 and MS2 pass independent snr and min length.
  /// @p map is CLEARED as it is consumed, so it is not resident through the valley-splitting peak.
  vector<Trace> detectTraces_(PeakMap& map, TraceStore& store, double delta_im, double noise,
                              double snr, double min_len, double split_valleys, int bands,
                              double split_scan_time,   // the run's cycle time for the splitter, 0 = per trace
                              int in_cat = -1,          // [ledger] the caller's category for `map`, -1 = uncharged
                              int part_cat = LC_BAND,   // [ledger] where this call's band partition belongs
                              Ms1Prune* ms1_prune = nullptr)   // [ms1-prune] the MS1 call only: the pre-prune edge sample
  {
    MassTraceDetection mtd;
    mtd.setLogType(ProgressLogger::NONE); // quiet + thread-safe (called from the parallel window loop)
    Param p = mtd.getParameters();
    p.setValue("mass_error_ppm", getDoubleOption_("trace:mass_error_ppm"));
    p.setValue("noise_threshold_int", noise);
    p.setValue("chrom_peak_snr", snr);
    p.setValue("ion_mobility_tolerance", delta_im);        // IM-aware tracing [C-1]
    p.setValue("min_trace_length", min_len);
    mtd.setParameters(p);

    vector<MassTrace> mts;
    // [outer] The traces as consecutive segments (the bands, then the valley-split chunks): their concatenation is the
    // sequence the serial gather and the serial re-concatenation used to build, and every later step walks it in that
    // order, so neither copy is needed.
    vector<vector<MassTrace>> seg;
    // [tstage] Sub-stage seconds of this phase, logged as one line.
    double tt_edge = 0, tt_dist = 0, tt_mtd = 0, tt_gather = 0, tt_split = 0;
    auto tt_now = []() { return std::chrono::steady_clock::now(); };
    auto tt_s = [](auto a, auto b) { return std::chrono::duration<double>(b - a).count(); };
    if (bands > 1)
    {
      // [bands] Partition the peaks by m/z and trace the bands concurrently. A trace cannot span more
      // than mass_error_ppm, so a halo of kBandHaloPpmMul tolerances makes the partition exact; a trace
      // is kept only by the band whose CORE holds its centroid.
      const double ppm = getDoubleOption_("trace:mass_error_ppm");
      // Band edges from m/z QUANTILES, not uniform m/z: peak density varies by orders of
      // magnitude across a window, and uniform edges would leave one band holding most of the work.
      auto _te0 = tt_now();
      vector<double> sample;
      if (ms1_prune && ms1_prune->on) sample.swap(ms1_prune->sample);   // [ms1-prune] taken before the prune; the pruned map's own would move the edges
      else
        for (const auto& s : map)
          for (Size i = 0; i < s.size(); i += 97) sample.push_back(s[i].getMZ());
      sort(sample.begin(), sample.end());
      vector<double> edge(bands + 1);
      if (sample.empty()) { edge.assign(bands + 1, 0.0); }
      else
        for (int b = 0; b <= bands; ++b)
          edge[b] = sample[std::min(sample.size() - 1, (size_t)((double)b / bands * sample.size()))];
      edge.front() = 0.0; edge.back() = std::numeric_limits<double>::max();
      if (ms1_prune && detOn_())   // [det] condition (i) of the MS1 prune, gated on its own
      { vector<uint64_t> eb; for (double e : edge) pushBits_(eb, e);
        writeLogInfo_("[ms1-edges] n=" + String(sample.size()) + " bands=" + String(bands) + " fnv=" + String(detDigest_(std::move(eb)))); }
      vector<double> ().swap(sample);

      // Partition ONCE (peaks are copied into exactly one core band plus any halo it falls in),
      // then release the source map so the halo duplicates are the only overhead.
      tt_edge = tt_s(_te0, tt_now());
      auto _td0 = tt_now();
      // [dist-par] Distribute chunks of spectra in parallel; chunks concatenate in chunk order, so every
      // band gets the same spectra in the same order (byte-identical), and each source spectrum is
      // released as it is copied (transient ~1x the map). docs/BASELINE.md, "`perf:ms1_trace_bands`:
      // the old sweep was measuring the serial distribute".
      const int n_thr = omp_get_max_threads();
      const int nchunk = std::max(1, std::min((int)map.size(), 4 * n_thr));
      vector<vector<PeakMap>> csub((size_t)nchunk, vector<PeakMap>(bands));
      std::exception_ptr dist_err;   // [task-err]
      #pragma omp taskloop grainsize(1) default(shared)
      for (int ci = 0; ci < nchunk; ++ci)
      { try {
      const size_t sp0 = (size_t)ci * map.size() / nchunk, sp1 = (size_t)(ci + 1) * map.size() / nchunk;
      for (size_t spi = sp0; spi < sp1; ++spi)
      {
        auto& s = map[spi];        // non-const: each spectrum is released once distributed
        const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
        // isSorted() is a full traversal of the spectrum, and the spectrum does not change between
        // bands -- evaluating it inside the band loop walked every peak once per band.
        const bool src_sorted = s.isSorted();
        for (int b = 0; b < bands; ++b)
        {
          // Size each halo from ITS OWN edge: one h derived from edge[b] under-sizes the upper
          // side by the ppm growth across the band.
          const double lo = edge[b] - edge[b] * ppm * 1e-6 * kBandHaloPpmMul;
          const double hi = edge[b + 1] + edge[b + 1] * ppm * 1e-6 * kBandHaloPpmMul;
          MSSpectrum t; t.setRT(s.getRT()); t.setMSLevel(1);
          OpenMS::DataArrays::FloatDataArray ia;
          ia.setName(Constants::UserParam::ION_MOBILITY);
          // [R3] sorted spectra: the band's range is two MZBegin lookups; the scan is the unsorted fallback.
          Size i0 = 0, i1 = s.size();
          if (src_sorted)
          {
            i0 = (Size)(s.MZBegin(lo) - s.begin());
            i1 = (Size)(s.MZBegin(hi) - s.begin());
            t.reserve(i1 - i0); ia.reserve(i1 - i0);
            for (Size i = i0; i < i1; ++i) { t.push_back(s[i]); if (ima && i < ima->size()) ia.push_back((*ima)[i]); }
          }
          else
            for (Size i = 0; i < s.size(); ++i)
              if (s[i].getMZ() >= lo && s[i].getMZ() < hi)
              { t.push_back(s[i]); if (ima && i < ima->size()) ia.push_back((*ima)[i]); }
          if (!t.empty()) { t.getFloatDataArrays().push_back(std::move(ia)); csub[ci][b].addSpectrum(std::move(t)); }
        }
        s.clear(true);   // distributed: release it now, so the source and the copies do not both
      }                  // stand at full size (the transient was ~2x the map)
      } catch (...) { taskErr_(dist_err); } }
      if (dist_err) std::rethrow_exception(dist_err);
      map.clear(true);
      vector<PeakMap> sub(bands);
      for (int b = 0; b < bands; ++b)
      {
        size_t nsp = 0; for (int ci = 0; ci < nchunk; ++ci) nsp += csub[ci][b].size();
        sub[b].reserveSpaceSpectra(nsp);
        for (int ci = 0; ci < nchunk; ++ci)
        { for (auto& sp : csub[ci][b]) sub[b].addSpectrum(std::move(sp)); csub[ci][b].clear(true); }
      }
      if (ms1_prune && detOn_())   // [det] which (frame, band) spectra exist: condition (ii) of the MS1 prune
      { String cnt; for (int b = 0; b < bands; ++b) cnt += (b ? "," : "") + String(sub[b].size());
        writeLogInfo_("[det] MS1 band spectra " + cnt); }
      // [ledger] the partition is a COPY of the whole map (plus a 1-5% halo) and lives from here
      // until each band is detected -- charged as ms1_map until now, which named the wrong
      // structure for the largest allocation of the MS1 phase. Hand the charge over at the swap.
      vector<long long> band_b((size_t)bands, 0);
      { long long tb = 0;
        for (int b = 0; b < bands; ++b) { band_b[(size_t)b] = ledMapBytes_(sub[b]); tb += band_b[(size_t)b]; }
        ledAdd_(part_cat, tb); if (in_cat >= 0) ledSet_(in_cat, 0); }
      vector<vector<PeakMap>>().swap(csub);
      tt_dist = tt_s(_td0, tt_now());
      // [ms1-trim] DIASPEXTRACT_MS1_TRIM=1 (MS1 call only): return the released MS1 map's free pages before
      // detection. Measured: docs/BASELINE.md, "TNBC 009: the pass-1 levers take 6.3 GB off".
      static const bool ms1_trim = []{ const char* e = std::getenv("DIASPEXTRACT_MS1_TRIM"); return e && *e && *e != '0'; }();
      if (ms1_trim && in_cat >= 0)
      {
#ifdef __GLIBC__
        const double _tt = phase_clock_(); const long _r0 = rss_mb_();
        malloc_trim(0);
        writeLogInfo_("[trim] MS1 band distribution -> detection: " + String(_r0) + " -> " + String(rss_mb_()) + " MB RSS in "
                      + String(phase_clock_() - _tt) + " s");
#endif
      }
      auto _tm0 = tt_now();
      vector<vector<MassTrace>> per(bands);
      std::exception_ptr band_err;   // [task-err]
      #pragma omp taskloop grainsize(1) default(shared)
      for (int b = 0; b < bands; ++b)
      { try {
        taskFault_("ms1_band");
        MassTraceDetection m2; m2.setLogType(ProgressLogger::NONE); m2.setParameters(p);
        vector<MassTrace> t;
        sub[b].sortSpectra();
        // NOTE: MassTraceDetection builds its OWN copy of this band internally; that copy is the
        // one structure of the MS1 phase the tool cannot charge without patching OpenMS.
        if (sub[b].size() < 3) { sub[b].clear(true); ledAdd_(part_cat, -band_b[(size_t)b]); continue; }
        m2.run(sub[b], t);
        sub[b].clear(true); ledAdd_(part_cat, -band_b[(size_t)b]);
        // keep only traces whose centroid lies in this band's CORE -> no duplicates from the halo
        for (auto& mt : t)
        {
          const double c = mt.getCentroidMZ();
          if (c >= edge[b] && c < edge[b + 1]) per[b].push_back(std::move(mt));
        }
      } catch (...) { taskErr_(band_err); } }
      if (band_err) std::rethrow_exception(band_err);
      tt_mtd = tt_s(_tm0, tt_now());
      auto _tg0 = tt_now();
      seg = std::move(per);   // [outer] the bands, in band order, are the segments: no serial gather
      tt_gather = tt_s(_tg0, tt_now());
    }
    else
    {
      if (map.size() >= 3) mtd.run(map, mts);   // MTD aborts below 3 spectra; a tiny input is not an error
      map.clear(true); // [compact] dead from here; ~1.1 GB/window freed before the EPD peak
      if (in_cat >= 0) ledSet_(in_cat, 0);
      seg.resize(1); seg[0].swap(mts);
    }
    size_t n_tr = 0; for (const auto& s : seg) n_tr += s.size();

    // [way-2] Split mass traces at chromatographic valleys: MassTraceDetection never splits at a local
    // minimum, so co-eluting neighbours within the m/z tolerance merge. width_filtering is forced off
    // (a deletion filter). Per-trace and deterministic, so it runs as a balanced parallel loop (4 chunks
    // per band); both consumers re-sort (MS2 after detectTraces_, MS1 before inferPrecursors_).
    auto _ts0 = tt_now();
    if (split_valleys > 0.0 && n_tr > 1)
    {
      const int nchunk = std::max(1, bands * 4);      // 4 chunks per band -> dynamic load balance
      const size_t n_in = n_tr;
      vector<size_t> at(seg.size() + 1, 0);           // [outer] segment s holds positions [at[s], at[s+1]) of the sequence
      for (size_t s = 0; s < seg.size(); ++s) at[s + 1] = at[s] + seg[s].size();
      vector<vector<MassTrace>> parts(nchunk);
      std::exception_ptr split_err;   // [task-err]
      #pragma omp taskloop grainsize(1) default(shared)
      for (int ci = 0; ci < nchunk; ++ci)
      { try {
        const size_t lo = (size_t)ci * n_in / nchunk, hi = (size_t)(ci + 1) * n_in / nchunk;
        if (lo >= hi) continue;
        // [outer] positions [lo, hi) of the segment sequence, moved out segment by segment: the chunk the gathered
        // vector gave (distinct chunks move distinct elements)
        vector<MassTrace> chunk; chunk.reserve(hi - lo);
        size_t s = (size_t)(std::upper_bound(at.begin(), at.end(), lo) - at.begin()) - 1;
        for (size_t g = lo; g < hi; ++g) { while (g >= at[s + 1]) ++s; chunk.push_back(std::move(seg[s][g - at[s]])); }
        ElutionPeakDetection epd; setSplitter(epd, split_valleys, split_scan_time);
        parts[ci].reserve(chunk.size());
        epd.detectPeaks(chunk, parts[ci]);
      } catch (...) { taskErr_(split_err); } }
      if (split_err) std::rethrow_exception(split_err);
      seg.swap(parts);   // [outer] the chunks' outputs, in chunk order, are the segments now: no serial re-concatenation
      n_tr = 0; for (const auto& s : seg) n_tr += s.size();
    }
    tt_split = tt_s(_ts0, tt_now());
    writeLogInfo_("[tstage] " + String(bands) + " bands -> " + String(n_tr) + " traces: edges "
                  + String(tt_edge) + " s | distribute " + String(tt_dist) + " s | detect "
                  + String(tt_mtd) + " s | gather " + String(tt_gather) + " s (serial) | valley-split "
                  + String(tt_split) + " s");
    auto _tc0 = tt_now();
    vector<Trace> out;
    // [outer] The MS1 call converts in parallel: pass 1 sizes every span (makeSpan's frame range), the offsets are the
    // running sum from the arena's end -- exactly where the serial appends put them -- and pass 2 fills each span at
    // its offset with toTrace's own statements. Every trace writes only its own span and its own slot, so the arena and
    // the trace list are byte-identical to the serial loop's. Contiguous ranges per task keep the walk sequential. The
    // two earlier parallel versions kept the arena append serial (docs/BASELINE.md, "Parallelising the MS1 trace
    // conversion"); DIASPEXTRACT_MS1_CONVERT_SERIAL=1 runs the serial loop below.
    static const bool conv_serial = []{ const char* e = std::getenv("DIASPEXTRACT_MS1_CONVERT_SERIAL"); return e && *e && *e != '0'; }();
    bool conv_par = ms1_prune && !conv_serial && n_tr >= 4096 && omp_in_parallel() && omp_get_num_threads() > 1;
    const char* conv_why = "";
    vector<size_t> at, off;   // segment s holds positions [at[s], at[s+1]) of the sequence; off: span length, then offset
    auto segOf = [&at](size_t g) { return (size_t)(std::upper_bound(at.begin(), at.end(), g) - at.begin()) - 1; };
    const int ntask = conv_par ? std::max(1, std::min(4 * omp_get_num_threads(), (int)(n_tr / 1024))) : 1;   // bounded: no task-count cliff
    std::exception_ptr conv_err;   // [task-err]
    if (conv_par)
    {
      at.assign(seg.size() + 1, 0);
      for (size_t s = 0; s < seg.size(); ++s) at[s + 1] = at[s] + seg[s].size();
      off.assign(n_tr + 1, 0);
      #pragma omp taskloop num_tasks(ntask) default(shared)
      for (int c = 0; c < ntask; ++c)
      { try {
        const size_t lo = (size_t)c * n_tr / ntask, hi = (size_t)(c + 1) * n_tr / ntask;
        if (lo >= hi) continue;
        size_t s = segOf(lo);
        for (size_t g = lo; g < hi; ++g) { while (g >= at[s + 1]) ++s; off[g + 1] = spanLength(seg[s][g - at[s]], store); }
      } catch (...) { taskErr_(conv_err); } }
      if (conv_err) std::rethrow_exception(conv_err);
      off[0] = store.inten.size();
      for (size_t g = 0; g < n_tr; ++g) off[g + 1] += off[g];   // off[g] = the arena size before trace g's append
      // [review] Trace::off is 32-bit and makeSpan narrows the arena size without a check. Past 2^32 - 1 floats two
      // traces' narrowed offsets would alias and their tasks would race on the same span; the serial loop below keeps
      // the old (deterministic) result there instead.
      if (off[n_tr] > (size_t)std::numeric_limits<uint32_t>::max()) { conv_par = false; conv_why = ", the arena passes 2^32 floats"; }
    }
    if (conv_par)
    {
      store.inten.resize(off[n_tr], 0.0f);
      out.resize(n_tr);
      #pragma omp taskloop num_tasks(ntask) default(shared)
      for (int c = 0; c < ntask; ++c)
      { try {
        const size_t lo = (size_t)c * n_tr / ntask, hi = (size_t)(c + 1) * n_tr / ntask;
        if (lo >= hi) continue;
        vector<pair<uint32_t, float>> pts;
        size_t s = segOf(lo);
        for (size_t g = lo; g < hi; ++g)
        {
          while (g >= at[s + 1]) ++s;
          MassTrace& mt = seg[s][g - at[s]];
          out[g] = toTraceAt(mt, store, off[g], pts);
          mt = MassTrace();   // the payload goes now, as in the serial loop
        }
      } catch (...) { taskErr_(conv_err); } }
      if (conv_err) std::rethrow_exception(conv_err);
      vector<vector<MassTrace>>().swap(seg);
    }
    else
    {
      out.reserve(n_tr);
      // Each MassTrace is released as soon as it is converted.
      for (auto& sv : seg)   // [outer] the segments in order: the sequence the gathered vector held
      {
        for (size_t i = 0; i < sv.size(); ++i) { out.push_back(toTrace(sv[i], store)); sv[i] = MassTrace(); }   // move-assign an empty one: the payload goes now, not at return
        vector<MassTrace>().swap(sv);
      }
    }
    if (ms1_prune) writeLogInfo_("[outer] MS1 trace conversion " + String(tt_s(_tc0, tt_now())) + " s (" + (conv_par ? "parallel" : "serial") + conv_why + ")");
    return out;
  }

  /// [route-1] Averagine isotope envelope, Poisson approximation: for neutral mass M the expected
  /// number of 13C is lambda ~ 0.000594*M, giving relative intensities [1, l, l^2/2, l^3/6, ...].
  /// This is what lets us pick the MONOISOTOPE: the mono is NOT always the most intense peak
  /// (above ~1800 Da M+1 exceeds M), so intensity-ranking alone mis-assigns it.
  static void averagineRatios(double neutral_mass, int n, vector<double>& out)
  {
    const double lam = std::max(0.10, 0.000594 * neutral_mass);
    out.assign(n, 0.0);
    double term = 1.0;
    for (int k = 0; k < n; ++k) { out[k] = term; term *= lam / (double)(k + 1); }
    double mx = 0.0; for (double v : out) mx = std::max(mx, v);
    if (mx > 0) for (double& v : out) v /= mx;
  }

  static double cosineSim(const vector<double>& a, const vector<double>& b)
  {
    double dot = 0, na = 0, nb = 0;
    for (size_t i = 0; i < a.size() && i < b.size(); ++i) { dot += a[i]*b[i]; na += a[i]*a[i]; nb += b[i]*b[i]; }
    return (na > 0 && nb > 0) ? dot / (sqrt(na) * sqrt(nb)) : 0.0;
  }

  /// [route-1] Pearson of two MS1 XICs over the frames where BOTH have a real point (mean-centred on
  /// that overlap, not zero-padded); -2 below 3 shared points or at zero variance.
  static double xicCorr(const TraceStore& s, const Trace& a, const Trace& b)
  {
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0; int m = 0;
    const uint32_t lo = std::max(a.frame0, b.frame0);
    const uint32_t hi = std::min(a.frame0 + a.len, b.frame0 + b.len);
    for (uint32_t f = lo; f < hi; ++f)
    {
      const double x = xv(s, a, f - a.frame0), y = xv(s, b, f - b.frame0);
      if (x > 0.0 && y > 0.0) { sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; ++m; }
    }
    if (m < 3) return -2.0;
    const double cov = sxy - sx * sy / m, vx = sxx - sx * sx / m, vy = syy - sy * sy / m;
    return (vx > 0 && vy > 0) ? cov / std::sqrt(vx * vy) : -2.0;
  }

  /// [isotol] charge:iso_im_tolerance 0 (auto): the run's isotope-partner 1/K0 tolerance, measured on its own MS1 traces
  /// before the greedy walk. `ms1` is the content-sorted list (m/z ascending), so the sample and every search are the same
  /// at any thread count and tiling. Every ceil(N / 200000)-th trace is a probe. At each charge z <= max_charge it is
  /// probed at its +1 isotope position (m/z + 1.0033548 / z) and at two controls 50 ppm either side, where no isotope of the
  /// probe at that charge can sit (another charge's can: near m/z 1,000 a 4+ ion's M+1 lies 0.2 ppm from the 5+ control,
  /// so a true partner there counts as a chance one), in the walk's partner window (probe m/z x mass_ppm, gate:delta_rt), and the nearest other trace's
  /// |d 1/K0| is kept. A tolerance t then admits P_iso(t) traces per isotope probe, of which about P_ctl(t) (the mean of
  /// the two controls) are chance traces: t* maximises P_iso(t) - 2 P_ctl(t), the true partners admitted minus the chance
  /// partners admitted, on a 0.001 grid over [0.002, 0.050] (ties to the wider), then clamped to [0.01, 0.05] (0.01: the
  /// tightest value screened without loss on dense runs; 0.05: the fixed value of 1.5.0). `info` gets the stamp text.
  double autoIsoImTol_(const vector<Trace>& ms1, const TraceStore& ms1st, int max_charge, double delta_rt, double mass_ppm, String& info)
  {
    const double ISO = 1.0033548, CTL = 50e-6, LO = 0.01, HI = 0.05;
    const int NG = 49;                                    // grid t_x = 0.002 + 0.001 x, x = 0..48
    auto tg = [](int x) { return 0.002 + 0.001 * x; };
    const size_t N = ms1.size();
    if (N < 2 || max_charge < 1)
    {
      info = "auto 0.05 (too few MS1 traces to measure)";
      writeLogInfo_("[isotol] auto: " + String(N) + " MS1 traces, too few to measure; 1/K0 tolerance 0.05");
      return HI;
    }
    vector<double> mz(N), rt(N), im(N);
    for (size_t i = 0; i < N; ++i) { mz[i] = ms1[i].mz; rt[i] = rtOf(ms1st, ms1[i]); im[i] = ms1[i].im; }
    if (!std::is_sorted(mz.begin(), mz.end()))
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[isotol] the MS1 traces are not sorted by m/z");
    const size_t stride = std::max<size_t>(1, (N + 199999) / 200000), ns = (N + stride - 1) / stride;
    const int Z = max_charge;
    vector<int8_t> gx(ns * (size_t)Z * 3, -1);            // per (probe, z, kind): first grid index x with gap <= t_x, -1 none
    #pragma omp parallel for schedule(dynamic, 256)
    for (long q = 0; q < (long)ns; ++q)
    {
      const size_t s = (size_t)q * stride;
      const double tol = mz[s] * mass_ppm * 1e-6;
      for (int z = 1; z <= Z; ++z)
        for (int c = 0; c < 3; ++c)
        {
          const double base = mz[s] + ISO / z;
          const double target = c == 0 ? base : c == 1 ? base * (1.0 + CTL) : base * (1.0 - CTL);
          const size_t klo = (size_t)(std::lower_bound(mz.begin(), mz.end(), target - tol) - mz.begin());
          const size_t khi = (size_t)(std::upper_bound(mz.begin(), mz.end(), target + tol) - mz.begin());
          double g = 1.0;
          for (size_t k = klo; k < khi; ++k)
          {
            if (k == s || fabs(rt[k] - rt[s]) > delta_rt) continue;
            const double d = fabs(im[k] - im[s]);
            if (d < g) g = d;
          }
          int8_t b = -1;
          for (int x = 0; x < NG; ++x) if (g <= tg(x)) { b = (int8_t)x; break; }
          gx[((size_t)q * (size_t)Z + (size_t)(z - 1)) * 3 + (size_t)c] = b;
        }
    }
    vector<long long> cum(3 * NG, 0);                     // kind c, grid x: probes whose nearest gap is <= t_x
    for (size_t i = 0; i < gx.size(); ++i) if (gx[i] >= 0) ++cum[(i % 3) * NG + (size_t)gx[i]];
    for (int c = 0; c < 3; ++c) for (int x = 1; x < NG; ++x) cum[c * NG + x] += cum[c * NG + x - 1];
    int best = 0; long long bj = 0;
    for (int x = 0; x < NG; ++x)
    {
      const long long J = cum[x] - cum[NG + x] - cum[2 * NG + x];
      if (x == 0 || J >= bj) { best = x; bj = J; }
    }
    const double tstar = tg(best), tol = std::min(HI, std::max(LO, tstar));
    const int xc = (int)std::lround((tol - 0.002) / 0.001);
    const double np = (double)ns * Z;
    auto rate = [&](int c, int x) { return (double)cum[c * NG + x] / np; };
    char b[400];
    std::snprintf(b, sizeof b, "auto %.3f (t* %.3f of isotope minus 2 x control traces, clamped to [0.01, 0.05]; per probe within %.3f: isotope %.3f, control %.3f; %zu probes x %d charges)",
                  tol, tstar, tol, rate(0, xc), 0.5 * (rate(1, xc) + rate(2, xc)), ns, Z);
    info = String(b);
    std::snprintf(b, sizeof b, "[isotol] auto: 1/K0 tolerance %.3f for isotope partners (t* %.3f, clamped to [0.01, 0.05]). %zu probes (stride %zu over %zu MS1 traces) x %d charges: "
                  "the nearest trace lies within the tolerance at the +1 isotope position for %.3f of the probes, at the +-50 ppm controls for %.3f (within 0.05: %.3f, %.3f)",
                  tol, tstar, ns, stride, N, Z, rate(0, xc), 0.5 * (rate(1, xc) + rate(2, xc)), rate(0, NG - 1), 0.5 * (rate(1, NG - 1) + rate(2, NG - 1)));
    writeLogInfo_(String(b));
    return tol;
  }

  /// [outer] The m/z order inferPrecursors_ and monoPositionPass_ share: the indices of @p ms1 by std::sort on m/z alone.
  /// The MS1 traces arrive content-sorted (m/z is the first key), so unless two of them share an m/z the m/z values
  /// strictly increase in container order, the identity is then the only sorted order, and it is taken without sorting.
  /// With any equal (or unordered) neighbour -- a 60-min run has 677,976, a 2-h HeLa run 19,200: quantised m/z -- the order of equal
  /// m/z is std::sort's own, so it is std::sort's sequence, computed by stdSortPar_. Returns the neighbours that are not
  /// strictly increasing (0 = the identity was taken).
  static size_t mzOrder_(const vector<Trace>& ms1, vector<size_t>& by_mz)
  {
    const size_t N = ms1.size();
    by_mz.resize(N);
    #pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)N; ++i) by_mz[(size_t)i] = (size_t)i;
    size_t nonincr = 0;
    #pragma omp parallel for schedule(static) reduction(+ : nonincr)
    for (long i = 1; i < (long)N; ++i) if (!(ms1[(size_t)i - 1].mz < ms1[(size_t)i].mz)) ++nonincr;
    if (nonincr) stdSortPar_(by_mz, [&ms1](size_t a, size_t b) { return ms1[a].mz < ms1[b].mz; });
    return nonincr;
  }

  /// [route-1] Isotope/charge inference on MS1 traces. Scores each (charge, monoisotope) hypothesis
  /// by averagine-shape agreement x isotope-XIC co-elution, instead of counting partners.
  vector<Precursor_> inferPrecursors_(const vector<Trace>& ms1, const TraceStore& ms1st, int max_charge,
                                      double delta_rt, double iso_im_tol, double mass_ppm,
                                      bool envelope_scoring, vector<MonoV2In_>* v2_out = nullptr)
  {
    const double ISO = 1.0033548;
    const double PROTON_ = 1.00727646;
    const size_t N = ms1.size();

    // m/z-sorted index for O(log N) isotope-partner lookup (avoids O(N^2) over all MS1 traces).
    auto _pi_t0 = std::chrono::steady_clock::now();
    vector<size_t> by_mz;
    const size_t mz_nonincr = mzOrder_(ms1, by_mz);   // [outer]
    vector<double> sorted_mz(N);
    // [soa] The candidate scan's RT and IM, laid out in m/z order beside sorted_mz (sequential reads).
    vector<double> sorted_rt(N), sorted_im(N);
    #pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long)N; ++ii)
    { const size_t i = (size_t)ii; const Trace& t = ms1[by_mz[i]]; sorted_mz[i] = t.mz; sorted_rt[i] = rtOf(ms1st, t); sorted_im[i] = t.im; }
    // [mzidx] Bucket index over sorted_mz keyed on the top bits of the IEEE-754 pattern: monotone for
    // positive doubles, so the bucket width is relative (ppm-shaped) and the std search on the bucket's
    // sub-range returns the same position as on the whole array. docs/BASELINE.md, "Computed indices
    // instead of binary searches".
    const double mz_lo = N ? sorted_mz.front() : 0.0, mz_hi = N ? sorted_mz.back() : 1.0;
    unsigned kshift = 32;
    auto rawKey = [](double x, unsigned sh) { uint64_t b; std::memcpy(&b, &x, sizeof b); return b >> sh; };
    while (kshift < 52 && (rawKey(mz_hi, kshift) - rawKey(mz_lo, kshift)) > (1u << 23)) ++kshift;
    const uint64_t klo = N ? rawKey(mz_lo, kshift) : 0;
    const size_t kspan = N ? (size_t)(rawKey(mz_hi, kshift) - klo) : 0;
    auto bucketOf = [&](double x) -> size_t {
      if (!(x > mz_lo)) return 0;
      const uint64_t k = rawKey(x, kshift);
      return k >= klo + kspan ? kspan : (size_t)(k - klo);
    };
    vector<uint32_t> bstart(kspan + 2, 0);
    {
      // bstart[b] = the first i with bucketOf(sorted_mz[i]) >= b, else N. [outer] bucketOf is monotone over sorted_mz, so
      // the buckets (bucketOf(sorted_mz[i-1]), bucketOf(sorted_mz[i])] start at i: every b is written by exactly one i.
      #pragma omp parallel for schedule(static)
      for (long ii = 0; ii < (long)N; ++ii)
      {
        const size_t i = (size_t)ii, b_hi = std::min(bucketOf(sorted_mz[i]), kspan);
        for (size_t b = (i == 0 ? 0 : bucketOf(sorted_mz[i - 1]) + 1); b <= b_hi; ++b) bstart[b] = (uint32_t)i;
      }
      for (size_t b = (N ? std::min(bucketOf(sorted_mz[N - 1]), kspan) + 1 : 0); b <= kspan; ++b) bstart[b] = (uint32_t)N;
      bstart[kspan + 1] = (uint32_t)N;
    }
    auto lbmz = [&](double x) -> size_t {
      if (N == 0 || x <= mz_lo) return 0;
      if (x > mz_hi) return N;
      const size_t b = bucketOf(x);
      return (size_t)(std::lower_bound(sorted_mz.begin() + bstart[b], sorted_mz.begin() + bstart[b + 1], x) - sorted_mz.begin());
    };
    auto ubmz = [&](double x) -> size_t {
      if (N == 0 || x < mz_lo) return 0;
      if (x >= mz_hi) return N;
      const size_t b = bucketOf(x);
      return (size_t)(std::upper_bound(sorted_mz.begin() + bstart[b], sorted_mz.begin() + bstart[b + 1], x) - sorted_mz.begin());
    };
    auto _pi_t1 = std::chrono::steady_clock::now();

    // Apex intensities in one contiguous array: the sort and every isotope walk read them.
    vector<float> inten(N);
    vector<size_t> order(N);
    #pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long)N; ++ii) { const size_t i = (size_t)ii; inten[i] = (float)intensityOf(ms1st, ms1[i]); order[i] = i; }
    // [outer] a total order (the index breaks every tie; the coordinates are finite, non-finite traces were dropped before
    // the content sort), so the sorted sequence is unique and the parallel sort returns std::sort's
    parallelSortTotal_(order, [&](size_t a, size_t b) {
      if (inten[a] != inten[b]) return inten[a] > inten[b];
      if (ms1[a].mz != ms1[b].mz) return ms1[a].mz < ms1[b].mz;
      return a < b; // total order so greedy isotope assignment is deterministic
    });

    auto _pi_t2 = std::chrono::steady_clock::now();
    vector<bool> used(N, false);
    vector<Precursor_> out;
    // [spec] Ordered speculation, exact: the only mutable state a seed reads is `used`, which is monotone.
    // A batch of seeds is evaluated concurrently against a frozen `used`, recording the partners each
    // search returned; the commit walks the batch in the original order, skips a seed that is now used,
    // and re-evaluates one whose recorded partner has since been claimed. A search that returned -1 saw
    // only gate failures or used candidates, both stable under a monotone bitmap, so it needs no
    // dependency. docs/BASELINE.md, "Ordered speculation on the greedy claim loop".
    // mc: [monocheck] v1 (true): walks checked, monos moved, rejected by profile, by 1/K0, steps passed with < 3
    // shared frames (or a flat profile), partners released, moved monos left with no claimed partner, steps
    // accepted by a profile that could be judged. v2: walks tested, monos moved right, steps rejected, extensions
    // tested, monos moved left, -, monos moved back onto the seed, -. Summed at commit, not in speculation.
    // [monocheck-v2] v2: the count walk's result for monoPositionPass_ (MonoV2In_), which runs after the IM charge veto.
    struct SeedRes { vector<Precursor_> emit; vector<size_t> claims; vector<size_t> deps; std::array<size_t, 8> mc{}; MonoV2In_ v2; };
    auto evalSeed = [&](size_t seed, SeedRes& R)
    {
      R.emit.clear(); R.claims.clear(); R.deps.clear(); R.mc.fill(0); R.v2 = MonoV2In_();
      // do/while(0): the body below is the sequential loop body VERBATIM, and its seed-level
      // `continue` statements exit this wrapper exactly as they used to exit one loop iteration,
      // while `continue` inside the body's own inner loops keeps binding to those loops.
      do
      {
      if (used[seed]) continue;
      const Trace& s = ms1[seed];
      const double s_rt = rtOf(ms1st, s);     // [soa-rt] hoisted: the candidate scan reads it per probe
      const double tol = s.mz * mass_ppm * 1e-6;

      // Find an unused isotope partner near `target` m/z, co-localized in RT and IM (binary-search
      // the m/z index, scan only the ppm-window candidates); -1 if none.
      // [isotol] charge:iso_partner closest: the unused candidate nearest the seed in (dRT / delta_rt)^2 + (d1/K0 / iso_im_tol)^2,
      // ties to the first in m/z. Recording only the chosen partner as the dependency stays exact: `used` is monotone, so a
      // candidate claimed later that was not chosen cannot change the minimum while the chosen one stays unused.
      auto findPartner = [&](double target) -> long {
        size_t klo = lbmz(target - tol), khi = ubmz(target + tol);   // [mzidx] same positions, local search
        long best = -1; double bd = 0.0;
        for (size_t k = klo; k < khi; ++k)
        {
          if (fabs(sorted_rt[k] - s_rt) > delta_rt || fabs(sorted_im[k] - s.im) > iso_im_tol) continue;   // [soa]
          size_t j = by_mz[k];
          if (used[j] || j == seed) continue;
          if (!iso_closest_) { R.deps.push_back((size_t)j); return (long)j; }
          const double dr = delta_rt > 0.0 ? (sorted_rt[k] - s_rt) / delta_rt : 0.0;
          const double di = iso_im_tol > 0.0 ? (sorted_im[k] - s.im) / iso_im_tol : 0.0;
          const double d = dr * dr + di * di;
          if (best < 0 || d < bd) { best = (long)j; bd = d; }
        }
        if (best >= 0) R.deps.push_back((size_t)best);
        return best;
      };

      // charge:scoring=count end to end, sharing nothing with the envelope path below (no gapped walk, no
      // lead offsets; mixing the two was a third, worse algorithm): contiguous walk from the seed with
      // break-on-first-miss, monoisotope = furthest-left partner reached, z by partner count, ties to lowest z.
      if (!envelope_scoring)
      {
        int best_z = 0, best_n = 1, best_nl = 0;
        bool best_gstop = false;                // [monocheck-v2] the winning walk stopped at the guard, not on a miss
        double best_mono = s.mz;
        vector<size_t> best_partners;
        for (int z = 1; z <= max_charge; ++z)
        {
          vector<size_t> partners;
          double mono = s.mz;
          double prev_int = inten[seed];
          vector<long> lefts;
          bool gstop = false;
          for (int k = 1; k <= 5; ++k)          // lighter partners define the monoisotope
          {
            long j = findPartner(s.mz - k * ISO / z);
            if (j < 0) break;                   // BREAK on first miss (contiguous only)
            // [mono-guard] findPartner never checks intensity, so a leftward step can latch onto noise one
            // isotope below. Under averagine the lighter of isotopes m, m+1 is (m+1)/lambda times the heavier,
            // minimised at m=0 -> 1/lambda (see charge:mono_averagine_guard).
            if (mono_guard_ > 0.0)
            {
              const double neutral = (s.mz - k * ISO / z) * z - z * PROTON_;
              const double lam = std::max(0.10, 0.000594 * neutral);
              if (inten[j] < mono_guard_ * prev_int / lam) { gstop = true; break; }
            }
            lefts.push_back(j);
            partners.push_back((size_t)j); mono = ms1[j].mz; prev_int = inten[j];
          }
          vector<long> rights;
          for (int k = 1; k <= 5; ++k)          // heavier partners add confidence
          {
            long j = findPartner(s.mz + k * ISO / z);
            if (j < 0) break;
            rights.push_back(j);
            partners.push_back((size_t)j);
          }
          // [mono-select] Pick the mono by averagine fit over the whole run (see charge:mono_averagine_select);
          // every found peak is still marked used, so only the reported mono m/z changes.
          if (mono_select_ && !lefts.empty())
          {
            const int nL = (int)lefts.size();
            vector<double> runI;                                  // ascending m/z: lefts, seed, rights
            for (int i = nL - 1; i >= 0; --i) runI.push_back(inten[lefts[i]]);
            runI.push_back(inten[seed]);
            for (long j : rights) runI.push_back(inten[j]);
            double best_cos = -1.0; int best_c = nL;              // default: seed is the mono
            vector<double> theo, obs;
            for (int c = 0; c <= nL; ++c)                         // mono cannot be right of the seed
            {
              const int len = std::min(4, (int)runI.size() - c);
              if (len < 2) break;                                 // need >=2 peaks to score a shape
              const double neutral = (s.mz - (double)(nL - c) * ISO / z) * z - z * PROTON_;
              averagineRatios(neutral, len, theo);
              obs.assign(runI.begin() + c, runI.begin() + c + len);
              double mx = 0.0; for (double v : obs) mx = std::max(mx, v);
              if (mx <= 0.0) continue;
              for (double& v : obs) v /= mx;
              // Weighted by evidence as in the envelope scorer, (n-1)/3 capped at 1: raw cosine of the
              // shorter vector a right-shifted candidate gets favours a too-heavy mono.
              const double cs = cosineSim(obs, theo) * std::min(1.0, (double)(len - 1) / 3.0);
              if (cs > best_cos) { best_cos = cs; best_c = c; }
            }
            mono = (best_c < nL) ? ms1[lefts[nL - 1 - best_c]].mz : s.mz;
          }
          const int n = 1 + (int)partners.size();
          if (n > best_n) { best_n = n; best_z = z; best_mono = mono; best_partners = partners; best_nl = (int)lefts.size(); best_gstop = gstop; }
        }
        // [monocheck] Re-walk the winning charge's leftward partners (best_partners[0, best_nl), nearest first):
        // each must share the seed's elution profile and mobility, else the walk stops there and the
        // rejected trace and those beyond it are released unclaimed. Charge and n_isotopes stay as counted.
        if (mono_check_ == 1 && best_nl > 0)
        {
          ++R.mc[0];
          int acc = 0; double acc_mz = s.mz;
          for (; acc < best_nl; ++acc)
          {
            const Trace& t = ms1[best_partners[acc]];
            if (fabs(t.im - s.im) > mono_check_max_dim_) { ++R.mc[3]; break; }
            const double r = xicCorr(ms1st, t, s);   // -2: < 3 shared frames (or flat), too little to judge
            if (r <= -2.0) ++R.mc[4];
            else if (r < mono_check_min_corr_) { ++R.mc[2]; break; }
            else ++R.mc[7];
            acc_mz = t.mz;
          }
          if (acc < best_nl)
          {
            R.mc[5] += (size_t)(best_nl - acc);
            if (acc_mz > best_mono) { best_mono = acc_mz; ++R.mc[1]; }   // mono_averagine_select may already sit further right
            best_partners.erase(best_partners.begin() + acc, best_partners.begin() + best_nl);
            if (best_partners.empty()) ++R.mc[6];
          }
        }
        // [monocheck-v2] Record what monoPositionPass_ needs (v2 reads only the immutable arrays and changes only
        // mono_mz, so deciding it after the loop gives the same result as deciding it here, lets an automatic im_box
        // use the run's own isotope scatter, and lets the decision follow the IM charge veto).
        if (mono_check_ == 2 && best_z > 0)
        {
          R.v2.z = best_z; R.v2.nl = best_nl; R.v2.mono = best_mono; R.v2.guard_stop = best_gstop;
          for (int i = 0; i < best_nl; ++i) R.v2.left_mz[(size_t)i] = ms1[best_partners[(size_t)i]].mz;
        }
        Precursor_ pc0;
        pc0.trace_idx = seed;
        pc0.rt = s_rt; pc0.im = s.im;
        pc0.mono_mz = best_mono;
        pc0.charge = best_z;                    // 0 => unknown, default_charge applied later
        pc0.n_isotopes = best_n;                // [Q4] envelope confidence for the mono call
        R.emit.push_back(pc0);
        for (size_t j : best_partners) R.claims.push_back(j);   // de-isotope
        continue;
      }

      struct Hyp { int z; double mono; double score; vector<size_t> partners; };
      vector<Hyp> hyps;
      vector<double> theo, obs(4), obsn(4);

      for (int z = 1; z <= max_charge; ++z)
      {
        // [route-1] GAPPED envelopes: no break-on-first-miss (a sub-threshold M+1 with a present
        // M+2 previously killed the whole walk). Try the seed AND up to 3 isotopes below it as the
        // monoisotope candidate; averagine decides which one really is the mono.
        for (int lead = 0; lead <= 3; ++lead)
        {
          const double mono_mz = s.mz - lead * ISO / z;
          long mono_idx = (lead == 0) ? (long)seed : findPartner(mono_mz);
          if (lead > 0 && mono_idx < 0) continue;
          const double neutral = mono_mz * z - z * PROTON_;
          if (neutral < 400.0 || neutral > 6000.0) continue;

          vector<size_t> members;
          vector<long> iso_idx(4, -1);
          obs.assign(4, 0.0);
          obs[0] = inten[mono_idx];
          iso_idx[0] = mono_idx;
          if ((size_t)mono_idx != seed) members.push_back((size_t)mono_idx);
          int found = 1;
          for (int k = 1; k <= 3; ++k)
          {
            long j = findPartner(mono_mz + k * ISO / z);
            if (j >= 0) { obs[k] = inten[j]; iso_idx[k] = j; ++found;
                          if ((size_t)j != seed) members.push_back((size_t)j); }
          }
          if (found < 2) continue;   // need >=2 isotopes before a shape can be scored

          averagineRatios(neutral, 4, theo);
          double mx = 0.0; for (double v : obs) mx = std::max(mx, v);
          if (mx <= 0) continue;
          for (int k = 0; k < 4; ++k) obsn[k] = obs[k] / mx;
          const double cos_sim = cosineSim(obsn, theo);

          double csum = 0.0; int cn = 0;
          for (int k = 1; k <= 3; ++k)
            if (iso_idx[k] >= 0) { csum += xicCorr(ms1st, ms1[iso_idx[0]], ms1[iso_idx[k]]); ++cn; }
          const double mean_corr = cn ? std::max(0.0, csum / cn) : 0.0;

          // Shape AND co-elution, weighted by evidence (n-1)/3 capped at 1: the cosine of a short vector is
          // nearly always high, so an unweighted 2-isotope hypothesis would beat 4 real isotopes.
          const double evidence = std::min(1.0, (found - 1) / 3.0);
          hyps.push_back({z, mono_mz, cos_sim * mean_corr * evidence, members});
        }
      }

      Precursor_ pc;
      pc.trace_idx = seed;
      pc.rt = s_rt; pc.im = s.im;

      if (hyps.empty())
      {
        pc.mono_mz = s.mz; pc.charge = 0;   // no envelope evidence at all -> unknown charge
        R.emit.push_back(pc);
        continue;
      }

      sort(hyps.begin(), hyps.end(), [](const Hyp& a, const Hyp& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.z != b.z) return a.z > b.z;   // [route-1] ties now favour HIGHER z (old code favoured low z)
        return a.mono < b.mono;
      });
      const Hyp& best = hyps[0];
      pc.mono_mz = best.mono; pc.charge = best.z;
      // partners + the seed = envelope size. `members` excludes the seed, so this is the whole
      // envelope only because the seed is always one of its peaks.
      pc.n_isotopes = 1 + (int)best.partners.size();   // [Q4] envelope confidence for the mono call
      R.emit.push_back(pc);
      for (size_t j : best.partners) R.claims.push_back(j);
      } while (0);
    };

    // [outer] 2048 (was 8192): the commit is serial and mostly fixed cost per seed, but each re-evaluation in it costs a
    // whole seed (~10 us); a quarter of the batch re-evaluates a third of the seeds (a 60-min run: 17,853 vs 51,891) for 4x the
    // fork/joins. Its greedy walk 1.45 -> 1.19 s at 64 threads, 1.55 -> 1.48 s at 192 (timed 2026-10-01). Committing
    // beside the next batch's speculation was exact but slower (more re-evaluations, a contended commit): not kept.
    size_t spec_batch = 2048;
    if (const char* e = std::getenv("DIASPEXTRACT_SPEC_BATCH")) { long v = atol(e); if (v >= 1) spec_batch = (size_t)v; }
    // [outer] about four speculation chunks per thread (at most 16 seeds, as before), so a small batch feeds every thread
    const int spec_chunk = (int)std::max<size_t>(1, std::min<size_t>(16, spec_batch / (4 * (size_t)std::max(1, omp_get_max_threads()))));
    vector<SeedRes> res(spec_batch);
    size_t n_committed = 0, n_conflict = 0, n_skipped = 0;
    std::array<size_t, 8> mc_tot{};   // [monocheck] SeedRes::mc over the committed seeds
    vector<MonoV2In_> v2in;           // [monocheck-v2] parallel to `out` when mono_check_ == 2
    double t_spec = 0.0, t_commit = 0.0, t_reeval = 0.0;   // [outer] the walk's sub-steps, log only
    auto now_s = []() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
    for (size_t b0 = 0; b0 < order.size(); b0 += spec_batch)
    {
      const size_t b1 = std::min(b0 + spec_batch, order.size());
      const double ts0 = now_s();
      // Phase 1: speculate. Nothing here writes `used`, so the concurrent reads are safe.
      #pragma omp parallel for schedule(dynamic, spec_chunk)
      for (long ii = (long)b0; ii < (long)b1; ++ii)
      {
        SeedRes& R = res[(size_t)ii - b0];
        const size_t seed = order[(size_t)ii];
        if (used[seed]) { R.emit.clear(); R.claims.clear(); R.deps.clear(); R.mc.fill(0); R.v2 = MonoV2In_(); continue; }
        evalSeed(seed, R);
      }
      const double tc0 = now_s();
      t_spec += tc0 - ts0;
      // Phase 2: commit in the original seed order.
      for (size_t ii = b0; ii < b1; ++ii)
      {
        const size_t seed = order[ii];
        if (used[seed]) { ++n_skipped; continue; }
        SeedRes& R = res[ii - b0];
        bool stale = false;
        for (size_t d : R.deps) if (used[d]) { stale = true; break; }
        if (stale) { const double tr0 = now_s(); evalSeed(seed, R); ++n_conflict; t_reeval += now_s() - tr0; }
        for (const auto& pc : R.emit) { out.push_back(pc); if (mono_check_ == 2) v2in.push_back(R.v2); }
        if (mono_check_) for (size_t q = 0; q < mc_tot.size(); ++q) mc_tot[q] += R.mc[q];
        used[seed] = true;
        for (size_t j : R.claims) used[j] = true;
        ++n_committed;
      }
      t_commit += now_s() - tc0;
    }
    // [z1-diag] DIASPEXTRACT_Z1_DIAG[=file]: for every z=1 call, is there a co-eluting MS1 trace at +ISO/2
    // (a halved 2+ ion's M+1) or +ISO/3, and was it already claimed.
    if (const char* dz = std::getenv("DIASPEXTRACT_Z1_DIAG"))
    {
      size_t n1 = 0, has2 = 0, has3 = 0, has2_free = 0;
      FILE* fo = (*dz && std::strcmp(dz, "1") != 0) ? std::fopen(dz, "w") : nullptr;
      if (fo) std::fprintf(fo, "mono_mz\trt\tim\tn_iso\tp2\tp2_free\tp3\n");
      for (const auto& pc : out)
      {
        if (pc.charge != 1) continue;
        ++n1;
        const Trace& sd = ms1[pc.trace_idx];
        const double sd_rt = rtOf(ms1st, sd);   // [soa-rt]
        const double tol = sd.mz * mass_ppm * 1e-6;
        auto has = [&](double target, bool& free_) -> bool {
          const size_t klo = lbmz(target - tol), khi = ubmz(target + tol);
          for (size_t k = klo; k < khi; ++k)
          {
            if (fabs(sorted_rt[k] - sd_rt) > delta_rt || fabs(sorted_im[k] - sd.im) > iso_im_tol) continue;
            const size_t j = by_mz[k]; if (j == pc.trace_idx) continue;
            free_ = !used[j]; return true;
          }
          return false;
        };
        bool f2 = false, f3 = false;
        const bool p2 = has(pc.mono_mz + ISO / 2.0, f2), p3 = has(pc.mono_mz + ISO / 3.0, f3);
        has2 += p2; has3 += p3; has2_free += (p2 && f2);
        if (fo) std::fprintf(fo, "%.5f\t%.3f\t%.5f\t%d\t%d\t%d\t%d\n", pc.mono_mz, pc.rt, pc.im, pc.n_isotopes, (int)p2, (int)f2, (int)p3);
      }
      if (fo) std::fclose(fo);
      writeLogInfo_("[z1-diag] " + String(n1) + " z=1 precursors: " + String(has2) + " (" + String((int)(1000.0 * has2 / std::max<size_t>(n1, 1)) / 10.0)
                    + "%) have a co-eluting MS1 trace at +ISO/2, i.e. where a 2+ ion's M+1 would be (" + String(has2_free)
                    + " of those unclaimed); " + String(has3) + " at +ISO/3");
    }
    {
      auto _pi_t3 = std::chrono::steady_clock::now();
      auto secs_ = [](auto d) { return std::chrono::duration<double>(d).count(); };
      writeLogInfo_("[infer] " + String(N) + " MS1 traces -> " + String(out.size()) + " precursors: mz-index "
                    + String(secs_(_pi_t1 - _pi_t0)) + " s | intensity-sort " + String(secs_(_pi_t2 - _pi_t1))
                    + " s | greedy walk " + String(secs_(_pi_t3 - _pi_t2)) + " s (speculative, batch "
                    + String(spec_batch) + ": " + String(n_committed) + " committed, " + String(n_conflict)
                    + " re-evaluated (" + String((int)(1000.0 * n_conflict / std::max<size_t>(n_committed, 1)) / 10.0)
                    + "%), " + String(n_skipped) + " already claimed)");
      writeLogInfo_("[outer] infer m/z order: " + (mz_nonincr ? String(mz_nonincr) + " neighbours not strictly increasing, std::sort's order (stdSortPar_)"
                                                              : String("strictly increasing, identity"))
                    + "; greedy walk: " + String((order.size() + spec_batch - 1) / spec_batch) + " batches, chunk " + String(spec_chunk)
                    + " | speculate " + String(t_spec) + " s | commit " + String(t_commit) + " s (re-evaluation " + String(t_reeval) + " s)");
    }
    if (mono_check_ == 1)
    {
      mono_check_stamp_ = String(envelope_scoring ? "n/a (charge:scoring envelope)" : "on") + ", min_corr " + String(mono_check_min_corr_)
                          + ", max_dim " + String(mono_check_max_dim_) + ": " + String(mc_tot[1]) + " of " + String(mc_tot[0]) + " monos moved right";
      writeLogInfo_("[monocheck] " + String(mc_tot[0]) + " count walks checked (winning charge with a leftward partner): "
                    + String(mc_tot[1]) + " monoisotopes moved right; steps rejected by profile (r < " + String(mono_check_min_corr_) + ") "
                    + String(mc_tot[2]) + ", by 1/K0 (> " + String(mono_check_max_dim_) + ") " + String(mc_tot[3]) + "; steps accepted "
                    + String(mc_tot[7]) + " by profile, " + String(mc_tot[4]) + " unchecked (< 3 shared frames or flat); " + String(mc_tot[5])
                    + " leftward partners released unclaimed; " + String(mc_tot[6]) + " precursors kept no claimed partner"
                    + (envelope_scoring ? " (charge:scoring envelope: not applied)" : ""));
    }
    if (v2_out) v2_out->swap(v2in);   // [monocheck-v2] parallel to `out`, for monoPositionPass_
    return out;
  }

  /// [monocheck-v2] charge:mono_profile_check v2, after the IM charge veto: each count walk is judged at its
  /// precursor's FINAL charge. A z=1 walk the veto re-calls 2+ stepped 1.003 m/z, two isotopes at 2+; judged at z=1
  /// spacing its steps and extension land two isotopes off (on a tuning run, identified 1.5.0 keys from such walks:
  /// 54.5 % not monoisotopic in 1.5.0, 55.7 % after a z=1-spaced v2), so its walk is read at the final charge:
  /// its leftmost partner is position -2 nl there, and the positions between are filled by whatever traces sit on
  /// them. PS(k) sums the apex intensity of EVERY MS1 trace (claimed or not, the seed included) at the seed's isotope
  /// position k within im_box 1/K0, in m/z order (deterministic): findPartner takes the first box candidate in m/z
  /// order and IM slicing splits one isotope over several traces, so single-trace ratios are noisy. A: an identified
  /// walk carries at most one spurious step and it is the leftmost (engines search -1..2), so step j = the leftmost is
  /// rejected when PS(-j) / PS(-j+1) x lambda < min_ratio and the mono moves one isotope right (iterate: repeat).
  /// B: when A rejected nothing and the walk ended on a miss (< 5 steps, not a charge:mono_averagine_guard stop), the
  /// true mono may be claimed by a brighter precursor, so the next lighter position becomes the mono, unclaimed, on a
  /// ratio >= extend. Only mono_mz changes: claims, charge and n_isotopes are the greedy loop's and the veto's.
  /// im_box 0 (auto): 3.5 x the 90th percentile of |d 1/K0| between the brightest decile of the walks (commit order is
  /// seed-intensity order) and their nearest trace one isotope heavier at the walk's charge (partner m/z tolerance,
  /// delta_rt, iso_im_tol), clamped to [0.01, 0.03]. `pcs` must be inferPrecursors_'s output in its order (the veto
  /// edits charges in place); `v2in` is the walk record inferPrecursors_ returned beside it.
  void monoPositionPass_(const vector<Trace>& ms1, const TraceStore& ms1st, vector<Precursor_>& pcs,
                         const vector<MonoV2In_>& v2in, double delta_rt, double iso_im_tol, double mass_ppm,
                         bool envelope_scoring)
  {
    const double ISO = 1.0033548;
    const double PROTON_ = 1.00727646;
    // walks tested, moved right, steps rejected, extensions tested, moved left, walks read at a veto-changed charge,
    // moved back onto the seed
    std::array<size_t, 8> mc_tot{};
    double v2_box = mono_pos_im_box_, v2_p90 = -1.0;
    size_t v2_nb = 0, v2_nhit = 0;
    if (!envelope_scoring)
    {
      if (v2in.size() != pcs.size())
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[monocheck-v2] the walk record does not match the precursor list");
      // The m/z-sorted index of inferPrecursors_ (same comparator on the same container, so the same order; the plain
      // std searches return the positions its bucket index returns).
      const size_t N = ms1.size();
      vector<size_t> by_mz;
      mzOrder_(ms1, by_mz);   // [outer] the same helper, so the same order
      vector<double> sorted_mz(N), sorted_rt(N), sorted_im(N);
      vector<float> inten(N);
      #pragma omp parallel for schedule(static)
      for (long ii = 0; ii < (long)N; ++ii)
      { const size_t i = (size_t)ii; const Trace& t = ms1[by_mz[i]]; sorted_mz[i] = t.mz; sorted_rt[i] = rtOf(ms1st, t); sorted_im[i] = t.im;
        inten[i] = (float)intensityOf(ms1st, ms1[i]); }
      auto lbmz = [&](double x) -> size_t { return (size_t)(std::lower_bound(sorted_mz.begin(), sorted_mz.end(), x) - sorted_mz.begin()); };
      auto ubmz = [&](double x) -> size_t { return (size_t)(std::upper_bound(sorted_mz.begin(), sorted_mz.end(), x) - sorted_mz.begin()); };
      if (v2_box <= 0.0)
      {
        size_t nz = 0;
        for (const MonoV2In_& v : v2in) if (v.z > 0) ++nz;
        v2_nb = (nz + 9) / 10;
        vector<size_t> idx;
        idx.reserve(v2_nb);
        for (size_t i = 0; i < pcs.size() && idx.size() < v2_nb; ++i) if (v2in[i].z > 0) idx.push_back(i);
        vector<double> gaps(idx.size(), -1.0);
        #pragma omp parallel for schedule(dynamic, 256)
        for (long ii = 0; ii < (long)idx.size(); ++ii)
        {
          const size_t seed = pcs[idx[(size_t)ii]].trace_idx;
          const Trace& s = ms1[seed];
          const double s_rt = rtOf(ms1st, s), tol = s.mz * mass_ppm * 1e-6;
          const double target = s.mz + ISO / v2in[idx[(size_t)ii]].z;
          const size_t klo = lbmz(target - tol), khi = ubmz(target + tol);
          double g = -1.0;
          for (size_t q = klo; q < khi; ++q)
          {
            if (fabs(sorted_rt[q] - s_rt) > delta_rt || by_mz[q] == seed) continue;
            const double d = fabs(sorted_im[q] - s.im);
            if (d <= iso_im_tol && (g < 0.0 || d < g)) g = d;
          }
          gaps[(size_t)ii] = g;
        }
        gaps.erase(std::remove_if(gaps.begin(), gaps.end(), [](double g) { return g < 0.0; }), gaps.end());
        v2_nhit = gaps.size();
        if (!gaps.empty())
        {
          const size_t k = (size_t)(0.9 * (double)(gaps.size() - 1));
          std::nth_element(gaps.begin(), gaps.begin() + (long)k, gaps.end());
          v2_p90 = gaps[k];
        }
        v2_box = v2_p90 < 0.0 ? 0.01 : std::min(0.03, std::max(0.01, 3.5 * v2_p90));
      }
      // Per precursor: bit 0 walk tested, 1 moved right, 2 moved back onto the seed, 3 extension tested, 4 moved left,
      // 5 read at a veto-changed charge; rej = steps rejected. Summed below in index order.
      vector<uint8_t> flag(pcs.size(), 0), rej(pcs.size(), 0);
      #pragma omp parallel for schedule(dynamic, 1024)
      for (long ii = 0; ii < (long)pcs.size(); ++ii)
      {
        const MonoV2In_& v = v2in[(size_t)ii];
        Precursor_& pc = pcs[(size_t)ii];
        // The final charge: the veto re-calls a z=1 walk 2+ (a multiple of the walk's charge) or unsets it (0: the
        // precursor becomes guessed and is dropped or defaulted later, untested).
        if (v.z <= 0 || pc.charge <= 0 || pc.charge % v.z != 0) continue;
        const int z = pc.charge, mul = pc.charge / v.z, nl = v.nl * mul;
        const size_t seed = pc.trace_idx;
        const Trace& s = ms1[seed];
        const double s_rt = rtOf(ms1st, s), tol = s.mz * mass_ppm * 1e-6;
        auto PS = [&](int k, long* best) -> double {
          const double target = s.mz + k * ISO / z;
          const size_t klo = lbmz(target - tol), khi = ubmz(target + tol);
          double sum = 0.0; float bi = -1.0f;
          for (size_t q = klo; q < khi; ++q)
          {
            if (fabs(sorted_rt[q] - s_rt) > delta_rt || fabs(sorted_im[q] - s.im) > v2_box) continue;
            const size_t j = by_mz[q];
            sum += inten[j];
            if (best && inten[j] > bi) { bi = inten[j]; *best = (long)j; }
          }
          return sum;
        };
        // q of the step onto position -j: PS(-j) / PS(-j+1) x lambda(neutral mass of the mono after that step). A walk
        // re-read at a veto charge never stood on its in-between positions: an empty one there is no evidence against
        // the step, so the step is compared over two positions, PS(-j) / PS(-j+2) x lambda^2 / 2 (Poisson M+2/M0).
        auto ratio = [&](int j, long* best) -> double {
          const double neutral = (s.mz - j * ISO / z) * z - z * PROTON_;
          const double lam = std::max(0.10, 0.000594 * neutral);
          const double den = PS(-j + 1, nullptr);
          if (den > 0.0) return PS(-j, best) / den * lam;
          if (mul > 1 && (j - 1) % mul != 0)
          {
            const double den2 = PS(-j + 2, nullptr);
            return den2 > 0.0 ? PS(-j, best) / den2 * lam * lam / 2.0 : 0.0;
          }
          return 0.0;
        };
        // the m/z standing at position -i: the walk's own partner where the walk stepped, else the most intense trace
        // there in the box, else the position itself
        auto leftMz = [&](int i) -> double {
          if (i % mul == 0) return v.left_mz[(size_t)(i / mul - 1)];
          long bt = -1;
          PS(-i, &bt);
          return bt >= 0 ? ms1[(size_t)bt].mz : s.mz - i * ISO / z;
        };
        double mono = v.mono;
        uint8_t f = mul > 1 ? 32 : 0, r = 0;
        int acc = nl;
        if (nl > 0)
        {
          f |= 1;
          while (acc > 0 && ratio(acc, nullptr) < mono_pos_min_ratio_)
          {
            --acc; ++r;
            if (!mono_pos_iterate_)
            {
              // Leftmost step only, at a veto charge: an EMPTY in-between position cannot be the monoisotope (leftMz
              // would report a bare position), so the mono moves on to the walk's next own position. mul == 1: no-op.
              while (acc > 0 && acc % mul != 0 && PS(-acc, nullptr) <= 0.0) --acc;
              break;
            }
          }
          if (acc < nl)
          {
            const double acc_mz = acc > 0 ? leftMz(acc) : s.mz;
            if (acc_mz > mono) { mono = acc_mz; f |= 2; if (acc == 0) f |= 4; }   // mono_averagine_select may already sit further right
          }
        }
        // B: only when A rejected nothing, the walk ended on a miss (< 5 steps, not a guard stop) and the mono is the
        // walk's leftmost peak (mono_averagine_select may have moved it right).
        const double leftmost_mz = v.nl > 0 ? v.left_mz[(size_t)(v.nl - 1)] : s.mz;
        if (mono_pos_extend_ > 0.0 && acc == nl && v.nl < 5 && !v.guard_stop && mono == leftmost_mz)
        {
          f |= 8;
          long bt = -1;
          const double q = ratio(nl + 1, &bt);
          if (q >= mono_pos_extend_ && q > 0.0 && bt >= 0) { mono = ms1[(size_t)bt].mz; f |= 16; pc.mono_ext = true; }
        }
        pc.mono_mz = mono;
        flag[(size_t)ii] = f; rej[(size_t)ii] = r;
      }
      for (size_t i = 0; i < pcs.size(); ++i)
      {
        const uint8_t f = flag[i];
        mc_tot[0] += f & 1; mc_tot[1] += (f >> 1) & 1; mc_tot[6] += (f >> 2) & 1; mc_tot[3] += (f >> 3) & 1; mc_tot[4] += (f >> 4) & 1;
        mc_tot[5] += (f >> 5) & 1; mc_tot[2] += rej[i];
      }
    }
    String box = String(mono_pos_im_box_);
    if (mono_pos_im_box_ <= 0.0)
    {
      char b[160];
      if (envelope_scoring) std::snprintf(b, sizeof b, "auto");
      else std::snprintf(b, sizeof b, "auto %.4f (3.5 x p90 %.4f of the +1-isotope 1/K0 gap, %zu of %zu brightest seeds)", v2_box,
                         v2_p90, v2_nhit, v2_nb);
      box = String(b);
    }
    const String cfg = "v2, min_ratio " + String(mono_pos_min_ratio_) + ", im_box " + box + ", extend "
                       + String(mono_pos_extend_) + (mono_pos_iterate_ ? ", iterate" : ", leftmost step only") + ", at the final charge";
    mono_check_stamp_ = (envelope_scoring ? String("n/a (charge:scoring envelope), ") : String("")) + cfg + ": " + String(mc_tot[1]) + " of "
                        + String(mc_tot[0]) + " monos moved right, " + String(mc_tot[4]) + " of " + String(mc_tot[3]) + " moved left";
    writeLogInfo_("[monocheck] " + cfg + ": " + String(mc_tot[0]) + " count walks with a leftward partner tested, " + String(mc_tot[2])
                  + " steps rejected, " + String(mc_tot[1]) + " monoisotopes moved right (" + String(mc_tot[6]) + " back onto the seed); "
                  + String(mc_tot[3]) + " walks tested for extension, " + String(mc_tot[4]) + " monoisotopes moved left (unclaimed); "
                  + String(mc_tot[5]) + " walks read at the charge the IM veto gave them; claims unchanged"
                  + (envelope_scoring ? " (charge:scoring envelope: not applied)" : ""));
  }

  /// [soa2] The fragment RT gate field. 99.4% of scoring visits fail the RT gate, so `q` quantises RT into
  /// buckets never narrower than delta_rt and a 2-byte test precedes the unchanged exact one. The reject
  /// admits +/-2 buckets where +/-1 suffices in exact arithmetic: insurance against a floor() rounding
  /// across a boundary, which would silently DROP a candidate.
  struct FragRt
  {
    vector<double> rt;
    vector<uint16_t> q;
    double rt0 = 0.0, inv = 0.0;
    double operator[](size_t i) const { return rt[i]; }
    /// Signed at the low end and clamped at both: a precursor may sit outside the window's
    /// fragment RT range. The test is written `!(x > 0.0)` rather than `x <= 0.0` so that a NaN
    /// coordinate lands in bucket 0 instead of reaching (int)NaN, which is undefined behaviour.
    int bucketOf(double v) const
    { const double x = (v - rt0) * inv; return !(x > 0.0) ? 0 : (x >= 65535.0 ? 65535 : (int)x); }
    /// inv == 0 is the DISABLED state, not an error: every bucket is then 0, every dq is 0, and
    /// the exact test alone decides. That is the fallback for every degenerate parameter below,
    /// and it costs no branch in the scan.
    void build(const vector<Trace>& fr, const TraceStore& st, double delta_rt, long nt = 1)
    {
      const size_t n = fr.size();
      rt.resize(n); q.assign(n, 0);
      rt0 = 0.0; inv = 0.0;
      // [serial] As nt tasks. The serial fold takes the first value, then any value strictly below (above) it, so
      // chunk 0 folds serially and every later chunk contributes its smallest (largest) non-NaN value, earliest on
      // ties, folded in chunk order: the same lo and hi for any values (a NaN first value sticks, a later NaN never
      // replaces).
      vector<double> clo(nt, 0.0), chi(nt, 0.0);
      vector<char> has(nt, 0);
      forTasks_(nt, [&](long c) {
        double l = 0.0, h = 0.0; bool any = false;
        for (size_t i = (size_t)c * n / (size_t)nt; i < (size_t)(c + 1) * n / (size_t)nt; ++i)
        {
          rt[i] = rtOf(st, fr[i]);
          if (c == 0) { if (!i || rt[i] < l) l = rt[i]; if (!i || rt[i] > h) h = rt[i]; }
          else if (rt[i] == rt[i]) { if (!any || rt[i] < l) l = rt[i]; if (!any || rt[i] > h) h = rt[i]; any = true; }
        }
        clo[c] = l; chi[c] = h; has[c] = any;
      });
      double lo = clo[0], hi = chi[0];
      for (long c = 1; c < nt; ++c) if (has[c]) { if (clo[c] < lo) lo = clo[c]; if (chi[c] > hi) hi = chi[c]; }
      // A denormal delta_rt makes 1/delta_rt infinite, and an RT range spanning +/-DBL_MAX makes
      // the span infinite; both then produce inf*0 = NaN inside bucketOf. Neither survives real
      // acquisition, but both are reachable and both were UB, so every one degrades to disabled.
      if (!(delta_rt > 0.0) || !std::isfinite(delta_rt)) return;
      if (!std::isfinite(lo) || !std::isfinite(hi)) return;
      const double span = hi - lo;
      if (!std::isfinite(span)) return;
      double iv = 1.0 / delta_rt;
      // Buckets are delta_rt wide, but widened if the window's RT span would overflow uint16.
      // Wider is always safe: it can only put more fragments in reach of the exact test.
      if (span > 0.0 && span * iv > 65000.0) iv = 65000.0 / span;
      if (!std::isfinite(iv) || !(iv > 0.0)) return;
      rt0 = lo; inv = iv;
      forTasks_(nt, [&](long c) {
        for (size_t i = (size_t)c * n / (size_t)nt; i < (size_t)(c + 1) * n / (size_t)nt; ++i) q[i] = (uint16_t)bucketOf(rt[i]); });
    }
  };

  /// [foreign] One window's precursors, the candidate foreign envelopes of its spectra: every precursor whose mono
  /// m/z lies in the window, with a finite apex RT and 1/K0 and an isotope-supported charge (a guessed precursor's
  /// charge is a default, so its envelope would be fictional; assembly:require_isotope_support drops them anyway),
  /// sorted by apex RT (ties by index). Built from the run's precursor list, not from the tile's owned subset, so it
  /// is the same in any tiling.
  struct ForeignIdx_
  {
    const vector<Precursor_>* p = nullptr;
    vector<double> rt;       ///< apex RT, ascending
    vector<uint32_t> idx;    ///< into *p, in rt order
    void build(const vector<Precursor_>& v, size_t lo, size_t hi)
    {
      p = &v;
      vector<pair<double, uint32_t>> o;
      o.reserve(hi > lo ? hi - lo : 0);
      for (size_t i = lo; i < hi; ++i)
        if (std::isfinite(v[i].rt) && std::isfinite(v[i].im) && v[i].charge > 0 && !v[i].guessed) o.emplace_back(v[i].rt, (uint32_t)i);
      std::sort(o.begin(), o.end());
      rt.resize(o.size()); idx.resize(o.size());
      for (size_t i = 0; i < o.size(); ++i) { rt[i] = o[i].first; idx[i] = o[i].second; }
    }
  };

  /// [foreign] assembly:foreign_precursor_weight: attenuate, among one precursor's admitted fragments, the residual
  /// (unfragmented) envelopes of the OTHER precursors isolated with it -- apex within drt of this one's, 1/K0 within
  /// foreign_dim_, neither a twin (mono within 0.01) nor this ion called again (either mono on the other's isotope
  /// grid, up to 3 steps, within max(foreign_ppm_, 30 ppm)). A candidate's positions are mono + k ISO/z, k = 0..3; a
  /// position is present when an admitted fragment lies within foreign_ppm_ of it. With at least foreign_min_iso_
  /// positions present, those fragments are multiplied by foreign_w_ (0 removes them) and their cap key is recomputed
  /// from the attenuated IM-weighted apex (@p frag_ai, attenuated with it). A fragment within max(foreign_ppm_, 30 ppm) of this
  /// precursor's own envelope (mono to M+3) neither counts nor is touched: stripping the own isotopes lost 2.3 % of
  /// peptides in tuning. Two matched positions are weak evidence for a 1+ (or same-charge) fragment and its own
  /// isotope, which line up with a co-isolated envelope of that charge as soon as one does (open). Returns the fragments
  /// attenuated. @p resc: gate:anchor_rescue's rescued flags, @p frag_kq: the [fuse] key factors; each parallel to @p frags
  /// when not empty, compacted with it.
  size_t attenuateForeign_(const Precursor_& pc, double drt, const ForeignIdx_& fx, vector<pair<double, double>>& frags,
                           vector<double>& frag_scores, vector<uint32_t>& frag_bin, vector<pair<int, int>>& frag_span,
                           vector<double>& frag_ai, vector<double>& frag_kq, vector<char>& resc) const
  {
    // The two own-ion tests take a wider tolerance than the envelope match: an ion's residual peaks scatter tens of
    // ppm around the MS1 mono grid (a 10 ng HeLa run: 31 % of the own-grid peaks within 60 ppm lie 20-60 ppm off), and with
    // both tests at foreign_ppm_ (10) 15 % of the removed fragments sat 10-60 ppm off the precursor's own M..M+3
    // (chance 0.06 %): the same ion called again with a mono just outside the sibling test took them.
    const double ISO = 1.0033548, ppm = foreign_ppm_ * 1e-6, same = std::max(ppm, 30e-6);
    auto& cnt = foreignCount_();
    ++cnt[8];
    auto on_grid = [&](double a, double b, int z) {                 // a within `same` of b + j ISO/z, |j| <= 3
      for (int j = -3; j <= 3; ++j) if (std::fabs(a - (b + j * ISO / z)) <= a * same) return true;
      return false; };
    vector<const Precursor_*> cand;
    long long n_skip = 0;                                            // added to the shared counter once per call
    const size_t lo = std::lower_bound(fx.rt.begin(), fx.rt.end(), pc.rt - drt) - fx.rt.begin();
    for (size_t j = lo; j < fx.rt.size() && fx.rt[j] <= pc.rt + drt; ++j)
    {
      const Precursor_& f = (*fx.p)[fx.idx[j]];
      if (&f == &pc || std::fabs(f.rt - pc.rt) > drt || std::fabs(f.im - pc.im) > foreign_dim_) continue;
      if (std::fabs(f.mono_mz - pc.mono_mz) < 0.01 || (pc.charge > 0 && on_grid(f.mono_mz, pc.mono_mz, pc.charge))
          || on_grid(pc.mono_mz, f.mono_mz, f.charge)) { ++n_skip; continue; }
      cand.push_back(&f);
    }
    if (n_skip) cnt[4] += n_skip;
    if (cand.empty()) return 0;
    ++cnt[0];
    vector<uint32_t> om(frags.size());                               // the fragments by m/z, ties by index
    for (size_t i = 0; i < om.size(); ++i) om[i] = (uint32_t)i;
    std::sort(om.begin(), om.end(), [&](uint32_t a, uint32_t b) { return frags[a].first != frags[b].first ? frags[a].first < frags[b].first : a < b; });
    auto own = [&](double mz) {
      if (pc.charge <= 0) return false;
      for (int k = 0; k <= 3; ++k) { const double m = pc.mono_mz + k * ISO / pc.charge; if (std::fabs(mz - m) <= m * same) return true; }
      return false; };
    vector<char> hit(frags.size(), 0);
    vector<uint32_t> at;
    long long n_env = 0, n_match = 0, err = 0;
    for (const Precursor_* f : cand)
    {
      at.clear();
      int present = 0; long long e = 0;
      for (int k = 0; k <= 3; ++k)
      {
        const double m = f->mono_mz + k * ISO / f->charge, tol = m * ppm;
        bool any = false;
        auto it = std::lower_bound(om.begin(), om.end(), m - tol, [&](uint32_t i, double v) { return frags[i].first < v; });
        for (; it != om.end() && frags[*it].first <= m + tol; ++it)
          if (!own(frags[*it].first)) { at.push_back(*it); e += std::llround((frags[*it].first - m) / m * 1e9); any = true; }
        present += any;
      }
      if (present < foreign_min_iso_) continue;
      ++n_env; n_match += (long long)at.size(); err += e;
      for (uint32_t i : at) hit[i] = 1;
    }
    size_t n = 0;
    for (size_t i = 0; i < hit.size(); ++i)
    {
      if (!hit[i]) continue;
      ++n;
      if (foreign_w_ <= 0.0) continue;
      frags[i].second *= foreign_w_;
      const double L = std::log1p(frag_ai[i]);                       // corr_inten keys are linear in it, the others in its log
      frag_scores[i] = cap_rank_ == 0 ? frag_scores[i] * foreign_w_ : L > 0.0 ? frag_scores[i] * std::log1p(frag_ai[i] * foreign_w_) / L : frag_scores[i];
      frag_ai[i] *= foreign_w_;                                      // [fuse] a fused key sums the attenuated apex
    }
    if (n && foreign_w_ <= 0.0)                                      // remove: compact every parallel array
    {
      size_t w = 0;
      for (size_t i = 0; i < hit.size(); ++i)
        if (!hit[i])
        {
          frags[w] = frags[i]; frag_scores[w] = frag_scores[i]; frag_bin[w] = frag_bin[i]; frag_span[w] = frag_span[i]; frag_ai[w] = frag_ai[i];
          if (!frag_kq.empty()) frag_kq[w] = frag_kq[i];
          if (!resc.empty()) resc[w] = resc[i];
          ++w;
        }
      frags.resize(w); frag_scores.resize(w); frag_bin.resize(w); frag_span.resize(w); frag_ai.resize(w);
      if (!frag_kq.empty()) frag_kq.resize(w);
      if (!resc.empty()) resc.resize(w);
    }
    cnt[1] += n_env; cnt[3] += (long long)n; cnt[6] += err; cnt[7] += n_match;
    if (n) ++cnt[2];
    return n;
  }

  /// [fuse] assembly:fragment_fuse_ppm: one peak per fragment ion. The IM picker keeps 1-2 flight-time bins of an ion
  /// that spans 4-7, so one ion becomes several traces and several peaks 5-20 ppm apart; the same-bin merge reaches only
  /// exact bins. Groups the fragments within @p ppm (anchored: strongest first, each unfused fragment an anchor taking
  /// every unfused one within ppm of its own m/z; @p linkage: m/z order, joined to the previous member within ppm and to
  /// the group's first within 2 x ppm). A fused peak: m/z = intensity-weighted mean of the members' (emitted) intensities,
  /// intensity = their sum, rescued flag = the strongest member's, key = the best-keyed member's factor @p kq (its score,
  /// and shared frames under evidence) times the summed IM-weighted apex @p ai (ln(1 + sum) when @p log_key), the formula
  /// of an unfused fragment's key applied to the whole ion; it is never below the best member's key (non-negative keys, up to rounding). It sits at its
  /// members' lowest index, so the relative order of the peaks is kept; a lone fragment is kept bit for bit. @p keys and
  /// @p resc (when not empty) are parallel to @p frags and compacted with it; @p ai and @p kq are parallel to it on entry
  /// and stale on return. Deterministic: every order has a total tie-break and the sums run in m/z order. Returns the
  /// groups of >= 2 fused and adds their members to @p n_members.
  static long long fuseFragments_(vector<pair<double, double>>& frags, vector<double>& keys, vector<char>& resc,
                                  const vector<double>& ai, const vector<double>& kq, bool log_key,
                                  double ppm, bool linkage, long long& n_members)
  {
    const size_t n = frags.size();
    if (n < 2 || !(ppm > 0.0)) return 0;
    const double tol = ppm * 1e-6;
    constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> om(n);                                          // by m/z, ties by index
    for (size_t i = 0; i < n; ++i) om[i] = (uint32_t)i;
    std::sort(om.begin(), om.end(), [&](uint32_t a, uint32_t b) { return frags[a].first != frags[b].first ? frags[a].first < frags[b].first : a < b; });
    vector<uint32_t> grp(n, NONE);
    uint32_t ng = 0;
    if (linkage)
    {
      size_t s = 0;                                                  // the current group's first member in om
      for (size_t k = 0; k < n; ++k)
      {
        if (k > 0)
        {
          const double m = frags[om[k]].first, mp = frags[om[k - 1]].first, ms = frags[om[s]].first;
          if (!(m - mp <= mp * tol && m - ms <= 2.0 * ms * tol)) { ++ng; s = k; }
        }
        grp[om[k]] = ng;
      }
      ++ng;
    }
    else
    {
      vector<uint32_t> oi(om);                                       // by intensity desc, then m/z asc, then index
      std::sort(oi.begin(), oi.end(), [&](uint32_t a, uint32_t b) {
        if (frags[a].second != frags[b].second) return frags[a].second > frags[b].second;
        if (frags[a].first != frags[b].first) return frags[a].first < frags[b].first;
        return a < b;
      });
      for (uint32_t a : oi)
      {
        if (grp[a] != NONE) continue;
        const double m = frags[a].first, d = m * tol;
        auto it = std::lower_bound(om.begin(), om.end(), m - d, [&](uint32_t i, double v) { return frags[i].first < v; });
        for (; it != om.end() && frags[*it].first <= m + d; ++it) if (grp[*it] == NONE) grp[*it] = ng;
        ++ng;
      }
    }
    if (ng == n) return 0;                                           // nothing within ppm: unchanged
    struct Acc { double w = 0.0, wmz = 0.0, ai = 0.0; uint32_t head = NONE, top = NONE, best = NONE, n = 0; };
    vector<Acc> acc(ng);
    for (uint32_t i : om)                                            // m/z order: a fixed summation order
    {
      Acc& g = acc[grp[i]];
      g.w += frags[i].second; g.wmz += frags[i].second * frags[i].first;
      g.ai += ai[i];
      if (g.best == NONE || keys[i] > keys[g.best]) g.best = i;      // the best key; ties: lower m/z
      g.head = std::min(g.head, i);
      if (g.top == NONE || frags[i].second > frags[g.top].second) g.top = i;   // the strongest; ties: lower m/z
      ++g.n;
    }
    long long groups = 0;
    size_t w = 0;
    for (size_t i = 0; i < n; ++i)                                   // write position w <= i <= every unread member
    {
      const Acc& g = acc[grp[i]];
      if (g.head != i) continue;
      if (g.n == 1) { frags[w] = frags[i]; keys[w] = keys[i]; if (!resc.empty()) resc[w] = resc[i]; }
      else
      {
        const double mz = g.w > 0.0 ? g.wmz / g.w : frags[g.top].first;
        const char r = resc.empty() ? 0 : resc[g.top];
        frags[w] = {mz, g.w}; keys[w] = kq[g.best] * (log_key ? std::log1p(g.ai) : g.ai);
        if (!resc.empty()) resc[w] = r;
        ++groups; n_members += g.n;
      }
      ++w;
    }
    frags.resize(w); keys.resize(w);
    if (!resc.empty()) resc.resize(w);
    return groups;
  }

  /// Per-precursor assembly: this precursor claims every fragment passing its gate (a fragment may
  /// be shared across precursors -> chimeric). Leaves @p out empty on failure.
  /// pdense is caller-owned scratch, so the caller decides the per-thread copy. Keeps the top max_frags
  /// by score, drops if < min_frags, annotates the synthetic precursor. Deterministic.
  void assembleOne_(const Precursor_& pc, double win_lo, double win_hi,
                    const vector<Trace>& frag_traces, const TraceStore& wst, const vector<double>& frag_im,
                    const FragRt& frag_rt,
                    const vector<Trace>& ms1_traces, const TraceStore& ms1st, const FragStats& fg,
                    double delta_im, double delta_rt, double min_corr, int min_corr_pts,
                    Size min_frags, Size max_frags, vector<float>& pdense, MSSpectrum& out,
                    const ForeignIdx_* fx = nullptr) const
  {
    vector<pair<double, double>> frags;
    vector<double> frag_scores;
    vector<double> frag_ai;      // [foreign]/[fuse] each fragment's IM-weighted apex, its cap key's intensity (only with fx or fusion)
    vector<double> frag_kq;      // [fuse] each fragment's cap key factor besides that intensity term (only with fusion)
    const bool fuse_on = fuse_ppm_ > 0.0;
    vector<uint32_t> frag_bin;   // [merge] the TOF bin of each fragment (0: OpenMS detector, no bin; a real bin 0, m/z just
                                 // below the usual 100 m/z floor, is compared by exact m/z like the OpenMS path)
    vector<pair<int, int>> frag_span;   // [merge] its detected window-frame span
    const double G = (double)fg.G;
    const Trace& p_tr = ms1_traces[pc.trace_idx];
    vector<int> touched;
    touched.reserve(p_tr.np());
    for (size_t k = 0; k < p_tr.span(); ++k)
    {
      if (!real(ms1st, p_tr, k)) continue;
      const double pval = (double)xv(ms1st, p_tr, k);
      const int gi = (*fg.nearest_local)[p_tr.frame0 + k];   // MS1 frame -> window frame, precomputed
      if (gi >= 0) { if (pdense[gi] == 0.0f) touched.push_back(gi); pdense[gi] += (float)pval; }
    }
    double psum = 0, psumsq = 0, plogsum = 0;
    // plogsum normalises the precursor profile for gate:coelution=logoverlap. Constant across
    // every fragment of this precursor, so hoist it out of the inner loop.
    for (int gi : touched)
    {
      double v = pdense[gi];
      psum += v; psumsq += v * v;
      if (log_overlap_) plogsum += std::log1p(v);
    }
    int p_lo = std::numeric_limits<int>::max(), p_hi = -1;              // [short] the precursor's span in window frames
    int p_apex = -1;                                                    // [coloc] the precursor apex frame
    for (int gi : touched)
    {
      p_lo = std::min(p_lo, gi); p_hi = std::max(p_hi, gi);
      if (p_apex < 0 || pdense[gi] > pdense[p_apex]) p_apex = gi;
    }
    // [fragint] the frames a fragment's amplitude is summed over: core, the K frames centred on the precursor apex
    // (ties to the earlier frame); sum, the precursor's span with its gaps: a faint MS1 trace has missing points and the
    // concentric gate keeps fragments with no point where pdense != 0 (summed only there, 29% of HeLa 10 ng's had none)
    int fint_apex = -1, fint_lo = 0, fint_hi = -1;
    long long n_fint = 0, n_fint_fb = 0;
    if (frag_int_ != 0)
    {
      for (int gi : touched) if (fint_apex < 0 || pdense[gi] > pdense[fint_apex] || (pdense[gi] == pdense[fint_apex] && gi < fint_apex)) fint_apex = gi;
      fint_lo = frag_int_ == 1 ? fint_apex - frag_core_half_ : p_lo;
      fint_hi = frag_int_ == 1 ? fint_apex + frag_core_half_ : p_hi;
    }
    auto fintAmp = [&](const Trace& f, double dim)                     // [fragint] a gate-passing or rescued fragment's amplitude
    {
      const int f0 = (int)f.frame0, f1 = f0 + (int)f.span() - 1;
      double b = 0.0;                                                  // the same frames for every fragment of this precursor
      for (int g = std::max(fint_lo, f0); g <= std::min(fint_hi, f1); ++g) b += xv(wst, f, (size_t)(g - f0));
      if (!(b > 0.0))                                                  // no point there: the point nearest the precursor apex (ties earlier),
      {                                                                // so the peak set stays the apex mode's and no peak is written at 0
        int kb = -1;
        for (int k = 0; k < (int)f.span(); ++k)
          if (real(wst, f, (size_t)k) && (kb < 0 || std::abs(f0 + k - fint_apex) < std::abs(f0 + kb - fint_apex))) kb = k;
        if (kb >= 0) b = xv(wst, f, (size_t)kb);
        ++n_fint_fb;
      }
      if (im_weight_sigma_ > 0.0) b *= std::exp(-(dim * dim) / (2.0 * im_weight_sigma_ * im_weight_sigma_));
      ++n_fint;
      return b;
    };
    // [adaptive] with gate:adaptive faint, the relaxed modes apply only to faint precursors (short precursor traces)
    const bool faint = adaptive_ && (int)touched.size() <= faint_frames_, all_p = !adaptive_ || faint;
    { const size_t n = touched.size(); pspanHist_()[n <= 4 ? (n ? n - 1 : 0) : n <= 6 ? 4 : n <= 10 ? 5 : n <= 20 ? 6 : 7]++; if (faint) pspanHist_()[8]++; }
    const bool g_shrink = shrink_ && all_p, g_coloc = coloc_ && all_p, g_conc = short_conc_ && all_p;
    const double g_drt = faint ? faint_drt_ : delta_rt, g_cp = faint ? faint_cp_ : corr_power_;
    const double p_ctr = 0.5 * (p_lo + p_hi), p_half = 0.5 * (p_hi - p_lo);
    double pmean = psum / G, pvar = psumsq - G * pmean * pmean;
    double pinv = pvar > 0 ? 1.0 / sqrt(pvar) : 0.0;
    // [anchor] gate:anchor_rescue, well-sampled precursors (signal in more than gate:faint_max_frames window frames, the
    // faint test of gate:adaptive, whatever that mode is set to): the accepted fragments that can anchor a consensus, and
    // the ones the gate rejected only for support (point or overlap floor) or by the correlation cut, with at least 2 points
    const bool g_anchor = anchor_rescue_ && (int)touched.size() > faint_frames_;
    vector<uint32_t> anchors, resc_fi;                                  // resc_fi: the rescued indices, sorted, for the gate dump
    vector<size_t> cands;                                               // cands: index << 2 | rejection (0 points, 1 overlap, 2 corr)
    vector<char> resc;                                                  // [anchor] 1 = rescued, parallel to frags (empty: none)
    if (pinv > 0.0)
    {
      size_t lo = lower_bound(frag_im.begin(), frag_im.end(), pc.im - delta_im) - frag_im.begin();
      size_t hi = upper_bound(frag_im.begin(), frag_im.end(), pc.im + delta_im) - frag_im.begin();
      const int pq = frag_rt.bucketOf(pc.rt);
      for (size_t fi = lo; fi < hi; ++fi)
      {
        const int dq = (int)frag_rt.q[fi] - pq;
        if (dq < -2 || dq > 2) continue;                                // 2 B reject, see FragRt
        // [soa] the gate reads a parallel array; the trace record is touched only after it passes.
        if (fabs(frag_rt[fi] - pc.rt) > g_drt) continue;               // RT gate ([adaptive] faint window)
        if (fg.invnorm[fi] == 0.0) continue;                           // degenerate fragment
        const Trace& f = frag_traces[fi];
        if (fabs(f.mz - pc.mono_mz) < 0.01) continue;                  // exclude precursor peak [M-7]
        double dot = 0; int overlap = 0;
        // [R2-late, 2026-09-02] overlap <= support, so a fragment with fewer support points than
        // min_corr_pts can never pass the overlap guard below: skip it before the dot product.
        // Output-identical by construction (the guard's continue is the only thing it would reach).
        const bool short_np = f.np() < (size_t)std::max(min_corr_pts, 0);
        if (short_np && !g_conc && !(g_shrink && (int)f.np() >= shrink_min_pts_))
        { if (g_anchor && f.np() >= 2) cands.push_back(fi << 2); continue; }
        int f_lo = -1, f_hi = -1, f_apex = -1;                         // [short]/[coloc] the fragment's detected span and apex
        float f_max = 0.0f;
        // two sequential streams: the fragment's span in the arena and pdense over window frames
        for (size_t k = 0; k < f.span(); ++k)
        { const float fv = xv(wst, f, k); if (fv == 0.0f) continue;
          if (f_lo < 0) f_lo = (int)(f.frame0 + k);
          f_hi = (int)(f.frame0 + k);
          if (fv > f_max) { f_max = fv; f_apex = (int)(f.frame0 + k); }
          const float pv = pdense[f.frame0 + k]; if (pv != 0.0f) { dot += (double)pv * fv; ++overlap; } }
        if (g_coloc)                                                    // [coloc] DIA-Umpire-style co-location
        {
          if (f_lo < 0 || p_hi < 0) continue;
          if (f_apex < p_lo - 1 || f_apex > p_hi + 1 || p_apex < f_lo - 1 || p_apex > f_hi + 1) continue;
          const int inter = std::min(f_hi, p_hi) - std::max(f_lo, p_lo) + 1, shorter = std::min(f_hi - f_lo, p_hi - p_lo) + 1;
          if (inter < coloc_overlap_ * shorter) continue;
        }
        double c, c_raw = 0.0;                                          // c_raw: the unshrunk score, for assembly:cap_rank
        const bool no_corr = g_shrink ? overlap < shrink_min_pts_ : (short_np || overlap < min_corr_pts);   // [short]/[shrink] precedence
        if (g_conc && no_corr)
        {
          // [short] too few points for a correlation: faint ions are the tips of their elution peaks, so their
          // detected spans shrink with intensity but stay centred on the apex. Keep a fragment whose span lies inside
          // the precursor's (one cycle of slack) with its centre near the precursor's; the score ranks by centring.
          if ((int)f.np() < short_min_pts_ || f_lo < 0 || p_hi < 0) continue;
          if (f_lo < p_lo - 1 || f_hi > p_hi + 1) continue;
          const double d = std::fabs(0.5 * (f_lo + f_hi) - p_ctr), tol = std::max(1.0, short_center_tol_ * p_half);
          if (d > tol) continue;
          c = min_corr * (1.0 - 0.5 * d / (tol + 1.0));                // at the gate: ranks and weighs as a borderline fragment
          c_raw = c;
        }
        else
        {
          if (overlap < (g_shrink ? shrink_min_pts_ : min_corr_pts))    // overlap guard [H-4]; [shrink] its own floor
          { if (g_anchor && f.np() >= 2) cands.push_back(fi << 2 | 1u); continue; }
          // [B0] Plain Pearson takes variance over the full grid G, so a sparse fragment is scored against ~G
          // implicit zeros. gate:variance_support uses the union support instead; gate:coelution=logoverlap is
          // the AlphaDIA-style log-space profile overlap, summed where the precursor has signal, in [0,1].
          if (var_support_ && !log_overlap_)                             // [Q1] pearson-only, per the flag doc
          {
            // [Q1] n = |precursor support U fragment support|; fsum/fsumsq over the fragment's real points.
            double fsum = 0.0, fsumsq = 0.0;
            for (size_t k = 0; k < f.span(); ++k)
            { const double fv = xv(wst, f, k); if (fv > 0.0) { fsum += fv; fsumsq += fv * fv; } }
            double n = (double)(touched.size() + f.np() - (size_t)overlap);
            double pm = psum / n, fm = fsum / n;
            double pv2 = psumsq - n * pm * pm, fv2 = fsumsq - n * fm * fm;
            if (pv2 <= 0.0 || fv2 <= 0.0) continue;                      // degenerate on union support
            c = (dot - n * pm * fm) / std::sqrt(pv2 * fv2);
          }
          else if (log_overlap_)
          {
            // Both log-profiles normalised to unit sum, then intersected: a shape statistic. Normalising by the
            // precursor's log-sum alone bounded the score by the intensity ratio (-18.9% peptides at gate 0.5).
            double fsum = 0.0;
            for (size_t k = 0; k < f.span(); ++k)
            { const double fv = xv(wst, f, k); if (fv > 0.0) fsum += std::log1p(fv); }
            if (fsum <= 0.0 || plogsum <= 0.0) continue;
            double inter = 0.0;
            for (size_t k = 0; k < f.span(); ++k)
            {
              const double fv = xv(wst, f, k); if (fv <= 0.0) continue;
              inter += std::min(std::log1p(fv) / fsum,
                                std::log1p((double)pdense[f.frame0 + k]) / plogsum);
            }
            c = inter;
          }
          else c = (dot - G * pmean * fg.mean[fi]) * pinv * fg.invnorm[fi]; // Pearson
          c_raw = c;
          if (g_shrink)                                                  // [shrink] Fisher-z posterior mean
          {
            const double z = std::atanh(std::max(-0.999, std::min(0.999, c))), s2 = 1.0 / std::max(overlap - 3, 1);
            const double t2 = shrink_sd_ * shrink_sd_, mu = std::atanh(std::max(-0.999, std::min(0.999, shrink_r0_)));
            c = std::tanh((z / s2 + mu / t2) / (1.0 / s2 + 1.0 / t2));
          }
          if (!(c >= min_corr))                                          // correlation gate (a NaN score is rejected too)
          { if (g_anchor && f.np() >= 2) cands.push_back(fi << 2 | 2u); continue; }
        }
        // [E5] The IM weight feeds both the cap key and the emitted intensity; corr^corr_power feeds only the emitted one.
        double inten = intensityOf(wst, f);
        const double dim = frag_im[fi] - pc.im;
        if (im_weight_sigma_ > 0.0) inten *= std::exp(-(dim * dim) / (2.0 * im_weight_sigma_ * im_weight_sigma_));
        double emit_inten = inten;
        if (frag_int_ != 0) emit_inten = fintAmp(f, dim);               // [fragint] base amplitude summed over frames; the cap key keeps the apex
        if (g_cp > 0.0 && c > 0.0) emit_inten *= std::pow(c, g_cp);
        frags.emplace_back(exportMz_(f.tof, bOf(wst, f), f.mz), emit_inten);
        frag_bin.push_back(f.tof);
        frag_span.emplace_back(f_lo, f_hi);
        if (fx || fuse_on) frag_ai.push_back(inten);
        frag_scores.push_back(cap_rank_ == 1 ? std::max(c_raw, 0.0) * std::max(c_raw, 0.0) * std::log1p(inten)
                            : cap_rank_ == 2 ? std::max(c_raw, 0.0) * std::log1p(inten) * (1.0 + overlap)
                            : c * inten);                              // [cap] the max_fragments ranking key
        if (fuse_on)                                                   // [fuse] the same key without its intensity term
          frag_kq.push_back(cap_rank_ == 1 ? std::max(c_raw, 0.0) * std::max(c_raw, 0.0) : cap_rank_ == 2 ? std::max(c_raw, 0.0) * (1.0 + overlap) : c);
        if (g_anchor && c >= anchor_min_corr_) anchors.push_back((uint32_t)fi);   // [anchor] a well-correlated accepted fragment
      }
      // [anchor] gate:anchor_rescue. A weak fragment of a strong precursor is a 2-3 point tip that fails the plain gate's
      // point floor or its full-grid Pearson, while the precursor's well-correlated fragments (the anchors) outline the
      // elution peak. Their profiles, each at unit maximum, are summed into a consensus over the window frames; a candidate
      // joins when its apex is within one cycle of the consensus apex and its own points follow the consensus (local
      // Pearson; with two points: both in the core, >= half the consensus maximum, and in its direction, scored at the
      // threshold like a borderline fragment). The score, capped at anchor_min_corr, weighs it like any other fragment, and
      // its cap key counts the frames it shares with the precursor, as every other fragment's does.
      // Deterministic: anchors and candidates come in index order, the budget sorts by (score, m/z, index).
      if (g_anchor)
      {
        auto& ac = anchorCount_();
        ++ac[0];
        if ((int)anchors.size() >= anchor_min_count_)
        {
          ++ac[1]; ac[2] += (long long)cands.size();
          int c_lo = std::numeric_limits<int>::max(), c_hi = -1;
          for (uint32_t ai : anchors)
          { const Trace& a = frag_traces[ai]; c_lo = std::min(c_lo, (int)a.frame0); c_hi = std::max(c_hi, (int)(a.frame0 + a.span()) - 1); }
          vector<double> cons((size_t)(c_hi - c_lo + 1), 0.0);
          for (uint32_t ai : anchors)
          {
            const Trace& a = frag_traces[ai];
            float amax = 0.0f;
            for (size_t k = 0; k < a.span(); ++k) amax = std::max(amax, xv(wst, a, k));
            if (!(amax > 0.0f)) continue;
            for (size_t k = 0; k < a.span(); ++k) cons[a.frame0 + k - (size_t)c_lo] += (double)xv(wst, a, k) / amax;
          }
          size_t ca = 0;
          for (size_t j = 1; j < cons.size(); ++j) if (cons[j] > cons[ca]) ca = j;   // the consensus apex (first maximum)
          const int c_apex = c_lo + (int)ca;
          const double core = 0.5 * cons[ca];
          struct Resc { double r, mz; uint32_t fi; int why, ov, t_lo, t_hi; };
          vector<Resc> rs;
          for (size_t cw : cands)
          {
            const uint32_t fi = (uint32_t)(cw >> 2);
            const Trace& f = frag_traces[fi];
            if (std::abs((int)(f.frame0 + f.apex) - c_apex) > 1) continue;   // the raw apex, within one cycle
            // sums shifted by the first point, so a constant profile has exactly zero variance
            double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, x0 = 0, y0 = 0, dx1 = 0, dy1 = 0;
            int n = 0, ov = 0, t_lo = -1, t_hi = -1;
            bool in_core = true;
            for (size_t k = 0; k < f.span(); ++k)
            {
              const double x = xv(wst, f, k); if (x <= 0.0) continue;
              const int g = (int)(f.frame0 + k);
              const double y = g < c_lo || g > c_hi ? 0.0 : cons[(size_t)(g - c_lo)];
              if (t_lo < 0) { t_lo = g; x0 = x; y0 = y; }
              t_hi = g;
              if (pdense[g] != 0.0f) ++ov;                               // shared with the precursor (pdense is still live)
              if (y < core) in_core = false;
              const double dx = x - x0, dy = y - y0;
              if (n == 1) { dx1 = dx; dy1 = dy; }
              sx += dx; sy += dy; sxx += dx * dx; syy += dy * dy; sxy += dx * dy; ++n;
            }
            double r = -1.0;
            if (n == 2) { if (in_core && dx1 * dy1 > 0.0) r = anchor_rescue_corr_; }
            else if (n > 2)
            {
              const double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
              if (vx > 0.0 && vy > 0.0) r = (sxy - sx * sy / n) / std::sqrt(vx * vy);
            }
            if (!(r >= anchor_rescue_corr_)) continue;
            rs.push_back({r, exportMz_(f.tof, bOf(wst, f), f.mz), fi, (int)(cw & 3), ov, t_lo, t_hi});
          }
          std::sort(rs.begin(), rs.end(), [](const Resc& a, const Resc& b) {
            if (a.r != b.r) return a.r > b.r;
            if (a.mz != b.mz) return a.mz < b.mz;
            return a.fi < b.fi;
          });
          const size_t nk = std::min(rs.size(), (size_t)anchor_max_rescued_);
          ac[4] += (long long)(rs.size() - nk);
          if (nk) { ac[3] += (long long)nk; ++ac[5]; resc.assign(frags.size(), 0); }
          long long n_why[3] = {0, 0, 0};
          for (size_t j = 0; j < nk; ++j)
          {
            const Resc& q = rs[j];
            const Trace& f = frag_traces[q.fi];
            const double c = std::min(q.r, anchor_min_corr_);
            double inten = intensityOf(wst, f);                          // as above: IM weight, then c^corr_power
            const double dim = frag_im[q.fi] - pc.im;
            if (im_weight_sigma_ > 0.0) inten *= std::exp(-(dim * dim) / (2.0 * im_weight_sigma_ * im_weight_sigma_));
            double emit_inten = frag_int_ != 0 ? fintAmp(f, dim) : inten;  // [fragint] the same amplitude as a gate-passing fragment
            if (g_cp > 0.0 && c > 0.0) emit_inten *= std::pow(c, g_cp);
            frags.emplace_back(q.mz, emit_inten);
            frag_bin.push_back(f.tof);
            frag_span.emplace_back(q.t_lo, q.t_hi);
            if (fx || fuse_on) frag_ai.push_back(inten);                 // [foreign]/[fuse] parallel to frags, as above
            frag_scores.push_back(cap_rank_ == 1 ? c * c * std::log1p(inten) : cap_rank_ == 2 ? c * std::log1p(inten) * (1.0 + q.ov) : c * inten);
            if (fuse_on) frag_kq.push_back(cap_rank_ == 1 ? c * c : cap_rank_ == 2 ? c * (1.0 + q.ov) : c);
            resc.push_back(1);
            resc_fi.push_back(q.fi);
            ++n_why[q.why];
          }
          for (int w = 0; w < 3; ++w) if (n_why[w]) ac[8 + w] += n_why[w];
          std::sort(resc_fi.begin(), resc_fi.end());
        }
      }
    }
    if (n_fint) { fragIntCount_()[0] += n_fint; fragIntCount_()[1] += n_fint_fb; }   // [fragint] once per precursor x window
    if (GateDump_& gd = GateDump_::get(); gd.on)                       // [gate-dump] read-only; before the reset
    {
      std::ostringstream rows;
      rows.imbue(std::locale::classic());
      auto it = std::lower_bound(gd.t.begin(), gd.t.end(), pc.mono_mz - 1.2, [](const GateDump_::T& a, double v) { return a.mz < v; });
      for (; it != gd.t.end() && it->mz <= pc.mono_mz + 1.2; ++it)
      {
        if (std::fabs(it->rt - pc.rt) > 20.0 || std::fabs(it->im - pc.im) > 0.05 || it->z <= 0) continue;
        int kk = 9;
        for (int k = -1; k <= 1; ++k)
          if (std::fabs(pc.mono_mz - (it->mz + k * 1.003355 / it->z)) <= it->mz * 20e-6) { kk = k; break; }
        if (kk == 9) continue;
        auto head = [&]() { rows << it->id << '\t' << kk << '\t' << pc.mono_mz << '\t' << pc.charge << '\t' << pc.rt << '\t' << pc.im
                                 << '\t' << pc.n_isotopes << '\t' << (int)pc.guessed << '\t'; };
        head(); rows << "-\t-\t-\t-\t-\t-\t" << frags.size() << (pinv > 0.0 ? "\tPRECURSOR\n" : "\tPRECURSOR_FLAT\n");
        if (!(pinv > 0.0)) continue;
        const size_t lo = lower_bound(frag_im.begin(), frag_im.end(), pc.im - 0.05) - frag_im.begin();
        const size_t hi = upper_bound(frag_im.begin(), frag_im.end(), pc.im + 0.05) - frag_im.begin();
        for (size_t fi = lo; fi < hi; ++fi)
        {
          const double drt = frag_rt[fi] - pc.rt, dim = frag_im[fi] - pc.im;
          if (std::fabs(drt) > 15.0) continue;
          const Trace& f = frag_traces[fi];
          double dot = 0; int overlap = 0;
          for (size_t k = 0; k < f.span(); ++k)
          { const float fv = xv(wst, f, k); if (fv == 0.0f) continue;
            const float pv = pdense[f.frame0 + k]; if (pv != 0.0f) { dot += (double)pv * fv; ++overlap; } }
          const double c = fg.invnorm[fi] == 0.0 ? 0.0 : (dot - G * pmean * fg.mean[fi]) * pinv * fg.invnorm[fi];
          const char* fate = std::fabs(dim) > delta_im ? "im" : std::fabs(drt) > g_drt ? "rt" : fg.invnorm[fi] == 0.0 ? "degenerate"
                           : std::fabs(f.mz - pc.mono_mz) < 0.01 ? "precursor_peak"
                           : !resc_fi.empty() && std::binary_search(resc_fi.begin(), resc_fi.end(), (uint32_t)fi) ? "rescued"   // [anchor]
                           : f.np() < (size_t)std::max(min_corr_pts, 0) ? "np"
                           : overlap < min_corr_pts ? "overlap" : c < min_corr ? "corr" : "kept";
          head(); rows << exportMz_(f.tof, bOf(wst, f), f.mz) << '\t' << dim << '\t' << drt << '\t' << f.np() << '\t' << overlap << '\t'
                       << c << '\t' << intensityOf(wst, f) << '\t' << fate << '\n';
        }
      }
      if (rows.tellp() > 0) { std::lock_guard<std::mutex> lk(gd.mu); gd.out << rows.str(); gd.out.flush(); }
    }
    for (int gi : touched) pdense[gi] = 0.0f;                          // reset scratch
    // [foreign] other precursors' residual envelopes, before the merge and the cap: a removed fragment frees its slot
    const size_t n_before = frags.size();                              // [foreign] for the fell-below count
    const size_t n_foreign = fx && !frags.empty() ? attenuateForeign_(pc, g_drt, *fx, frags, frag_scores, frag_bin, frag_span, frag_ai, frag_kq, resc) : 0;
    // [merge] guard: one TOF bin is one m/z. What toTof's same-frame merge cannot reach (slices beyond its span, valley-
    // split children, band edges, tile carries) is merged here, before the cap, into the strongest member: members whose
    // RT span overlaps its span are summed (parallel slices of one ion), RT-disjoint pieces are dropped (the stronger
    // piece stands, as the twin fold's max); the merged peak takes the strongest member's m/z and the best cap key of
    // the summed members (not recomputed). Two peaks at one m/z stall OpenMS' MatchedIterator.
    if (merge_same_bin_ && frags.size() > 1)
    {
      vector<size_t> ord(frags.size());
      for (size_t i = 0; i < ord.size(); ++i) ord[i] = i;
      auto key = [&](size_t i) { return std::make_tuple(frag_bin[i] == 0, frag_bin[i], frag_bin[i] ? 0.0 : frags[i].first, i); };
      std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return key(a) < key(b); });
      auto same = [&](size_t a, size_t b) { return frag_bin[a] ? frag_bin[a] == frag_bin[b] : (!frag_bin[b] && frags[a].first == frags[b].first); };
      vector<char> gone(frags.size(), 0);
      long long nm = 0;
      for (size_t g = 0; g < ord.size();)
      {
        size_t e = g + 1;
        while (e < ord.size() && same(ord[g], ord[e])) ++e;
        if (e - g > 1)
        {
          const size_t head = ord[g];                                    // lowest index: the output order is unchanged
          size_t top = head;
          for (size_t k = g; k < e; ++k) if (frags[ord[k]].second > frags[top].second) top = ord[k];
          double sum = 0.0, cap = -std::numeric_limits<double>::infinity(), asum = 0.0;
          size_t kb = top;                                               // [fuse] the member whose key the peak takes
          for (size_t k = g; k < e; ++k)
          {
            const size_t i = ord[k];
            if (!(frag_span[i].second < frag_span[top].first || frag_span[top].second < frag_span[i].first))
            {
              sum += frags[i].second;
              if (frag_scores[i] > cap) kb = i;
              cap = std::max(cap, frag_scores[i]);
              if (fuse_on) asum += frag_ai[i];
            }
            if (k > g) gone[i] = 1;
          }
          frags[head] = {frags[top].first, sum}; frag_scores[head] = cap;
          if (fuse_on) { frag_ai[head] = asum; frag_kq[head] = frag_kq[kb]; }   // [fuse] the summed apex, the key's member's factor
          if (!resc.empty()) resc[head] = resc[top];                    // [anchor] rescued if its m/z is a rescued one's
          nm += (long long)(e - g - 1);
        }
        g = e;
      }
      if (nm)
      {
        size_t w = 0;
        for (size_t i = 0; i < frags.size(); ++i)
          if (!gone[i])
          {
            frags[w] = frags[i]; frag_scores[w] = frag_scores[i];
            if (fuse_on) { frag_ai[w] = frag_ai[i]; frag_kq[w] = frag_kq[i]; }
            if (!resc.empty()) resc[w] = resc[i];
            ++w;
          }
        frags.resize(w); frag_scores.resize(w);
        if (fuse_on) { frag_ai.resize(w); frag_kq.resize(w); }
        if (!resc.empty()) resc.resize(w);
        mergeCount_()[2] += nm;
      }
    }
    // [fuse] one peak per fragment ion, after the merge (frag_bin, frag_span, frag_ai and frag_kq are stale from here on)
    // and before min_fragments and the cap: the cap ranks fused peaks by a key recomputed from the summed apex, so the
    // slots of an ion's pieces refill
    const size_t n_merged = frags.size();                              // [foreign]/[fuse] the fell-below counts
    if (fuse_on && n_merged > 0)
    {
      long long n_mem = 0;
      const long long n_grp = fuseFragments_(frags, frag_scores, resc, frag_ai, frag_kq, cap_rank_ != 0, fuse_ppm_, fuse_linkage_, n_mem);
      auto& fc = fuseCount_();
      ++fc[0]; fc[1] += (long long)n_merged; fc[2] += (long long)frags.size();
      if (n_grp) { fc[3] += n_grp; fc[4] += n_mem; }
      if (n_merged >= min_frags && frags.size() < min_frags) ++fc[5];
    }
    if (frags.size() < min_frags)
    {   // [foreign] counted only when the removal (weight 0) took a spectrum that had enough fragments below the floor
      if (n_foreign && foreign_w_ <= 0.0 && n_before >= min_frags && n_merged < min_frags) ++foreignCount_()[5];
      return;
    }
    long long n_resc = 0;                                               // [anchor] rescued peaks in this spectrum
    for (char r : resc) n_resc += r;
    if (frags.size() > max_frags)
    {
      vector<size_t> idx(frags.size());
      for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
      // Rank by correlation. Ties break on m/z, then index, so the output is deterministic.
      partial_sort(idx.begin(), idx.begin() + max_frags, idx.end(), [&](size_t a, size_t b) {
        const double ka = frag_scores[a];
        const double kb = frag_scores[b];
        if (ka != kb) return ka > kb;
        if (frags[a].first != frags[b].first) return frags[a].first < frags[b].first;
        return a < b;
      });
      vector<pair<double, double>> keep; keep.reserve(max_frags);
      for (size_t i = 0; i < max_frags; ++i) keep.push_back(frags[idx[i]]);
      if (!resc.empty()) { n_resc = 0; for (size_t i = 0; i < max_frags; ++i) n_resc += resc[idx[i]]; }
      frags.swap(keep);
      if (post_fusion_cap_)                                            // [recap] the kept peaks' keys, in the same order
      {
        vector<double> keep_k(max_frags);
        for (size_t i = 0; i < max_frags; ++i) keep_k[i] = frag_scores[idx[i]];
        frag_scores.swap(keep_k);
      }
    }
    if (n_resc) { anchorCount_()[6] += n_resc; ++anchorCount_()[7]; }
    out.setMSLevel(2);
    out.setType(SpectrumSettings::SpectrumType::CENTROID);
    out.setRT(pc.rt);
    // [mem] exact reserve: no doubling slack held from emit until the write.
    out.reserve(frags.size());
    for (const auto& fr : frags) out.emplace_back(fr.first, fr.second);
    if (post_fusion_cap_)   // [recap] each peak's cap key rides along as a float data array (sorted with the peaks), dropped before the write
    {
      OpenMS::DataArrays::FloatDataArray k;
      k.setName("spx_cap_key");
      k.assign(frag_scores.begin(), frag_scores.end());
      out.getFloatDataArrays().push_back(std::move(k));
    }
    out.sortByPosition();
    OpenMS::Precursor prec;
    prec.setMZ(pc.mono_mz);   // MS1 traces come from the OpenMS detector: no bin to calibrate
    if (pc.charge > 0) prec.setCharge(pc.charge);
    prec.setIsolationWindowLowerOffset(max(0.0, pc.mono_mz - win_lo));   // [C-5]
    prec.setIsolationWindowUpperOffset(max(0.0, win_hi - pc.mono_mz));
    prec.setDriftTime(pc.im);
    prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    out.setPrecursors({prec});
    // [provenance] userParams: guessed (no isotope support) vs supported, and the isotope peaks behind the
    // mono call (0/1: the mono may be off by +-k*1.00335 Da, a hint for open-search correction).
    out.setMetaValue("spx_guessed", (int)(pc.guessed ? 1 : 0));
    out.setMetaValue("spx_n_isotopes", pc.n_isotopes);
    if (twin_ext_mark_ && pc.mono_ext) out.setMetaValue("spx_mono_ext", 1);   // [twin-rt] read and removed by mergeTwins_
  }

  /// [arena] Assertions on the SHIPPED compactUnreferenced(): the compaction must survive traces
  /// whose container order differs from their arena-offset order, which the previous walk did not.
  ExitCodes selftestArena_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (!ok) { writeLogError_(String("FAIL: ") + what); ++fail; }
      else writeLogInfo_(String("  ok  ") + what);
    };
    // four spans A(4) B(3) C(5) D(2) at offsets 0, 4, 7, 12, every value distinct; span k of trace T
    // holds T*100 + k + 1, with a gap (zero) inside D so interior zeros are exercised
    auto build = [](TraceStore& st, vector<Trace>& tr) {
      st = TraceStore(); tr.clear();
      const vector<pair<uint32_t, int>> spec = {{4, 1}, {3, 2}, {5, 3}, {3, 4}};
      for (const auto& [L, tag] : spec)
      {
        vector<pair<uint32_t, float>> pts;
        for (uint32_t k = 0; k < L; ++k) if (!(tag == 4 && k == 1)) pts.emplace_back(k, (float)(tag * 100 + k + 1));
        tr.push_back(makeSpan(pts, st));
      }
    };
    auto bytes = [](const TraceStore& st, const Trace& t) {
      return vector<float>(st.inten.begin() + t.off, st.inten.begin() + t.off + t.len);
    };
    TraceStore st0; vector<Trace> tr0; build(st0, tr0);
    vector<vector<float>> ref; for (const auto& t : tr0) ref.push_back(bytes(st0, t));
    chk(st0.inten.size() == 15 && tr0[3].len == 3 && tr0[3].npts == 2, "fixture: 15 floats, D spans 3 frames with a gap");

    // 1. container order D,B,A,C (offsets 12,4,0,7), keep D,B,A, drop C. The old container-order walk
    //    moved D to 0 and B to 2 before reading A from [0,4): A came back as D's and B's bytes.
    { vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]}; TraceStore s = st0;
      const Size freed = compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      chk(freed == 1, "1: one profile released");
      chk(s.inten.size() == 10, "1: arena holds exactly the kept spans (4+3+3)");
      chk(bytes(s, v[0]) == ref[3], "1: D survives (interior zero kept)");
      chk(bytes(s, v[1]) == ref[1], "1: B survives");
      chk(bytes(s, v[2]) == ref[0], "1: A survives -- the old walk overwrote it");
      chk(v[3].len == 0 && v[3].npts == 0, "1: C freed (span cleared)"); }
    // 2. offset order, everything needed: unchanged bytes, no move
    { vector<Trace> v = tr0; TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, true)) == 0, "2: nothing released");
      bool ok = s.inten.size() == 15; for (size_t i = 0; i < 4; ++i) ok = ok && bytes(s, v[i]) == ref[i];
      chk(ok, "2: offset-order input is a no-op"); }
    // 3. the freed trace first, then a trimmed span: the trim moved off/len/apex (as trimToSpan does)
    { vector<Trace> v = {tr0[2], tr0[0], tr0[1], tr0[3]}; TraceStore s = st0;
      v[1].off += 1; v[1].len -= 2; v[1].frame0 += 1; v[1].npts = 2;          // A trimmed to [1,3)
      const vector<float> a_trim(ref[0].begin() + 1, ref[0].begin() + 3);
      chk(compactUnreferenced(v, s, vector<bool>{false, true, true, true}) == 1, "3: C released first");
      chk(bytes(s, v[1]) == a_trim && bytes(s, v[2]) == ref[1] && bytes(s, v[3]) == ref[3], "3: trimmed A, B, D survive");
      chk(s.inten.size() == 2 + 3 + 3, "3: arena sized to the trimmed spans"); }
    // 4. nobody needed / a needed trace whose span was already released (len == 0) / duplicate references
    { vector<Trace> v = tr0; TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, false)) == 4 && s.inten.empty(), "4a: none needed -> empty arena");
      vector<Trace> w = tr0; TraceStore s2 = st0; w[0].freeProfile();
      chk(compactUnreferenced(w, s2, vector<bool>(4, true)) == 0 && s2.inten.size() == 11 && bytes(s2, w[1]) == ref[1], "4b: a needed len==0 trace owns no bytes"); }
    // 5. compaction is idempotent
    { vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]}; TraceStore s = st0;
      compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      const vector<float> once = s.inten;
      compactUnreferenced(v, s, vector<bool>{true, true, true, false});
      chk(s.inten == once && bytes(s, v[2]) == ref[0], "5: second compaction is a no-op"); }
    // 6. overlapping spans are refused before any byte moves
    { vector<Trace> v = tr0; TraceStore s = st0; v[1].off = 2;                    // B now overlaps A
      bool threw = false; try { compactUnreferenced(v, s, vector<bool>(4, true)); } catch (const OpenMS::Exception::Precondition&) { threw = true; }
      chk(threw, "6: overlapping spans throw Precondition"); }

    writeLogInfo_(fail ? "[arena] SELFTEST FAILED" : "[arena] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /// [twins] Assertions on mergeTwins_ with hand-built spectra: each side of every clause of the relation,
  /// the cluster closure, the survivor order and the fusion.
  ExitCodes selftestTwins_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (!ok) { writeLogError_(String("FAIL: ") + what); ++fail; }
      else writeLogInfo_(String("  ok  ") + what);
    };
    auto sp = [](double rt, double mz, int z, double im, int niso, const vector<pair<double, float>>& pk) {
      MSSpectrum s; s.setMSLevel(2); s.setRT(rt);
      OpenMS::Precursor p; p.setMZ(mz); if (z > 0) p.setCharge(z); p.setDriftTime(im); s.setPrecursors({p});
      s.setMetaValue("spx_n_isotopes", niso);
      for (const auto& q : pk) s.emplace_back(q.first, q.second);
      return s;
    };
    const vector<pair<double, float>> P{{100.0, 10.0f}, {200.0, 20.0f}};
    auto pair_folds = [&](MSSpectrum b) { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), std::move(b)}; return mergeTwins_(v, 0.005) == 1 && v.size() == 1; };

    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, {{100.0, 10.0f}, {200.0, 20.0f}}), sp(10.0, 500.001, 2, 1.0, 3, {{200.0005, 5.0f}, {300.0, 30.0f}})};
      const size_t gone = mergeTwins_(v, 0.005);
      chk(gone == 1 && v.size() == 1, "a twin pair folds into one spectrum");
      chk(v.size() == 1 && (int)v[0].getMetaValue("spx_n_isotopes") == 3 && v[0].getPrecursors()[0].getMZ() == 500.001, "the member with more isotopes survives, with its own precursor");
      chk(v.size() == 1 && v[0].size() == 3 && v[0][0].getMZ() == 100.0 && v[0][2].getMZ() == 300.0, "the survivor takes the union of the fragments");
      chk(v.size() == 1 && v[0].size() == 3 && fabs(v[0][1].getMZ() - (200.0 * 20.0 + 200.0005 * 5.0) / 25.0) < 1e-9 && v[0][1].getIntensity() == 20.0f,
          "peaks within 10 ppm fuse to the intensity-weighted m/z at the maximum intensity"); }
    chk(!pair_folds(sp(10.0, 500.001, 2, 1.006, 2, P)), "1/K0 beyond the tolerance: untouched");
    chk(pair_folds(sp(10.0, 500.001, 2, 1.004, 2, P)), "1/K0 inside the tolerance: folded");
    chk(!pair_folds(sp(11.1, 500.0, 2, 1.0, 2, P)), "the next MS1 frame: untouched");
    // [twin-rt] the near-frame relation (assembly:twin_rt_frames): |dRT| <= rt_tol, inside one RT cell
    { auto nearFold = [&](double rt_b, double tol, const vector<double>* cells, std::array<size_t, 5>* st = nullptr) {
        vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(rt_b, 500.001, 2, 1.0, 3, P)};
        return mergeTwins_(v, 0.005, false, tol, cells, 3.0, st) == 1 && v.size() == 1; };
      const vector<double> edges{20.0};   // cell 0 owns rt < 17, cell 1 rt >= 17 (the first edge T with (T - rt) > 3)
      std::array<size_t, 5> st{};
      chk(nearFold(11.1, 1.5, nullptr, &st) && st[0] == 0 && st[1] == 1 && st[3] == 0, "the next MS1 frame inside rt_tol: folded, counted as a near-frame link");
      chk(nearFold(11.5, 1.5, nullptr), "|dRT| exactly rt_tol: folded (the bound is inclusive)");
      chk(!nearFold(11.75, 1.5, nullptr), "|dRT| beyond rt_tol: untouched");
      chk(!nearFold(11.1, 0.0, nullptr), "rt_tol 0: the same-frame relation, the next frame untouched");
      chk(nearFold(16.9, 7.0, &edges), "both inside one cell: dRT 6.9 within a 7-s tolerance folds");
      { std::array<size_t, 5> se{};
        vector<MSSpectrum> v{sp(16.5, 500.0, 2, 1.0, 2, P), sp(17.5, 500.001, 2, 1.0, 3, P)};
        chk(mergeTwins_(v, 0.005, false, 1.65, &edges, 3.0, &se) == 0 && v.size() == 2 && se[2] == 1, "a pair across a cell edge stays two spectra (counted)"); }
      { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.002, 2, 1.0, 2, P), sp(11.1, 500.001, 2, 1.004, 2, P), sp(12.2, 500.0, 2, 1.0, 2, P)};
        std::array<size_t, 5> sc{};
        chk(mergeTwins_(v, 0.005, false, 1.2, nullptr, 3.0, &sc) == 3 && v.size() == 1 && sc[0] == 1 && sc[3] == 1 && sc[4] == 1,
            "frames 10, 11.1, 12.2 chain into one cluster (span 2.2 > rt_tol: counted as chained); the same frame alone folds one"); }
      { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(11.1, 500.0, 2, 1.006, 2, P), sp(11.1, 500.0, 3, 1.0, 2, P), sp(11.1, 500.0125, 2, 1.0, 2, P)};
        chk(mergeTwins_(v, 0.005, false, 1.65) == 0 && v.size() == 4, "near frames: 1/K0, charge and 20 ppm still apply"); }
      { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(11.1, 500.0002, 2, 1.0, 2, P), sp(11.6, 500.00002, 2, 1.0, 2, P), sp(10.0, 500.0005, 2, 1.0, 2, P)};
        chk(mergeTwins_(v, 0.005, false, 1.65, nullptr, 3.0, nullptr, 1) == 2 && v.size() == 2,
            "near_bin: across frames only m/z within 0.05 ppm links (10/500.0 with 11.6/500.00002 at 0.04 ppm, not 11.1/500.0002 at 0.4 ppm); one frame keeps the 20-ppm relation"); }
      { auto ext = [&](MSSpectrum s) { s.setMetaValue("spx_mono_ext", 1); return s; };
        vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(11.1, 500.00002, 2, 1.0, 2, P), sp(20.0, 600.0, 2, 1.0, 2, P), ext(sp(21.1, 600.00002, 2, 1.0, 2, P))};
        chk(mergeTwins_(v, 0.005, false, 1.65, nullptr, 3.0, nullptr, 2) == 1 && v.size() == 3 && v[0].getRT() == 10.0 && v[1].getRT() == 11.1,
            "near_mode 2: a one-bin pair across frames folds only when one member's mono came from the v2 extension");
        bool clean = true; for (const auto& q : v) clean = clean && !q.metaValueExists("spx_mono_ext");
        chk(clean, "near_mode 2: the spx_mono_ext marks are removed (never written)"); }
      { vector<MSSpectrum> v{sp(12.0, 500.001, 2, 1.0, 2, P), sp(10.0, 500.0, 2, 1.0, 3, P)};
        mergeTwins_(v, 0.005, false, 2.5);
        chk(v.size() == 1 && v[0].getRT() == 10.0 && v[0].getPrecursors()[0].getMZ() == 500.0, "near frames: the survivor rule is the fold's (more isotopes), with its own RT"); } }
    chk(!pair_folds(sp(10.0, 500.0, 3, 1.0, 2, P)), "another charge: untouched");
    chk(!pair_folds(sp(10.0, 500.0125, 2, 1.0, 2, P)), "25 ppm apart: untouched");
    { vector<MSSpectrum> v{sp(10.0, 500.0, 0, 1.0, 2, P), sp(10.0, 500.0, 0, 1.0, 2, P)};
      chk(mergeTwins_(v, 0.005) == 0 && v.size() == 2, "charge 0 (no charge call): untouched"); }
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.0075, 2, 1.0, 2, P), sp(10.0, 500.015, 2, 1.0, 2, P)};
      chk(mergeTwins_(v, 0.005) == 2 && v.size() == 1, "a chain 15 + 15 ppm is one cluster although its ends are 30 ppm apart"); }
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.001, 2, 1.0, 2, P)};
      mergeTwins_(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.0, "equal isotopes and peaks: the earliest survives"); }
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.001, 2, 1.0, 2, {{100.0, 1.0f}, {150.0, 1.0f}, {250.0, 1.0f}})};
      mergeTwins_(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.001, "equal isotopes: more peaks survives"); }
    { vector<MSSpectrum> v{sp(5.0, 300.0, 2, 0.8, 2, P), sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 700.0, 3, 1.1, 2, P), sp(10.0, 500.001, 2, 1.0, 3, P), sp(12.0, 400.0, 2, 0.9, 2, P)};
      mergeTwins_(v, 0.005);
      chk(v.size() == 4 && v[0].getRT() == 5.0 && v[1].getPrecursors()[0].getMZ() == 700.0 && v[2].getPrecursors()[0].getMZ() == 500.001 && v[3].getRT() == 12.0,
          "every survivor keeps its own place (the later twin here), the folded member leaves no gap"); }
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.0, 2, 1.001, 2, P)};
      chk(mergeTwins_(v, 0.0) == 1 && v.size() == 2, "tolerance 0 folds identical-1/K0 twins only"); }
    chk(pair_folds(sp(10.0, 500.00975, 2, 1.0, 2, P)), "19.5 ppm apart: folded (the 20-ppm relation, pinned from below)");
    chk(!pair_folds(sp(10.0, 500.01025, 2, 1.0, 2, P)), "20.5 ppm apart: untouched (pinned from above)");
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 2, P), sp(10.0, 500.001, 2, 1.25, 2, P)};
      chk(mergeTwins_(v, 0.25) == 1, "|d1/K0| exactly the tolerance: folded (the bound is inclusive)"); }
    { vector<MSSpectrum> v{sp(10.0, 500.001, 2, 1.0, 2, P), sp(10.0, 500.0, 2, 1.0, 2, P)};
      chk(mergeTwins_(v, 0.005) == 1 && v.size() == 1, "the higher m/z first in the input: still folded (the backward scan)"); }
    { vector<MSSpectrum> v{sp(10.0, 500.0, 2, 1.0, 3, {{100.0, 1.0f}, {200.0, 1.0f}}), sp(10.0, 500.001, 2, 1.0, 2, {{110.0, 1.0f}, {210.0, 1.0f}, {310.0, 1.0f}})};
      mergeTwins_(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.0, "isotopes outrank peaks: 3 isotopes / 2 peaks beats 2 isotopes / 3 peaks"); }
    { vector<MSSpectrum> v{sp(20.0, 400.0, 2, 0.9, 2, {{400.0, 1.0f}, {400.003, 2.0f}})};
      chk(mergeTwins_(v, 0.005) == 0 && v.size() == 1 && v[0].size() == 2 && v[0][0].getMZ() == 400.0 && v[0][1].getMZ() == 400.003
          && v[0][0].getIntensity() == 1.0f && v[0][1].getIntensity() == 2.0f, "a spectrum without a twin keeps its peaks bit for bit, 7.5 ppm apart or not"); }
    // the fusion anchor is the RUNNING intensity-weighted mean, at 10 ppm of it
    auto fused = [&](const vector<pair<double, float>>& a, const vector<pair<double, float>>& b) {
      vector<MSSpectrum> v{sp(30.0, 600.0, 2, 1.0, 3, a), sp(30.0, 600.0, 2, 1.0, 2, b)};
      mergeTwins_(v, 0.005); return v.size() == 1 ? v[0].size() : (size_t)0; };
    chk(fused({{200.0, 1.0f}}, {{200.0021, 1.0f}}) == 2, "fusion: 10.5 ppm apart stays two peaks");
    chk(fused({{200.0, 1.0f}, {200.0038, 1.0f}}, {{200.0019, 1.0f}}) == 2, "fusion: equal weights, the third peak is 14.25 ppm from the running mean: two peaks");
    chk(fused({{200.0, 1.0f}, {200.0038, 1.0f}}, {{200.0019, 100.0f}}) == 1, "fusion: a heavy middle peak pulls the mean to 9.6 ppm of the third: one peak (not anchored on the first)");
    { MSSpectrum a = sp(40.0, 700.0, 3, 1.1, 3, P), b = sp(40.0, 700.002, 3, 1.101, 2, P);
      OpenMS::Precursor pa = a.getPrecursors()[0]; pa.setIsolationWindowLowerOffset(7.0); pa.setIsolationWindowUpperOffset(18.0); a.setPrecursors({pa});
      a.setMetaValue("spx_guessed", 0); b.setMetaValue("spx_guessed", 1);
      vector<MSSpectrum> v{a, b};
      mergeTwins_(v, 0.005);
      const auto& q = v[0].getPrecursors()[0];
      chk(v.size() == 1 && v[0].getRT() == 40.0 && q.getMZ() == 700.0 && q.getCharge() == 3 && q.getDriftTime() == 1.1
          && q.getIsolationWindowLowerOffset() == 7.0 && q.getIsolationWindowUpperOffset() == 18.0 && (int)v[0].getMetaValue("spx_guessed") == 0,
          "the survivor keeps its own RT, precursor, 1/K0, isolation offsets and userParams"); }
    // [recap] the keyed fold (assembly:post_fusion_cap) and the cap after it
    auto keyed = [](MSSpectrum s, const vector<float>& k) {
      OpenMS::DataArrays::FloatDataArray a; a.setName("spx_cap_key"); a.assign(k.begin(), k.end()); s.getFloatDataArrays().push_back(std::move(a)); return s; };
    { const vector<pair<double, float>> A{{100.0, 10.0f}, {200.0, 20.0f}}, B{{200.0005, 5.0f}, {300.0, 30.0f}};
      vector<MSSpectrum> u{sp(10.0, 500.0, 2, 1.0, 2, A), sp(10.0, 500.001, 2, 1.0, 3, B)};
      vector<MSSpectrum> v{keyed(sp(10.0, 500.0, 2, 1.0, 2, A), {1.0f, 2.0f}), keyed(sp(10.0, 500.001, 2, 1.0, 3, B), {7.0f, 0.5f})};
      mergeTwins_(u, 0.005); mergeTwins_(v, 0.005, true);
      bool same = u.size() == 1 && v.size() == 1 && u[0].size() == v[0].size();
      for (Size i = 0; same && i < u[0].size(); ++i) same = u[0][i].getMZ() == v[0][i].getMZ() && u[0][i].getIntensity() == v[0][i].getIntensity();
      chk(same, "keyed fold: the peaks of the unkeyed fold, bit for bit");
      const auto& k = v[0].getFloatDataArrays();
      chk(k.size() == 1 && k[0].size() == 3 && k[0][0] == 1.0f && k[0][1] == 7.0f && k[0][2] == 0.5f, "keyed fold: a fused peak keeps its members' max key");
      const std::array<size_t, 3> r = capAfterFusion_(v, 2);
      chk(r[0] == 1 && r[1] == 1 && r[2] == 3 && v[0].size() == 2 && v[0][0].getMZ() == 100.0 && v[0][1].getMZ() == u[0][1].getMZ()
          && v[0].getFloatDataArrays().empty(), "post-fusion cap: the top keys survive in m/z order, the key array is dropped"); }
    { vector<MSSpectrum> v{keyed(sp(10.0, 500.0, 2, 1.0, 2, {{100.0, 5.0f}, {200.0, 9.0f}, {300.0, 9.0f}}), {1.0f, 1.0f, 1.0f}),
                           keyed(sp(12.0, 500.0, 2, 1.0, 2, {{100.0, 5.0f}}), {3.0f})};
      const std::array<size_t, 3> r = capAfterFusion_(v, 1);
      chk(r[0] == 1 && r[1] == 2 && v[0].size() == 1 && v[0][0].getMZ() == 200.0, "post-fusion cap ties: the higher intensity, then the lower m/z");
      chk(v[1].size() == 1 && v[1][0].getMZ() == 100.0 && v[1].getFloatDataArrays().empty(), "post-fusion cap: a spectrum within the cap keeps its peaks, without the key array"); }
    writeLogInfo_(fail ? "[twins] SELFTEST FAILED" : "[twins] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /// Export-time calibration: the reported m/z is computed once, here, from the bin and its apex frame's
  /// factor. tof 0 (the OpenMS detector has no bin) keeps the cached m/z.
  static double exportMz_(uint32_t tof, double b, double cached)
  {
    if (tof == 0 || !(b > 0.0) || !tofAxis().ok) return cached;
    const double v = tofAxis().mzOf(tof, b);
    return std::isfinite(v) && v > 0.0 ? v : cached;
  }

  /// [trim] malloc_trim(0), logged with the RSS before and after and how long it took.
  void trimLogged_(const String& what) const
  {
    const double t0 = phase_clock_(); const long r0 = rss_mb_();
#ifdef __GLIBC__
    malloc_trim(0);
#endif
    writeLogInfo_("[trim] " + what + ": " + String(r0) + " -> " + String(rss_mb_()) + " MB RSS in " + String(phase_clock_() - t0) + " s" + mem_());
  }

  ExitCodes main_(int argc, const char** argv) override
  {
    phase_clock_();   // [perf-instr] fix the epoch (a first-call static) before loading
    LedRun _ledger;   // [ledger] DIASPEXTRACT_LEDGER=<file>: the per-structure memory profile

    // -threads not given: use every core. The command line is scanned, so an explicit -threads 1 stays serial.
    int n_threads_req = getIntOption_("threads");
    {
      bool given = false;
      for (int i = 1; i < argc && !given; ++i) given = (String(argv[i]) == "-threads");
      if (!given && n_threads_req == 1)
      {
        n_threads_req = std::max(1, (int)std::thread::hardware_concurrency());
        omp_set_num_threads(n_threads_req);
        writeLogInfo_("-threads not given: using all " + String(n_threads_req) + " cores. Pass "
                      "-threads <n> to limit it, which is what you want on a shared machine.");
      }
    }
    const String in = getStringOption_("in");
    const String out = getStringOption_("out");
    const double delta_im = getDoubleOption_("gate:delta_im");
    // [frag-im] the gate's fragment-precursor 1/K0 window, apart from the tracing tolerance delta_im (negative = delta_im)
    const double gate_dim = getDoubleOption_("gate:fragment_delta_im") >= 0.0 ? getDoubleOption_("gate:fragment_delta_im") : delta_im;
    const double delta_rt = getDoubleOption_("gate:delta_rt");
    if (getFlag_("diag:selftest_arena")) return selftestArena_();
    if (getFlag_("diag:selftest_twins")) return selftestTwins_();
    log_overlap_ = (getStringOption_("gate:coelution") == "logoverlap");
    im_weight_sigma_ = getDoubleOption_("assembly:im_weight_sigma");
    mono_guard_ = getDoubleOption_("charge:mono_averagine_guard");
    mono_select_ = getFlag_("charge:mono_averagine_select");
    {
      const String mc = getStringOption_("charge:mono_profile_check");
      mono_check_ = (mc == "true") ? 1 : (mc == "v2") ? 2 : 0;
    }
    mono_check_min_corr_ = getDoubleOption_("charge:mono_profile_min_corr");
    mono_check_max_dim_ = getDoubleOption_("charge:mono_profile_max_dim");
    mono_pos_min_ratio_ = getDoubleOption_("charge:mono_position_min_ratio");
    mono_pos_im_box_ = getDoubleOption_("charge:mono_position_im_box");
    mono_pos_extend_ = getDoubleOption_("charge:mono_position_extend");
    mono_pos_iterate_ = getStringOption_("charge:mono_position_iterate") == "true";
    iso_closest_ = getStringOption_("charge:iso_partner") == "closest";
    var_support_ = getFlag_("gate:variance_support");
    short_conc_ = (getStringOption_("gate:short_traces") == "concentric");
    short_center_tol_ = getDoubleOption_("gate:short_center_tol");
    short_min_pts_ = getIntOption_("gate:short_min_points");
    coloc_ = (getStringOption_("gate:colocation") == "ranges");
    coloc_overlap_ = getDoubleOption_("gate:colocation_overlap");
    shrink_ = (getStringOption_("gate:shrinkage") == "fisher");
    shrink_r0_ = getDoubleOption_("gate:shrink_prior_r");
    shrink_sd_ = getDoubleOption_("gate:shrink_prior_sd");
    shrink_min_pts_ = getIntOption_("gate:shrink_min_points");
    adaptive_ = (getStringOption_("gate:adaptive") == "faint");
    adaptive_auto_ = (getStringOption_("gate:adaptive") == "auto");
    auto_faint_frac_ = getDoubleOption_("gate:auto_faint_fraction");
    cap_rank_ = getStringOption_("assembly:cap_rank") == "umpire" ? 1 : getStringOption_("assembly:cap_rank") == "evidence" ? 2 : 0;
    im_slice_merge_ = getDoubleOption_("trace:im_slice_merge");
    merge_same_bin_ = getStringOption_("assembly:merge_same_bin") == "true";
    frag_int_ = getStringOption_("assembly:frag_intensity") == "core" ? 1 : getStringOption_("assembly:frag_intensity") == "sum" ? 2 : 0;
    if (frag_int_ == 1 && getIntOption_("assembly:frag_core_frames") % 2 == 0)   // used by 'core' only
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "assembly:frag_core_frames must be odd: the core is centred on the precursor's apex frame");
    frag_core_half_ = getIntOption_("assembly:frag_core_frames") / 2;
    post_fusion_cap_ = getStringOption_("assembly:post_fusion_cap") == "true";
    foreign_w_ = getDoubleOption_("assembly:foreign_precursor_weight");
    foreign_ppm_ = getDoubleOption_("assembly:foreign_precursor_ppm");
    foreign_dim_ = getDoubleOption_("assembly:foreign_precursor_dim");
    foreign_min_iso_ = getIntOption_("assembly:foreign_precursor_min_iso");
    fuse_ppm_ = getDoubleOption_("assembly:fragment_fuse_ppm");
    fuse_linkage_ = getStringOption_("assembly:fragment_fuse_mode") == "linkage";
    faint_frames_ = getIntOption_("gate:faint_max_frames");
    faint_drt_ = getDoubleOption_("gate:faint_delta_rt");
    faint_cp_ = getDoubleOption_("gate:faint_corr_power");
    anchor_rescue_ = getStringOption_("gate:anchor_rescue") == "true";
    anchor_min_count_ = getIntOption_("gate:anchor_min_count");
    anchor_min_corr_ = getDoubleOption_("gate:anchor_min_corr");
    anchor_rescue_corr_ = getDoubleOption_("gate:anchor_rescue_corr");
    anchor_max_rescued_ = getIntOption_("gate:anchor_max_rescued");
    if (shrink_ && log_overlap_)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "gate:shrinkage fisher needs a correlation: it does not apply to gate:coelution logoverlap");
    corr_power_ = getDoubleOption_("assembly:corr_power");
    const double min_corr = getDoubleOption_("gate:min_correlation");
    const int min_corr_pts = getIntOption_("gate:min_correlation_points");
    const Size min_frags = (Size)getIntOption_("assembly:min_fragments");
    const Size max_frags = (Size)getIntOption_("assembly:max_fragments");
    const double twin_im = getDoubleOption_("assembly:twin_im_tolerance");
    const int twin_rt_frames = getIntOption_("assembly:twin_rt_frames");   // [twin-rt] -1 off, 0 auto, n frames
    const String twin_rt_match = getStringOption_("assembly:twin_rt_match");
    const int twin_rt_mode = twin_rt_match == "ppm" ? 0 : twin_rt_match == "bin" ? 1 : 2;
    const int max_charge = getIntOption_("max_charge");
    mzEstimator() = getStringOption_("trace:mz_estimator");   // before any toTrace() call
    const double mass_ppm = getDoubleOption_("trace:mass_error_ppm");

    PeakMap exp;
    // Accept mzML, mzPeak, or a Bruker .d directory; FileHandler auto-detects the format
    // (mzPeak -> MzPeakFile; .d -> BrukerTimsFile, requires WITH_OPENTIMS).
    PeakMap ms1_map;
    // [ms1-prune] one run-level prune (see Ms1Prune), used by the streaming consumer and the resident dispatch alike
    Ms1Prune ms1_prune;
    ms1_prune.on = (getStringOption_("perf:ms1_prune") == "true");
    ms1_prune.noise = getDoubleOption_("trace:noise_threshold_int");
    ms1_prune.hop = mass_ppm * 1e-6 * kBandHaloPpmMul;
    { const char* e = std::getenv("DIASPEXTRACT_MS1_PRUNE_NO_WITNESS"); ms1_prune.no_witness = e && *e && *e != '0'; }
    { const char* e = std::getenv("DIASPEXTRACT_MS1_PRUNE_NO_CHAIN"); ms1_prune.no_chain = e && *e && *e != '0'; }   // its own switch: check 19b
    map<WinKey, PeakSlab> ms2_by_window;
    CompactStats cstat;
    const bool stream_load = (getStringOption_("perf:stream_load") == "true");
    const bool trim_on = (getStringOption_("perf:malloc_trim") == "true");
    bool streamed = false;
    // Output format, resolved before the load: it decides the tiling below (mzPeak output runs as one tile), and a
    // plain build refuses -out_type mzpeak here in milliseconds. The extension is authoritative, because that is what
    // a user typing `-out pseudo.mzML` means; -out_type overrides it; mzPeak is the default when neither says.
    FileTypes::Type out_type = FileTypes::MZML;
    const String out_type_opt = getStringOption_("out_type");
#ifdef DIASPEXTRACT_WITH_MZPEAK
    if (!out_type_opt.empty())
      out_type = (out_type_opt == "mzML") ? FileTypes::MZML : FileTypes::MZPEAK;
    else
    {
      const FileTypes::Type by_name = FileHandler::getTypeByFileName(out);
      out_type = (by_name == FileTypes::MZML || by_name == FileTypes::MZPEAK) ? by_name
                                                                             : FileTypes::MZPEAK;
    }
    if (out_type == FileTypes::MZPEAK)
      writeLogInfo_("Writing mzPeak. NOTE: DDA search engines read mzML, not mzPeak -- if this file "
                    "is going straight into a search, write .mzML instead (or -out_type mzML).");
#else
    // Built without mzPeak support: mzML is the only container this binary can write. Refuse an
    // explicit request for mzPeak rather than silently writing something else.
    if (!out_type_opt.empty() && out_type_opt != "mzML")
      throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "-out_type mzpeak needs a build with mzPeak support (-DMZPEAK_ROOT=...); this binary "
            "can only write mzML.");
#endif
    // [3a] trace:band_edges=acquisition needs the tdf's GlobalMetadata bounds: refuse BEFORE the
    // load when no candidate tdf provides them (an mzPeak archive is checked after its embedded
    // tdf is recovered, below). Milliseconds against minutes of loading.
    const bool band_edges_acq = (getStringOption_("trace:band_edges") == "acquisition");
    const bool integer_req = (getStringOption_("trace:detector") == "integer");
    const char* tr_env = std::getenv("DIASPEXTRACT_TILE_RESIDENT");
    const bool force_resident = tr_env && String(tr_env) != "0" && String(tr_env) != "";
    if (!band_edges_acq && integer_req && stream_load && FileHandler::getTypeByFileName(in) == FileTypes::BRUKER_TDF && !force_resident)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "trace:band_edges=slab needs the whole window resident; the .d streaming source reads it per tile. Use trace:band_edges=acquisition, or DIASPEXTRACT_TILE_RESIDENT=1 for the resident source.");
    std::vector<std::string> tdf_cand{in + "/analysis.tdf"};   // the calibration sources, in the order the load tries them
    if (const char* sc = std::getenv("DIASPEXTRACT_MZPEAK_TDF")) tdf_cand.push_back(sc);
    if (band_edges_acq && integer_req)
    {
      // [3a] The axis is selected HERE, in the order the load below would use it (the first
      // candidate whose calibration loads is THE source), so the bounds that are checked belong to
      // the tdf that will actually calibrate the run -- and the check costs a two-row query rather
      // than the whole load. No candidate at all is NOT an error: the integer detector then falls
      // back to the OpenMS detector further down, loudly, and traces no bands to partition.
      for (const std::string& c : tdf_cand)
      {
        String w_axis;
        if (!loadTofAxis(c, w_axis)) continue;   // not a calibration source; try the next
        if (!(tofAxis().n_bins > 1))
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "trace:band_edges=acquisition needs MzAcqRangeLower/Upper and DigitizerNumSamples from the tdf's GlobalMetadata; the calibration source " + c + " does not provide them (" + tofAxis().bounds_why + "). Use -trace:band_edges slab, or -trace:detector openms.");
        writeLogInfo_("flight-time axis from " + c + " (preflight): " + String(tofAxis().b_by_frame.size()) + " frame factors; digitizer bins " + String(tofAxis().n_bins) + ", acquisition m/z " + String(tofAxis().acq_lo) + "-" + String(tofAxis().acq_hi));
        break;
      }
    }
    // the .d loader's configuration; set here because the dnoise setup below predicts the loader's m/z model from it
    BrukerTimsFile::Config cfg;
    cfg.ms1_n_neighbors     = getIntOption_("trace:native_ms1_neighbors");
    // [dnoise] dnoise v0.1.0's MS1 path on the raw frames (DnoiseRun), set up once for the run: every refusal the options
    // and the tdf can decide comes here, in milliseconds, before the load
    DnoiseRun dnoise;
    const bool dnoise_req = (getStringOption_("dnoise:ms1") == "true");
    if (dnoise_req && FileHandler::getTypeByFileName(in) != FileTypes::BRUKER_TDF)
      writeLogInfo_("[dnoise] MS1 denoising skipped: the input is not a Bruker .d");
    else if (dnoise_req)
    {
      auto refuse = [](const String& why) {
        return Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "dnoise:ms1 " + why + " -dnoise:ms1 false runs without MS1 denoising."); };
      if (cfg.ms1_n_neighbors > 0)
        throw refuse("denoises each raw MS1 frame, and trace:native_ms1_neighbors > 0 sums neighbouring frames before DIAspeXtract sees them.");
      // The loader's m/z model, predicted: under the default AUTO strategy an SDK path (the Config's, else
      // OPENMS_BRUKER_SDK_PATH) makes openTimsDataHandle calibrate m/z through the Bruker SDK, which the TOF recovery
      // cannot invert. An SDK that then fails to load would leave the table model, so this can refuse a run that would
      // have recovered: fail-closed. checkLoader() confirms the loader's actual choice during the load.
      { const char* env_sdk = std::getenv("OPENMS_BRUKER_SDK_PATH");
        const std::string sdk = !cfg.bruker_sdk_path.empty() ? cfg.bruker_sdk_path : std::string(env_sdk ? env_sdk : "");
        if (!sdk.empty())
          throw refuse("recovers each point's TOF bin through the tdf's MzCalibration table model, but with a Bruker SDK path set (" + String(sdk)
                       + ") the loader calibrates m/z through the Bruker SDK: unset OPENMS_BRUKER_SDK_PATH to denoise."); }
      String why_axis;
      if (!tofAxis().ok && !loadTofAxis(in + "/analysis.tdf", why_axis))   // the TOF recovery runs during the load
        throw refuse("recovers each point's TOF bin through the tdf's MzCalibration table model, which is not usable here (" + why_axis + ").");
      spx::dnoise::Params& p = dnoise.p;
      p.mz_half_width = (uint32_t)getIntOption_("dnoise:mz_half_width");
      p.min_feature_length = (size_t)getIntOption_("dnoise:min_feature_length");
      p.max_internal_gap = (size_t)getIntOption_("dnoise:max_internal_gap");
      p.iterations = (size_t)getIntOption_("dnoise:iterations");
      p.min_window_intensity = (uint64_t)getIntOption_("dnoise:min_window_intensity");
      p.min_feature_intensity = (uint64_t)getIntOption_("dnoise:min_feature_intensity");
      p.halo = (getStringOption_("dnoise:halo") == "true");
      p.halo_peak_fraction = getDoubleOption_("dnoise:halo_peak_fraction");
      p.halo_mz_idx_half_width = (uint32_t)getIntOption_("dnoise:halo_mz_idx_half_width");
      p.halo_scan_half_width = (size_t)getIntOption_("dnoise:halo_scan_half_width");
      p.dia_ms1_window = (getStringOption_("dnoise:dia_ms1_window") == "true");
      p.dia_ms1_mz_pad = getDoubleOption_("dnoise:dia_ms1_mz_pad");
      p.dia_ms1_im_pad = getDoubleOption_("dnoise:dia_ms1_im_pad");
      const std::string why = dnoise.setup(in + "/analysis.tdf");
      if (!why.empty()) throw refuse("cannot be exact on this input: " + why + ".");
      writeLogInfo_("[dnoise] " + String(spx::dnoise::kPortedFrom) + " MS1 path on the raw frames, before the pick: " + dnoise.params()
                    + "; " + String(dnoise.im_asc.size()) + " TimsCalibration row(s) over " + String(dnoise.tdf.meta.max_num_scans)
                    + " scans (MAX(Frames.NumScans), the gate's); scan tables per row over its MS1 frames' scans: " + dnoise.scanTables());
    }
    String frame_table_src = "slab";   // [frames] where the frozen frame tables came from: loader | mzpeak | slab
    // [3b] the .d streaming SOURCE: pass 1 reads MS1 only (the loader exports the whole frame
    // table regardless), pass 2 reads each tile's MS2 frames through the loader's RT range. The
    // resident source (mzPeak, mzML, perf:stream_load=false, or DIASPEXTRACT_TILE_RESIDENT=1 for an
    // A/B) holds the whole run and slices it, as before.
    bool tile_source = false;
    PickCompactConsumer::FrameTables run_tables;   // the reader's frozen tables, kept for pass 2
    PeakPickerIM spicker;
    {
      Param sp = spicker.getParameters();
      sp.setValue("pickIMCluster:im_tolerance_cluster", delta_im);
      sp.setValue("pickIMCluster:ppm_tolerance_cluster", mass_ppm);
      spicker.setParameters(sp);
    }

    if (stream_load)
    {
      // [stream] frame-by-frame: pick + compact on arrival, never hold the whole run.
      // Also the ONLY place the reader's NATIVE pre-centroid frame aggregation is reachable.
      PickCompactConsumer consumer(spicker, ms2_by_window, ms1_map, cstat);
      consumer.ms1_prune = &ms1_prune;
      consumer.dnoise = &dnoise;
#ifdef DIASPEXTRACT_WITH_MZPEAK
      const bool in_mzpeak = (FileHandler::getTypeByFileName(in) == FileTypes::MZPEAK);
#else
      const bool in_mzpeak = false;
#endif
      // [3b] the streaming source needs the Bruker loader: a .d only, never an mzML or an mzPeak
      // archive (the resident path reads those whole). DIASPEXTRACT_TILE_RESIDENT set to any non-empty
      // value other than 0 (parsed once, above) forces the resident source for a same-build A/B.
      tile_source = !in_mzpeak && !force_resident && FileHandler::getTypeByFileName(in) == FileTypes::BRUKER_TDF;
      if (tile_source)
      {
        // pass 1: MS1 only; lo > hi selects no MS2 frame, the loader exports the table anyway
        cfg.dia_ms2_rt_lo = std::numeric_limits<double>::infinity();
        cfg.dia_ms2_rt_hi = -std::numeric_limits<double>::infinity();
      }
      // [frames] The frozen frame tables come from the READER's metadata, keyed with the reader's
      // own precursor expressions so the keys equal the consumer's: the .d loader exports its
      // frame table and windows (patch: lastDIAFrameTable/lastDIAWindows; frames with no raw peaks
      // are skipped for every window, as the loop skips them), the mzPeak reader its cached sweep
      // (a spectrum with number_of_peaks == 0 yields only empty window spectra). An mzML is
      // resident and untiled: no table, the slab is its own.
      String table_src_if_ok = "slab";   // becomes the stamp only after the streamed load and the invariant succeed
      if (in_mzpeak)
      {
#ifdef DIASPEXTRACT_WITH_MZPEAK
        table_src_if_ok = "mzpeak";
        if (band_edges_acq && integer_req)
        {
          // [3a] the sweep recovers the embedded tdf's calibration and bounds; it is cached, so
          // the load below reuses it -- the refusal comes before any peak is decoded
          const auto m = spx::mzPeakRunMeta(in);
          if (!m->exact || !(m->exact->n_bins > 1))
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "trace:band_edges=acquisition needs DigitizerNumSamples and MzAcqRangeLower/Upper from the archive's (or sidecar) analysis.tdf; not provided (" + String(m->exact ? m->exact->bounds_why : "no exact calibration") + ").");
        }
        consumer.tablefn = [&]() {
          PickCompactConsumer::FrameTables t;
          const auto m = spx::lastMzPeakMeta();
          if (!m) return t;
          for (size_t i : m->ms2)
          {
            if (m->n_peaks[i] == 0) continue;
            for (const auto& w : m->wins[i])
            {
              const double c = 0.5 * (w.lo + w.hi), off = 0.5 * (w.hi - w.lo);   // frameToSpectra_'s precursor
              auto& tab = t[winKey(c - off, c + off, 0)];
              tab.rt.push_back(m->rt[i]); tab.id.push_back(m->frame_id[i] < 0 ? 0u : (uint32_t)m->frame_id[i]);
            }
          }
          return t; };
#endif
      }
      else
      {
        table_src_if_ok = "loader";
        consumer.tablefn = [&]() {
          PickCompactConsumer::FrameTables t;
          const auto& ft = BrukerTimsFile::lastDIAFrameTable();
          for (const auto& w : BrukerTimsFile::lastDIAWindows())
          {
            const auto g = ft.find(w.window_group);
            if (g == ft.end()) continue;
            const double c = w.mz_center, off = w.mz_width / 2.0;                 // BrukerTimsFile's precursor
            const WinKey k = winKey(c - off, c + off, (uint32_t)w.window_group);
            if (t.count(k))   // two table rows with one (m/z, group) key -- the consumer would merge their frames into one slab
              throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[frames] two DIA windows share one (m/z lo, m/z hi, group) key; DIAspeXtract keys windows by isolation window and group", String(c - off) + "-" + String(c + off) + " group " + String(w.window_group));
            auto& tab = t[k];
            for (const auto& r : g->second) if (r.num_peaks > 0) { tab.rt.push_back(r.time); tab.id.push_back(r.frame_id); }
          }
          return t; };
      }
      writeLogInfo_("[outer] start-up before the stream reader: " + String(phase_clock_()) + " s");
      try
      {
        if (in_mzpeak)
        {
#ifdef DIASPEXTRACT_WITH_MZPEAK
          spx::loadMzPeakStreaming(in, consumer, n_threads_req);   // consumer counts frames itself
#else
          throw Exception::NotImplemented(__FILE__, __LINE__, "streaming .mzpeak input needs a build with -DMZPEAK_ROOT");
#endif
        }
        else BrukerTimsFile().loadDIAStreaming(in, consumer, cfg);
        consumer.finish();                                          // the last partial batch
        streamed = true;
        if (ms1HistOn_())
        { const auto& h = ms1Hist_(); long long tot = 0; for (const auto& v : h) tot += v.load();
          auto pc = [&](long long v) { return String(tot > 0 ? 100.0 * (double)v / (double)tot : 0.0); };
          writeLogInfo_("[ms1-hist] picked MS1 peaks " + String(tot) + " by intensity: <100 " + String(h[0].load()) + " (" + pc(h[0].load()) + "%), 100-300 "
                        + String(h[1].load()) + " (" + pc(h[1].load()) + "%), 300-1000 " + String(h[2].load()) + " (" + pc(h[2].load()) + "%), 1000-10000 "
                        + String(h[3].load()) + " (" + pc(h[3].load()) + "%), >=10000 " + String(h[4].load()) + " (" + pc(h[4].load()) + "%)"); }
        phaseAdd_("LOAD(stream)", 0.0, 0.0, false);                // [perf-instr] start-to-here
        writeLogInfo_("[perf-load] pick(flush) wall=" + String(flush_wall_(), 1) + " s inside LOAD(stream)");
        writeLogInfo_("[outer] MS1 pick " + consumer.ms1_pick_mode_);
        if (std::getenv("DIASPEXTRACT_LOAD_ONLY"))
        { if (tile_source) writeLogInfo_("[perf-load] DIASPEXTRACT_LOAD_ONLY under the .d streaming source measures PASS 1 (MS1) only; the MS2 frames are read per tile in pass 2.");
          report_phases_(phase_clock_()); return EXECUTION_OK; }   // [perf-load] decode/hand-off/pick split only
        writeLogInfo_("[stream] frame-by-frame load done: " + String(consumer.frames_seen)
                      + " frames, MS1 " + String(ms1_map.size()) + ", windows "
                      + String(ms2_by_window.size()) + clk_() + rss_() + mem_());
      }
      catch (const InputRefusal&) { throw; }   // a refusal of the input, not a reader failure: never fall back onto partial slabs
      catch (const Exception::BaseException& e)
      {
        writeLogWarn_(String("[stream] streaming load failed (") + e.getName()
                      + "); falling back to loadExperiment. " + e.getMessage());
      }
      // [3b] the source is COMMITTED only by a successful pass 1: a fallback to the resident load
      // reads the whole run, so the tile loop must slice it rather than call the loader again
      if (!streamed && tile_source)
      { tile_source = false; writeLogWarn_("[3b] the streaming source did not read pass 1; the resident source takes over (the whole run is loaded)."); }
      // [frames] the invariant the frozen tables rest on, checked OUTSIDE the recoverable catch: a
      // delivered window's slab is its table (rt and id, frame for frame) after padding
      if (streamed && tile_source && consumer.tables().empty())
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[3b] the .d loader exported no frame table; the tiled source cannot read this input", in);
      if (streamed && !consumer.tables().empty())
      {
        size_t n_pad = 0, n_del = 0, n_absent = 0;
        for (const auto& kv : consumer.tables())
        {
          const String wl = String(kv.first[0] / 100.0) + "-" + String(kv.first[1] / 100.0);
          const auto w = ms2_by_window.find(kv.first);
          if (w == ms2_by_window.end()) { ++n_absent; if (!tile_source) writeLogInfo_("[frames] window " + wl + ": absent (no frame delivered), " + String(kv.second.rt.size()) + " in the table"); continue; }
          if (w->second.rt != kv.second.rt || w->second.frame_id != kv.second.id)
            throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[frames] a window's slab is not its frozen frame table after padding", wl);
          n_pad += kv.second.padded; n_del += kv.second.delivered;
          writeLogInfo_("[frames] window " + wl + ": delivered " + String(kv.second.delivered) + ", padded " + String(kv.second.padded) + " of " + String(kv.second.rt.size()));
        }
        frame_table_src = table_src_if_ok;
        run_tables = consumer.tables();
        writeLogInfo_("[frames] source=" + frame_table_src + " windows=" + String(consumer.tables().size()) + " delivered=" + String(n_del) + " padded=" + String(n_pad) + " absent=" + String(n_absent)
                      + (tile_source ? " (pass 1: MS1 only; the windows' frames are read per tile)" : ""));
#ifdef DIASPEXTRACT_WITH_MZPEAK
        if (in_mzpeak && spx::lastMzPeakMeta())
          writeLogInfo_("[frames] mzPeak number_of_peaks: " + String(spx::lastMzPeakMeta()->n_peaks_missing_ms2) + " MS2 spectra lack it and count as present"
                        + (spx::lastMzPeakMeta()->n_peaks_missing_ms2 ? " (a frame with no peaks among them would be padded)" : "; a recorded zero is skipped like the .d loop's"));
#endif
      }
    }

    if (!streamed)
    {
      // A failed stream attempt may already have delivered part of the run -- MS1 spectra (sampled, pruned and
      // denoised), window slabs, compaction counts and their ledger charges. The resident load below reads the WHOLE
      // run again and only appends, so discard all of it first: kept, the delivered frames were loaded twice (silently
      // wrong with the prune off, a refusal with it on, "not in retention-time order" once partial slabs existed).
      ms1_map = PeakMap();
      ms2_by_window.clear();
      cstat = CompactStats();
      ms1_prune.reset();
      for (auto& v : ms1Hist_()) v.store(0, std::memory_order_relaxed);   // [ms1-hist] the failed attempt's picked peaks
      dnoise.reset();   // [dnoise] counts and frame marks: the resident pick below filters every frame again
      ledSet_(LC_MS1MAP, 0); ledSet_(LC_SLAB, 0);
      Phase _ph("LOAD(full)");
      FileHandler().loadExperiment(in, exp,
#ifdef DIASPEXTRACT_WITH_MZPEAK
                                   {FileTypes::MZML, FileTypes::MZPEAK, FileTypes::BRUKER_TDF},
#else
                                   {FileTypes::MZML, FileTypes::BRUKER_TDF},
#endif
                                   log_type_);
    }

    //-------------------------------------------------------------
    // Input normalization / validation [C-4]
    //-------------------------------------------------------------
    // The IM check applies to the resident load only: a streamed run leaves `exp` empty (ms1_map.empty() guards it below).
    if (!streamed)
    {
      IMFormat imf = IMTypes::determineIMFormat(exp, 1);
      if (imf == IMFormat::NONE)
      {
        writeLogError_("Error: input has no ion-mobility data at MS1. DIAspeXtract requires diaPASEF (IM DIA) input.");
        return ILLEGAL_PARAMETERS;
      }
    }
    // deliberate: assumes 1/K0 (VSSC). ms/CCS conversion and FAIMS rejection are Phase 2 [M-2].

    // Centroid IM frames that are still profile, keeping per-peak IM at a tolerance no coarser
    // than the trace IM tolerance [C-2]. Split MS1 vs MS2 into separate working maps.
    if (!streamed)
    {
    PeakPickerIM picker = spicker;

    // Pick all IM frames in parallel. firstprivate gives each thread a deep copy (Param is a value type), so the
    // mutable ccs_warning_shown_ flag is not shared. [par-Crit-1]
    writeLogInfo_("Loaded " + String(exp.size()) + " frames." + clk_() + rss_()); // raw residency peak? [mem]
    writeLogInfo_("Peak-picking " + String(exp.size()) + " frames (OpenMP)...");
    // An exception escaping the OMP region would terminate the process: captured, rethrown serially.
    std::exception_ptr pick_err;
    // [dnoise] as in flushMS1_: the raw MS1 frames of a .d, before their pick (a frame it empties keeps its IM array,
    // so the dispatch below keeps it too)
    if (dnoise.on) dnoise.checkLoader();
    std::vector<DnoiseRun::Work> dn_work(dnoise.on ? (size_t)std::max(1, omp_get_max_threads()) : 0);
    std::vector<DnoiseRun::Counts> dn_cnt(dnoise.on ? exp.size() : 0);
    #pragma omp parallel for firstprivate(picker) schedule(dynamic, 16)
    for (long i = 0; i < (long)exp.size(); ++i)
    {
      try
      {
        MSSpectrum& s = exp[i];
        // [dnoise] before the ion-mobility skip, as the streaming consumer filters every MS1 frame: frame() refuses a
        // non-empty MS1 frame without 1/K0 and counts an empty one, so no MS1 frame bypasses it and both paths count alike
        if (dnoise.on && s.getMSLevel() == 1) dnoise.frame(s, dn_work.at((size_t)omp_get_thread_num()), dn_cnt[(size_t)i]);
        if (!s.containsIMData()) continue;
        if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) picker.pickIMCluster(s);
        ensureIMArrayName(s); // so MassTraceDetection reads per-peak IM (IM-aware tracing) [C-1]
        if (s.getMSLevel() == 1) requireFiniteMs1Mz(s);   // [ms1-finite] exactly the spectra the prune and the MS1 map take
      }
      catch (...)
      {
        #pragma omp critical
        if (!pick_err) pick_err = std::current_exception();
      }
    }
    if (pick_err) std::rethrow_exception(pick_err);
    dnoise.add(dn_cnt);
    // [dnoise] each thread's scratch keeps the capacity of its largest frame (up to ~4 GiB at 100 threads): free it with
    // the pick, not at the end of this block after the prune and the dispatch
    std::vector<DnoiseRun::Work>().swap(dn_work); std::vector<DnoiseRun::Counts>().swap(dn_cnt);
    // [ms1-prune] exactly the MS1 spectra the dispatch below keeps, pruned as flushMS1_ prunes a batch
    if (ms1_prune.on) ms1_prune.batch(exp.getSpectra(), [](const MSSpectrum& s) { return s.containsIMData() && s.getMSLevel() == 1; });

    // Serial dispatch: MOVE (not copy) frames into the MS1 map / per-window MS2 maps to avoid
    // transient memory doubling; container inserts aren't thread-safe so this stays serial.
    vector<uint32_t> cmzq; vector<float> cint; vector<uint16_t> cimq;   // one MS2 frame's compaction region, reused
    for (MSSpectrum& s : exp)
    {
      if (!s.containsIMData()) continue;
      if (s.getMSLevel() == 1)
      {
        ms1_map.addSpectrum(std::move(s));
      }
      else if (s.getMSLevel() == 2 && !s.getPrecursors().empty())
      {
        const OpenMS::Precursor& pr = s.getPrecursors()[0];
        double lo = pr.getMZ() - pr.getIsolationWindowLowerOffset();
        double hi = pr.getMZ() + pr.getIsolationWindowUpperOffset();
        const uint32_t fid = frameIdOf(s.getNativeID());
        if (fid == 0) ++cstat.unmapped_frame;
        if (cmzq.size() < s.size()) { cmzq.resize(s.size()); cint.resize(s.size()); cimq.resize(s.size()); }
        const size_t kept = compactifyInto(s, cstat, cmzq.data(), cint.data(), cimq.data(), cmzq.size());
        ms2_by_window[winKey(lo, hi, windowGroupOf(s.getNativeID()))].append(s.getRT(), fid, cmzq.data(), cint.data(), cimq.data(), kept, ms1_map.size());
        s.clear(true);   // [mem] the picked frame is dead once compacted; waiting for exp.clear()
      }                  // below held the whole picked run (20 B/peak) beside the compact store
    }
    }
    exp.clear(true); // frames moved out; release the container
    {
      size_t cb = 0, cp = 0;
      for (const auto& kv : ms2_by_window) { cb += kv.second.bytes(); cp += kv.second.peaks(); }
      const size_t lost = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
      if (lost)
      {
        writeLogWarn_("[compact] DROPPED " + String(lost) + " peaks: no_im_array=" + String(cstat.no_im_array)
                      + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz)
                      + " bad_im=" + String(cstat.bad_im) + " (kept " + String(cstat.kept) + ") -- if this is"
                      + " nonzero the compact store is NOT lossless and results differ from the old path.");
      }
      size_t ms1_pk = 0; for (const auto& sp : ms1_map) ms1_pk += sp.size();
      { long long sb = 0; for (const auto& kv : ms2_by_window) sb += (long long)kv.second.bytes();
        ledSet_(LC_SLAB, sb); ledSet_(LC_MS1MAP, (long long)ms1_pk * 20); }   // touched bytes; the 2.5x reserve is untouched pages
      writeLogInfo_("Split into MS1 + " + String(ms2_by_window.size()) + " windows (picked); compact store "
                    + String(cb / (1024ULL * 1024ULL)) + " MB for " + String(cp) + " peaks (~"
                    + String(cp ? (double)cb / cp : 0.0) + " B/peak); MS1 " + String(ms1_map.size()) + " frames / "
                    + String(ms1_pk) + " peaks (" + String(ms1_pk * 20 / (1024ULL * 1024ULL)) + " MB as PeakMap)." + clk_() + rss_() + mem_()); // [mem]
      if (dnoise.on)   // [dnoise] after the MS1 load, whichever path picked it: every MS1 frame with points exactly once
      { dnoise.requireComplete(); writeLogInfo_(dnoise.summary()); }
      if (ms1_prune.on)
      { const Ms1Prune::Counts& t = ms1_prune.total;
        writeLogInfo_("[ms1-prune] picked " + String(t.picked) + " survivors " + String(t.survivors) + " witnesses " + String(t.kept - t.survivors)
                      + " at_noise " + String(t.at_noise) + " nan " + String(t.nan) + " bad_mz " + String(t.bad_mz) + ": kept " + String(t.kept) + " peaks = "
                      + String(t.kept * 20 / (1024LL * 1024LL)) + " MB as PeakMap of " + String(t.picked * 20 / (1024LL * 1024LL)) + " MB picked; sample "
                      + String(ms1_prune.sample.size()) + " m/z from " + String(ms1_prune.frames) + " spectra"
                      + (ms1_prune.no_witness ? String(" -- DIASPEXTRACT_MS1_PRUNE_NO_WITNESS: survivors only, a negative control, output NOT identical") : String())
                      + (ms1_prune.no_chain ? String(" -- DIASPEXTRACT_MS1_PRUNE_NO_CHAIN: interior chain links dropped, a negative control, output NOT identical") : String())); }
    }

    if (ms1_map.empty())
    {
      writeLogError_("Error: no MS1 ion-mobility frames found; cannot seed precursors.");
      return INCOMPATIBLE_INPUT_DATA;
    }

    // [cell] the frozen frame tables (see FrozenWin), copied before indexRt() frees the slabs' frame
    // times and toTof() their frame ids
    map<WinKey, FrozenWin> frozen;
    if (tile_source)
      for (const auto& kv : run_tables)   // [3b] every window of the reader's table, no slab needed
      { FrozenWin& fw = frozen[kv.first]; fw.frame_rt = kv.second.rt; fw.frame_id = kv.second.id; }
    else
      for (const auto& kv : ms2_by_window)
      { FrozenWin& fw = frozen[kv.first]; fw.frame_rt = kv.second.rt; fw.frame_id = kv.second.frame_id; }
    // The global RT axis (see rtAxis()): one entry per distinct frame time, MS1 and MS2 alike.
    {
      vector<double>& ax = rtAxis();
      ax.clear();
      { size_t nf = ms1_map.size(); for (const auto& kv : ms2_by_window) nf += kv.second.frames();
        ax.reserve(nf + 1024); }   // MS1 *and* every window's MS2 frames land here [mem]
      for (const auto& sp : ms1_map) ax.push_back(sp.getRT());
      for (const auto& kv : frozen) ax.insert(ax.end(), kv.second.frame_rt.begin(), kv.second.frame_rt.end());   // the tables' frame times (== the slabs' after padding)
      sort(ax.begin(), ax.end());
      ax.erase(unique(ax.begin(), ax.end()), ax.end());
      writeLogInfo_("RT axis: " + String(ax.size()) + " distinct frame retention times"
                    + (ax.empty() ? "" : " spanning " + String(ax.front()) + "-" + String(ax.back()) + " s"));
      if (ax.empty())
      {
        writeLogError_("Error: no frame retention times; cannot build the RT axis.");
        return INCOMPATIBLE_INPUT_DATA;
      }
      // [tile-2d] frame times -> axis indices, and the load-time growth slack returned, window by
      // window (one window's copy in flight at a time, not all of them at once).
      const double _ts = phase_clock_(); const long _r0 = rss_mb_();
      size_t reserved = 0;
      for (auto& kv : ms2_by_window) { indexRt(kv.second); reserved += kv.second.reserved(); }
      for (auto& kv : frozen)
      { FrozenWin& fw = kv.second; fw.rt_index.resize(fw.frames()); for (size_t f = 0; f < fw.frames(); ++f) fw.rt_index[f] = rtIndex(fw.frame_rt[f]); }
#ifdef __GLIBC__
      if (trim_on) malloc_trim(0);   // the load's arena churn, returned before the MS1 phase
#endif
      writeLogInfo_("[slab] " + String(ms2_by_window.size()) + " windows indexed on the RT axis in " + String(phase_clock_() - _ts)
                    + " s; " + String(reserved / (1024ULL * 1024ULL)) + " MB reserved but untouched (virtual only) [RSS "
                    + String(_r0) + " -> " + String(rss_mb_()) + " MB]" + mem_());
    }

    //-------------------------------------------------------------
    // MS1 traces + precursor/charge inference (shared, once) [H-8,H-10]
    //-------------------------------------------------------------
    // [outer] sortSpectra() = the RT sort, then each spectrum's own m/z sort (a no-op check on a picked spectrum). The
    // per-spectrum sorts are independent, so they run in parallel after the same RT sort.
    {
      const double _tss = phase_clock_();
      ms1_map.sortSpectra(false);
      std::exception_ptr ss_err;
      #pragma omp parallel for schedule(dynamic, 16)
      for (long i = 0; i < (long)ms1_map.size(); ++i)
      { try { ms1_map[(Size)i].sortByPosition(); } catch (...) { taskErr_(ss_err); } }
      if (ss_err) std::rethrow_exception(ss_err);
      writeLogInfo_("[outer] MS1 spectra sorted in " + String(phase_clock_() - _tss) + " s");
    }
    const double ms1_snr = getDoubleOption_("trace:ms1_chrom_peak_snr");
    const double ms2_snr = getDoubleOption_("trace:ms2_chrom_peak_snr");
    const double ms2_minlen = getDoubleOption_("trace:ms2_min_length_sec");
    const double ms2_split = getDoubleOption_("trace:ms2_split_valleys");
    // [split-scan-time] ONE scalar for the whole run, from the MS1 frame times: one MS1 frame per
    // cycle, so it is the MS2 cadence too. Per window it would differ in the last bits (and per
    // tile later), and a sparse window would not be scanned differently from its neighbours.
    // [cell] The fixed RT cells (tile:rt_sec): interior edges T_1..T_{K-1} on MS1 frame times, so
    // every window's frames split at the same instants; cell k = [T_k, T_{k+1}), T_0 = -inf.
    vector<double> cell_edges;
    String cell_boundaries;
    double cell_pitch = getDoubleOption_("tile:rt_sec");
    vector<double> ms1_rt; ms1_rt.reserve(ms1_map.size());
    for (const auto& sp : ms1_map) ms1_rt.push_back(sp.getRT());
    std::sort(ms1_rt.begin(), ms1_rt.end());
    const double cycle = frameSpacing(ms1_rt);
    if (cell_pitch > 0.0)
    {
      const size_t n = cycle > 0.0 ? std::max<size_t>(1, (size_t)llround(cell_pitch / cycle)) : 0;
      if (n > 0)
      {
        for (size_t i = n; i < ms1_rt.size(); i += n) cell_edges.push_back(ms1_rt[i]);
        if (!cell_edges.empty() && 2 * (ms1_rt.size() - cell_edges.size() * n) < n) cell_edges.pop_back();   // remainder < pitch/2 joins the last cell
      }
      for (double e : cell_edges) { char buf[32]; std::snprintf(buf, sizeof buf, "%.17g", e); cell_boundaries += (cell_boundaries.empty() ? "" : ",") + String(buf); }   // round-trip exact
      writeLogInfo_("[cell] pitch " + String(cell_pitch) + " s = " + String(n) + " MS1 frames (cycle " + String(cycle) + " s): "
                    + String(cell_edges.size() + 1) + " cells, edges " + (cell_boundaries.empty() ? String("none") : cell_boundaries));
    }
    double split_scan_time = 0.0;
    if (getStringOption_("trace:split_scan_time") == "frame")
    {
      split_scan_time = cycle;
      if (split_scan_time > 0.0)
        writeLogInfo_("[split] scan time " + String(split_scan_time) + " s from " + String(ms1_rt.size()) + " MS1 frames -> chrom_fwhm / scan time is the splitter's window");
      else writeLogWarn_("[split] no positive gap between MS1 frame times; the valley splitter falls back to each trace's own average cycle time (trace:split_scan_time=trace)");
    }
    const double max_span = getDoubleOption_("trace:max_span_sec");
    int ms1_bands = getIntOption_("perf:ms1_trace_bands");
    if (ms1_bands < 1) ms1_bands = 1;
    double _t_ms1 = phase_clock_(), _c_ms1 = cpu_seconds_();   // [perf-instr]
    // detectTraces_ spawns its bands as tasks, so this phase opens its own parallel region.
    vector<Trace> ms1_traces;
    // The MS1 traces' frame table is the MS1 map's own spectrum sequence (RT-sorted above); their
    // spans live in this store, which must outlive scoring.
    TraceStore ms1_store;
    { vector<uint32_t> ri; ri.reserve(ms1_map.size());
      for (const auto& sp : ms1_map) ri.push_back(rtIndex(sp.getRT()));
      ms1_store.setFrames(ri, nullptr); }
    // [minlen] trace:min_length_sec < 0 (auto): the run's minimum MS1 trace length from its own peak widths. A fixed
    // sample of the MS1 map -- every peak in the 1-Th m/z slices [25k, 25k + 1), plus the band halo -- is traced like
    // the run (tolerances, snr, valley splitter, bands) but with no minimum length, so the estimate does not depend on
    // the threshold it sets; a whole Th per slice samples every mass-defect phase. The median FWHM W of the brightest
    // tenth of the traces centred in a slice (ranked by apex intensity, then FWHM: deterministic) sets the minimum to
    // trace:min_length_fwhm x W, at least half a cycle (2 MS1 points). Fewer than 100 measurable traces, or no cycle,
    // fall back to 3 s, the value before 1.6.0. A fixed value skips all of it and leaves the header unchanged.
    double ms1_min_len = getDoubleOption_("trace:min_length_sec");
    String minlen_stamp;   // spx:ms1_min_length, auto only
    if (ms1_min_len < 0.0)
    {
      const double _tm = phase_clock_();
      const double frac = getDoubleOption_("trace:min_length_fwhm"), ppm = getDoubleOption_("trace:mass_error_ppm");
      const double kPitch = 25.0, kWidth = 1.0;
      auto slicePos = [&](double x) { return x - kPitch * std::floor(x / kPitch); };   // position within its 25-Th period
      PeakMap smp;
      TraceStore smp_store;
      smp_store.setFrames(ms1_store.rt_index, nullptr);
      vector<Trace> st;
      size_t smp_pk = 0;
      const size_t og0 = rtOffGrid().load(), sb0 = rtStrictlyBetween().load();   // [soa-rt] counts the run's traces, not the sample's
      std::exception_ptr smp_err;
      #pragma omp parallel
      #pragma omp single
      try
      {
        // the cut: chunks of spectra in parallel, each spectrum kept in place (empty ones too: the frame sequence)
        vector<MSSpectrum> part(ms1_map.size());
        const int nchunk = std::max(1, std::min((int)ms1_map.size(), 4 * omp_get_max_threads()));
        std::exception_ptr cut_err;   // [task-err]
        #pragma omp taskloop grainsize(1) default(shared)
        for (int ci = 0; ci < nchunk; ++ci)
        { try {
          const size_t a = (size_t)ci * ms1_map.size() / nchunk, b = (size_t)(ci + 1) * ms1_map.size() / nchunk;
          for (size_t i = a; i < b; ++i)
          {
            const MSSpectrum& s = ms1_map[i];
            const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
            MSSpectrum& t = part[i]; t.setRT(s.getRT()); t.setMSLevel(1);
            OpenMS::DataArrays::FloatDataArray ia;
            ia.setName(Constants::UserParam::ION_MOBILITY);
            for (Size k = 0; k < s.size(); ++k)
            {
              const double x = s[k].getMZ(), h = x * ppm * 1e-6 * kBandHaloPpmMul, r = slicePos(x);
              if (r < kWidth + h || r >= kPitch - h) { t.push_back(s[k]); if (ima && k < ima->size()) ia.push_back((*ima)[k]); }
            }
            t.getFloatDataArrays().push_back(std::move(ia));
          }
        } catch (...) { taskErr_(cut_err); } }
        if (cut_err) std::rethrow_exception(cut_err);
        smp.reserveSpaceSpectra(part.size());
        for (auto& sp : part) { smp_pk += sp.size(); smp.addSpectrum(std::move(sp)); }
        vector<MSSpectrum>().swap(part);
        st = detectTraces_(smp, smp_store, delta_im, getDoubleOption_("trace:noise_threshold_int"), ms1_snr, 0.0,
                           getDoubleOption_("trace:ms1_split_valleys"), ms1_bands, split_scan_time, -1, LC_MS1BAND);
      }
      catch (...) { smp_err = std::current_exception(); }
      if (smp_err) std::rethrow_exception(smp_err);
      rtOffGrid().store(og0); rtStrictlyBetween().store(sb0);
      vector<pair<double, double>> iw;   // (apex intensity, FWHM s) of the traces centred in a slice
      for (const Trace& t : st)
        if (std::isfinite(t.mz) && slicePos(t.mz) < kWidth) iw.emplace_back(intensityOf(smp_store, t), xicFwhm(smp_store, t));
      std::sort(iw.begin(), iw.end(), [](const pair<double, double>& a, const pair<double, double>& b)
                { return a.first != b.first ? a.first > b.first : a.second < b.second; });
      const size_t top = iw.size() / 10;
      vector<double> w;
      for (size_t i = 0; i < top; ++i) if (iw[i].second > 0.0) w.push_back(iw[i].second);
      std::sort(w.begin(), w.end());
      auto q = [&w](double p) { return w[(size_t)(p * (double)(w.size() - 1) + 0.5)]; };
      const String head = "sample " + String(smp_pk) + " MS1 peaks in 1-Th slices every 25 Th -> " + String(st.size()) + " traces at min length 0, "
                          + String(iw.size()) + " centred in a slice; brightest tenth " + String(top) + ", " + String(w.size()) + " with a FWHM (>= 3 points)";
      if (w.size() >= 100 && cycle > 0.0)
      {
        const double w50 = q(0.5);
        ms1_min_len = std::max(frac * w50, 0.5 * cycle);
        const int pts = 1 + (int)std::min(std::ceil(ms1_min_len / cycle - 1e-9), 1e9);   // [review] capped before the int conversion
        writeLogInfo_("[minlen] auto: " + head + "; FWHM median " + String(w50) + " s = " + String(w50 / cycle) + " cycles (p10 " + String(q(0.1))
                      + " s = " + String(q(0.1) / cycle) + ", p90 " + String(q(0.9)) + " s = " + String(q(0.9) / cycle) + "), cycle " + String(cycle)
                      + " s -> trace:min_length_sec " + String(ms1_min_len) + " = max(" + String(frac) + " x median, half a cycle): >= "
                      + String(pts) + " MS1 points [" + String(phase_clock_() - _tm) + " s]");
        minlen_stamp = "auto " + String(ms1_min_len) + " s = max(" + String(frac) + " x FWHM median " + String(w50) + " s, half the cycle " + String(cycle)
                       + " s), >= " + String(pts) + " MS1 points; FWHM of " + String(w.size()) + " bright sampled traces";
      }
      else
      {
        ms1_min_len = 3.0;
        writeLogWarn_("[minlen] auto: " + head + ", cycle " + String(cycle) + " s: too few to estimate (100 needed, and a positive cycle); trace:min_length_sec 3, the value before 1.6.0");
        minlen_stamp = "auto fallback 3 s (" + String(w.size()) + " bright sampled traces with a FWHM, cycle " + String(cycle) + " s)";
      }
    }
    // [ms1-prune] every MS1 spectrum must have left its m/z sample before losing a peak, or the band edges
    // would come from a pruned map. Checked here, outside the parallel region an exception cannot leave.
    if (ms1_prune.on && (ms1_prune.frames != ms1_map.size() || ms1_prune.sample.size() != ms1_prune.sample_n))
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[ms1-prune] the pre-prune m/z sample covers " + String(ms1_prune.frames)
            + " spectra (" + String(ms1_prune.sample.size()) + " of " + String(ms1_prune.sample_n) + " m/z) but the MS1 map holds " + String(ms1_map.size())
            + "; -perf:ms1_prune false keeps every peak");
    // An exception must not leave the parallel region (terminate): captured inside the single, rethrown after it. This
    // covers detectTraces_'s own code, the [det] diagnostics included, and its taskloop tasks, which rethrow here ([task-err]).
    std::exception_ptr ms1_err;
    #pragma omp parallel
    #pragma omp single
    try
    {
      ms1_traces = detectTraces_(ms1_map, ms1_store, delta_im, getDoubleOption_("trace:noise_threshold_int"),
                                               ms1_snr, ms1_min_len, getDoubleOption_("trace:ms1_split_valleys"), ms1_bands,
                                               split_scan_time, LC_MS1MAP, LC_MS1BAND, &ms1_prune);   // [ledger] the partition takes the charge over
    }
    catch (...) { ms1_err = std::current_exception(); }
    if (ms1_err) std::rethrow_exception(ms1_err);
    { const double cap = max_span;
      // [outer] per trace, reading the store and writing only its own span fields: in parallel
      std::exception_ptr trim_err;
      if (cap > 0.0)
      {
        #pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)ms1_traces.size(); ++i)
        { try { trimToSpan(ms1_store, ms1_traces[(size_t)i], cap); } catch (...) { taskErr_(trim_err); } }
      }
      if (trim_err) std::rethrow_exception(trim_err); }
    vector<double>().swap(ms1_prune.sample);   // [ms1-prune] consumed by the band edges, or unused at one band
    ledSet_(LC_MS1TRACE, (long long)(ms1_traces.capacity() * sizeof(Trace)) + ms1_store.bytes());
    if (detOn_()) writeLogInfo_("[det] MS1 traces n=" + String(ms1_traces.size()) + " digest=" + String(traceDigest_(ms1_store, ms1_traces)));
    phaseAdd_("MS1_TRACE", _t_ms1, _c_ms1, false);

    double iso_im_tol = getDoubleOption_("charge:iso_im_tolerance");   // [isotol] 0 = auto, resolved after the content sort
    String iso_tol_stamp = String(iso_im_tol);
    const double ms2_noise = getDoubleOption_("trace:ms2_noise_threshold_int");
    writeLogInfo_("MS2 trace gate: apex > " + String(ms2_snr * ms2_noise) + " (snr " + String(ms2_snr)
                  + " x noise " + String(ms2_noise) + "), membership > " + String(ms2_noise)
                  + ", min_length " + String(ms2_minlen) + " s");
    const bool env_scoring = (getStringOption_("charge:scoring") == "envelope");
    // [determinism] inferPrecursors_ is greedy, so its output depends on ms1_traces' container order, which
    // EPD's parallel split makes thread-dependent: sort by content first. A non-finite coordinate would make
    // the comparator non-transitive (UB in sort), so those traces are dropped before it.
    {
      const size_t before = ms1_traces.size();
      ms1_traces.erase(remove_if(ms1_traces.begin(), ms1_traces.end(), [&ms1_store](const Trace& t) {
        return !std::isfinite(t.mz) || !std::isfinite(rtOf(ms1_store, t)) || !std::isfinite(t.im)
               || !std::isfinite(intensityOf(ms1_store, t));
      }), ms1_traces.end());
      if (before != ms1_traces.size())
        writeLogInfo_("dropped " + String(before - ms1_traces.size()) + " MS1 traces with non-finite coordinates");
    }
    // [tile-2g] Content-only key, never a container index or arena offset, so the order does not depend on
    // which tile detected a trace. Equal-key neighbours are counted: std::sort is not stable.
    auto traceLess = [&ms1_store](const Trace& a, const Trace& b) {
      if (a.mz != b.mz) return a.mz < b.mz;
      { const double ra = rtOf(ms1_store, a), rb = rtOf(ms1_store, b); if (ra != rb) return ra < rb; }
      if (a.im != b.im) return a.im < b.im;
      { const double ia = intensityOf(ms1_store, a), ib = intensityOf(ms1_store, b); if (ia != ib) return ia < ib; }
      if (a.npts != b.npts) return a.npts < b.npts;
      if (a.len != b.len) return a.len < b.len;
      return ms1_store.frame_rt[a.frame0] < ms1_store.frame_rt[b.frame0];
    };
    {
      // [outer] in parallel when that is provably std::sort's order (no equal-key neighbours), else std::sort itself
      const double _tcs = phase_clock_();
      const size_t ties = sortExact_(ms1_traces, traceLess);
      writeLogInfo_("[det] MS1 traces sorted by content: " + String(ms1_traces.size()) + " traces, "
                    + String(ties) + " equal-key neighbours");
      writeLogInfo_("[outer] MS1 content sort " + String(phase_clock_() - _tcs) + " s (" + (ties ? "equal keys: std::sort" : "unique order") + ")");
    }
    if (iso_im_tol <= 0.0) iso_im_tol = autoIsoImTol_(ms1_traces, ms1_store, max_charge, delta_rt, mass_ppm, iso_tol_stamp);   // [isotol]
    else writeLogInfo_("[isotol] 1/K0 tolerance for isotope partners " + String(iso_im_tol) + " (fixed)");
    writeLogInfo_(String("[isotol] partner rule: ") + (iso_closest_ ? "closest to the seed in apex RT and 1/K0 (charge:iso_partner closest)" : "first in m/z order"));
    std::optional<Phase> win_phase;                            // [perf-instr] the window loop (function scope)
    double _t_inf = phase_clock_(), _c_inf = cpu_seconds_();   // [perf-instr]
    vector<MonoV2In_> mono_v2in;   // [monocheck-v2] the count walks, parallel to `precursors` until monoPositionPass_
    vector<Precursor_> precursors = inferPrecursors_(ms1_traces, ms1_store, max_charge, delta_rt, iso_im_tol,
                                                     mass_ppm, env_scoring, mono_check_ == 2 ? &mono_v2in : nullptr);
    const double _t_inf_fn = phase_clock_();   // [outer] sub-step clock of PRECURSOR_INFER, log only
    if (adaptive_auto_)   // [adaptive] auto: decide the run's gate from its precursor span distribution (MS1 points)
    {
      size_t nf = 0;
      for (const Precursor_& pc : precursors) if ((int)ms1_traces[pc.trace_idx].np() <= faint_frames_) ++nf;
      const double frac = precursors.empty() ? 0.0 : (double)nf / (double)precursors.size();
      adaptive_ = frac >= auto_faint_frac_;
      if (!adaptive_) shrink_ = coloc_ = short_conc_ = false;   // auto:off = the plain gate for every precursor, not 'off'
      writeLogInfo_("[adaptive] auto: " + String(nf) + " of " + String(precursors.size()) + " precursors span <= " + String(faint_frames_)
                    + " MS1 points (" + String(frac) + ", threshold " + String(auto_faint_frac_) + "): faint gate " + (adaptive_ ? "ON" : "OFF"));
    }
    if (detOn_())
    {
      vector<uint64_t> b;
      for (const auto& pc : precursors) { pushBits_(b, pc.rt); pushBits_(b, pc.im); }
      writeLogInfo_("[det] precursors n=" + String(precursors.size()) + " digest=" + String(detDigest_(std::move(b))));
    }
    // [im-veto 2026-09-08] Charge halving, corrected by ion mobility. See charge:im_charge_veto.
    if (getStringOption_("charge:im_charge_veto") == "true")
    {
      auto fitLine = [&](int z, double& a, double& b, double& sig, size_t& n) -> bool {
        std::map<int, vector<double>> bins;                           // 25-Th bins of 1/K0, confident calls only
        for (const auto& pc : precursors)
          if (pc.charge == z && pc.n_isotopes >= 3 && std::isfinite(pc.im) && pc.mono_mz > 0.0)
            bins[(int)(pc.mono_mz / 25.0)].push_back(pc.im);
        vector<double> xs, ys; n = 0;
        for (auto& kv : bins)
        {
          auto& v = kv.second; if (v.size() < 30) continue;
          std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
          xs.push_back((kv.first + 0.5) * 25.0); ys.push_back(v[v.size() / 2]); n += v.size();
        }
        if (xs.size() < 4) return false;
        double sx = 0, sy = 0, sxx = 0, sxy = 0; const double m = (double)xs.size();
        for (size_t i = 0; i < xs.size(); ++i) { sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i]; }
        const double den = m * sxx - sx * sx; if (den <= 0.0) return false;
        b = (m * sxy - sx * sy) / den; a = (sy - b * sx) / m;
        vector<double> res;
        for (const auto& pc : precursors)
          if (pc.charge == z && pc.n_isotopes >= 3 && std::isfinite(pc.im) && pc.mono_mz > 0.0)
            res.push_back(fabs(pc.im - (a + b * pc.mono_mz)));
        if (res.empty()) return false;
        std::nth_element(res.begin(), res.begin() + res.size() / 2, res.end());
        sig = 1.4826 * res[res.size() / 2];
        return sig > 0.0;
      };
      double a2 = 0, b2 = 0, s2 = 0, a3 = 0, b3 = 0, s3 = 0; size_t n2 = 0, n3 = 0;
      const bool ok2 = fitLine(2, a2, b2, s2, n2), ok3 = fitLine(3, a3, b3, s3, n3);
      const double kband = getDoubleOption_("charge:im_veto_band");
      size_t n1 = 0, to2 = 0, to3 = 0;
      if (ok2)
      {
        for (auto& pc : precursors)
        {
          if (pc.charge != 1 || !std::isfinite(pc.im)) continue;
          ++n1;
          const double r2 = fabs(pc.im - (a2 + b2 * pc.mono_mz)) / s2;
          const double r3 = ok3 ? fabs(pc.im - (a3 + b3 * pc.mono_mz)) / s3 : 1e9;
          if (r2 <= kband && r2 <= r3) { pc.charge = 2; ++to2; }
          // On the 3+ band: unset the charge (re-calling it 3+ identified nothing; docs/BASELINE.md, "The
          // veto's 3+ arm dropped"). It becomes a guessed precursor, which require_isotope_support drops.
          else if (ok3 && r3 <= kband) { pc.charge = 0; ++to3; }
        }
        writeLogInfo_("[im-veto] 2+ line 1/K0 = " + String(a2) + " + " + String(b2) + " * m/z (n=" + String(n2) + ", sigma "
                      + String(s2) + ")" + (ok3 ? "; 3+ line " + String(a3) + " + " + String(b3) + " * m/z (n=" + String(n3) + ", sigma " + String(s3) + ")" : String("; no 3+ line"))
                      + ". Of " + String(n1) + " z=1 calls, " + String(to2) + " (" + String((int)(1000.0 * to2 / std::max<size_t>(n1, 1)) / 10.0)
                      + "%) sit on the 2+ band and are re-called 2+, " + String(to3) + " on the 3+ band are dropped, " + String(n1 - to2 - to3) + " kept as 1+");
      }
      else writeLogInfo_("[im-veto] too few confident 2+ calls to fit a mobility line; no z=1 re-calls");
    }
    // [monocheck-v2] At the final charge: after the veto (which reads the walks' own monoisotopes, as without v2),
    // before the charge-less precursors are defaulted or dropped (the veto edits in place, so the order still matches).
    const double _t_inf_v2a = phase_clock_();   // [outer]
    if (mono_check_ == 2)
    {
      monoPositionPass_(ms1_traces, ms1_store, precursors, mono_v2in, delta_rt, iso_im_tol, mass_ppm, env_scoring);
      vector<MonoV2In_>().swap(mono_v2in);
    }
    const double _t_inf_v2b = phase_clock_();   // [outer]
    const int default_charge = getIntOption_("assembly:default_charge");
    Size n_default = 0;
    // [way-4] A GUESSED charge and an isotope-SUPPORTED charge must not be emitted identically.
    // Closed search hides the difference (engine enumerates); open search does not.
    for (auto& pc : precursors)
      if (pc.charge == 0) { pc.guessed = true; pc.charge = default_charge; ++n_default; }
    // [emission] assembly:require_isotope_support: drop the guessed singletons.
    const bool require_iso = (getStringOption_("assembly:require_isotope_support") == "true");
    if (require_iso)
    {
      const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [](const Precursor_& pc){ return pc.guessed; }), precursors.end());
      writeLogInfo_("require_isotope_support: dropped " + String(before - precursors.size())
                    + " guessed precursors, " + String(precursors.size()) + " isotope-supported remain");
    }
    {
      const int zmin = getIntOption_("charge:min_charge"); const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [zmin](const Precursor_& pc){ return pc.charge > 0 && pc.charge < zmin; }), precursors.end());
      if (before != precursors.size())
        writeLogInfo_("charge:min_charge=" + String(zmin) + ": dropped " + String(before - precursors.size())
                      + " precursors below that charge, " + String(precursors.size()) + " remain");
    }
    // [step-02 emission-controlled arm, 2026-09-02] DIASPEXTRACT_MIN_ISOTOPES=k keeps only precursors whose
    // envelope has >= k isotope peaks (n_isotopes counts the mono): a precursor QUALITY gate that cuts
    // emission without touching fragment sharing. Pre-registered falsifier in BASELINE.md.
    if (const char* mi = std::getenv("DIASPEXTRACT_MIN_ISOTOPES"))
    {
      const int k = std::atoi(mi); const size_t before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(),
                       [k](const Precursor_& pc){ return pc.n_isotopes < k; }), precursors.end());
      writeLogInfo_("DIASPEXTRACT_MIN_ISOTOPES=" + String(k) + ": dropped " + String(before - precursors.size())
                    + " precursors with fewer isotope peaks, " + String(precursors.size()) + " remain");
    }
    // [ms1-funnel] diag:dump_ms1_tsv: MS1 traces and precursors, to attribute where a precursor is lost.
    const String ms1_dump = getStringOption_("diag:dump_ms1_tsv");
    if (!ms1_dump.empty())
    {
      auto fwhm_of = [&ms1_store](const Trace& tr) -> double { return xicFwhm(ms1_store, tr); };   // [fwhm]
      ofstream ft((ms1_dump + ".traces.tsv").c_str());
      // [model] 12 significant digits: the default 6 quantises m/z 765 to ~1.3 ppm.
      ft.precision(12);
      ft << "mz\trt\tim\tintensity\tfwhm_s\tn_xic\trt_first\trt_last\n";   // [tile-0d] span endpoints
      for (const auto& t : ms1_traces)
        ft << t.mz << '\t' << rtOf(ms1_store, t) << '\t' << t.im << '\t' << intensityOf(ms1_store, t) << '\t'
           << fwhm_of(t) << '\t' << t.np() << '\t'
           << (t.span() ? rtAtSpan(ms1_store, t, 0) : 0.0) << '\t' << (t.span() ? rtAtSpan(ms1_store, t, t.span() - 1) : 0.0) << '\n';
      ft.close();
      ofstream fp((ms1_dump + ".precursors.tsv").c_str());
      fp.precision(12);
      // [model] trace_idx makes the hypothesis -> seed-trace link explicit, so an offline replica
      // of inferPrecursors_ can be checked row-for-row rather than assumed to agree.
      fp << "mono_mz\tcharge\trt\tim\tguessed\ttrace_idx\n";
      for (const auto& pc : precursors)
        fp << pc.mono_mz << '\t' << pc.charge << '\t' << pc.rt << '\t' << pc.im << '\t'
           << (pc.guessed ? 1 : 0) << '\t' << pc.trace_idx << '\n';
      fp.close();
      writeLogInfo_("MS1 funnel dump: " + String(ms1_traces.size()) + " traces, "
                    + String(precursors.size()) + " precursors -> " + ms1_dump + ".{traces,precursors}.tsv");
    }
    // [twin-rt] assembly:twin_rt_frames: the twin fold's frame tolerance, |dRT| <= (n + 0.5) cycles (0 s = one frame).
    // Auto (0): n = the median FWHM of the final precursors' seed traces (>= 3 points), in cycles, rounded, >= 1. The
    // median of a multiset does not depend on the container order, so neither does n.
    double twin_rt_tol = 0.0;
    String twin_rt_desc;
    // [review] `ext` folds across frames only through v2's extension: without v2 no spectrum carries the mark, so the
    // relation is the same-frame one, and it is logged and stamped as such (the 1.5.0 header under the options reset)
    const bool twin_rt_inert = twin_rt_mode == 2 && mono_check_ != 2;
    if (twin_im >= 0.0 && twin_rt_frames >= 0 && cycle > 0.0 && twin_rt_inert)
      writeLogInfo_("[twins] assembly:twin_rt_frames: the same MS1 frame only (twin_rt_match ext needs charge:mono_profile_check v2)");
    if (twin_im >= 0.0 && twin_rt_frames >= 0 && cycle > 0.0 && !twin_rt_inert)
    {
      int n = twin_rt_frames;
      String how = String(n) + " (set)";
      if (n == 0)
      {
        vector<double> fw;
        fw.reserve(precursors.size());
        for (const Precursor_& pc : precursors) { const double f = xicFwhm(ms1_store, ms1_traces[pc.trace_idx]); if (f > 0.0) fw.push_back(f); }
        double med = -1.0;
        if (!fw.empty()) { const size_t m = fw.size() / 2; std::nth_element(fw.begin(), fw.begin() + (long)m, fw.end()); med = fw[m]; }
        // [review] rounded in double and capped before the int conversion (a near-zero cycle overflowed int)
        n = med > 0.0 ? (int)std::max(1.0, std::min(std::round(med / cycle), 1e6)) : 1;
        char b[200];
        std::snprintf(b, sizeof b, "%d (auto: median MS1 FWHM %.3f s over %zu seed traces = %.2f cycles)", n, med, fw.size(), med / cycle);
        how = String(b);
      }
      twin_rt_tol = (n + 0.5) * cycle;
      char b[96];
      std::snprintf(b, sizeof b, "|dRT| <= %.4f s", twin_rt_tol);
      twin_rt_desc = "MS1 frames up to " + how + " cycles apart (" + String(b) + ", cycle " + String(cycle) + " s"
                     + (twin_rt_mode == 0 ? ", precursor m/z within 20 ppm across frames" : twin_rt_mode == 1 ? ", precursor m/z within 0.05 ppm across frames"
                        : ", across frames precursor m/z within 0.05 ppm and one monoisotope placed by the v2 extension") + "), one RT cell";
      twin_ext_mark_ = twin_rt_mode == 2;
      writeLogInfo_("[twins] assembly:twin_rt_frames: the fold joins " + twin_rt_desc);
    }
    // sort precursors by mono m/z so each window can binary-search its members instead of
    // rescanning all ~1M precursors per window (total order for determinism).
    // [outer] trace_idx is unique (one precursor per committed seed), so the order is total and sortExact_ sorts in
    // parallel; an equal-key pair (a non-finite key) would send it to std::sort itself.
    const double _t_inf_srt = phase_clock_();
    const size_t prec_ties = sortExact_(precursors, [](const Precursor_& a, const Precursor_& b) {
      if (a.mono_mz != b.mono_mz) return a.mono_mz < b.mono_mz;
      if (a.rt != b.rt) return a.rt < b.rt;
      if (a.im != b.im) return a.im < b.im;
      return a.trace_idx < b.trace_idx;
    });
    const double _t_inf_cmp = phase_clock_();
    ledSet_(LC_PRECURSOR, (long long)(precursors.capacity() * sizeof(Precursor_) + precursors.size() * 8));
    vector<double> prec_mz(precursors.size());
    for (size_t i = 0; i < precursors.size(); ++i) prec_mz[i] = precursors[i].mono_mz;
    double _t_inf_trim = 0.0;
    // [mem] Scoring reads only the profiles of traces a precursor points at: release the rest.
    // (ms1_map is already empty: detectTraces_ consumed it.)
    {
      vector<bool> needed(ms1_traces.size(), false);
      for (const auto& pc : precursors)
        if (pc.trace_idx < ms1_traces.size()) needed[pc.trace_idx] = true;
      // Compacted in place, in offset order (see compactUnreferenced()).
      const Size freed = compactUnreferenced(ms1_traces, ms1_store, needed);
      ledSet_(LC_MS1TRACE, (long long)(ms1_traces.capacity() * sizeof(Trace)) + ms1_store.bytes());
      writeLogInfo_("Released XICs of " + String(freed) + " unreferenced MS1 traces." + rss_() + mem_());
      _t_inf_trim = phase_clock_();
#ifdef __GLIBC__
    if (trim_on) trimLogged_("MS1 -> window loop");
#endif
    }
    writeLogInfo_("[outer] PRECURSOR_INFER steps: inferPrecursors_ " + String(_t_inf_fn - _t_inf) + " s | veto " + String(_t_inf_v2a - _t_inf_fn)
                  + " s | monocheck v2 " + String(_t_inf_v2b - _t_inf_v2a) + " s | filters " + String(_t_inf_srt - _t_inf_v2b) + " s | precursor sort "
                  + String(_t_inf_cmp - _t_inf_srt) + " s (" + (prec_ties ? "equal keys: std::sort" : "unique order") + ") | compaction "
                  + String(_t_inf_trim - _t_inf_cmp) + " s | trim " + String(phase_clock_() - _t_inf_trim) + " s");
    writeLogInfo_(clk_() + " Detected " + String(ms1_traces.size()) + " MS1 traces -> " + String(precursors.size())
                  + " precursor hypotheses (" + String(n_default) + " assigned default charge)." + rss_());
    std::exception_ptr worker_err;

    //-------------------------------------------------------------
    // Window loop, tile by tile: the master admits windows, their work runs as tasks, results land in
    // index-addressed slots and each tile is canonically sorted, so output is thread-count-independent.
    // Exceptions are captured and rethrown serially.
    //-------------------------------------------------------------
    vector<pair<WinKey, PeakSlab*>> window_list;
    window_list.reserve(ms2_by_window.size());
    for (auto& kv : ms2_by_window) window_list.emplace_back(kv.first, &kv.second);
    const size_t n_win = frozen.size();   // the run's windows (the resident list is empty under the tiled source)
    vector<WinKey> fkeys; fkeys.reserve(n_win); for (const auto& kv : frozen) fkeys.push_back(kv.first);
    // [heavy-first] Biggest windows (compact bytes) first, so the largest is not left alone in the tail; output is sorted per tile.
    auto heavyFirst = [](const auto& a, const auto& b) { return a.second->bytes() > b.second->bytes(); };
    phaseAdd_("PRECURSOR_INFER", _t_inf, _c_inf, false);
    win_phase.emplace("WINDOW_LOOP");
    std::stable_sort(window_list.begin(), window_list.end(), heavyFirst);
    // MEMORY vs SPEED: every concurrently-processed window holds its own frames + traces + grid, so
    // peak RSS scales with the number of windows in flight, not with the window count. Capping
    // concurrency trades wall time for RAM. 0 = unlimited (previous behaviour). [mem]
    const int max_conc = getIntOption_("perf:max_concurrent_windows");
    const int n_threads = omp_get_max_threads();
    // [bands] perf:trace_bands splits each window's m/z range into bands traced concurrently (the window
    // count is fixed by the method); 0 = auto, enough bands that windows x bands covers every thread.
    const int bands_opt = getIntOption_("perf:trace_bands");
    int n_bands = bands_opt;
    if (n_bands <= 0)  // auto: enough bands that windows x bands covers every thread
      n_bands = (int)std::ceil((double)n_threads / std::max<size_t>(n_win, 1));
    n_bands = std::max(1, n_bands);
    // [master/worker] n_conc only caps windows IN FLIGHT (perf:max_concurrent_windows, 0 = uncapped); the
    // free-RAM admission gate is the real bound. Uncapped measured 7:04 vs 7:15 capped at the same 157 GB
    // peak on dataset D (context: docs/BASELINE.md, "WINDOW-LOOP PERFORMANCE LINE (2026-09-03, evening)").
    int n_conc = (int)n_win;
    if (max_conc > 0) n_conc = std::min(n_conc, max_conc);
    if (n_conc < 1) n_conc = 1;
    // outer loop over windows + inner loop over bands = two live levels
    if (n_bands > 1) omp_set_max_active_levels(2);
    writeLogInfo_("Parallelism: " + String(n_win) + " windows x " + String(n_bands)
                  + " m/z bands = " + String((int)n_win * n_bands) + " trace units; "
                  + String(n_threads) + " worker threads over one task pool; up to "
                  + String(n_conc) + " windows in flight (memory-bounded)");
    bool integer_detector = integer_req;
    String why = "not attempted";
    if (integer_detector && !tofAxis().ok)
    {
      // The flight-time axis comes from the input's own analysis.tdf: a .d directory has one, and an
      // mzPeak archive can be accompanied by one (DIASPEXTRACT_MZPEAK_TDF, as the exact-calibration
      // path already uses).
      why = "no candidate tdf";
      for (const std::string& c : tdf_cand)
        if (loadTofAxis(c, why)) { writeLogInfo_("flight-time axis from " + c + ": "
              + String(tofAxis().b_by_frame.size()) + " frame factors"); break; }
#ifdef DIASPEXTRACT_WITH_MZPEAK
      // An mzPeak archive carries its own vendor/analysis.tdf.gz, and the exact-m/z path has already
      // read MzCalibration + Frames.T1 from it: the same calibration the .d path loads.
      if (const auto mzp = spx::lastMzPeakMeta(); !tofAxis().ok && mzp && mzp->exact)
      {
        const spx::MzPeakExactMz& ex = *mzp->exact;
        setTofAxis(ex.cal, ex.t1_by_frame);
        tofAxis().acq_lo = ex.acq_lo; tofAxis().acq_hi = ex.acq_hi;
        tofAxis().n_bins = ex.n_bins; tofAxis().bounds_why = ex.bounds_why;
        writeLogInfo_("flight-time axis from the mzPeak archive's embedded analysis.tdf: "
                      + String(tofAxis().b_by_frame.size()) + " frame factors");
      }
#endif
    }
    // FAIL CLOSED: a frame without a parseable "frame=<Id>" cannot be calibrated per frame, and the
    // integer detector must not calibrate it at reference temperature -- fall back to openms instead.
    if (integer_detector && tofAxis().ok && cstat.unmapped_frame > 0)
    {
      writeLogWarn_("trace:detector=integer needs a per-frame calibration factor for EVERY frame, and "
                    + String(cstat.unmapped_frame) + " frame(s) carry a nativeID with no parseable "
                    "\"frame=<Id>\" key, so they cannot be addressed in the tdf's Frames table. Falling "
                    "back to trace:detector=openms rather than calibrating them at reference "
                    "temperature. A Bruker .d always names its frames; an mzML converted from one may "
                    "not. Pass -trace:detector openms to silence this.");
      integer_detector = false;
    }
    if (integer_detector && !tofAxis().ok)
    {
      // Fall back rather than refuse: the OpenMS detector needs no calibration, and spx:detector records which ran.
      writeLogWarn_("trace:detector=integer needs the vendor flight-time calibration (MzCalibration "
                    "+ Frames.T1 from analysis.tdf) and it was not available for this input ("
                    + why + "). Falling back to trace:detector=openms; the two detectors agree on about "
                    "85% of the union of identified peptides, so this file is NOT comparable peptide-for-peptide "
                    "with one traced on the integer axis. Pass -trace:detector openms to silence this.");
      integer_detector = false;
    }
    const double mem_frac = getDoubleOption_("perf:mem_fraction");
    std::atomic<size_t> inflight{0};
    std::atomic<int> adm_hwm{0};   // high-water mark of admitted windows, for the log and the e2e
    // [Fix A] One malloc_trim when the last window has materialised returns the load's bookkeeping left on
    // glibc's free lists (one tile only; docs/BASELINE.md, "2026-09-09: Fix A").
    std::atomic<size_t> n_materialised{0};
    // Projected footprint per window: ~44 B per compact peak sustained (measured over 28 windows of a 2-h
    // acquisition) + 256 MiB. docs/BASELINE.md, "Step 6: the admission gate now has bounds it can hold".
    auto project = [&](const PeakSlab& sl) { return sl.peaks() * (size_t)44 + (size_t)256 * 1024 * 1024; };
    for (auto& kv : frozen)
    {
      FrozenWin& fw = kv.second;
      fw.b.resize(fw.frames());
      for (size_t f = 0; f < fw.frames(); ++f) fw.b[f] = tofAxis().factor(fw.frame_id[f]);   // the value toTof() gives the slab
      fw.nearest_local = nearestLocal(fw.frame_rt, ms1_store.frame_rt, delta_rt);
    }
    { long long fb = 0;
      for (const auto& kv : frozen)
        fb += (long long)(kv.second.frame_rt.capacity() * 8 + kv.second.frame_id.capacity() * 4
                          + kv.second.rt_index.capacity() * 4 + kv.second.b.capacity() * 8 + kv.second.nearest_local.capacity() * 4);
      ledSet_(LC_FROZEN, fb); }
    if (cell_pitch > 0.0 && !integer_detector)
    {
      writeLogWarn_("tile:rt_sec applies to the integer detector; trace:detector=openms traces whole windows (no cells) -- the header records no grid.");
      cell_edges.clear(); cell_boundaries.clear(); cell_pitch = -1.0;
    }
    if (cell_pitch > 0.0 && bands_opt <= 0)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "tile:rt_sec needs a fixed perf:trace_bands: the auto setting (0) ties the band count to the thread count, and each band has its own `visited`, so the output would depend on the threads.");
    if (tile_source && integer_detector && bands_opt <= 0)
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "the .d streaming source needs a fixed perf:trace_bands: the auto setting derives the band count from the WINDOW COUNT, and the streaming source counts every window of the acquisition while the resident source counts only the windows that delivered a frame -- the two would band differently (review 3b).");
    // [tile] The cells grouped into TILES. A tile owns the precursors whose MS1 frame time lies in
    // [lo - tile_drt, hi - tile_drt) -- the half-open ranges partition RT (lo = -inf for the first,
    // hi = +inf for the last) -- and traces its windows on the tile's frames only; a precursor's
    // candidates (|rtOf - rt| <= tile_drt) then lie in [lo - 2 tile_drt, hi): the tile's own
    // fragments plus those CARRIED from the previous tile (rtOf >= lo - 2 tile_drt, re-carried
    // transitively). One tile = today's loop, statement for statement. tile_drt is the WIDEST RT
    // window any precursor's gate may use: gate:faint_delta_rt under gate:adaptive (a faint precursor
    // 3-6 s from a cut line otherwise sees different candidates per tile grid -- gate D 2026-09-26).
    const double tile_drt = adaptive_ ? std::max(delta_rt, faint_drt_) : delta_rt;
    struct Tile { double lo, hi; bool last; };
    vector<Tile> tiles;
    const int cells_per_tile = getIntOption_("tile:cells_per_tile");
    {
      const size_t n_cells_run = cell_edges.size() + 1;
      const bool one_piece = (out_type != FileTypes::MZML);   // mzPeak output is bulk-only: one tile
      const size_t per = (cells_per_tile <= 0 || cell_pitch <= 0.0 || one_piece) ? n_cells_run : (size_t)cells_per_tile;
      if (one_piece && cells_per_tile > 0 && cell_pitch > 0.0 && n_cells_run > 1)
        writeLogWarn_("[tile] mzPeak output is written in one piece: the " + String(n_cells_run) + " cells run as one tile (the whole run's memory); write mzML (-out_type mzML) for tiled memory.");
      for (size_t c0 = 0; c0 < n_cells_run; c0 += per)
      {
        const size_t c1 = std::min(n_cells_run, c0 + per);
        tiles.push_back(Tile{c0 == 0 ? -std::numeric_limits<double>::infinity() : cell_edges[c0 - 1],
                             c1 == n_cells_run ? std::numeric_limits<double>::infinity() : cell_edges[c1 - 1], c1 == n_cells_run});
      }
      String tl; for (const Tile& t : tiles) tl += (tl.empty() ? "" : " ") + String("[") + String(t.lo) + "," + String(t.hi) + ")";
      writeLogInfo_("[tile] " + String(tiles.size()) + " tile(s) of up to " + String(per) + " cell(s): " + tl);
    }
    // [tile] Run-level per-window band edges: the whole resident slab's flight-time extremes
    // (exactly the values the body computed on the whole window), so every tile traces in the
    // same bands. The bins are converted here, once, in place (toTof is idempotent).
    map<WinKey, pair<TofIdx, TofIdx>> wbands;
    if (integer_detector)
    {
      if (band_edges_acq && !(tofAxis().n_bins > 1))
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "trace:band_edges=acquisition needs DigitizerNumSamples and MzAcqRangeLower/Upper from the tdf's GlobalMetadata, which this input's calibration source did not provide (" + tofAxis().bounds_why + ").");
      if (!(tofAxis().n_bins > 1)) writeLogInfo_("[tof] no run-level axis bounds from the tdf (" + tofAxis().bounds_why + "); the acquisition edge set is not on record for this run");
      if (tile_source && !band_edges_acq)   // unreachable: refused before the load
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "trace:band_edges=slab needs the whole window resident; the .d streaming source reads it per tile.");
      for (const WinKey& k : fkeys) wbands[k] = {0xFFFFFFFFu, 0};
      #pragma omp parallel for schedule(dynamic, 1)
      for (long i = 0; i < (long)fkeys.size(); ++i)
      {
        const WinKey& k = fkeys[(size_t)i];
        TofIdx tlo = 0xFFFFFFFFu, thi = 0;
        const auto rs = ms2_by_window.find(k);   // a resident slab: converted here, once, in place
        if (rs != ms2_by_window.end())
        {
          PeakSlab& sl = rs->second;
          toTof(sl);
          for (TofIdx t : sl.tof) { tlo = std::min(tlo, t); thi = std::max(thi, t); }
        }
        // [3a] the run-level alternative: tlo 0, thi the larger of the digitizer's last bin
        // (DigitizerNumSamples - 1: no bin index can exceed it) and the acquisition's upper m/z
        // through the calibration at every frame of the frozen table; logged in both modes so
        // the two sets are on record side by side
        TofIdx thi_mz = 0, thi_acq = 0;
        if (tofAxis().n_bins > 1)
        {
          const TofIdx last = (TofIdx)(tofAxis().n_bins - 1);
          TofIdx thi_q = last;
          for (double b : frozen.at(k).b)
          {
            if (tofAxis().acq_hi > 0.0) thi_mz = std::max(thi_mz, tofAxis().tofOf(tofAxis().acq_hi, b));
            // the last bin's m/z, through the store's 1e-5 quantum, back to a bin: on a coarse
            // calibration the round trip can land one bin higher (review 3a, round 2)
            thi_q = std::max(thi_q, tofAxis().tofOf(tofAxis().mzOf(last, b) + 1.0e-5, b));
          }
          thi_acq = std::max(thi_q, thi_mz);
          if (thi_acq >= 4294967294u)
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[tof] the acquisition's upper bin is not representable (the partition would wrap)", String(thi_acq));
        }
        #pragma omp critical(winlog)
        writeLogInfo_("[tof] window " + String(k[0] / 100.0) + "-" + String(k[1] / 100.0) + " slab tlo=" + String(tlo) + " thi=" + String(thi)
                      + " | acquisition tlo=0 thi=" + String(thi_acq) + " (digitizer " + String(tofAxis().n_bins > 1 ? tofAxis().n_bins - 1 : 0L) + ", MzAcqRangeUpper " + String(thi_mz) + ")"
                      + (band_edges_acq ? " (in use)" : " (slab in use)"));
        wbands.at(k) = band_edges_acq ? pair<TofIdx, TofIdx>{0, thi_acq} : pair<TofIdx, TofIdx>{tlo, thi};
      }
    }
    // [tile] carried boundary fragments (records + spans), in for tile k = out of tile k - 1; and
    // the once-only ownership check per precursor (a window's m/z range is tile-independent)
    struct Carry { vector<Trace> rec; vector<float> inten;
                   long long bytes() const { return (long long)(rec.capacity() * sizeof(Trace) + inten.capacity() * 4); } };
    // [pipe] One slot per window: `c` is carry_out(w, k - 1), and `ready` counts the window's tiles that have
    // published theirs, so tile k may take `c` once ready == k. Under the per-tile barrier that always holds;
    // with pipelined tiles it is the one dependency between a window's tiles -- tile k + 1 may slab, prepare,
    // detect, merge and split while tile k still runs; it waits only where it takes the carry (assembly).
    struct CarrySlot { Carry c; std::atomic<size_t> ready{0}; };
    // [pipe] A window's handle on its slot for one tile. take() waits until the window's previous tile has
    // published, then moves that carry out; publish() hands this tile's carry (put() before it, or none) to
    // the next tile. The destructor publishes on every exit of the window body, waiting first if the body left
    // before take(), so a slot is never published out of order and the next tile never waits forever.
    struct CarryGate
    {
      CarrySlot& s; const size_t tk; bool took = false, done = false;
      void wait() const   // [pipe] idle until the carry arrives: help the tile reader or writer meanwhile
      { while (s.ready.load(std::memory_order_acquire) < tk) if (!teamAssist_(4)) std::this_thread::sleep_for(std::chrono::microseconds(200)); }
      Carry take() { wait(); took = true; Carry c = std::move(s.c); s.c = Carry(); ledAdd_(LC_CARRY, -c.bytes()); return c; }
      void put(Carry&& c) { ledAdd_(LC_CARRY, c.bytes()); s.c = std::move(c); }
      void publish()
      {
        if (done) return;
        if (!took) { (void)take(); }   // the body left early: drop the previous tile's carry, in order
        done = true;
        s.ready.store(tk + 1, std::memory_order_release);
      }
      ~CarryGate() { publish(); }
    };
    map<WinKey, CarrySlot> carry;
    map<WinKey, uint64_t> slab_chain;   // [det] the per-window slab digest, chained across tiles
    map<WinKey, vector<uint8_t>> visited;
    for (const WinKey& k : fkeys)
    {
      carry[k];
      const double lo = k[0] / 100.0, hi = k[1] / 100.0;
      const size_t plo = lower_bound(prec_mz.begin(), prec_mz.end(), lo) - prec_mz.begin();
      const size_t phi = upper_bound(prec_mz.begin(), prec_mz.end(), hi) - prec_mz.begin();
      visited[k].assign(phi > plo ? phi - plo : 0, 0);
    }
    // [tile] The output: settings and stamps are fixed before the first tile, the mzML writer is
    // opened once and fed a sorted block per tile (with one tile the count is exact and the bytes
    // are the bulk writer's); mzPeak output is bulk-only, so it runs as one tile (the tiling above).
    PeakMap out_exp;
    std::shared_ptr<DataProcessing> out_dp = std::make_shared<DataProcessing>(getProcessingInfo_(DataProcessing::DATA_PROCESSING));
    // [provenance] Header stamps that make an archived file attributable: calibration path (exported from
    // libOpenMS, so no header static is duplicated across the DSO boundary), detector, grid, intensity weights.
#ifdef DIASPEXTRACT_WITH_MZPEAK
    if (FileHandler::getTypeByFileName(in) == FileTypes::MZPEAK)
      out_exp.setMetaValue("spx:mz_calibration", spx::lastMzPeakCalibration());   // not BrukerTimsFile's state
    else
#endif
    out_exp.setMetaValue("spx:mz_calibration", BrukerTimsFile::lastMzCalibration());
    out_exp.setMetaValue("spx:detector", String(integer_detector ? "integer" : "openms"));
    out_exp.setMetaValue("spx:tile_rt_sec", cell_pitch);
    out_exp.setMetaValue("spx:tile_cells_per_tile", cells_per_tile);
    out_exp.setMetaValue("spx:band_edges", integer_detector ? getStringOption_("trace:band_edges") : String("none"));   // what was applied
    out_exp.setMetaValue("spx:tile_source", tile_source ? "stream" : "resident");
    out_exp.setMetaValue("spx:tiles", (int)tiles.size());
    out_exp.setMetaValue("spx:frame_table", frame_table_src);   // the frozen tables' source (loader | mzpeak | slab)
    out_exp.setMetaValue("spx:pearson_G", "frames");
    out_exp.setMetaValue("spx:tile_boundaries", cell_boundaries.empty() ? String("none") : cell_boundaries);
    out_exp.setMetaValue("spx:corr_power", corr_power_);
    out_exp.setMetaValue("spx:gate_adaptive", adaptive_auto_ ? String(adaptive_ ? "auto:on" : "auto:off") : String(adaptive_ ? "faint" : "off"));
    out_exp.setMetaValue("spx:gate_faint", String(faint_frames_) + " frames, " + String(faint_drt_) + " s, power " + String(faint_cp_));
    out_exp.setMetaValue("spx:gate_shrinkage", shrink_ ? String("fisher r0 ") + String(shrink_r0_) + " sd " + String(shrink_sd_) + " from " + String(shrink_min_pts_) : String("off"));
    out_exp.setMetaValue("spx:gate_colocation", coloc_ ? String("ranges ") + String(coloc_overlap_) : String("off"));
    out_exp.setMetaValue("spx:gate_short_traces", short_conc_ ? String("concentric tol ") + String(short_center_tol_) + " from " + String(short_min_pts_) : String("off"));
    out_exp.setMetaValue("spx:cap_rank", String(cap_rank_ == 2 ? "evidence" : cap_rank_ == 1 ? "umpire" : "corr_inten"));
    out_exp.setMetaValue("spx:same_bin_merge", String("im_slices ") + (!integer_detector ? String("n/a (openms detector)") : im_slice_merge_ >= 0.0 ? String(im_slice_merge_) + " (span " + String(2.0 * im_slice_merge_) + ")" : String("off"))
                         + ", assembly " + (merge_same_bin_ ? (integer_detector ? "on" : "on (exact m/z)") : "off"));
    if (anchor_rescue_)   // [anchor] stamped only when on: the default header is unchanged
      out_exp.setMetaValue("spx:anchor", String("anchors ") + String(anchor_min_count_) + " at r " + String(anchor_min_corr_) + ", rescue r "
                                          + String(anchor_rescue_corr_) + ", max " + String(anchor_max_rescued_) + " per spectrum");
    out_exp.setMetaValue("spx:im_weight_sigma", im_weight_sigma_);
    // [dnoise] what was applied: the port, "off", or skipped because the input is no .d
    out_exp.setMetaValue("spx:dnoise_ms1", dnoise.on ? String(String(spx::dnoise::kPortedFrom) + " MS1") : String(dnoise_req ? "skipped (input is not a Bruker .d)" : "off"));
    out_exp.setMetaValue("spx:dnoise_ms1_params", dnoise.on ? dnoise.params() : String("none"));
    out_exp.setMetaValue("spx:dnoise_ms1_points", dnoise.on ? String("raw " + String(dnoise.total.raw) + " kept " + String(dnoise.total.kept)) : String("none"));
    out_exp.setMetaValue("spx:require_isotope_support", (int)require_iso);
    if (!minlen_stamp.empty()) out_exp.setMetaValue("spx:ms1_min_length", minlen_stamp);   // [minlen] auto only: a fixed value leaves the header unchanged
    if (mono_check_) out_exp.setMetaValue("spx:monocheck", mono_check_stamp_);   // [monocheck] only when on
    if (getDoubleOption_("charge:iso_im_tolerance") <= 0.0 || iso_closest_)       // [isotol] only when not the 1.5.0 behaviour
      out_exp.setMetaValue("spx:iso_im_tolerance", iso_tol_stamp + (iso_closest_ ? ", closest partner in RT and 1/K0" : ", first partner in m/z"));
    if (post_fusion_cap_)   // [recap] only when on: the default header is unchanged
      out_exp.setMetaValue("spx:post_fusion_cap", String("max_fragments ") + String(max_frags) + " after the twin fold, ranked by the "
                           + (cap_rank_ == 2 ? "evidence" : cap_rank_ == 1 ? "umpire" : "corr_inten") + " cap key");
    out_exp.setMetaValue("spx:twin_merge", twin_im < 0.0 ? String("off") : twin_rt_desc.empty() ? String("same frame, same charge, 20 ppm, 1/K0 <= " + String(twin_im) + ", fused 10 ppm")
                         : String(twin_rt_desc + ", same charge, 20 ppm, 1/K0 <= " + String(twin_im) + ", fused 10 ppm"));   // [twin-rt] the 1.5.0 text when off
    if (frag_int_ != 0)   // [fragint] stamped only when on: the default header is unchanged
      out_exp.setMetaValue("spx:fragint", frag_int_ == 1 ? String("core ") + String(2 * frag_core_half_ + 1) + " frames" : String("sum"));
    if (foreign_w_ < 1.0)   // [foreign] stamped only when on (on by default since 1.6.0; 1 leaves the 1.5.0 header)
      out_exp.setMetaValue("spx:foreign", String("weight " + String(foreign_w_) + ", " + String(foreign_ppm_) + " ppm, 1/K0 " + String(foreign_dim_)
                                                 + ", min_iso " + String(foreign_min_iso_) + ", positions mono..M+3, own ion kept ("
                                                 + (foreign_ppm_ > 30.0 ? String(foreign_ppm_) : String("30")) + " ppm)"));   // the own-ion tests' max(ppm, 30)
    if (fuse_ppm_ > 0.0)   // [fuse] stamped only when on (on by default since 1.6.0; 0 leaves the 1.5.0 header)
    {
      String s = String(fuse_linkage_ ? "linkage " : "anchored ") + String(fuse_ppm_) + " ppm";
      if (fuse_linkage_) s += String(" (group span <= ") + String(2.0 * fuse_ppm_) + " ppm)";
      s += ", summed intensity, intensity-weighted m/z, cap key from the summed apex; after the same-bin merge, before min_fragments and max_fragments";
      out_exp.setMetaValue("spx:fragment_fuse", s);
    }
    std::unique_ptr<TileWriter> writer;
    // written as `<out>.part` and renamed on completion: a tile that throws must not leave a
    // structurally valid file holding the earlier tiles only (review, step 2)
    const String out_part = out + ".part";
    if (out_type == FileTypes::MZML) { writer.reset(new TileWriter(out_part)); writer->setExperimentalSettings(out_exp); }
    vector<MSSpectrum> all_out;   // the bulk (mzPeak) path only
    Size n_out = 0;
    Size n_twins = 0;   // [twins] spectra folded into a twin, over every tile
    std::array<size_t, 5> twin_st{};   // [twin-rt] mergeTwins_'s counters, over every tile
    std::array<size_t, 4> recap{};   // [recap] spectra re-capped, peaks removed, largest spectrum before, spectra seen
    writeLogInfo_("Processing " + String(n_win) + " isolation windows in " + String(tiles.size()) + " tile(s) from the " + (tile_source ? "streaming" : "resident") + " source (OpenMP over windows, <="
                  + String(n_conc) + " concurrent, admission bounded by " + String((int)(mem_frac * 100))
                  + "% of free RAM)..." + rss_() + mem_());

    // [master/worker] The master walks the windows heaviest-first and admits them against memory BEFORE
    // spawning, so the byte budget and the window cap hold by construction; every unit of work inside a
    // window is an OpenMP task, so idle threads steal across windows. Results go to index-addressed slots
    // (win_out[wi], pslot[j], per[ci][b]), never appended in completion order.
    std::atomic<int> active{0};
    // [pipe] Pipelined tiles. With the per-tile barrier, every tile waits for its slowest window chain and the
    // read, sort, merge, write and trim between two tiles overlap nothing. Pipelined, ONE parallel region
    // spans the whole loop and its single thread drives it: it reads tile k + 1 and admits its windows while
    // tile k's still run, and finishes (sorts, merges, writes, trims) the tiles strictly in tile order as each
    // completes. A window's tiles depend on each other only through the carry (CarryGate); tiles own disjoint
    // precursors and are sorted and written per tile as before, so the output is the barrier's, bit for bit.
    // DIASPEXTRACT_PIPE_TILES=<n> (default 2, 1 = the per-tile barrier) bounds the tiles whose windows may run
    // at once: tile j is read only after every window of tile j - n has finished, which bounds the slabs
    // resident and the chain of carry waits; the window cap and the memory admission below still apply to
    // every window in flight, over all tiles. The reader's and writer's loops run as tasks of the window team
    // (teamFor_); only the vendor reader's own decode loop (libOpenMS) opens a nested team.
    const int pipe_depth = []() { const char* e = std::getenv("DIASPEXTRACT_PIPE_TILES");
                                  const long v = e && *e ? std::atol(e) : 2L; return (int)std::max(1L, std::min(v, 16L)); }();
    const bool pipe = pipe_depth > 1 && tiles.size() > 1 && n_threads > 1;
    std::atomic<bool> pipe_abort{false};   // a window or the driver failed: admit nothing more, skip unstarted windows
    writeLogInfo_("[pipe] " + (pipe ? "tiles pipelined, up to " + String(pipe_depth) + " with windows in flight"
                                    : String("off (") + (pipe_depth <= 1 ? "DIASPEXTRACT_PIPE_TILES=1" : tiles.size() <= 1 ? "one tile" : "one thread")
                                      + "): a barrier after every tile"));
    // [pipe] One tile in flight: its slabs, its windows heaviest-first, their output slots, and how many of its
    // windows have not finished yet (the driver finishes the tile once that is 0).
    struct TileRun
    {
      size_t tk = 0;
      map<WinKey, PeakSlab> tile_slabs;
      vector<pair<WinKey, PeakSlab*>> wl;
      vector<vector<MSSpectrum>> win_out;
      Size win_done = 0;
      std::atomic<size_t> left{0};
    };
    auto loadTile = [&](TileRun& R)
    {
    const size_t tk = R.tk;
    const Tile& tile = tiles[tk];
    ledPhaseTile_("tile", tk + 1);
    // [tile] the tile's slabs: one tile = the resident slabs themselves (moved out by the body, as
    // today); several = a slice of each window's frames in [lo, hi) (the resident source keeps the
    // whole run: this arm is the exactness gate, not the memory configuration)
    map<WinKey, PeakSlab>& tile_slabs = R.tile_slabs;
    vector<pair<WinKey, PeakSlab*>>& wl = R.wl;
    if (tile_source)
    {
      // [3b] read the tile: a fresh consumer over the loader's RT range, lockstep-padded against
      // the frozen tables restricted to the tile, so a window's slab is exactly the table's frames
      // in [lo, hi); a window with none gets an empty slab (the body's empty-slab branch)
      Phase _ph("LOAD(tile)");
      const double _tr0 = phase_clock_();
      PeakMap ms1_none;
      PickCompactConsumer tc(spicker, tile_slabs, ms1_none, cstat);
      tc.tablefn = [&]() {
        PickCompactConsumer::FrameTables t;
        for (const auto& kv : frozen)
        {
          PickCompactConsumer::FrameTab& tab = t[kv.first];
          for (size_t f = 0; f < kv.second.frames(); ++f)
            if (kv.second.frame_rt[f] >= tile.lo && kv.second.frame_rt[f] < tile.hi) { tab.rt.push_back(kv.second.frame_rt[f]); tab.id.push_back(kv.second.frame_id[f]); }
        }
        return t; };
      BrukerTimsFile::Config c2 = cfg;
      c2.load_ms1 = false; c2.dia_ms2_rt_lo = tile.lo; c2.dia_ms2_rt_hi = tile.hi;
      BrukerTimsFile().loadDIAStreaming(in, tc, c2);   // a failure here is fatal: no fallback onto a half-read tile
      tc.finish();
      size_t n_del = 0, n_pad = 0, n_absent = 0, n_peaks = 0;
      for (const auto& kv : tc.tables())
      {
        PeakSlab& sl = tile_slabs[kv.first];   // a window that delivered nothing still owns its table's frames
        if (sl.frames() == 0 && !kv.second.rt.empty())
        {
          // every frame of the slice was skipped (an empty ion-mobility slice in each): pad them
          // all, exactly as the resident path does, so the frame sequence is the table's
          ++n_absent;
          for (size_t f = 0; f < kv.second.rt.size(); ++f) sl.append(kv.second.rt[f], kv.second.id[f], nullptr, nullptr, nullptr, 0, kv.second.rt.size());
          n_pad += kv.second.rt.size();
        }
        else { n_del += kv.second.delivered; n_pad += kv.second.padded; }
        if (sl.rt != kv.second.rt || sl.frame_id != kv.second.id)
          throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[3b] a tile's slab is not the frozen table's frames for the tile", String(kv.first[0] / 100.0) + "-" + String(kv.first[1] / 100.0));
        n_peaks += sl.peaks();
      }
      for (auto& kv : tile_slabs) indexRt(kv.second);
      { long long sb = 0; for (const auto& kv : tile_slabs) sb += (long long)kv.second.bytes();
        if (pipe) ledAdd_(LC_SLAB, sb - tc.slab_charged_); else ledSet_(LC_SLAB, sb); }   // [pipe] other tiles' slabs may still be live: replace only this read's batch charges
      if (detOn_() && integer_detector)   // the OpenMS path needs the slabs in m/z quanta (materializeWindow)
        // [det] the reader's exactness: the per-(window, tile) slab digest, chained over the tiles.
        // The tiles partition a window's frames in RT order and the digest is a sequential fold, so
        // the chain of a window's tiles EQUALS the whole-window digest a resident run prints -- the
        // check that the tiled reads are the resident peaks, frame for frame. The conversion to
        // flight-time bins has to happen first (it fills `b` and consumes `frame_id`); it is
        // idempotent, so the window body's own call becomes a no-op.
        for (auto& kv : tile_slabs)
        { toTof(kv.second);
          // the chain must START from the same FNV offset basis the whole-window digest starts
          // from, or it cannot equal it (a zero-initialised map entry would seed the fold with 0)
          uint64_t& h = slab_chain.try_emplace(kv.first, 1469598103934665603ULL).first->second;
          h = slabDigest_(kv.second, h);
          writeLogInfo_("[det] window " + String(kv.first[0] / 100.0) + "-" + String(kv.first[1] / 100.0) + " tile " + String(tk + 1)
                        + " slab n=" + String(kv.second.peaks()) + " digest=" + String(slabDigest_(kv.second)) + " chain=" + String(h)); }
      ledPhaseTile_("loop", tk + 1);
      writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + " [" + String(tile.lo) + "," + String(tile.hi) + ") read: "
                    + String(tile_slabs.size()) + " windows, frames delivered " + String(n_del) + " padded " + String(n_pad) + " absent-windows " + String(n_absent)
                    + ", " + String(n_peaks) + " peaks, " + String(phase_clock_() - _tr0, 1) + " s" + rss_() + mem_());
    }
    else if (tiles.size() > 1)
    {
      size_t tile_peaks = 0;
      for (const auto& w : window_list)
      {
        const PeakSlab& sl = *w.second;
        auto frt = [&](size_t f) { return rtAxis()[sl.rt_index[f]]; };
        size_t f0 = 0, f1 = sl.frames();
        { size_t lo = 0, hi = sl.frames(); while (lo < hi) { const size_t m = (lo + hi) / 2; if (frt(m) < tile.lo) lo = m + 1; else hi = m; } f0 = lo; }
        { size_t lo = f0, hi = sl.frames(); while (lo < hi) { const size_t m = (lo + hi) / 2; if (frt(m) < tile.hi) lo = m + 1; else hi = m; } f1 = lo; }
        PeakSlab& ts = tile_slabs[w.first];
        ts = sliceFrames(sl, f0, f1);
        tile_peaks += ts.peaks();
      }
      { long long sb = 0; for (const auto& kv : tile_slabs) sb += (long long)kv.second.bytes(); ledAdd_(LC_SLAB, sb); }
      writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + " [" + String(tile.lo) + "," + String(tile.hi) + "): "
                    + String(tile_slabs.size()) + " windows, " + String(tile_peaks) + " peaks resident" + rss_() + mem_());
    }
    if (tile_source || tiles.size() > 1)
    {
      for (auto& kv : tile_slabs) wl.emplace_back(kv.first, &kv.second);
      std::stable_sort(wl.begin(), wl.end(), heavyFirst);
    }
    else wl = window_list;
    R.win_out.resize(wl.size());
    R.left.store(wl.size());
    };   // loadTile
    // [pipe] the window body, defined below the admission that spawns it
    std::function<void(TileRun*, long, size_t)> runWindow;
    auto admitTile = [&](TileRun* R)
    {
    auto* rw = &runWindow;
    for (long wi = 0; wi < (long)R->wl.size() && !pipe_abort.load(); ++wi)
    {
      const size_t need = project(*R->wl[wi].second);
      for (;;)
      {
        const size_t budget = (size_t)(availableBytes_() * mem_frac);
        // [dyn-mem] admit while fewer than n_conc windows are in flight and the projection fits mem_fraction
        // of MemAvailable (re-read at every admission); an empty pipeline always admits, however large the window
        if ((active.load() == 0) || (active.load() < n_conc && inflight.load() + need <= budget)) break;
        // [review] A team of ONE (-threads 1, OMP_THREAD_LIMIT=1, OMP_DYNAMIC) has nobody else to run the windows admitted
        // so far, and libgomp's taskyield is a no-op: once the window cap or the memory budget held, this loop spun forever
        // (the barrier path since before the pipeline, and the pipeline's solo fallback). In a team of one every queued
        // window is a child of this thread, so taskwait runs them, which releases their admission charges.
        if (omp_get_num_threads() < 2)
        {
          #pragma omp taskwait
          continue;
        }
        if (!pipe)   // [pipe] the driver must not be handed a window here (a no-op in libgomp, not in every runtime)
        {
          #pragma omp taskyield
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
      }
      inflight += need; ++active;
      if (active.load() > adm_hwm.load()) adm_hwm.store(active.load());
      #pragma omp task default(none) firstprivate(rw, R, wi, need)
      (*rw)(R, wi, need);
    }
    };   // admitTile
    runWindow = [&](TileRun* R, const long wi, const size_t need)
    {
      const size_t tk = R->tk;
      const Tile& tile = tiles[tk];
      vector<pair<WinKey, PeakSlab*>>& wl = R->wl;
      vector<vector<MSSpectrum>>& win_out = R->win_out;
      Size& win_done = R->win_done;
      // A lambda, because `continue` is not a legal exit from a task's structured block and the
      // window body has one early exit (a window whose fragments all fail the finiteness filter).
      [&]() {
      try
      {
        const double win_lo = wl[wi].first[0] / 100.0;
        const double win_hi = wl[wi].first[1] / 100.0;
        // admitted by the master before this task existed; release both charges however we leave
        struct Release { std::atomic<size_t>& b; std::atomic<int>& a; size_t n;
                         ~Release() { b -= n; --a; } } rel{inflight, active, need};
        // [pipe] this window's carry slot: its destructor publishes this tile's carry on every exit
        CarryGate gate{carry.at(wl[wi].first), tk};
        if (pipe_abort.load()) return;   // [pipe] a window or the driver failed: the run is ending
        // Per-window stage seconds and the completion line, printed on EVERY exit of the body
        // (empty window, any assembly arm, exception): what the profile could only infer. [mem]
        double w_prep = 0, w_trace = 0, w_split = 0, w_score = 0, w_emit = 0;
        // [tsub] w_trace bundles more than the band taskloop -- the slab build, the serial arena
        // merge, and a full sort of every trace record in the window. Time them apart, so the next
        // decision about this stage is measured rather than attributed to the parallel part.
        double w_slab = 0, w_detect = 0, w_merge = 0, w_sort = 0;
        auto secs = [](auto d) { return std::chrono::duration<double>(d).count(); };
        struct WinDone { std::function<void()> f; ~WinDone() { f(); } } win_done_guard{[&]() {
          #pragma omp critical(winlog)
          writeLogInfo_("  window " + String(++win_done) + "/" + String(wl.size()) + " " + String(win_lo) + "-" + String(win_hi)
                        + ": prep " + String(w_prep) + " s | trace " + String(w_trace) + " s | split " + String(w_split)
                        + " s | score " + String(w_score) + " s | emit " + String(w_emit) + " s | spectra " + String(win_out[wi].size()) + clk_()
                        + "\n    [tsub] " + String(win_lo) + "-" + String(win_hi) + " trace = slab " + String(w_slab)
                        + " + detect " + String(w_detect) + " + merge " + String(w_merge) + " + sort " + String(w_sort) + " s");
        }};
        // [tile-2d] take ONLY this window's slab (a move: the store is the slab) and convert its
        // m/z quanta to flight-time bins in place. Only the OpenMS path needs the peaks converted
        // back into Peak1D doubles, which is the largest per-window allocation.
        PeakMap wmap;
        PeakSlab wslab;
        if (integer_detector)
        {
          { auto _tsl = std::chrono::steady_clock::now();
            wslab = std::move(*wl[wi].second);
            toTof(wslab);
            w_slab = secs(std::chrono::steady_clock::now() - _tsl); }
          if (detOn_())
          {
            #pragma omp critical
            writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " slab n=" + String(wslab.peaks()) + " digest=" + String(slabDigest_(wslab)));
          }
        }
        else
        {
          wmap = materializeWindow(*wl[wi].second);
          if (detOn_())   // [det] compact store digest BEFORE tracing: splits "pick" from "trace+assemble"
          {
            // materializeWindow consumed the frames; digest the materialised map instead (same content)
            vector<uint64_t> b;
            for (const auto& sp : wmap) for (const auto& pk : sp) { pushBits_(b, pk.getMZ()); float f = pk.getIntensity(); uint32_t w; std::memcpy(&w, &f, 4); b.push_back(w); }
            #pragma omp critical
            writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " compact n=" + String(b.size() / 2) + " digest=" + String(detDigest_(std::move(b))));
          }
        }
        if (++n_materialised == wl.size() && trim_on && tiles.size() == 1) trimLogged_("last window materialised -> store bookkeeping returned");
        wmap.sortSpectra();
        auto _t2 = std::chrono::steady_clock::now();
        vector<Trace> frag_traces;
        TraceStore wst;                                   // this window's frame table + arenas; outlives scoring
        // [ledger] the two structures that dominate a window: the store's arena and the trace
        // records. Scoped, so every exit path returns them to the counter.
        LedCharge _wst(LC_WSTORE, 0), _recs(LC_PARENT, 0);
        auto ledStore = [&]() { _wst.change(wst.bytes());
                                _recs.change((long long)(frag_traces.capacity() * sizeof(Trace))); };
        const FrozenWin& fwin = frozen.at(wl[wi].first);
        if (integer_detector && wslab.frames() == 0)
        {
          // [tile] no frame of this window lies in the tile: nothing to trace, the carried-in
          // fragments (below) are all it scores with
          wst.setFrames(fwin.rt_index, &fwin.b);
          #pragma omp critical(winlog)
          writeLogInfo_("[tile] window " + String(win_lo) + "-" + String(win_hi) + " tile " + String(tk + 1) + ": no frames in the tile");
        }
        else if (integer_detector)
        {
          wst.setFrames(fwin.rt_index, &fwin.b);             // the FROZEN table, not the resident slab's
          // Banded like the OpenMS path. A band SEEDS only in its own flight-time core but may
          // extend anywhere, so no halo is needed and no trace is cut; the one approximation is
          // that two bands can take the same peak (each has its own `visited`), which is the same
          // class of error the OpenMS halo carries -- measured at 0.1% of peptides. Bands are
          // TASKS, so they fill the pool alongside every other window's work.
          // [tile] the window's RUN-LEVEL band edges (the whole slab's extremes), the same in every tile
          const TofIdx tlo = wbands.at(wl[wi].first).first, thi = wbands.at(wl[wi].first).second;
          const int nb = (thi <= tlo) ? 1 : n_bands;
          size_t n_seeds = 0, n_parents = 0;
          const size_t n_peaks = wslab.peaks();
          {
          // The band edges are computed ONCE and shared by the seed segmentation and the detector,
          // so the two can never disagree about which band owns a seed.
          vector<TofIdx> bnd;
          if (nb > 1)
          {
            const double wdt = (double)(thi - tlo + 1) / nb;
            bnd.resize(nb + 1);
            for (int b = 0; b < nb; ++b) bnd[b] = tlo + (TofIdx)(b * wdt);
            bnd[nb] = thi + 1;
          }
          else bnd = {0, 0};
          // [cell] Detection per fixed RT cell: the window's frames are cut at the run-level cell
          // edges (tile:rt_sec), each cell is traced on a slice of its own -- seeds, `visited`,
          // frame_live and the tof index are the cell's -- and every cell's band stores are absorbed
          // into the window store afterwards, in (cell, band) order, with the frame indices rebased
          // by the cell's first frame. One cell (tile:rt_sec=-1) is the whole window on the slab
          // itself: the same statements, no copy, byte-identical to the uncelled code.
          // the cut is SLAB-local (sliceFrames indexes the slab); the rebase below is a VALUE
          // lookup of each cell's first frame time in the frozen table
          vector<double> slab_rt(wslab.frames());
          for (size_t f = 0; f < wslab.frames(); ++f) slab_rt[f] = rtAxis()[wslab.rt_index[f]];
          vector<size_t> cut{0};
          for (double e : cell_edges)
          { const size_t f = (size_t)(lower_bound(slab_rt.begin(), slab_rt.end(), e) - slab_rt.begin());
            if (f > cut.back() && f < slab_rt.size()) cut.push_back(f); }
          cut.push_back(slab_rt.size());
          const size_t n_cells = cut.size() - 1;
          const bool one_cell = (n_cells == 1);
          // Each band of each cell appends spans to a PRIVATE store; after the taskloops join, the
          // window store absorbs them in (cell, band) order and every offset is rebased (R4 sequencing).
          vector<vector<TraceStore>> bst(n_cells, vector<TraceStore>(nb));
          vector<vector<vector<Trace>>> per(n_cells, vector<vector<Trace>>(nb));
          long long band_bytes = 0, parent_bytes = 0;   // [ledger] returned to the counter at the merge
          for (size_t ci = 0; ci < n_cells; ++ci)
          {
            const size_t f0 = cut[ci], f1 = cut[ci + 1];
            PeakSlab cslab;
            if (!one_cell) cslab = sliceFrames(wslab, f0, f1);
            const PeakSlab& cs = one_cell ? wslab : cslab;
            if (!cs.tof.empty())
            {
              // the window's band edges must enclose the cell's bins: bandOf() and the tof index
              // wrap for a bin below tlo (review M2). A run-level freeze must be a superset. In
              // every mode, one cell included (review 3a: the acquisition edges are metadata).
              // [serial] the extremes as tasks: min and max are exact in any order
              const size_t np = cs.tof.size();
              const long nr = parTasks_(np, 1 << 18);
              vector<TofIdx> clos(nr, 0xFFFFFFFFu), chis(nr, 0);
              forTasks_(nr, [&](long c) {
                for (size_t k = (size_t)c * np / (size_t)nr; k < (size_t)(c + 1) * np / (size_t)nr; ++k)
                { clos[c] = std::min(clos[c], cs.tof[k]); chis[c] = std::max(chis[c], cs.tof[k]); } });
              const TofIdx clo = *std::min_element(clos.begin(), clos.end()), chi = *std::max_element(chis.begin(), chis.end());
              if (clo < tlo || chi > thi) throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "cell bins outside the window's band edges", String(clo) + "-" + String(chi) + " vs " + String(tlo) + "-" + String(thi));
            }
            auto _tp = std::chrono::steady_clock::now();
            TracePrep tprep = prepareTracing(cs, ms2_noise, ms2_snr, bnd, tlo, thi);   // once per cell
            LedCharge _prep(LC_PREP, (long long)(tprep.order.capacity() * 4 + tprep.frame_live.capacity()
                                                 + tprep.seg.capacity() * sizeof(size_t) + tprep.fidx.capacity() * 4));
            w_prep += secs(std::chrono::steady_clock::now() - _tp);
            n_seeds += tprep.order.size();
            auto _td0 = std::chrono::steady_clock::now();
            if (nb == 1)
              per[ci][0] = detectTracesInteger_(cs, tprep, bst[ci][0], mass_ppm, delta_im, ms2_noise, ms2_minlen);
            else
            {
              vector<vector<Trace>>& pc = per[ci]; vector<TraceStore>& bc = bst[ci];
              std::exception_ptr band_err;   // [task-err] rethrown into the window body's catch
              #pragma omp taskloop grainsize(1) default(shared)
              for (int b = 0; b < nb; ++b)
              { try {
                teamAssist_(1);   // [pipe] a chunk of the pipelined tile reader or writer, if one is live
                pc[b] = detectTracesInteger_(cs, tprep, bc[b], mass_ppm, delta_im, ms2_noise, ms2_minlen,
                                             bnd[b], bnd[b + 1], b);
              } catch (...) { taskErr_(band_err); } }
              if (band_err) std::rethrow_exception(band_err);
            }
            w_detect += secs(std::chrono::steady_clock::now() - _td0);
            { long long ab = 0, pb = 0;
              for (int b = 0; b < nb; ++b) { ab += (long long)(bst[ci][b].inten.capacity() * 4 + bst[ci][b].bins.capacity() * 2); pb += (long long)(per[ci][b].capacity() * sizeof(Trace)); }
              ledAdd_(LC_BAND, ab); ledAdd_(LC_PARENT, pb); band_bytes += ab; parent_bytes += pb; }
          }
          // Each cell's first frame in the FROZEN table, by value (the slab-local cut is an index
          // into the resident slab, which may hold only some of the table's frames later); the
          // resident frames must be, frame for frame, the table's -- asserted per cell.
          vector<size_t> f0z(n_cells);
          for (size_t ci = 0; ci < n_cells; ++ci)
          {
            const double t0 = slab_rt[cut[ci]];
            f0z[ci] = (size_t)(lower_bound(wst.frame_rt.begin(), wst.frame_rt.end(), t0) - wst.frame_rt.begin());
            for (size_t j = 0; j < cut[ci + 1] - cut[ci]; ++j)
              if (f0z[ci] + j >= wst.frames() || wst.rt_index[f0z[ci] + j] != wslab.rt_index[cut[ci] + j])
                throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                      "resident frames are not the frozen table's at cell " + String(ci) + " frame", String(j));
          }
          // The slab is dead after the frame check above: release it BEFORE the merge, the window's peak allocation.
          ledAdd_(LC_SLAB, -(long long)wslab.bytes());
          vector<TofIdx>().swap(wslab.tof); vector<float>().swap(wslab.inten);
          vector<uint16_t>().swap(wslab.imq); vector<uint32_t>().swap(wslab.frame_off);
          vector<double>().swap(wslab.b); vector<uint32_t>().swap(wslab.rt_index);
          auto _tm0 = std::chrono::steady_clock::now();
          size_t tot = 0, arena = 0; bool any_bins = false;
          for (size_t ci = 0; ci < n_cells; ++ci) for (int b = 0; b < nb; ++b)
          { tot += per[ci][b].size(); arena += bst[ci][b].inten.size(); any_bins |= !bst[ci][b].bins.empty(); }
          frag_traces.reserve(tot);
          wst.inten.reserve(arena); if (any_bins) wst.bins.reserve(arena);   // as in splitIntegerTraces
          if (nb > 1)
            writeLogInfo_("[mem] window " + String(win_lo) + "-" + String(win_hi) + " merge: " + String(n_cells) + " cell(s) x " + String(nb) + " band arenas, "
                          + String(arena) + " floats (" + String(arena * (any_bins ? 6 : 4) / (1024ULL * 1024ULL)) + " MB) x2 while absorbed, " + String(tot) + " traces");
          { // [serial] absorbed in (cell, band) order, rebased by the cell's first frame in the frozen table; copies as tasks
            vector<std::pair<TraceStore*, vector<Trace>*>> src; vector<uint32_t> fsh;
            for (size_t ci = 0; ci < n_cells; ++ci)
              for (int b = 0; b < nb; ++b) { src.push_back({&bst[ci][b], &per[ci][b]}); fsh.push_back((uint32_t)f0z[ci]); }
            absorbAll_(wst, frag_traces, src, fsh); }
          w_merge = secs(std::chrono::steady_clock::now() - _tm0);
          ledAdd_(LC_BAND, -band_bytes); ledAdd_(LC_PARENT, -parent_bytes); ledStore();
          }                                                  // the band stores die here
          // Valley splitting, as on the OpenMS path (skipping it lost 92% of peptides); it rebuilds the arena, and the per-frame bins are dead after it.
          n_parents = frag_traces.size();
          auto _ts = std::chrono::steady_clock::now();
          frag_traces = splitIntegerTraces(frag_traces, ms2_split, wst, n_bands * 4, split_scan_time);
          w_split = secs(std::chrono::steady_clock::now() - _ts); ledStore();
          if (max_span > 0.0)   // [serial] per trace, independent: as tasks over trace chunks
          { const size_t n = frag_traces.size(); const long nr = parTasks_(n, 1 << 16);
            forTasks_(nr, [&](long c) { for (size_t i = (size_t)c * n / (size_t)nr; i < (size_t)(c + 1) * n / (size_t)nr; ++i)
                                          trimToSpan(wst, frag_traces[i], max_span); }); }
          vector<int16_t>().swap(wst.bins); ledStore();
          writeLogInfo_("[mem] window " + String(win_lo) + "-" + String(win_hi) + ": seeds " + String(n_seeds) + " of " + String(n_peaks)
                        + " peaks (" + String((int)(1000.0 * n_seeds / std::max<size_t>(n_peaks, 1)) / 10.0) + "%); traces " + String(n_parents)
                        + " -> " + String(frag_traces.size()) + " (cap " + String(frag_traces.capacity()) + ", "
                        + String(frag_traces.capacity() * sizeof(Trace) / (1024ULL * 1024ULL)) + " MB); arena " + String(wst.inten.size()) + "/"
                        + String(wst.inten.capacity()) + " floats (" + String(wst.inten.capacity() * 4 / (1024ULL * 1024ULL)) + " MB)" + rss_());
        }
        else
        {
          { vector<uint32_t> ri; ri.reserve(wmap.size());          // wmap is RT-sorted above
            for (const auto& sp : wmap) ri.push_back(rtIndex(sp.getRT()));
            wst.setFrames(ri, nullptr); }
          frag_traces = detectTraces_(wmap, wst, delta_im, ms2_noise, ms2_snr, ms2_minlen, ms2_split, n_bands, split_scan_time);
          if (max_span > 0.0) for (auto& t : frag_traces) trimToSpan(wst, t, max_span);
        }
        if (detOn_())
        {
          #pragma omp critical
          writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " frag traces n=" + String(frag_traces.size()) + " digest=" + String(traceDigest_(wst, frag_traces)));
        }
        // (compact frames already released inside materializeWindow) [par-Crit-3]

        // [tile] the previous tile's boundary fragments join this window's list, spans copied
        // onto its arena; their frame indices are the frozen table's already (the same table in
        // every tile), and they enter BEFORE the content sort, so they take their canonical places.
        // [pipe] take() waits for this window's previous tile to publish them: the one point where a
        // window's tiles meet (everything above ran on this tile's slab alone).
        size_t n_carried_in = 0;
        { const Carry ci = gate.take();
          if (!ci.rec.empty())
          {
            frag_traces.reserve(frag_traces.size() + ci.rec.size());
            for (const Trace& r : ci.rec)
            {
              if (wst.inten.size() + r.len >= (size_t)std::numeric_limits<uint32_t>::max())
                throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);   // as absorb()
              Trace t = r; t.off = (uint32_t)wst.inten.size();
              wst.inten.insert(wst.inten.end(), ci.inten.begin() + r.off, ci.inten.begin() + r.off + r.len);
              frag_traces.push_back(t);
            }
            n_carried_in = ci.rec.size();
          } }
        // [tile] the precursors this tile OWNS, decided with the SCORER's own arithmetic: a
        // fragment at frame time x is a candidate of rt iff fabs(x - rt) <= tile_drt (doubles), so
        // rt has no candidate at or beyond an edge T iff (T - rt) > tile_drt -- the same
        // subtraction, monotone in T. The first tile k with (T_{k+1} - rt) > tile_drt owns rt:
        // exactly one does (-inf never, +inf always). Marked here, before any early exit, so the
        // coverage check after the last tile sees every precursor exactly once. Only this
        // window's precursors (by mono m/z); NaN sorts to the end, excluded; charge 0 kept.
        const size_t plo = lower_bound(prec_mz.begin(), prec_mz.end(), win_lo) - prec_mz.begin();
        const size_t phi = upper_bound(prec_mz.begin(), prec_mz.end(), win_hi) - prec_mz.begin();
        vector<uint32_t> pidx; pidx.reserve(phi > plo ? phi - plo : 0);
        for (size_t pi = plo; pi < phi; ++pi)
        {
          const Precursor_& pc = precursors[pi];
          if (!std::isfinite(pc.rt) || !std::isfinite(pc.im)) continue;
          if (!(tk == 0 || !(tile.lo - pc.rt > tile_drt))) continue;   // an earlier tile owns it
          if (!(tile.last || (tile.hi - pc.rt > tile_drt))) continue;  // a later tile owns it
          pidx.push_back((uint32_t)pi);
        }
        { vector<uint8_t>& vis = visited.at(wl[wi].first);
          for (uint32_t pi : pidx)
            if (vis[pi - plo]++ != 0)
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[tile] a precursor is owned by two tiles",
                                            String(precursors[pi].rt, 6) + " s in window " + String(win_lo) + "-" + String(win_hi)); }
        auto _tso = std::chrono::steady_clock::now();
        frag_traces.erase(remove_if(frag_traces.begin(), frag_traces.end(), [&wst](const Trace& t) {
          return !std::isfinite(t.im) || !std::isfinite(t.mz) || !std::isfinite(rtOf(wst, t));
        }), frag_traces.end());
        if (frag_traces.empty()) return;   // nothing traceable in this window

        // Content-only and TOTAL: traces equal in (im, mz, rt, apex intensity) fall through to the remaining
        // fields, so their order never depends on what else is resident; exact duplicates are interchangeable.
        // [serial] Sorted as tasks (parallelSort_): it may order exact duplicates differently from std::sort, which
        // the comparator already allows -- the order of what else is resident moves them too.
        parallelSort_(frag_traces.data(), frag_traces.size(), [&wst](const Trace& a, const Trace& b) {
          if (a.im != b.im) return a.im < b.im;
          if (a.mz != b.mz) return a.mz < b.mz;
          { const double ra = rtOf(wst, a), rb = rtOf(wst, b); if (ra != rb) return ra < rb; }
          { const double ia = intensityOf(wst, a), ib = intensityOf(wst, b); if (ia != ib) return ia < ib; }
          if (a.npts != b.npts) return a.npts < b.npts;
          if (a.len != b.len) return a.len < b.len;
          { const double fa = wst.frame_rt[a.frame0], fb = wst.frame_rt[b.frame0]; if (fa != fb) return fa < fb; }
          // ...then the profile and the bin fields (two bands can share an apex with different flanks)
          for (size_t k = 0; k < a.len; ++k) { const float xa = xv(wst, a, k), xb = xv(wst, b, k); if (xa != xb) return xa < xb; }
          if (a.tof != b.tof) return a.tof < b.tof;
          if (a.bframe != b.bframe) return a.bframe < b.bframe;
          if (a.rt_at != b.rt_at) return a.rt_at < b.rt_at;
          return false;
        }, parTasks_(frag_traces.size(), 1 << 15));
        // [tile] carry the boundary fragments to the next tile: every record (own or carried in)
        // whose RT is within 2 tile_drt of this tile's upper edge -- the candidates of a precursor
        // the next tile owns. Over ALL of frag_traces, so the carry is transitive by induction.
        // A precursor a later tile owns has (T - rt) <= tile_drt, its candidates x have
        // (rt - x) <= tile_drt: x >= T - 2 tile_drt up to two roundings of ~1e-12 s. The slack
        // over-carries by a hair (a carried fragment no precursor accepts changes nothing), never
        // under-carries (review, step 2: the exact form dropped a fragment at 1.2000000000000002).
        // [pipe] Taken and published HERE, after the content sort: scoring below only reads frag_traces
        // and the arena, so this is the carry the end of the body would produce, and the window's next
        // tile may assemble while this one scores.
        static const bool no_carry = std::getenv("DIASPEXTRACT_TILE_NO_CARRY") != nullptr;   // e2e falsifier only
        if (!tile.last && !no_carry)
        {
          Carry co;
          const double keep_from = tile.hi - 2.0 * tile_drt - 1e-6;
          for (const Trace& t : frag_traces)
            if (rtOf(wst, t) >= keep_from)
            {
              Trace c = t; c.off = (uint32_t)co.inten.size();
              co.inten.insert(co.inten.end(), wst.inten.begin() + t.off, wst.inten.begin() + t.off + t.len);
              co.rec.push_back(c);
            }
          if (n_carried_in || !co.rec.empty())
          {
            #pragma omp critical(winlog)
            writeLogInfo_("[tile] window " + String(win_lo) + "-" + String(win_hi) + " tile " + String(tk + 1) + ": carried in "
                          + String(n_carried_in) + ", carrying out " + String(co.rec.size()) + " fragments (" + String((co.rec.size() * sizeof(Trace) + co.inten.size() * 4) / 1024) + " KB)");
          }
          gate.put(std::move(co));
        }
        gate.publish();
        vector<double> frag_im(frag_traces.size());
        { const size_t n = frag_traces.size(); const long nr = parTasks_(n, 1 << 17);   // [serial] as tasks
          forTasks_(nr, [&](long c) { for (size_t i = (size_t)c * n / (size_t)nr; i < (size_t)(c + 1) * n / (size_t)nr; ++i) frag_im[i] = frag_traces[i].im; }); }
        FragRt frag_rt; frag_rt.build(frag_traces, wst, tile_drt, parTasks_(frag_traces.size(), 1 << 17));   // [adaptive] buckets cover the wider window
        w_sort = secs(std::chrono::steady_clock::now() - _tso);

        // Precompute the common correlation grid once per window; pdense is reused across the
        // window's precursors (this window runs on one thread). [perf]
        auto _t3 = std::chrono::steady_clock::now();
        w_trace = secs(_t3 - _t2) - w_prep - w_split;
        FragStats fg = buildFragStats(frag_traces, wst, fwin);
        ledStore();   // the finiteness filter and the content sort have run
        LedCharge _frag(LC_FRAG, (long long)(frag_im.capacity() * 8 + frag_rt.rt.capacity() * 8 + frag_rt.q.capacity() * 2
                                             + fg.mean.capacity() * 8 + fg.invnorm.capacity() * 8));
        vector<float> pdense(wst.frames(), 0.0f);   // indexed by window-local frame
        // [foreign] the window's precursors by apex RT, owned by this tile or not: the same candidates in any tiling
        ForeignIdx_ fx;
        if (foreign_w_ < 1.0) fx.build(precursors, plo, phi);

        vector<MSSpectrum>& bucket = win_out[wi];
        {
          // Share-all (default, current best): each fragment goes to EVERY gated precursor.
          //
          // [perf-3] Precursors are the unit of parallelism (one window can hold ~24% of them). Each result goes to
          // its own pre-sized slot (an earlier append to the shared bucket lost spectra), so task granularity cannot
          // drop or reorder an emission; pdense is per-task scratch and must stay firstprivate.
          auto _t5 = std::chrono::steady_clock::now();
          vector<MSSpectrum> pslot(pidx.size());   // [tile] the owned precursors, in index order as before
          // [ledger] ~700 B of empty header per OWNED precursor, most of which never emit, plus the
          // payloads as they are produced -- charged HERE and not at the end of the stage, so the
          // scoring stage's real footprint is visible instead of appearing all at once after it.
          LedCharge _pslot(LC_SPECTRA, (long long)pidx.size() * (long long)sizeof(MSSpectrum));
          // [gomp-cliff] Bound the TASK COUNT: libgomp's GOMP_taskloop runs the whole loop inline in the
          // creating thread once team task_count + num_tasks > 64 x nthreads, silently serialising it
          // (grainsize does not bound the count). So <= n_threads worker tasks pull 32-precursor chunks from
          // a shared cursor; keep the team-wide total well under the threshold when adding task layers.
          // docs/BASELINE.md, "Scheduling: the window loop was bounded by one window's serial chain".
          const long n_prec_w = (long)pidx.size();
          const long n_task_w = std::max(1L, std::min((n_prec_w + 31) / 32, (long)n_threads));
          std::atomic<long> pcur{0};
          std::exception_ptr asm_err;   // [task-err] rethrown into the window body's catch
          #pragma omp taskloop num_tasks(n_task_w) firstprivate(pdense) default(shared)
          for (long w = 0; w < n_task_w; ++w)
          { try {
            taskFault_("assemble");
            for (;;)
            {
              teamAssist_(1);   // [pipe] a chunk of the pipelined tile reader or writer, if one is live
              const long cs = pcur.fetch_add(32);
              if (cs >= n_prec_w) break;
              const long ce = std::min(cs + 32, n_prec_w);
              for (long j = cs; j < ce; ++j)
              {
                const Precursor_& pc = precursors[pidx[(size_t)j]];
                MSSpectrum ms2;
                assembleOne_(pc, win_lo, win_hi, frag_traces, wst, frag_im, frag_rt, ms1_traces, ms1_store, fg,
                             gate_dim, delta_rt, min_corr, min_corr_pts, min_frags, max_frags, pdense, ms2,
                             foreign_w_ < 1.0 ? &fx : nullptr);
                if (!ms2.empty())
                { ledAdd_(LC_SPECTRA, (long long)(sizeof(MSSpectrum) + ms2.size() * (post_fusion_cap_ ? 16 : 12))); pslot[(size_t)j] = std::move(ms2); }   // [recap] +4 B/peak key
              }
            }
          } catch (...) { taskErr_(asm_err); } }
          if (asm_err) std::rethrow_exception(asm_err);
          auto _t6 = std::chrono::steady_clock::now();
          w_score = secs(_t6 - _t5);
          for (auto& s : pslot) if (!s.empty()) bucket.push_back(std::move(s));
          // (the payloads were charged at emission; moving them into the bucket changes no bytes)
          const auto _t7 = std::chrono::steady_clock::now();
          w_emit = secs(_t7 - _t6);
        }
      }
      catch (...)
      {
        #pragma omp critical
        if (!worker_err) worker_err = std::current_exception();
        pipe_abort.store(true);
      }
      }();   // window body
      R->left.fetch_sub(1, std::memory_order_acq_rel);   // [pipe] last: the driver may finish the tile after this
    };   // runWindow
    // [pipe] After every window of the tile has finished: the coverage check (last tile), the canonical sort,
    // the twin fold, the post-fusion cap and the write of the tile's block, then the between-tile trim.
    auto finishTile = [&](TileRun& R)
    {
    const size_t tk = R.tk;
    const Tile& tile = tiles[tk];
    vector<vector<MSSpectrum>>& win_out = R.win_out;
    if (tile.last)
    {
      // [tile] coverage: every finite precursor of every window was owned by exactly one tile
      // (the per-tile check catches a second owner; this catches none)
      size_t n_cov = 0;
      for (const WinKey& k : fkeys)
      {
        const double lo = k[0] / 100.0, hi = k[1] / 100.0;
        const size_t plo = lower_bound(prec_mz.begin(), prec_mz.end(), lo) - prec_mz.begin();
        const size_t phi = upper_bound(prec_mz.begin(), prec_mz.end(), hi) - prec_mz.begin();
        const vector<uint8_t>& vis = visited.at(k);
        for (size_t pi = plo; pi < phi; ++pi)
        {
          if (!std::isfinite(precursors[pi].rt) || !std::isfinite(precursors[pi].im)) continue;
          if (vis[pi - plo] != 1)
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[tile] a precursor was owned by no tile",
                                          String(precursors[pi].rt, 6) + " s in window " + String(lo) + "-" + String(hi));
          ++n_cov;
        }
      }
      writeLogInfo_("[tile] coverage: " + String(n_cov) + " precursor-window pairs owned exactly once over " + String(tiles.size()) + " tile(s)");
      if (tile_source)
      {
        // the tile reads share `cstat` with pass 1; the pre-pass-2 report could not see them
        const size_t lost2 = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
        if (lost2) writeLogWarn_("[compact] DROPPED " + String(lost2) + " peaks over the whole run (pass 1 + every tile): no_im_array="
                                 + String(cstat.no_im_array) + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz)
                                 + " bad_im=" + String(cstat.bad_im) + " (kept " + String(cstat.kept) + ")");
        else writeLogInfo_("[compact] no peak dropped over the whole run (kept " + String(cstat.kept) + ")");
        if (detOn_()) for (const auto& kv : slab_chain)
          writeLogInfo_("[det] window " + String(kv.first[0] / 100.0) + "-" + String(kv.first[1] / 100.0) + " slab chain=" + String(kv.second));
      }
      // MEMORY: every MS1 structure is dead once the last tile's window loop has joined -- scoring
      // was its only reader (assembleOne_ dereferences ms1_traces[pc.trace_idx]
      // and its span). Released BEFORE this tile's output is assembled and written. [mem]
      vector<Trace>().swap(ms1_traces);
      ms1_store = TraceStore();
      vector<Precursor_>().swap(precursors);
      vector<double>().swap(prec_mz);
      ledSet_(LC_MS1TRACE, 0); ledSet_(LC_PRECURSOR, 0);
      writeLogInfo_("[mem] MS1 traces, spans, precursors released after the window loop." + rss_() + mem_());
    }
    {
      // [tile] this tile's spectra in the canonical order: written as one block (mzML), or kept for
      // the bulk writer (mzPeak, one tile). Tiles are RT-disjoint, so the blocks concatenate to the
      // whole run's order.
      vector<MSSpectrum> tile_out;
      { size_t n_tot = 0; for (const auto& bucket : win_out) n_tot += bucket.size(); tile_out.reserve(n_tot); }
      for (auto& bucket : win_out) { for (auto& s : bucket) tile_out.push_back(std::move(s)); bucket.clear(); bucket.shrink_to_fit(); }
      long long tile_spec_bytes = 0;
      for (const auto& sp : tile_out) tile_spec_bytes += (long long)(sizeof(MSSpectrum) + sp.size() * (post_fusion_cap_ ? 16 : 12));   // as charged at emission
      { Phase _ph("SORT(canonical)"); canonicalSort_(tile_out); }
      if (twin_im >= 0.0)
      {
        Phase _ph("MERGE(twins)");
        const size_t before = tile_out.size(), gone = mergeTwins_(tile_out, twin_im, post_fusion_cap_, twin_rt_tol, &cell_edges, tile_drt, &twin_st, twin_rt_mode);
        n_twins += gone;
        writeLogInfo_("[twins] tile " + String(tk + 1) + "/" + String(tiles.size()) + ": " + String(gone) + " of " + String(before) + " spectra folded into a twin");
      }
      if (post_fusion_cap_)
      {
        Phase _ph("RECAP(post-fusion)");
        const std::array<size_t, 3> r = capAfterFusion_(tile_out, max_frags);
        recap[0] += r[0]; recap[1] += r[1]; recap[2] = std::max(recap[2], r[2]); recap[3] += tile_out.size();
      }
      if (writer)
      {
        for (auto& sp : tile_out) sp.getDataProcessing().push_back(out_dp);   // ONE entry, by pointer, for every spectrum
        if (tiles.size() == 1) writer->setExpectedSize(tile_out.size(), 0);   // exact count: the bulk writer's bytes
        ledPhaseTile_("write", tk + 1);
        { Phase _ph("WRITE(mzML)"); writer->writeBlock(tile_out); }
        n_out += tile_out.size();
        ledAdd_(LC_SPECTRA, -tile_spec_bytes); ledSet_(LC_WRITER, (long long)n_out * 40);
        writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + ": " + String(tile_out.size()) + " spectra written" + clk_() + rss_() + mem_());
      }
      else
      {
        all_out.reserve(all_out.size() + tile_out.size());
        for (auto& sp : tile_out) all_out.push_back(std::move(sp));   // [ledger] the same bytes, still charged
      }
    }
#ifdef __GLIBC__
    // [trim] Between tiles, under perf:malloc_trim: nothing else returns the allocator's free pages here -- the in-loop
    // trim above runs only when there is one tile -- so a tiled loop's RSS would ratchet toward the arena pass 1 sized.
    // After tile k is written, before tile k+1 is read; about 4 GB off TNBC 009 and 0.5 GB off dataset D. Output-identical.
    // [pipe] Pipelined, the later tiles in flight keep their pages: the trim returns what tile k left behind.
    if (trim_on && tk + 1 < tiles.size())
    { const double _tt = phase_clock_(); const long _r0 = rss_mb_();
      malloc_trim(0);
      writeLogInfo_("[trim] tile " + String(tk + 1) + "/" + String(tiles.size()) + ": " + String(_r0) + " -> " + String(rss_mb_()) + " MB RSS in "
                    + String(phase_clock_() - _tt) + " s" + mem_()); }
#endif
    };   // finishTile

    if (!pipe)
      for (size_t tk = 0; tk < tiles.size(); ++tk)
      {
        TileRun R;
        R.tk = tk;
        loadTile(R);
        #pragma omp parallel num_threads(n_threads)
        #pragma omp single
        admitTile(&R);
        if (worker_err) std::rethrow_exception(worker_err);
        finishTile(R);
      }
    else
    {
      // [pipe] The driver is thread 0 of ONE parallel region over all tiles: it reads tile fin + pipe_depth - 1
      // while tile fin's windows run, admits its windows as the cap and memory allow, and once tile fin's windows
      // are done hands the tile to the WRITER, a task that finishes the tiles strictly in tile order. So the read
      // of the next tile never waits behind a write (a first build ran both on the driver: on a 60-min run at 100 threads
      // read + write per heavy tile outlasted the tile's windows, which then sat idle at three boundaries).
      // The driver does not run windows: it reaches no taskgroup end or taskwait (teamFor_, and parTasks_ keeps its
      // passes inline), so its n_threads - 1 team mates run the windows, the writer, and the reader's and writer's helper
      // tasks when they are free. One exception remains: past libgomp's task cliff (team task count > 64 x nthreads) a
      // spawned window or writer runs undeferred in the driver -- a stall, not a deadlock (review 2026-10-01).
      // Thread 0 rather than `single`'s first arrival: the vendor reader's decode loop opens a nested team next
      // to the driver's place, so under OMP_PROC_BIND=close it shares the window team's cores, never others.
      // Deadlock-free: a window waits only on its own previous tile (CarryGate), which was admitted earlier and
      // is dequeued earlier (the team queue is FIFO), and a tile is read only after tile j - pipe_depth finished,
      // so that previous tile never waits itself; the writer of tile f waits only for the writer of tile f - 1,
      // spawned (and dequeued) before it.
      omp_set_max_active_levels(std::max(2, omp_get_max_active_levels()));
      vector<std::unique_ptr<TileRun>> runs(tiles.size());
      std::exception_ptr drv_err;
      auto drvFail = [&]() {
        #pragma omp critical(pipe_err)
        if (!drv_err) drv_err = std::current_exception();
        pipe_abort.store(true);
      };
      std::atomic<size_t> written{0};   // tiles the writer has finished, in tile order
      auto writeTile = [&](const size_t f)   // the writer: outside the region, so a queued task never outlives it
      {
        while (written.load(std::memory_order_acquire) < f) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        try { if (!pipe_abort.load()) finishTile(*runs[f]); } catch (...) { drvFail(); }
        runs[f].reset();
        written.store(f + 1, std::memory_order_release);   // on every exit: the next writer never waits forever
      };
      #pragma omp parallel num_threads(n_threads)
      if (omp_get_thread_num() == 0)
      {
        // A team of one (OMP_DYNAMIC, a thread limit) has nobody to run the windows while the driver polls: then
        // one tile at a time, run by the driver's own taskwait and written by the driver -- the barrier's schedule.
        const bool solo = omp_get_num_threads() < 2;
        const size_t depth = solo ? 1 : (size_t)pipe_depth;
        size_t next = 0;   // the next tile to read
        auto issue = [&](const size_t upto)
        {
          for (; next < tiles.size() && next < upto && !pipe_abort.load(); ++next)
          {
            runs[next].reset(new TileRun);
            runs[next]->tk = next;
            loadTile(*runs[next]);
            admitTile(runs[next].get());
          }
        };
        onPipeDriver_() = !solo;   // [review] no taskloop on the driver (parTasks_); solo runs its windows itself anyway
        try
        {
          for (size_t fin = 0; fin < tiles.size() && !pipe_abort.load(); ++fin)
          {
            // [review] Backpressure on the writer: at most `depth` + 1 tiles whose windows are done wait for their write. Without
            // it the driver kept reading and admitting while the writer lagged (slow output storage), so finished spectra and
            // writer tasks spinning on their predecessor piled up without a bound. Never binds while the writer keeps pace:
            // by the read/write order in their logs it would not have bound in any of 38 pipelined timing
            // runs. Deadlock-free: writers 0..fin-1 are spawned, and they wait only on each other, in spawn order.
            while (written.load(std::memory_order_acquire) + depth < fin && !pipe_abort.load())
              if (!teamAssist_(8)) std::this_thread::sleep_for(std::chrono::milliseconds(1));   // help the writer meanwhile
            if (pipe_abort.load()) break;
            issue(fin + depth);
            if (!runs[fin]) break;   // the reader stopped (a failure)
            TileRun& R = *runs[fin];
            if (solo)
            {
              #pragma omp taskwait
            }
            while (R.left.load(std::memory_order_acquire) != 0 && !pipe_abort.load())
              if (!teamAssist_(8)) std::this_thread::sleep_for(std::chrono::milliseconds(1));   // help the writer meanwhile
            if (pipe_abort.load()) break;
            if (solo) writeTile(fin);
            else
            {
              auto* wt = &writeTile;
              #pragma omp task default(none) firstprivate(wt, fin)
              (*wt)(fin);
            }
          }
        }
        catch (...) { drvFail(); }
        onPipeDriver_() = false;
      }   // the region's barrier: every admitted window and every writer has finished
      if (drv_err) std::rethrow_exception(drv_err);
      if (worker_err) std::rethrow_exception(worker_err);
    }

    win_phase.reset();   // [perf-instr] WINDOW_LOOP ends here: the tile loop's own phases were enclosed in it, the final write is not
    writeLogInfo_("[soa-rt] off-grid " + String(rtOffGrid().load())
                  + " (of which strictly between two frames " + String(rtStrictlyBetween().load()) + ")");
    if (rtOffGrid().load() != 0)
      writeLogWarn_("[soa-rt] " + String(rtOffGrid().load()) + " traces had an RT that is not one of "
                    "the store's frame times (ElutionPeakDetection skipped, or a span trimmed away from "
                    "its smoothed apex); their RT was snapped to the nearest frame. Zero is the "
                    "output-identical case -- a nonzero count here means the digest will differ.");
    else writeLogInfo_("[soa-rt] every trace RT was exactly a frame time (0 snapped).");
    writeLogInfo_("All windows done (admission high-water " + String(adm_hwm.load()) + " of "
                  + String(n_conc) + " allowed concurrent windows)." + clk_() + rss_() + mem_()); // [mem]
#ifdef __GLIBC__
    if (trim_on) trimLogged_("window loop -> write");
#endif

    if (writer)
    {
      { Phase _ph("WRITE(mzML)"); writer.reset(); }   // the count patched, the index and footer written
      if (std::rename(out_part.c_str(), out.c_str()) != 0)
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out, "cannot rename the completed " + out_part);
    }
    else
    {
      // the bulk path (mzPeak output, one tile): as before, whole
      { Phase _ph("ASSEMBLE(PeakMap)"); out_exp.setSpectra(std::move(all_out)); }
      n_out = out_exp.size();
      for (Size i = 0; i < out_exp.size(); ++i) out_exp[i].getDataProcessing().push_back(out_dp);
      { Phase _ph("WRITE(mzPeak)"); FileHandler().storeExperiment(out, out_exp, {out_type}, log_type_); }
    }
    writeLogInfo_("Wrote " + String(n_out) + " pseudo-MS2 spectra to " + out
                  + " (" + String(FileTypes::typeToName(out_type)) + ")");
    if (twin_im >= 0.0) writeLogInfo_("[twins] " + String(n_twins) + " spectra folded into a twin (assembly:twin_im_tolerance " + String(twin_im) + ")");
    if (twin_im >= 0.0 && twin_rt_tol > 0.0)   // [twin-rt] the before/after of the near-frame relation
      writeLogInfo_("[twins] assembly:twin_rt_frames: " + String(twin_st[0]) + " twin pairs at one MS1 frame, " + String(twin_st[1])
                    + " same-monoisotope pairs at nearby frames joined, " + String(twin_st[2]) + " left apart at an RT cell edge; "
                    + String(n_twins - std::min<size_t>(n_twins, twin_st[3])) + " of the " + String(n_twins) + " folded spectra need a near-frame link (the same frame alone folds "
                    + String(twin_st[3]) + "); " + String(twin_st[4]) + " clusters span more than the tolerance (chained)");
    if (post_fusion_cap_)
      writeLogInfo_("[recap] assembly:post_fusion_cap: " + String(recap[0]) + " of " + String(recap[3]) + " spectra re-capped to " + String(max_frags)
                    + " peaks after the twin fold, " + String(recap[1]) + " peaks removed (largest spectrum before: " + String(recap[2]) + " peaks)");
    {   // [merge] how much splitting the two merges removed
      auto& m = mergeCount_();
      writeLogInfo_("[merge] im-slice merge (trace:im_slice_merge " + (!integer_detector ? String("n/a, openms detector") : im_slice_merge_ >= 0.0 ? String(im_slice_merge_) : String("off")) + "): " + String(m[0].load())
                    + " picked MS2 peaks merged away in " + String(m[1].load()) + " frames; assembly same-bin merge (" + (merge_same_bin_ ? "on" : "off") + "): "
                    + String(m[2].load()) + " fragments merged away");
    }
    if (frag_int_ != 0)   // [fragint] what the fragment-intensity mode did (gate-passing and rescued fragments, before the merge and the cap)
      writeLogInfo_("[fragint] assembly:frag_intensity " + (frag_int_ == 1 ? "core (" + String(2 * frag_core_half_ + 1) + " frames around the precursor apex)" : String("sum (the precursor's span)"))
                    + ": " + String(fragIntCount_()[0].load()) + " fragments took the summed amplitude, " + String(fragIntCount_()[1].load())
                    + " of them had no point there and took their point nearest the precursor apex");
    if (foreign_w_ < 1.0)
    {   // [foreign] what assembly:foreign_precursor_weight did
      auto& c = foreignCount_();
      const long long nm = c[7].load();
      writeLogInfo_("[foreign] assembly:foreign_precursor_weight " + String(foreign_w_) + " (" + String(foreign_ppm_) + " ppm, 1/K0 " + String(foreign_dim_)
                    + ", min_iso " + String(foreign_min_iso_) + "): of " + String(c[8].load()) + " precursor-window assemblies with fragments, "
                    + String(c[0].load()) + " had a foreign candidate (" + String(c[4].load()) + " candidates skipped as the precursor's own ion); "
                    + String(c[1].load()) + " foreign envelopes matched; "
                    + String(c[3].load()) + " fragments " + (foreign_w_ > 0.0 ? "attenuated" : "removed") + " in " + String(c[2].load()) + " assemblies, "
                    + String(c[5].load()) + " of which then fell below assembly:min_fragments; mean m/z error of the " + String(nm)
                    + " matches " + String(nm ? (double)c[6].load() / 1000.0 / (double)nm : 0.0, false) + " ppm");
    }
    if (fuse_ppm_ > 0.0)
    {   // [fuse] what assembly:fragment_fuse_ppm did (per precursor-window assembly, before min_fragments and the cap)
      auto& c = fuseCount_();
      const long long g = c[3].load();
      writeLogInfo_("[fuse] assembly:fragment_fuse_ppm " + String(fuse_ppm_) + " (" + (fuse_linkage_ ? "linkage" : "anchored") + "): "
                    + String(c[1].load()) + " fragments in, " + String(c[2].load()) + " peaks out, in " + String(c[0].load())
                    + " precursor-window assemblies; " + String(g) + " groups of >= 2 fused, mean size "
                    + String(g ? (double)c[4].load() / (double)g : 0.0, false) + "; " + String(c[5].load())
                    + " assemblies fell below assembly:min_fragments");
    }
    {   // [pspan] the precursor MS1-span histogram, for choosing gate:faint_max_frames per data set
      auto& h = pspanHist_(); long long tot = 0; for (int i = 0; i < 8; ++i) tot += h[i].load();
      writeLogInfo_("[pspan] precursors by MS1 span (window frames): 1:" + String(h[0].load()) + " 2:" + String(h[1].load()) + " 3:" + String(h[2].load())
                    + " 4:" + String(h[3].load()) + " 5-6:" + String(h[4].load()) + " 7-10:" + String(h[5].load()) + " 11-20:" + String(h[6].load())
                    + " >20:" + String(h[7].load()) + " total " + String(tot) + (adaptive_ ? " faint(<=" + String(faint_frames_) + "): " + String(h[8].load()) : String("")));
    }
    if (anchor_rescue_)
    {   // [anchor] what the rescue did (counted at assembly: before the twin fold)
      auto& a = anchorCount_();
      char b[640];
      snprintf(b, sizeof b, "[anchor] %lld of %lld non-faint precursors anchored (>= %d fragments at r >= %g), %lld candidates; %lld fragments "
                            "rescued (r >= %g; rejected at the point floor %lld, overlap guard %lld, correlation cut %lld) in %lld precursors, "
                            "%.2f per anchored precursor, %lld more beyond %d per spectrum; %lld rescued peaks in %lld written spectra after the "
                            "same-bin merge and the cap (%.2f per such spectrum)",
               a[1].load(), a[0].load(), anchor_min_count_, anchor_min_corr_, a[2].load(), a[3].load(), anchor_rescue_corr_, a[8].load(),
               a[9].load(), a[10].load(), a[5].load(), a[1].load() ? (double)a[3].load() / a[1].load() : 0.0, a[4].load(),
               anchor_max_rescued_, a[6].load(), a[7].load(), a[7].load() ? (double)a[6].load() / a[7].load() : 0.0);
      writeLogInfo_(b);
    }
    report_phases_(phase_clock_());
    return EXECUTION_OK;
  }
};

/// @endcond

extern char** environ;

int main(int argc, const char** argv)
{
  // The tool was SpeXtractor (variables SPEXTRACTOR_*) and then DIAspeXtractor (DIASPEXTRACTOR_*); it is DIAspeXtract
  // now, and every variable is DIASPEXTRACT_* (so are the three the patched libOpenMS reads). An old name would be
  // ignored without a word -- a chain's A/B would measure the default -- so any variable under either old prefix refuses
  // the run before anything else happens, naming the variable and its DIASPEXTRACT_ name. There is no dual read.
  static const struct { const char* prefix; size_t len; } old_prefixes[] = {{"SPEXTRACTOR_", 12}, {"DIASPEXTRACTOR_", 15}};
  bool old_name = false;
  for (char** e = environ; e != nullptr && *e != nullptr; ++e)
  {
    for (const auto& op : old_prefixes)
    {
      if (std::strncmp(*e, op.prefix, op.len) != 0) continue;
      const int n = (int)std::strcspn(*e, "=");                 // the variable name
      const int m = n - (int)op.len;                            // the part after the old prefix
      std::fprintf(stderr, "Error: %.*s is set, but the tool is now DIAspeXtract and reads only DIASPEXTRACT_* variables "
                           "(SpeXtractor's SPEXTRACTOR_* and DIAspeXtractor's DIASPEXTRACTOR_* are refused, not ignored): "
                           "rename it to DIASPEXTRACT_%.*s or unset it. Nothing was run.\n", n, *e, m, *e + op.len);
      old_name = true;
    }
  }
  if (old_name) return TOPPBase::ILLEGAL_PARAMETERS;
  // TOPPBase's update check sends this tool's name and version to an OpenMS REST endpoint on every run and,
  // offline, prints Qt QIODevice errors on the stderr the harness parses (build.yml asserts --help stays clean).
  // Set only if unset, so a user's export still decides.
  ::setenv("OPENMS_DISABLE_UPDATE_CHECK", "ON", 0);
  TOPPDIAspeXtract tool;
  return tool.main(argc, argv);
}
