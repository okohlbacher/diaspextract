// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// Copyright (c) 2026, Oliver Kohlbacher and contributors
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include "MzPeakStreamLoad.h" // [mzpeak] streaming .mzpeak input (DIASPEXTRACT_WITH_MZPEAK)
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include "OmsString.h" // String(x) on top of std::string (OpenMS 3.7 has no OpenMS::String)

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/S3InputSource.h>   // [s3] -in s3://...
#include <OpenMS/FORMAT/DATAACCESS/HTTPInputSource.h> // redactURL: a URI may carry credentials
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SysInfo.h> // getAvailablePhysicalCores: -threads default and cap
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include "CgroupMem.h"    // [dyn-mem] cgroup headroom net of page cache
#include "TdfLoad.h"      // loadTdfCalibration: sqlite3, tool-only
#include "NoImsPath.h"    // [noims] the separate extraction path for input without ion mobility
#include "DnoiseMs1.h"    // [dnoise] the dnoise v0.1.0 MS1 port (header-only; tests/test_dnoise_ms1.cpp compiles the same code)
#include "PortableMath.h" // [xplat] exp/log/log1p/pow with the same bits on every platform (the libm's differ)
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>     // [faims]
#include <OpenMS/IONMOBILITY/IMDataConverter.h> // [faims] splitByFAIMSCV
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/VersionInfo.h> // verboseVersion_ carries the OpenMS build alongside ours

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
#include <deque>
#include <mutex>
#include <thread>
#include <map>
#include <string>
#include <vector>
#include <cstdlib> // getenv
#include <cstring> // memcpy
#include <optional>
#include <memory> // [outer] unique_ptr
#include <time.h> // clock_gettime(CLOCK_THREAD_CPUTIME_ID): [dnoise] filter CPU
#ifdef __GLIBC__
  #include <malloc.h> // malloc_trim [trim]
#endif
#include <filesystem> // rename: replaces an existing -out on every platform (std::rename does not on Windows)
#include <system_error>
#include <cerrno>
#ifndef _WIN32
  #include <fcntl.h> // open/fsync of the completed .part before its rename
  #include <unistd.h>
#endif

// The eight OpenMS patches in patches/ are build prerequisites. This is the compile-time check for
// openms-masstrace-move.patch: without defaulted moves an rvalue MassTrace binds to the copy constructor (not noexcept), and
// every "move" of a trace in the band/split gathering is a deep copy of its points.
static_assert(std::is_nothrow_move_constructible_v<OpenMS::MassTrace> && std::is_nothrow_move_assignable_v<OpenMS::MassTrace>,
              "OpenMS MassTrace has no noexcept move operations: apply patches/openms-masstrace-move.patch to the OpenMS tree");
#include <exception> // exception_ptr
#include <omp.h>
#include <iomanip>
#include <sstream>

// [win] Windows portability: Win32 stand-ins for the POSIX and /proc probes (cpuSeconds, procKb, threadCpu_).
// Included after every other header, so the Win32 macros reach this file only. The Windows compiler is clang-cl (see
// CMakeLists.txt: cl's OpenMP has no taskloop).
#ifdef _WIN32
  #define NOMINMAX
  #define WIN32_LEAN_AND_MEAN
  #define NOGDI
  #include <windows.h>
  #include <psapi.h>  // GetProcessMemoryInfo (kernel32's K32 entry point: no extra library)
  #include <stdlib.h> // _environ, _putenv_s
/**
  @brief Sum of two Windows FILETIME values (100 ns ticks), in seconds.

  @param[in] a First time
  @param[in] b Second time
  @return (a + b) in seconds
*/
static inline double spxFileTimeSec(FILETIME a, FILETIME b) // two FILETIMEs (100 ns ticks) summed, in seconds
{
  return ((double)(((unsigned long long)a.dwHighDateTime << 32) | a.dwLowDateTime)
          + (double)(((unsigned long long)b.dwHighDateTime << 32) | b.dwLowDateTime))
         * 1e-7;
}
#else
  #include <sys/resource.h>
#endif
#ifdef __APPLE__
  #include <mach/mach.h> // host_statistics64, task_info: procKb_'s stand-ins for /proc on macOS
#endif

// The release version. The standalone CMake build passes it from project(VERSION); the in-tree
// OpenMS build does not, so the fallback here is the same number.
#ifndef DIASPEXTRACT_VERSION
  #define DIASPEXTRACT_VERSION "1.9.4"
#endif

using namespace OpenMS;
using namespace std;
using DIAspeXtract::String;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_DIAspeXtract DIAspeXtract

  @brief Builds pseudo-MS/MS ("pseudo-DDA") spectra from DIA data for DDA search engines.

  Reconstructs one DDA-like MS2 spectrum per precursor from data-independent acquisition, so that any DDA search
  engine, including open and blind searches, can search DIA data without a spectral library. It does not perform the
  database search, FDR control, library building or quantification.

  Two kinds of input take two separate paths:
  - <b>diaPASEF</b> (Bruker timsTOF, ion-mobility DIA; a .d directory, mzML or mzPeak). MS1 frames are denoised and
    centroided along ion mobility; precursor and fragment mass traces are detected on the instrument's flight-time
    bins; isotope traces at the same ion mobility form precursors (charge 1-5, checked by ion mobility and by the
    isotope pattern). A fragment joins a precursor when the two co-localize in ion mobility and retention time and
    their elution profiles correlate. The run is processed in retention-time tiles, the isolation windows in parallel.
  - <b>DIA without ion mobility</b> (Orbitrap, SWATH; centroided mzML). Precursors are @ref Biosaur2Algorithm
    isotope-pattern features of the MS1 scans; fragments are the mass traces of each isolation window whose elution
    profiles correlate with the precursor's. None of the diaPASEF options apply to this path.

  The output is mzML (or mzPeak in a build with mzPeak support), written as &lt;out&gt;.part and renamed when complete.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_DIAspeXtract.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_DIAspeXtract.html
*/

/// @cond TOPPCLASSES

namespace
{
/**
  @brief Seconds since tool start.

  [perf]

  @return Seconds since the first call
*/
static double phaseClock()
{
  static const auto t0 = std::chrono::steady_clock::now();
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
/**
  @brief The current phaseClock() as a log suffix " [t=<s>s]".

  @return The stamp text
*/
static std::string clockStamp()
{
  char b[32];
  snprintf(b, sizeof b, " [t=%.0fs]", phaseClock());
  return b;
}
/**
  @brief A number for a log line: four significant digits (OpenMS's String(double) prints all 17).

  Header stamps keep
  full precision, so a value can be passed back exactly.

  @param[in] x The value
  @return x printed with "%.4g"
*/
static String numStr(double x)
{
  char b[32];
  snprintf(b, sizeof b, "%.4g", x);
  return b;
}

/// @brief [perf-instr] Per-phase wall/CPU/RSS accounting, printed unconditionally at the end.
struct PhaseStat
{
  double wall = 0;
  double cpu = 0;
  long rss_end_mb = 0;
  int n = 0;
  bool nested = false;
};
/**
  @brief [perf-instr] A phase opened inside another (the tile loop's read/sort/write inside WINDOW_LOOP)
  is marked enclosed and left out of the TOTAL.

  Phases are opened on the master thread, and with pipelined
  tiles also by the writer task while the master reads (atomic depth, locked phaseAdd).

  @return The process-wide nesting depth of open phases
*/
inline std::atomic<int>& phaseDepth()
{
  static std::atomic<int> d {0};
  return d;
}
/**
  @brief [task-err] An exception escaping an OpenMP task terminates the process with no message.

  Task bodies catch
  everything and hand the first exception here; the code after the loop rethrows it once the tasks have joined.

  @param[in,out] first Receives the current exception if it is still empty (set under a critical section)
*/
inline void taskErr(std::exception_ptr& first)
{
#pragma omp critical(task_err)
  if (! first) { first = std::current_exception(); }
}
/**
  @brief [task-err] e2e falsifier only: DIASPEXTRACT_TASK_THROW=<site> throws inside that task, so the suite can show the run
  ends with an error message and a non-zero exit instead of an abort.

  @param[in] site Name of the task site
  @throws Exception::InvalidValue if DIASPEXTRACT_TASK_THROW names @p site
*/
inline void taskFault(const char* site)
{
  static const char* want = std::getenv("DIASPEXTRACT_TASK_THROW");
  if (want && String(want) == site)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("[task-err] injected failure in task ") + site, site);
  }
}
/**
  @brief [pipe] Helper tasks of teamFor created and not yet started, over every call: capped at the team size, so
  stale helpers (started after their loop ran out) cannot pile up towards the libgomp cliff (64 x nthreads).

  @return The process-wide counter of queued helper tasks
*/
inline std::atomic<int>& teamHelpersQueued()
{
  static std::atomic<int> q {0};
  return q;
}
/**
  @brief [pipe] A loop that the pipelined tile loop runs on the window team's driver (the tile reader and writer).

  Outside a parallel region the call sites keep their `omp parallel for`. Inside the team a `parallel for`
  would open a NESTED team -- fresh threads on the cores the windows use -- so the loop runs here instead:
  the caller pulls `chunk` indices at a time from a shared cursor, and up to one helper task per team
  thread pulls with it once a thread is free. The caller never reaches a task scheduling point, so it is
  never handed a window task (libgomp lets a thread waiting at a taskgroup end run ANY child of its task),
  and it never waits for a helper to START: a helper that starts after the cursor ran out touches nothing
  but its own reference to the cursor, so only helpers that claimed a chunk are waited for (seq_cst: a
  helper raises `busy` before it claims). f(i, slot): `slot` = the puller's omp_get_thread_num(), unique
  among the concurrent pullers (a thread pulls one chunk at a time and never nests a pull), for per-slot
  scratch where a thread number was used. Every caller writes index-addressed results only, so the output
  is the parallel for's. The first exception is rethrown here, after every claimed chunk has finished.
  Helper tasks queue behind the band and score tasks and barely run while windows hold the threads, so the loop is
  also published for teamAssist: window code pulls a chunk at its break points instead of opening a nested team.
*/
struct PullJob
{
  std::atomic<long> cur {0};
  std::atomic<int> busy {0};
  std::exception_ptr err;
  long n = 0, chunk = 1;
  void (*call)(void*, long, int) = nullptr;
  void* fp = nullptr;
  /**
    @brief Pull chunks from the shared cursor and run the loop body on them.

    Up to max_chunks chunks (all with max_chunks < 0), returns how many it claimed; `busy` is raised BEFORE
    the first claim (see above). An exception of the body is stored in `err` and ends the loop (the cursor is moved to the end).

    @param[in] max_chunks Maximum number of chunks to claim (negative: all that remain)
    @return Number of chunks claimed
  */
  long pull(const long max_chunks)
  {
    const int slot = omp_get_thread_num();
    busy.fetch_add(1);
    long a = 0, k = 0;
    for (; (max_chunks < 0 || k < max_chunks) && (a = cur.fetch_add(chunk)) < n; ++k)
    {
      try
      {
        for (long i = a, e = std::min(a + chunk, n); i < e; ++i)
        {
          call(fp, i, slot);
        }
      }
      catch (...)
      {
        taskErr(err);
        cur.store(n);
      }
    }
    busy.fetch_sub(1);
    return k;
  }
};
/**
  @brief Mutex guarding pullJobs().

  @return The process-wide mutex
*/
inline std::mutex& pullMu()
{
  static std::mutex m;
  return m;
}
/**
  @brief The live teamFor loops published for teamAssist().

  @return The process-wide list (guarded by pullMu())
*/
inline std::vector<std::shared_ptr<PullJob>>& pullJobs()
{
  static std::vector<std::shared_ptr<PullJob>> v;
  return v;
}
/**
  @brief Number of live teamFor loops; teamAssist() reads it before taking the lock.

  @return The process-wide counter
*/
inline std::atomic<int>& pullLive()
{
  static std::atomic<int> n {0};
  return n;
}
/**
  @brief [pipe] Help a live teamFor loop (the pipelined tile reader or writer) with up to max_chunks chunks; true if
  it did any.

  One relaxed load when none is live. Called inside the window team only (slot = thread number).

  @param[in] max_chunks Maximum number of chunks to run (negative: all that remain)
  @return true if any chunk was run
*/
inline bool teamAssist(const long max_chunks)
{
  if (pullLive().load(std::memory_order_relaxed) == 0) { return false; }
  std::shared_ptr<PullJob> j; // the copy keeps the job alive; its loop body is the caller's until it returns
  {
    std::lock_guard<std::mutex> g(pullMu());
    for (const auto& p : pullJobs())
    {
      if (p->cur.load(std::memory_order_relaxed) < p->n)
      {
        j = p;
        break;
      }
    }
  }
  return j && j->pull(max_chunks) > 0;
}
/**
  @brief [overlap] teamFor in two halves: teamForStart publishes the loop and queues its helpers and returns at once (the
  caller may do other work meanwhile, e.g. read the next frames); teamForFinish pulls what is left, unpublishes the
  loop and waits for every claimed chunk.

  `f` must outlive teamForFinish (the helpers call it through the job).

  @param[in] n Number of indices
  @param[in] chunk Indices per chunk
  @param[in] f Loop body f(i, slot); must outlive teamForFinish
  @return The published job, to be passed to teamForFinish()
*/
template<class F>
static std::shared_ptr<PullJob> teamForStart(long n, long chunk, F& f)
{
  auto st = std::make_shared<PullJob>();
  st->n = std::max(0L, n);
  st->chunk = chunk;
  st->fp = (void*)&f;
  st->call = [](void* p, long i, int slot) { (*static_cast<F*>(p))(i, slot); };
  if (n <= 0) { return st; }
  {
    std::lock_guard<std::mutex> g(pullMu());
    pullJobs().push_back(st);
  }
  ++pullLive();
  const int team = std::max(1, omp_get_num_threads());
  const int nt = (int)std::min<long>((n + chunk - 1) / chunk, (long)team);
  for (int t = 1; t < nt && teamHelpersQueued().load() < team; ++t)
  {
    ++teamHelpersQueued();
#pragma omp task default(none) firstprivate(st)
    {
      --teamHelpersQueued();
      st->pull(-1);
      while (teamAssist(-1)) {}
    } // [overlap] then any other live loop (the decode beside a pick): the cap held its helpers back
  }
  return st;
}
/**
  @brief Second half of teamFor: pull the remaining chunks, unpublish the loop and wait for every claimed chunk.

  @param[in] st The job returned by teamForStart()
  @par Exceptions
    Rethrows the first exception thrown by the loop body, of whatever type it was.
*/
static void teamForFinish(const std::shared_ptr<PullJob>& st)
{
  if (st->n <= 0) { return; }
  st->pull(-1);
  {
    std::lock_guard<std::mutex> g(pullMu());
    auto& v = pullJobs();
    v.erase(std::remove(v.begin(), v.end(), st), v.end());
  }
  --pullLive();
  while (st->busy.load() != 0)
  {
    std::this_thread::yield();
  }
  if (st->err) { std::rethrow_exception(st->err); }
}
/**
  @brief Run f(i, slot) for i in [0, n) on the caller and helper tasks of the current team (see PullJob).

  @param[in] n Number of indices
  @param[in] chunk Indices per chunk
  @param[in] f Loop body f(i, slot)
  @par Exceptions
    Rethrows the first exception thrown by the loop body, of whatever type it was.
*/
template<class F>
static void teamFor(long n, long chunk, F&& f)
{
  if (n <= 0) { return; }
  teamForFinish(teamForStart(n, chunk, f));
}
/**
  @brief True on the pipelined tile loop's driver while it drives.

  Every window task is a child of the driver's
  implicit task, so a taskloop's taskgroup end on the driver (libgomp then runs ANY queued child of the waiting task)
  could hand it a whole window, stalling the reads and admissions.
  parTasks keeps every pass the driver reaches inline instead (nt = 1 is forTasks's serial branch: the same result).

  @return Reference to this thread's flag
*/
inline bool& onPipeDriver()
{
  static thread_local bool v = false;
  return v;
}
/**
  @brief [serial] Task count for one of a window's data-parallel passes over n items: at most 16 and none smaller than
  min_per; 1 (run inline) outside a parallel region and on the pipe driver.

  A bounded COUNT, never a grainsize over
  data: libgomp runs a whole taskloop inline once the team's task count passes 64 x nthreads ([gomp-cliff]).

  @param[in] n Number of items
  @param[in] min_per Minimum number of items per task
  @return The task count (>= 1)
*/
inline long parTasks(Size n, Size min_per)
{
  const long cap = omp_in_parallel() && ! onPipeDriver() ? std::min(16, omp_get_num_threads()) : 1;
  return std::max(1L, std::min(cap, (long)(n / std::max<Size>(min_per, 1))));
}
/**
  @brief [serial] body(c) for c in [0, nt): nt tasks of the current team (the creating thread waits in the taskloop's
  taskgroup and runs its own children, so a busy pool degrades to the serial loop), or inline for nt <= 1.

  The first exception is rethrown once every task has joined.

  @param[in] nt Number of tasks
  @param[in] body Called as body(c) for each task index c
  @par Exceptions
    Rethrows the first exception thrown by @p body, of whatever type it was.
*/
template<class F>
void forTasks(long nt, F&& body)
{
  if (nt <= 1)
  {
    body(0L);
    return;
  }
  std::exception_ptr err;
#pragma omp taskloop num_tasks(nt) default(shared)
  for (long c = 0; c < nt; ++c)
  {
    try
    {
      body(c);
    }
    catch (...)
    {
      taskErr(err);
    }
  }
  if (err) { std::rethrow_exception(err); }
}
/**
  @brief [serial] Merge-path co-rank: how many of the first r elements of the STABLE merge of a[0, m) and b[0, k) come
  from a (a first on ties, as std::merge takes them).

  @param[in] r Output rank
  @param[in] a First sorted range
  @param[in] m Length of @p a
  @param[in] b Second sorted range
  @param[in] k Length of @p b
  @param[in] less Comparator
  @return The number of elements taken from @p a
*/
template<class T, class Less>
Size coRank(Size r, const T* a, Size m, const T* b, Size k, const Less& less)
{
  Size lo = r > k ? r - k : 0, hi = std::min(r, m);
  while (lo < hi)
  {
    const Size i = lo + (hi - lo) / 2, j = r - i; // i < m and 1 <= j <= k here
    if (! less(b[j - 1], a[i])) { lo = i + 1; }
    else
    {
      hi = i;
    }
  }
  return lo;
}
/**
  @brief [serial] Sort a[0, n) as nt tasks: nt chunks std::sort-ed, then pairwise std::merge levels, each pair's output
  cut at merge-path co-ranks into pieces of about n / nt, ping-ponging with one scratch buffer of n elements.

  The result equals std::sort's whenever `less` is a TOTAL order (no two elements equivalent: the sorted sequence
  is then unique). Under a weaker order equivalent elements may come out in another order.

  @param[in,out] a The elements to sort
  @param[in] n Number of elements
  @param[in] less Comparator
  @param[in] nt Number of tasks (<= 1: std::sort)
  @throws std::bad_alloc if the scratch buffer cannot be allocated
*/
template<class T, class Less>
void parallelSort(T* a, Size n, const Less& less, long nt)
{
  static_assert(std::is_trivially_copyable<T>::value, "parallelSort_ moves elements through raw scratch memory");
  if (nt <= 1 || n < 2 * (Size)nt)
  {
    std::sort(a, a + n, less);
    return;
  }
  vector<Size> run(nt + 1);
  for (long c = 0; c <= nt; ++c)
  {
    run[c] = (Size)c * n / (Size)nt;
  }
  forTasks(nt, [&](long c) { std::sort(a + run[c], a + run[c + 1], less); });
  // uninitialised scratch, first touched by the merge tasks (an implicit-lifetime type: malloc creates the objects)
  std::unique_ptr<T, void (*)(T*)> buf(static_cast<T*>(std::malloc(n * sizeof(T))), [](T* p) { std::free(p); });
  if (! buf) { throw std::bad_alloc(); }
  T* src = a;
  T* dst = buf.get();
  struct Piece
  { Size s, mid, e, r0, r1; }; // merge src[s, mid) with src[mid, e); output ranks [r0, r1)
  while (run.size() > 2)
  {
    const Size nr = run.size() - 1;
    vector<Size> nxt {0};
    vector<Piece> pcs;
    for (Size p = 0; p < nr; p += 2)
    {
      const Size s = run[p], mid = run[p + 1], e = p + 1 < nr ? run[p + 2] : mid; // an odd last run is copied
      const Size len = e - s, np = std::max<Size>(1, (len * (Size)nt + n - 1) / n);
      for (Size q = 0; q < np; ++q)
      {
        pcs.push_back({s, mid, e, q * len / np, (q + 1) * len / np});
      }
      nxt.push_back(e);
    }
    forTasks((long)pcs.size(), [&](long c) {
      const Piece& p = pcs[(Size)c];
      const T* run_a = src + p.s;
      const T* run_b = src + p.mid;
      const Size m = p.mid - p.s, k = p.e - p.mid;
      const Size i0 = coRank(p.r0, run_a, m, run_b, k, less), i1 = coRank(p.r1, run_a, m, run_b, k, less);
      std::merge(run_a + i0, run_a + i1, run_b + (p.r0 - i0), run_b + (p.r1 - i1), dst + p.s + p.r0, less);
    });
    std::swap(src, dst);
    run.swap(nxt);
  }
  if (src != a)
  {
    forTasks(nt, [&](long c) {
      const Size lo = (Size)c * n / (Size)nt, hi = (Size)(c + 1) * n / (Size)nt;
      std::copy(src + lo, src + hi, a + lo);
    });
  }
}

/**
  @brief [outer] std::sort's result, computed in parallel, for a comparator that is a strict TOTAL order on the elements
  present (no two equivalent): the sorted sequence is then unique, so any correct sort returns exactly what std::sort
  returns.

  Chunks are sorted concurrently and merged pairwise; each merge is cut into equal output pieces by merge path
  (exactly std::merge's sequence). Serial std::sort when small, on one thread, or inside a parallel region. Callers
  guarantee the total order by construction (a unique last key) or check the result (sortExact).

  @param[in,out] v The elements to sort
  @param[in] less A strict total order on the elements
*/
template<typename T, typename Less>
void parallelSortTotal(std::vector<T>& v, Less less)
{
  const Size n = v.size();
  const int nt = omp_get_max_threads();
  if (n < ((Size)1 << 15) || nt < 2 || omp_in_parallel())
  {
    std::sort(v.begin(), v.end(), less);
    return;
  }
  Size nc = 1; // a power of two, about one chunk per thread
  while (nc < (Size)nt && n / (2 * nc) >= 4096)
  {
    nc *= 2;
  }
  auto bound = [n, nc](Size c) { return c * n / nc; }; // chunk c is [bound(c), bound(c + 1))
  std::vector<T> buf(n);
  T* src = v.data();
  T* dst = buf.data();
#pragma omp parallel for schedule(dynamic, 1)
  for (long c = 0; c < (long)nc; ++c)
  {
    std::sort(src + bound((Size)c), src + bound((Size)c + 1), less);
  }
  for (Size w = 1; w < nc; w *= 2) // runs of w chunks, merged pairwise
  {
    const Size nm = nc / (2 * w), np = std::max<Size>(1, (2 * (Size)nt + nm - 1) / nm); // merges, pieces per merge
#pragma omp parallel for schedule(dynamic, 1)
    for (long q = 0; q < (long)(nm * np); ++q)
    {
      const Size m = (Size)q / np, p = (Size)q % np;
      const Size a0 = bound(2 * m * w), b0 = bound((2 * m + 1) * w), b1 = bound((2 * m + 2) * w);
      const T* a = src + a0;
      const T* b = src + b0;
      const Size na = b0 - a0, nb = b1 - b0, k0 = p * (na + nb) / np, k1 = (p + 1) * (na + nb) / np;
      const Size i0 = coRank(k0, a, na, b, nb, less), i1 = coRank(k1, a, na, b, nb, less);
      std::merge(a + i0, a + i1, b + (k0 - i0), b + (k1 - i1), dst + a0 + k0, less);
    }
    std::swap(src, dst);
  }
  if (src != v.data()) { v.swap(buf); }
}

/**
  @brief [outer] Sorts @p v exactly as std::sort(v.begin(), v.end(), less) does, in parallel whenever that is provably the
  same result: a permutation is sorted in parallel (index as the last key, so a total order) and, when no two
  neighbours are equivalent under @p less (every element then has one place, and std::sort's sequence is that one),
  gathered.

  With an equivalent pair, @p v -- still untouched -- goes through std::sort, so the order of equal keys is
  std::sort's own. Returns the number of equal-key neighbours of the result (the same count either way).

  @param[in,out] v The elements to sort
  @param[in] less Comparator
  @return The number of equal-key neighbours in the sorted result
*/
template<typename T, typename Less>
Size sortExact(std::vector<T>& v, Less less)
{
  const Size n = v.size();
  auto serial = [&]() {
    std::sort(v.begin(), v.end(), less);
    Size t = 0;
    for (Size i = 1; i < n; ++i)
    {
      if (! less(v[i - 1], v[i])) { ++t; }
    }
    return t;
  };
  if (n < ((Size)1 << 15) || n > (Size)std::numeric_limits<uint32_t>::max() || omp_get_max_threads() < 2 || omp_in_parallel()) { return serial(); }
  std::vector<uint32_t> perm(n);
#pragma omp parallel for schedule(static)
  for (long i = 0; i < (long)n; ++i)
  {
    perm[(Size)i] = (uint32_t)i;
  }
  parallelSortTotal(perm, [&v, &less](uint32_t a, uint32_t b) {
    if (less(v[a], v[b])) { return true; }
    if (less(v[b], v[a])) { return false; }
    return a < b;
  });
  Size ties = 0;
#pragma omp parallel for schedule(static) reduction(+ : ties)
  for (long i = 1; i < (long)n; ++i)
  {
    if (! less(v[perm[(Size)i - 1]], v[perm[(Size)i]])) { ++ties; }
  }
  if (ties) { return serial(); }
  std::vector<T> out(n);
#pragma omp parallel for schedule(static)
  for (long i = 0; i < (long)n; ++i)
  {
    out[(Size)i] = v[perm[(Size)i]];
  }
  v.swap(out);
  return 0;
}

/**
  @brief The accumulated statistics per phase name.

  @return The process-wide map
*/
inline std::map<std::string, PhaseStat>& phaseStats()
{
  static std::map<std::string, PhaseStat> m;
  return m;
}
/**
  @brief Phase names in the order they were first closed, for reportPhases().

  @return The process-wide list
*/
inline std::vector<std::string>& phaseOrder()
{
  static std::vector<std::string> v;
  return v;
}
/**
  @brief CPU time (user + system) of the whole process, over all threads.

  @return Seconds, 0 if unavailable
*/
static double cpuSeconds()
{
#ifdef _WIN32
  FILETIME c, e, k, u; // kernel + user over all threads, as RUSAGE_SELF
  return GetProcessTimes(GetCurrentProcess(), &c, &e, &k, &u) ? spxFileTimeSec(k, u) : 0.0;
#else
  struct rusage ru; // RUSAGE_SELF sums ALL threads -> parallel efficiency
  if (getrusage(RUSAGE_SELF, &ru) != 0) { return 0.0; }
  return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6 + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
#endif
}
/**
  @brief The kB value of `key` in a /proc file; -1 where it is absent (Linux; Windows answers its keys from Win32; else -1).

  @param[in] path The /proc file (Linux)
  @param[in] key The key including its colon, e.g. "VmRSS:"
  @return The value in kB, or -1
*/
static long long procKb(const char* path, const char* key)
{
#ifdef _WIN32
  // [win] the three keys read here, from Win32: VmRSS/VmHWM = working set and its peak, MemAvailable = free physical
  (void)path;
  const std::string k(key);
  if (k == "MemAvailable:")
  {
    MEMORYSTATUSEX m {};
    m.dwLength = sizeof m;
    return GlobalMemoryStatusEx(&m) ? (long long)(m.ullAvailPhys / 1024) : -1;
  }
  PROCESS_MEMORY_COUNTERS pm {};
  if (! GetProcessMemoryInfo(GetCurrentProcess(), &pm, sizeof pm)) { return -1; }
  return k == "VmRSS:" ? (long long)(pm.WorkingSetSize / 1024) : k == "VmHWM:" ? (long long)(pm.PeakWorkingSetSize / 1024) : -1;
#elif defined(__APPLE__)
  // [mac] no /proc: answered from Mach. Without this MemAvailable read as 0 and the window loop admitted one window
  // at a time (admission high-water 1 of 64). MemAvailable = free + inactive pages (psutil's "available" on macOS);
  // VmRSS = resident size; VmHWM = getrusage's ru_maxrss, in bytes on Darwin.
  (void)path;
  const std::string k(key);
  if (k == "MemAvailable:")
  {
    vm_statistics64_data_t vs {};
    mach_msg_type_number_t n = HOST_VM_INFO64_COUNT;
    vm_size_t page = 0;
    if (host_page_size(mach_host_self(), &page) != KERN_SUCCESS
        || host_statistics64(mach_host_self(), HOST_VM_INFO64, (host_info64_t)&vs, &n) != KERN_SUCCESS)
    {
      return -1;
    }
    return (long long)(((unsigned long long)vs.free_count + vs.inactive_count) * page / 1024);
  }
  if (k == "VmRSS:")
  {
    mach_task_basic_info_data_t ti {};
    mach_msg_type_number_t n = MACH_TASK_BASIC_INFO_COUNT;
    return task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&ti, &n) == KERN_SUCCESS ? (long long)(ti.resident_size / 1024) : -1;
  }
  if (k == "VmHWM:")
  {
    struct rusage ru;
    return getrusage(RUSAGE_SELF, &ru) == 0 ? (long long)(ru.ru_maxrss / 1024) : -1;
  }
  return -1;
#else
  std::ifstream f(path);
  std::string k;
  while (f >> k)
  {
    if (k == key)
    {
      long long v = -1;
      f >> v;
      return v;
    }
    f.ignore(1 << 20, '\n');
  }
  return -1;
#endif
}
/**
  @brief Resident set size of the process.

  @return RSS in MB (0 if unreadable)
*/
static long rssMb()
{ return (long)(std::max(0LL, procKb("/proc/self/status", "VmRSS:")) / 1024); }
/**
  @brief Add one closed phase to its statistics (thread-safe).

  @param[in] name Phase name
  @param[in] w0 phaseClock() at the start
  @param[in] c0 cpuSeconds() at the start
  @param[in] nested Whether the phase was opened inside another
*/
static void phaseAdd(const std::string& name, double w0, double c0, bool nested)
{
  // [pipe] the pipelined tile loop times the tile read (driver) and the tile write (writer task) concurrently
  static std::mutex m;
  std::lock_guard<std::mutex> g(m);
  auto& st = phaseStats()[name];
  if (st.n++ == 0) { phaseOrder().push_back(name); }
  st.wall += phaseClock() - w0;
  st.cpu += cpuSeconds() - c0;
  if (nested) { st.nested = true; }
  st.rss_end_mb = rssMb();
}
/// @brief RAII timer of one named phase: records wall and CPU time from construction to destruction via phaseAdd().
struct Phase
{
  std::string name;
  double w0, c0;
  int d;
  /**
    @brief Open a phase.

    @param[in] n Phase name
  */
  explicit Phase(std::string n): name(std::move(n)), w0(phaseClock()), c0(cpuSeconds()), d(phaseDepth()++)
  {
  }
  /// @brief Close the phase and add it to its statistics.
  ~Phase()
  {
    --phaseDepth();
    phaseAdd(name, w0, c0, d > 0);
  }
};
/**
  @brief Log the per-phase breakdown (wall, share of the total, CPU, parallelism, RSS) and the measured total.

  @param[in] total_wall Wall time of the whole run, in seconds
*/
static void reportPhases(double total_wall)
{
  if (phaseOrder().empty()) { return; }
  OPENMS_LOG_INFO << "\n[perf] phase breakdown (wall s | % of total | CPU s | par x | RSS MB at end)\n";
  double acc = 0;
  for (const auto& n : phaseOrder())
  {
    const auto& p = phaseStats()[n];
    if (! p.nested)
    {
      acc += p.wall; // enclosed rows are already inside their parent's wall
    }
    char b[220];
    const char* tag = p.nested ? (p.n > 1 ? "  (summed, enclosed)" : "  (enclosed)") : (p.n > 1 ? "  (summed)" : "");
    snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  %9.1f  %5.1fx  %8ld%s\n", n.c_str(), p.wall,
             total_wall > 0 ? 100.0 * p.wall / total_wall : 0.0, p.cpu, p.wall > 0 ? p.cpu / p.wall : 0.0, p.rss_end_mb, tag);
    OPENMS_LOG_INFO << b;
  }
  char b[160];
  snprintf(b, sizeof b, "[perf] %-22s %8.1f  %5.1f%%  (unattributed %.1f s)\n", "TOTAL(measured)", acc,
           total_wall > 0 ? 100.0 * acc / total_wall : 0.0, total_wall - acc);
  OPENMS_LOG_INFO << b;
}

// The option mirrors below are assigned from their options at the start of main_; the initial values are the defaults.
double cfg_im_weight_sigma = 0.005; ///< [recall] assembly:im_weight_sigma
int cfg_mono_check = 2;             ///< [monocheck] charge:mono_profile_check: 0 false, 2 v2
double cfg_mono_pos_min_ratio = 0.2, cfg_mono_pos_im_box = 0.0,
       cfg_mono_pos_extend = 0.6;                ///< [monocheck-v2] charge:mono_position_* (im_box 0 = auto)
bool cfg_mono_pos_iterate = true;                ///< [monocheck-v2] charge:mono_position_iterate
String cfg_mono_check_stamp;                     ///< [monocheck] the spx:monocheck run stamp, filled by inferPrecursors_
bool cfg_short_conc = true;                      ///< [short] gate:short_traces == "concentric": span test for traces too short for Pearson
double cfg_short_center_tol = 0.5;               ///< [short] gate:short_center_tol, fraction of the precursor half-span
int cfg_short_min_pts = 1;                       ///< [short] gate:short_min_points
bool cfg_coloc = true;                           ///< [coloc] gate:colocation == "ranges": each apex inside the other span, span overlap
double cfg_coloc_overlap = 0.3;                  ///< [coloc] gate:colocation_overlap, fraction of the shorter span
bool cfg_shrink = true;                          ///< [shrink] gate:shrinkage == "fisher": Fisher-z posterior mean instead of a hard point cut
double cfg_shrink_r0 = 0.6, cfg_shrink_sd = 0.3; ///< [shrink] prior N(atanh r0, sd^2) on the Fisher z of the true correlation
int cfg_shrink_min_pts = 1;                      ///< [shrink] gate:shrink_min_points, overlapping frames needed to compute r at all
bool cfg_adaptive = true;                        ///< [adaptive] gate:adaptive == "faint": the relaxed modes only for faint precursors
bool cfg_adaptive_auto = false;    ///< [adaptive] gate:adaptive == "auto": faint mode only when the run's faint fraction >= cfg_auto_faint_frac
double cfg_auto_faint_frac = 0.45; ///< [adaptive] gate:auto_faint_fraction
int cfg_faint_frames = 4;          ///< [adaptive] gate:faint_max_frames
double cfg_faint_drt = 6.0, cfg_faint_cp = 1.0; ///< [adaptive] gate:faint_delta_rt, gate:faint_corr_power
double cfg_corr_power = 2.0;                    ///< assembly:corr_power: emitted fragment intensity *= corr^corr_power
double cfg_im_slice_merge = 0.015;              ///< [merge] trace:im_slice_merge, 1/K0 link (group span <= 2x); negative = off
bool cfg_merge_same_bin = true;                 ///< [merge] assembly:merge_same_bin
bool cfg_spectrum_meta = false;                 ///< [diag] diag:spectrum_meta: per-spectrum precursor support userParams
bool cfg_twin_ext_mark
  = false; ///< [twin-rt] mark spectra whose mono v2's extension placed (spx_mono_ext, never written); set only when the cross-frame fold is active
double cfg_foreign_w = 0.0;     ///< [foreign] assembly:foreign_precursor_weight (1 = off, 0 removes)
double cfg_foreign_ppm = 10.0;  ///< [foreign] assembly:foreign_precursor_ppm
double cfg_foreign_dim = 0.015; ///< [foreign] assembly:foreign_precursor_dim
int cfg_foreign_min_iso = 2;    ///< [foreign] assembly:foreign_precursor_min_iso
double cfg_fuse_ppm = 25.0;     ///< [fuse] assembly:fragment_fuse_ppm (0 = off)

/**
  @brief [merge] [0] picked MS2 peaks merged away in toTof (one frame, one bin, 1/K0 within trace:im_slice_merge);
  [1] frames that had such a merge; [2] fragments merged away at assembly (one bin in one spectrum)

  @return The process-wide counters
*/
inline std::array<std::atomic<long long>, 3>& mergeCount()
{
  static std::array<std::atomic<long long>, 3> c {};
  return c;
}
/**
  @brief [foreign] per assembly (one precursor in one window): [6] assemblies checked (with fragments), [0] of them with a
  foreign candidate, [1] foreign envelopes matched, [2] assemblies with a fragment attenuated, [3] fragments
  attenuated, [4] candidates skipped as this precursor's own ion, [5] assemblies of [2] with at least
  assembly:min_fragments before the removal (weight 0) and fewer after it and the same-bin merge

  @return The process-wide counters
*/
inline std::array<std::atomic<long long>, 7>& foreignCount()
{
  static std::array<std::atomic<long long>, 7> c {};
  return c;
}
/**
  @brief [fuse] assembly:fragment_fuse_ppm, per assembly (one precursor in one window) with at least one fragment: [0]
  assemblies, [1] fragments in, [2] peaks out, [3] groups of >= 2 fused, [4] fragments in those groups, [5] assemblies
  with at least assembly:min_fragments before the fusion and fewer after it

  @return The process-wide counters
*/
inline std::array<std::atomic<long long>, 6>& fuseCount()
{
  static std::array<std::atomic<long long>, 6> c {};
  return c;
}
/**
  @brief " [RSS cur, peak hwm]" from /proc/self/status; "" where unsupported.

  [mem]

  @return The stamp text, or an empty string
*/
static std::string rss()
{
  const long long cur = procKb("/proc/self/status", "VmRSS:"), peak = procKb("/proc/self/status", "VmHWM:");
  if (cur < 0 && peak < 0) { return ""; }
  auto kb = [](long long v) { return v < 0 ? std::string("?") : std::to_string(v) + " kB"; };
  return " [RSS " + kb(cur) + ", peak " + kb(peak) + "]";
}

/**
  @brief The headroom under the cgroup limits of this process (src/CgroupMem.h).

  Limit minus usage net of reclaimable page
  cache, the least over the leaf and its ancestors. Unlimited or unreadable: SIZE_MAX.

  @return Bytes of headroom, SIZE_MAX when unlimited or unreadable
*/
static Size cgroupHeadroom()
{ return DIAspeXtract::Cgroup::headroom(); }
/**
  @brief [dyn-mem] MemAvailable in bytes, capped by cgroupHeadroom().

  [dyn-mem] MemAvailable in bytes: the kernel's estimate of what can be allocated without swapping on a
  shared node. Our own RSS is NOT added back: after perf:malloc_trim it is the windows in flight, memory
  already spoken for. 0 if unreadable (the caller then admits one window at a time).

  @return Bytes available to this process
*/
static Size availableBytes()
{ return std::min((Size)std::max(0LL, procKb("/proc/meminfo", "MemAvailable:")) * (Size)1024, cgroupHeadroom()); }

/**
  @brief The global RT axis: every distinct frame retention time (MS1 and MS2), sorted.

  No RT is stored per
  point: each TraceStore maps its frames to indices on this axis (rt_index). Built once after loading,
  read-only thereafter.

  @return The process-wide axis
*/
inline vector<double>& rtAxis()
{
  static vector<double> v;
  return v;
}

/**
  @brief Index of the element of sorted @p a nearest to @p x (a tie keeps the later one).

  @param[in] a Sorted, non-empty values
  @param[in] x The value to look up
  @return Index of the nearest element
*/
inline uint32_t nearestIndex(const vector<double>& a, double x)
{
  auto it = lower_bound(a.begin(), a.end(), x);
  if (it == a.end()) { return (uint32_t)(a.size() - 1); }
  if (it != a.begin() && (*it - x) > (x - *(it - 1))) { --it; }
  return (uint32_t)(it - a.begin());
}

/**
  @brief Frame index of an RT.

  Used only while BUILDING traces: the value comes from a frame, so the
  match is exact, and the nearest-neighbour fallback only guards a last-ulp difference.

  @param[in] rt Retention time
  @return Index on rtAxis() (0 if the axis is empty)
*/
inline uint32_t rtIndex(double rt)
{ return rtAxis().empty() ? 0 : nearestIndex(rtAxis(), rt); }

/// @brief One mass trace: cached m/z, ion mobility, flight-time bin and a profile span in its TraceStore's arena (40 bytes).
struct Trace
{
  double mz = 0.0;
  double im = 0.0; ///< centroid ion mobility (1/K0)
  /// The flight-time bin this trace was measured on, and the calibration factor of its apex
  /// frame. The pair is the AUTHORITATIVE m/z: `mz` below is a cached convenience for the gates
  /// that still compare in m/z space, and the value actually written out is recomputed from these
  /// two at export. 0 means the trace came from the OpenMS detector, which has no bin.
  uint32_t tof = 0;
  uint32_t bframe = 0; ///< store-local frame whose factor calibrates `tof` (see bOf())
  /// The profile is a span of the owning TraceStore's arena, [off, off+len), one slot per frame
  /// from `frame0` (store-local: a window's frames are not contiguous on the global RT axis).
  /// ZERO MEANS MISSING, valid because every observed intensity is > noise_threshold >= 0;
  /// consumers walk the real points in frame order.
  uint32_t frame0 = 0; ///< first frame of the span, store-local
  uint32_t off = 0;    ///< offset of the span in the store's arena
  uint16_t len = 0;    ///< frames spanned (a gradient is ~1,400; asserted < 65536)
  uint16_t npts = 0;   ///< REAL points in the span
  uint16_t apex = 0;   ///< position of the apex within the span
  /// [soa-rt] The RT is a frame index. At the shipped defaults every producer reports one of the store's
  /// frame times (the integer detector its apex frame, ElutionPeakDetection its smoothed-maximum point);
  /// an off-grid RT (EPD skipped) is snapped and counted, see rtOffGrid(). Read it with rtOf(); never
  /// derive it from `apex` (the RAW maximum, a different frame). SIGNED: trimToSpan moves `frame0` by
  /// +lo and this by -lo, so the reported RT survives any trim.
  int16_t rt_at = 0; ///< frame whose time is this trace's RT, relative to frame0 (see rtOf())
  /**
    @brief Number of real (non-zero) points in the span.

    @return npts
  */
  Size np() const
  { return npts; }
  /**
    @brief Number of frames the span covers.

    @return len
  */
  Size span() const
  { return len; }
  /// @brief Release the profile: the span becomes empty, the scalars stay.
  void freeProfile()
  {
    len = 0;
    npts = 0;
  } ///< the arena is the store's; nothing to free
};

/**
  @brief The frame table and the arenas that a set of traces point into: one per window, one for the
  MS1 map.

  `rt_index`/`b` are per FRAME (store-local index); `inten` holds every trace's span
  back to back (zero = missing); `bins` holds the per-frame flight-time bin for the integer
  path only and is freed once valley splitting has assigned each child its apex bin.
*/
struct TraceStore
{
  vector<uint32_t> rt_index; ///< frame -> global RT axis
  vector<double> b;          ///< frame -> calibration factor (empty on the OpenMS path)
  vector<double> frame_rt;   ///< frame -> RT value, for toTrace()'s frame lookup
  vector<float> inten;       ///< arena of spans
  /// [bin-delta] A point's flight-time bin as the SIGNED offset from its trace's own bin
  /// (`Trace::tof`): a trace spans a few bins, int16 leaves ~30,000 of headroom, and the writer
  /// throws rather than truncate. Reconstruct with `t.tof + bins[t.off + k]`; a binless trace
  /// (tof 0, zero-filled slots) reconstructs to 0.
  vector<int16_t> bins; ///< arena of per-frame bin OFFSETS (integer path, until EPD)
  /**
    @brief Number of frames in the store's frame table.

    @return rt_index.size()
  */
  Size frames() const
  { return rt_index.size(); }
  /**
    @brief Set the frame table and fill frame_rt from rtAxis().

    @param[in] ri Global RT-axis index per frame
    @param[in] bb Calibration factor per frame, or nullptr for none
  */
  void setFrames(const vector<uint32_t>& ri, const vector<double>* bb)
  {
    rt_index = ri;
    if (bb) { b = *bb; }
    else
    {
      b.clear();
    }
    frame_rt.resize(ri.size());
    for (Size f = 0; f < ri.size(); ++f)
    {
      frame_rt[f] = rtAxis()[ri[f]];
    }
  }
  /**
    @brief Store-local frame of an RT value (exact by construction; nearest guards a last-ulp slip).

    @param[in] rt Retention time
    @return Store-local frame index
  */
  uint32_t frameOf(double rt) const
  { return nearestIndex(frame_rt, rt); }
};
static_assert(sizeof(Trace) == 40,
              "Trace must stay 40 bytes: it is the largest structure in the tool "
              "(2.1e9 records on a 2-h file), and a silent 8-byte regrowth costs ~16 GiB");

/**
  @brief [serial] Appends every source's arenas to `dst` in order and its records, rebased (off += its base; frame0 and
  bframe += its frame shift), to `out`, with the copies run as bounded tasks at prefix-summed offsets.

  Needs an
  empty `dst` and `out` and every source's bins index-parallel to its inten (or absent), as every caller has. The
  sources are emptied.

  @param[in,out] dst Destination store (must be empty)
  @param[in,out] out Destination records (must be empty)
  @param[in,out] src Pairs of source store and records; emptied
  @param[in] fshift Frame shift per source
  @throws OpenMS::Exception::Precondition if @p dst or @p out is not empty or a source's bins are not index-parallel to its inten
  @throws OpenMS::Exception::OutOfRange if the combined arena reaches 2^32 - 1 floats
*/
inline void absorbAll(TraceStore& dst, vector<Trace>& out, const vector<std::pair<TraceStore*, vector<Trace>*>>& src, const vector<uint32_t>& fshift)
{
  const Size ns = src.size();
  vector<Size> ib(ns + 1, 0), tb(ns + 1, 0); // arena / record prefix sums
  bool par = dst.inten.empty() && dst.bins.empty() && out.empty(), any_bins = false;
  for (Size s = 0; s < ns; ++s)
  {
    const TraceStore& o = *src[s].first;
    if (! o.bins.empty())
    {
      any_bins = true;
      par &= (o.bins.size() == o.inten.size());
    }
    ib[s + 1] = ib[s] + o.inten.size();
    tb[s + 1] = tb[s] + src[s].second->size();
  }
  if (! par)
  {
    throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "absorbAll_ needs an empty destination and index-parallel bins");
  }
  // Trace::off is 32-bit: refuse an arena that would reach 2^32 - 1 floats
  if (ib[ns] >= (Size)std::numeric_limits<uint32_t>::max()) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  // once any source carries bins they are index-parallel to inten, zero where a source has none
  dst.inten.resize(ib[ns]);
  if (any_bins) { dst.bins.resize(ib[ns], 0); }
  out.resize(tb[ns]);
  // call f(s, a, e) for the part [a, e) of source s that falls in slice c of nt of the concatenation `pre`
  auto slices = [ns](const vector<Size>& pre, long c, long nt, auto&& f) {
    const Size n_total = pre[ns], lo = (Size)c * n_total / (Size)nt, hi = (Size)(c + 1) * n_total / (Size)nt;
    Size s = (Size)(std::upper_bound(pre.begin(), pre.end(), lo) - pre.begin()) - 1;
    for (; s < ns && pre[s] < hi; ++s)
    {
      const Size a = std::max(lo, pre[s]) - pre[s], e = std::min(hi, pre[s + 1]) - pre[s];
      if (a < e) { f(s, a, e); }
    }
  };
  const long nt = parTasks(ib[ns] + 8 * tb[ns], 1 << 20); // a record is ~8 floats of copying
  forTasks(nt, [&](long c) {
    slices(ib, c, nt, [&](Size s, Size a, Size e) {
      const TraceStore& o = *src[s].first;
      std::copy(o.inten.begin() + a, o.inten.begin() + e, dst.inten.begin() + ib[s] + a);
      if (! o.bins.empty()) { std::copy(o.bins.begin() + a, o.bins.begin() + e, dst.bins.begin() + ib[s] + a); }
    });
    slices(tb, c, nt, [&](Size s, Size a, Size e) {
      const uint32_t base = (uint32_t)ib[s], fs = fshift[s];
      const vector<Trace>& r = *src[s].second;
      for (Size q = a; q < e; ++q)
      {
        Trace t = r[q];
        t.off += base;
        t.frame0 += fs;
        t.bframe += fs;
        out[tb[s] + q] = t;
      }
    });
  });
  for (Size s = 0; s < ns; ++s)
  {
    vector<float>().swap(src[s].first->inten);
    vector<int16_t>().swap(src[s].first->bins);
    vector<Trace>().swap(*src[s].second);
  }
}

/**
  @brief Span accessors take the owning store explicitly (no per-record back-pointer); bounds asserted in debug builds.

  @param[in] s Owning store
  @param[in] t The trace
  @param[in] k Position in the span
  @return The intensity at @p k (0 = missing)
*/
inline float xv(const TraceStore& s, const Trace& t, Size k)
{
  assert(t.off + k < s.inten.size());
  return s.inten[t.off + k];
}
/**
  @brief Whether span position @p k holds an observed point.

  @param[in] s Owning store
  @param[in] t The trace
  @param[in] k Position in the span
  @return true if the intensity is positive
*/
inline bool real(const TraceStore& s, const Trace& t, Size k)
{ return xv(s, t, k) > 0.0f; }
/**
  @brief Retention time of span position @p k.

  @param[in] s Owning store
  @param[in] t The trace
  @param[in] k Position in the span
  @return The frame's retention time
*/
inline double rtAtSpan(const TraceStore& s, const Trace& t, Size k)
{
  assert(t.frame0 + k < s.rt_index.size());
  return rtAxis()[s.rt_index[t.frame0 + k]];
}
/**
  @brief [soa-rt] The trace's reported RT.

  `frame_rt` is one double per frame of the owning store
  (~38 KB for a window), so this is a single load from a table that stays hot, against the
  8 bytes per record it replaces. Identical by construction: frame_rt[f] IS rtAxis()[ri[f]],
  the same array element the producers read.

  @param[in] s Owning store
  @param[in] t The trace
  @return The trace's retention time
*/
inline double rtOf(const TraceStore& s, const Trace& t)
{
  const int64_t f = (int64_t)t.frame0 + t.rt_at;
  assert(f >= 0 && (Size)f < s.frame_rt.size());
  return s.frame_rt[(Size)f];
}
/**
  @brief Traces whose RT is NOT one of the store's frame times: only reachable when
  ElutionPeakDetection is skipped (trace:ms{1,2}_split_valleys = 0, or a single trace), where
  MassTraceDetection reports an intensity-weighted mean instead.

  Counted, not silent.

  @return The process-wide counter
*/
inline std::atomic<Size>& rtOffGrid()
{
  static std::atomic<Size> n {0};
  return n;
}
/**
  @brief Record `rt` as a frame offset.

  Snaps to the nearest frame if the value is off-grid (counted).

  @param[in] s Owning store
  @param[in,out] t The trace; its rt_at is set
  @param[in] rt Retention time
  @throws OpenMS::Exception::OutOfRange if the frame offset does not fit into 16 bits
*/
inline void setRtFrame(const TraceStore& s, Trace& t, double rt)
{
  if (s.frame_rt.empty())
  {
    t.rt_at = 0;
    return;
  }
  const uint32_t f = s.frameOf(rt);
  if (s.frame_rt[f] != rt) { rtOffGrid().fetch_add(1, std::memory_order_relaxed); }
  const int64_t d = (int64_t)f - (int64_t)t.frame0;
  if (d < -32768 || d > 32767) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  t.rt_at = (int16_t)d;
}
/**
  @brief The calibration factor of the trace's own apex/measurement frame.

  Empty table (the OpenMS
  detector has no per-frame factors) yields 0.0, and exportMz_ returns the cached m/z for those
  traces because their `tof` is 0.

  @param[in] s Owning store
  @param[in] t The trace
  @return The calibration factor, or 0.0
*/
inline double bOf(const TraceStore& s, const Trace& t)
{ return s.b.empty() ? 0.0 : s.b[t.bframe]; }
/**
  @brief The trace's apex intensity: the arena value at `apex` (makeSpan's maximum; trimming moves `apex`
  with the data).

  0 for a released profile (len == 0), which is never scored.

  @param[in] s Owning store
  @param[in] t The trace
  @return The apex intensity
*/
inline double intensityOf(const TraceStore& s, const Trace& t)
{ return t.len ? (double)s.inten[t.off + t.apex] : 0.0; }

/**
  @brief [fwhm] FWHM (seconds) from the trace's own XIC: apex, then outwards to half maximum with linear interpolation;
  a side that never reaches half maximum ends at the first/last point.

  -1 below 3 points.

  @param[in] s Owning store
  @param[in] tr The trace
  @return FWHM in seconds, or -1
*/
inline double xicFwhm(const TraceStore& s, const Trace& tr)
{
  if (tr.np() < 3) { return -1.0; }
  vector<double> prt, pv; // the real points, in frame order
  for (Size k = 0; k < tr.span(); ++k)
  {
    if (real(s, tr, k))
    {
      prt.push_back(rtAtSpan(s, tr, k));
      pv.push_back((double)xv(s, tr, k));
    }
  }
  Size ap = 0;
  for (Size i = 1; i < pv.size(); ++i)
  {
    if (pv[i] > pv[ap]) { ap = i; }
  }
  const double half = pv[ap] * 0.5;
  double lo = prt[0], hi = prt[pv.size() - 1];
  for (Size i = ap; i > 0; --i)
  {
    if (pv[i - 1] <= half)
    {
      const double d = pv[i] - pv[i - 1];
      const double f = d > 0 ? (pv[i] - half) / d : 0.0;
      lo = prt[i] - f * (prt[i] - prt[i - 1]);
      break;
    }
  }
  for (Size i = ap; i + 1 < pv.size(); ++i)
  {
    if (pv[i + 1] <= half)
    {
      const double d = pv[i] - pv[i + 1];
      const double f = d > 0 ? (pv[i] - half) / d : 0.0;
      hi = prt[i] + f * (prt[i + 1] - prt[i]);
      break;
    }
  }
  return hi - lo;
}

/**
  @brief Trim a trace to at most `cap` seconds around its apex, AFTER detection and valley splitting, so
  peak ownership is unchanged (rationale: trace:max_span_sec).

  Only the span fields change and no
  bytes are copied; the arena keeps the trimmed bytes until a compaction pass.

  @param[in] s Owning store
  @param[in,out] t The trace to trim
  @param[in] cap Maximum span in seconds (<= 0: no trim)
  @throws OpenMS::Exception::OutOfRange if the RT frame offset does not fit into 16 bits
*/
inline void trimToSpan(const TraceStore& s, Trace& t, double cap)
{
  if (cap <= 0.0 || t.span() == 0) { return; }
  const Size a = t.apex;
  Size lo = a, hi = a;
  // grow symmetrically outward from the apex while the span fits, preferring the earlier frame on
  // a tie so the result does not depend on which side is tested first
  for (;;)
  {
    const bool can_lo = lo > 0, can_hi = hi + 1 < t.span();
    if (! can_lo && ! can_hi) { break; }
    const double d_lo = can_lo ? rtAtSpan(s, t, hi) - rtAtSpan(s, t, lo - 1) : 1e30;
    const double d_hi = can_hi ? rtAtSpan(s, t, hi + 1) - rtAtSpan(s, t, lo) : 1e30;
    if (d_lo <= d_hi)
    {
      if (! can_lo || d_lo > cap)
      {
        if (! can_hi || d_hi > cap) { break; }
        ++hi;
      }
      else
      {
        --lo;
      }
    }
    else
    {
      if (! can_hi || d_hi > cap)
      {
        if (! can_lo || d_lo > cap) { break; }
        --lo;
      }
      else
      {
        ++hi;
      }
    }
  }
  // Zero-trim the ends. makeSpan starts and ends a span on a real point by construction, but
  // growing outward from the apex can stop on a gap, and a leading or trailing zero is pure
  // padding: it costs 4 B, it is skipped by every consumer, and it makes `span()` overstate the
  // trace. Interior zeros are the gaps and MUST stay -- they are the presence encoding.
  while (lo < hi && xv(s, t, lo) == 0.0f)
  {
    ++lo;
  }
  while (hi > lo && xv(s, t, hi) == 0.0f)
  {
    --hi;
  }
  if (lo == 0 && hi + 1 == t.span())
  {
    return; // already within the cap, ends real
  }
  uint16_t np = 0;
  for (Size k = lo; k <= hi; ++k)
  {
    if (xv(s, t, k) > 0.0f) { ++np; }
  }
  // [soa-rt] rt_at names EPD's smoothed apex, not `apex`, and may lie outside the kept span. It is
  // not clamped: frame0 and rt_at move by opposite amounts, so the reported RT is unchanged.
  const long rr = (long)t.rt_at - (long)lo;
  if (rr < -32768L || rr > 32767L) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  t.frame0 += (uint32_t)lo;
  t.off += (uint32_t)lo;
  t.len = (uint16_t)(hi - lo + 1);
  t.apex = (uint16_t)(a - lo);
  t.npts = np;
  t.rt_at = (int16_t)rr;
}

/**
  @brief Drop the spans nobody references and compact the arena IN PLACE in OFFSET order (ms1_traces is
  m/z-sorted, so container order would overwrite spans not yet copied).

  Spans must be disjoint --
  validated before a byte moves. A freed trace keeps its scalars; returns the profiles released.

  @param[in,out] tr The traces; unneeded ones are freed and the rest rebased
  @param[in,out] store Their store; its arena is compacted
  @param[in] needed Per trace, whether its profile is kept
  @return Number of profiles released
  @throws OpenMS::Exception::Precondition if kept spans overlap or exceed the arena
*/
inline Size compactUnreferenced(vector<Trace>& tr, TraceStore& store, const vector<bool>& needed)
{
  vector<uint32_t> ids;
  Size freed = 0;
  for (Size i = 0; i < tr.size(); ++i)
  {
    if (! needed[i])
    {
      if (tr[i].np() > 0) { ++freed; }
      tr[i].freeProfile();
      continue;
    }
    if (tr[i].len > 0)
    {
      ids.push_back((uint32_t)i); // a zero-length span owns no bytes
    }
  }
  // [outer] valid spans have distinct offsets (equal ones overlap and throw below either way), so with the index as a
  // last key the order is total and the parallel sort returns std::sort's
  parallelSortTotal(ids, [&](uint32_t a, uint32_t b) { return tr[a].off != tr[b].off ? tr[a].off < tr[b].off : a < b; });
  uint64_t prev_end = 0;
  for (uint32_t i : ids)
  {
    const uint64_t beg = tr[i].off, end = beg + (uint64_t)tr[i].len;
    if (end > store.inten.size() || beg < prev_end)
    {
      throw OpenMS::Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MS1 spans overlap or exceed the arena");
    }
    prev_end = end;
  }
  Size w = 0;
  for (uint32_t i : ids)
  {
    Trace& t = tr[i];
    if (w != t.off) { std::memmove(store.inten.data() + w, store.inten.data() + t.off, (Size)t.len * sizeof(float)); }
    t.off = (uint32_t)w;
    w += t.len;
  }
  store.inten.resize(w);
  store.inten.shrink_to_fit();
  return freed;
}

/**
  @brief Build a span from (store-local frame, intensity) points into `store`, returning the trace
  with `st` UNSET: the caller sets it once the store is immutable.

  Points may arrive unsorted.
  The span is appended to the arena, or placed at @p off when given (that part of the arena
  already sized and zero).

  @param[in,out] pts (frame, intensity) points; sorted in place
  @param[in,out] store Receives the span
  @param[in] off Arena offset, or SIZE_MAX to append
  @return The trace with its span fields and apex set
  @throws OpenMS::Exception::OutOfRange if the span covers 65536 frames or more
  @throws OpenMS::Exception::InvalidValue if the arena would reach 2^32 points
*/
inline Trace makeSpan(vector<pair<uint32_t, float>>& pts, TraceStore& store, Size off = SIZE_MAX)
{
  Trace t;
  if (pts.empty()) { return t; }
  sort(pts.begin(), pts.end());
  const uint32_t f0 = pts.front().first, f1 = pts.back().first;
  const Size trace_len = (Size)f1 - f0 + 1;
  if (trace_len >= 65536) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  if (off == SIZE_MAX) { off = store.inten.size(); }
  if (off + trace_len > (Size)std::numeric_limits<uint32_t>::max()) // Trace::off is 32-bit: a wrapped offset would overwrite another trace
  {
    throw OpenMS::Exception::InvalidValue(
      __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "one trace store holds 2^32 or more profile "
      "points; for an isolation window, run in smaller tiles (-tile:cells_per_tile 1, a smaller -tile:rt_sec); the MS1 store "
      "is built before tiling, so a run this large has to be split by retention time (bench/slice_d.py)",
      String(off + trace_len));
  }
  if (off + trace_len > store.inten.size()) { store.inten.resize(off + trace_len, 0.0f); }
  t.frame0 = f0;
  t.len = (uint16_t)trace_len;
  t.off = (uint32_t)off;
  float best = -1.0f;
  uint16_t np = 0;
  for (const auto& q : pts)
  {
    const Size k = q.first - f0;
    if (store.inten[t.off + k] == 0.0f)
    {
      ++np; // one point per frame; a repeat overwrites
    }
    store.inten[t.off + k] = q.second;
    if (q.second > best)
    {
      best = q.second;
      t.apex = (uint16_t)k;
    }
  }
  t.npts = np;
  return t;
}

/// @brief A precursor hypothesis after isotope/charge inference.
struct PrecursorCand
{
  double mono_mz = 0.0;
  int charge = 0; ///< 0 == unknown
  double rt = 0.0;
  double im = 0.0;
  Size trace_idx = 0;    ///< index into ms1 traces
  bool mono_ext = false; ///< [twin-rt] the monoisotope was placed by monocheck v2's extension (B)
  int n_isotopes = 0;    ///< isotope peaks behind the mono call (0/1 = weak, mono ±1.00335 uncertain)
};

/**
  @brief [monocheck-v2] What charge:mono_profile_check v2 needs of a count walk, one per precursor of inferPrecursors_
  (z = 0: nothing to test): the winning charge, its leftward partners' m/z (nearest first) and the monoisotope
  before v2.
*/
struct MonoV2In
{
  int z = 0, nl = 0;
  double mono = 0.0;
  std::array<double, 5> left_mz {};
};

/**
  @brief Give the spectrum's ion-mobility float data array the name MassTraceDetection looks for.

  MassTraceDetection locates the per-peak IM array by the exact name
  Constants::UserParam::ION_MOBILITY ("Ion Mobility"). Converters name it via CV term and
  PeakPickerIM names it "Ion Mobility Centroid"; rename whichever exists so tracing is IM-aware.

  @param[in,out] s The spectrum
*/
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

/**
  @brief [outer] The span length makeSpan gives @p mt: its points' frame range, 0 without points.

  Throws as makeSpan does.

  @param[in] mt The mass trace
  @param[in] store Store whose frame table maps its RTs
  @return Span length in frames
  @throws OpenMS::Exception::OutOfRange if the span covers 65536 frames or more
*/
inline Size spanLength(const MassTrace& mt, const TraceStore& store)
{
  const Size np = mt.getSize();
  if (np == 0) { return 0; }
  uint32_t f0 = std::numeric_limits<uint32_t>::max(), f1 = 0;
  for (Size i = 0; i < np; ++i)
  {
    const uint32_t f = store.frameOf(mt[i].getRT());
    f0 = std::min(f0, f);
    f1 = std::max(f1, f);
  }
  const Size trace_len = (Size)f1 - f0 + 1;
  if (trace_len >= 65536) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
  return trace_len;
}

/**
  @brief A MassTrace as a Trace, its span appended to `store` (@p off SIZE_MAX) or placed at arena offset @p off (already
  sized and zero, see makeSpan).

  @p pts is the caller's scratch.

  @param[in] mt The mass trace
  @param[in,out] store Receives the span
  @param[in] off Arena offset, or SIZE_MAX to append
  @param[in,out] pts Scratch for the points
  @return The trace
  @throws OpenMS::Exception::OutOfRange, OpenMS::Exception::InvalidValue as makeSpan()
*/
inline Trace toTraceAt(const MassTrace& mt, TraceStore& store, Size off, vector<pair<uint32_t, float>>& pts)
{
  // Reported m/z of a mass trace: the m/z of its most intense peak (the centroid when it has none).
  const double mz = mt.getSize() > 0 ? mt[mt.findMaxByIntPeak(false)].getMZ() : mt.getCentroidMZ();
  // One peak per frame (MassTraceDetection takes at most one candidate per spectrum), so the
  // profile is a span over the store's frame table with zeros where frames were missed.
  const Size np = mt.getSize();
  pts.clear();
  pts.reserve(np);
  for (Size i = 0; i < np; ++i)
  {
    pts.emplace_back(store.frameOf(mt[i].getRT()), (float)mt[i].getIntensity());
  }
  Trace t = makeSpan(pts, store, off);
  t.mz = mz;
  t.im = mt.getCentroidIM();
  setRtFrame(store, t, mt.getCentroidRT()); // [soa-rt] a frame time, stored as an offset
  return t;
}
/**
  @brief A MassTrace as a Trace, its span appended to @p store.

  @param[in] mt The mass trace
  @param[in,out] store Receives the span
  @return The trace
*/
inline Trace toTrace(const MassTrace& mt, TraceStore& store)
{
  vector<pair<uint32_t, float>> pts;
  return toTraceAt(mt, store, SIZE_MAX, pts);
}


/// [compact] Quantised per-peak store, 10 B/peak, far finer than every downstream tolerance:
///   m/z : uint32 at 1e-5 Da     (0.007 ppm at m/z 1400)
///   1/K0: uint16 over [0.4,1.8] (2.1e-5)
///   intensity: float32, unchanged
constexpr double MZ_Q = 1e5; ///< m/z quantum = 1e-5 Da

/// The m/z axis is the instrument's flight-time bin: lossless, the same bin for the same ion in every
/// frame, and a ppm tolerance is a handful of integer bins. The per-frame factor `b` (digitizer
/// temperature) turns a bin into m/z at report time.
using TofIdx = uint32_t;

/**
  @brief The flight-time axis for one run: the calibration model plus each frame's temperature factor.

  Fail-closed: without a calibration there is no TOF axis and the caller must not pretend there
  is one.
*/
struct TofAxis
{
  diaspextract::TdfMzCalibration cal;
  vector<double> b_by_frame; ///< frame factor per FRAME ID (index 0 unused, as in the tdf)
  bool ok = false;
  double acq_lo = std::numeric_limits<double>::quiet_NaN(); ///< the acquisition's m/z range (GlobalMetadata)
  double acq_hi = std::numeric_limits<double>::quiet_NaN();
  long n_bins = 0;   ///< DigitizerNumSamples: every bin index is below it
  String bounds_why; ///< why the bounds are missing, if they are

  /**
    @brief Calibration factor of a frame; the reference factor (at T1_ref) for an unknown frame.

    @param[in] frame_id Vendor Frames.Id
    @return The factor
  */
  double factor(Size frame_id) const
  { return frame_id < b_by_frame.size() && ! std::isnan(b_by_frame[frame_id]) ? b_by_frame[frame_id] : cal.frameFactor(cal.T1_ref); }
  /**
    @brief m/z of a flight-time bin.

    @param[in] tof Flight-time bin
    @param[in] b Frame calibration factor
    @return The m/z
  */
  double mzOf(TofIdx tof, double b) const
  { return cal.tofToMz((double)tof, b); }
  /**
    @brief Nearest flight-time bin of an m/z, saturated at 4294967294.

    @param[in] mz The m/z
    @param[in] b Frame calibration factor
    @return The bin (0 if the model gives no positive time)
  */
  TofIdx tofOf(double mz, double b) const
  {
    const double t = cal.mzToTof(mz, b);
    if (! (t > 0.0)) { return 0u; }
    return t >= 4294967294.0 ? 4294967294u : (TofIdx)llround(t); // saturate: a wrap here would collapse a band partition silently
  }
};

/**
  @brief The run's flight-time axis.

  @return The process-wide axis
*/
inline TofAxis& tofAxis()
{
  static TofAxis a;
  return a;
}

/**
  @brief Build the flight-time axis from a calibration and the per-frame T1 (index = Frames.Id).

  @param[in] cal The calibration model
  @param[in] t1 T1 per frame id (NaN for an id the table lacks)
*/
inline void setTofAxis(const diaspextract::TdfMzCalibration& cal, const std::vector<double>& t1)
{
  TofAxis& ax = tofAxis();
  ax.cal = cal;
  // A frame without a positive factor gets the reference one. Index 0 is the "unknown frame"
  // sentinel (Frames.Id starts at 1; frameIdOf() returns 0 when it cannot parse) and stays 0.0. An id the table lacks
  // (t1 NaN: a gap) stays NaN, so frameCalibrated() counts it as unmapped and factor() gives the reference one.
  const double bref = ax.cal.frameFactor(ax.cal.T1_ref);
  ax.b_by_frame.assign(t1.size(), 0.0);
  for (Size i = 1; i < t1.size(); ++i)
  {
    if (std::isnan(t1[i]))
    {
      ax.b_by_frame[i] = t1[i];
      continue;
    }
    const double f = ax.cal.frameFactor(t1[i]);
    ax.b_by_frame[i] = f > 0.0 ? f : bref;
  }
  ax.ok = true;
}
/**
  @brief Load the axis from an analysis.tdf; false (axis unusable) when anything is missing.

  @param[in] tdf_path Path of the analysis.tdf
  @param[out] why Reason when false is returned
  @return true if the axis was set
*/
inline bool loadTofAxis(const std::string& tdf_path, String& why)
{
  diaspextract::TdfMzCalibration cal;
  std::vector<double> t1;
  if (! DIAspeXtract::loadTdfCalibration(tdf_path, cal, t1, why)) { return false; }
  setTofAxis(cal, t1);
  {
    std::string w2;
    DIAspeXtract::TdfAxisBounds bd;
    if (DIAspeXtract::loadTdfAxisBounds(tdf_path, bd, w2))
    {
      tofAxis().acq_lo = bd.mz_lo;
      tofAxis().acq_hi = bd.mz_hi;
      tofAxis().n_bins = bd.n_bins;
    }
    else
    {
      tofAxis().bounds_why = w2;
    }
  }
  return true;
}
constexpr double IM_LO = 0.4, IM_HI = 1.8;
constexpr double ISO_STEP = 1.0033548; ///< the isotope spacing (Da) of the isotope walks and envelopes
constexpr double IM_Q = 65535.0 / (IM_HI - IM_LO);

/**
  @brief Whether a frame has a per-frame calibration factor.

  The frame has a per-frame calibration factor: a parseable id the calibration table holds. Before the table is
  loaded (an openms-detector run) only the id is checked.

  @param[in] fid Vendor frame id
  @return true if the frame is calibrated
*/
inline bool frameCalibrated(uint32_t fid)
{ return fid != 0 && (! tofAxis().ok || (fid < tofAxis().b_by_frame.size() && ! std::isnan(tofAxis().b_by_frame[fid]))); }
/**
  @brief Counts of the compaction of picked peaks.

  Peaks that cannot be represented exactly are DROPPED and COUNTED: clamping 1/K0 would stack
  outliers on two artificial mobilities, and NaN into llround is UB.
*/
struct CompactStats
{
  Size no_im_array = 0, size_mismatch = 0, bad_mz = 0, bad_im = 0, kept = 0;
  Size unmapped_frame = 0; ///< no parseable "frame=" in the nativeID, or a frame the calibration table lacks
  /**
    @brief merge per-frame stats after a parallel batch [perf-load]

    @param[in] o The stats to add
    @return *this
  */
  CompactStats& operator+=(const CompactStats& o)
  {
    no_im_array += o.no_im_array;
    size_mismatch += o.size_mismatch;
    bad_mz += o.bad_mz;
    bad_im += o.bad_im;
    kept += o.kept;
    unmapped_frame += o.unmapped_frame;
    return *this;
  }
};

/**
  @brief The unsigned integer after `key` in a native ID; 0 when absent.

  @param[in] nid The native ID
  @param[in] key The key including its '='
  @return The value, or 0
*/
inline uint32_t nidField(const String& nid, const char* key)
{
  const Size i = nid.find(key);
  if (i == std::string::npos) { return 0; }
  uint32_t v = 0;
  Size k = i + strlen(key);
  for (; k < nid.size() && isdigit((unsigned char)nid[k]); ++k)
  {
    v = v * 10 + (uint32_t)(nid[k] - '0');
  }
  return v;
}
/**
  @brief The vendor Frames.Id ("frame=", written by the .d loader and the mzPeak reader; "mzpeak=" is an archive index).

  0 = unknown.

  @param[in] nid The native ID
  @return The frame id, or 0
*/
inline uint32_t frameIdOf(const String& nid)
{ return nidField(nid, "frame="); }
/**
  @brief "windowGroup=<WindowGroup>": which acquisition group this frame belongs to.

  A group holds each
  isolation m/z at most once, so (m/z, group) names one ion-mobility slice of the window; the
  scan range alone does not (two slices can share ScanNumBegin = 0 and differ in ScanNumEnd).
  0 when absent (mzPeak input), which keeps that path keyed by m/z alone.

  @param[in] nid The native ID
  @return The window group, or 0
*/
inline uint32_t windowGroupOf(const String& nid)
{ return nidField(nid, "windowGroup="); }

/**
  @brief [tile] One frame's picked peaks written into a caller-provided region sized by the raw count
  (a pick never emits more points).

  No per-frame allocation: vectors allocated on picking threads
  and freed on the main thread strand free-list in their arenas.

  @param[in] s The picked frame
  @param[in,out] st Counts of dropped and kept peaks
  @param[out] mzq m/z quanta
  @param[out] inten Intensities
  @param[out] imq Quantised 1/K0
  @param[in] room Size of the region
  @return Number of peaks written
  @throws OpenMS::Exception::InvalidValue if the picked frame is larger than @p room
*/
Size compactifyInto(const MSSpectrum& s, CompactStats& st, uint32_t* mzq, float* inten, uint16_t* imq, Size room)
{
  if (! s.containsIMData())
  {
    st.no_im_array += s.size();
    return 0;
  }
  const auto& im = s.getFloatDataArrays()[s.getIMData().first];
  const Size n = s.size();
  if (im.size() != n)
  {
    st.size_mismatch += n;
    return 0;
  } // no out-of-bounds read
  if (n > room)
  {
    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "picked frame larger than its raw region", std::to_string(n));
  }
  Size k = 0;
  for (Size i = 0; i < n; ++i)
  {
    const double mz = s[i].getMZ();
    const double imv = (double)im[i];
    if (! std::isfinite(mz) || mz <= 0.0 || mz * MZ_Q > 4.2e9)
    {
      ++st.bad_mz;
      continue;
    }
    if (! std::isfinite(imv) || imv < IM_LO || imv > IM_HI)
    {
      ++st.bad_im;
      continue;
    }
    mzq[k] = (uint32_t)llround(mz * MZ_Q);
    inten[k] = s[i].getIntensity();
    imq[k] = (uint16_t)llround((imv - IM_LO) * IM_Q);
    ++k;
    ++st.kept;
  }
  return k;
}

/**
  @brief [det] Order-insensitive 64-bit digest of a stage's numeric output, for localising
  run-to-run nondeterminism: sort the bit patterns, then FNV-1a.

  Logged under DIASPEXTRACT_DET=1.

  @param[in] bits The bit patterns (sorted here)
  @return The digest
*/
static uint64_t detDigest(vector<uint64_t> bits)
{
  std::sort(bits.begin(), bits.end());
  uint64_t h = 1469598103934665603ULL;
  for (uint64_t b : bits)
  {
    h ^= b;
    h *= 1099511628211ULL;
  }
  return h;
}
/**
  @brief Append the bit pattern of a double.

  @param[in,out] b The bit patterns
  @param[in] v The value
*/
static void pushBits(vector<uint64_t>& b, double v)
{
  uint64_t u;
  std::memcpy(&u, &v, 8);
  b.push_back(u);
}
/**
  @brief Whether DIASPEXTRACT_DET is set (determinism digests are logged).

  @return true if set
*/
static bool detOn()
{
  static const bool on = std::getenv("DIASPEXTRACT_DET") != nullptr;
  return on;
}
/**
  @brief Order-insensitive digest of traces: m/z, RT, 1/K0, apex intensity and flight-time bin of each.

  @param[in] st Owning store
  @param[in] ts The traces
  @return The digest
*/
static uint64_t traceDigest(const TraceStore& st, const vector<Trace>& ts)
{
  // im (the fragment sort key and scoring's search axis) and tof (the exported m/z) are hashed too.
  vector<uint64_t> b;
  b.reserve(ts.size() * 5);
  for (const auto& t : ts)
  {
    pushBits(b, t.mz);
    pushBits(b, rtOf(st, t));
    pushBits(b, t.im);
    pushBits(b, intensityOf(st, t));
    b.push_back((uint64_t)t.tof);
  }
  return detDigest(std::move(b));
}

/**
  @brief A window's peaks in index space, frame-major, appended by the loader one frame at a time and
  moved into the detector.

  Until toTof() runs, `tof` holds the 1e-5 Da m/z quantum (`mzq`), which
  the OpenMS-detector path reads as is.
*/
struct PeakSlab
{
  vector<uint32_t> frame_off; ///< size F+1: frame f owns [frame_off[f], frame_off[f+1])
  vector<uint32_t> rt_index;  ///< size F: the frame's index on the GLOBAL RT axis (indexRt)
  vector<TofIdx> tof;         ///< flight-time bin, ascending within a frame (m/z quanta while mzq)
  vector<double> b;           ///< size F: the frame's calibration factor (bin -> m/z) (toTof)
  vector<float> inten;
  vector<uint16_t> imq;      ///< the store's quantised 1/K0, kept for export
  vector<double> rt;         ///< size F until indexRt: the frame's retention time
  vector<uint32_t> frame_id; ///< size F until toTof: vendor frame Id, the calibration key
  bool mzq = true;           ///< tof[] holds m/z quanta, not bins, until toTof()
  /**
    @brief Number of frames.

    @return frame_off.size() - 1 (0 when empty)
  */
  Size frames() const
  { return frame_off.empty() ? 0 : frame_off.size() - 1; }
  /**
    @brief Number of peaks.

    @return tof.size()
  */
  Size peaks() const
  { return tof.size(); }
  /**
    @brief Bytes held by the per-peak arrays.

    @return The size in bytes
  */
  Size bytes() const
  { return tof.size() * 4 + inten.size() * 4 + imq.size() * 2; }
  /**
    @brief Append one picked frame.

    Capacity is projected from the running peaks-per-frame x
    `expected_frames` (the MS1 frame count: one MS2 frame per window per cycle, and MS1 arrives
    first on both readers), over-reserved 2.5x: untouched reserve costs no RSS, a reallocation
    copies the window. Never shrunk.

    @param[in] frame_rt Retention time
    @param[in] fid Vendor frame id
    @param[in] mzq_ m/z quanta
    @param[in] inten_ Intensities
    @param[in] imq_ Quantised 1/K0
    @param[in] n Number of peaks
    @param[in] expected_frames Projected number of frames
    @throws OpenMS::Exception::InvalidValue if the window would hold more than 2^32 peaks
  */
  void append(double frame_rt, uint32_t fid, const uint32_t* mzq_, const float* inten_, const uint16_t* imq_, Size n, Size expected_frames)
  {
    if (frame_off.empty()) { frame_off.push_back(0); }
    const Size n_frames = frames();
    auto grow = [&](auto& v, Size add) {
      if (v.size() + add <= v.capacity()) { return; }
      const double per_frame = n_frames ? (double)v.size() / (double)n_frames : (double)add;
      const Size want = (Size)(per_frame * (double)std::max(expected_frames, n_frames + 1) * 2.5) + add;
      v.reserve(std::max(v.size() + add, want));
    };
    if (n) // an EMPTY frame (a padded table frame) carries no peak pointers
    {
      grow(tof, n);
      grow(inten, n);
      grow(imq, n);
      tof.insert(tof.end(), mzq_, mzq_ + n);
      inten.insert(inten.end(), inten_, inten_ + n);
      imq.insert(imq.end(), imq_, imq_ + n);
    }
    rt.push_back(frame_rt);
    frame_id.push_back(fid);
    if (tof.size() > (Size)std::numeric_limits<uint32_t>::max()) // frame_off is 32-bit: a wrapped bound would mix frames
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "one isolation window holds more than 2^32 peaks; "
                                            "run in smaller tiles (-tile:cells_per_tile 1, a smaller -tile:rt_sec)",
                                            String(tof.size()));
    }
    frame_off.push_back((uint32_t)tof.size());
  }
};

/**
  @brief [cell] A window's FROZEN frame table, in RT order: frame times, vendor ids, RT-axis indices,
  calibration factors, and the MS1-frame -> window-frame map the Pearson uses.

  Store and scoring
  support come from here in every mode, so a tile holding only some frames scores exactly as the
  whole run does. Filled from the reader's tables (.d streaming source) or from the resident slabs.
*/
struct FrozenWin
{
  vector<double> frame_rt;
  vector<uint32_t> frame_id;
  vector<uint32_t> rt_index;
  vector<double> b;
  vector<int> nearest_local; ///< MS1-store frame -> nearest window frame within delta_rt, or -1
  /**
    @brief Number of frames in the table.

    @return frame_rt.size()
  */
  Size frames() const
  { return frame_rt.size(); }
};

/**
  @brief The MS1-frame -> window-frame map over ALL window frames (every frame carries fragment
  support), ties to the LATER frame.

  @param[in] win_rt Window frame times, sorted
  @param[in] ms1_rt MS1 frame times
  @param[in] delta_rt Maximum RT distance
  @return Per MS1 frame the nearest window frame, or -1 when none lies within @p delta_rt
*/
inline vector<int> nearestLocal(const vector<double>& win_rt, const vector<double>& ms1_rt, double delta_rt)
{
  vector<int> nl(ms1_rt.size(), -1);
  Size hi = 0;
  for (Size m = 0; m < ms1_rt.size(); ++m)
  {
    const double t = ms1_rt[m];
    while (hi < win_rt.size() && win_rt[hi] < t)
    {
      ++hi;
    }
    double bd = delta_rt;
    int gi = -1;
    if (hi > 0)
    {
      const double d = fabs(win_rt[hi - 1] - t);
      if (d <= bd)
      {
        bd = d;
        gi = (int)(hi - 1);
      }
    }
    if (hi < win_rt.size())
    {
      const double d = fabs(win_rt[hi] - t);
      if (d <= bd)
      {
        bd = d;
        gi = (int)hi;
      }
    }
    nl[m] = gi;
  }
  return nl;
}

/**
  @brief Ordered, incremental digest of a window's slab: frame table (rt index, factor bits, count) then
  every (tof, intensity bits, imq) tuple in frame order.

  Unlike detDigest it is NOT sorted, so a
  reordered equal-intensity peak or a tuple swapped between frames changes it -- which is what a
  representation change of the store must not do. [det]

  @param[in] sl The slab
  @param[in] seed Initial FNV-1a state
  @return The digest
*/
static uint64_t slabDigest(const PeakSlab& sl, uint64_t seed = 1469598103934665603ULL)
{
  uint64_t h = seed;
  auto mix = [&](uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
  };
  for (Size f = 0; f < sl.frames(); ++f)
  {
    uint64_t bb;
    std::memcpy(&bb, &sl.b[f], 8);
    mix(sl.rt_index[f]);
    mix(bb);
    mix(sl.frame_off[f + 1] - sl.frame_off[f]);
    for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
    {
      uint32_t ib;
      std::memcpy(&ib, &sl.inten[k], 4);
      mix(sl.tof[k]);
      mix(ib);
      mix(sl.imq[k]);
    }
  }
  return h;
}

/**
  @brief Frame times -> indices on the GLOBAL RT axis, once the axis exists; the frames must be in RT
  order (the extension loop steps frames by index and takes the RT span as last minus first --
  the OpenMS path sorts its spectra for the same reason).

  The loader delivers them ordered; this
  checks rather than assumes.

  @param[in,out] sl The slab; its rt_index is filled
  @throws OpenMS::Exception::InvalidValue if the frames are not in retention-time order
*/
void indexRt(PeakSlab& sl)
{
  const Size n_frames = sl.frames();
  sl.rt_index.resize(n_frames);
  for (Size i = 0; i < n_frames; ++i)
  {
    sl.rt_index[i] = rtIndex(sl.rt[i]);
  }
  for (Size i = 1; i < n_frames; ++i)
  {
    if (sl.rt_index[i] < sl.rt_index[i - 1])
    {
      throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "window frames are not in retention-time order at frame",
                                            std::to_string(i));
    }
  }
  vector<double>().swap(sl.rt);
}

/**
  @brief [cell] One RT cell of a window's slab -- frames [f0, f1) -- as a slab of its own, so the
  detector's seeds, `visited`, frame_live and tof index are the cell's alone.

  Traces found on
  it carry cell-local frame indices; the caller rebases them by f0 into the window store.

  @param[in] sl The window's slab
  @param[in] f0 First frame
  @param[in] f1 End frame (exclusive)
  @return The cell's slab
*/
PeakSlab sliceFrames(const PeakSlab& sl, Size f0, Size f1)
{
  PeakSlab s;
  const uint32_t p0 = sl.frame_off[f0], p1 = sl.frame_off[f1];
  s.frame_off.assign(sl.frame_off.begin() + f0, sl.frame_off.begin() + f1 + 1);
  for (auto& o : s.frame_off)
  {
    o -= p0;
  }
  s.rt_index.assign(sl.rt_index.begin() + f0, sl.rt_index.begin() + f1);
  s.b.assign(sl.b.begin() + f0, sl.b.begin() + f1);
  s.tof.assign(sl.tof.begin() + p0, sl.tof.begin() + p1);
  s.inten.assign(sl.inten.begin() + p0, sl.inten.begin() + p1);
  s.imq.assign(sl.imq.begin() + p0, sl.imq.begin() + p1);
  s.mzq = sl.mzq;
  return s;
}

/**
  @brief m/z quanta -> flight-time bins, in place.

  The store keeps m/z at a 1e-5 Da quantum, ~400x finer
  than one TOF bin at m/z 600, so inverting the calibration recovers the ORIGINAL bin exactly --
  no approximation is introduced here, and the loader does not have to know the axis.

  [merge] The IM picker's 1/K0 box (gate:delta_im) cuts an ion whose mobility peak is wider into several centroids per
  frame that land in ONE bin; the tracer takes one peak per frame, so each slice would become its own trace and a
  spectrum could carry the same m/z twice. Here, where the bin first exists, a frame's peaks that share a bin and
  chain within trace:im_slice_merge in 1/K0 become one peak: intensity summed, 1/K0 intensity-weighted. Frame-local,
  so thread and tile invariant; peaks stay ascending by bin within a frame.

  @param[in,out] sl The slab; converted from m/z quanta to flight-time bins
*/
void toTof(PeakSlab& sl)
{
  if (! sl.mzq) { return; }
  const TofAxis& ax = tofAxis();
  const Size n_frames = sl.frames();
  sl.b.resize(n_frames);
  // single linkage along 1/K0 with a link of cfg_im_slice_merge, the whole group within 2x that (one ion's 2-3 slices)
  const long tol_q = cfg_im_slice_merge >= 0.0 ? (long)llround(cfg_im_slice_merge * IM_Q) : -1, span_q = 2 * tol_q;
  // [serial] Frames convert and merge independently, so they run as tasks over frame chunks, each frame compacted
  // within its own range (its cursor starts at its first peak and never passes the read position, as the global
  // cursor never did); the frames are then closed up in order -- the one-cursor loop's layout, byte for byte.
  const long nt = parTasks(sl.tof.size(), 1 << 17);
  vector<uint32_t> nkept(n_frames);        // the frame's peaks after the merge
  vector<long long> mc(nt, 0), fmc(nt, 0); // per task: peaks merged away, frames with a merge
  forTasks(nt, [&](long ch) {
    vector<pair<uint16_t, uint32_t>> run;             // (imq, slab index) of one same-bin run
    vector<std::tuple<uint32_t, float, uint16_t>> cl; // (first slab index, intensity, imq) of the run's groups
    long long merged = 0, frames_merged = 0;
    for (Size i = (Size)ch * n_frames / (Size)nt; i < (Size)(ch + 1) * n_frames / (Size)nt; ++i)
    {
      const double b = ax.factor(sl.frame_id[i]);
      sl.b[i] = b;
      const uint32_t lo = sl.frame_off[i], hi = sl.frame_off[i + 1];
      for (uint32_t k = lo; k < hi; ++k)
      {
        sl.tof[k] = ax.tofOf((double)sl.tof[k] / MZ_Q, b);
      }
      uint32_t w = lo;
      bool any = false;
      for (uint32_t k = lo; k < hi;)
      {
        uint32_t e = k + 1;
        while (e < hi && sl.tof[e] == sl.tof[k])
        {
          ++e;
        }
        if (e - k == 1 || tol_q < 0)
        {
          for (uint32_t j = k; j < e; ++j)
          {
            sl.tof[w] = sl.tof[j];
            sl.inten[w] = sl.inten[j];
            sl.imq[w] = sl.imq[j];
            ++w;
          }
          k = e;
          continue;
        }
        run.clear();
        for (uint32_t j = k; j < e; ++j)
        {
          run.emplace_back(sl.imq[j], j);
        }
        std::sort(run.begin(), run.end()); // by 1/K0, then slab order: deterministic
        cl.clear();                        // every read of the run happens before any write
        for (Size r = 0; r < run.size();)
        {
          Size q = r + 1;
          while (q < run.size() && (long)run[q].first - (long)run[q - 1].first <= tol_q && (long)run[q].first - (long)run[r].first <= span_q)
          {
            ++q;
          }
          if (q - r == 1)
          {
            cl.emplace_back(run[r].second, sl.inten[run[r].second], run[r].first);
            r = q;
            continue;
          } // untouched
          double isum = 0.0, wim = 0.0;
          uint32_t first = run[r].second;
          for (Size u = r; u < q; ++u)
          {
            first = std::min(first, run[u].second);
            const float v = sl.inten[run[u].second];
            if (v > 0.0f && std::isfinite(v))
            {
              isum += v;
              wim += (double)v * run[u].first;
            }
          }
          cl.emplace_back(first, (float)isum, isum > 0.0 ? (uint16_t)std::min(65535.0, std::max(0.0, std::round(wim / isum))) : run[r].first);
          merged += (long long)(q - r - 1);
          any = true;
          r = q;
        }
        // groups written at their first member's slab position: with nothing merged the slab is byte-identical,
        // so the tracer's same-bin tie-break (slab order) is unchanged wherever no merge happened
        std::sort(cl.begin(), cl.end(), [](const auto& x, const auto& y) { return std::get<0>(x) < std::get<0>(y); });
        const TofIdx t = sl.tof[k];
        for (const auto& c : cl)
        {
          sl.tof[w] = t;
          sl.inten[w] = std::get<1>(c);
          sl.imq[w] = std::get<2>(c);
          ++w;
        }
        k = e;
      }
      nkept[i] = w - lo;
      frames_merged += any;
    }
    mc[ch] = merged;
    fmc[ch] = frames_merged;
  });
  long long merged = 0, frames_merged = 0;
  for (long c = 0; c < nt; ++c)
  {
    merged += mc[c];
    frames_merged += fmc[c];
  }
  uint32_t w = 0;
  for (Size i = 0; i < n_frames; ++i) // close up: every destination lies at or before its source
  {
    const uint32_t lo = sl.frame_off[i], n = nkept[i];
    if (w != lo)
    {
      std::copy(sl.tof.begin() + lo, sl.tof.begin() + lo + n, sl.tof.begin() + w);
      std::copy(sl.inten.begin() + lo, sl.inten.begin() + lo + n, sl.inten.begin() + w);
      std::copy(sl.imq.begin() + lo, sl.imq.begin() + lo + n, sl.imq.begin() + w);
    }
    sl.frame_off[i] = w;
    w += n;
  }
  if (n_frames) { sl.frame_off[n_frames] = w; }
  sl.tof.resize(w);
  sl.inten.resize(w);
  sl.imq.resize(w);
  mergeCount()[0] += merged;
  mergeCount()[1] += frames_merged;
  vector<uint32_t>().swap(sl.frame_id);
  sl.mzq = false;
}


/**
  @brief Per-cell, read-only preparation shared by every band task of detectTracesInteger: per-frame
  "has an above-noise peak" (OpenMS deletes sub-noise peaks before asking whether a frame is
  empty), the intensity-ordered seed list segmented by band, and (banded) the tof index.
*/
struct TracePrep
{
  vector<char> frame_live; ///< F: has at least one member
  vector<uint32_t> order;  ///< seeds, intensity desc, later-index-first among equals
  vector<Size> seg;        ///< nb+1 offsets into `order`: band b owns [seg[b], seg[b+1])
  /// [tofidx] Per-frame direct-address index over the flight-time bins: bucket `(tof - tlo) >> tsh`
  /// of frame f starts at peak `fidx[f * (nbk+1) + bucket]`, so best() narrows its search by arithmetic.
  vector<uint32_t> fidx;
  TofIdx tlo = 0;
  unsigned tsh = 0;
  uint32_t nbk = 0;
  /**
    @brief The sub-range of frame @p f's peaks [a, b) that can hold the first peak with tof >= @p lo_t (from the tof index).

    @param[in] f Frame
    @param[in] lo_t Lowest flight-time bin of interest
    @param[in] a Frame's first peak
    @param[in] b Frame's end peak
    @param[out] s0 Start of the sub-range
    @param[out] s1 End of the sub-range
  */
  inline void range(Size f, TofIdx lo_t, uint32_t a, uint32_t b, uint32_t& s0, uint32_t& s1) const
  {
    if (fidx.empty())
    {
      s0 = a;
      s1 = b;
      return;
    }
    const uint64_t d = lo_t <= tlo ? 0 : ((uint64_t)(lo_t - tlo) >> tsh);
    const uint32_t k = d >= nbk ? nbk : (uint32_t)d;
    const uint32_t* row = &fidx[f * ((Size)nbk + 1)];
    s0 = row[k];
    s1 = (k + 1 <= nbk) ? row[k + 1] : b;
  }
};
/**
  @brief [tofidx] Bucket every frame's flight-time bins so `best()` can address its scan start instead
  of binary-searching for it.

  ~32 peaks per bucket; the shift is global so bucket boundaries
  line up across frames. Exact: the bucket start is the first peak with `tof >= bucket_lo`, so
  the sub-range handed to lower_bound provably contains the same position it would have found.

  @param[in] sl The slab (in flight-time bins)
  @param[in,out] tp Receives the index
  @param[in] tlo Lowest bin
  @param[in] thi Highest bin
*/
void buildTofIndex(const PeakSlab& sl, TracePrep& tp, TofIdx tlo, TofIdx thi)
{
  const Size n_frames = sl.frames(), n_points = sl.tof.size();
  if (n_frames == 0 || n_points == 0 || thi < tlo) { return; }
  const uint64_t span = (uint64_t)thi - tlo + 1;
  const Size ppf = std::max<Size>(1, n_points / n_frames);
  uint32_t want = 64;
  while (want < 4096 && (Size)want * 32 < ppf)
  {
    want <<= 1;
  }
  unsigned sh = 0;
  while ((span >> sh) > want)
  {
    ++sh;
  }
  const uint32_t nbk = (uint32_t)std::min<uint64_t>(want, (span >> sh) + 1);
  // one row of nbk+1 offsets per frame; the extra entry is the row's end
  if ((double)n_frames * (nbk + 1) > 3e8)
  {
    return; // pathological geometry: skip the index
  }
  tp.tlo = tlo;
  tp.tsh = sh;
  tp.nbk = nbk;
  tp.fidx.assign(n_frames * ((Size)nbk + 1), 0);
#pragma omp taskloop default(shared)
  for (long fi = 0; fi < (long)n_frames; ++fi)
  {
    const Size f = (Size)fi;
    const uint32_t a = sl.frame_off[f], e = sl.frame_off[f + 1];
    uint32_t* row = &tp.fidx[f * ((Size)nbk + 1)];
    uint32_t k = a;
    for (uint32_t bk = 0; bk <= nbk; ++bk)
    {
      const uint64_t blo = (uint64_t)bk << sh;
      while (k < e && (uint64_t)(sl.tof[k] - tlo) < blo)
      {
        ++k;
      }
      row[bk] = k;
    }
  }
}

/**
  @brief Build the read-only tracing preparation of one cell.

  `bnd` are the nb+1 flight-time band edges the detector will be given (bnd[0] = tlo,
  bnd[nb] = thi+1); pass {0, 0} for the unbanded path, which wants one globally sorted list.
  [seg] The seed list is segmented by band and each segment sorted in parallel. Exact because the
  comparator is a total order: filter_b(sort(all)) == sort(filter_b(all)).

  @param[in] sl The cell's slab
  @param[in] noise Noise threshold
  @param[in] snr Seed gate: apex > snr * noise
  @param[in] bnd Band edges, or {0, 0}
  @param[in] tlo Lowest bin
  @param[in] thi Highest bin
  @return The preparation
  @par Exceptions
    Rethrows the first exception of a parallel part, of whatever type it was.
*/
TracePrep prepareTracing(const PeakSlab& sl, double noise, double snr, const vector<TofIdx>& bnd, TofIdx tlo, TofIdx thi)
{
  TracePrep tp;
  const Size n_frames = sl.frames(), n_points = sl.tof.size();
  // Membership IS inten > noise (no per-peak map); seeds are counted first so `order` is sized exactly.
  tp.frame_live.assign(n_frames, 0);
  const double min_apex = snr * noise;
  // OpenMS: stable_sort ascending by intensity, then iterate in reverse -> among equal
  // intensities the LATER peak seeds first.
  auto cmp = [&](uint32_t a, uint32_t b) {
    if (sl.inten[a] != sl.inten[b]) { return sl.inten[a] > sl.inten[b]; }
    return a > b;
  };
  const int nb = (int)bnd.size() - 1;
  const bool banded = nb >= 2; // unbanded: one globally sorted list
  const Size nseg = banded ? (Size)nb : 1;
  // Band of a peak: the same half-open [bnd[b], bnd[b+1]) partition the detector applies, found
  // by upper_bound so that empty bands (bnd[b] == bnd[b+1], reachable when nb exceeds the tof
  // range) are skipped exactly as the detector's own range test skips them.
  // Clamped at both ends: a bin outside [bnd[0], bnd[nb]) would otherwise index cnt/cur out of
  // bounds. The window's edges are asserted to enclose every cell's bins before this runs (a
  // loud throw); the clamp is the memory-safety floor under that assert.
  auto band_of = [&](TofIdx t) {
    if (! banded) { return (Size)0; }
    const auto it = std::upper_bound(bnd.begin(), bnd.end(), t);
    if (it == bnd.begin()) { return (Size)0; }
    return std::min((Size)(it - bnd.begin() - 1), (Size)(nb - 1));
  };
  // [serial] Count, then distribute, as tasks over peak-balanced FRAME chunks (frame_live stays chunk-private):
  // chunk c counts its seeds per band, the prefix over (band, chunk) is its write cursor, so every band receives
  // its seeds in ascending peak index -- what the one-pass loops wrote. Chunk c owns peaks [kb(c), kb(c + 1)).
  const long nt = parTasks(n_points, 1 << 17);
  vector<Size> fc(nt + 1, n_frames);
  fc[0] = 0;
  for (long c = 1; c < nt; ++c)
  {
    fc[c] = (Size)(std::lower_bound(sl.frame_off.begin(), sl.frame_off.begin() + n_frames, (uint32_t)((Size)c * n_points / (Size)nt))
                   - sl.frame_off.begin());
  }
  auto kb = [&](long c) { return fc[c] >= n_frames ? (uint32_t)n_points : sl.frame_off[fc[c]]; };
  vector<Size> cnt((Size)nt * nseg, 0);
  forTasks(nt, [&](long c) {
    Size* cc = &cnt[(Size)c * nseg];
    for (Size f = fc[c]; f < fc[c + 1]; ++f)
    {
      for (uint32_t k = sl.frame_off[f]; k < sl.frame_off[f + 1]; ++k)
      {
        if ((double)sl.inten[k] > noise)
        {
          tp.frame_live[f] = 1;
          if ((double)sl.inten[k] > min_apex) { ++cc[band_of(sl.tof[k])]; }
        }
      }
    }
  });
  tp.seg.assign(nseg + 1, 0);
  vector<Size> cur((Size)nt * nseg);
  Size n_seed = 0;
  for (Size b = 0; b < nseg; ++b)
  {
    tp.seg[b] = n_seed;
    for (long c = 0; c < nt; ++c)
    {
      cur[(Size)c * nseg + b] = n_seed;
      n_seed += cnt[(Size)c * nseg + b];
    }
  }
  tp.seg[nseg] = n_seed;
  tp.order.resize(n_seed);
  forTasks(nt, [&](long c) {
    Size* cc = &cur[(Size)c * nseg];
    for (uint32_t k = kb(c), ke = kb(c + 1); k < ke; ++k)
    {
      if ((double)sl.inten[k] > noise && (double)sl.inten[k] > min_apex) { tp.order[cc[band_of(sl.tof[k])]++] = k; }
    }
  });
  // [serial] A segment holding many seeds (the dense band) is itself sorted as tasks. The comparator is a TOTAL
  // order (the index breaks every tie), so the sorted segment is unique: parallelSort returns std::sort's result.
  auto sort_seg = [&](Size b) {
    const Size n = tp.seg[b + 1] - tp.seg[b];
    parallelSort(tp.order.data() + tp.seg[b], n, cmp, parTasks(n, 1 << 18));
  };
  if (! banded)
  {
    sort_seg(0);
    return tp;
  }
  std::exception_ptr sort_err; // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
  for (int b = 0; b < nb; ++b)
  {
    try
    {
      sort_seg((Size)b);
    }
    catch (...)
    {
      taskErr(sort_err);
    }
  }
  if (sort_err) { std::rethrow_exception(sort_err); }
  buildTofIndex(sl, tp, tlo, thi);
  return tp;
}

/**
  @brief Mass-trace detection on integer arrays with the semantics of OpenMS MassTraceDetection::run_
  (written from that source):
    * seeds in intensity order, ties later-index-first;
    * the two directions extend interleaved, one frame down then one up, sharing one centroid;
    * accept the CLOSEST above-noise peak (in m/z) within +-3 sd of the running centroid and inside
      the mobility tolerance; if it is taken the frame is a MISS (the next-closest is not tried);
    * sd starts at centroid*ppm and is re-estimated from accepted peaks;
    * a direction stops after more than 5 consecutive misses; an empty frame is not a miss;
    * valid iff RT span >= min length and size / (visited - trailing misses) >= min_sample_rate;
      peaks are marked taken only for a valid trace; the seed counts twice in the centroid.

  It does not reproduce the OpenMS detector peak for peak: judge it on identified peptides with both
  engines, never on a digest. Valley splitting runs afterwards.

  @param[in] sl The cell's slab (in flight-time bins)
  @param[in] tp Its preparation
  @param[in,out] store Receives the trace spans
  @param[in] mass_ppm Mass tolerance (ppm)
  @param[in] im_tol Ion-mobility tolerance (1/K0)
  @param[in] noise Noise threshold
  @param[in] min_len_sec Minimum trace length (seconds)
  @param[in] band Band index into the preparation's seed segments, -1 for the unbanded path
  @return The traces
  @throws OpenMS::Exception::OutOfRange if a bin or frame offset does not fit into 16 bits
*/
vector<Trace> detectTracesInteger(const PeakSlab& sl,
                                  const TracePrep& tp,
                                  TraceStore& store,
                                  double mass_ppm,
                                  double im_tol,
                                  double noise,
                                  double min_len_sec,
                                  int band = -1)
{
  vector<Trace> out;
  const Size n_frames = sl.frames(), n_points = sl.tof.size();
  if (n_frames == 0 || n_points == 0) { return out; }
  const TofAxis& ax = tofAxis();
  const double msr = 0.5; // MassTraceDetection's default min_sample_rate
  const double imq_tol = im_tol * IM_Q;
  const Size max_missed = 5; // trace_termination_outliers

  vector<uint64_t> visited((n_points + 63) / 64, 0); // one BIT per peak
  auto seen = [&](Size k) { return (visited[k >> 6] >> (k & 63)) & 1u; };
  auto mark = [&](Size k) { visited[k >> 6] |= (uint64_t)1 << (k & 63); };
  auto mz_at = [&](uint32_t k, Size f) { return ax.mzOf(sl.tof[k], sl.b[f]); };
  auto iwm = [](double v, double w, double& c, double& cnt, double& den) { // updateIterativeWeightedMean_
    const double ct = 1.0 + (w * v) / cnt, dt = 1.0 + w / den;
    c *= (ct / dt);
    cnt += w * v;
    den += w;
  };
  auto sd_robust = [](double mz, double w, double mean, double& sd, double& wsum) { // updateWeightedSDEstimateRobust
    const double d1 = pmath::log(wsum) + 2.0 * pmath::log(sd);
    const double d2 = pmath::log(w) + 2.0 * pmath::log(std::abs(mz - mean));
    const double denom = std::sqrt(pmath::exp(d1) + pmath::exp(d2));
    const double ws = wsum + w;
    const double t = denom / std::sqrt(ws);
    if (t > std::numeric_limits<double>::epsilon()) { sd = t; }
    wsum = ws;
  };
  /// findBestPeak_ + isPeakAcceptable_: the closest peak in m/z within +-3 sd of the running
  /// centroid whose mobility is inside the gate; -1 for "no candidate", and `taken` set when the
  /// closest one is already used, which OpenMS counts as a miss rather than trying the next.
  /// The +-3 sd window is converted to a BIN range once per frame step, so the walk is bounded
  /// by two calibration calls and a binary search; only candidates inside that window are
  /// converted back to m/z for the comparison.
  auto best = [&](Size f, double cmz, double sd, double cim, bool& taken) -> long {
    taken = false;
    const uint32_t a = sl.frame_off[f], b = sl.frame_off[f + 1];
    if (b <= a) { return -1; }
    const double bf = sl.b[f];
    const double lo_mz = cmz - 3.0 * sd, hi_mz = cmz + 3.0 * sd;
    const TofIdx lo_t = ax.tofOf(lo_mz > 0.0 ? lo_mz : 1e-9, bf);
    const TofIdx hi_t = ax.tofOf(hi_mz, bf);
    uint32_t r0 = a, r1 = b;
    tp.range(f, lo_t, a, b, r0, r1); // [tofidx] computed, not searched
    const uint32_t s0 = (uint32_t)(lower_bound(sl.tof.begin() + r0, sl.tof.begin() + r1, lo_t) - sl.tof.begin());
    long bk = -1;
    double bd = std::numeric_limits<double>::infinity();
    for (uint32_t k = s0; k < b && sl.tof[k] <= hi_t; ++k)
    {
      if (! ((double)sl.inten[k] > noise)) { continue; }
      if (std::abs((double)sl.imq[k] - cim) > imq_tol) { continue; }
      // "Closest" in m/z, as OpenMS does, not in bins: across the centroid the calibration's
      // curvature can flip which is nearer. Only candidates inside the bin window pay the conversion.
      const double d = std::abs(ax.mzOf(sl.tof[k], bf) - cmz);
      if (d < bd)
      {
        bd = d;
        bk = (long)k;
      } // strict <: first of equals wins
    }
    if (bk >= 0 && seen((Size)bk))
    {
      taken = true;
      return -1;
    }
    return bk;
  };

  struct Dir
  {
    Size missed = 0, scans = 0;
    bool active = true;
  };
  vector<pair<uint32_t, uint32_t>> got; // (frame, peak); marked only if the trace is valid
  // [seg] Only this band's own segment of the seed list (band -1: the whole list). The segment holds exactly the
  // band's seeds: prepareTracing partitions by the same edges, and every cell bin lies inside them (checked).
  const bool segd = band >= 0;
  const Size s0 = segd ? tp.seg[band] : 0, s1 = segd ? tp.seg[band + 1] : tp.order.size();
  for (Size si = s0; si < s1; ++si)
  {
    const uint32_t seed = tp.order[si];
    if (seen(seed)) { continue; }
    const Size f0 = (uint32_t)(std::upper_bound(sl.frame_off.begin(), sl.frame_off.end(), seed) - sl.frame_off.begin() - 1);
    const double w0 = sl.inten[seed], m0 = mz_at(seed, f0);
    double cmz = m0, cnt = w0 * m0, den = w0;
    iwm(m0, w0, cmz, cnt, den); // seed counted twice, as run_ does
    double cim = sl.imq[seed], cnt_im = w0 * cim, den_im = w0;
    iwm((double)sl.imq[seed], w0, cim, cnt_im, den_im);
    double sd = (cmz / 1e6) * mass_ppm, wsum = w0;
    got.clear();
    got.emplace_back((uint32_t)f0, seed);
    Dir dn, up;
    long fd = (long)f0, fu = (long)f0;
    while ((fd > 0 && dn.active) || (fu + 1 < (long)n_frames && up.active))
    {
      for (int side = 0; side < 2; ++side)
      {
        Dir& d = side == 0 ? dn : up;
        long& fi = side == 0 ? fd : fu;
        if (! d.active) { continue; }
        if (side == 0 ? ! (fi > 0) : ! (fi + 1 < (long)n_frames)) { continue; }
        const Size f = (Size)(side == 0 ? fi - 1 : fi + 1);
        // "empty" AFTER sub-noise deletion, as OpenMS tests it; an empty frame is not a miss.
        if (tp.frame_live[f])
        {
          bool taken = false;
          const long k = best(f, cmz, sd, cim, taken);
          if (k >= 0)
          {
            const double mz = mz_at((uint32_t)k, f), w = sl.inten[k];
            iwm(mz, w, cmz, cnt, den);
            iwm((double)sl.imq[k], w, cim, cnt_im, den_im);
            sd_robust(mz, w, cmz, sd, wsum); // reestimate_mt_sd: the tool never turns it off
            got.emplace_back((uint32_t)f, (uint32_t)k);
            d.missed = 0;
          }
          else
          {
            ++d.missed;
          }
        }
        fi += (side == 0 ? -1 : 1);
        ++d.scans;
        if (d.missed > max_missed) { d.active = false; }
      }
    }
    // isTraceValid_
    uint32_t fmin = got[0].first, fmax = fmin;
    for (const auto& g : got)
    {
      fmin = std::min(fmin, g.first);
      fmax = std::max(fmax, g.first);
    }
    const double span_s = rtAxis()[sl.rt_index[fmax]] - rtAxis()[sl.rt_index[fmin]];
    if (span_s < min_len_sec) { continue; }
    const Size total = dn.scans + up.scans + 1;
    const Size adjusted = total > dn.missed + up.missed ? total - dn.missed - up.missed : 1;
    if ((double)got.size() / (double)adjusted < msr) { continue; }
    for (const auto& g : got)
    {
      mark(g.second); // valid: now they are taken
    }

    // the span: (store-local frame, intensity) points, plus the per-frame bin arena the valley
    // splitter needs to give each child ITS measured apex bin
    vector<pair<uint32_t, float>> pts;
    pts.reserve(got.size());
    uint32_t kap = got[0].second, fap = got[0].first;
    for (const auto& g : got)
    {
      pts.emplace_back(g.first, sl.inten[g.second]);
      if (sl.inten[g.second] > sl.inten[kap])
      {
        kap = g.second;
        fap = g.first;
      }
    }
    Trace t = makeSpan(pts, store);
    const uint32_t t_tof = sl.tof[kap]; // this trace's own bin: the offsets' origin
    store.bins.resize(store.inten.size(), 0);
    for (const auto& g : got)
    {
      const int64_t d = (int64_t)sl.tof[g.second] - (int64_t)t_tof; // [bin-delta]
      if (d < -32768 || d > 32767) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
      store.bins[t.off + (g.first - t.frame0)] = (int16_t)d;
    }
    t.tof = t_tof;
    t.bframe = fap; // got-order first max: the frame the bin was measured on
    t.mz = cmz;     // the weighted centroid, as OpenMS reports
    {
      const int64_t d = (int64_t)fap - (int64_t)t.frame0; // [soa-rt] the apex frame, 16-bit as in setRtFrame
      if (d < -32768 || d > 32767) { throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
      t.rt_at = (int16_t)d;
    }
    t.im = IM_LO + cim / IM_Q;
    out.push_back(std::move(t));
  }
  return out;
}

/**
  @brief The acquisition's cycle time: the smallest POSITIVE gap between consecutive frame times.

  Every gap
  is a multiple of the cycle, so skipped frames cannot inflate it (a median can). It is the valley
  splitter's scan time and the tile:rt_sec cell unit; 0 when no two
  frames differ, and the splitter then uses each trace's own average.

  @param[in] frame_rt Frame times
  @return The cycle time, or 0
*/
inline double frameSpacing(const vector<double>& frame_rt)
{
  double m = 0.0;
  for (Size i = 1; i < frame_rt.size(); ++i)
  {
    const double d = frame_rt[i] - frame_rt[i - 1];
    if (d > 0.0 && (m == 0.0 || d < m)) { m = d; }
  }
  return m;
}

/**
  @brief ElutionPeakDetection as a valley splitter: it splits, never deletes (width_filtering off, no SNR filter).

  @param[in,out] epd The detector to configure
  @param[in] chrom_fwhm Expected chromatographic FWHM
  @param[in] scan_time Scan time, 0 = per trace
*/
inline void setSplitter(ElutionPeakDetection& epd, double chrom_fwhm, double scan_time)
{
  epd.setLogType(ProgressLogger::NONE);
  Param ep = epd.getParameters();
  ep.setValue("chrom_fwhm", chrom_fwhm);
  ep.setValue("width_filtering", "off");
  ep.setValue("masstrace_snr_filtering", "false");
  ep.setValue("scan_time", scan_time); // 0 = the trace's own average cycle time
  epd.setParameters(ep);
}

/**
  @brief Integer traces -> MassTrace -> ElutionPeakDetection -> Trace, each child's flight-time bin
  re-attached from its apex frame's calibration.

  Chunked and parallel like the OpenMS path's split;
  only the traces become OpenMS objects, never the peak arrays.

  @param[in,out] in The traces to split (consumed)
  @param[in] split_valleys Valley-splitting parameter
  @param[in,out] wst The window store
  @param[in] nchunk Number of chunks
  @param[in] scan_time Scan time, 0 = per trace
  @return The split traces
  @par Exceptions
    Rethrows the first exception of a chunk, of whatever type it was.
*/
vector<Trace> splitIntegerTraces(vector<Trace>& in, double split_valleys, TraceStore& wst, int nchunk, double scan_time)
{
  // No `in.size() < 2` guard: a lone parent must split as it would among others (tile invariance).
  if (split_valleys <= 0.0 || in.empty()) { return std::move(in); }
  const TofAxis& ax = tofAxis();
  // [serial] Every parent splits on its own (below), so neither the chunk count nor the cut points can change the
  // output. A large window (>= 8192 parents per chunk) gets up to one chunk per team thread instead of `nchunk`,
  // and the cuts sit at equal cumulative cost (real points + 8 per trace) rather than equal counts: a window's
  // parents come band by band, most intense (longest) first, so equal counts put the longest traces together.
  const Size n_in = in.size();
  const int team = omp_in_parallel() ? omp_get_num_threads() : 1;
  nchunk = std::max(nchunk, std::min<int>(team, (int)std::min<Size>(n_in / 8192, 1 << 20)));
  nchunk = std::max(1, std::min<int>(nchunk, (int)n_in));
  vector<Size> cut(nchunk + 1, n_in);
  cut[0] = 0;
  {
    uint64_t tot = 0, acc = 0;
    int c = 1;
    for (const Trace& t : in)
    {
      tot += (uint64_t)t.npts + 8;
    }
    for (Size q = 0; q < n_in && c < nchunk; ++q)
    {
      acc += (uint64_t)in[q].npts + 8;
      while (c < nchunk && acc * (uint64_t)nchunk >= tot * (uint64_t)c)
      {
        cut[c++] = q + 1;
      }
    }
  }
  vector<vector<Trace>> parts(nchunk);
  vector<TraceStore> stores(nchunk); // chunk-private arenas, merged after
  std::exception_ptr task_err;       // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
  for (int ci = 0; ci < nchunk; ++ci)
  {
    try
    {
      const Size lo = cut[ci], hi = cut[ci + 1];
      TraceStore& cst = stores[ci];
      cst.setFrames(wst.rt_index, &wst.b);
      // One MassTrace in flight at a time. Per-trace detectPeaks equals the batch call: with
      // width_filtering=off and masstrace_snr_filtering=false no cross-trace routine runs, and the
      // patched batch path, nested as it is here, is itself a per-trace loop appending in the same
      // order. epd is configured once per chunk: per trace that would be ~4e4 Param round-trips.
      ElutionPeakDetection epd;
      setSplitter(epd, split_valleys, scan_time);
      vector<Trace>& out = parts[ci];
      out.reserve((hi - lo) + (hi - lo) / 3); // headroom for the children the split adds
      vector<Peak2D> pk;                      // reused: a vector, not a per-point list
      vector<MassTrace> split;                // reused: detectPeaks clears it itself
      for (Size q = lo; q < hi; ++q)
      {
        const Trace& t = in[q];
        pk.clear();
        for (Size k = 0; k < t.span(); ++k) // REAL points only: EPD smooths what it is given
        {
          if (! real(wst, t, k)) { continue; }
          Peak2D pt;
          pt.setRT(rtAtSpan(wst, t, k));
          pt.setIntensity(xv(wst, t, k));
          // [bin-delta] reconstruct the absolute bin
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
          if (c.np() == 0) { continue; }
          const uint32_t cf = c.frame0 + c.apex; // apex frame, window-local
          c.bframe = cf;
          const double cb = wst.b.empty() ? ax.factor(0) : wst.b[cf];
          const double amz = cm.getSize() ? cm[cm.findMaxByIntPeak(false)].getMZ() : c.mz;
          c.tof = ax.tofOf(amz, cb);
          out.push_back(std::move(c));
        }
      }
    }
    catch (...)
    {
      taskErr(task_err);
    }
  }
  if (task_err) { std::rethrow_exception(task_err); }
  // the parents' spans are dead: rebuild the window arena from the chunk arenas
  vector<Trace>().swap(in);
  vector<float>().swap(wst.inten);
  vector<int16_t>().swap(wst.bins);
  Size tot = 0, arena = 0;
  bool any_bins = false;
  for (int ci = 0; ci < nchunk; ++ci)
  {
    tot += parts[ci].size();
    arena += stores[ci].inten.size();
    any_bins |= ! stores[ci].bins.empty();
  }
  // Reserve the merged arena once; the chunk arenas stay live until absorbed.
  wst.inten.reserve(arena);
  if (any_bins) { wst.bins.reserve(arena); }
  vector<Trace> out;
  out.reserve(tot);
  vector<std::pair<TraceStore*, vector<Trace>*>> src(nchunk); // [serial] absorbed in chunk order, copies as tasks
  for (int ci = 0; ci < nchunk; ++ci)
  {
    src[ci] = {&stores[ci], &parts[ci]};
  }
  absorbAll(wst, out, src, vector<uint32_t>(nchunk, 0));
  return out;
}

/**
  @brief Rebuild a PeakMap for ONE window from its slab (still in m/z quanta), for MassTraceDetection.

  The slab is consumed.

  @param[in,out] sl The slab (consumed)
  @return One spectrum per frame
  @throws OpenMS::Exception::InvalidValue if the slab is no longer in m/z quanta
*/
PeakMap materializeWindow(PeakSlab& sl)
{
  if (! sl.mzq)
  {
    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "materializeWindow needs a slab in m/z quanta", "tof");
  }
  PeakMap m;
  const Size n_frames = sl.frames();
  m.reserve(n_frames);
  for (Size f = 0; f < n_frames; ++f)
  {
    MSSpectrum s;
    s.setRT(rtAxis()[sl.rt_index[f]]);
    s.setMSLevel(1); // MassTraceDetection traces MS1-level only
    const Size n = sl.frame_off[f + 1] - sl.frame_off[f];
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
  sl = PeakSlab(); // consumed
  return m;
}

/// Isolation window identity: (lo*100, hi*100, WindowGroup). The group is part of it because one
/// diaPASEF scheme (PXD017703 "py3") acquires the SAME m/z window in two groups with shifted,
/// overlapping mobility slices; keyed by m/z alone the halves would merge into one non-monotone
/// window. With one slice per m/z the key orders and partitions exactly as (lo, hi) did.
using WinKey = std::array<int, 3>;
/**
  @brief The key of an isolation window.

  @param[in] lo Lower m/z bound
  @param[in] hi Upper m/z bound
  @param[in] group Window group
  @return The key
*/
inline WinKey winKey(double lo, double hi, uint32_t group)
{ return {(int)llround(lo * 100.0), (int)llround(hi * 100.0), (int)group}; }
/**
  @brief The isolation bounds precursors are routed on: the key's own, i.e. rounded to 0.01 Th (1.7.0's rule; exact bounds
  would change the spectra of acquisitions whose edges lie off that grid).

  @param[in] k The key
  @return (lo, hi)
*/
inline std::pair<double, double> winBounds(const WinKey& k)
{ return {k[0] / 100.0, k[1] / 100.0}; }
/**
  @brief A window's m/z range as the logs and refusals name it: "lo-hi".

  @param[in] k The key
  @return "lo-hi"
*/
inline String winLabel(const WinKey& k)
{ return String(k[0] / 100.0) + "-" + String(k[1] / 100.0); }

/**
  @brief A refusal of the input: structural frame tables [frames], a non-finite MS1 m/z [ms1-finite], the dnoise MS1 port
  [dnoise].

  Rethrown past the stream load's fallback, which would otherwise load the whole run again onto the partial
  slabs and MS1 map it leaves behind, only to refuse it there.
*/
struct InputRefusal : Exception::InvalidValue
{ using Exception::InvalidValue::InvalidValue; };

/// The MS1 band partition's halo in units of trace:mass_error_ppm: detectTraces_ widens every band by
/// this many tolerances on each side, and the load-time MS1 prune is exact only because its hop is
/// derived from the same number. There is deliberately only this one. [ms1-prune]
constexpr double BAND_HALO_PPM_MUL = 20.0;

/**
  @brief MSSpectrum::select right-sizes the peaks but reserves every data array at its OLD size (OpenMS MSSpectrum.cpp:49,
  71, 93): without this a selected spectrum keeps its Ion Mobility array's capacity for every peak it lost, 4 B
  each.

  [ms1-prune] [dnoise]

  @param[in,out] s The spectrum
*/
inline void shrinkDataArrays(MSSpectrum& s)
{
  for (auto& a : s.getFloatDataArrays())
  {
    a.shrink_to_fit();
  }
  for (auto& a : s.getIntegerDataArrays())
  {
    a.shrink_to_fit();
  }
  for (auto& a : s.getStringDataArrays())
  {
    a.shrink_to_fit();
  }
}

/**
  @brief [ms1-finite] Refuses an input with a non-finite MS1 m/z.

  Nothing after the pick is defined on one, prune on or off:
  sortByPosition leaves finite m/z out of order around a NaN while isSorted() stays true, a NaN band edge leaves two
  band cores owning nothing, and MassTraceDetection's m/z search matches anything.

  @param[in] s The MS1 spectrum
  @throws InputRefusal if an m/z is not finite
*/
inline void requireFiniteMs1Mz(const MSSpectrum& s)
{
  for (Size i = 0; i < s.size(); ++i)
  {
    if (! std::isfinite(s[i].getMZ()))
    {
      throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                         "[ms1] MS1 spectrum has a non-finite m/z (" + std::to_string(s[i].getMZ()) + ") at peak " + std::to_string(i) + " of "
                           + std::to_string(s.size()) + " after peak picking; the input is malformed and the run is refused",
                         s.getNativeID());
    }
  }
}

/**
  @brief [ms1-prune] Picked MS1 centroids at or below trace:noise_threshold_int, dropped at load: MassTraceDetection never
  traces them (its work map keeps (double)intensity > noise only).

  Two things still read them, and both stay exact:
   (i)  the band EDGES are quantiles of a stride-97 m/z sample of the whole map, so the sample is taken here, from the
        sorted spectrum BEFORE a peak is removed, and detectTraces_ uses it;
   (ii) a band spectrum EXISTS iff any peak lies in the band's range, and one whose peaks are all sub-threshold is still
        a SCAN to MassTraceDetection (it counts against min_sample_rate). So a spectrum keeps its first and last peak plus
        the fewest sub-threshold witnesses that leave no two consecutive kept peaks more than hop (relative m/z) apart.
        A band range with its halo is at least 2 hop / (1 - hop) of its lower bound wide, so a range that holds any peak
        holds a kept one, with a factor-2 margin over floating-point rounding.
  Never removes a spectrum or strips the IM array.
*/
struct Ms1Prune
{
  /// @brief Per-spectrum and total prune counts.
  struct Counts
  { long long picked = 0, survivors = 0, kept = 0, at_noise = 0, nan = 0; };
  bool on = false;               ///< perf:ms1_prune
  bool no_witness = false;       ///< DIASPEXTRACT_MS1_PRUNE_NO_WITNESS=1: survivors only -- a negative control, NOT exact
  bool no_chain = false;         ///< DIASPEXTRACT_MS1_PRUNE_NO_CHAIN=1: interior chain links dropped, anchors kept -- a negative control, NOT exact
  double noise = 0.0;            ///< trace:noise_threshold_int, the very double MassTraceDetection is given
  double hop = 0.0;              ///< trace:mass_error_ppm x 1e-6 x BAND_HALO_PPM_MUL
  std::vector<double> sample;    ///< the stride-97 m/z sample of every spectrum pruned so far, taken unpruned
  Size frames = 0, sample_n = 0; ///< spectra pruned; the sample size their picked counts imply
  Counts total;

  /// @brief Back to no spectra pruned: the sample and the counts are cleared.
  void reset()
  {
    std::vector<double>().swap(sample);
    frames = sample_n = 0;
    total = Counts();
  }

  /**
    @brief One spectrum: sort by m/z (a no-op on a picked one), write the m/z values detectTraces_ would
    sample to @p samp, then keep the survivors and the witnesses.

    The chain runs once to count, and a second
    time to collect only when it removes a peak: a spectrum it would not shrink costs no index vector and no select.

    @param[in,out] s The spectrum
    @param[out] samp Receives the sample (every 97th m/z)
    @param[in,out] c Counts
  */
  void frame(MSSpectrum& s, double* samp, Counts& c) const
  {
    s.sortByPosition();
    const Size n = s.size();
    for (Size i = 0; i < n; i += 97)
    {
      *samp++ = s[i].getMZ();
    }
    const Size none = (Size)-1;
    auto reach = [&](Size k, double x) {
      const double m = s[k].getMZ();
      return x <= m + hop * m;
    };
    auto chain = [&](bool count, auto&& take) { // take(i) for every kept peak, in ascending order
      Size last = none, cand = none;            // the last kept peak; the farthest sub-threshold peak within its reach
      for (Size i = 0; i < n; ++i)
      {
        const double v = s[i].getIntensity(); // float -> double, as MassTraceDetection reads it
        const bool surv = v > noise;          // its own test; the complement also drops NaN and == noise
        const double x = s[i].getMZ();
        if (count)
        {
          c.survivors += surv;
          c.at_noise += (v == noise);
          c.nan += (v != v);
        }
        if (no_witness)
        {
          if (surv) { take(i); }
          continue;
        }
        if (last != none && reach(last, x))
        {
          if (surv)
          {
            take(i);
            last = i;
            cand = none;
          }
          else
          {
            cand = i;
          }
          continue;
        }
        if (cand != none) // out of reach: the candidate becomes a witness first
        {
          if (! no_chain)
          {
            take(cand); // (the NO_CHAIN control walks on from it without keeping it)
          }
          last = cand;
          cand = none;
          if (reach(last, x))
          {
            if (surv)
            {
              take(i);
              last = i;
            }
            else
            {
              cand = i;
            }
            continue;
          }
        }
        take(i);
        last = i; // the first peak, or the far side of a gap nothing bridges
      }
      if (cand != none)
      {
        take(cand); // the last peak
      }
    };
    Size kept = 0;
    chain(true, [&kept](Size) { ++kept; });
    c.picked += (long long)n; // [ms1-finite] both load paths refuse a non-finite m/z before the prune
    c.kept += (long long)kept;
    if (kept == n) { return; }
    std::vector<Size> keep;
    keep.reserve(kept);
    chain(false, [&keep](Size i) { keep.push_back(i); });
    s.select(keep);      // peaks and every data array in lockstep
    shrinkDataArrays(s); // select right-sizes only the peaks
  }

  /**
    @brief Prune every spectrum of @p v that @p take selects, in parallel.

    A spectrum's sample size does not
    depend on the prune, so the offsets are fixed first; the order is irrelevant (detectTraces_ sorts).

    @param[in,out] v The spectra
    @param[in] take Selects the spectra to prune
    @par Exceptions
      Rethrows the first exception of a spectrum, of whatever type it was.
  */
  void batch(std::vector<MSSpectrum>& v, bool (*take)(const MSSpectrum&))
  {
    std::vector<char> sel(v.size());
    std::vector<Size> off(v.size());
    Size at = sample.size();
    for (Size i = 0; i < v.size(); ++i)
    {
      sel[i] = take(v[i]);
      off[i] = at;
      if (sel[i]) { at += (v[i].size() + 96) / 97; }
    }
    sample.resize(at);
    std::vector<Counts> cs(v.size());
    std::exception_ptr err;
#pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < (long)v.size(); ++i)
    {
      if (! sel[(Size)i]) { continue; }
      try
      {
        frame(v[(Size)i], sample.data() + off[(Size)i], cs[(Size)i]);
      }
      catch (...)
      {
        taskErr(err);
      }
    }
    if (err) { std::rethrow_exception(err); }
    for (Size i = 0; i < v.size(); ++i)
    {
      if (! sel[i]) { continue; }
      const Counts& c = cs[i];
      ++frames;
      sample_n += (Size)(c.picked + 96) / 97;
      total.picked += c.picked;
      total.survivors += c.survivors;
      total.kept += c.kept;
      total.at_noise += c.at_noise;
      total.nan += c.nan;
    }
  }
};

/**
  @brief [dnoise] dnoise v0.1.0's default MS1 path (DnoiseMs1.h) on every MS1 frame of a Bruker .d, on the raw
  points, before the pick.

  The loader delivers m/z, float32 1/K0 and a corrected intensity, so the
  (TOF bin, scan, raw count) dnoise filters are recovered per point by DIAspeXtract::Dnoise::recover and the run
  is REFUSED on anything that does not recover exactly -- as it is when the loader did not use the
  MzCalibration table model, a spectrum is not a whole raw frame, or an MS1 frame with points does
  not reach the filter exactly once. Kept points are selected in place; an emptied frame stays.
*/
struct DnoiseRun
{
  /// @brief Per-frame and total dnoise counts and filter CPU time.
  struct Counts
  {
    long long frames = 0, emptied = 0, raw = 0, streak = 0, halo = 0, kept = 0, ambiguous = 0;
    double cpu = 0.0;
    uint32_t fid = 0;
  };
  /// @brief one per thread for one batch: dnoise's Scratch plus the recovered integers (~45 B per raw point in flight)
  struct Work
  {
    DIAspeXtract::Dnoise::Scratch w;
    std::vector<uint32_t> scan, tof, raw;
    std::vector<uint8_t> keep;
    std::vector<Size> idx;
  };
  bool on = false; ///< dnoise:ms1 true on a Bruker .d, set up
  DIAspeXtract::Dnoise::Params p;
  std::optional<DIAspeXtract::Dnoise::Ms1WindowGate> gate;
  DIAspeXtract::TdfDnoiseInputs tdf;
  std::vector<std::vector<float>> im_asc; ///< per TimsCalibration row: the float32 1/K0 of each scan its MS1 frames have, ascending
  std::vector<char> im_falls;             ///< the row's 1/K0 falls with the scan, so its im_asc is reversed
  Counts total;
  DIAspeXtract::Dnoise::FrameCoverage seen; ///< the frame ids add() merged, each once

#ifdef _WIN32
  /**
    @brief CPU time of the calling thread (Windows).

    @return Seconds, 0 if unavailable
  */
  static double threadCpu_()
  {
    FILETIME c, e, k, u;
    return GetThreadTimes(GetCurrentThread(), &c, &e, &k, &u) ? spxFileTimeSec(k, u) : 0.0;
  }
#else
  /**
    @brief CPU time of the calling thread.

    @return Seconds, 0 if unavailable
  */
  static double threadCpu_()
  {
    timespec ts;
    return clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) == 0 ? (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9 : 0.0;
  }
#endif

  /**
    @brief Read the tdf, build the gate and the scan tables.

    Empty = ready; anything else is why the run is refused.

    @param[in] tdf_path Path of the analysis.tdf
    @return Empty when ready, else the reason for the refusal
  */
  std::string setup(const std::string& tdf_path)
  {
    std::string why;
    if (! DIAspeXtract::loadTdfDnoise(tdf_path, tdf, why)) { return why; }
    DIAspeXtract::Dnoise::GateBuild g = DIAspeXtract::Dnoise::buildMs1WindowGate(tdf.meta, tdf.windows, p);
    if (! g.error.empty()) { return "the dnoise MS1 port cannot build its window gate from this file: " + g.error; }
    gate = std::move(g.gate);
    for (Size id = 1; id < tdf.frames.size(); ++id) // above 100 ms raw counts are approximate (recoverPoint), not refused
    {
      if (tdf.frames[id].msms_type == 0)
      {
        if (std::string r = DIAspeXtract::Dnoise::Recover::correctionRefusal(id, tdf.frames[id].corr); ! r.empty()) { return r; }
      }
    }
    // one table per TimsCalibration row, over the most scans an MS1 frame using it has (tdf.ms1_max_scans): a row no MS1
    // frame uses stays empty and matches nothing. The gate keeps the run's MAX(NumScans), as dnoise does.
    im_asc.assign(tdf.tims_cal.size(), {});
    im_falls.assign(tdf.tims_cal.size(), 0);
    for (Size ci = 0; ci < tdf.tims_cal.size(); ++ci)
    {
      bool falls = false;
      const std::string r = DIAspeXtract::Dnoise::Recover::buildScanTable(tdf.tims_cal[ci], tdf.ms1_max_scans[ci], im_asc[ci], falls);
      if (! r.empty()) { return r; }
      im_falls[ci] = falls;
    }
    reset();
    on = true;
    return std::string();
  }

  /**
    @brief The scans each row's table covers, for the setup log: "944", or "944,0" with a row no MS1 frame uses.

    @return The scan counts, comma-separated
  */
  String scanTables() const
  {
    String out;
    for (Size ci = 0; ci < tdf.ms1_max_scans.size(); ++ci)
    {
      out += (ci ? "," : "") + String(tdf.ms1_max_scans[ci]);
    }
    return out;
  }

  /**
    @brief The loader must have calibrated m/z with the model the TOF recovery inverts.

    Serially, once per batch.

    @throws InputRefusal if the loader did not use the MzCalibration table model
  */
  void checkLoader() const
  {
    if (BrukerTimsFile::lastMzCalibration() != "tdf_table_modeltype1")
    {
      throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                         "[dnoise] MS1 denoising recovers each point's TOF bin through the MzCalibration "
                         "table model, and the loader calibrated m/z otherwise; -dnoise:ms1 false runs without MS1 denoising",
                         BrukerTimsFile::lastMzCalibration());
    }
  }

  /**
    @brief One MS1 spectrum as the loader delivered it, filtered in place.

    @p k is the calling thread's, @p c the frame's.

    @param[in,out] s The MS1 spectrum
    @param[in,out] k The calling thread's scratch
    @param[out] c The frame's counts
    @throws InputRefusal if a point cannot be recovered exactly or the frame does not match the tdf
  */
  void frame(MSSpectrum& s, Work& k, Counts& c) const
  {
    using DIAspeXtract::Dnoise::Recover::Refusal;
    const double cpu0 = threadCpu_();
    const uint32_t fid = frameIdOf(s.getNativeID());
    const Size n = s.size();
    auto refuse = [&](const String& what) {
      return InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[dnoise] " + what + "; -dnoise:ms1 false runs without MS1 denoising",
                          s.getNativeID());
    };
    const bool has_im = s.containsIMData();
    switch (DIAspeXtract::Dnoise::Recover::frameCheck(tdf.frames, fid, n, has_im, has_im ? s.getFloatDataArrays()[s.getIMData().first].size() : 0))
    {
      case Refusal::NONE:
        break;
      case Refusal::FRAME:
        throw refuse("an MS1 spectrum is not an MS1 frame of the input's analysis.tdf");
      case Refusal::PEAKS:
        throw refuse("frame " + String(fid) + " arrived with " + String(n) + " of its " + String(tdf.frames[fid].num_peaks)
                     + " raw points: only a whole raw frame can be denoised (no frame aggregation, no loader centroiding)");
      default:
        throw refuse("frame " + String(fid) + " carries no 1/K0 per point"); // Refusal::IM
    }
    const DIAspeXtract::TdfDnoiseInputs::Frame& fr = tdf.frames[fid];
    ++c.frames;
    c.raw += (long long)n;
    c.fid = fid;
    if (n == 0)
    {
      c.cpu += threadCpu_() - cpu0;
      return;
    }
    const auto& im = s.getFloatDataArrays()[s.getIMData().first];
    const TofAxis& ax = tofAxis();
    const double b = ax.factor(fid), corr = fr.corr;
    const std::vector<float>& asc = im_asc[fr.cal];
    const bool falls = im_falls[fr.cal] != 0;
    k.scan.resize(n);
    k.tof.resize(n);
    k.raw.resize(n);
    for (Size i = 0; i < n; ++i)
    {
      const float v = s[i].getIntensity();
      double t = 0.0;
      int hits = 0;
      switch (DIAspeXtract::Dnoise::Recover::recoverPoint(ax.cal, b, s[i].getMZ(), im[i], v, corr, asc, falls, fr.num_scans, k.tof[i], k.scan[i],
                                                          k.raw[i], t, hits))
      {
        case Refusal::NONE:
          c.ambiguous += (hits > 1);
          break; // above 100 ms accumulation: the nearest of several raw counts
        case Refusal::TOF:
          throw refuse("m/z " + String(s[i].getMZ()) + " in frame " + String(fid) + " is not the MzCalibration table model's value of a TOF bin ("
                       + String(t) + "): the loader calibrated m/z otherwise");
        case Refusal::SCAN:
          throw refuse("1/K0 " + String(im[i]) + " in frame " + String(fid) + " is not the TimsCalibration model's value of one of its "
                       + String(fr.num_scans) + " scans: the loader converted 1/K0 otherwise");
        case Refusal::INTENSITY_RANGE:
          throw refuse("intensity " + String(v) + " in frame " + String(fid)
                       + " is >= 2^24 = 16777216: the loader's float32 intensity no longer identifies one raw count");
        default: // Refusal::INTENSITY
          throw refuse("intensity " + String(v) + " in frame " + String(fid) + " is not opentims' correction of any raw count");
      }
    }
    DIAspeXtract::Dnoise::FrameCounts fc;
    try
    {
      fc = DIAspeXtract::Dnoise::denoiseMs1Frame(k.scan.data(), k.tof.data(), k.raw.data(), n, fr.num_scans, p, gate ? &*gate : nullptr, k.w, k.keep);
    }
    catch (const std::invalid_argument& e) // the OpenMS-free port throws standard exceptions; the tool reports OpenMS ones
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, e.what(), String(fid));
    }
    c.streak += (long long)fc.after_streak;
    c.halo += (long long)fc.after_halo;
    c.kept += (long long)fc.kept;
    c.emptied += (fc.kept == 0);
    if (fc.kept != n)
    {
      k.idx.clear();
      for (Size i = 0; i < n; ++i)
      {
        if (k.keep[i]) { k.idx.push_back(i); }
      }
      s.select(k.idx);     // the peaks and every data array in lockstep
      shrinkDataArrays(s); // select leaves every data array at its raw capacity
    }
    c.cpu += threadCpu_() - cpu0;
  }

  /**
    @brief A batch's per-spectrum counts, serially: each denoised frame (a slot with frames 0 is an MS2 spectrum's) marked
    once -- a frame delivered twice is refused, it would be filtered and traced twice.

    @param[in] v Counts per spectrum
    @throws InputRefusal if a frame was delivered twice
  */
  void add(const std::vector<Counts>& v)
  {
    for (const Counts& c : v)
    {
      if (c.frames == 0) { continue; }
      if (! seen.mark(c.fid))
      {
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           "[dnoise] MS1 frame " + String(c.fid)
                             + " reached the filter twice: the loader delivered it more than once; -dnoise:ms1 false runs without MS1 denoising",
                           "frame=" + String(c.fid));
      }
      total.frames += c.frames;
      total.emptied += c.emptied;
      total.raw += c.raw;
      total.streak += c.streak;
      total.halo += c.halo;
      total.kept += c.kept;
      total.ambiguous += c.ambiguous;
      total.cpu += c.cpu;
    }
  }

  /**
    @brief After the load, serially: every MS1 frame with points (Frames.NumPeaks > 0) reached the filter.

    Counts cannot
    tell -- a frame delivered twice and a missing one with equal NumPeaks cancel out -- so this reads the marks.

    @throws InputRefusal if an MS1 frame with points never reached the filter
  */
  void requireComplete() const
  {
    if (const Size miss = seen.firstMissing(tdf.frames))
    {
      throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                         "[dnoise] MS1 frame " + String(miss) + " (" + String(tdf.frames[miss].num_peaks)
                           + " raw points in Frames.NumPeaks) never reached the filter: the loader did not deliver every MS1 frame; -dnoise:ms1 "
                             "false runs without MS1 denoising",
                         "frame=" + String(miss));
    }
  }

  /**
    @brief The counts and marks back to none: at setup, and when a failed stream load falls back to the resident one, which
    delivers every frame again.
  */
  void reset()
  {
    total = Counts();
    seen.reset();
  }

  /**
    @brief The applied parameters and the run's gate size, for the log and spx:dnoise_ms1_params.

    @return The parameter text
  */
  String params() const
  {
    return "mz_half_width=" + String(p.mz_half_width) + " min_feature_length=" + String(p.min_feature_length)
           + " max_internal_gap=" + String(p.max_internal_gap) + " iterations=" + String(p.iterations)
           + " min_window_intensity=" + String(p.min_window_intensity) + " min_feature_intensity=" + String(p.min_feature_intensity)
           + " halo=" + (p.halo ? "true" : "false") + " halo_peak_fraction=" + String(p.halo_peak_fraction)
           + " halo_mz_idx_half_width=" + String(p.halo_mz_idx_half_width) + " halo_scan_half_width=" + String(p.halo_scan_half_width)
           + " dia_ms1_window=" + (p.dia_ms1_window ? "true" : "false") + " dia_ms1_mz_pad=" + String(p.dia_ms1_mz_pad)
           + " dia_ms1_im_pad=" + String(p.dia_ms1_im_pad) + " gate_boxes=" + String(gate ? gate->boxes().size() : (Size)0);
  }

  /**
    @brief The run's dnoise counts as one log line.

    @return The summary text
  */
  String summary() const
  {
    char b[420];
    snprintf(
      b, sizeof b,
      "[dnoise] %s MS1: %lld frames (%lld emptied), %lld raw points -> streak %lld -> halo %lld -> window gate: kept %lld MS1 points (%.2f%%); "
      "%lld raw counts approximate (accumulation above 100 ms: the nearest of several counts); filter CPU %.1f s",
      DIAspeXtract::Dnoise::PORTED_FROM, total.frames, total.emptied, total.raw, total.streak, total.halo, total.kept,
      total.raw > 0 ? 100.0 * (double)total.kept / (double)total.raw : 0.0, total.ambiguous, total.cpu);
    return String(b);
  }
};

/**
  @brief [stream] Pick-and-compact consumer: each frame is picked and compacted on arrival, so the raw run
  is never resident.

  Windows are keyed from the spectrum's own precursor, not swath_nr.
*/
class PickCompactConsumer : public FullSwathFileConsumer
{
public:
  /**
    @brief Construct the consumer over the caller's picker, window slabs, MS1 map and stats.

    @param[in] picker The ion-mobility peak picker (copied per thread)
    @param[in,out] windows Receives one slab per isolation window
    @param[in,out] ms1 Receives the picked MS1 spectra
    @param[in,out] stats Receives the compaction counts
  */
  PickCompactConsumer(PeakPickerIM& picker, map<WinKey, PeakSlab>& windows, PeakMap& ms1, CompactStats& stats):
      picker_(picker),
      windows_(windows),
      ms1_(ms1),
      stats_(stats)
  {
  }
  /// @brief [overlap] A reader that failed mid-tile leaves a pick job queued: wait for it, its helpers use pend_.
  ~PickCompactConsumer()
  {
    for (PickBatch* b : {pend_.get(), next_.get()})
    {
      if (b->job)
      {
        try
        {
          teamForFinish(b->job);
        }
        catch (...)
        {
        }
      }
    }
  }

  Size frames_seen = 0;

  /**
    @brief [frames] A window's FROZEN frame table (rt, vendor id) with the lockstep cursor: every table
    frame the reader skipped (no raw peaks, an empty mobility slice) is appended EMPTY with the
    table's rt and id, so the slab's frame sequence IS the table's in every mode (an empty frame
    is a scan, not a miss).
  */
  struct FrameTab
  {
    std::vector<double> rt;
    std::vector<uint32_t> id;
    Size cur = 0, delivered = 0, padded = 0;
    /**
      @brief Pad the slab up to the delivered frame and advance the cursor past it.

      Pad up to the delivered frame (its rt AND vendor id equal the table row's -- an rt alone is
      ambiguous when two rows share one, e.g. an mzPeak archive without retention times); false
      = not in the table

      @param[in,out] sl The window's slab
      @param[in] t Retention time of the delivered frame
      @param[in] fid Vendor id of the delivered frame
      @return false if the frame is not in the table
    */
    bool step(PeakSlab& sl, double t, uint32_t fid)
    {
      while (cur < rt.size() && (rt[cur] != t || id[cur] != fid))
      {
        sl.append(rt[cur], id[cur], nullptr, nullptr, nullptr, 0, rt.size());
        ++cur;
        ++padded;
      }
      if (cur == rt.size()) { return false; }
      ++cur;
      ++delivered;
      return true;
    }
    /**
      @brief Append every remaining table frame to the slab as an empty frame.

      @param[in,out] sl The window's slab
    */
    void tail(PeakSlab& sl)
    {
      for (; cur < rt.size(); ++cur, ++padded)
      {
        sl.append(rt[cur], id[cur], nullptr, nullptr, nullptr, 0, rt.size());
      }
    }
  };
  using FrameTables = std::map<WinKey, FrameTab>;
  /// The reader's frozen tables, fetched at the first MS2 flush: both streaming readers publish
  /// their metadata before their first MS2 frame (the .d loader exports its frame table and
  /// windows before its MS2 loop, the mzPeak reader sweeps the archive before decoding). Unset
  /// (an mzML) = no tables, no padding: the slab is its own table.
  std::function<FrameTables()> tablefn;
  /**
    @brief The reader's frozen frame tables (empty for an mzML).

    @return The tables
  */
  const FrameTables& tables() const
  { return tables_; }
  /// [ms1-prune] pass 1's run-level MS1 prune and its edge sample; null = every picked peak is kept
  Ms1Prune* ms1_prune = nullptr;
  /// [dnoise] the run's MS1 denoising of the raw frames before their pick; null or off = none
  DnoiseRun* dnoise = nullptr;

  /**
    @brief [perf-load] The readers hand frames over serially: pick+compact a buffered BATCH in parallel, then
    append strictly by buffer index (arrival order).

    One pick batch: the frames (moved out of buf_), their windows, and everything the pick writes. Two of them
    alternate on the pipelined reader, so one batch is picked while the reader decodes the next.
  */
  struct PickBatch
  {
    std::vector<MapType::SpectrumType> buf;
    std::vector<WinKey> key;
    std::vector<Size> off; ///< per-frame regions of the scratch (raw sizes)
    std::vector<uint32_t> mzq;
    std::vector<float> inten;
    std::vector<uint16_t> imq;
    std::vector<Size> kept;
    std::vector<double> frt;
    std::vector<uint32_t> fid;
    std::vector<CompactStats> st;                  ///< per-frame stats, merged serially at the append
    std::vector<std::unique_ptr<PeakPickerIM>> pk; ///< [pipe] one picker copy per puller slot
    std::exception_ptr err;
    std::function<void(long, int)> body; ///< the pick loop body the job calls (pipelined)
    std::shared_ptr<PullJob> job;        ///< the live pick job, null when none
  };
  /**
    @brief Take buf_/key_ into b and size its outputs (the scratch keeps its capacity between batches).

    @param[out] b The batch
  */
  void takeBatch_(PickBatch& b)
  {
    b.buf.swap(buf_);
    b.key.swap(key_);
    buf_.clear();
    key_.clear();
    const Size n = b.buf.size();
    b.off.resize(n + 1);
    Size total = 0;
    for (Size i = 0; i < n; ++i)
    {
      b.off[i] = total;
      total += b.buf[i].size();
    }
    b.off[n] = total;
    if (b.mzq.size() < total)
    {
      b.mzq.resize(total);
      b.inten.resize(total);
      b.imq.resize(total);
    }
    b.kept.assign(n, 0);
    b.frt.assign(n, 0.0);
    b.fid.assign(n, 0);
    b.st.assign(n, CompactStats());
    b.err = nullptr;
  }
  /**
    @brief Pick and compact frame @p i of the batch; an exception is stored in the batch.

    @param[in,out] b The batch
    @param[in] i Frame index in the batch
    @param[in,out] pk The calling thread's picker
  */
  void pickOne_(PickBatch& b, const Size i, PeakPickerIM& pk)
  {
    try
    {
      auto& s = b.buf[i];
      ensureIMArrayName(s);
      if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) { pk.pickIMCluster(s); }
      ensureIMArrayName(s);
      b.frt[i] = s.getRT();
      b.fid[i] = frameIdOf(s.getNativeID());
      if (! frameCalibrated(b.fid[i])) { ++b.st[i].unmapped_frame; }
      const Size o = b.off[i];
      b.kept[i] = compactifyInto(s, b.st[i], &b.mzq[o], &b.inten[o], &b.imq[o], b.off[i + 1] - o);
      s.clear(true); // raw frame dies here, inside the loop
    }
    catch (...)
    {
      taskErr(b.err);
    }
  }
  /**
    @brief [pipe][overlap] Queue the batch's pick on the window team and return: the reader decodes the next frames while
    free threads (and windows at their break points) pick these.

    @param[in,out] b The batch
  */
  void startPick_(PickBatch& b)
  {
    b.pk.resize((Size)std::max(1, omp_get_num_threads()));
    b.body = [this, &b](long i, int slot) {
      std::unique_ptr<PeakPickerIM>& p = b.pk[(Size)slot];
      if (! p)
      {
        p.reset(new PeakPickerIM(picker_)); // Param is a deep-copy value type
      }
      pickOne_(b, (Size)i, *p);
    };
    b.job = teamForStart((long)b.buf.size(), 1, b.body);
  }
  /**
    @brief Wait for the batch's queued pick to finish.

    @param[in,out] b The batch
    @par Exceptions
      Rethrows the first exception of the pick, of whatever type it was.
  */
  void finishPick_(PickBatch& b)
  {
    if (! b.job) { return; }
    std::shared_ptr<PullJob> j;
    j.swap(b.job);
    teamForFinish(j);
    if (b.err) { std::rethrow_exception(b.err); }
  }
  /**
    @brief Append in arrival order per window, the windows of this batch in parallel: a window's
    frames stay in order (one thread walks them), distinct windows never touch each other's
    slab, and the copy of ~1 GB per batch into cold pages runs on many threads.

    @param[in,out] b The batch; its frames are cleared
    @throws InputRefusal if a delivered MS2 window is not in the reader's frame table
  */
  void appendBatch_(PickBatch& b)
  {
    const int n = (int)b.buf.size();
    struct Batch
    {
      PeakSlab* slab;
      FrameTab* tab;
      std::vector<int> idx;
    };
    std::vector<Batch> per_win;
    {
      std::map<WinKey, Size> at;
      for (int i = 0; i < n; ++i)
      {
        stats_ += b.st[(Size)i];
        const WinKey& k = b.key[(Size)i];
        auto it = at.find(k);
        if (it == at.end())
        {
          it = at.emplace(k, per_win.size()).first;
          auto tt = tables_.find(k);
          if (! tables_.empty() && tt == tables_.end())
          {
            throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[frames] a delivered MS2 window is not in the reader's frame table",
                               winLabel(k) + " group " + String(k[2]));
          }
          per_win.push_back(Batch {&windows_[k], tt == tables_.end() ? nullptr : &tt->second, {}});
        }
        per_win[it->second].idx.push_back(i);
      }
    }
    const Size n_ms1 = ms1_.size();
    std::exception_ptr pad_err;
    auto append_win = [&](const long w) {
      Batch& bw = per_win[(Size)w];
      for (int i : bw.idx)
      {
        const Size exp_f = bw.tab ? bw.tab->rt.size() : n_ms1;                  // the window's frame count: exact when a table exists
        if (bw.tab && ! bw.tab->step(*bw.slab, b.frt[(Size)i], b.fid[(Size)i])) // [frames] lockstep with the frozen table
        {
#pragma omp critical(pad_err)
          if (! pad_err)
          {
            pad_err = std::make_exception_ptr(
              InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           "[frames] a delivered MS2 frame (rt, frame id) is not in its window's frozen frame table (window "
                             + winLabel(b.key[(Size)i]) + ", frame " + String(b.fid[(Size)i]) + ")",
                           String(b.frt[(Size)i], 6)));
          }
          break;
        }
        const Size o = b.off[(Size)i];
        bw.slab->append(b.frt[(Size)i], b.fid[(Size)i], &b.mzq[o], &b.inten[o], &b.imq[o], b.kept[(Size)i], exp_f);
      }
    };
    if (! omp_in_parallel())
    {
#pragma omp parallel for schedule(dynamic, 1)
      for (long w = 0; w < (long)per_win.size(); ++w)
      {
        try
        {
          append_win(w);
        }
        catch (...)
        {
          taskErr(pad_err);
        }
      } // [task-err] bad_alloc, the slab bounds
    }
    else
    {
      teamFor((long)per_win.size(), 1, [&](long w, int) { append_win(w); }); // [pipe] as the pick
    }
    if (pad_err) { std::rethrow_exception(pad_err); }
    b.buf.clear();
    b.key.clear();
  }
  /// @brief [overlap] The pending batch (picked or being picked, not yet appended): finish its pick, append it.
  void drainPending_()
  {
    if (! pending_) { return; }
    pending_ = false;
    finishPick_(*pend_);
    appendBatch_(*pend_);
  }

  /**
    @brief Pick, compact and append the buffered MS2 frames.

    [perf-load] The readers hand frames over serially: pick+compact a buffered BATCH in parallel, then
    append strictly by buffer index (arrival order).
    [pipe][overlap] On the pipelined tile reader the pick of a batch runs while the reader decodes the next one:
    this call queues the new batch's pick, then finishes and appends the PREVIOUS batch (still in arrival order,
    before the new one), and returns. finish() appends the last one. The output is the serial order's.

    @par Exceptions
      Rethrows the first exception of the pick, of whatever type it was.
  */
  void flush_()
  {
    if (! tables_init_)
    {
      tables_init_ = true;
      if (tablefn) { tables_ = tablefn(); }
    }
    if (buf_.empty()) { return; }
    if (! omp_in_parallel())
    {
      drainPending_();
      PickBatch& b = *next_;
      takeBatch_(b);
      PeakPickerIM picker = picker_; // Param is a deep-copy value type
      const int n = (int)b.buf.size();
// chunk 1: one frame per chunk, so a batch spreads over every worker
#pragma omp parallel for firstprivate(picker) schedule(dynamic, 1)
      for (int i = 0; i < n; ++i)
      {
        pickOne_(b, (Size)i, picker);
      }
      if (b.err) { std::rethrow_exception(b.err); }
      appendBatch_(b);
      return;
    }
    PickBatch& b = *next_; // free: its last append is done
    takeBatch_(b);
    startPick_(b);
    drainPending_(); // the previous batch, while the helpers pick this one
    std::swap(pend_, next_);
    pending_ = true;
  }

protected:
  /**
    @brief Buffer one MS1 spectrum; a full buffer is picked and appended to the MS1 map.

    @param[in,out] s The spectrum (moved from)
    @throws InputRefusal if dnoise is on and the loader did not use the MzCalibration table model
  */
  void consumeMS1Spectrum_(MapType::SpectrumType& s) override
  {
    // MS1 stays a PeakMap (MassTraceDetection needs it whole); buffered and picked in parallel like MS2, order kept by index.
    // [dnoise] the loader chose its m/z model before its first frame: refuse at that frame, not after a whole batch
    // (pickMS1Batch_ checks every batch again and stays the authority)
    if (dnoise && dnoise->on && ! dn_loader_checked_)
    {
      dnoise->checkLoader();
      dn_loader_checked_ = true;
    }
    ms1_buf_.push_back(std::move(s));
    ++frames_seen;
    if (ms1_buf_.size() >= BATCH_SIZE) { pickMS1Batch_(ms1_buf_); }
  }

  /**
    @brief One MS1 batch, picked and appended to the MS1 map, then emptied: the order of the MS1 map is the arrival order.

    @param[in,out] buf The MS1 frames; emptied
    @par Exceptions
      Rethrows the first exception of the pick, of whatever type it was.
  */
  void pickMS1Batch_(std::vector<MapType::SpectrumType>& buf)
  {
    const int n = (int)buf.size();
    if (n == 0) { return; }
    PeakPickerIM picker = picker_;
    std::exception_ptr pick_err;
    // [dnoise] per thread for this batch only (freed with it), per frame counts merged below
    const bool dn = dnoise && dnoise->on;
    if (dn) { dnoise->checkLoader(); }
    std::vector<DnoiseRun::Work> dn_work(dn ? (Size)std::max(1, omp_get_max_threads()) : 0);
    std::vector<DnoiseRun::Counts> dn_cnt(dn ? (Size)n : 0);
#pragma omp parallel for firstprivate(picker) schedule(dynamic, 1)
    for (int i = 0; i < n; ++i)
    {
      try
      {
        auto& s = buf[(Size)i];
        ensureIMArrayName(s);
        if (dn)
        {
          dnoise->frame(s, dn_work.at((Size)omp_get_thread_num()), dn_cnt[(Size)i]); // [dnoise] the raw points, before the pick
        }
        if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) { picker.pickIMCluster(s); }
        ensureIMArrayName(s);
        requireFiniteMs1Mz(s); // [ms1-finite] before the prune, prune on or off
      }
      catch (...)
      {
        taskErr(pick_err);
      }
    }
    if (pick_err) { std::rethrow_exception(pick_err); }
    if (dn) { dnoise->add(dn_cnt); }
    if (ms1_prune && ms1_prune->on)
    {
      ms1_prune->batch(buf, [](const MSSpectrum&) { return true; }); // [ms1-prune]
    }
    for (auto& s : buf)
    {
      ms1_.addSpectrum(std::move(s));
    }
    buf.clear();
  }

  /**
    @brief Buffer one MS2 spectrum under its window key (from its own precursor); a full buffer is flushed.

    @param[in,out] s The spectrum (moved from); one without a precursor is ignored
  */
  void consumeSwathSpectrum_(MapType::SpectrumType& s, Size /*swath_nr*/) override
  {
    const auto& prec = s.getPrecursors();
    if (prec.empty())
    {
      return; // no isolation window -> not routable
    }
    const double c = prec[0].getMZ();
    const double lo = c - prec[0].getIsolationWindowLowerOffset();
    const double hi = c + prec[0].getIsolationWindowUpperOffset();
    pickMS1Batch_(ms1_buf_); // MS1 arrives first; pick it before any MS2 work
    key_.push_back(winKey(lo, hi, windowGroupOf(s.getNativeID())));
    buf_.push_back(std::move(s)); // pick+compact deferred to flush_()
    ++frames_seen;
    if (buf_.size() >= BATCH_SIZE) { flush_(); }
  }

  /// @brief Calls finish().
  void ensureMapsAreFilled_() override
  { finish(); }

public:
  /**
    @brief Flush the buffered tails.

    Callers MUST call this after the load: FullSwathFileConsumer reaches
    ensureMapsAreFilled_ only via retrieveSwathMaps, which this consumer never uses, so without it the
    last `n mod BATCH_SIZE` spectra are dropped.
  */
  void finish()
  {
    pickMS1Batch_(ms1_buf_);
    flush_();
    drainPending_();
    // [frames] the tail of every delivered window's table; a table window that delivered nothing
    // gets no slab (it is counted as absent by the caller)
    for (auto& kv : tables_)
    {
      auto w = windows_.find(kv.first);
      if (w != windows_.end()) { kv.second.tail(w->second); }
    }
  }

private:
  /// Frames per pick batch (~1 GB raw buffer), MS1 and MS2. Batch composition does not change the output.
  static constexpr Size BATCH_SIZE = 256;
  std::vector<MapType::SpectrumType> buf_;
  std::vector<WinKey> key_;
  FrameTables tables_;
  bool tables_init_ = false;                                               ///< [frames] the reader's frozen tables + cursors
  std::unique_ptr<PickBatch> next_ {new PickBatch}, pend_ {new PickBatch}; ///< [overlap] the free batch, the pending one
  bool pending_ = false;                                                   ///< [overlap] pend_ holds a batch whose append is due
  PeakPickerIM& picker_;
  std::vector<MapType::SpectrumType> ms1_buf_; // MS1 frames awaiting the parallel pick
  bool dn_loader_checked_ = false;             ///< [dnoise] checkLoader ran on the first MS1 frame
  map<WinKey, PeakSlab>& windows_;
  PeakMap& ms1_;
  CompactStats& stats_;
};

/**
  @brief What the scorer needs per window: the Pearson denominator G (the frozen table's frame count),
  per-fragment mean and 1/norm over G, and the MS1 -> window nearest-frame map.
*/
struct FragStats
{
  vector<double> mean, invnorm;               ///< per fragment, over the full grid of grid_frames frames
  const vector<int>* nearest_local = nullptr; ///< MS1-store frame -> nearest window frame within delta_rt, or -1 (the FrozenWin's)
  Size grid_frames = 0;                       ///< frames in the window's grid
};

/**
  @brief The Pearson support is the window's FROZEN frame table: G = F (every frame), and the
  MS1 -> window map is the table's, so neither depends on which frames are resident.

  @param[in] frags The window's fragment traces
  @param[in] wst Their store
  @param[in] fw The window's frozen frame table
  @return The statistics
*/
FragStats buildFragStats(const vector<Trace>& frags, const TraceStore& wst, const FrozenWin& fw)
{
  FragStats g;
  g.grid_frames = fw.frames();
  const double n_grid = (double)g.grid_frames;
  g.nearest_local = &fw.nearest_local;
  g.mean.resize(frags.size());
  g.invnorm.resize(frags.size());
#pragma omp taskloop default(shared)
  for (long ii = 0; ii < (long)frags.size(); ++ii)
  {
    const Trace& f = frags[(Size)ii];
    double sum = 0, sumsq = 0;
    for (Size k = 0; k < f.span(); ++k)
    {
      const double v = xv(wst, f, k);
      if (v > 0.0)
      {
        sum += v;
        sumsq += v * v;
      }
    }
    g.mean[ii] = sum / n_grid;
    const double var = sumsq - n_grid * g.mean[ii] * g.mean[ii];
    g.invnorm[ii] = var > 0 ? 1.0 / sqrt(var) : 0.0; // 0 => constant/degenerate -> never correlates
  }
  return g;
}

/**
  @brief Canonical total-order sort (RT, precursor m/z, charge, then the fragment sequence):
  thread-count-independent and free of RT-tie ambiguity (RT-only sortSpectra is NOT canonical).

  Sorts a compact key array and permutes once: std::sort moves what it is handed, and an
  MSSpectrum is ~750 bytes against 32 for a key. The keys start in the spectra's order, so
  introsort makes the same decisions as sorting the spectra would, ties included.
  Tiles are RT-disjoint and RT is the first key, so per-tile sorted blocks concatenate to the
  order one sort of everything gives.

  @param[in,out] v The spectra
*/
static void canonicalSort(vector<MSSpectrum>& v)
{
  struct SortKey
  {
    double rt, mz;
    int charge;
    uint32_t n, idx;
  };
  vector<SortKey> keys(v.size());
  auto key = [&](const long i) {
    const MSSpectrum& a = v[i];
    keys[i].rt = a.getRT();
    keys[i].mz = a.getPrecursors().empty() ? 0.0 : a.getPrecursors()[0].getMZ();
    keys[i].charge = a.getPrecursors().empty() ? 0 : a.getPrecursors()[0].getCharge();
    keys[i].n = (uint32_t)a.size();
    keys[i].idx = (uint32_t)i;
  };
  if (! omp_in_parallel())
  {
#pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)v.size(); ++i)
    {
      key(i);
    }
  }
  else
  {
    teamFor((long)v.size(), 4096, [&](long i, int) { key(i); }); // [pipe] the pipelined tile writer
  }
  sort(keys.begin(), keys.end(), [&v](const SortKey& ka, const SortKey& kb) {
    if (ka.rt != kb.rt) { return ka.rt < kb.rt; }
    if (ka.mz != kb.mz) { return ka.mz < kb.mz; }
    if (ka.charge != kb.charge) { return ka.charge < kb.charge; }
    // fragment-sequence tiebreaks so equal-precursor spectra have a canonical order
    if (ka.n != kb.n) { return ka.n < kb.n; }
    const MSSpectrum& a = v[ka.idx];
    const MSSpectrum& b = v[kb.idx];
    for (Size i = 0; i < a.size(); ++i)
    {
      if (a[i].getMZ() != b[i].getMZ()) { return a[i].getMZ() < b[i].getMZ(); }
      if (a[i].getIntensity() != b[i].getIntensity()) { return a[i].getIntensity() < b[i].getIntensity(); }
    }
    // equal peaks: the fields that still serialise differently -- two overlapping isolation
    // windows can emit one precursor with the same peaks and different offsets
    const auto& pa = a.getPrecursors();
    const auto& pb = b.getPrecursors();
    if (pa.empty() != pb.empty()) { return pa.empty(); }
    if (! pa.empty())
    {
      if (pa[0].getIsolationWindowLowerOffset() != pb[0].getIsolationWindowLowerOffset())
      {
        return pa[0].getIsolationWindowLowerOffset() < pb[0].getIsolationWindowLowerOffset();
      }
      if (pa[0].getIsolationWindowUpperOffset() != pb[0].getIsolationWindowUpperOffset())
      {
        return pa[0].getIsolationWindowUpperOffset() < pb[0].getIsolationWindowUpperOffset();
      }
      if (pa[0].getDriftTime() != pb[0].getDriftTime()) { return pa[0].getDriftTime() < pb[0].getDriftTime(); }
    }
    for (const char* key : {"spx_guessed", "spx_n_isotopes"})
    {
      const int ia = a.metaValueExists(key) ? (int)a.getMetaValue(key) : -1, ib = b.metaValueExists(key) ? (int)b.getMetaValue(key) : -1;
      if (ia != ib) { return ia < ib; }
    }
    return false;
  });
  vector<MSSpectrum> ord;
  ord.reserve(v.size());
  for (const SortKey& k : keys)
  {
    ord.push_back(std::move(v[k.idx]));
  }
  v.swap(ord);
}

/**
  @brief [twins] One precursor claimed twice is emitted as two spectra with complementary fragments: two seeds that
  resolved to the same monoisotope (the claim loop's `used[]` owns traces, not mono positions), or one
  precursor in the 1-m/z overlap of two isolation windows.

  Folds each such cluster into one spectrum.
  Relation: same MS1 frame (RT bit-equal; every precursor RT is a frame time), same charge > 0, precursor
  m/z within 20 ppm of the member earlier in `v`, |d1/K0| <= im_tol; clusters are its connected
  components (a chain of links can span more than the tolerances). Survivor: most isotopes, then most peaks, then the
  earliest. Its peaks become the union of the cluster's, fused within 10 ppm (m/z intensity-weighted, intensity = max); every other field is its
  own. Same frame => same RT cell => same tile under every grouping, so the output stays grouping- and
  thread-invariant.
  `v` in canonical order; survivors keep their places. Returns the number of spectra removed.
  [twin-rt] `rt_tol` > 0 widens "same MS1 frame" to |dRT| <= rt_tol (assembly:twin_rt_frames): two calls of one
  monoisotope at nearby frames (seeds of one ion whose apexes differ, or v2's extension onto another precursor's
  monoisotope) fold like same-frame twins. Such a pair must also lie in one RT cell: the cell that owns a precursor
  (`cells`: the cell edges, `cell_drt`: tile_drt; the tile loop's own test, (T - rt) > cell_drt), so a pair never
  straddles a tile under any grouping and the output stays grouping- and thread-invariant. rt_tol 0 is the
  same-frame relation (|dRT| <= 0 and one frame is one cell), bit for bit. `st` (optional) adds {pairs linked at one
  frame, pairs linked across frames, near pairs refused at a cell edge}. A pair from two frames must have its precursor
  m/z within BIN_PPM (one flight-time bin: one MS1 trace, or one bin at two frames, is the monoisotope of both), and
  one of the two must carry the userParam spx_mono_ext (its mono was placed by monocheck v2's extension); the mark
  is removed from every spectrum.

  @param[in,out] v The spectra, in canonical order
  @param[in] im_tol Maximum |d 1/K0|
  @param[in] rt_tol Maximum |dRT| (0: same frame)
  @param[in] cells RT cell edges, or nullptr
  @param[in] cell_drt The tile loop's cell drt
  @param[out] st Optional link counts
  @return Number of spectra removed
  @par Exceptions
    Rethrows the first exception of a parallel part, of whatever type it was.
*/
static Size mergeTwins(vector<MSSpectrum>& v,
                       double im_tol,
                       double rt_tol = 0.0,
                       const vector<double>* cells = nullptr,
                       double cell_drt = 0.0,
                       std::array<Size, 3>* st = nullptr)
{
  const double TWIN_PPM = 20.0, FUSE_PPM = 10.0; // fixed tolerances
  // [twin-rt] one flight-time bin across frames: the per-frame calibration moves one bin's m/z by < 0.01 ppm, distinct
  // monoisotope traces of one charge and 1/K0 at nearby frames sit >= 0.3 ppm apart
  const double BIN_PPM = 0.05;
  struct K
  {
    double rt, mz, im;
    int z;
    uint32_t i, cell;
    bool ext;
  };
  auto cell_of = [&](double rt) -> uint32_t { // the first edge T with (T - rt) > cell_drt; edges ascend, so the test is monotone
    if (! cells || cells->empty()) { return 0; }
    return (uint32_t)(std::partition_point(cells->begin(), cells->end(), [&](double cell_rt) { return ! (cell_rt - rt > cell_drt); })
                      - cells->begin());
  };
  vector<K> ks;
  ks.reserve(v.size());
  for (Size i = 0; i < v.size(); ++i)
  {
    const auto& p = v[i].getPrecursors();
    if (! p.empty() && p[0].getCharge() > 0)
    {
      ks.push_back({v[i].getRT(), p[0].getMZ(), p[0].getDriftTime(), p[0].getCharge(), (uint32_t)i, cell_of(v[i].getRT()),
                    v[i].metaValueExists("spx_mono_ext")});
    }
  }
  // one charge, then m/z: a pair is found from its lower-m/z member's 20-ppm window whatever its frames
  sort(ks.begin(), ks.end(), [](const K& a, const K& b) {
    if (a.z != b.z) { return a.z < b.z; }
    if (a.mz != b.mz) { return a.mz < b.mz; }
    if (a.rt != b.rt) { return a.rt < b.rt; }
    return a.i < b.i;
  });
  // the twin relation, a earlier in v than b, its m/z window a's: 0 none, 1 one frame, 2 nearby frames, 3 refused at a cell edge
  auto linked = [&](const K& a, const K& b) -> int {
    if (! (fabs(b.mz - a.mz) <= a.mz * TWIN_PPM * 1e-6) || ! (fabs(b.im - a.im) <= im_tol) || ! (fabs(b.rt - a.rt) <= rt_tol)) { return 0; }
    const bool same = b.rt == a.rt;
    if (! same && fabs(b.mz - a.mz) > a.mz * BIN_PPM * 1e-6) { return 0; }
    if (! same && ! a.ext && ! b.ext) { return 0; }
    if (b.cell != a.cell) { return 3; }
    return same ? 1 : 2;
  };
  vector<uint32_t> up(v.size());
  for (Size i = 0; i < up.size(); ++i)
  {
    up[i] = (uint32_t)i;
  }
  auto root = [&](uint32_t a) {
    while (up[a] != a)
    {
      up[a] = up[up[a]];
      a = up[a];
    }
    return a;
  };
  Size n_same = 0, n_near = 0, n_edge = 0;
  for (Size g0 = 0, g1; g0 < ks.size(); g0 = g1)
  {
    for (g1 = g0 + 1; g1 < ks.size() && ks[g1].z == ks[g0].z; ++g1) {} // one charge
    for (Size a = g0; a < g1; ++a)
    {
      // the probe is the EARLIER member of a pair, with its own tolerance
      const double tol = ks[a].mz * TWIN_PPM * 1e-6, lo = ks[a].mz - tol, hi = ks[a].mz + tol;
      auto link = [&](Size b) {
        if (ks[b].i <= ks[a].i) { return; }
        const int l = linked(ks[a], ks[b]);
        if (l == 3) { ++n_edge; }
        if (l != 1 && l != 2) { return; }
        ++(l == 1 ? n_same : n_near);
        {
          const uint32_t ra = root(ks[a].i), rb = root(ks[b].i);
          if (ra != rb) { up[ra] = rb; }
        }
      };
      for (Size b = a + 1; b < g1 && ks[b].mz <= hi; ++b)
      {
        link(b);
      }
      for (Size b = a; b-- > g0 && ks[b].mz >= lo;)
      {
        link(b);
      }
    }
  }
  if (st)
  {
    (*st)[0] += n_same;
    (*st)[1] += n_near;
    (*st)[2] += n_edge;
  }
  // survivor per cluster: (isotopes, peaks) strictly greater replaces, so a tie keeps the earliest
  const uint32_t NONE = UINT32_MAX;
  vector<uint32_t> n_in(v.size(), 0), best(v.size(), NONE);
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    ++n_in[root(i)];
  }
  auto isotopes = [&v](uint32_t i) { return v[i].metaValueExists("spx_n_isotopes") ? (int)v[i].getMetaValue("spx_n_isotopes") : 0; };
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    const uint32_t r = root(i);
    if (n_in[r] < 2) { continue; }
    const uint32_t b = best[r];
    if (b == NONE || std::make_pair(isotopes(i), v[i].size()) > std::make_pair(isotopes(b), v[b].size())) { best[r] = i; }
  }
  // A component is a chain of pairwise links and can span more than the tolerances (500.0000, 500.0095 and 500.0190
  // link pairwise within 20 ppm, the ends 38 ppm apart). The whole component folds: on real runs such chains are one
  // ion called at drifting m/z; folding only the survivor's direct links wrote 10.7 % more spectra on a public
  // cell-line run for the same peptides (Sage +0.02 %, MSFragger -0.19 %).
  vector<uint32_t> grp(v.size(), NONE); // the survivor a spectrum folds into (itself for a survivor)
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    if (n_in[root(i)] >= 2) { grp[i] = best[root(i)]; }
  }
  vector<uint32_t> surv;
  vector<int> slot(v.size(), -1);
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    if (grp[i] == i)
    {
      slot[i] = (int)surv.size();
      surv.push_back(i);
    }
  }
  vector<vector<uint32_t>> members(surv.size());
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    if (grp[i] != NONE) { members[(Size)slot[grp[i]]].push_back(i); }
  }
  // fuse: each thread reads and writes one cluster's spectra only
  auto fuse = [&](const long c) {
    vector<pair<double, double>> ps;
    for (uint32_t i : members[(Size)c])
    {
      for (const auto& pk : v[i])
      {
        ps.emplace_back(pk.getMZ(), (double)pk.getIntensity());
      }
    }
    sort(ps.begin(), ps.end());
    vector<pair<double, double>> fused;
    double gm = 0.0, gw = 0.0, gi = 0.0;
    for (const auto& pk : ps)
    {
      if (gw != 0.0 && pk.first - gm / gw > (gm / gw) * (FUSE_PPM * 1e-6))
      {
        fused.emplace_back(gm / gw, gi);
        gm = gw = gi = 0.0;
      }
      gm += pk.first * pk.second;
      gw += pk.second;
      gi = std::max(gi, pk.second);
    }
    if (gw != 0.0) { fused.emplace_back(gm / gw, gi); }
    MSSpectrum& s = v[surv[(Size)c]];
    s.clear(false); // peaks only: RT, precursor and userParams are the survivor's
    s.reserve(fused.size());
    for (const auto& f : fused)
    {
      s.emplace_back(f.first, (float)f.second);
    }
  };
  if (! omp_in_parallel())
  {
    std::exception_ptr fold_err; // [task-err] an exception escaping the region (bad_alloc) would abort without a message
#pragma omp parallel for schedule(dynamic, 64)
    for (long c = 0; c < (long)surv.size(); ++c)
    {
      try
      {
        fuse(c);
      }
      catch (...)
      {
        taskErr(fold_err);
      }
    }
    if (fold_err) { std::rethrow_exception(fold_err); }
  }
  else
  {
    teamFor((long)surv.size(), 64, [&](long c, int) { fuse(c); }); // [pipe] the pipelined tile writer; teamFor rethrows the first exception
  }
  Size w = 0;
  for (uint32_t i = 0; i < (uint32_t)v.size(); ++i)
  {
    if (grp[i] != NONE && grp[i] != i) { continue; }
    if (w != i) { v[w] = std::move(v[i]); }
    ++w;
  }
  const Size removed = v.size() - w;
  v.erase(v.begin() + (long)w, v.end());
  for (auto& sp : v)
  {
    if (sp.metaValueExists("spx_mono_ext"))
    {
      sp.removeMetaValue("spx_mono_ext"); // [twin-rt] never written
    }
  }
  return removed;
}

#ifdef _WIN32
static constexpr const char* NULL_DEVICE = "NUL"; // the encoders' own files: never written, only opened
#else
static constexpr const char* NULL_DEVICE = "/dev/null";
#endif
/**
  @brief [tile] Streams spectra to an mzML in BLOCKS: the header once, then every block encoded in
  parallel chunks (the tool-side twin of patches/openms-mzml-parallel-write.patch) and appended in
  order, the index offsets rebased as each chunk lands.

  Blocks arrive in the canonical order (a
  tile's sorted spectra; tiles are RT-disjoint), so the file is in the order the bulk writer
  produced. With the total known up front (setExpectedSize) the count is written exactly and the
  bytes from <spectrumList> on are the bulk writer's; without it a zero-padded placeholder is
  overwritten at the end -- same width, so the indexedmzML offsets stay valid.
*/
class TileWriter : public PlainMSDataWritingConsumer
{
public:
  /**
    @brief Open @p filename for writing.

    @param[in] filename Output path
  */
  explicit TileWriter(const String& filename): PlainMSDataWritingConsumer(filename), path_(filename)
  {
  }
  /// @brief Patch the spectrum count if it is still a placeholder; the base destructor writes what is left.
  ~TileWriter() override
  {
    if (count_off_ >= 0 && ofs_.is_open()) { patchCount_(); }
    // the base destructor closes <spectrumList> and writes the footer from spectra_offsets_ (nothing once finish() ran)
  }

  /**
    @brief Completes and closes the file: with no spectrum at all, the header and an empty <spectrumList> (a valid mzML,
    not a 0-byte file); then the count, the list end, the index and footer.

    Throws if any write failed (a full disk,
    a quota), so a cut-off .part is never renamed into place. The base destructor then has nothing left to write.

    @throws Exception::UnableToCreateFile if a write failed
  */
  void finish()
  {
    if (! started_writing_)
    {
      MapType empty;
      empty = settings_;
      Internal::MzMLHandler::writeHeader_(ofs_, empty, dps_, *validator_);
      started_writing_ = true;
    }
    if (! writing_spectra_)
    {
      ofs_ << "\t\t<spectrumList count=\"0\" defaultDataProcessingRef=\"dp_sp_0\">\n";
      writing_spectra_ = true;
    }
    if (count_off_ >= 0) { patchCount_(); }
    ofs_ << "\t\t</spectrumList>\n";
    writing_spectra_ = false;
    Internal::MzMLHandlerHelper::writeFooter_(ofs_, options_, spectra_offsets_, chromatograms_offsets_);
    started_writing_ = false;
    ofs_.close();
    if (ofs_.fail())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path_, "writing the mzML failed (disk full or quota?)");
    }
  }

  /**
    @brief Append a block of spectra, encoded in parallel chunks; the header is written before the first block.

    @param[in] block Spectra in canonical order
    @throws Exception::UnableToCreateFile if a write failed
  */
  void writeBlock(const vector<MSSpectrum>& block)
  {
    if (block.empty()) { return; }
    if (! started_writing_)
    {
      MapType dummy;
      dummy = settings_;
      dummy.addSpectrum(block[0]);
      Internal::MzMLHandler::writeHeader_(ofs_, dummy, dps_, *validator_);
      started_writing_ = true;
    }
    if (! writing_spectra_)
    {
      ofs_ << "\t\t<spectrumList count=\"";
      if (spectra_expected_ > 0) { ofs_ << spectra_expected_; }
      else
      {
        count_off_ = ofs_.tellp();
        ofs_ << "0000000000";
      }
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
    // teamFor, at most the team size; the chunks and their order are the same either way.
    const bool in_team = omp_in_parallel();
    const Size CHUNK = 256;
    const Size n_thr = (Size)std::max(1, in_team ? omp_get_num_threads() : omp_get_max_threads());
    while (enc_.size() < n_thr)
    {
      enc_.emplace_back(new TileWriter(NULL_DEVICE));
      enc_.back()->options_ = options_;
    }
    const Size wave = CHUNK * n_thr;
    for (Size w0 = 0; w0 < block.size(); w0 += wave)
    {
      const Size w1 = std::min(w0 + wave, block.size());
      const Size n_chunk = (w1 - w0 + CHUNK - 1) / CHUNK;
      vector<std::string> text(n_chunk);
      vector<vector<pair<std::string, Int64>>> offs(n_chunk);
      std::exception_ptr err;
      auto encode = [&](const long c, TileWriter& h) {
        try
        {
          std::ostringstream buf;
          buf.copyfmt(ofs_); // identical formatting state, or the numbers differ
          const Size lo = w0 + (Size)c * CHUNK, hi = std::min(lo + CHUNK, w1);
          const Size o0 = h.spectra_offsets_.size();
          for (Size i = lo; i < hi; ++i)
          {
            h.writeSpectrum_(buf, block[i], spectra_written_ + i, *h.validator_, /*renew_native_ids=*/true, dps_);
          }
          text[(Size)c] = buf.str();
          offs[(Size)c].assign(h.spectra_offsets_.begin() + o0, h.spectra_offsets_.end()); // buffer-relative
          h.spectra_offsets_.resize(o0);
        }
        catch (...)
        {
          taskErr(err);
        }
      };
      if (! in_team)
      {
#pragma omp parallel for schedule(dynamic, 1)
        for (long c = 0; c < (long)n_chunk; ++c)
        {
          encode(c, *enc_[(Size)omp_get_thread_num()]);
        }
      }
      else
      {
        teamFor((long)n_chunk, 1, [&](long c, int slot) { encode(c, *enc_[(Size)slot]); });
      }
      if (err) { std::rethrow_exception(err); }
      for (Size c = 0; c < n_chunk; ++c)
      {
        const Int64 base = (Int64)ofs_.tellp(); // where this chunk lands in the FILE
        ofs_ << text[c];
        for (const auto& p : offs[c])
        {
          spectra_offsets_.emplace_back(p.first, base + p.second);
        }
        std::string().swap(text[c]);
      }
    }
    spectra_written_ += block.size();
    if (! ofs_)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path_, "writing the mzML failed (disk full or quota?)");
    }
  }

private:
  /// @brief the placeholder count, overwritten in place: same width, so the index offsets stay valid
  void patchCount_()
  {
    const std::streampos end = ofs_.tellp();
    ofs_.seekp(count_off_);
    const char fill = ofs_.fill('0');
    ofs_ << std::setw(10) << spectra_written_;
    ofs_.fill(fill);
    ofs_.seekp(end);
  }
  String path_;
  std::streamoff count_off_ = -1;
  vector<std::unique_ptr<TileWriter>> enc_; ///< per-thread encoders (never write their null device)
};
} // namespace

/// @brief [s3] An S3 object fetched to a temporary local file under the same name, removed again when the run ends.
struct S3Copy
{
  std::filesystem::path dir;
  /// @brief Remove the temporary directory, if any.
  ~S3Copy()
  {
    std::error_code ec;
    if (! dir.empty()) { std::filesystem::remove_all(dir, ec); }
  }

  /**
    @brief Copy an S3 object to a local file named like it, in a fresh temporary directory.

    @param[in] uri The s3:// URI (query parameters configure the connection)
    @param[in] log Receives the log line of the copy
    @return Path of the local copy
    @throws OpenMS::Exception::InvalidParameter if this build cannot read S3 or the object is not an mzML, mzPeak or .raw file
    @throws OpenMS::Exception::UnableToCreateFile if the copy is incomplete
  */
  String fetch(const std::string& uri, const std::function<void(const String&)>& log)
  {
    if (! S3InputSource::isAvailable())
    {
      throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "this build of DIAspeXtract cannot read S3 (its OpenMS has no Arrow S3 support): " + uri);
    }
    const std::string key = uri.substr(0, uri.find('?')); // query parameters configure the connection, not the name
    const std::string name = key.substr(key.find_last_of('/') + 1);
    const FileTypes::Type t = name.empty() ? FileTypes::UNKNOWN : FileHandler::getTypeByFileName(name);
    if (t != FileTypes::MZML && t != FileTypes::MZPEAK && t != FileTypes::RAW)
    {
      throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "an S3 input must be one file (mzML, mzPeak or .raw); a Bruker .d folder has to be local: " + uri);
    }
    const Int64 size = S3InputSource::objectSize(uri); // a missing object: OpenMS logs Arrow's reason
    dir = std::filesystem::temp_directory_path() / ("diaspextract-s3-" + File::getUniqueName(false));
    std::filesystem::create_directories(dir);
    const std::filesystem::path local = dir / name;
    std::ofstream os(local, std::ios::binary);
    const auto t0 = std::chrono::steady_clock::now();
    const Int64 chunk = Int64(64) << 20;
    for (Int64 off = 0; off < size; off += chunk)
    {
      const std::string part = S3InputSource::readRange(uri, off, std::min(chunk, size - off));
      os.write(part.data(), (std::streamsize)part.size());
    }
    os.close();
    if (! os || (Int64)std::filesystem::file_size(local) != size)
    {
      throw OpenMS::Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, local.string(),
                                                  "the copy of " + uri + " is incomplete (disk full?)");
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    log("Fetched " + HTTPInputSource::redactURL(uri) + " (" + String(size / 1e6, false) + " MB in " + String(secs, false) + " s) to "
        + local.string());
    return local.string();
  }
};

// -noims:* options named on the command line, each followed by a space: an explicit value wins over the TOF defaults
static std::string g_noims_cli;

/// @brief The DIAspeXtract TOPP tool: builds pseudo-MS/MS spectra from DIA data for DDA search engines.
class TOPPDIAspeXtract : public TOPPBase
{
public:
  /// @brief Construct the tool and set its own version (with the OpenMS version in the verbose form).
  TOPPDIAspeXtract():
      // toolhandler_test=false: a standalone tool, absent from OpenMS's tool registry (TOPPBase rejects an unregistered name).
      TOPPBase("DIAspeXtract", "Builds pseudo-MS/MS spectra from DIA data (diaPASEF, Orbitrap, SWATH) for DDA search engines.", {}, false)
  {
    // Otherwise TOPPBase reports OpenMS's version as ours: version_ feeds the ini and the mzML software entry,
    // verboseVersion_ the --help Version line, which carries both numbers because a bug report needs the OpenMS build too.
    version_ = DIASPEXTRACT_VERSION;
    verboseVersion_ = String(DIASPEXTRACT_VERSION) + " (OpenMS " + VersionInfo::getVersion() + ")";
  }

  /**
    @brief [docurl] patches/openms-topp-external.patch makes this virtual.

    TOPPBase's default is an OpenMS doxygen page that
    does not exist for a tool outside the release; the first line of --help and every CTD's docurl point here instead.

    @return The tool's GitHub URL
  */
  std::string getDocumentationURL() const override
  { return "https://github.com/okohlbacher/diaspextract"; }

protected:
  /**
    @brief OpenMS's log streams are not thread-safe: two threads logging at once lose or repeat a line, or leave stray bytes.

    The window workers and the pipelined tile reader and writer log concurrently, so every line goes through one lock.

    @return The lock
  */
  static std::mutex& logMu_()
  {
    static std::mutex m;
    return m;
  }
  /**
    @brief TOPPBase::writeLogInfo_ under logMu_().

    @param[in] t The message
  */
  void writeLogInfo_(const String& t) const
  {
    std::lock_guard<std::mutex> g(logMu_());
    TOPPBase::writeLogInfo_(t);
  }
  /**
    @brief TOPPBase::writeLogWarn_ under logMu_().

    @param[in] t The message
  */
  void writeLogWarn_(const String& t) const
  {
    std::lock_guard<std::mutex> g(logMu_());
    TOPPBase::writeLogWarn_(t);
  }
  /**
    @brief TOPPBase::writeLogError_ under logMu_().

    @param[in] t The message
  */
  void writeLogError_(const String& t) const
  {
    std::lock_guard<std::mutex> g(logMu_());
    TOPPBase::writeLogError_(t);
  }

  /// @brief Register the tool's input/output, options and flags with their defaults and ranges.
  void registerOptionsAndFlags_() override
  {
    // a true/false option: a string option restricted to "true" and "false"
    auto reg_bool = [this](const String& name, const String& def, const String& desc, bool advanced = false) {
      registerStringOption_(name, "<true/false>", def, desc, false, advanced);
      setValidStrings_(name, {"true", "false"});
    };
    // skipexists, and no setValidFormats_: an s3:// URI is no local file, and with query parameters it has no extension
    // TOPPBase recognises (it then sniffs the content as a local file). main_ checks the input itself.
    registerInputFile_("in", "<file>", "",
                       "Input DIA run: a Bruker diaPASEF .d directory, an mzML (with ion mobility, or centroided Orbitrap/SWATH "
                       "DIA without it, which takes the separate noims path), a Thermo .raw file (noims path), or mzPeak. mzML and mzPeak "
                       "can also be read from S3: s3://bucket/key, with region or endpoint_override as URI query parameters "
                       "and the credentials from the usual AWS environment. (valid formats: 'mzML', 'mzpeak', 'd', 'raw')",
                       true, false, {"skipexists"});
    registerOutputFile_(
      "out", "<file>", "",
      "Output spectra. The format follows the extension: .mzML or .mzpeak. "
      "Search engines read mzML, so use .mzML if the next step is a search. mzPeak output runs as one tile and needs the whole run's memory.");
    setValidFormats_("out", {"mzpeak", "mzML"});
    registerStringOption_("out_type", "<type>", "",
                          "Set the output format instead of taking it from the extension of -out. Empty (default) = follow "
                          "the extension; an unknown extension gives mzPeak.",
                          false);
    setValidStrings_("out_type", {"", "mzpeak", "mzML"});

    registerTOPPSubsection_("gate", "Which fragments join a precursor's spectrum (diaPASEF input only)");
    registerDoubleOption_("gate:delta_im", "<1/K0>", 0.01,
                          "Largest ion-mobility (1/K0) difference between a precursor and a fragment. It also sets "
                          "the 1/K0 tolerance of peak picking and tracing; to widen only the fragment match, use gate:fragment_delta_im.",
                          false);
    registerDoubleOption_("gate:delta_rt", "<sec>", 4.0, "Largest apex retention-time difference (seconds) between a precursor and a fragment.",
                          false);
    registerDoubleOption_("gate:min_correlation", "<-1..1>", 0.3,
                          "Least correlation between the elution profiles of a precursor and a fragment. Raise "
                          "it for cleaner but sparser spectra.",
                          false);
    registerIntOption_("gate:min_correlation_points", "<n>", 3,
                       "Fewest shared points a fragment's correlation needs. Where gate:shrinkage applies, "
                       "gate:shrink_min_points is used instead.",
                       false, true);
    setMinInt_("gate:min_correlation_points", 0);

    registerTOPPSubsection_("assembly", "Building the spectra (diaPASEF input only)");
    registerIntOption_("assembly:min_fragments", "<n>", 3, "Write a spectrum only if it has at least this many fragments.", false);
    registerIntOption_("assembly:max_fragments", "<n>", 500,
                       "Keep at most this many fragments per spectrum, best-ranked first. The twin fold can "
                       "exceed it. At least assembly:min_fragments.",
                       false);
    registerDoubleOption_("assembly:im_weight_sigma", "<1/K0>", 0.005,
                          "Weight each fragment's intensity by a Gaussian of its 1/K0 distance to the "
                          "precursor (sigma in 1/K0), before the fragment cap. 0 = off.",
                          false);
    setMinFloat_("assembly:im_weight_sigma", 0.0);
    registerDoubleOption_("assembly:corr_power", "<k>", 2.0,
                          "Weight each fragment's intensity by corr^k, corr its elution correlation with the precursor. 0 = off.", false);
    setMinFloat_("assembly:corr_power", 0.0);
    registerDoubleOption_(
      "assembly:twin_im_tolerance", "<1/K0>", 0.005,
      "Fold spectra that claim the same precursor twice (same charge, m/z within "
      "20 ppm, 1/K0 within this tolerance) into one spectrum with the fragments of both. Negative = off; 0 folds only twins of identical 1/K0.",
      false);
    registerIntOption_("assembly:twin_rt_frames", "<n>", 0,
                       "Also fold twins whose apexes are up to n + 0.5 MS1 cycles apart and in the same "
                       "retention-time cell, when their precursor m/z lie within 0.05 ppm (one flight-time bin) and charge:mono_profile_check v2 "
                       "moved one of the two monoisotopes. 0 (default) = auto, from the run's MS1 peak width; -1 = the same MS1 frame only.",
                       false);

    registerTOPPSubsection_("trace", "Mass-trace detection (diaPASEF input only)");
    registerDoubleOption_("trace:mass_error_ppm", "<ppm>", 15.0, "m/z tolerance (ppm) of mass-trace detection.", false);
    registerDoubleOption_("trace:noise_threshold_int", "<int>", 100.0, "MS1 (precursor) noise level: peaks at or below it are not traced.", false);
    registerDoubleOption_("trace:ms2_noise_threshold_int", "<int>", 10.0, "MS2 (fragment) noise level: peaks at or below it are not traced.", false);
    registerDoubleOption_("trace:min_length_sec", "<sec>", -1.0,
                          "Shortest MS1 trace kept (seconds, first to last point). Negative (default) = auto: "
                          "trace:min_length_fwhm times the run's median MS1 peak width, at least half an MS1 cycle; 3 s when the run has too few "
                          "bright traces to measure the width. 0 = no filter.",
                          false);
    registerDoubleOption_("trace:min_length_fwhm", "<x>", 0.4,
                          "Auto trace:min_length_sec: the shortest MS1 trace as this fraction of the run's median MS1 peak width.", false, true);
    registerDoubleOption_("trace:ms1_chrom_peak_snr", "<x>", 3.0, "An MS1 trace's apex must exceed this times trace:noise_threshold_int.", false,
                          true);
    registerDoubleOption_("trace:ms2_chrom_peak_snr", "<x>", 1.0,
                          "A fragment trace's apex must exceed this times trace:ms2_noise_threshold_int. "
                          "Lower it to keep weaker fragments.",
                          false);
    registerDoubleOption_("trace:ms2_min_length_sec", "<sec>", 0.0, "Shortest fragment trace kept (seconds). 0 = no filter.", false);
    setMinFloat_("trace:ms2_min_length_sec", 0.0);
    registerTOPPSubsection_("tile", "Retention-time cells and tiles (diaPASEF input only)"); // registered, or --helphelp throws ElementNotFound
    registerDoubleOption_(
      "tile:rt_sec", "<sec>", 600.0,
      "Width (seconds) of the fixed retention-time cells the integer detector traces in; a trace "
      "that crosses a cell edge is cut there; 60 to a few thousand seconds. 0 or negative = the whole run as one cell. Changes the output.",
      false, true);
    registerIntOption_("tile:cells_per_tile", "<n>", 1,
                       "Process the run in tiles of this many retention-time cells to bound memory; fewer cells per "
                       "tile use less memory. 0 = the whole run as one tile, at several times the memory. mzPeak output always runs as one tile.",
                       false, true);
    setMinInt_("tile:cells_per_tile", 0);
    registerStringOption_("trace:detector", "<mode>", "integer",
                          "'integer' (default) traces on the instrument's flight-time bins and needs the "
                          "Bruker calibration; without it the run falls back to 'openms' with a warning. 'openms' uses OpenMS MassTraceDetection. "
                          "The two find partly different peptides; integer uses much less memory.",
                          false);
    setValidStrings_("trace:detector", {"openms", "integer"});
    registerDoubleOption_("trace:ms2_split_valleys", "<chrom_fwhm>", 7.0,
                          "Split fragment traces at chromatographic valleys; the value is the "
                          "expected peak width (seconds). 0 = off.",
                          false);
    setMinFloat_("trace:ms2_split_valleys", 0.0);
    registerFlag_("diag:selftest_twins", "Run the built-in twin-fold self-test and exit.", true);
    registerFlag_("diag:selftest_arena", "Run the built-in MS1 storage self-test and exit.", true);
    registerFlag_("diag:spectrum_meta",
                  "Write per-spectrum userParams for diagnosis: spx_p_frames (window frames with precursor signal), "
                  "spx_faint (1: at most gate:faint_max_frames of them), spx_p_apex (precursor apex intensity).",
                  true);
    registerDoubleOption_("trace:ms1_split_valleys", "<chrom_fwhm>", 7.0,
                          "Split precursor traces at chromatographic valleys, so two peptides that "
                          "elute apart at the same m/z stay separate; the value is the expected peak width (seconds). 0 = off.",
                          false);
    setMinFloat_("trace:ms1_split_valleys", 0.0);
    registerDoubleOption_("trace:max_span_sec", "<sec>", 120.0, "Trim each trace to at most this many seconds around its apex. 0 = no limit.", false);
    setMinFloat_("trace:max_span_sec", 0.0);

    // NB: "threads" is reserved by TOPPBase (the standard -threads option), so this lives under "perf".

    reg_bool("perf:malloc_trim", "true",
             "Return freed memory to the system between phases, and between tiles when free memory is below the tool's own "
             "footprint (Linux/glibc only). Lowers peak memory at some extra run time; set false when run time matters more.");
    reg_bool("perf:stream_load", "true",
             "Read the input frame by frame (Bruker .d, or mzPeak in a build that supports it; mzML is always read whole). false loads "
             "the whole run first: much more memory and a longer run. Changes the output.");
    registerIntOption_("trace:native_ms1_neighbors", "<n>", 0,
                       "Add n neighbouring MS1 frames on each side to every MS1 frame before peak picking "
                       "(Bruker .d with perf:stream_load true). 0 = off. Needs -dnoise:ms1 false.",
                       false, true);
    setMinInt_("trace:native_ms1_neighbors", 0);

    registerTOPPSubsection_("perf", "Speed and memory, diaPASEF input only (the output does not change unless an option says so)");
    registerStringOption_(
      "gate:short_traces", "<mode>", "concentric",
      "Fragment traces too short for a correlation: 'concentric' (default) keeps one "
      "whose retention-time span lies inside the precursor's (one frame of slack) and is centred on it; 'off' drops them. Applies to the "
      "precursors gate:adaptive selects.",
      false, true);
    setValidStrings_("gate:short_traces", {"off", "concentric"});
    registerDoubleOption_("gate:short_center_tol", "<frac>", 0.5,
                          "gate:short_traces concentric: how far the fragment's span centre may lie from the "
                          "precursor's, as a fraction of the precursor's half-span (at least one cycle).",
                          false, true);
    setMinFloat_("gate:short_center_tol", 0.0);
    registerIntOption_("gate:short_min_points", "<n>", 1, "gate:short_traces concentric: fewest points a short fragment trace needs.", false, true);
    setMinInt_("gate:short_min_points", 1);
    registerDoubleOption_("gate:fragment_delta_im", "<1/K0>", -1.0,
                          "Largest precursor-fragment 1/K0 difference for the fragment match only, leaving "
                          "the tracing tolerance unchanged. Negative (default) = gate:delta_im; otherwise 0 to about 0.1. Experimental.",
                          false, true);
    registerDoubleOption_("trace:im_slice_merge", "<1/K0>", 0.015,
                          "Integer detector: before tracing, merge a frame's fragment peaks in one "
                          "flight-time bin whose 1/K0 values lie within this distance, so one ion does not become several traces. Negative = off.",
                          false, true);
    setMaxFloat_("trace:im_slice_merge", 1.0); // 1/K0 units: anything wider than the mobility axis is a typo
    reg_bool("assembly:merge_same_bin", "true",
             "Merge a spectrum's fragments at the same m/z (one flight-time bin) into one peak before the fragment "
             "cap. Two peaks at one m/z break the fragment matching of OpenMS-based search engines.",
             true);
    registerDoubleOption_("assembly:foreign_precursor_weight", "<w>", 0.0,
                          "Multiply fragments that lie on the isotope envelope of another precursor "
                          "isolated in the same window by this weight. 0 (default) removes them; 1 = off. The precursor's own "
                          "envelope is never touched.",
                          false, false);
    setMinFloat_("assembly:foreign_precursor_weight", 0.0);
    setMaxFloat_("assembly:foreign_precursor_weight", 1.0);
    registerDoubleOption_("assembly:foreign_precursor_ppm", "<ppm>", 10.0,
                          "assembly:foreign_precursor_weight: m/z tolerance (ppm) between a foreign "
                          "envelope position and a fragment.",
                          false, true);
    setMinFloat_("assembly:foreign_precursor_ppm", 0.0);
    registerDoubleOption_("assembly:foreign_precursor_dim", "<1/K0>", 0.015,
                          "assembly:foreign_precursor_weight: largest 1/K0 difference between this "
                          "precursor and a foreign one.",
                          false, true);
    setMinFloat_("assembly:foreign_precursor_dim", 0.0);
    registerIntOption_("assembly:foreign_precursor_min_iso", "<n>", 2,
                       "assembly:foreign_precursor_weight: how many of a foreign envelope's four "
                       "positions (M to M+3) must carry a fragment.",
                       false, true);
    setMinInt_("assembly:foreign_precursor_min_iso", 1);
    setMaxInt_("assembly:foreign_precursor_min_iso", 4);
    registerDoubleOption_("assembly:fragment_fuse_ppm", "<ppm>", 25.0,
                          "Fuse a spectrum's fragments within this many ppm into one peak (summed "
                          "intensity, intensity-weighted m/z) before the fragment cap: one fragment ion often shows up as several peaks a few ppm "
                          "apart. 0 = off.",
                          false, false);
    setMinFloat_("assembly:fragment_fuse_ppm", 0.0);
    setMaxFloat_("assembly:fragment_fuse_ppm", 50.0);
    registerStringOption_("gate:adaptive", "<mode>", "faint",
                          "Which precursors get the extra gate rules (gate:shrinkage, gate:colocation, "
                          "gate:short_traces). 'faint' (default): only faint precursors, which also use gate:faint_delta_rt and "
                          "gate:faint_corr_power. 'auto': 'faint' on runs with many faint precursors (gate:auto_faint_fraction), otherwise no "
                          "extra rules; 'off': the rules apply to every precursor.",
                          false, true);
    setValidStrings_("gate:adaptive", {"off", "faint", "auto"});
    registerDoubleOption_("gate:auto_faint_fraction", "<frac>", 0.45,
                          "gate:adaptive auto: use the faint gate when at least this fraction of the "
                          "run's precursors are faint.",
                          false, true);
    setMinFloat_("gate:auto_faint_fraction", 0.0);
    setMaxFloat_("gate:auto_faint_fraction", 1.0);
    registerIntOption_("gate:faint_max_frames", "<n>", 4, "gate:adaptive: a precursor with signal in at most this many frames counts as faint.",
                       false, true);
    setMinInt_("gate:faint_max_frames", 1);
    registerDoubleOption_("gate:faint_delta_rt", "<sec>", 6.0, "gate:adaptive faint: gate:delta_rt for faint precursors.", false, true);
    setMinFloat_("gate:faint_delta_rt", 0.0);
    registerDoubleOption_("gate:faint_corr_power", "<k>", 1.0, "gate:adaptive faint: assembly:corr_power for faint precursors.", false, true);
    setMinFloat_("gate:faint_corr_power", 0.0);
    registerStringOption_(
      "gate:colocation", "<mode>", "ranges",
      "'ranges' (default): keep a fragment only if its apex lies inside the precursor's "
      "retention-time span, the precursor's apex inside the fragment's (each with one frame of slack), and the spans overlap by at least "
      "gate:colocation_overlap. 'off' skips this test.",
      false, true);
    setValidStrings_("gate:colocation", {"off", "ranges"});
    registerDoubleOption_("gate:colocation_overlap", "<frac>", 0.3,
                          "gate:colocation ranges: least overlap of the two spans, as a fraction of the shorter one.", false, true);
    setMinFloat_("gate:colocation_overlap", 0.0);
    setMaxFloat_("gate:colocation_overlap", 1.0);
    registerStringOption_("gate:shrinkage", "<mode>", "fisher",
                          "'fisher' (default): pull the correlation of a fragment with few shared points toward "
                          "a prior (gate:shrink_prior_r, gate:shrink_prior_sd) before the gate:min_correlation test, instead of requiring "
                          "gate:min_correlation_points points. 'off' uses the plain point count. At the defaults a fragment with at most 4 shared "
                          "points passes whatever its correlation (docs/options.md, Faint precursors).",
                          false, true);
    setValidStrings_("gate:shrinkage", {"off", "fisher"});
    registerDoubleOption_("gate:shrink_prior_r", "<r>", 0.6, "gate:shrinkage fisher: prior mean correlation (-1 < r < 1).", false, true);
    registerDoubleOption_("gate:shrink_prior_sd", "<sd>", 0.3,
                          "gate:shrinkage fisher: prior standard deviation in Fisher-z units; smaller pulls harder.", false, true);
    setMinFloat_("gate:shrink_prior_sd", 0.01);
    registerIntOption_("gate:shrink_min_points", "<n>", 1,
                       "gate:shrinkage fisher: fewest shared points to compute a correlation (replaces gate:min_correlation_points).", false, true);
    setMinInt_("gate:shrink_min_points", 1);
    registerIntOption_("perf:ms1_trace_bands", "<n>", 48,
                       "Trace MS1 in this many m/z bands in parallel. Changes the output slightly. 1 = off; 0 is read as 1.", false);
    reg_bool("perf:ms1_prune", "true",
             "Drop MS1 peaks at or below trace:noise_threshold_int while loading, which saves memory. false keeps every peak.");
    registerIntOption_("perf:trace_bands", "<n>", 12,
                       "Trace each window's fragments in this many m/z bands in parallel, for more parallelism than "
                       "the window count gives. Changes the output. 1 = off; 0 = auto (the output then depends on the thread count; the integer "
                       "detector refuses it with tile:rt_sec > 0 or a streamed .d).",
                       false);
    registerDoubleOption_(
      "perf:mem_fraction", "<f>", 0.75,
      "Fraction of free memory the window loop may use, re-checked for every window: MemAvailable, "
      "or the headroom under a cgroup limit when lower (macOS: free plus inactive memory). One window always runs, however large.",
      false);
    registerIntOption_("perf:max_concurrent_windows", "<n>", 0,
                       "Most isolation windows processed at once; fewer use less memory but take longer. 0 = "
                       "no cap (free memory still limits them).",
                       false);
    setMinInt_("perf:max_concurrent_windows", 0);
    registerIntOption_("perf:decode_threads", "<n>", 0,
                       "Threads that decode the Bruker frames (MS1 on the first read; MS2, sliced into the windows, "
                       "when the run or a tile is read; diaPASEF). 0 = all threads on the first read and the window team's free threads "
                       "on the tile reads; n = a team of n threads for both; 1 = serial. The spectra do not depend on it.",
                       false);
    setMinInt_("perf:decode_threads", 0);

    registerTOPPSubsection_("diag", "Diagnostics");
    registerStringOption_("diag:dump_ms1_tsv", "<prefix>", "",
                          "Write the MS1 traces and the precursors inferred from them to <prefix>.traces.tsv and "
                          "<prefix>.precursors.tsv, to see where a precursor was lost. Empty = off.",
                          false, true);

    registerIntOption_("max_charge", "<n>", 5, "Highest precursor charge considered (diaPASEF input only).", false);

    // NOTE: the subsection MUST be registered or printUsage_() throws ElementNotFound — which turns
    // any user parameter typo into a FATAL crash instead of a usage message.
    registerTOPPSubsection_("charge", "Isotopes, charge and monoisotope (diaPASEF input only)");
    registerIntOption_("charge:min_charge", "<n>", 1,
                       "Lowest precursor charge written. Default 1: singly charged precursors are real 1+ ions but "
                       "carry a higher FDR, so control FDR per charge downstream. 2 writes only 2+ and higher, the same spectra as at 1.",
                       false);
    reg_bool("charge:im_charge_veto", "true",
             "Use ion mobility to correct charge calls: a 1+ call on the run's 2+ mobility trend is re-called 2+, one "
             "on the 3+ trend is dropped. Genuine 1+ ions keep their charge.");
    registerDoubleOption_("charge:im_veto_band", "<sigma>", 2.5,
                          "charge:im_charge_veto: half-width of the 2+ and 3+ mobility bands, in robust "
                          "standard deviations of the fit.",
                          false);
    setMinFloat_("charge:im_veto_band", 0.0);
    registerStringOption_(
      "charge:mono_profile_check", "<false/v2>", "v2",
      "Check each precursor's monoisotope against the summed MS1 intensity at "
      "its isotope positions and move it by one isotope when averagine disagrees. 'v2' (default) = on; 'false' = off; 'v3' (experimental) = "
      "the monoisotope position with the best averagine cosine x observed fraction (Dinosaur/Monocle shift test).",
      false, false);
    setValidStrings_("charge:mono_profile_check", {"false", "v2", "v3"});
    registerDoubleOption_("charge:mono_position_min_ratio", "<q>", 0.2,
                          "charge:mono_profile_check v2: reject the lightest isotope when its "
                          "averagine-scaled intensity ratio is below this; the monoisotope then moves one isotope heavier.",
                          false, true);
    registerDoubleOption_("charge:mono_position_im_box", "<1/K0>", 0.0,
                          "charge:mono_profile_check v2: 1/K0 half-width of the isotope intensity sums. "
                          "0 = auto (per run, 0.01 to 0.03).",
                          false, true);
    registerDoubleOption_("charge:mono_position_extend", "<q>", 0.6,
                          "charge:mono_profile_check v2: when nothing was rejected, move the monoisotope "
                          "one isotope lighter if that position's ratio is at least this. 0 = off.",
                          false, true);
    reg_bool("charge:mono_position_iterate", "true",
             "charge:mono_profile_check v2: after a rejection, test the new lightest isotope too, and repeat. false tests once.", true);
    registerDoubleOption_("charge:iso_im_tolerance", "<1/K0>", 0.0,
                          "1/K0 tolerance for matching isotope partners. 0 (default) = measured per run "
                          "(0.01 to 0.05); a value above 0 fixes it.",
                          false);

    registerTOPPSubsection_("dnoise", "MS1 denoising before peak picking (Bruker .d only)");
    reg_bool("dnoise:ms1", "true",
             "Denoise every MS1 frame before peak picking, with a port of dnoise v0.1.0 (Garrett, Diedrich & Yates III, bioRxiv "
             "2026.08.27.747603). Bruker .d input only; other inputs pass through. Lowers peak memory; false turns it off.");

    registerTOPPSubsection_("noims",
                            "Input without ion mobility (Orbitrap, SWATH): a separate extraction path; the other sections do not apply to it");
    registerDoubleOption_("noims:min_corr", "<r>", 0.3, "Non-IMS path: minimum fragment-precursor profile correlation.", false, true);
    setMinFloat_("noims:min_corr", 0.0);
    setMaxFloat_("noims:min_corr", 1.0);
    registerDoubleOption_("noims:ppm", "<ppm>", 10.0, "Non-IMS path: mass-trace and isotope m/z tolerance.", false, true);
    registerStringOption_(
      "noims:instrument", "<type>", "auto",
      "Non-IMS path: mass analyzer the defaults are tuned for; auto = from the input's instrument (TOF without Orbitrap or FTICR = tof).", false,
      true);
    setValidStrings_("noims:instrument", {"auto", "orbitrap", "tof"});
    registerStringOption_("noims:tof", "<mode>", "off", "Non-IMS path: Biosaur2's TOF noise filter on MS1.", false, true);
    setValidStrings_("noims:tof", {"off", "on"});
    registerDoubleOption_("noims:min_intensity", "<i>", 1.0, "Non-IMS path: least MS1 peak intensity for precursor detection (Biosaur2 mini).", false,
                          true);
    registerIntOption_("noims:min_hill_scans", "<n>", 2, "Non-IMS path: fewest MS1 scans per precursor hill (Biosaur2 minlh).", false, true);
    registerIntOption_("noims:min_charge", "<z>", 1, "Non-IMS path: least precursor charge (Biosaur2 cmin).", false, true);
    registerIntOption_("noims:min_fragments", "<n>", 3, "Non-IMS path: fewest fragments a spectrum needs.", false, true);
    registerDoubleOption_("noims:min_fragment_intensity", "<i>", 0.0, "Non-IMS path: least apex intensity of a fragment trace.", false, true);
    registerIntOption_("noims:min_fragment_points", "<n>", 1, "Non-IMS path: fewest MS2 scans of a fragment trace.", false, true);
    registerIntOption_("noims:max_fragments", "<n>", 200, "Non-IMS path: most fragments kept per spectrum.", false, true);
    registerDoubleOption_("noims:weight_power", "<p>", 24.0, "Non-IMS path: fragment intensity is weighted by r^p.", false, true);
    registerDoubleOption_("noims:ms2_ppm", "<ppm>", 0.0, "Non-IMS path: MS2 mass-trace tolerance; 0 = noims:ppm.", false, true);
    setMinFloat_("noims:ms2_ppm", 0.0);
    registerDoubleOption_("noims:split", "<f>", 0.0,
                          "Non-IMS path: split MS2 traces at elution valleys below f x the lower flanking maximum; 0 = off.", false, true);
    setMinFloat_("noims:split", 0.0);
    setMaxFloat_("noims:split", 1.0);
    registerStringOption_(
      "noims:envelope_charges", "<mode>", "first",
      "Non-IMS path: MS2 isotope envelopes: first = a -1 isotope at the highest matching charge drops the envelope; all = try the lower charges.",
      false, true);
    setValidStrings_("noims:envelope_charges", {"first", "all"});
    registerStringOption_("noims:precursor_peaks", "<mode>", "keep", "Non-IMS path: the precursor's own isotope peaks in its spectrum: keep or drop.",
                          false, true);
    setValidStrings_("noims:precursor_peaks", {"keep", "drop"});
    registerDoubleOption_("noims:complement_boost", "<f>", 0.0,
                          "Non-IMS path: weight x f for fragments with a complementary partner (sum = precursor mass + 2 protons); 0 = off.", false,
                          true);
    setMinFloat_("noims:complement_boost", 0.0);
    registerIntOption_("noims:span_pad", "<n>", 0,
                       "Non-IMS path: widen each precursor's RT span by n window scans per side for the fragment correlation; 0 = off.", false, true);
    setMinInt_("noims:span_pad", 0);
    registerStringOption_("noims:envelope_pattern", "<mode>", "off",
                          "Non-IMS path: MS2 isotope envelopes must fit the averagine +1/mono intensity ratio; the best-fitting charge wins.", false,
                          true);
    setValidStrings_("noims:envelope_pattern", {"off", "on"});
    registerIntOption_("noims:hill_seeds", "<n>", 0,
                       "Non-IMS path: MS1 hills outside every feature become charge 2 and 3 precursors, kept with at least this many complementary "
                       "fragment pairs (0 = off).",
                       false, true);
    setMinInt_("noims:hill_seeds", 0);
    registerIntOption_("noims:pair_top", "<n>", 50, "Non-IMS path: hill seeds count complementary fragment pairs among the n most intense peaks.",
                       false, true);
    setMinInt_("noims:pair_top", 2);
    registerStringOption_("noims:seed_mono", "<mode>", "off", "Non-IMS path: drop a hill seed at a charge when a hill one isotope step lighter overlaps it in RT (it is then no monoisotopic peak).",
                          false, true);
    setValidStrings_("noims:seed_mono", {"off", "on"});
    registerStringOption_("noims:ms2_shift", "<ppm|auto>", "0",
                          "Non-IMS path: shift of the fragment m/z in ppm; auto = estimated from the y1 ions of tryptic peptides.", false, true);
    registerStringOption_(
      "noims:tag_source", "<mode>", "off",
      "Non-IMS path: write each spectrum's precursor source as userParam precursor_source (0 MS1 feature, 1 MS2 envelope, 2 corrected envelope, "
      "3 hill seed); hill seeds also get hill_pairs and hill_span (s).",
      false, true);
    setValidStrings_("noims:tag_source", {"off", "on"});
    setMinInt_("noims:max_fragments", 1);
    setMinFloat_("noims:weight_power", 0.0);
    setMinFloat_("noims:min_fragment_intensity", 0.0);
    setMinInt_("noims:min_fragment_points", 1);
    setMinFloat_("noims:min_intensity", 0.0);
    setMinInt_("noims:min_hill_scans", 1);
    setMinInt_("noims:min_charge", 1);
    setMaxInt_("noims:min_charge", 4);
    setMinInt_("noims:min_fragments", 1);

    // Bounds, not documentation: without these a negative count wraps through `Size` (min_fragments
    // -1 becomes ~1.8e19 and nothing is ever emitted) and a negative tolerance silently gates
    // everything out. Both fail as an empty-but-valid result, the worst failure mode here.
    setMinInt_("assembly:min_fragments", 1);
    setMinInt_("assembly:max_fragments", 1);
    setMinFloat_("assembly:twin_im_tolerance", -1.0); // any negative value = off
    setMinInt_("assembly:twin_rt_frames", -1);        // -1 = off (same frame), 0 = auto
    setMinInt_("max_charge", 1);
    setMinInt_("charge:min_charge", 1);
    setMinFloat_("trace:mass_error_ppm", 0.0);
    setMinFloat_("noims:ppm", 0.0);
    setMinInt_("perf:trace_bands", 0); // 0 = auto (documented)
    setMinFloat_("gate:delta_im", 0.0);
    setMinFloat_("gate:delta_rt", 0.0);
    setMinFloat_("gate:min_correlation", -1.0);
    setMaxFloat_("gate:min_correlation", 1.0);
    setMinFloat_("charge:iso_im_tolerance", 0.0);
    setMinFloat_("charge:mono_position_min_ratio", 0.0);
    setMinFloat_("charge:mono_position_im_box", 0.0);
    setMinFloat_("trace:min_length_fwhm", 0.0); // [minlen] 0 = the 2-point floor
    setMinFloat_("charge:mono_position_extend", 0.0);
    setMinFloat_("trace:ms1_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_chrom_peak_snr", 0.0);
    setMinFloat_("trace:ms2_noise_threshold_int", 0.0);
    setMinFloat_("trace:noise_threshold_int", 0.0);
    setMinInt_("perf:ms1_trace_bands", 0); // 0 is read as 1 (off): MS1 has no auto
    setMinFloat_("perf:mem_fraction", 0.05);
    setMaxFloat_("perf:mem_fraction", 0.95);
  }

  /**
    @brief IM-aware MassTraceDetection (banded) + valley splitting on @p map.

    The apex gate is snr * noise
    while membership needs only > noise, so MS1 and MS2 pass independent snr and min length.
    @p map is CLEARED as it is consumed, so it is not resident through the valley-splitting peak.

    @param[in,out] map The spectra to trace; cleared
    @param[in,out] store Receives the trace spans
    @param[in] delta_im Ion-mobility tolerance
    @param[in] noise Noise threshold
    @param[in] snr Apex gate factor
    @param[in] min_len Minimum trace length (s)
    @param[in] split_valleys Valley-splitting parameter
    @param[in] bands Number of m/z bands
    @param[in] splitter_scan_time The run's cycle time for the splitter, 0 = per trace
    @param[in,out] ms1_prune The MS1 call only: the pre-prune edge sample
    @return The traces
    @par Exceptions
      Rethrows the first exception of a band task, of whatever type it was.
  */
  vector<Trace> detectTraces_(PeakMap& map,
                              TraceStore& store,
                              double delta_im,
                              double noise,
                              double snr,
                              double min_len,
                              double split_valleys,
                              int bands,
                              double splitter_scan_time,     // the run's cycle time for the splitter, 0 = per trace
                              Ms1Prune* ms1_prune = nullptr) // [ms1-prune] the MS1 call only: the pre-prune edge sample
  {
    MassTraceDetection mtd;
    mtd.setLogType(ProgressLogger::NONE); // quiet + thread-safe (called from the parallel window loop)
    Param p = mtd.getParameters();
    p.setValue("mass_error_ppm", getDoubleOption_("trace:mass_error_ppm"));
    p.setValue("noise_threshold_int", noise);
    p.setValue("chrom_peak_snr", snr);
    p.setValue("ion_mobility_tolerance", delta_im); // IM-aware tracing
    p.setValue("min_trace_length", min_len);
    mtd.setParameters(p);

    vector<MassTrace> mts;
    // [outer] The traces as consecutive segments (the bands, then the valley-split chunks): every later step walks their
    // concatenation in order, so no gathered copy is needed.
    vector<vector<MassTrace>> seg;
    if (bands > 1)
    {
      // [bands] Partition the peaks by m/z and trace the bands concurrently. A trace cannot span more
      // than mass_error_ppm, so a halo of BAND_HALO_PPM_MUL tolerances makes the partition exact; a trace
      // is kept only by the band whose CORE holds its centroid.
      const double ppm = getDoubleOption_("trace:mass_error_ppm");
      // Band edges from m/z QUANTILES, not uniform m/z: peak density varies by orders of
      // magnitude across a window, and uniform edges would leave one band holding most of the work.
      vector<double> sample;
      if (ms1_prune && ms1_prune->on)
      {
        sample.swap(ms1_prune->sample); // [ms1-prune] taken before the prune; the pruned map's own would move the edges
      }
      else
      {
        for (const auto& s : map)
        {
          for (Size i = 0; i < s.size(); i += 97)
          {
            sample.push_back(s[i].getMZ());
          }
        }
      }
      sort(sample.begin(), sample.end());
      vector<double> edge(bands + 1);
      if (! sample.empty()) // else every edge stays 0
      {
        for (int b = 0; b <= bands; ++b)
        {
          edge[b] = sample[std::min(sample.size() - 1, (Size)((double)b / bands * sample.size()))];
        }
      }
      edge.front() = 0.0;
      edge.back() = std::numeric_limits<double>::max();
      if (ms1_prune && detOn()) // [det] condition (i) of the MS1 prune, gated on its own
      {
        vector<uint64_t> eb;
        for (double e : edge)
        {
          pushBits(eb, e);
        }
        writeLogInfo_("[ms1-edges] n=" + String(sample.size()) + " bands=" + String(bands) + " fnv=" + String(detDigest(std::move(eb))));
      }
      vector<double>().swap(sample);

      // Partition ONCE (peaks are copied into exactly one core band plus any halo it falls in),
      // then release the source map so the halo duplicates are the only overhead.
      // [dist-par] Distribute chunks of spectra in parallel; chunks concatenate in chunk order, so every
      // band gets the same spectra in the same order (byte-identical), and each source spectrum is
      // released as it is copied (transient ~1x the map).
      const int n_thr = omp_get_max_threads();
      const int nchunk = std::max(1, std::min((int)map.size(), 4 * n_thr));
      vector<vector<PeakMap>> csub((Size)nchunk, vector<PeakMap>(bands));
      std::exception_ptr dist_err; // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
      for (int ci = 0; ci < nchunk; ++ci)
      {
        try
        {
          const Size sp0 = (Size)ci * map.size() / nchunk, sp1 = (Size)(ci + 1) * map.size() / nchunk;
          for (Size spi = sp0; spi < sp1; ++spi)
          {
            auto& s = map[spi]; // non-const: each spectrum is released once distributed
            const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
            // isSorted() is a full traversal of the spectrum, and the spectrum does not change between
            // bands -- evaluating it inside the band loop walked every peak once per band.
            const bool src_sorted = s.isSorted();
            for (int b = 0; b < bands; ++b)
            {
              // Size each halo from ITS OWN edge: one h derived from edge[b] under-sizes the upper
              // side by the ppm growth across the band.
              const double lo = edge[b] - edge[b] * ppm * 1e-6 * BAND_HALO_PPM_MUL;
              const double hi = edge[b + 1] + edge[b + 1] * ppm * 1e-6 * BAND_HALO_PPM_MUL;
              MSSpectrum t;
              t.setRT(s.getRT());
              t.setMSLevel(1);
              OpenMS::DataArrays::FloatDataArray ia;
              ia.setName(Constants::UserParam::ION_MOBILITY);
              // sorted spectra: the band's range is two MZBegin lookups; the scan is the unsorted fallback.
              Size i0 = 0, i1 = s.size();
              if (src_sorted)
              {
                i0 = (Size)(s.MZBegin(lo) - s.begin());
                i1 = (Size)(s.MZBegin(hi) - s.begin());
                t.reserve(i1 - i0);
                ia.reserve(i1 - i0);
                for (Size i = i0; i < i1; ++i)
                {
                  t.push_back(s[i]);
                  if (ima && i < ima->size()) { ia.push_back((*ima)[i]); }
                }
              }
              else
              {
                for (Size i = 0; i < s.size(); ++i)
                {
                  if (s[i].getMZ() >= lo && s[i].getMZ() < hi)
                  {
                    t.push_back(s[i]);
                    if (ima && i < ima->size()) { ia.push_back((*ima)[i]); }
                  }
                }
              }
              if (! t.empty())
              {
                t.getFloatDataArrays().push_back(std::move(ia));
                csub[ci][b].addSpectrum(std::move(t));
              }
            }
            s.clear(true); // distributed: release it now, so the source and the copies do not both
          } // stand at full size
        }
        catch (...)
        {
          taskErr(dist_err);
        }
      }
      if (dist_err) { std::rethrow_exception(dist_err); }
      map.clear(true);
      vector<PeakMap> sub(bands);
      for (int b = 0; b < bands; ++b)
      {
        Size nsp = 0;
        for (int ci = 0; ci < nchunk; ++ci)
        {
          nsp += csub[ci][b].size();
        }
        sub[b].reserveSpaceSpectra(nsp);
        for (int ci = 0; ci < nchunk; ++ci)
        {
          for (auto& sp : csub[ci][b])
          {
            sub[b].addSpectrum(std::move(sp));
          }
          csub[ci][b].clear(true);
        }
      }
      if (ms1_prune && detOn()) // [det] which (frame, band) spectra exist: condition (ii) of the MS1 prune
      {
        String cnt;
        for (int b = 0; b < bands; ++b)
        {
          cnt += (b ? "," : "") + String(sub[b].size());
        }
        writeLogInfo_("[det] MS1 band spectra " + cnt);
      }
      vector<vector<PeakMap>>().swap(csub);
      vector<vector<MassTrace>> per(bands);
      std::exception_ptr band_err; // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
      for (int b = 0; b < bands; ++b)
      {
        try
        {
          taskFault("ms1_band");
          MassTraceDetection m2;
          m2.setLogType(ProgressLogger::NONE);
          m2.setParameters(p);
          vector<MassTrace> t;
          sub[b].sortSpectra();
          if (sub[b].size() < 3)
          {
            sub[b].clear(true);
            continue;
          }
          m2.run(sub[b], t);
          sub[b].clear(true);
          // keep only traces whose centroid lies in this band's CORE -> no duplicates from the halo
          for (auto& mt : t)
          {
            const double c = mt.getCentroidMZ();
            if (c >= edge[b] && c < edge[b + 1]) { per[b].push_back(std::move(mt)); }
          }
        }
        catch (...)
        {
          taskErr(band_err);
        }
      }
      if (band_err) { std::rethrow_exception(band_err); }
      seg = std::move(per); // [outer] the bands, in band order, are the segments
    }
    else
    {
      if (map.size() >= 3)
      {
        mtd.run(map, mts); // MTD aborts below 3 spectra; a tiny input is not an error
      }
      map.clear(true); // [compact] dead from here; ~1.1 GB/window freed before the EPD peak
      seg.resize(1);
      seg[0].swap(mts);
    }
    Size n_tr = 0;
    for (const auto& s : seg)
    {
      n_tr += s.size();
    }

    // Split mass traces at chromatographic valleys: MassTraceDetection never splits at a local
    // minimum, so co-eluting neighbours within the m/z tolerance merge. width_filtering is forced off
    // (a deletion filter). Per-trace and deterministic, so it runs as a balanced parallel loop (4 chunks
    // per band); both consumers re-sort (MS2 after detectTraces_, MS1 before inferPrecursors_).
    if (split_valleys > 0.0 && n_tr > 1)
    {
      const int nchunk = std::max(1, bands * 4); // 4 chunks per band -> dynamic load balance
      const Size n_in = n_tr;
      vector<Size> at(seg.size() + 1, 0); // [outer] segment s holds positions [at[s], at[s+1]) of the sequence
      for (Size s = 0; s < seg.size(); ++s)
      {
        at[s + 1] = at[s] + seg[s].size();
      }
      vector<vector<MassTrace>> parts(nchunk);
      std::exception_ptr split_err; // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
      for (int ci = 0; ci < nchunk; ++ci)
      {
        try
        {
          const Size lo = (Size)ci * n_in / nchunk, hi = (Size)(ci + 1) * n_in / nchunk;
          if (lo >= hi) { continue; }
          // [outer] positions [lo, hi) of the segment sequence, moved out segment by segment (distinct chunks move
          // distinct elements)
          vector<MassTrace> chunk;
          chunk.reserve(hi - lo);
          Size s = (Size)(std::upper_bound(at.begin(), at.end(), lo) - at.begin()) - 1;
          for (Size g = lo; g < hi; ++g)
          {
            while (g >= at[s + 1])
            {
              ++s;
            }
            chunk.push_back(std::move(seg[s][g - at[s]]));
          }
          ElutionPeakDetection epd;
          setSplitter(epd, split_valleys, splitter_scan_time);
          parts[ci].reserve(chunk.size());
          epd.detectPeaks(chunk, parts[ci]);
        }
        catch (...)
        {
          taskErr(split_err);
        }
      }
      if (split_err) { std::rethrow_exception(split_err); }
      seg.swap(parts); // [outer] the chunks' outputs, in chunk order, are the segments now
      n_tr = 0;
      for (const auto& s : seg)
      {
        n_tr += s.size();
      }
    }
    vector<Trace> out;
    // [outer] The MS1 call converts in parallel: pass 1 sizes every span (makeSpan's frame range), the offsets are the
    // running sum from the arena's end -- exactly where the serial appends put them -- and pass 2 fills each span at
    // its offset with toTrace's own statements. Every trace writes only its own span and its own slot, so the arena and
    // the trace list are byte-identical to the serial loop's. Contiguous ranges per task keep the walk sequential.
    bool conv_par = ms1_prune && n_tr >= 4096 && omp_in_parallel() && omp_get_num_threads() > 1;
    vector<Size> at, off; // segment s holds positions [at[s], at[s+1]) of the sequence; off: span length, then offset
    auto seg_of = [&at](Size g) { return (Size)(std::upper_bound(at.begin(), at.end(), g) - at.begin()) - 1; };
    const int ntask = conv_par ? std::max(1, std::min(4 * omp_get_num_threads(), (int)(n_tr / 1024))) : 1; // bounded: no task-count cliff
    std::exception_ptr conv_err;                                                                           // [task-err]
    if (conv_par)
    {
      at.assign(seg.size() + 1, 0);
      for (Size s = 0; s < seg.size(); ++s)
      {
        at[s + 1] = at[s] + seg[s].size();
      }
      off.assign(n_tr + 1, 0);
#pragma omp taskloop num_tasks(ntask) default(shared)
      for (int c = 0; c < ntask; ++c)
      {
        try
        {
          const Size lo = (Size)c * n_tr / ntask, hi = (Size)(c + 1) * n_tr / ntask;
          if (lo >= hi) { continue; }
          Size s = seg_of(lo);
          for (Size g = lo; g < hi; ++g)
          {
            while (g >= at[s + 1])
            {
              ++s;
            }
            off[g + 1] = spanLength(seg[s][g - at[s]], store);
          }
        }
        catch (...)
        {
          taskErr(conv_err);
        }
      }
      if (conv_err) { std::rethrow_exception(conv_err); }
      off[0] = store.inten.size();
      for (Size g = 0; g < n_tr; ++g)
      {
        off[g + 1] += off[g]; // off[g] = the arena size before trace g's append
      }
      // Trace::off is 32-bit: past 2^32 - 1 floats the serial loop below runs instead, and makeSpan refuses there
      if (off[n_tr] > (Size)std::numeric_limits<uint32_t>::max()) { conv_par = false; }
    }
    if (conv_par)
    {
      store.inten.resize(off[n_tr], 0.0f);
      out.resize(n_tr);
#pragma omp taskloop num_tasks(ntask) default(shared)
      for (int c = 0; c < ntask; ++c)
      {
        try
        {
          const Size lo = (Size)c * n_tr / ntask, hi = (Size)(c + 1) * n_tr / ntask;
          if (lo >= hi) { continue; }
          vector<pair<uint32_t, float>> pts;
          Size s = seg_of(lo);
          for (Size g = lo; g < hi; ++g)
          {
            while (g >= at[s + 1])
            {
              ++s;
            }
            MassTrace& mt = seg[s][g - at[s]];
            out[g] = toTraceAt(mt, store, off[g], pts);
            mt = MassTrace(); // the payload goes now, as in the serial loop
          }
        }
        catch (...)
        {
          taskErr(conv_err);
        }
      }
      if (conv_err) { std::rethrow_exception(conv_err); }
      vector<vector<MassTrace>>().swap(seg);
    }
    else
    {
      out.reserve(n_tr);
      // Each MassTrace is released as soon as it is converted.
      for (auto& sv : seg) // [outer] the segments in order
      {
        for (Size i = 0; i < sv.size(); ++i)
        {
          out.push_back(toTrace(sv[i], store));
          sv[i] = MassTrace();
        } // move-assign an empty one: the payload goes now, not at return
        vector<MassTrace>().swap(sv);
      }
    }
    return out;
  }

  /**
    @brief [isotol] charge:iso_im_tolerance 0 (auto): the run's isotope-partner 1/K0 tolerance, measured on its own MS1 traces
    before the greedy walk.

    `ms1` is the content-sorted list (m/z ascending), so the sample and every search are the same
    at any thread count and tiling. Every ceil(N / 200000)-th trace is a probe. At each charge z <= max_charge it is
    probed at its +1 isotope position (m/z + 1.0033548 / z) and at two controls 50 ppm either side, where no isotope of the
    probe at that charge can sit (another charge's can: near m/z 1,000 a 4+ ion's M+1 lies 0.2 ppm from the 5+ control,
    so a true partner there counts as a chance one), in the walk's partner window (probe m/z x mass_ppm, gate:delta_rt), and the nearest other
    trace's |d 1/K0| is kept. A tolerance t then admits P_iso(t) traces per isotope probe, of which about P_ctl(t) (the mean of the two controls)
    are chance traces: t* maximises P_iso(t) - 2 P_ctl(t), the true partners admitted minus the chance partners admitted, on a 0.001 grid over
    [0.002, 0.050] (ties to the wider), then clamped to [0.01, 0.05]. `info` gets the stamp text.

    @param[in] ms1 The MS1 traces, sorted by m/z
    @param[in] ms1st Their store
    @param[in] max_charge Highest charge probed
    @param[in] delta_rt RT tolerance
    @param[in] mass_ppm Mass tolerance (ppm)
    @param[out] info The stamp text
    @return The 1/K0 tolerance
    @throws Exception::Precondition if the MS1 traces are not sorted by m/z
  */
  double autoIsoImTol_(const vector<Trace>& ms1, const TraceStore& ms1st, int max_charge, double delta_rt, double mass_ppm, String& info)
  {
    const double CTL = 50e-6, LO = 0.01, HI = 0.05;
    const int NG = 49; // grid t_x = 0.002 + 0.001 x, x = 0..48
    auto tg = [](int x) { return 0.002 + 0.001 * x; };
    const Size n_ms1_scans = ms1.size();
    if (n_ms1_scans < 2) // max_charge >= 1 (setMinInt_)
    {
      info = "auto 0.05 (too few MS1 traces to measure)";
      writeLogInfo_("[isotol] auto: " + String(n_ms1_scans) + " MS1 traces, too few to measure; 1/K0 tolerance 0.05");
      return HI;
    }
    vector<double> mz(n_ms1_scans), rt(n_ms1_scans), im(n_ms1_scans);
    for (Size i = 0; i < n_ms1_scans; ++i)
    {
      mz[i] = ms1[i].mz;
      rt[i] = rtOf(ms1st, ms1[i]);
      im[i] = ms1[i].im;
    }
    if (! std::is_sorted(mz.begin(), mz.end()))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[isotol] the MS1 traces are not sorted by m/z");
    }
    const Size stride = std::max<Size>(1, (n_ms1_scans + 199999) / 200000), ns = (n_ms1_scans + stride - 1) / stride;
    const int z_max = max_charge;
    vector<int8_t> gx(ns * (Size)z_max * 3, -1); // per (probe, z, kind): first grid index x with gap <= t_x, -1 none
#pragma omp parallel for schedule(dynamic, 256)
    for (long q = 0; q < (long)ns; ++q)
    {
      const Size s = (Size)q * stride;
      const double tol = mz[s] * mass_ppm * 1e-6;
      for (int z = 1; z <= z_max; ++z)
      {
        for (int c = 0; c < 3; ++c)
        {
          const double base = mz[s] + ISO_STEP / z;
          const double target = c == 0 ? base : c == 1 ? base * (1.0 + CTL) : base * (1.0 - CTL);
          const Size klo = (Size)(std::lower_bound(mz.begin(), mz.end(), target - tol) - mz.begin());
          const Size khi = (Size)(std::upper_bound(mz.begin(), mz.end(), target + tol) - mz.begin());
          double g = 1.0;
          for (Size k = klo; k < khi; ++k)
          {
            if (k == s || fabs(rt[k] - rt[s]) > delta_rt) { continue; }
            const double d = fabs(im[k] - im[s]);
            if (d < g) { g = d; }
          }
          int8_t b = -1;
          for (int x = 0; x < NG; ++x)
          {
            if (g <= tg(x))
            {
              b = (int8_t)x;
              break;
            }
          }
          gx[((Size)q * (Size)z_max + (Size)(z - 1)) * 3 + (Size)c] = b;
        }
      }
    }
    vector<long long> cum(3 * NG, 0); // kind c, grid x: probes whose nearest gap is <= t_x
    for (Size i = 0; i < gx.size(); ++i)
    {
      if (gx[i] >= 0) { ++cum[(i % 3) * NG + (Size)gx[i]]; }
    }
    for (int c = 0; c < 3; ++c)
    {
      for (int x = 1; x < NG; ++x)
      {
        cum[c * NG + x] += cum[c * NG + x - 1];
      }
    }
    int best = 0;
    long long bj = 0;
    for (int x = 0; x < NG; ++x)
    {
      const long long j_score = cum[x] - cum[NG + x] - cum[2 * NG + x];
      if (x == 0 || j_score >= bj)
      {
        best = x;
        bj = j_score;
      }
    }
    const double tstar = tg(best), tol = std::min(HI, std::max(LO, tstar));
    const int xc = (int)std::lround((tol - 0.002) / 0.001);
    const double np = (double)ns * z_max;
    auto rate = [&](int c, int x) { return (double)cum[c * NG + x] / np; };
    char b[400];
    std::snprintf(b, sizeof b,
                  "auto %.3f (t* %.3f of isotope minus 2 x control traces, clamped to [0.01, 0.05]; per probe within %.3f: isotope %.3f, control "
                  "%.3f; %zu probes x %d charges)",
                  tol, tstar, tol, rate(0, xc), 0.5 * (rate(1, xc) + rate(2, xc)), ns, z_max);
    info = String(b);
    std::snprintf(b, sizeof b,
                  "[isotol] auto: 1/K0 tolerance %.3f for isotope partners (t* %.3f, clamped to [0.01, 0.05]). %zu probes (stride %zu over %zu MS1 "
                  "traces) x %d charges: "
                  "the nearest trace lies within the tolerance at the +1 isotope position for %.3f of the probes, at the +-50 ppm controls for %.3f "
                  "(within 0.05: %.3f, %.3f)",
                  tol, tstar, ns, stride, n_ms1_scans, z_max, rate(0, xc), 0.5 * (rate(1, xc) + rate(2, xc)), rate(0, NG - 1),
                  0.5 * (rate(1, NG - 1) + rate(2, NG - 1)));
    writeLogInfo_(String(b));
    return tol;
  }

  /**
    @brief [outer] The m/z order inferPrecursors_ and monoPositionPass_ share: the indices of @p ms1 by m/z, equal m/z by index.

    The MS1 traces arrive content-sorted (m/z is the first key), so unless two of them share an m/z the m/z values
    strictly increase in container order, the identity is then the only sorted order, and it is taken without sorting.
    m/z is quantised, so equal values are common; [xplat] the index breaks those ties -- std::sort's own order of equal
    keys differs between libstdc++, libc++ and the MS STL, and moved Linux's precursors against macOS and Windows.

    @param[in] ms1 The MS1 traces
    @param[out] by_mz Their indices in m/z order
  */
  static void mzOrder_(const vector<Trace>& ms1, vector<Size>& by_mz)
  {
    const Size n_ms1_scans = ms1.size();
    by_mz.resize(n_ms1_scans);
#pragma omp parallel for schedule(static)
    for (long i = 0; i < (long)n_ms1_scans; ++i)
    {
      by_mz[(Size)i] = (Size)i;
    }
    Size nonincr = 0;
#pragma omp parallel for schedule(static) reduction(+ : nonincr)
    for (long i = 1; i < (long)n_ms1_scans; ++i)
    {
      if (! (ms1[(Size)i - 1].mz < ms1[(Size)i].mz)) { ++nonincr; }
    }
    if (nonincr)
    {
      parallelSortTotal(by_mz, [&ms1](Size a, Size b) { return ms1[a].mz < ms1[b].mz || (ms1[a].mz == ms1[b].mz && a < b); });
    }
  }

  /**
    @brief Isotope/charge inference on MS1 traces: a greedy partner-count walk, brightest seed first, with exact ordered
    speculation.

    @param[in] ms1 The MS1 traces
    @param[in] ms1st Their store
    @param[in] max_charge Highest charge tried
    @param[in] delta_rt RT tolerance
    @param[in] iso_im_tol Isotope-partner 1/K0 tolerance
    @param[in] mass_ppm Mass tolerance (ppm)
    @param[out] v2_out Optional walk record for monoPositionPass_, parallel to the result
    @return The precursor candidates
    @par Exceptions
      Rethrows the first exception of a parallel part, of whatever type it was.
  */
  vector<PrecursorCand> inferPrecursors_(const vector<Trace>& ms1,
                                         const TraceStore& ms1st,
                                         int max_charge,
                                         double delta_rt,
                                         double iso_im_tol,
                                         double mass_ppm,
                                         vector<MonoV2In>* v2_out = nullptr)
  {
    const Size n_ms1_scans = ms1.size();

    // m/z-sorted index for O(log N) isotope-partner lookup (avoids O(N^2) over all MS1 traces).
    vector<Size> by_mz;
    mzOrder_(ms1, by_mz); // [outer]
    vector<double> sorted_mz(n_ms1_scans);
    // [soa] The candidate scan's RT and IM, laid out in m/z order beside sorted_mz (sequential reads).
    vector<double> sorted_rt(n_ms1_scans), sorted_im(n_ms1_scans);
#pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long)n_ms1_scans; ++ii)
    {
      const Size i = (Size)ii;
      const Trace& t = ms1[by_mz[i]];
      sorted_mz[i] = t.mz;
      sorted_rt[i] = rtOf(ms1st, t);
      sorted_im[i] = t.im;
    }
    // [mzidx] Bucket index over sorted_mz keyed on the top bits of the IEEE-754 pattern: monotone for
    // positive doubles, so the bucket width is relative (ppm-shaped) and the std search on the bucket's
    // sub-range returns the same position as on the whole array.
    const double mz_lo = n_ms1_scans ? sorted_mz.front() : 0.0, mz_hi = n_ms1_scans ? sorted_mz.back() : 1.0;
    unsigned kshift = 32;
    auto raw_key = [](double x, unsigned sh) {
      uint64_t b;
      std::memcpy(&b, &x, sizeof b);
      return b >> sh;
    };
    while (kshift < 52 && (raw_key(mz_hi, kshift) - raw_key(mz_lo, kshift)) > (1u << 23))
    {
      ++kshift;
    }
    const uint64_t klo = n_ms1_scans ? raw_key(mz_lo, kshift) : 0;
    const Size kspan = n_ms1_scans ? (Size)(raw_key(mz_hi, kshift) - klo) : 0;
    auto bucket_of = [&](double x) -> Size {
      if (! (x > mz_lo)) { return 0; }
      const uint64_t k = raw_key(x, kshift);
      return k >= klo + kspan ? kspan : (Size)(k - klo);
    };
    vector<uint32_t> bstart(kspan + 2, 0);
    {
// bstart[b] = the first i with bucket_of(sorted_mz[i]) >= b, else N. [outer] bucket_of is monotone over sorted_mz, so
// the buckets (bucket_of(sorted_mz[i-1]), bucket_of(sorted_mz[i])] start at i: every b is written by exactly one i.
#pragma omp parallel for schedule(static)
      for (long ii = 0; ii < (long)n_ms1_scans; ++ii)
      {
        const Size i = (Size)ii, b_hi = std::min(bucket_of(sorted_mz[i]), kspan);
        for (Size b = (i == 0 ? 0 : bucket_of(sorted_mz[i - 1]) + 1); b <= b_hi; ++b)
        {
          bstart[b] = (uint32_t)i;
        }
      }
      for (Size b = (n_ms1_scans ? std::min(bucket_of(sorted_mz[n_ms1_scans - 1]), kspan) + 1 : 0); b <= kspan; ++b)
      {
        bstart[b] = (uint32_t)n_ms1_scans;
      }
      bstart[kspan + 1] = (uint32_t)n_ms1_scans;
    }
    auto lbmz = [&](double x) -> Size {
      if (n_ms1_scans == 0 || x <= mz_lo) { return 0; }
      if (x > mz_hi) { return n_ms1_scans; }
      const Size b = bucket_of(x);
      return (Size)(std::lower_bound(sorted_mz.begin() + bstart[b], sorted_mz.begin() + bstart[b + 1], x) - sorted_mz.begin());
    };
    auto ubmz = [&](double x) -> Size {
      if (n_ms1_scans == 0 || x < mz_lo) { return 0; }
      if (x >= mz_hi) { return n_ms1_scans; }
      const Size b = bucket_of(x);
      return (Size)(std::upper_bound(sorted_mz.begin() + bstart[b], sorted_mz.begin() + bstart[b + 1], x) - sorted_mz.begin());
    };

    // Apex intensities in one contiguous array: the sort and every isotope walk read them.
    vector<float> inten(n_ms1_scans);
    vector<Size> order(n_ms1_scans);
#pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long)n_ms1_scans; ++ii)
    {
      const Size i = (Size)ii;
      inten[i] = (float)intensityOf(ms1st, ms1[i]);
      order[i] = i;
    }
    // [outer] a total order (the index breaks every tie; the coordinates are finite, non-finite traces were dropped before
    // the content sort), so the sorted sequence is unique and the parallel sort returns std::sort's
    parallelSortTotal(order, [&](Size a, Size b) {
      if (inten[a] != inten[b]) { return inten[a] > inten[b]; }
      if (ms1[a].mz != ms1[b].mz) { return ms1[a].mz < ms1[b].mz; }
      return a < b; // total order so greedy isotope assignment is deterministic
    });

    vector<bool> used(n_ms1_scans, false);
    vector<PrecursorCand> out;
    // [spec] Ordered speculation, exact: the only mutable state a seed reads is `used`, which is monotone.
    // A batch of seeds is evaluated concurrently against a frozen `used`, recording the partners each
    // search returned; the commit walks the batch in the original order, skips a seed that is now used,
    // and re-evaluates one whose recorded partner has since been claimed. A search that returned -1 saw
    // only gate failures or used candidates, both stable under a monotone bitmap, so it needs no
    // dependency.
    // [monocheck-v2] the count walk's result for monoPositionPass_ (MonoV2In), which runs after the IM charge veto.
    struct SeedRes
    {
      vector<PrecursorCand> emit;
      vector<Size> claims;
      vector<Size> deps;
      MonoV2In v2;
    };
    auto eval_seed = [&](Size seed, SeedRes& seed_res) {
      seed_res.emit.clear();
      seed_res.claims.clear();
      seed_res.deps.clear();
      seed_res.v2 = MonoV2In();
      if (used[seed]) { return; }
      const Trace& s = ms1[seed];
      const double s_rt = rtOf(ms1st, s); // [soa-rt] hoisted: the candidate scan reads it per probe
      const double tol = s.mz * mass_ppm * 1e-6;

      // Find an unused isotope partner near `target` m/z, co-localized in RT and IM (binary-search
      // the m/z index, scan only the ppm-window candidates): the first in m/z order; -1 if none.
      auto find_partner = [&](double target) -> long {
        Size klo = lbmz(target - tol), khi = ubmz(target + tol); // [mzidx] same positions, local search
        for (Size k = klo; k < khi; ++k)
        {
          if (fabs(sorted_rt[k] - s_rt) > delta_rt || fabs(sorted_im[k] - s.im) > iso_im_tol)
          {
            continue; // [soa]
          }
          Size j = by_mz[k];
          if (used[j] || j == seed) { continue; }
          seed_res.deps.push_back((Size)j);
          return (long)j;
        }
        return -1;
      };

      // Contiguous walk from the seed with break-on-first-miss, monoisotope = furthest-left partner reached, z by
      // partner count, ties to lowest z.
      int best_z = 0, best_n = 1, best_nl = 0;
      double best_mono = s.mz;
      vector<Size> best_partners;
      for (int z = 1; z <= max_charge; ++z)
      {
        vector<Size> partners;
        double mono = s.mz;
        vector<long> lefts;
        for (int k = 1; k <= 5; ++k) // lighter partners define the monoisotope
        {
          long j = find_partner(s.mz - k * ISO_STEP / z);
          if (j < 0)
          {
            break; // BREAK on first miss (contiguous only)
          }
          lefts.push_back(j);
          partners.push_back((Size)j);
          mono = ms1[j].mz;
        }
        for (int k = 1; k <= 5; ++k) // heavier partners add confidence
        {
          long j = find_partner(s.mz + k * ISO_STEP / z);
          if (j < 0) { break; }
          partners.push_back((Size)j);
        }
        const int n = 1 + (int)partners.size();
        if (n > best_n)
        {
          best_n = n;
          best_z = z;
          best_mono = mono;
          best_partners = partners;
          best_nl = (int)lefts.size();
        }
      }
      // [monocheck-v2] Record what monoPositionPass_ needs (v2 reads only the immutable arrays and changes only
      // mono_mz, so deciding it after the loop gives the same result as deciding it here, lets an automatic im_box
      // use the run's own isotope scatter, and lets the decision follow the IM charge veto).
      if (cfg_mono_check >= 2 && best_z > 0)
      {
        seed_res.v2.z = best_z;
        seed_res.v2.nl = best_nl;
        seed_res.v2.mono = best_mono;
        for (int i = 0; i < best_nl; ++i)
        {
          seed_res.v2.left_mz[(Size)i] = ms1[best_partners[(Size)i]].mz;
        }
      }
      PrecursorCand pc0;
      pc0.trace_idx = seed;
      pc0.rt = s_rt;
      pc0.im = s.im;
      pc0.mono_mz = best_mono;
      pc0.charge = best_z;     // 0 => unknown: a guessed precursor, dropped later
      pc0.n_isotopes = best_n; // envelope confidence for the mono call
      seed_res.emit.push_back(pc0);
      for (Size j : best_partners)
      {
        seed_res.claims.push_back(j); // de-isotope
      }
    };

    // [outer] batch size: re-evaluations in the serial commit vs fork/join count
    constexpr Size spec_batch = 2048;
    // [outer] about four speculation chunks per thread (at most 16 seeds), so a small batch feeds every thread
    const int spec_chunk = (int)std::max<Size>(1, std::min<Size>(16, spec_batch / (4 * (Size)std::max(1, omp_get_max_threads()))));
    vector<SeedRes> res(spec_batch);
    vector<MonoV2In> v2in; // [monocheck-v2] parallel to `out` when cfg_mono_check == 2
    for (Size b0 = 0; b0 < order.size(); b0 += spec_batch)
    {
      const Size b1 = std::min(b0 + spec_batch, order.size());
      // Phase 1: speculate. Nothing here writes `used`, so the concurrent reads are safe.
      std::exception_ptr spec_err;
#pragma omp parallel for schedule(dynamic, spec_chunk)
      for (long ii = (long)b0; ii < (long)b1; ++ii)
      {
        SeedRes& seed_res = res[(Size)ii - b0];
        try
        {
          eval_seed(order[(Size)ii], seed_res);
        } // clears R; a used seed leaves it empty
        catch (...)
        {
          taskErr(spec_err);
        } // [task-err] bad_alloc would abort the run without a message
      }
      if (spec_err) { std::rethrow_exception(spec_err); }
      // Phase 2: commit in the original seed order.
      for (Size ii = b0; ii < b1; ++ii)
      {
        const Size seed = order[ii];
        if (used[seed]) { continue; }
        SeedRes& seed_res = res[ii - b0];
        bool stale = false;
        for (Size d : seed_res.deps)
        {
          if (used[d])
          {
            stale = true;
            break;
          }
        }
        if (stale) { eval_seed(seed, seed_res); }
        for (const auto& pc : seed_res.emit)
        {
          out.push_back(pc);
          if (cfg_mono_check >= 2) { v2in.push_back(seed_res.v2); }
        }
        used[seed] = true;
        for (Size j : seed_res.claims)
        {
          used[j] = true;
        }
      }
    }
    if (v2_out)
    {
      v2_out->swap(v2in); // [monocheck-v2] parallel to `out`, for monoPositionPass_
    }
    return out;
  }

  /**
    @brief [monocheck-v2] charge:mono_profile_check v2, after the IM charge veto: each count walk is judged at its
    precursor's FINAL charge.

    A z=1 walk the veto re-calls 2+ stepped 1.003 m/z, two isotopes at 2+; judged at z=1
    spacing its steps and extension land two isotopes off, so its walk is read at the final charge:
    its leftmost partner is position -2 nl there, and the positions between are filled by whatever traces sit on
    them. PS(k) sums the apex intensity of EVERY MS1 trace (claimed or not, the seed included) at the seed's isotope
    position k within im_box 1/K0, in m/z order (deterministic): find_partner takes the first box candidate in m/z
    order and IM slicing splits one isotope over several traces, so single-trace ratios are noisy. A: an identified
    walk carries at most one spurious step and it is the leftmost (engines search -1..2), so step j = the leftmost is
    rejected when PS(-j) / PS(-j+1) x lambda < min_ratio and the mono moves one isotope right (iterate: repeat).
    B: when A rejected nothing and the walk ended on a miss (< 5 steps), the
    true mono may be claimed by a brighter precursor, so the next lighter position becomes the mono, unclaimed, on a
    ratio >= extend. Only mono_mz changes: claims, charge and n_isotopes are the greedy loop's and the veto's.
    im_box 0 (auto): 3.5 x the 90th percentile of |d 1/K0| between the brightest decile of the walks (commit order is
    seed-intensity order) and their nearest trace one isotope heavier at the walk's charge (partner m/z tolerance,
    delta_rt, iso_im_tol), clamped to [0.01, 0.03]. `pcs` must be inferPrecursors_'s output in its order (the veto
    edits charges in place); `v2in` is the walk record inferPrecursors_ returned beside it.

    @param[in] ms1 The MS1 traces
    @param[in] ms1st Their store
    @param[in,out] pcs The precursors; mono_mz may change
    @param[in] v2in The walk record
    @param[in] delta_rt RT tolerance
    @param[in] iso_im_tol Isotope-partner 1/K0 tolerance
    @param[in] mass_ppm Mass tolerance (ppm)
    @throws Exception::Precondition if the walk record does not match the precursor list
  */
  void monoPositionPass_(const vector<Trace>& ms1,
                         const TraceStore& ms1st,
                         vector<PrecursorCand>& pcs,
                         const vector<MonoV2In>& v2in,
                         double delta_rt,
                         double iso_im_tol,
                         double mass_ppm)
  {
    const double PROTON_ = 1.00727646;
    // walks tested, moved right, steps rejected, extensions tested, moved left, walks read at a veto-changed charge,
    // moved back onto the seed
    std::array<Size, 7> mc_tot {};
    double v2_box = cfg_mono_pos_im_box, v2_p90 = -1.0;
    Size v2_nb = 0, v2_nhit = 0;
    if (v2in.size() != pcs.size())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[monocheck-v2] the walk record does not match the precursor list");
    }
    // The m/z-sorted index of inferPrecursors_ (same comparator on the same container, so the same order; the plain
    // std searches return the positions its bucket index returns).
    const Size n_total = ms1.size();
    vector<Size> by_mz;
    mzOrder_(ms1, by_mz); // [outer] the same helper, so the same order
    vector<double> sorted_mz(n_total), sorted_rt(n_total), sorted_im(n_total);
    vector<float> inten(n_total);
#pragma omp parallel for schedule(static)
    for (long ii = 0; ii < (long)n_total; ++ii)
    {
      const Size i = (Size)ii;
      const Trace& t = ms1[by_mz[i]];
      sorted_mz[i] = t.mz;
      sorted_rt[i] = rtOf(ms1st, t);
      sorted_im[i] = t.im;
      inten[i] = (float)intensityOf(ms1st, ms1[i]);
    }
    auto lbmz = [&](double x) -> Size { return (Size)(std::lower_bound(sorted_mz.begin(), sorted_mz.end(), x) - sorted_mz.begin()); };
    auto ubmz = [&](double x) -> Size { return (Size)(std::upper_bound(sorted_mz.begin(), sorted_mz.end(), x) - sorted_mz.begin()); };
    if (v2_box <= 0.0)
    {
      Size nz = 0;
      for (const MonoV2In& v : v2in)
      {
        if (v.z > 0) { ++nz; }
      }
      v2_nb = (nz + 9) / 10;
      vector<Size> idx;
      idx.reserve(v2_nb);
      for (Size i = 0; i < pcs.size() && idx.size() < v2_nb; ++i)
      {
        if (v2in[i].z > 0) { idx.push_back(i); }
      }
      vector<double> gaps(idx.size(), -1.0);
#pragma omp parallel for schedule(dynamic, 256)
      for (long ii = 0; ii < (long)idx.size(); ++ii)
      {
        const Size seed = pcs[idx[(Size)ii]].trace_idx;
        const Trace& s = ms1[seed];
        const double s_rt = rtOf(ms1st, s), tol = s.mz * mass_ppm * 1e-6;
        const double target = s.mz + ISO_STEP / v2in[idx[(Size)ii]].z;
        const Size klo = lbmz(target - tol), khi = ubmz(target + tol);
        double g = -1.0;
        for (Size q = klo; q < khi; ++q)
        {
          if (fabs(sorted_rt[q] - s_rt) > delta_rt || by_mz[q] == seed) { continue; }
          const double d = fabs(sorted_im[q] - s.im);
          if (d <= iso_im_tol && (g < 0.0 || d < g)) { g = d; }
        }
        gaps[(Size)ii] = g;
      }
      gaps.erase(std::remove_if(gaps.begin(), gaps.end(), [](double g) { return g < 0.0; }), gaps.end());
      v2_nhit = gaps.size();
      if (! gaps.empty())
      {
        const Size k = (Size)(0.9 * (double)(gaps.size() - 1));
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
      const MonoV2In& v = v2in[(Size)ii];
      PrecursorCand& pc = pcs[(Size)ii];
      // The final charge: the veto re-calls a z=1 walk 2+ (a multiple of the walk's charge) or unsets it (0: the
      // precursor becomes guessed and is dropped later, untested).
      if (v.z <= 0 || pc.charge <= 0 || pc.charge % v.z != 0) { continue; }
      const int z = pc.charge, mul = pc.charge / v.z, nl = v.nl * mul;
      const Size seed = pc.trace_idx;
      const Trace& s = ms1[seed];
      const double s_rt = rtOf(ms1st, s), tol = s.mz * mass_ppm * 1e-6;
      auto iso_sum = [&](int k, long* best) -> double {
        const double target = s.mz + k * ISO_STEP / z;
        const Size klo = lbmz(target - tol), khi = ubmz(target + tol);
        double sum = 0.0;
        float bi = -1.0f;
        for (Size q = klo; q < khi; ++q)
        {
          if (fabs(sorted_rt[q] - s_rt) > delta_rt || fabs(sorted_im[q] - s.im) > v2_box) { continue; }
          const Size j = by_mz[q];
          sum += inten[j];
          if (best && inten[j] > bi)
          {
            bi = inten[j];
            *best = (long)j;
          }
        }
        return sum;
      };
      // q of the step onto position -j: PS(-j) / PS(-j+1) x lambda(neutral mass of the mono after that step). A walk
      // re-read at a veto charge never stood on its in-between positions: an empty one there is no evidence against
      // the step, so the step is compared over two positions, PS(-j) / PS(-j+2) x lambda^2 / 2 (Poisson M+2/M0).
      auto ratio = [&](int j, long* best) -> double {
        const double neutral = (s.mz - j * ISO_STEP / z) * z - z * PROTON_;
        const double lam = std::max(0.10, 0.000594 * neutral);
        const double den = iso_sum(-j + 1, nullptr);
        if (den > 0.0) { return iso_sum(-j, best) / den * lam; }
        if (mul > 1 && (j - 1) % mul != 0)
        {
          const double den2 = iso_sum(-j + 2, nullptr);
          return den2 > 0.0 ? iso_sum(-j, best) / den2 * lam * lam / 2.0 : 0.0;
        }
        return 0.0;
      };
      // the m/z standing at position -i: the walk's own partner where the walk stepped, else the most intense trace
      // there in the box, else the position itself
      auto mz_left = [&](int i) -> double {
        if (i % mul == 0) { return v.left_mz[(Size)(i / mul - 1)]; }
        long bt = -1;
        iso_sum(-i, &bt);
        return bt >= 0 ? ms1[(Size)bt].mz : s.mz - i * ISO_STEP / z;
      };
      if (cfg_mono_check == 3)
      {
        // [monocheck-v3] every position c in [-(nl+1), 0] with intensity (the seed is isotope -c): cosine of PS(c-1 .. c+5)
        // against [0, A0 .. A5] (Poisson averagine at the mono's mass) x the averagine fraction observed; best c wins
        double best_s = -1.0, best_mz = v.mono;
        int best_c = 1;
        for (int c = -(nl + 1); c <= 0; ++c)
        {
          long bt = -1;
          const double pc0 = c == 0 ? 1.0 : iso_sum(c, &bt);
          if (pc0 <= 0.0) { continue; }
          const double neutral = (s.mz + c * ISO_STEP / z) * z - z * PROTON_;
          const double lam = std::max(0.10, 0.000594 * neutral);
          double sab = 0, saa = 0, sbb = 0, seen = 0, all = 0, a = pmath::exp(-lam);
          for (int i = -1; i <= 5; ++i)
          {
            const double t = i < 0 ? 0.0 : a, b = iso_sum(c + i, nullptr);
            if (i >= 0)
            {
              all += t;
              if (b > 0.0) { seen += t; }
              a *= lam / (i + 1);
            }
            sab += t * b;
            saa += t * t;
            sbb += b * b;
          }
          const double sc = sbb > 0.0 ? sab / std::sqrt(saa * sbb) * seen / all : 0.0;
          if (sc > best_s)
          {
            best_s = sc;
            best_c = c;
            best_mz = c == 0                        ? s.mz
                      : (-c <= nl && -c % mul == 0) ? v.left_mz[(Size)(-c / mul - 1)]
                                                    : (bt >= 0 ? ms1[(Size)bt].mz : s.mz + c * ISO_STEP / z);
          }
        }
        uint8_t f = (mul > 1 ? 32 : 0) | 1;
        if (best_c <= 0 && best_mz != v.mono)
        {
          if (best_mz > v.mono)
          {
            f |= 2;
            if (best_c == 0) { f |= 4; }
          }
          else
          {
            f |= 8 | 16;
            pc.mono_ext = -best_c > nl;
          }
          pc.mono_mz = best_mz;
        }
        flag[(Size)ii] = f;
        continue;
      }
      double mono = v.mono;
      uint8_t f = mul > 1 ? 32 : 0, r = 0;
      int acc = nl;
      if (nl > 0)
      {
        f |= 1;
        while (acc > 0 && ratio(acc, nullptr) < cfg_mono_pos_min_ratio)
        {
          --acc;
          ++r;
          if (! cfg_mono_pos_iterate)
          {
            // Leftmost step only, at a veto charge: an EMPTY in-between position cannot be the monoisotope (mz_left
            // would report a bare position), so the mono moves on to the walk's next own position. mul == 1: no-op.
            while (acc > 0 && acc % mul != 0 && iso_sum(-acc, nullptr) <= 0.0)
            {
              --acc;
            }
            break;
          }
        }
        if (acc < nl)
        {
          const double acc_mz = acc > 0 ? mz_left(acc) : s.mz;
          if (acc_mz > mono)
          {
            mono = acc_mz;
            f |= 2;
            if (acc == 0) { f |= 4; }
          }
        }
      }
      // B: only when A rejected nothing and the walk ended on a miss (< 5 steps).
      if (cfg_mono_pos_extend > 0.0 && acc == nl && v.nl < 5)
      {
        f |= 8;
        long bt = -1;
        const double q = ratio(nl + 1, &bt);
        if (q >= cfg_mono_pos_extend && q > 0.0 && bt >= 0)
        {
          mono = ms1[(Size)bt].mz;
          f |= 16;
          pc.mono_ext = true;
        }
      }
      pc.mono_mz = mono;
      flag[(Size)ii] = f;
      rej[(Size)ii] = r;
    }
    for (Size i = 0; i < pcs.size(); ++i)
    {
      const uint8_t f = flag[i];
      mc_tot[0] += f & 1;
      mc_tot[1] += (f >> 1) & 1;
      mc_tot[6] += (f >> 2) & 1;
      mc_tot[3] += (f >> 3) & 1;
      mc_tot[4] += (f >> 4) & 1;
      mc_tot[5] += (f >> 5) & 1;
      mc_tot[2] += rej[i];
    }
    String box = String(cfg_mono_pos_im_box);
    if (cfg_mono_pos_im_box <= 0.0)
    {
      char b[160];
      std::snprintf(b, sizeof b, "auto %.4f (3.5 x p90 %.4f of the +1-isotope 1/K0 gap, %zu of %zu brightest seeds)", v2_box, v2_p90, v2_nhit, v2_nb);
      box = String(b);
    }
    const String cfg = cfg_mono_check == 3
                         ? String("v3 (averagine cosine x observed fraction, positions -(walk+1) .. 0), im_box ") + box + ", at the final charge"
                         : "v2, min_ratio " + String(cfg_mono_pos_min_ratio) + ", im_box " + box + ", extend " + String(cfg_mono_pos_extend)
                             + (cfg_mono_pos_iterate ? ", iterate" : ", leftmost step only") + ", at the final charge";
    cfg_mono_check_stamp = cfg + ": " + String(mc_tot[1]) + " of " + String(mc_tot[0]) + " monos moved right, " + String(mc_tot[4]) + " of "
                           + String(mc_tot[3]) + " moved left";
    writeLogInfo_("[monocheck] " + cfg + ": " + String(mc_tot[0]) + " count walks with a leftward partner tested, " + String(mc_tot[2])
                  + " steps rejected, " + String(mc_tot[1]) + " monoisotopes moved right (" + String(mc_tot[6]) + " back onto the seed); "
                  + String(mc_tot[3]) + " walks tested for extension, " + String(mc_tot[4]) + " monoisotopes moved left (unclaimed); "
                  + String(mc_tot[5]) + " walks read at the charge the IM veto gave them; claims unchanged");
  }

  /**
    @brief [soa2] The fragment RT gate field.

    99.4% of scoring visits fail the RT gate, so `q` quantises RT into
    buckets never narrower than delta_rt and a 2-byte test precedes the unchanged exact one. The reject
    admits +/-2 buckets where +/-1 suffices in exact arithmetic: insurance against a floor() rounding
    across a boundary, which would silently DROP a candidate.
  */
  struct FragRt
  {
    vector<double> rt;
    vector<uint16_t> q;
    double rt0 = 0.0, inv = 0.0;
    /**
      @brief Exact RT of fragment @p i.

      @param[in] i Fragment index
      @return Its RT
    */
    double operator[](Size i) const
    { return rt[i]; }
    /**
      @brief The RT bucket of a value.

      Signed at the low end and clamped at both: a precursor may sit outside the window's
      fragment RT range. The test is written `!(x > 0.0)` rather than `x <= 0.0` so that a NaN
      coordinate lands in bucket 0 instead of reaching (int)NaN, which is undefined behaviour.

      @param[in] v Retention time
      @return Bucket in [0, 65535]
    */
    int bucket_of(double v) const
    {
      const double x = (v - rt0) * inv;
      return ! (x > 0.0) ? 0 : (x >= 65535.0 ? 65535 : (int)x);
    }
    /**
      @brief Fill the exact RTs and their buckets for a window's fragments.

      inv == 0 is the DISABLED state, not an error: every bucket is then 0, every dq is 0, and
      the exact test alone decides. That is the fallback for every degenerate parameter below,
      and it costs no branch in the scan.

      @param[in] fr The fragment traces
      @param[in] st Their store
      @param[in] delta_rt RT tolerance (smallest bucket width)
      @param[in] nt Number of tasks
    */
    void build(const vector<Trace>& fr, const TraceStore& st, double delta_rt, long nt)
    {
      const Size n = fr.size();
      rt.resize(n);
      q.assign(n, 0);
      rt0 = 0.0;
      inv = 0.0;
      // [serial] As nt tasks. The serial fold takes the first value, then any value strictly below (above) it, so
      // chunk 0 folds serially and every later chunk contributes its smallest (largest) non-NaN value, earliest on
      // ties, folded in chunk order: the same lo and hi for any values (a NaN first value sticks, a later NaN never
      // replaces).
      vector<double> clo(nt, 0.0), chi(nt, 0.0);
      vector<char> has(nt, 0);
      forTasks(nt, [&](long c) {
        double l = 0.0, h = 0.0;
        bool any = false;
        for (Size i = (Size)c * n / (Size)nt; i < (Size)(c + 1) * n / (Size)nt; ++i)
        {
          rt[i] = rtOf(st, fr[i]);
          if (c == 0)
          {
            if (! i || rt[i] < l) { l = rt[i]; }
            if (! i || rt[i] > h) { h = rt[i]; }
          }
          else if (rt[i] == rt[i])
          {
            if (! any || rt[i] < l) { l = rt[i]; }
            if (! any || rt[i] > h) { h = rt[i]; }
            any = true;
          }
        }
        clo[c] = l;
        chi[c] = h;
        has[c] = any;
      });
      double lo = clo[0], hi = chi[0];
      for (long c = 1; c < nt; ++c)
      {
        if (has[c])
        {
          if (clo[c] < lo) { lo = clo[c]; }
          if (chi[c] > hi) { hi = chi[c]; }
        }
      }
      // A denormal delta_rt makes 1/delta_rt infinite, and an RT range spanning +/-DBL_MAX makes
      // the span infinite; both then produce inf*0 = NaN inside bucket_of. Neither survives real
      // acquisition, but both are reachable, so each degrades to disabled.
      if (! (delta_rt > 0.0) || ! std::isfinite(delta_rt)) { return; }
      if (! std::isfinite(lo) || ! std::isfinite(hi)) { return; }
      const double span = hi - lo;
      if (! std::isfinite(span)) { return; }
      double iv = 1.0 / delta_rt;
      // Buckets are delta_rt wide, but widened if the window's RT span would overflow uint16.
      // Wider is always safe: it can only put more fragments in reach of the exact test.
      if (span > 0.0 && span * iv > 65000.0) { iv = 65000.0 / span; }
      if (! std::isfinite(iv) || ! (iv > 0.0)) { return; }
      rt0 = lo;
      inv = iv;
      forTasks(nt, [&](long c) {
        for (Size i = (Size)c * n / (Size)nt; i < (Size)(c + 1) * n / (Size)nt; ++i)
        {
          q[i] = (uint16_t)bucket_of(rt[i]);
        }
      });
    }
  };

  /**
    @brief [foreign] One window's precursors, the candidate foreign envelopes of its spectra: every precursor whose mono
    m/z lies in the window (a precursor whose mono lies just below it is not a candidate: extending the candidates there
    moved no peptide on the screen panel), with a finite apex RT and 1/K0 and an isotope-supported charge (a guessed precursor's
    charge is a default, so its envelope would be fictional; they are dropped anyway),
    sorted by apex RT (ties by index).

    Built from the run's precursor list, not from the tile's owned subset, so it
    is the same in any tiling.
  */
  struct ForeignIdx_
  {
    const vector<PrecursorCand>* p = nullptr;
    vector<double> rt;    ///< apex RT, ascending
    vector<uint32_t> idx; ///< into *p, in rt order
    /**
      @brief Index the candidates among v[lo, hi) with a finite RT and 1/K0 and a charge > 0, by apex RT.

      @param[in] v The run's precursors (referenced, must outlive the index)
      @param[in] lo First index
      @param[in] hi End index (exclusive)
    */
    void build(const vector<PrecursorCand>& v, Size lo, Size hi)
    {
      p = &v;
      vector<pair<double, uint32_t>> o;
      o.reserve(hi > lo ? hi - lo : 0);
      for (Size i = lo; i < hi; ++i)
      {
        if (std::isfinite(v[i].rt) && std::isfinite(v[i].im) && v[i].charge > 0) { o.emplace_back(v[i].rt, (uint32_t)i); }
      }
      std::sort(o.begin(), o.end());
      rt.resize(o.size());
      idx.resize(o.size());
      for (Size i = 0; i < o.size(); ++i)
      {
        rt[i] = o[i].first;
        idx[i] = o[i].second;
      }
    }
  };

  /**
    @brief [foreign] assembly:foreign_precursor_weight: attenuate, among one precursor's admitted fragments, the residual
    (unfragmented) envelopes of the OTHER precursors isolated with it -- apex within drt of this one's, 1/K0 within
    cfg_foreign_dim, neither a twin (mono within 0.01) nor this ion called again (either mono on the other's isotope
    grid, up to 3 steps, within max(cfg_foreign_ppm, 30 ppm)).

    A candidate's positions are mono + k ISO/z, k = 0..3; a
    position is present when an admitted fragment lies within cfg_foreign_ppm of it. With at least cfg_foreign_min_iso
    positions present, those fragments are multiplied by cfg_foreign_w (0 removes them) and their cap key is recomputed
    from the attenuated IM-weighted apex (@p frag_ai, attenuated with it). A fragment within max(cfg_foreign_ppm, 30 ppm) of this
    precursor's own envelope (mono to M+3) neither counts nor is touched. Two matched positions are weak evidence for a 1+
    (or same-charge) fragment and its own isotope, which line up with a co-isolated envelope of that charge as soon as one
    does (open). Returns the fragments
    attenuated. @p frag_kq: the [fuse] key factors, parallel to @p frags when not empty, compacted with it.

    @param[in] pc This precursor
    @param[in] drt Apex RT tolerance
    @param[in] fx The window's foreign index
    @param[in,out] frags Admitted fragments (m/z, intensity)
    @param[in,out] frag_scores Their cap keys
    @param[in,out] frag_bin Their flight-time bins
    @param[in,out] frag_span Their window-frame spans
    @param[in,out] frag_ai Their IM-weighted apex intensities
    @param[in,out] frag_kq Their [fuse] key factors, or empty
    @return Number of fragments attenuated
  */
  Size attenuateForeign_(const PrecursorCand& pc,
                         double drt,
                         const ForeignIdx_& fx,
                         vector<pair<double, double>>& frags,
                         vector<double>& frag_scores,
                         vector<uint32_t>& frag_bin,
                         vector<pair<int, int>>& frag_span,
                         vector<double>& frag_ai,
                         vector<double>& frag_kq) const
  {
    // The two own-ion tests take a wider tolerance than the envelope match: an ion's residual peaks scatter tens of
    // ppm around the MS1 mono grid, so at cfg_foreign_ppm the same ion called again with a mono just outside the sibling
    // test would take the precursor's own M..M+3 peaks.
    const double ppm = cfg_foreign_ppm * 1e-6, same = std::max(ppm, 30e-6);
    auto& cnt = foreignCount();
    ++cnt[6];
    auto on_grid = [&](double a, double b, int z) { // a within `same` of b + j ISO/z, |j| <= 3
      for (int j = -3; j <= 3; ++j)
      {
        if (std::fabs(a - (b + j * ISO_STEP / z)) <= a * same) { return true; }
      }
      return false;
    };
    vector<const PrecursorCand*> cand;
    long long n_skip = 0; // added to the shared counter once per call
    const Size lo = std::lower_bound(fx.rt.begin(), fx.rt.end(), pc.rt - drt) - fx.rt.begin();
    for (Size j = lo; j < fx.rt.size() && fx.rt[j] <= pc.rt + drt; ++j)
    {
      const PrecursorCand& f = (*fx.p)[fx.idx[j]];
      if (&f == &pc || std::fabs(f.rt - pc.rt) > drt || std::fabs(f.im - pc.im) > cfg_foreign_dim) { continue; }
      if (std::fabs(f.mono_mz - pc.mono_mz) < 0.01 || (pc.charge > 0 && on_grid(f.mono_mz, pc.mono_mz, pc.charge))
          || on_grid(pc.mono_mz, f.mono_mz, f.charge))
      {
        ++n_skip;
        continue;
      }
      cand.push_back(&f);
    }
    if (n_skip) { cnt[4] += n_skip; }
    if (cand.empty()) { return 0; }
    ++cnt[0];
    vector<uint32_t> om(frags.size()); // the fragments by m/z, ties by index
    for (Size i = 0; i < om.size(); ++i)
    {
      om[i] = (uint32_t)i;
    }
    std::sort(om.begin(), om.end(),
              [&](uint32_t a, uint32_t b) { return frags[a].first != frags[b].first ? frags[a].first < frags[b].first : a < b; });
    auto own = [&](double mz) {
      if (pc.charge <= 0) { return false; }
      for (int k = 0; k <= 3; ++k)
      {
        const double m = pc.mono_mz + k * ISO_STEP / pc.charge;
        if (std::fabs(mz - m) <= m * same) { return true; }
      }
      return false;
    };
    vector<char> hit(frags.size(), 0);
    vector<uint32_t> at;
    long long n_env = 0;
    for (const PrecursorCand* f : cand)
    {
      at.clear();
      int present = 0;
      for (int k = 0; k <= 3; ++k)
      {
        const double m = f->mono_mz + k * ISO_STEP / f->charge, tol = m * ppm;
        bool any = false;
        auto it = std::lower_bound(om.begin(), om.end(), m - tol, [&](uint32_t i, double v) { return frags[i].first < v; });
        for (; it != om.end() && frags[*it].first <= m + tol; ++it)
        {
          if (! own(frags[*it].first))
          {
            at.push_back(*it);
            any = true;
          }
        }
        present += any;
      }
      if (present < cfg_foreign_min_iso) { continue; }
      ++n_env;
      for (uint32_t i : at)
      {
        hit[i] = 1;
      }
    }
    Size n = 0;
    for (Size i = 0; i < hit.size(); ++i)
    {
      if (! hit[i]) { continue; }
      ++n;
      if (cfg_foreign_w <= 0.0) { continue; }
      frags[i].second *= cfg_foreign_w;
      const double log_key = pmath::log1p(frag_ai[i]); // the key is linear in its log
      frag_scores[i] = log_key > 0.0 ? frag_scores[i] * pmath::log1p(frag_ai[i] * cfg_foreign_w) / log_key : frag_scores[i];
      frag_ai[i] *= cfg_foreign_w; // [fuse] a fused key sums the attenuated apex
    }
    if (n && cfg_foreign_w <= 0.0) // remove: compact every parallel array
    {
      Size w = 0;
      for (Size i = 0; i < hit.size(); ++i)
      {
        if (! hit[i])
        {
          frags[w] = frags[i];
          frag_scores[w] = frag_scores[i];
          frag_bin[w] = frag_bin[i];
          frag_span[w] = frag_span[i];
          frag_ai[w] = frag_ai[i];
          if (! frag_kq.empty()) { frag_kq[w] = frag_kq[i]; }
          ++w;
        }
      }
      frags.resize(w);
      frag_scores.resize(w);
      frag_bin.resize(w);
      frag_span.resize(w);
      frag_ai.resize(w);
      if (! frag_kq.empty()) { frag_kq.resize(w); }
    }
    cnt[1] += n_env;
    cnt[3] += (long long)n;
    if (n) { ++cnt[2]; }
    return n;
  }

  /**
    @brief [fuse] assembly:fragment_fuse_ppm: one peak per fragment ion.

    The IM picker keeps 1-2 flight-time bins of an ion
    that spans 4-7, so one ion becomes several traces and several peaks 5-20 ppm apart; the same-bin merge reaches only
    exact bins. Groups the fragments within @p ppm (anchored: strongest first, each unfused fragment an anchor taking
    every unfused one within ppm of its own m/z). A fused peak: m/z = intensity-weighted mean of the members' (emitted) intensities,
    intensity = their sum, key = the best-keyed member's factor @p kq (its score
    and shared frames) times ln(1 + the summed IM-weighted apex @p ai), the formula
    of an unfused fragment's key applied to the whole ion; it is never below the best member's key (non-negative keys, up to rounding). It sits at
    its members' lowest index, so the relative order of the peaks is kept; a lone fragment is kept bit for bit. @p keys and are parallel to @p frags
    and compacted with it; @p ai and @p kq are parallel to it on entry and stale on return. Deterministic: every order has a total tie-break and the
    sums run in m/z order. Returns the groups of >= 2 fused and adds their members to @p n_members.

    @param[in,out] frags Fragments (m/z, intensity); compacted
    @param[in,out] keys Their keys; compacted
    @param[in] ai Their IM-weighted apex intensities
    @param[in] kq Their key factors
    @param[in] ppm Fusion tolerance (ppm)
    @param[in,out] n_members Incremented by the members of the fused groups
    @return Number of groups of >= 2 fused
  */
  static long long fuseFragments_(vector<pair<double, double>>& frags,
                                  vector<double>& keys,
                                  const vector<double>& ai,
                                  const vector<double>& kq,
                                  double ppm,
                                  long long& n_members)
  {
    const Size n = frags.size();
    if (n < 2 || ! (ppm > 0.0)) { return 0; }
    const double tol = ppm * 1e-6;
    constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max();
    vector<uint32_t> om(n); // by m/z, ties by index
    for (Size i = 0; i < n; ++i)
    {
      om[i] = (uint32_t)i;
    }
    std::sort(om.begin(), om.end(),
              [&](uint32_t a, uint32_t b) { return frags[a].first != frags[b].first ? frags[a].first < frags[b].first : a < b; });
    vector<uint32_t> grp(n, NONE);
    uint32_t ng = 0;
    vector<uint32_t> oi(om); // by intensity desc, then m/z asc, then index
    std::sort(oi.begin(), oi.end(), [&](uint32_t a, uint32_t b) {
      if (frags[a].second != frags[b].second) { return frags[a].second > frags[b].second; }
      if (frags[a].first != frags[b].first) { return frags[a].first < frags[b].first; }
      return a < b;
    });
    for (uint32_t a : oi)
    {
      if (grp[a] != NONE) { continue; }
      const double m = frags[a].first, d = m * tol;
      auto it = std::lower_bound(om.begin(), om.end(), m - d, [&](uint32_t i, double v) { return frags[i].first < v; });
      for (; it != om.end() && frags[*it].first <= m + d; ++it)
      {
        if (grp[*it] == NONE) { grp[*it] = ng; }
      }
      ++ng;
    }
    if (ng == n)
    {
      return 0; // nothing within ppm: unchanged
    }
    struct Acc
    {
      double w = 0.0, wmz = 0.0, ai = 0.0;
      uint32_t head = NONE, top = NONE, best = NONE, n = 0;
    };
    vector<Acc> acc(ng);
    for (uint32_t i : om) // m/z order: a fixed summation order
    {
      Acc& g = acc[grp[i]];
      g.w += frags[i].second;
      g.wmz += frags[i].second * frags[i].first;
      g.ai += ai[i];
      if (g.best == NONE || keys[i] > keys[g.best])
      {
        g.best = i; // the best key; ties: lower m/z
      }
      g.head = std::min(g.head, i);
      if (g.top == NONE || frags[i].second > frags[g.top].second)
      {
        g.top = i; // the strongest; ties: lower m/z
      }
      ++g.n;
    }
    long long groups = 0;
    Size w = 0;
    for (Size i = 0; i < n; ++i) // write position w <= i <= every unread member
    {
      const Acc& g = acc[grp[i]];
      if (g.head != i) { continue; }
      if (g.n == 1)
      {
        frags[w] = frags[i];
        keys[w] = keys[i];
      }
      else
      {
        const double mz = g.w > 0.0 ? g.wmz / g.w : frags[g.top].first;
        frags[w] = {mz, g.w};
        keys[w] = kq[g.best] * pmath::log1p(g.ai);
        ++groups;
        n_members += g.n;
      }
      ++w;
    }
    frags.resize(w);
    keys.resize(w);
    return groups;
  }

  /**
    @brief Per-precursor assembly: this precursor claims every fragment passing its gate (a fragment may
    be shared across precursors -> chimeric).

    Leaves @p out empty on failure.
    pdense is caller-owned scratch, so the caller decides the per-thread copy. Keeps the top max_frags
    by score, drops if < min_frags, annotates the synthetic precursor. Deterministic.

    @param[in] pc The precursor
    @param[in] win_lo Window lower m/z bound
    @param[in] win_hi Window upper m/z bound
    @param[in] frag_traces The window's fragment traces
    @param[in] wst Their store
    @param[in] frag_im Their 1/K0 values
    @param[in] frag_rt Their RT gate field
    @param[in] ms1_traces The MS1 traces
    @param[in] ms1st Their store
    @param[in] fg The window's fragment statistics
    @param[in] delta_im 1/K0 tolerance
    @param[in] delta_rt RT tolerance
    @param[in] min_corr Minimum correlation
    @param[in] min_corr_pts Minimum points for a correlation
    @param[in] min_frags Minimum number of fragments
    @param[in] max_frags Maximum number of fragments
    @param[in,out] pdense Scratch over window frames (zero on entry and on return)
    @param[out] out The pseudo-MS/MS spectrum, empty on failure
    @param[in] fx The window's foreign index, or nullptr
  */
  void assembleOne_(const PrecursorCand& pc,
                    double win_lo,
                    double win_hi,
                    const vector<Trace>& frag_traces,
                    const TraceStore& wst,
                    const vector<double>& frag_im,
                    const FragRt& frag_rt,
                    const vector<Trace>& ms1_traces,
                    const TraceStore& ms1st,
                    const FragStats& fg,
                    double delta_im,
                    double delta_rt,
                    double min_corr,
                    int min_corr_pts,
                    Size min_frags,
                    Size max_frags,
                    vector<float>& pdense,
                    MSSpectrum& out,
                    const ForeignIdx_* fx) const
  {
    vector<pair<double, double>> frags;
    vector<double> frag_scores;
    vector<double> frag_ai; // [foreign]/[fuse] each fragment's IM-weighted apex, its cap key's intensity (only with fx or fusion)
    vector<double> frag_kq; // [fuse] each fragment's cap key factor besides that intensity term (only with fusion)
    const bool fuse_on = cfg_fuse_ppm > 0.0;
    vector<uint32_t> frag_bin;        // [merge] the TOF bin of each fragment (0: OpenMS detector, no bin; a real bin 0, m/z just
                                      // below the usual 100 m/z floor, is compared by exact m/z like the OpenMS path)
    vector<pair<int, int>> frag_span; // [merge] its detected window-frame span
    const double n_grid = (double)fg.grid_frames;
    const Trace& p_tr = ms1_traces[pc.trace_idx];
    vector<int> touched;
    touched.reserve(p_tr.np());
    for (Size k = 0; k < p_tr.span(); ++k)
    {
      if (! real(ms1st, p_tr, k)) { continue; }
      const double pval = (double)xv(ms1st, p_tr, k);
      const int gi = (*fg.nearest_local)[p_tr.frame0 + k]; // MS1 frame -> window frame, precomputed
      if (gi >= 0)
      {
        if (pdense[gi] == 0.0f) { touched.push_back(gi); }
        pdense[gi] += (float)pval;
      }
    }
    double psum = 0, psumsq = 0;
    for (int gi : touched)
    {
      double v = pdense[gi];
      psum += v;
      psumsq += v * v;
    }
    int p_lo = std::numeric_limits<int>::max(), p_hi = -1; // [short] the precursor's span in window frames
    int p_apex = -1;                                       // [coloc] the precursor apex frame
    for (int gi : touched)
    {
      p_lo = std::min(p_lo, gi);
      p_hi = std::max(p_hi, gi);
      if (p_apex < 0 || pdense[gi] > pdense[p_apex]) { p_apex = gi; }
    }
    // [recall] assembly:im_weight_sigma: a Gaussian of the fragment's 1/K0 distance @p dim to the precursor
    auto im_weigh = [this](double& x, double dim) {
      if (cfg_im_weight_sigma > 0.0) { x *= pmath::exp(-(dim * dim) / (2.0 * cfg_im_weight_sigma * cfg_im_weight_sigma)); }
    };
    // [adaptive] with gate:adaptive faint, the relaxed modes apply only to faint precursors (short precursor traces)
    const bool faint = cfg_adaptive && (int)touched.size() <= cfg_faint_frames, all_p = ! cfg_adaptive || faint;
    const Size n_touched = touched.size(); // [diag] spectrum_meta, after pdense is reset
    const double p_top = p_apex >= 0 ? (double)pdense[p_apex] : 0.0;
    const bool g_shrink = cfg_shrink && all_p, g_coloc = cfg_coloc && all_p, g_conc = cfg_short_conc && all_p;
    const double g_drt = faint ? cfg_faint_drt : delta_rt, g_cp = faint ? cfg_faint_cp : cfg_corr_power;
    const double p_ctr = 0.5 * (p_lo + p_hi), p_half = 0.5 * (p_hi - p_lo);
    double pmean = psum / n_grid, pvar = psumsq - n_grid * pmean * pmean;
    double pinv = pvar > 0 ? 1.0 / sqrt(pvar) : 0.0;
    if (pinv > 0.0)
    {
      Size lo = lower_bound(frag_im.begin(), frag_im.end(), pc.im - delta_im) - frag_im.begin();
      Size hi = upper_bound(frag_im.begin(), frag_im.end(), pc.im + delta_im) - frag_im.begin();
      const int pq = frag_rt.bucket_of(pc.rt);
      for (Size fi = lo; fi < hi; ++fi)
      {
        const int dq = (int)frag_rt.q[fi] - pq;
        if (dq < -2 || dq > 2)
        {
          continue; // 2 B reject, see FragRt
        }
        // [soa] the gate reads a parallel array; the trace record is touched only after it passes.
        if (fabs(frag_rt[fi] - pc.rt) > g_drt)
        {
          continue; // RT gate ([adaptive] faint window)
        }
        if (fg.invnorm[fi] == 0.0)
        {
          continue; // degenerate fragment
        }
        const Trace& f = frag_traces[fi];
        if (fabs(f.mz - pc.mono_mz) < 0.01)
        {
          continue; // exclude precursor peak
        }
        double dot = 0;
        int overlap = 0;
        // overlap <= support, so a fragment with fewer support points than
        // min_corr_pts can never pass the overlap guard below: skip it before the dot product.
        // Output-identical by construction (the guard's continue is the only thing it would reach).
        const bool short_np = f.np() < (Size)std::max(min_corr_pts, 0);
        if (short_np && ! g_conc && ! (g_shrink && (int)f.np() >= cfg_shrink_min_pts)) { continue; }
        int f_lo = -1, f_hi = -1, f_apex = -1; // [short]/[coloc] the fragment's detected span and apex
        float f_max = 0.0f;
        // two sequential streams: the fragment's span in the arena and pdense over window frames
        for (Size k = 0; k < f.span(); ++k)
        {
          const float fv = xv(wst, f, k);
          if (fv == 0.0f) { continue; }
          if (f_lo < 0) { f_lo = (int)(f.frame0 + k); }
          f_hi = (int)(f.frame0 + k);
          if (fv > f_max)
          {
            f_max = fv;
            f_apex = (int)(f.frame0 + k);
          }
          const float pv = pdense[f.frame0 + k];
          if (pv != 0.0f)
          {
            dot += (double)pv * fv;
            ++overlap;
          }
        }
        if (g_coloc) // [coloc] DIA-Umpire-style co-location
        {
          if (f_lo < 0 || p_hi < 0) { continue; }
          if (f_apex < p_lo - 1 || f_apex > p_hi + 1 || p_apex < f_lo - 1 || p_apex > f_hi + 1) { continue; }
          const int inter = std::min(f_hi, p_hi) - std::max(f_lo, p_lo) + 1, shorter = std::min(f_hi - f_lo, p_hi - p_lo) + 1;
          if (inter < cfg_coloc_overlap * shorter) { continue; }
        }
        double c, c_raw = 0.0; // c_raw: the unshrunk score, for the cap key
        const bool no_corr = g_shrink ? overlap < cfg_shrink_min_pts : (short_np || overlap < min_corr_pts); // [short]/[shrink] precedence
        if (g_conc && no_corr)
        {
          // [short] too few points for a correlation: faint ions are the tips of their elution peaks, so their
          // detected spans shrink with intensity but stay centred on the apex. Keep a fragment whose span lies inside
          // the precursor's (one cycle of slack) with its centre near the precursor's; the score ranks by centring.
          if ((int)f.np() < cfg_short_min_pts || f_lo < 0 || p_hi < 0) { continue; }
          if (f_lo < p_lo - 1 || f_hi > p_hi + 1) { continue; }
          const double d = std::fabs(0.5 * (f_lo + f_hi) - p_ctr), tol = std::max(1.0, cfg_short_center_tol * p_half);
          if (d > tol) { continue; }
          c = min_corr * (1.0 - 0.5 * d / (tol + 1.0)); // at the gate: ranks and weighs as a borderline fragment
          c_raw = c;
        }
        else
        {
          if (overlap < (g_shrink ? cfg_shrink_min_pts : min_corr_pts))
          {
            continue; // overlap guard; [shrink] its own floor
          }
          c = (dot - n_grid * pmean * fg.mean[fi]) * pinv * fg.invnorm[fi]; // Pearson over the full RT grid
          c_raw = c;
          if (g_shrink) // [shrink] Fisher-z posterior mean
          {
            const double z = pmath::atanh(std::max(-0.999, std::min(0.999, c))), s2 = 1.0 / std::max(overlap - 3, 1);
            const double t2 = cfg_shrink_sd * cfg_shrink_sd, mu = pmath::atanh(std::max(-0.999, std::min(0.999, cfg_shrink_r0)));
            c = pmath::tanh((z / s2 + mu / t2) / (1.0 / s2 + 1.0 / t2));
          }
          if (! (c >= min_corr))
          {
            continue; // correlation gate (a NaN score is rejected too)
          }
        }
        // The IM weight feeds both the cap key and the emitted intensity; corr^corr_power feeds only the emitted one.
        double inten = intensityOf(wst, f);
        const double dim = frag_im[fi] - pc.im;
        im_weigh(inten, dim);
        double emit_inten = inten;
        if (g_cp > 0.0 && c > 0.0) { emit_inten *= pmath::pow(c, g_cp); }
        frags.emplace_back(exportMz_(f.tof, bOf(wst, f), f.mz), emit_inten);
        frag_bin.push_back(f.tof);
        frag_span.emplace_back(f_lo, f_hi);
        if (fx || fuse_on) { frag_ai.push_back(inten); }
        frag_scores.push_back(std::max(c_raw, 0.0) * pmath::log1p(inten) * (1.0 + overlap)); // [cap] the max_fragments ranking key
        if (fuse_on)                                                                         // [fuse] the same key without its intensity term
        {
          frag_kq.push_back(std::max(c_raw, 0.0) * (1.0 + overlap));
        }
      }
    }
    for (int gi : touched)
    {
      pdense[gi] = 0.0f; // reset scratch
    }
    // [foreign] other precursors' residual envelopes, before the merge and the cap: a removed fragment frees its slot
    const Size n_before = frags.size(); // [foreign] for the fell-below count
    const Size n_foreign = fx && ! frags.empty() ? attenuateForeign_(pc, g_drt, *fx, frags, frag_scores, frag_bin, frag_span, frag_ai, frag_kq) : 0;
    // [merge] guard: one TOF bin is one m/z. What toTof's same-frame merge cannot reach (slices beyond its span, valley-
    // split children, band edges, tile carries) is merged here, before the cap, into the strongest member: members whose
    // RT span overlaps its span are summed (parallel slices of one ion), RT-disjoint pieces are dropped (the stronger
    // piece stands, as the twin fold's max); the merged peak takes the strongest member's m/z and the best cap key of
    // the summed members (not recomputed). Two peaks at one m/z stall OpenMS' MatchedIterator.
    if (cfg_merge_same_bin && frags.size() > 1)
    {
      vector<Size> ord(frags.size());
      for (Size i = 0; i < ord.size(); ++i)
      {
        ord[i] = i;
      }
      auto key = [&](Size i) { return std::make_tuple(frag_bin[i] == 0, frag_bin[i], frag_bin[i] ? 0.0 : frags[i].first, i); };
      std::sort(ord.begin(), ord.end(), [&](Size a, Size b) { return key(a) < key(b); });
      auto same = [&](Size a, Size b) { return frag_bin[a] ? frag_bin[a] == frag_bin[b] : (! frag_bin[b] && frags[a].first == frags[b].first); };
      vector<char> gone(frags.size(), 0);
      long long nm = 0;
      for (Size g = 0; g < ord.size();)
      {
        Size e = g + 1;
        while (e < ord.size() && same(ord[g], ord[e]))
        {
          ++e;
        }
        if (e - g > 1)
        {
          const Size head = ord[g]; // lowest index: the output order is unchanged
          Size top = head;
          for (Size k = g; k < e; ++k)
          {
            if (frags[ord[k]].second > frags[top].second) { top = ord[k]; }
          }
          double sum = 0.0, cap = -std::numeric_limits<double>::infinity(), asum = 0.0;
          Size kb = top; // [fuse] the member whose key the peak takes
          for (Size k = g; k < e; ++k)
          {
            const Size i = ord[k];
            if (! (frag_span[i].second < frag_span[top].first || frag_span[top].second < frag_span[i].first))
            {
              sum += frags[i].second;
              if (frag_scores[i] > cap) { kb = i; }
              cap = std::max(cap, frag_scores[i]);
              if (fuse_on) { asum += frag_ai[i]; }
            }
            if (k > g) { gone[i] = 1; }
          }
          frags[head] = {frags[top].first, sum};
          frag_scores[head] = cap;
          if (fuse_on)
          {
            frag_ai[head] = asum;
            frag_kq[head] = frag_kq[kb];
          } // [fuse] the summed apex, the key's member's factor
          nm += (long long)(e - g - 1);
        }
        g = e;
      }
      if (nm)
      {
        Size w = 0;
        for (Size i = 0; i < frags.size(); ++i)
        {
          if (! gone[i])
          {
            frags[w] = frags[i];
            frag_scores[w] = frag_scores[i];
            if (fuse_on)
            {
              frag_ai[w] = frag_ai[i];
              frag_kq[w] = frag_kq[i];
            }
            ++w;
          }
        }
        frags.resize(w);
        frag_scores.resize(w);
        if (fuse_on)
        {
          frag_ai.resize(w);
          frag_kq.resize(w);
        }
        mergeCount()[2] += nm;
      }
    }
    // [fuse] one peak per fragment ion, after the merge (frag_bin, frag_span, frag_ai and frag_kq are stale from here on)
    // and before min_fragments and the cap: the cap ranks fused peaks by a key recomputed from the summed apex, so the
    // slots of an ion's pieces refill
    const Size n_merged = frags.size(); // [foreign]/[fuse] the fell-below counts
    if (fuse_on && n_merged > 0)
    {
      long long n_mem = 0;
      const long long n_grp = fuseFragments_(frags, frag_scores, frag_ai, frag_kq, cfg_fuse_ppm, n_mem);
      auto& fc = fuseCount();
      ++fc[0];
      fc[1] += (long long)n_merged;
      fc[2] += (long long)frags.size();
      if (n_grp)
      {
        fc[3] += n_grp;
        fc[4] += n_mem;
      }
      if (n_merged >= min_frags && frags.size() < min_frags) { ++fc[5]; }
    }
    if (frags.size() < min_frags)
    { // [foreign] counted only when the removal (weight 0) took a spectrum that had enough fragments below the floor
      if (n_foreign && cfg_foreign_w <= 0.0 && n_before >= min_frags && n_merged < min_frags) { ++foreignCount()[5]; }
      return;
    }
    if (frags.size() > max_frags)
    {
      vector<Size> idx(frags.size());
      for (Size i = 0; i < idx.size(); ++i)
      {
        idx[i] = i;
      }
      // Rank by correlation. Ties break on m/z, then index, so the output is deterministic.
      partial_sort(idx.begin(), idx.begin() + max_frags, idx.end(), [&](Size a, Size b) {
        const double ka = frag_scores[a];
        const double kb = frag_scores[b];
        if (ka != kb) { return ka > kb; }
        if (frags[a].first != frags[b].first) { return frags[a].first < frags[b].first; }
        return a < b;
      });
      vector<pair<double, double>> keep;
      keep.reserve(max_frags);
      for (Size i = 0; i < max_frags; ++i)
      {
        keep.push_back(frags[idx[i]]);
      }
      frags.swap(keep);
    }
    out.setMSLevel(2);
    out.setType(SpectrumSettings::SpectrumType::CENTROID);
    out.setRT(pc.rt);
    // [mem] exact reserve: no doubling slack held from emit until the write.
    out.reserve(frags.size());
    for (const auto& fr : frags)
    {
      out.emplace_back(fr.first, fr.second);
    }
    out.sortByPosition();
    OpenMS::Precursor prec;
    prec.setMZ(pc.mono_mz); // MS1 traces come from the OpenMS detector: no bin to calibrate
    if (pc.charge > 0) { prec.setCharge(pc.charge); }
    prec.setIsolationWindowLowerOffset(max(0.0, pc.mono_mz - win_lo));
    prec.setIsolationWindowUpperOffset(max(0.0, win_hi - pc.mono_mz));
    prec.setDriftTime(pc.im);
    prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    out.setPrecursors({prec});
    // [provenance] userParams: guessed (no isotope support; always 0, guessed precursors are dropped before assembly) vs supported, and the isotope
    // peaks behind the mono call (0/1: the mono may be off by +-k*1.00335 Da, a hint for open-search correction).
    out.setMetaValue("spx_guessed", 0);
    out.setMetaValue("spx_n_isotopes", pc.n_isotopes);
    if (cfg_twin_ext_mark && pc.mono_ext)
    {
      out.setMetaValue("spx_mono_ext", 1); // [twin-rt] read and removed by mergeTwins
    }
    if (cfg_spectrum_meta) // [diag] diag:spectrum_meta
    {
      out.setMetaValue("spx_p_frames", (int)n_touched);
      out.setMetaValue("spx_faint", (int)n_touched <= cfg_faint_frames ? 1 : 0);
      out.setMetaValue("spx_p_apex", p_top);
    }
  }

  /**
    @brief [arena] Assertions on the SHIPPED compactUnreferenced(): the compaction must survive traces
    whose container order differs from their arena-offset order.

    @return EXECUTION_OK if every assertion holds, UNEXPECTED_RESULT otherwise
  */
  ExitCodes selftestArena_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (! ok)
      {
        writeLogError_(String("FAIL: ") + what);
        ++fail;
      }
      else
      {
        writeLogInfo_(String("  ok  ") + what);
      }
    };
    // four spans A(4) B(3) C(5) D(2) at offsets 0, 4, 7, 12, every value distinct; span k of trace T
    // holds T*100 + k + 1, with a gap (zero) inside D so interior zeros are exercised
    auto build = [](TraceStore& st, vector<Trace>& tr) {
      st = TraceStore();
      tr.clear();
      const vector<pair<uint32_t, int>> spec = {{4, 1}, {3, 2}, {5, 3}, {3, 4}};
      for (const auto& [n_pts, tag] : spec)
      {
        vector<pair<uint32_t, float>> pts;
        for (uint32_t k = 0; k < n_pts; ++k)
        {
          if (! (tag == 4 && k == 1)) { pts.emplace_back(k, (float)(tag * 100 + k + 1)); }
        }
        tr.push_back(makeSpan(pts, st));
      }
    };
    auto bytes = [](const TraceStore& st, const Trace& t) { return vector<float>(st.inten.begin() + t.off, st.inten.begin() + t.off + t.len); };
    TraceStore st0;
    vector<Trace> tr0;
    build(st0, tr0);
    vector<vector<float>> ref;
    for (const auto& t : tr0)
    {
      ref.push_back(bytes(st0, t));
    }
    chk(st0.inten.size() == 15 && tr0[3].len == 3 && tr0[3].npts == 2, "fixture: 15 floats, D spans 3 frames with a gap");

    // 1. container order D,B,A,C (offsets 12,4,0,7), keep D,B,A, drop C. A container-order walk would
    //    move D to 0 and B to 2 before reading A from [0,4): A would come back as D's and B's bytes.
    {
      vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]};
      TraceStore s = st0;
      const Size freed = compactUnreferenced(v, s, vector<bool> {true, true, true, false});
      chk(freed == 1, "1: one profile released");
      chk(s.inten.size() == 10, "1: arena holds exactly the kept spans (4+3+3)");
      chk(bytes(s, v[0]) == ref[3], "1: D survives (interior zero kept)");
      chk(bytes(s, v[1]) == ref[1], "1: B survives");
      chk(bytes(s, v[2]) == ref[0], "1: A survives (a container-order walk overwrites it)");
      chk(v[3].len == 0 && v[3].npts == 0, "1: C freed (span cleared)");
    }
    // 2. offset order, everything needed: unchanged bytes, no move
    {
      vector<Trace> v = tr0;
      TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, true)) == 0, "2: nothing released");
      bool ok = s.inten.size() == 15;
      for (Size i = 0; i < 4; ++i)
      {
        ok = ok && bytes(s, v[i]) == ref[i];
      }
      chk(ok, "2: offset-order input is a no-op");
    }
    // 3. the freed trace first, then a trimmed span: the trim moved off/len/apex (as trimToSpan does)
    {
      vector<Trace> v = {tr0[2], tr0[0], tr0[1], tr0[3]};
      TraceStore s = st0;
      v[1].off += 1;
      v[1].len -= 2;
      v[1].frame0 += 1;
      v[1].npts = 2; // A trimmed to [1,3)
      const vector<float> a_trim(ref[0].begin() + 1, ref[0].begin() + 3);
      chk(compactUnreferenced(v, s, vector<bool> {false, true, true, true}) == 1, "3: C released first");
      chk(bytes(s, v[1]) == a_trim && bytes(s, v[2]) == ref[1] && bytes(s, v[3]) == ref[3], "3: trimmed A, B, D survive");
      chk(s.inten.size() == 2 + 3 + 3, "3: arena sized to the trimmed spans");
    }
    // 4. nobody needed / a needed trace whose span was already released (len == 0) / duplicate references
    {
      vector<Trace> v = tr0;
      TraceStore s = st0;
      chk(compactUnreferenced(v, s, vector<bool>(4, false)) == 4 && s.inten.empty(), "4a: none needed -> empty arena");
      vector<Trace> w = tr0;
      TraceStore s2 = st0;
      w[0].freeProfile();
      chk(compactUnreferenced(w, s2, vector<bool>(4, true)) == 0 && s2.inten.size() == 11 && bytes(s2, w[1]) == ref[1],
          "4b: a needed len==0 trace owns no bytes");
    }
    // 5. compaction is idempotent
    {
      vector<Trace> v = {tr0[3], tr0[1], tr0[0], tr0[2]};
      TraceStore s = st0;
      compactUnreferenced(v, s, vector<bool> {true, true, true, false});
      const vector<float> once = s.inten;
      compactUnreferenced(v, s, vector<bool> {true, true, true, false});
      chk(s.inten == once && bytes(s, v[2]) == ref[0], "5: second compaction is a no-op");
    }
    // 6. overlapping spans are refused before any byte moves
    {
      vector<Trace> v = tr0;
      TraceStore s = st0;
      v[1].off = 2; // B now overlaps A
      bool threw = false;
      try
      {
        compactUnreferenced(v, s, vector<bool>(4, true));
      }
      catch (const OpenMS::Exception::Precondition&)
      {
        threw = true;
      }
      chk(threw, "6: overlapping spans throw Precondition");
    }

    writeLogInfo_(fail ? "[arena] SELFTEST FAILED" : "[arena] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /**
    @brief [twins] Assertions on mergeTwins with hand-built spectra: each side of every clause of the relation,
    the cluster closure, the survivor order and the fusion.

    @return EXECUTION_OK if every assertion holds, UNEXPECTED_RESULT otherwise
  */
  ExitCodes selftestTwins_()
  {
    int fail = 0;
    auto chk = [&](bool ok, const char* what) {
      if (! ok)
      {
        writeLogError_(String("FAIL: ") + what);
        ++fail;
      }
      else
      {
        writeLogInfo_(String("  ok  ") + what);
      }
    };
    auto sp = [](double rt, double mz, int z, double im, int niso, const vector<pair<double, float>>& pk) {
      MSSpectrum s;
      s.setMSLevel(2);
      s.setRT(rt);
      OpenMS::Precursor p;
      p.setMZ(mz);
      if (z > 0) { p.setCharge(z); }
      p.setDriftTime(im);
      s.setPrecursors({p});
      s.setMetaValue("spx_n_isotopes", niso);
      for (const auto& q : pk)
      {
        s.emplace_back(q.first, q.second);
      }
      return s;
    };
    const vector<pair<double, float>> peaks_in {{100.0, 10.0f}, {200.0, 20.0f}};
    auto pair_folds = [&](MSSpectrum b) {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), std::move(b)};
      return mergeTwins(v, 0.005) == 1 && v.size() == 1;
    };

    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, {{100.0, 10.0f}, {200.0, 20.0f}}),
                            sp(10.0, 500.001, 2, 1.0, 3, {{200.0005, 5.0f}, {300.0, 30.0f}})};
      const Size gone = mergeTwins(v, 0.005);
      chk(gone == 1 && v.size() == 1, "a twin pair folds into one spectrum");
      chk(v.size() == 1 && (int)v[0].getMetaValue("spx_n_isotopes") == 3 && v[0].getPrecursors()[0].getMZ() == 500.001,
          "the member with more isotopes survives, with its own precursor");
      chk(v.size() == 1 && v[0].size() == 3 && v[0][0].getMZ() == 100.0 && v[0][2].getMZ() == 300.0, "the survivor takes the union of the fragments");
      chk(v.size() == 1 && v[0].size() == 3 && fabs(v[0][1].getMZ() - (200.0 * 20.0 + 200.0005 * 5.0) / 25.0) < 1e-9
            && v[0][1].getIntensity() == 20.0f,
          "peaks within 10 ppm fuse to the intensity-weighted m/z at the maximum intensity");
    }
    chk(! pair_folds(sp(10.0, 500.001, 2, 1.006, 2, peaks_in)), "1/K0 beyond the tolerance: untouched");
    chk(pair_folds(sp(10.0, 500.001, 2, 1.004, 2, peaks_in)), "1/K0 inside the tolerance: folded");
    chk(! pair_folds(sp(11.1, 500.0, 2, 1.0, 2, peaks_in)), "the next MS1 frame: untouched");
    // [twin-rt] the near-frame relation (assembly:twin_rt_frames): |dRT| <= rt_tol, inside one RT cell
    // across frames a pair must sit in one flight-time bin (0.05 ppm) and one member must carry spx_mono_ext
    {
      auto ext = [&](MSSpectrum s) {
        s.setMetaValue("spx_mono_ext", 1);
        return s;
      };
      auto near_fold = [&](double rt_b, double tol, const vector<double>* cells, std::array<Size, 3>* st = nullptr) {
        vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), ext(sp(rt_b, 500.00002, 2, 1.0, 3, peaks_in))};
        return mergeTwins(v, 0.005, tol, cells, 3.0, st) == 1 && v.size() == 1;
      };
      const vector<double> edges {20.0}; // cell 0 owns rt < 17, cell 1 rt >= 17 (the first edge T with (T - rt) > 3)
      std::array<Size, 3> st {};
      chk(near_fold(11.1, 1.5, nullptr, &st) && st[0] == 0 && st[1] == 1, "the next MS1 frame inside rt_tol: folded, counted as a near-frame link");
      chk(near_fold(11.5, 1.5, nullptr), "|dRT| exactly rt_tol: folded (the bound is inclusive)");
      chk(! near_fold(11.75, 1.5, nullptr), "|dRT| beyond rt_tol: untouched");
      chk(! near_fold(11.1, 0.0, nullptr), "rt_tol 0: the same-frame relation, the next frame untouched");
      chk(near_fold(16.9, 7.0, &edges), "both inside one cell: dRT 6.9 within a 7-s tolerance folds");
      {
        std::array<Size, 3> se {};
        vector<MSSpectrum> v {sp(16.5, 500.0, 2, 1.0, 2, peaks_in), ext(sp(17.5, 500.00002, 2, 1.0, 3, peaks_in))};
        chk(mergeTwins(v, 0.005, 1.65, &edges, 3.0, &se) == 0 && v.size() == 2 && se[2] == 1,
            "a pair across a cell edge stays two spectra (counted)");
      }
      {
        vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.002, 2, 1.0, 2, peaks_in),
                              ext(sp(11.1, 500.00002, 2, 1.004, 2, peaks_in)), sp(12.2, 500.0, 2, 1.0, 2, peaks_in)};
        std::array<Size, 3> sc {};
        chk(mergeTwins(v, 0.005, 1.2, nullptr, 3.0, &sc) == 3 && v.size() == 1 && sc[0] == 1,
            "frames 10, 11.1, 12.2 chain into one cluster (span 2.2 > rt_tol); one pair is linked at one frame");
      }
      {
        vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), ext(sp(11.1, 500.0, 2, 1.006, 2, peaks_in)),
                              ext(sp(11.1, 500.0, 3, 1.0, 2, peaks_in)), ext(sp(11.1, 500.0125, 2, 1.0, 2, peaks_in))};
        chk(mergeTwins(v, 0.005, 1.65) == 0 && v.size() == 4, "near frames: 1/K0, charge and m/z still apply");
      }
      {
        vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), ext(sp(11.1, 500.0002, 2, 1.0, 2, peaks_in)),
                              sp(10.0, 500.0005, 2, 1.0, 2, peaks_in)};
        chk(mergeTwins(v, 0.005, 1.65) == 1 && v.size() == 2,
            "across frames only m/z within 0.05 ppm links (not 11.1/500.0002 at 0.4 ppm); one frame keeps the 20-ppm relation");
      }
      {
        vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(11.1, 500.00002, 2, 1.0, 2, peaks_in), sp(20.0, 600.0, 2, 1.0, 2, peaks_in),
                              ext(sp(21.1, 600.00002, 2, 1.0, 2, peaks_in))};
        chk(mergeTwins(v, 0.005, 1.65) == 1 && v.size() == 3 && v[0].getRT() == 10.0 && v[1].getRT() == 11.1,
            "a one-bin pair across frames folds only when one member's mono came from the v2 extension");
        bool clean = true;
        for (const auto& q : v)
        {
          clean = clean && ! q.metaValueExists("spx_mono_ext");
        }
        chk(clean, "the spx_mono_ext marks are removed (never written)");
      }
      {
        vector<MSSpectrum> v {ext(sp(12.0, 500.00002, 2, 1.0, 2, peaks_in)), sp(10.0, 500.0, 2, 1.0, 3, peaks_in)};
        mergeTwins(v, 0.005, 2.5);
        chk(v.size() == 1 && v[0].getRT() == 10.0 && v[0].getPrecursors()[0].getMZ() == 500.0,
            "near frames: the survivor rule is the fold's (more isotopes), with its own RT");
      }
    }
    chk(! pair_folds(sp(10.0, 500.0, 3, 1.0, 2, peaks_in)), "another charge: untouched");
    chk(! pair_folds(sp(10.0, 500.0125, 2, 1.0, 2, peaks_in)), "25 ppm apart: untouched");
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 0, 1.0, 2, peaks_in), sp(10.0, 500.0, 0, 1.0, 2, peaks_in)};
      chk(mergeTwins(v, 0.005) == 0 && v.size() == 2, "charge 0 (no charge call): untouched");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.0075, 2, 1.0, 2, peaks_in), sp(10.0, 500.015, 2, 1.0, 2, peaks_in)};
      chk(mergeTwins(v, 0.005) == 2 && v.size() == 1, "a chain 15 + 15 ppm is one cluster although its ends are 30 ppm apart");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.0075, 2, 1.0, 3, peaks_in), sp(10.0, 500.015, 2, 1.0, 2, peaks_in)};
      chk(mergeTwins(v, 0.005) == 2 && v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.0075,
          "a chain whose survivor is its middle: both ends are within 20 ppm of it and fold");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.001, 2, 1.0, 2, peaks_in)};
      mergeTwins(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.0, "equal isotopes and peaks: the earliest survives");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.001, 2, 1.0, 2, {{100.0, 1.0f}, {150.0, 1.0f}, {250.0, 1.0f}})};
      mergeTwins(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.001, "equal isotopes: more peaks survives");
    }
    {
      vector<MSSpectrum> v {sp(5.0, 300.0, 2, 0.8, 2, peaks_in), sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 700.0, 3, 1.1, 2, peaks_in),
                            sp(10.0, 500.001, 2, 1.0, 3, peaks_in), sp(12.0, 400.0, 2, 0.9, 2, peaks_in)};
      mergeTwins(v, 0.005);
      chk(v.size() == 4 && v[0].getRT() == 5.0 && v[1].getPrecursors()[0].getMZ() == 700.0 && v[2].getPrecursors()[0].getMZ() == 500.001
            && v[3].getRT() == 12.0,
          "every survivor keeps its own place (the later twin here), the folded member leaves no gap");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.0, 2, 1.001, 2, peaks_in)};
      chk(mergeTwins(v, 0.0) == 1 && v.size() == 2, "tolerance 0 folds identical-1/K0 twins only");
    }
    chk(pair_folds(sp(10.0, 500.00975, 2, 1.0, 2, peaks_in)), "19.5 ppm apart: folded (the 20-ppm relation, pinned from below)");
    chk(! pair_folds(sp(10.0, 500.01025, 2, 1.0, 2, peaks_in)), "20.5 ppm apart: untouched (pinned from above)");
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 2, peaks_in), sp(10.0, 500.001, 2, 1.25, 2, peaks_in)};
      chk(mergeTwins(v, 0.25) == 1, "|d1/K0| exactly the tolerance: folded (the bound is inclusive)");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.001, 2, 1.0, 2, peaks_in), sp(10.0, 500.0, 2, 1.0, 2, peaks_in)};
      chk(mergeTwins(v, 0.005) == 1 && v.size() == 1, "the higher m/z first in the input: still folded (the backward scan)");
    }
    {
      vector<MSSpectrum> v {sp(10.0, 500.0, 2, 1.0, 3, {{100.0, 1.0f}, {200.0, 1.0f}}),
                            sp(10.0, 500.001, 2, 1.0, 2, {{110.0, 1.0f}, {210.0, 1.0f}, {310.0, 1.0f}})};
      mergeTwins(v, 0.005);
      chk(v.size() == 1 && v[0].getPrecursors()[0].getMZ() == 500.0, "isotopes outrank peaks: 3 isotopes / 2 peaks beats 2 isotopes / 3 peaks");
    }
    {
      vector<MSSpectrum> v {sp(20.0, 400.0, 2, 0.9, 2, {{400.0, 1.0f}, {400.003, 2.0f}})};
      chk(mergeTwins(v, 0.005) == 0 && v.size() == 1 && v[0].size() == 2 && v[0][0].getMZ() == 400.0 && v[0][1].getMZ() == 400.003
            && v[0][0].getIntensity() == 1.0f && v[0][1].getIntensity() == 2.0f,
          "a spectrum without a twin keeps its peaks bit for bit, 7.5 ppm apart or not");
    }
    // the fusion anchor is the RUNNING intensity-weighted mean, at 10 ppm of it
    auto fused = [&](const vector<pair<double, float>>& a, const vector<pair<double, float>>& b) {
      vector<MSSpectrum> v {sp(30.0, 600.0, 2, 1.0, 3, a), sp(30.0, 600.0, 2, 1.0, 2, b)};
      mergeTwins(v, 0.005);
      return v.size() == 1 ? v[0].size() : (Size)0;
    };
    chk(fused({{200.0, 1.0f}}, {{200.0021, 1.0f}}) == 2, "fusion: 10.5 ppm apart stays two peaks");
    chk(fused({{200.0, 1.0f}, {200.0038, 1.0f}}, {{200.0019, 1.0f}}) == 2,
        "fusion: equal weights, the third peak is 14.25 ppm from the running mean: two peaks");
    chk(fused({{200.0, 1.0f}, {200.0038, 1.0f}}, {{200.0019, 100.0f}}) == 1,
        "fusion: a heavy middle peak pulls the mean to 9.6 ppm of the third: one peak (not anchored on the first)");
    {
      MSSpectrum a = sp(40.0, 700.0, 3, 1.1, 3, peaks_in), b = sp(40.0, 700.002, 3, 1.101, 2, peaks_in);
      OpenMS::Precursor pa = a.getPrecursors()[0];
      pa.setIsolationWindowLowerOffset(7.0);
      pa.setIsolationWindowUpperOffset(18.0);
      a.setPrecursors({pa});
      a.setMetaValue("spx_guessed", 0);
      b.setMetaValue("spx_guessed", 1);
      vector<MSSpectrum> v {a, b};
      mergeTwins(v, 0.005);
      const auto& q = v[0].getPrecursors()[0];
      chk(v.size() == 1 && v[0].getRT() == 40.0 && q.getMZ() == 700.0 && q.getCharge() == 3 && q.getDriftTime() == 1.1
            && q.getIsolationWindowLowerOffset() == 7.0 && q.getIsolationWindowUpperOffset() == 18.0 && (int)v[0].getMetaValue("spx_guessed") == 0,
          "the survivor keeps its own RT, precursor, 1/K0, isolation offsets and userParams");
    }
    writeLogInfo_(fail ? "[twins] SELFTEST FAILED" : "[twins] selftest passed");
    return fail ? UNEXPECTED_RESULT : EXECUTION_OK;
  }

  /**
    @brief Export-time calibration: the reported m/z is computed once, here, from the bin and its apex frame's
    factor.

    tof 0 (the OpenMS detector has no bin) keeps the cached m/z.

    @param[in] tof Flight-time bin
    @param[in] b Apex frame calibration factor
    @param[in] cached The cached m/z
    @return The m/z to write
  */
  static double exportMz_(uint32_t tof, double b, double cached)
  {
    if (tof == 0 || ! (b > 0.0) || ! tofAxis().ok) { return cached; }
    const double v = tofAxis().mzOf(tof, b);
    return std::isfinite(v) && v > 0.0 ? v : cached;
  }

  /**
    @brief [trim] malloc_trim(0), logged with the RSS before and after and how long it took.

    @param[in] what Label for the log line
  */
  void trimLogged_(const String& what) const
  {
    const double t0 = phaseClock();
    const long r0 = rssMb();
#ifdef __GLIBC__
    malloc_trim(0);
#endif
    writeLogInfo_("[trim] " + what + ": " + String(r0) + " -> " + String(rssMb()) + " MB RSS in " + numStr(phaseClock() - t0) + " s");
  }

  /**
    @brief The completed <out>.part becomes <out>: the data reach the disk before the name does, so a crash after the rename
    never leaves a short file under <out> (a file system without fsync, EINVAL/ENOTSUP, is not an error);
    std::filesystem::rename also replaces an existing <out> on Windows.

    @param[in] part The completed <out>.part
    @param[in] out The final path
    @throws Exception::UnableToCreateFile if the flush or the rename fails
  */
  void commitOutput_(const String& part, const String& out)
  {
#ifndef _WIN32
    const int fd = ::open(part.c_str(), O_RDONLY);
    const bool synced = fd >= 0 && (::fsync(fd) == 0 || errno == EINVAL || errno == ENOTSUP);
    if (fd >= 0) { ::close(fd); }
    if (! synced)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, part, "cannot flush the completed output to disk");
    }
#endif
    std::error_code ec;
    std::filesystem::rename(std::string(part), std::string(out), ec);
    if (ec)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, out,
                                          "cannot rename the completed " + part + ": " + ec.message());
    }
  }

  /**
    @brief Extract a FAIMS run per compensation voltage on the non-IMS path and write the merged result.

    [faims] A FAIMS run interleaves several compensation voltages; each is a run of its own for the non-IMS path, so each
    is extracted alone and its spectra carry their voltage. The merged result is in retention-time order.

    @param[in] exp The input (moved from)
    @param[in] out Output path
    @param[in] out_type Output file type
    @return EXECUTION_OK, or the error code of a part
  */
  ExitCodes runFaims_(PeakMap&& exp, const String& out, FileTypes::Type out_type)
  {
    PeakMap res;
    auto parts = IMDataConverter::splitByFAIMSCV(std::move(exp));
    writeLogInfo_("[faims] " + String(parts.size()) + " compensation voltages: each takes the non-IMS path on its own.");
    for (auto& [cv, sub] : parts)
    {
      writeLogInfo_("[faims] CV " + String(cv) + " V: " + String(sub.size()) + " spectra");
      PeakMap part;
      if (const ExitCodes rc = noImsSpectra_(sub, part); rc != EXECUTION_OK) { return rc; }
      for (MSSpectrum& sp : part)
      {
        sp.setDriftTime(cv);
        sp.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
        res.addSpectrum(std::move(sp));
      }
    }
    std::stable_sort(res.getSpectra().begin(), res.getSpectra().end(),
                     [](const MSSpectrum& a, const MSSpectrum& b) { return a.getRT() < b.getRT(); });
    for (Size i = 0; i < res.size(); ++i)
    {
      res[i].setNativeID("index=" + String(i));
    }
    return writeNoIms_(res, out, out_type);
  }

  /**
    @brief Run the non-IMS path on @p exp and write the result.

    [noims] input without ion mobility: the separate path of NoImsPath.h, written to <out>.part and committed like the
    IMS path's output. Centroided input only: profile spectra would make every profile sample a peak.

    @param[in,out] exp The input (its spectra are moved out)
    @param[in] out Output path
    @param[in] out_type Output file type
    @return EXECUTION_OK, or the error code
  */
  ExitCodes runNoIms_(PeakMap& exp, const String& out, FileTypes::Type out_type)
  {
    PeakMap res;
    if (const ExitCodes rc = noImsSpectra_(exp, res); rc != EXECUTION_OK) { return rc; }
    return writeNoIms_(res, out, out_type);
  }

  // TOF defaults of the non-IMS path (PXD002952 TripleTOF 6600 SWATH, Sage peptides at 1 %)
  static constexpr double TOF_MS2_PPM = 4.0, TOF_SPLIT = 0.5, TOF_PPM = 7.0, TOF_MIN_CORR = 0.5, TOF_WEIGHT_POW = 2.0, TOF_COMP_BOOST = 2.0;
  static constexpr int TOF_MIN_CHARGE = 2, TOF_MIN_HILL_SCANS = 3, TOF_SPAN_PAD = 5, TOF_HILL_SEEDS = 3,
                       TOF_PAIR_TOP = 200; // hill seeds: 3 pairs beat 1 on six TripleTOF runs with entrapment (+413 peptides, 45 % fewer spectra)

  /**
    @brief Validate the non-IMS input, apply the TOF defaults when the instrument is a TOF, and run NoIms::run.

    @param[in,out] exp The input (its spectra are moved out)
    @param[out] res The pseudo-MS/MS spectra
    @return EXECUTION_OK, or INCOMPATIBLE_INPUT_DATA for profile spectra or an MS2 spectrum without exactly one precursor
  */
  ExitCodes noImsSpectra_(PeakMap& exp, PeakMap& res)
  {
    for (const MSSpectrum& s : exp)
    {
      if (s.getType(true) == SpectrumSettings::SpectrumType::PROFILE)
      {
        writeLogError_("Error: the input has profile spectra (MS level " + String(s.getMSLevel()) + ", " + s.getNativeID()
                       + "); the non-IMS path needs centroided MS1 and MS2 (centroid them first, e.g. with PeakPickerHiRes).");
        return INCOMPATIBLE_INPUT_DATA;
      }
      if (s.getMSLevel() == 2 && s.getPrecursors().size() > 1)
      { // multiplexed DIA: one <precursor>, with its own isolation window, per window
        writeLogError_("Error: MS2 spectrum " + s.getNativeID() + " has " + String(s.getPrecursors().size())
                       + " isolation windows (multiplexed DIA); the non-IMS path reads one isolation window per MS2 scan.");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
    writeLogInfo_("[noims] the input has no ion mobility: using the separate non-IMS path.");
    // the path walks scans in file order and binary-searches their RTs: put an unordered file in RT order (stable, so
    // an ordered file is left exactly as it is)
    auto rt_less = [](const MSSpectrum& a, const MSSpectrum& b) { return a.getRT() < b.getRT(); };
    if (! std::is_sorted(exp.begin(), exp.end(), rt_less))
    {
      std::stable_sort(exp.getSpectra().begin(), exp.getSpectra().end(), rt_less);
      writeLogInfo_("[noims] the input spectra were not in retention-time order; sorted them.");
    }
    DIAspeXtract::NoIms::Params noims_opts;
    noims_opts.min_corr = getDoubleOption_("noims:min_corr");
    noims_opts.ppm = getDoubleOption_("noims:ppm");
    noims_opts.tof = getStringOption_("noims:tof") == "on";
    noims_opts.min_intensity = getDoubleOption_("noims:min_intensity");
    noims_opts.min_hill_scans = getIntOption_("noims:min_hill_scans");
    noims_opts.min_charge = getIntOption_("noims:min_charge");
    noims_opts.min_frags = getIntOption_("noims:min_fragments");
    noims_opts.min_frag_intensity = getDoubleOption_("noims:min_fragment_intensity");
    noims_opts.min_frag_points = getIntOption_("noims:min_fragment_points");
    noims_opts.max_frags = getIntOption_("noims:max_fragments");
    noims_opts.weight_pow = getDoubleOption_("noims:weight_power");
    noims_opts.ms2_ppm = getDoubleOption_("noims:ms2_ppm");
    noims_opts.split = getDoubleOption_("noims:split");
    noims_opts.env_all_charges = getStringOption_("noims:envelope_charges") == "all";
    noims_opts.drop_precursor = getStringOption_("noims:precursor_peaks") == "drop";
    noims_opts.comp_boost = getDoubleOption_("noims:complement_boost");
    noims_opts.span_pad = getIntOption_("noims:span_pad");
    noims_opts.env_pattern = getStringOption_("noims:envelope_pattern") == "on";
    noims_opts.hill_seeds = getIntOption_("noims:hill_seeds");
    noims_opts.tag_src = getStringOption_("noims:tag_source") == "on";
    noims_opts.pair_top = getIntOption_("noims:pair_top");
    noims_opts.seed_mono = getStringOption_("noims:seed_mono") == "on";
    {
      const String m = getStringOption_("noims:ms2_shift");
      noims_opts.ms2_shift_auto = m == "auto";
      noims_opts.ms2_shift = noims_opts.ms2_shift_auto ? 0 : std::atof(m.c_str());
    }
    // [tof] TOF analyzers (SWATH) record each fragment m/z with more scatter and many unrelated ions per scan: wide MS2
    // traces join neighbouring ions and recurring m/z into one long trace whose single apex reaches one precursor only.
    // TOF input therefore gets narrow MS2 traces split into elution peaks, and fewer precursors: charge >= 2, seen in
    // >= 3 MS1 scans (charge 1 gave 0.2 % PSMs per spectrum on SWATH). Only options left at their defaults change; Orbitrap and FTICR input keeps
    // every default.
    String inst = getStringOption_("noims:instrument");
    if (inst == "auto")
    {
      bool tof = false, ft = false;
      for (const auto& a : exp.getInstrument().getMassAnalyzers())
      {
        tof |= a.getType() == MassAnalyzer::AnalyzerType::TOF;
        ft |= a.getType() == MassAnalyzer::AnalyzerType::ORBITRAP || a.getType() == MassAnalyzer::AnalyzerType::FOURIERTRANSFORM;
      }
      inst = tof && ! ft ? "tof" : "orbitrap";
    }
    if (inst == "tof")
    {
      const DIAspeXtract::NoIms::Params noims_dflt;
      // ponytail: "explicit" = named on the command line; a value set only in an -ini file counts if it differs from the general default
      const auto dflt
        = [](const char* o, bool at_default) { return at_default && g_noims_cli.find(std::string("-noims:") + o + " ") == std::string::npos; };
      if (dflt("ms2_ppm", noims_opts.ms2_ppm == noims_dflt.ms2_ppm)) { noims_opts.ms2_ppm = TOF_MS2_PPM; }
      if (dflt("split", noims_opts.split == noims_dflt.split)) { noims_opts.split = TOF_SPLIT; }
      if (dflt("min_charge", noims_opts.min_charge == noims_dflt.min_charge)) { noims_opts.min_charge = TOF_MIN_CHARGE; }
      if (dflt("min_hill_scans", noims_opts.min_hill_scans == noims_dflt.min_hill_scans)) { noims_opts.min_hill_scans = TOF_MIN_HILL_SCANS; }
      // tuned on PXD002952 run 008 against Sage peptides at 1 %: tighter tolerances, flatter correlation weight, wider RT span,
      // complement boost, averagine-checked MS2 envelopes, hill seeds, and a per-run fragment m/z shift
      if (dflt("ppm", noims_opts.ppm == noims_dflt.ppm)) { noims_opts.ppm = TOF_PPM; }
      if (dflt("min_corr", noims_opts.min_corr == noims_dflt.min_corr)) { noims_opts.min_corr = TOF_MIN_CORR; }
      if (dflt("weight_power", noims_opts.weight_pow == noims_dflt.weight_pow)) { noims_opts.weight_pow = TOF_WEIGHT_POW; }
      if (dflt("complement_boost", noims_opts.comp_boost == noims_dflt.comp_boost)) { noims_opts.comp_boost = TOF_COMP_BOOST; }
      if (dflt("span_pad", noims_opts.span_pad == noims_dflt.span_pad)) { noims_opts.span_pad = TOF_SPAN_PAD; }
      if (dflt("hill_seeds", noims_opts.hill_seeds == noims_dflt.hill_seeds)) { noims_opts.hill_seeds = TOF_HILL_SEEDS; }
      if (dflt("pair_top", noims_opts.pair_top == noims_dflt.pair_top)) { noims_opts.pair_top = TOF_PAIR_TOP; }
      if (dflt("envelope_pattern", ! noims_opts.env_pattern)) { noims_opts.env_pattern = true; }
      if (dflt("seed_mono", ! noims_opts.seed_mono)) { noims_opts.seed_mono = true; }
      if (dflt("ms2_shift", ! noims_opts.ms2_shift_auto && noims_opts.ms2_shift == 0)) { noims_opts.ms2_shift_auto = true; }
      writeLogInfo_("[noims] TOF input: MS2 traces at " + String(noims_opts.ms2_ppm) + " ppm split at valleys below " + String(noims_opts.split)
                    + ", precursors of charge >= " + String(noims_opts.min_charge) + " seen in >= " + String(noims_opts.min_hill_scans)
                    + " MS1 scans");
    }
    ProgressLogger noims_progress;
    noims_progress.setLogType(log_type_);
    DIAspeXtract::NoIms::run(exp, noims_opts, res, [this](const String& m) { writeLogInfo_(m); }, noims_progress); // moves the spectra out of exp
    return EXECUTION_OK;
  }

  /**
    @brief Write the non-IMS result to <out>.part, check that an mzML is complete, and commit it.

    @param[in] res The spectra
    @param[in] out Output path
    @param[in] out_type Output file type
    @return EXECUTION_OK
    @throws Exception::UnableToCreateFile if the mzML is incomplete or the commit fails
  */
  ExitCodes writeNoIms_(PeakMap& res, const String& out, FileTypes::Type out_type)
  {
    const String part = out + ".part";
    FileHandler().storeExperiment(part, res, {out_type}, log_type_);
    if (out_type == FileTypes::MZML) // FileHandler does not report a failed write: the file must end in its footer
    {
      std::ifstream in(part.c_str(), std::ios::binary | std::ios::ate);
      const std::streamoff size = in ? (std::streamoff)in.tellg() : (std::streamoff)0;
      std::string tail(size > 64 ? 64 : (Size)std::max<std::streamoff>(size, 0), ' ');
      if (size > 0)
      {
        in.seekg(size - (std::streamoff)tail.size());
        in.read(&tail[0], (std::streamsize)tail.size());
      }
      if (! in || tail.find("</indexedmzML>") == std::string::npos)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, part, "writing the mzML failed (disk full or quota?)");
      }
    }
    commitOutput_(part, out);
    writeLogInfo_("Wrote " + String(res.size()) + " pseudo-MS2 spectra to " + out + " (" + String(FileTypes::typeToName(out_type)) + ")");
    return EXECUTION_OK;
  }

  /**
    @brief The tool's main: load the input, pick and trace MS1 and MS2, infer precursors, assemble and write the pseudo-MS/MS spectra.

    @param[in] argc Argument count
    @param[in] argv Arguments
    @return The tool's exit code
  */
  ExitCodes main_(int argc, const char** argv) override
  {
    phaseClock(); // [perf-instr] fix the epoch (a first-call static) before loading

    // -threads not given: use every core. The command line is scanned, so an explicit -threads 1 stays serial.
    int n_threads_req = getIntOption_("threads");
    {
      bool given = false;
      for (int i = 1; i < argc && ! given; ++i)
      {
        given = (String(argv[i]) == "-threads");
      }
      // the physical cores this process may use: affinity mask and cgroup quota (Linux), hw.physicalcpu (macOS),
      // the affinity mask (Windows); SMT siblings count once
      const int phys = (int)SysInfo::getAvailablePhysicalCores();
      if (! given && n_threads_req == 1)
      {
        n_threads_req = phys;
        omp_set_num_threads(n_threads_req);
        writeLogInfo_("-threads not given: using all " + String(n_threads_req)
                      + " cores. Pass -threads <n> to use fewer, "
                        "for example on a shared machine.");
      }
      // [threads] More threads than physical cores only add SMT siblings and oversubscription: with libomp (macOS,
      // Windows) an oversubscribed team stalls in lock hand-offs to descheduled threads (64 on 32 cores: 4 h instead of
      // 7 min). Clamped, with a warning; the output does not depend on the thread count.
      if (n_threads_req < 1)
      {
        n_threads_req = phys; // -threads 0 (the GUI's "all cores"): TOPPBase set the same count
      }
      if (n_threads_req > phys)
      {
        writeLogWarn_("-threads " + String(n_threads_req) + " is more than the " + String(phys) + " physical cores here: using " + String(phys)
                      + ".");
        n_threads_req = phys;
        omp_set_num_threads(n_threads_req);
      }
    }
    // [s3] an s3:// input is copied to a local file once, then read like any other: the readers (mzPeak's libzip archive,
    // the parallel mzML parser) all want a seekable local file. The copy costs the object's size in temporary disk.
    // ponytail: whole-object fetch; range reads straight from S3 would need the mzPeak library's zip I/O reworked.
    const String in_arg = getStringOption_("in");
    S3Copy s3_copy;
    const String in = S3InputSource::isS3URI(in_arg) ? s3_copy.fetch(in_arg, [this](const String& m) { writeLogInfo_(m); }) : in_arg;
    if (! S3InputSource::isS3URI(in_arg) && ! File::exists(in))
    {
      throw OpenMS::Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in);
    }
    {
      const FileTypes::Type t = FileHandler::getTypeByFileName(in);
      if (t != FileTypes::MZML && t != FileTypes::MZPEAK && t != FileTypes::BRUKER_TDF && t != FileTypes::RAW)
      {
        throw OpenMS::Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                  "-in must be an mzML, mzPeak or Thermo .raw file or a Bruker .d folder: " + in_arg);
      }
    }
    const String out = getStringOption_("out");
    const double delta_im = getDoubleOption_("gate:delta_im");
    // [frag-im] the gate's fragment-precursor 1/K0 window, apart from the tracing tolerance delta_im (negative = delta_im)
    const double gate_dim = getDoubleOption_("gate:fragment_delta_im") >= 0.0 ? getDoubleOption_("gate:fragment_delta_im") : delta_im;
    const double delta_rt = getDoubleOption_("gate:delta_rt");
    if (getFlag_("diag:selftest_arena")) { return selftestArena_(); }
    if (getFlag_("diag:selftest_twins")) { return selftestTwins_(); }
    cfg_im_weight_sigma = getDoubleOption_("assembly:im_weight_sigma");
    {
      const String m = getStringOption_("charge:mono_profile_check");
      cfg_mono_check = m == "v3" ? 3 : m == "v2" ? 2 : 0;
    }
    cfg_mono_pos_min_ratio = getDoubleOption_("charge:mono_position_min_ratio");
    cfg_mono_pos_im_box = getDoubleOption_("charge:mono_position_im_box");
    cfg_mono_pos_extend = getDoubleOption_("charge:mono_position_extend");
    cfg_mono_pos_iterate = getStringOption_("charge:mono_position_iterate") == "true";
    cfg_short_conc = (getStringOption_("gate:short_traces") == "concentric");
    cfg_short_center_tol = getDoubleOption_("gate:short_center_tol");
    cfg_short_min_pts = getIntOption_("gate:short_min_points");
    cfg_coloc = (getStringOption_("gate:colocation") == "ranges");
    cfg_coloc_overlap = getDoubleOption_("gate:colocation_overlap");
    cfg_shrink = (getStringOption_("gate:shrinkage") == "fisher");
    cfg_shrink_r0 = getDoubleOption_("gate:shrink_prior_r");
    cfg_shrink_sd = getDoubleOption_("gate:shrink_prior_sd");
    cfg_shrink_min_pts = getIntOption_("gate:shrink_min_points");
    cfg_adaptive = (getStringOption_("gate:adaptive") == "faint");
    cfg_adaptive_auto = (getStringOption_("gate:adaptive") == "auto");
    cfg_auto_faint_frac = getDoubleOption_("gate:auto_faint_fraction");
    cfg_im_slice_merge = getDoubleOption_("trace:im_slice_merge");
    cfg_merge_same_bin = getStringOption_("assembly:merge_same_bin") == "true";
    cfg_foreign_w = getDoubleOption_("assembly:foreign_precursor_weight");
    cfg_foreign_ppm = getDoubleOption_("assembly:foreign_precursor_ppm");
    cfg_foreign_dim = getDoubleOption_("assembly:foreign_precursor_dim");
    cfg_foreign_min_iso = getIntOption_("assembly:foreign_precursor_min_iso");
    cfg_fuse_ppm = getDoubleOption_("assembly:fragment_fuse_ppm");
    cfg_faint_frames = getIntOption_("gate:faint_max_frames");
    cfg_spectrum_meta = getFlag_("diag:spectrum_meta");
    cfg_faint_drt = getDoubleOption_("gate:faint_delta_rt");
    cfg_faint_cp = getDoubleOption_("gate:faint_corr_power");
    cfg_corr_power = getDoubleOption_("assembly:corr_power");
    const double min_corr = getDoubleOption_("gate:min_correlation");
    const int min_corr_pts = getIntOption_("gate:min_correlation_points");
    const Size min_frags = (Size)getIntOption_("assembly:min_fragments");
    const Size max_frags = (Size)getIntOption_("assembly:max_fragments");
    if (max_frags < min_frags)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "assembly:max_fragments (" + String(max_frags) + ") is below assembly:min_fragments (" + String(min_frags)
                                          + "): every written spectrum would hold fewer fragments than the minimum.");
    }
    const double twin_im = getDoubleOption_("assembly:twin_im_tolerance");
    const int twin_rt_frames = getIntOption_("assembly:twin_rt_frames"); // [twin-rt] -1 off, 0 auto, n frames
    const int max_charge = getIntOption_("max_charge");
    if (getIntOption_("charge:min_charge") > max_charge)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "charge:min_charge (" + String(getIntOption_("charge:min_charge")) + ") is above max_charge ("
                                          + String(max_charge) + "): no precursor charge could be written.");
    }
    const double mass_ppm = getDoubleOption_("trace:mass_error_ppm");

    PeakMap exp;
    // Accept mzML, mzPeak, or a Bruker .d directory; FileHandler auto-detects the format
    // (mzPeak -> MzPeakFile; .d -> BrukerTimsFile, requires WITH_OPENTIMS).
    PeakMap ms1_map;
    // [ms1-prune] one run-level prune (see Ms1Prune), used by the streaming consumer and the resident dispatch alike
    Ms1Prune ms1_prune;
    ms1_prune.on = (getStringOption_("perf:ms1_prune") == "true");
    ms1_prune.noise = getDoubleOption_("trace:noise_threshold_int");
    ms1_prune.hop = mass_ppm * 1e-6 * BAND_HALO_PPM_MUL;
    {
      const char* e = std::getenv("DIASPEXTRACT_MS1_PRUNE_NO_WITNESS");
      ms1_prune.no_witness = e && *e && *e != '0';
    }
    {
      const char* e = std::getenv("DIASPEXTRACT_MS1_PRUNE_NO_CHAIN");
      ms1_prune.no_chain = e && *e && *e != '0';
    } // its own switch: an e2e falsifier
    map<WinKey, PeakSlab> ms2_by_window;
    CompactStats cstat;
    // an mzML is always read whole: the frame-by-frame reader reads a .d (or an mzPeak archive) only
    const FileTypes::Type in_type = FileHandler::getTypeByFileName(in);
#ifdef DIASPEXTRACT_WITH_MZPEAK
    // [noims] an mzPeak archive without ion mobility is read whole, as an mzML is
    const bool streamable = in_type == FileTypes::BRUKER_TDF || (in_type == FileTypes::MZPEAK && DIAspeXtract::mzPeakHasIonMobility(in));
#else
    const bool streamable = in_type == FileTypes::BRUKER_TDF;
#endif
    const bool stream_load = getStringOption_("perf:stream_load") == "true" && streamable;
    const bool trim_on = (getStringOption_("perf:malloc_trim") == "true");
    bool streamed = false;
    // Output format, resolved before the load: it decides the tiling below (mzPeak output runs as one tile), and a
    // plain build refuses -out_type mzpeak here in milliseconds. The extension is authoritative, because that is what
    // a user typing `-out pseudo.mzML` means; -out_type overrides it; mzPeak is the default when neither says.
    FileTypes::Type out_type = FileTypes::MZML;
    const String out_type_opt = getStringOption_("out_type");
    if (! out_type_opt.empty()) { out_type = (out_type_opt == "mzML") ? FileTypes::MZML : FileTypes::MZPEAK; }
    else
    {
      const FileTypes::Type by_name = FileHandler::getTypeByFileName(out);
      out_type = (by_name == FileTypes::MZML || by_name == FileTypes::MZPEAK) ? by_name : FileTypes::MZPEAK;
    }
    if (out_type == FileTypes::MZPEAK)
    {
      writeLogInfo_("Writing mzPeak. Search engines read mzML: for a search, name the output .mzML (or pass -out_type mzML).");
    }
    // [preflight] what the .d loader would abort on, refused here with its cause and an exit code a caller can route on:
    // a ddaPASEF run (INCOMPATIBLE_INPUT_DATA) and a cut analysis.tdf_bin (INPUT_FILE_CORRUPT). First, before the
    // calibration and dnoise setups read the same tdf, so a broken file is named as such whatever it breaks there.
    if (FileHandler::getTypeByFileName(in) == FileTypes::BRUKER_TDF)
    {
      bool dda = false;
      if (const std::string w = DIAspeXtract::tdfNotDiaReason(in + "/analysis.tdf", dda); ! w.empty())
      {
        writeLogError_("Error: " + in + " " + w + ". DIAspeXtract needs diaPASEF data" + (dda ? "; search ddaPASEF runs directly." : "."));
        return dda ? INCOMPATIBLE_INPUT_DATA : INPUT_FILE_CORRUPT;
      }
      if (const std::string w = DIAspeXtract::tdfBinTruncationReason(in); ! w.empty())
      {
        writeLogError_("Error: " + in + ": " + w + ".");
        return INPUT_FILE_CORRUPT;
      }
    }
    // [tof] The integer detector's band edges come from the acquisition: they need the tdf's GlobalMetadata bounds, so
    // refuse BEFORE the load when no candidate tdf provides them (an mzPeak archive is checked after its embedded
    // tdf is recovered, below). Milliseconds against minutes of loading.
    const bool integer_req = (getStringOption_("trace:detector") == "integer");
    std::vector<std::string> tdf_cand {in + "/analysis.tdf"}; // the calibration sources, in the order the load tries them
    if (const char* sc = std::getenv("DIASPEXTRACT_MZPEAK_TDF")) { tdf_cand.push_back(sc); }
    String why = "no candidate tdf"; // why no candidate calibrates the run, for the fallback warning below
    if (integer_req)
    {
      // [tof] The axis is selected HERE, in the order the load below would use it (the first
      // candidate whose calibration loads is THE source), so the bounds that are checked belong to
      // the tdf that will actually calibrate the run -- and the check costs a two-row query rather
      // than the whole load. The candidates are the input's own analysis.tdf (a .d has one) and
      // DIASPEXTRACT_MZPEAK_TDF beside an mzPeak archive. No candidate at all is NOT an error: the
      // integer detector then falls back to the OpenMS detector further down, loudly.
      for (const std::string& c : tdf_cand)
      {
        if (! loadTofAxis(c, why))
        {
          continue; // not a calibration source; try the next
        }
        if (! (tofAxis().n_bins > 1))
        {
          throw Exception::InvalidParameter(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "trace:detector integer needs the acquisition m/z range and digitizer size (GlobalMetadata) from the calibration source " + c
              + ", which lacks them (" + tofAxis().bounds_why + "). Use -trace:detector openms.");
        }
        writeLogInfo_("flight-time axis from " + c + " (preflight): " + String(tofAxis().b_by_frame.size()) + " frame factors; digitizer bins "
                      + String(tofAxis().n_bins) + ", acquisition m/z " + String(tofAxis().acq_lo) + "-" + String(tofAxis().acq_hi));
        break;
      }
    }
    // the .d loader's configuration; set here because the dnoise setup below predicts the loader's m/z model from it
    BrukerTimsFile::Config cfg;
    cfg.ms1_n_neighbors = getIntOption_("trace:native_ms1_neighbors");
    // [pipe] the reader's own decode team: on the first pass nothing else runs, so all threads (perf:decode_threads 0);
    // the pipelined tile reads below hand it the window team instead (dia_ms2_parallel_for)
    cfg.dia_ms2_decode_threads = getIntOption_("perf:decode_threads") > 0 ? getIntOption_("perf:decode_threads") : n_threads_req;
    cfg.ms1_decode_threads = cfg.dia_ms2_decode_threads; // the first pass reads the MS1 frames; the reader's own default is 4
    // [dnoise] dnoise v0.1.0's MS1 path on the raw frames (DnoiseRun), set up once for the run: every refusal the options
    // and the tdf can decide comes here, in milliseconds, before the load
    DnoiseRun dnoise;
    const bool dnoise_req = (getStringOption_("dnoise:ms1") == "true");
    if (dnoise_req && FileHandler::getTypeByFileName(in) != FileTypes::BRUKER_TDF)
    {
      writeLogInfo_("[dnoise] MS1 denoising skipped: the input is not a Bruker .d");
    }
    else if (dnoise_req)
    {
      auto refuse = [](const String& why) {
        return Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "dnoise:ms1 " + why + " Pass -dnoise:ms1 false to run without MS1 denoising.");
      };
      if (cfg.ms1_n_neighbors > 0) { throw refuse("cannot be combined with trace:native_ms1_neighbors > 0."); }
      // The loader's m/z model, predicted: under the default AUTO strategy an SDK path (the Config's, else
      // OPENMS_BRUKER_SDK_PATH) makes openTimsDataHandle calibrate m/z through the Bruker SDK, which the TOF recovery
      // cannot invert. An SDK that then fails to load would leave the table model, so this can refuse a run that would
      // have recovered: fail-closed. checkLoader() confirms the loader's actual choice during the load.
      {
        const char* env_sdk = std::getenv("OPENMS_BRUKER_SDK_PATH");
        const std::string sdk = ! cfg.bruker_sdk_path.empty() ? cfg.bruker_sdk_path : std::string(env_sdk ? env_sdk : "");
        if (! sdk.empty())
        {
          throw refuse("needs the run's own m/z calibration table, but a Bruker SDK path is set (" + String(sdk)
                       + "): unset OPENMS_BRUKER_SDK_PATH to denoise.");
        }
      }
      String why_axis;
      if (! tofAxis().ok && ! loadTofAxis(in + "/analysis.tdf", why_axis)) // the TOF recovery runs during the load
      {
        if (DIAspeXtract::tdfReasonIsCorrupt(why_axis)) // values no instrument writes: the file is broken, not unsupported
        {
          writeLogError_("Error: " + in + ": the m/z calibration in analysis.tdf is corrupt (" + why_axis + ").");
          return INPUT_FILE_CORRUPT;
        }
        throw refuse("needs the run's m/z calibration table (MzCalibration in analysis.tdf), which is not usable here (" + why_axis + ").");
      }
      const std::string why = dnoise.setup(in + "/analysis.tdf");
      if (! why.empty() && DIAspeXtract::tdfReasonIsCorrupt(why))
      {
        writeLogError_("Error: " + in + ": analysis.tdf is corrupt (" + why + ").");
        return INPUT_FILE_CORRUPT;
      }
      if (! why.empty()) { throw refuse("cannot be exact on this input: " + why + "."); }
      writeLogInfo_("[dnoise] " + String(DIAspeXtract::Dnoise::PORTED_FROM) + " MS1 path on the raw frames, before the pick: " + dnoise.params()
                    + "; " + String(dnoise.im_asc.size()) + " TimsCalibration row(s) over " + String(dnoise.tdf.meta.max_num_scans)
                    + " scans (MAX(Frames.NumScans), the gate's); scan tables per row over its MS1 frames' scans: " + dnoise.scanTables());
    }
    String frame_table_src = "slab"; // [frames] where the frozen frame tables came from: loader | mzpeak | slab
    // [stream] the .d streaming SOURCE: pass 1 reads MS1 only (the loader exports the whole frame
    // table regardless), pass 2 reads each tile's MS2 frames through the loader's RT range. The
    // resident source (mzPeak, mzML, perf:stream_load=false) holds the whole run and slices it.
    bool tile_source = false;
    PickCompactConsumer::FrameTables run_tables; // the reader's frozen tables, kept for pass 2
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
      // [stream] the streaming source needs the Bruker loader: a .d only, never an mzML or an mzPeak
      // archive (the resident path reads those whole).
      tile_source = ! in_mzpeak && FileHandler::getTypeByFileName(in) == FileTypes::BRUKER_TDF;
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
      String table_src_if_ok = "slab"; // becomes the stamp only after the streamed load and the invariant succeed
      if (in_mzpeak)
      {
#ifdef DIASPEXTRACT_WITH_MZPEAK
        table_src_if_ok = "mzpeak";
        if (integer_req)
        {
          // [tof] the sweep recovers the embedded tdf's calibration and bounds; it is cached, so
          // the load below reuses it -- the refusal comes before any peak is decoded
          const auto m = DIAspeXtract::mzPeakRunMeta(in);
          if (! m->exact || ! (m->exact->n_bins > 1))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "trace:detector integer needs the acquisition m/z range and digitizer size from the analysis.tdf in or "
                                              "beside the mzPeak archive, which lacks them ("
                                                + String(m->exact ? m->exact->bounds_why : "no exact calibration")
                                                + "). Use -trace:detector openms.");
          }
        }
        consumer.tablefn = [&]() {
          PickCompactConsumer::FrameTables t;
          const auto m = DIAspeXtract::lastMzPeakMeta();
          if (! m) { return t; }
          for (Size i : m->ms2)
          {
            if (m->n_peaks[i] == 0) { continue; }
            for (const auto& w : m->wins[i])
            {
              const double c = 0.5 * (w.lo + w.hi), off = 0.5 * (w.hi - w.lo); // frameToSpectra's precursor
              auto& tab = t[winKey(c - off, c + off, 0)];
              tab.rt.push_back(m->rt[i]);
              tab.id.push_back(m->frame_id[i] < 0 ? 0u : (uint32_t)m->frame_id[i]);
            }
          }
          return t;
        };
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
            if (g == ft.end()) { continue; }
            const double c = w.mz_center, off = w.mz_width / 2.0; // BrukerTimsFile's precursor
            const WinKey k = winKey(c - off, c + off, (uint32_t)w.window_group);
            if (t.count(k)) // two table rows with one (m/z, group) key -- the consumer would merge their frames into one slab
            {
              throw InputRefusal(
                __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "[frames] two DIA windows share one (m/z lo, m/z hi, group) key; DIAspeXtract keys windows by isolation window and group",
                String(c - off) + "-" + String(c + off) + " group " + String(w.window_group));
            }
            auto& tab = t[k];
            for (const auto& r : g->second)
            {
              if (r.num_peaks > 0)
              {
                tab.rt.push_back(r.time);
                tab.id.push_back(r.frame_id);
              }
            }
          }
          return t;
        };
      }
      try
      {
#ifdef DIASPEXTRACT_WITH_MZPEAK
        if (in_mzpeak)
        {
          DIAspeXtract::loadMzPeakStreaming(in, consumer, n_threads_req); // consumer counts frames itself
        }
        else
#endif
        {
          BrukerTimsFile().loadDIAStreaming(in, consumer, cfg);
        }
        consumer.finish(); // the last partial batch
        streamed = true;
        phaseAdd("LOAD(stream)", 0.0, 0.0, false); // [perf-instr] start-to-here
        writeLogInfo_("[stream] frame-by-frame load done: " + String(consumer.frames_seen) + " frames, MS1 " + String(ms1_map.size()) + ", windows "
                      + String(ms2_by_window.size()) + clockStamp() + rss());
      }
      catch (const InputRefusal&)
      {
        throw;
      } // a refusal of the input, not a reader failure: never fall back onto partial slabs
#ifdef DIASPEXTRACT_WITH_MZPEAK
      catch (const DIAspeXtract::MzPeakError&)
      {
        throw;
      } // an unusable archive: a whole-run read of it fails the same way
#endif
      catch (const Exception::BaseException& e)
      {
        writeLogWarn_(String("[stream] streaming load failed (") + e.getName() + "); reading the whole run instead. " + e.getMessage());
      }
      // [stream] the source is COMMITTED only by a successful pass 1: a fallback to the resident load
      // reads the whole run, so the tile loop must slice it rather than call the loader again
      if (! streamed && tile_source)
      {
        tile_source = false;
        writeLogWarn_("[stream] the streaming read did not finish; reading the whole run instead.");
      }
      // [frames] the invariant the frozen tables rest on, checked OUTSIDE the recoverable catch: a
      // delivered window's slab is its table (rt and id, frame for frame) after padding
      if (streamed && tile_source && consumer.tables().empty())
      {
        throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                           "[3b] the .d loader exported no frame table; the tiled source cannot read this input", in);
      }
      if (streamed && ! consumer.tables().empty())
      {
        Size n_pad = 0, n_del = 0, n_absent = 0;
        for (const auto& kv : consumer.tables())
        {
          const String wl = winLabel(kv.first);
          const auto w = ms2_by_window.find(kv.first);
          if (w == ms2_by_window.end())
          {
            ++n_absent;
            if (! tile_source)
            {
              writeLogInfo_("[frames] window " + wl + ": absent (no frame delivered), " + String(kv.second.rt.size()) + " in the table");
            }
            continue;
          }
          if (w->second.rt != kv.second.rt || w->second.frame_id != kv.second.id)
          {
            throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[frames] a window's slab is not its frozen frame table after padding",
                               wl);
          }
          n_pad += kv.second.padded;
          n_del += kv.second.delivered;
          writeLogInfo_("[frames] window " + wl + ": delivered " + String(kv.second.delivered) + ", padded " + String(kv.second.padded) + " of "
                        + String(kv.second.rt.size()));
        }
        frame_table_src = table_src_if_ok;
        run_tables = consumer.tables();
        writeLogInfo_("[frames] source=" + frame_table_src + " windows=" + String(consumer.tables().size()) + " delivered=" + String(n_del)
                      + " padded=" + String(n_pad) + " absent=" + String(n_absent)
                      + (tile_source ? " (pass 1: MS1 only; the windows' frames are read per tile)" : ""));
#ifdef DIASPEXTRACT_WITH_MZPEAK
        if (in_mzpeak && DIAspeXtract::lastMzPeakMeta())
        {
          writeLogInfo_("[frames] mzPeak number_of_peaks: " + String(DIAspeXtract::lastMzPeakMeta()->n_peaks_missing_ms2)
                        + " MS2 spectra lack it and count as present"
                        + (DIAspeXtract::lastMzPeakMeta()->n_peaks_missing_ms2 ? " (a frame with no peaks among them would be padded)"
                                                                               : "; a recorded zero is skipped like the .d loop's"));
        }
#endif
      }
    }

    if (! streamed)
    {
      // A failed stream attempt may already have delivered part of the run -- MS1 spectra (sampled, pruned and
      // denoised), window slabs and compaction counts. The resident load below reads the WHOLE
      // run again and only appends, so discard all of it first: kept, the delivered frames were loaded twice (silently
      // wrong with the prune off, a refusal with it on, "not in retention-time order" once partial slabs existed).
      ms1_map = PeakMap();
      ms2_by_window.clear();
      cstat = CompactStats();
      ms1_prune.reset();
      dnoise.reset(); // [dnoise] counts and frame marks: the resident pick below filters every frame again
      Phase _ph("LOAD(full)");
      FileHandler().loadExperiment(in, exp, {FileTypes::MZML, FileTypes::MZPEAK, FileTypes::BRUKER_TDF, FileTypes::RAW}, log_type_);
    }

    // Input normalization / validation, resident load only (a streamed run leaves `exp` empty; ms1_map.empty() guards it
    // below). Assumes 1/K0 (VSSC). Centroid IM frames that are still profile, keeping per-peak IM at a tolerance no coarser
    // than the trace IM tolerance, and split MS1 vs MS2 into separate working maps.
    if (! streamed)
    {
      if (! FAIMSHelper::getCompensationVoltages(exp).empty())
      {
        return runFaims_(std::move(exp), out, out_type); // [faims]
      }
      IMFormat imf = IMTypes::determineIMFormat(exp, 1);
      if (imf == IMFormat::NONE)
      {
        if (IMTypes::determineIMFormat(exp, 2) != IMFormat::NONE)
        {
          writeLogError_("Error: the input has ion mobility at MS2 but not at MS1; neither path reads such a mix.");
          return INCOMPATIBLE_INPUT_DATA;
        }
        return runNoIms_(exp, out, out_type); // [noims] a separate path; the IMS path below is never entered
      }
      // Pick all IM frames in parallel. firstprivate gives each thread a deep copy (Param is a value type), so the
      // mutable ccs_warning_shown_ flag is not shared.
      writeLogInfo_("Loaded " + String(exp.size()) + " frames." + clockStamp() + rss()); // raw residency peak? [mem]
      writeLogInfo_("Peak-picking " + String(exp.size()) + " frames (OpenMP)...");
      // An exception escaping the OMP region would terminate the process: captured, rethrown serially.
      std::exception_ptr pick_err;
      // [dnoise] as in pickMS1Batch_: the raw MS1 frames of a .d, before their pick (a frame it empties keeps its IM array,
      // so the dispatch below keeps it too)
      if (dnoise.on) { dnoise.checkLoader(); }
      std::vector<DnoiseRun::Work> dn_work(dnoise.on ? (Size)std::max(1, omp_get_max_threads()) : 0);
      std::vector<DnoiseRun::Counts> dn_cnt(dnoise.on ? exp.size() : 0);
#pragma omp parallel for firstprivate(spicker) schedule(dynamic, 16)
      for (long i = 0; i < (long)exp.size(); ++i)
      {
        try
        {
          MSSpectrum& s = exp[i];
          // [dnoise] before the ion-mobility skip, as the streaming consumer filters every MS1 frame: frame() refuses a
          // non-empty MS1 frame without 1/K0 and counts an empty one, so no MS1 frame bypasses it and both paths count alike
          if (dnoise.on && s.getMSLevel() == 1) { dnoise.frame(s, dn_work.at((Size)omp_get_thread_num()), dn_cnt[(Size)i]); }
          if (! s.containsIMData()) { continue; }
          if (s.getIMPeakType() != IMPeakType::IM_CENTROIDED) { spicker.pickIMCluster(s); }
          ensureIMArrayName(s); // so MassTraceDetection reads per-peak IM (IM-aware tracing)
          if (s.getMSLevel() == 1)
          {
            requireFiniteMs1Mz(s); // [ms1-finite] exactly the spectra the prune and the MS1 map take
          }
        }
        catch (...)
        {
          taskErr(pick_err);
        }
      }
      if (pick_err) { std::rethrow_exception(pick_err); }
      dnoise.add(dn_cnt);
      // [dnoise] each thread's scratch keeps the capacity of its largest frame (up to ~4 GiB at 100 threads): free it with
      // the pick, not at the end of this block after the prune and the dispatch
      std::vector<DnoiseRun::Work>().swap(dn_work);
      std::vector<DnoiseRun::Counts>().swap(dn_cnt);
      // [ms1-prune] exactly the MS1 spectra the dispatch below keeps, pruned as pickMS1Batch_ prunes a batch
      if (ms1_prune.on)
      {
        ms1_prune.batch(exp.getSpectra(), [](const MSSpectrum& s) { return s.containsIMData() && s.getMSLevel() == 1; });
      }

      // Serial dispatch: MOVE (not copy) frames into the MS1 map / per-window MS2 maps to avoid
      // transient memory doubling; container inserts aren't thread-safe so this stays serial.
      vector<uint32_t> cmzq;
      vector<float> cint;
      vector<uint16_t> cimq; // one MS2 frame's compaction region, reused
      for (MSSpectrum& s : exp)
      {
        if (! s.containsIMData()) { continue; }
        if (s.getMSLevel() == 1) { ms1_map.addSpectrum(std::move(s)); }
        else if (s.getMSLevel() == 2 && ! s.getPrecursors().empty())
        {
          const OpenMS::Precursor& pr = s.getPrecursors()[0];
          double lo = pr.getMZ() - pr.getIsolationWindowLowerOffset();
          double hi = pr.getMZ() + pr.getIsolationWindowUpperOffset();
          const uint32_t fid = frameIdOf(s.getNativeID());
          if (! frameCalibrated(fid)) { ++cstat.unmapped_frame; }
          if (cmzq.size() < s.size())
          {
            cmzq.resize(s.size());
            cint.resize(s.size());
            cimq.resize(s.size());
          }
          const Size kept = compactifyInto(s, cstat, cmzq.data(), cint.data(), cimq.data(), cmzq.size());
          ms2_by_window[winKey(lo, hi, windowGroupOf(s.getNativeID()))].append(s.getRT(), fid, cmzq.data(), cint.data(), cimq.data(), kept,
                                                                               ms1_map.size());
          s.clear(true); // [mem] the picked frame is dead once compacted; waiting for exp.clear()
        } // below held the whole picked run (20 B/peak) beside the compact store
      }
    }
    exp.clear(true); // frames moved out; release the container
    {
      Size cb = 0, cp = 0;
      for (const auto& kv : ms2_by_window)
      {
        cb += kv.second.bytes();
        cp += kv.second.peaks();
      }
      const Size lost = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
      if (lost)
      {
        writeLogWarn_("[compact] DROPPED " + String(lost) + " peaks: no_im_array=" + String(cstat.no_im_array)
                      + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz) + " bad_im=" + String(cstat.bad_im)
                      + " (kept " + String(cstat.kept) + "); they are left out of the analysis.");
      }
      Size ms1_pk = 0;
      for (const auto& sp : ms1_map)
      {
        ms1_pk += sp.size();
      }
      writeLogInfo_("Split into MS1 + " + String(ms2_by_window.size()) + " windows (picked); compact store " + String(cb / (1024ULL * 1024ULL))
                    + " MB for " + String(cp) + " peaks (~" + String(cp ? (double)cb / cp : 0.0) + " B/peak); MS1 " + String(ms1_map.size())
                    + " frames / " + String(ms1_pk) + " peaks (" + String(ms1_pk * 20 / (1024ULL * 1024ULL)) + " MB as PeakMap)." + clockStamp()
                    + rss()); // [mem]
      if (dnoise.on)          // [dnoise] after the MS1 load, whichever path picked it: every MS1 frame with points exactly once
      {
        dnoise.requireComplete();
        writeLogInfo_(dnoise.summary());
      }
      if (ms1_prune.on)
      {
        const Ms1Prune::Counts& t = ms1_prune.total;
        writeLogInfo_("[ms1-prune] picked " + String(t.picked) + " survivors " + String(t.survivors) + " witnesses " + String(t.kept - t.survivors)
                      + " at_noise " + String(t.at_noise) + " nan " + String(t.nan) + ": kept " + String(t.kept)
                      + " peaks = " + String(t.kept * 20 / (1024LL * 1024LL)) + " MB as PeakMap of " + String(t.picked * 20 / (1024LL * 1024LL))
                      + " MB picked; sample " + String(ms1_prune.sample.size()) + " m/z from " + String(ms1_prune.frames) + " spectra"
                      + (ms1_prune.no_witness
                           ? String(" -- DIASPEXTRACT_MS1_PRUNE_NO_WITNESS: survivors only, a negative control, output NOT identical")
                           : String())
                      + (ms1_prune.no_chain
                           ? String(" -- DIASPEXTRACT_MS1_PRUNE_NO_CHAIN: interior chain links dropped, a negative control, output NOT identical")
                           : String()));
      }
    }

    if (ms1_map.empty())
    {
      writeLogError_("Error: the input has no MS1 frames with ion mobility, so there are no precursors.");
      return INCOMPATIBLE_INPUT_DATA;
    }

    // [cell] the frozen frame tables (see FrozenWin), copied before indexRt() frees the slabs' frame
    // times and toTof() their frame ids
    map<WinKey, FrozenWin> frozen;
    if (tile_source)
    {
      for (const auto& kv : run_tables) // [stream] every window of the reader's table, no slab needed
      {
        FrozenWin& fw = frozen[kv.first];
        fw.frame_rt = kv.second.rt;
        fw.frame_id = kv.second.id;
      }
    }
    else
    {
      for (const auto& kv : ms2_by_window)
      {
        FrozenWin& fw = frozen[kv.first];
        fw.frame_rt = kv.second.rt;
        fw.frame_id = kv.second.frame_id;
      }
    }
    // The global RT axis (see rtAxis()): one entry per distinct frame time, MS1 and MS2 alike.
    {
      vector<double>& ax = rtAxis();
      ax.clear();
      {
        Size nf = ms1_map.size();
        for (const auto& kv : ms2_by_window)
        {
          nf += kv.second.frames();
        }
        ax.reserve(nf + 1024);
      } // MS1 *and* every window's MS2 frames land here [mem]
      for (const auto& sp : ms1_map)
      {
        ax.push_back(sp.getRT());
      }
      for (const auto& kv : frozen)
      {
        ax.insert(ax.end(), kv.second.frame_rt.begin(), kv.second.frame_rt.end()); // the tables' frame times (== the slabs' after padding)
      }
      sort(ax.begin(), ax.end());
      ax.erase(unique(ax.begin(), ax.end()), ax.end());
      writeLogInfo_("RT axis: " + String(ax.size()) + " distinct frame retention times"
                    + (ax.empty() ? "" : " spanning " + numStr(ax.front()) + "-" + numStr(ax.back()) + " s"));
      if (ax.empty())
      {
        writeLogError_("Error: the input's frames carry no retention times.");
        return INCOMPATIBLE_INPUT_DATA;
      }
      // [tile] frame times -> axis indices, and the load-time growth slack returned, window by
      // window (one window's copy in flight at a time, not all of them at once).
      for (auto& kv : ms2_by_window)
      {
        indexRt(kv.second);
      }
      for (auto& kv : frozen)
      {
        FrozenWin& fw = kv.second;
        fw.rt_index.resize(fw.frames());
        for (Size f = 0; f < fw.frames(); ++f)
        {
          fw.rt_index[f] = rtIndex(fw.frame_rt[f]);
        }
      }
#ifdef __GLIBC__
      if (trim_on)
      {
        malloc_trim(0); // the load's arena churn, returned before the MS1 phase
      }
#endif
    }

    //-------------------------------------------------------------
    // MS1 traces + precursor/charge inference (shared, once)
    //-------------------------------------------------------------
    // [outer] sortSpectra() = the RT sort, then each spectrum's own m/z sort (a no-op check on a picked spectrum). The
    // per-spectrum sorts are independent, so they run in parallel after the same RT sort.
    {
      ms1_map.sortSpectra(false);
      std::exception_ptr ss_err;
#pragma omp parallel for schedule(dynamic, 16)
      for (long i = 0; i < (long)ms1_map.size(); ++i)
      {
        try
        {
          ms1_map[(Size)i].sortByPosition();
        }
        catch (...)
        {
          taskErr(ss_err);
        }
      }
      if (ss_err) { std::rethrow_exception(ss_err); }
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
    vector<double> ms1_rt;
    ms1_rt.reserve(ms1_map.size());
    for (const auto& sp : ms1_map)
    {
      ms1_rt.push_back(sp.getRT());
    }
    std::sort(ms1_rt.begin(), ms1_rt.end());
    const double cycle = frameSpacing(ms1_rt);
    if (cell_pitch > 0.0)
    {
      const Size n = cycle > 0.0 ? std::max<Size>(1, (Size)llround(cell_pitch / cycle)) : 0;
      if (n > 0)
      {
        for (Size i = n; i < ms1_rt.size(); i += n)
        {
          cell_edges.push_back(ms1_rt[i]);
        }
        if (! cell_edges.empty() && 2 * (ms1_rt.size() - cell_edges.size() * n) < n)
        {
          cell_edges.pop_back(); // remainder < pitch/2 joins the last cell
        }
      }
      for (double e : cell_edges)
      {
        char buf[32];
        std::snprintf(buf, sizeof buf, "%.17g", e);
        cell_boundaries += (cell_boundaries.empty() ? "" : ",") + String(buf);
      } // round-trip exact
      writeLogInfo_("[cell] pitch " + String(cell_pitch) + " s = " + String(n) + " MS1 frames (cycle " + numStr(cycle)
                    + " s): " + String(cell_edges.size() + 1) + " cells, edges " + (cell_boundaries.empty() ? String("none") : cell_boundaries));
    }
    const double splitter_scan_time = cycle;
    if (splitter_scan_time > 0.0)
    {
      writeLogInfo_("[split] scan time " + numStr(splitter_scan_time) + " s from " + String(ms1_rt.size())
                    + " MS1 frames -> chrom_fwhm / scan time is the splitter's window");
    }
    else
    {
      writeLogWarn_("[split] the MS1 frame times show no positive gap; valley splitting uses each trace's own cycle time.");
    }
    const double max_span = getDoubleOption_("trace:max_span_sec");
    int ms1_bands = getIntOption_("perf:ms1_trace_bands");
    if (ms1_bands < 1) { ms1_bands = 1; }
    double _t_ms1 = phaseClock(), _c_ms1 = cpuSeconds(); // [perf-instr]
    // detectTraces_ spawns its bands as tasks, so this phase opens its own parallel region.
    vector<Trace> ms1_traces;
    // The MS1 traces' frame table is the MS1 map's own spectrum sequence (RT-sorted above); their
    // spans live in this store, which must outlive scoring.
    TraceStore ms1_store;
    {
      vector<uint32_t> ri;
      ri.reserve(ms1_map.size());
      for (const auto& sp : ms1_map)
      {
        ri.push_back(rtIndex(sp.getRT()));
      }
      ms1_store.setFrames(ri, nullptr);
    }
    // [minlen] trace:min_length_sec < 0 (auto): the run's minimum MS1 trace length from its own peak widths. A fixed
    // sample of the MS1 map -- every peak in the 1-Th m/z slices [25k, 25k + 1), plus the band halo -- is traced like
    // the run (tolerances, snr, valley splitter, bands) but with no minimum length, so the estimate does not depend on
    // the threshold it sets; a whole Th per slice samples every mass-defect phase. The median FWHM W of the brightest
    // tenth of the traces centred in a slice (ranked by apex intensity, then FWHM: deterministic) sets the minimum to
    // trace:min_length_fwhm x W, at least half a cycle (2 MS1 points). Fewer than 100 measurable traces, or no cycle,
    // fall back to 3 s. A fixed value skips all of it and leaves the header unchanged.
    double ms1_min_len = getDoubleOption_("trace:min_length_sec");
    String minlen_stamp; // spx:ms1_min_length, auto only
    if (ms1_min_len < 0.0)
    {
      const double frac = getDoubleOption_("trace:min_length_fwhm"), ppm = getDoubleOption_("trace:mass_error_ppm");
      const double PITCH = 25.0, WIDTH = 1.0;
      auto slice_pos = [&](double x) { return x - PITCH * std::floor(x / PITCH); }; // position within its 25-Th period
      PeakMap smp;
      TraceStore smp_store;
      smp_store.setFrames(ms1_store.rt_index, nullptr);
      vector<Trace> st;
      Size smp_pk = 0;
      const Size og0 = rtOffGrid().load(); // [soa-rt] counts the run's traces, not the sample's
      std::exception_ptr smp_err;
#pragma omp parallel
#pragma omp single
      try
      {
        // the cut: chunks of spectra in parallel, each spectrum kept in place (empty ones too: the frame sequence)
        vector<MSSpectrum> part(ms1_map.size());
        const int nchunk = std::max(1, std::min((int)ms1_map.size(), 4 * omp_get_max_threads()));
        std::exception_ptr cut_err; // [task-err]
#pragma omp taskloop grainsize(1) default(shared)
        for (int ci = 0; ci < nchunk; ++ci)
        {
          try
          {
            const Size a = (Size)ci * ms1_map.size() / nchunk, b = (Size)(ci + 1) * ms1_map.size() / nchunk;
            for (Size i = a; i < b; ++i)
            {
              const MSSpectrum& s = ms1_map[i];
              const auto* ima = s.getFloatDataArrays().empty() ? nullptr : &s.getFloatDataArrays()[0];
              MSSpectrum& t = part[i];
              t.setRT(s.getRT());
              t.setMSLevel(1);
              OpenMS::DataArrays::FloatDataArray ia;
              ia.setName(Constants::UserParam::ION_MOBILITY);
              for (Size k = 0; k < s.size(); ++k)
              {
                const double x = s[k].getMZ(), h = x * ppm * 1e-6 * BAND_HALO_PPM_MUL, r = slice_pos(x);
                if (r < WIDTH + h || r >= PITCH - h)
                {
                  t.push_back(s[k]);
                  if (ima && k < ima->size()) { ia.push_back((*ima)[k]); }
                }
              }
              t.getFloatDataArrays().push_back(std::move(ia));
            }
          }
          catch (...)
          {
            taskErr(cut_err);
          }
        }
        if (cut_err) { std::rethrow_exception(cut_err); }
        smp.reserveSpaceSpectra(part.size());
        for (auto& sp : part)
        {
          smp_pk += sp.size();
          smp.addSpectrum(std::move(sp));
        }
        vector<MSSpectrum>().swap(part);
        st = detectTraces_(smp, smp_store, delta_im, getDoubleOption_("trace:noise_threshold_int"), ms1_snr, 0.0,
                           getDoubleOption_("trace:ms1_split_valleys"), ms1_bands, splitter_scan_time);
      }
      catch (...)
      {
        smp_err = std::current_exception();
      }
      if (smp_err) { std::rethrow_exception(smp_err); }
      rtOffGrid().store(og0);
      vector<pair<double, double>> iw; // (apex intensity, FWHM s) of the traces centred in a slice
      for (const Trace& t : st)
      {
        if (std::isfinite(t.mz) && slice_pos(t.mz) < WIDTH) { iw.emplace_back(intensityOf(smp_store, t), xicFwhm(smp_store, t)); }
      }
      std::sort(iw.begin(), iw.end(), [](const pair<double, double>& a, const pair<double, double>& b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;
      });
      const Size top = iw.size() / 10;
      vector<double> w;
      for (Size i = 0; i < top; ++i)
      {
        if (iw[i].second > 0.0) { w.push_back(iw[i].second); }
      }
      std::sort(w.begin(), w.end());
      auto q = [&w](double p) { return w[(Size)(p * (double)(w.size() - 1) + 0.5)]; };
      const String head = "sample " + String(smp_pk) + " MS1 peaks in 1-Th slices every 25 Th -> " + String(st.size()) + " traces at min length 0, "
                          + String(iw.size()) + " centred in a slice; brightest tenth " + String(top) + ", " + String(w.size())
                          + " with a FWHM (>= 3 points)";
      if (w.size() >= 100 && cycle > 0.0)
      {
        const double w50 = q(0.5);
        ms1_min_len = std::max(frac * w50, 0.5 * cycle);
        const int pts = 1 + (int)std::min(std::ceil(ms1_min_len / cycle - 1e-9), 1e9); // capped before the int conversion
        writeLogInfo_("[minlen] auto: " + head + "; FWHM median " + numStr(w50) + " s = " + numStr(w50 / cycle) + " cycles (p10 " + numStr(q(0.1))
                      + " s = " + numStr(q(0.1) / cycle) + ", p90 " + numStr(q(0.9)) + " s = " + numStr(q(0.9) / cycle) + "), cycle " + numStr(cycle)
                      + " s -> trace:min_length_sec " + numStr(ms1_min_len) + " = max(" + numStr(frac) + " x median, half a cycle): >= " + String(pts)
                      + " MS1 points");
        minlen_stamp = "auto " + String(ms1_min_len) + " s = max(" + String(frac) + " x FWHM median " + String(w50) + " s, half the cycle "
                       + String(cycle) + " s), >= " + String(pts) + " MS1 points; FWHM of " + String(w.size()) + " bright sampled traces";
      }
      else
      {
        ms1_min_len = 3.0;
        writeLogWarn_("[minlen] auto: " + head + ", cycle " + numStr(cycle)
                      + " s: too few traces to measure the peak width (100 and a positive cycle time needed); using trace:min_length_sec 3.");
        minlen_stamp = "auto fallback 3 s (" + String(w.size()) + " bright sampled traces with a FWHM, cycle " + String(cycle) + " s)";
      }
    }
    // [ms1-prune] every MS1 spectrum must have left its m/z sample before losing a peak, or the band edges
    // would come from a pruned map. Checked here, outside the parallel region an exception cannot leave.
    if (ms1_prune.on && (ms1_prune.frames != ms1_map.size() || ms1_prune.sample.size() != ms1_prune.sample_n))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "[ms1-prune] the pre-prune m/z sample covers " + String(ms1_prune.frames) + " spectra ("
                                      + String(ms1_prune.sample.size()) + " of " + String(ms1_prune.sample_n) + " m/z) but the MS1 map holds "
                                      + String(ms1_map.size()) + "; -perf:ms1_prune false keeps every peak");
    }
    // An exception must not leave the parallel region (terminate): captured inside the single, rethrown after it. This
    // covers detectTraces_'s own code, the [det] diagnostics included, and its taskloop tasks, which rethrow here ([task-err]).
    std::exception_ptr ms1_err;
#pragma omp parallel
#pragma omp single
    try
    {
      ms1_traces = detectTraces_(ms1_map, ms1_store, delta_im, getDoubleOption_("trace:noise_threshold_int"), ms1_snr, ms1_min_len,
                                 getDoubleOption_("trace:ms1_split_valleys"), ms1_bands, splitter_scan_time, &ms1_prune);
    }
    catch (...)
    {
      ms1_err = std::current_exception();
    }
    if (ms1_err) { std::rethrow_exception(ms1_err); }
    {
      const double cap = max_span;
      // [outer] per trace, reading the store and writing only its own span fields: in parallel
      std::exception_ptr trim_err;
      if (cap > 0.0)
      {
#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)ms1_traces.size(); ++i)
        {
          try
          {
            trimToSpan(ms1_store, ms1_traces[(Size)i], cap);
          }
          catch (...)
          {
            taskErr(trim_err);
          }
        }
      }
      if (trim_err) { std::rethrow_exception(trim_err); }
    }
    vector<double>().swap(ms1_prune.sample); // [ms1-prune] consumed by the band edges, or unused at one band
    if (detOn()) { writeLogInfo_("[det] MS1 traces n=" + String(ms1_traces.size()) + " digest=" + String(traceDigest(ms1_store, ms1_traces))); }
    phaseAdd("MS1_TRACE", _t_ms1, _c_ms1, false);

    double iso_im_tol = getDoubleOption_("charge:iso_im_tolerance"); // [isotol] 0 = auto, resolved after the content sort
    String iso_tol_stamp = String(iso_im_tol);
    const double ms2_noise = getDoubleOption_("trace:ms2_noise_threshold_int");
    writeLogInfo_("MS2 trace gate: apex > " + String(ms2_snr * ms2_noise) + " (snr " + String(ms2_snr) + " x noise " + String(ms2_noise)
                  + "), membership > " + String(ms2_noise) + ", min_length " + String(ms2_minlen) + " s");
    // [determinism] inferPrecursors_ is greedy, so its output depends on ms1_traces' container order, which
    // EPD's parallel split makes thread-dependent: sort by content first. A non-finite coordinate would make
    // the comparator non-transitive (UB in sort), so those traces are dropped before it.
    {
      const Size before = ms1_traces.size();
      ms1_traces.erase(remove_if(ms1_traces.begin(), ms1_traces.end(),
                                 [&ms1_store](const Trace& t) {
                                   return ! std::isfinite(t.mz) || ! std::isfinite(rtOf(ms1_store, t)) || ! std::isfinite(t.im)
                                          || ! std::isfinite(intensityOf(ms1_store, t));
                                 }),
                       ms1_traces.end());
      if (before != ms1_traces.size()) { writeLogInfo_("dropped " + String(before - ms1_traces.size()) + " MS1 traces with non-finite coordinates"); }
    }
    // [tile] Content-only key, never a container index or arena offset, so the order does not depend on
    // which tile detected a trace. Equal-key neighbours are counted: std::sort is not stable.
    auto trace_less = [&ms1_store](const Trace& a, const Trace& b) {
      if (a.mz != b.mz) { return a.mz < b.mz; }
      {
        const double ra = rtOf(ms1_store, a), rb = rtOf(ms1_store, b);
        if (ra != rb) { return ra < rb; }
      }
      if (a.im != b.im) { return a.im < b.im; }
      {
        const double ia = intensityOf(ms1_store, a), ib = intensityOf(ms1_store, b);
        if (ia != ib) { return ia < ib; }
      }
      if (a.npts != b.npts) { return a.npts < b.npts; }
      if (a.len != b.len) { return a.len < b.len; }
      if (ms1_store.frame_rt[a.frame0] != ms1_store.frame_rt[b.frame0]) { return ms1_store.frame_rt[a.frame0] < ms1_store.frame_rt[b.frame0]; }
      // equal summaries: the profile decides, so equivalent traces are identical ones and std::sort's tie order cannot matter
      if (a.tof != b.tof) { return a.tof < b.tof; }
      if (a.apex != b.apex) { return a.apex < b.apex; }
      const float* pa = ms1_store.inten.data() + a.off;
      const float* pb = ms1_store.inten.data() + b.off;
      return std::lexicographical_compare(pa, pa + a.len, pb, pb + b.len);
    };
    {
      // [outer] in parallel when that is provably std::sort's order (no equal-key neighbours), else std::sort itself
      const Size ties = sortExact(ms1_traces, trace_less);
      writeLogInfo_("[det] MS1 traces sorted by content: " + String(ms1_traces.size()) + " traces, " + String(ties) + " equal-key neighbours");
    }
    if (iso_im_tol <= 0.0)
    {
      iso_im_tol = autoIsoImTol_(ms1_traces, ms1_store, max_charge, delta_rt, mass_ppm, iso_tol_stamp); // [isotol]
    }
    else
    {
      writeLogInfo_("[isotol] 1/K0 tolerance for isotope partners " + String(iso_im_tol) + " (fixed)");
    }
    std::optional<Phase> win_phase;                      // [perf-instr] the window loop (function scope)
    double _t_inf = phaseClock(), _c_inf = cpuSeconds(); // [perf-instr]
    vector<MonoV2In> mono_v2in;                          // [monocheck-v2] the count walks, parallel to `precursors` until monoPositionPass_
    vector<PrecursorCand> precursors
      = inferPrecursors_(ms1_traces, ms1_store, max_charge, delta_rt, iso_im_tol, mass_ppm, cfg_mono_check >= 2 ? &mono_v2in : nullptr);
    if (cfg_adaptive_auto) // [adaptive] auto: decide the run's gate from its precursor span distribution (MS1 points)
    {
      Size nf = 0;
      for (const PrecursorCand& pc : precursors)
      {
        if ((int)ms1_traces[pc.trace_idx].np() <= cfg_faint_frames) { ++nf; }
      }
      const double frac = precursors.empty() ? 0.0 : (double)nf / (double)precursors.size();
      cfg_adaptive = frac >= cfg_auto_faint_frac;
      if (! cfg_adaptive)
      {
        cfg_shrink = cfg_coloc = cfg_short_conc = false; // auto:off = the plain gate for every precursor, not 'off'
      }
      writeLogInfo_("[adaptive] auto: " + String(nf) + " of " + String(precursors.size()) + " precursors span <= " + String(cfg_faint_frames)
                    + " MS1 points (" + String(frac) + ", threshold " + String(cfg_auto_faint_frac) + "): faint gate "
                    + (cfg_adaptive ? "ON" : "OFF"));
    }
    if (detOn())
    {
      vector<uint64_t> b;
      for (const auto& pc : precursors)
      {
        pushBits(b, pc.rt);
        pushBits(b, pc.im);
      }
      writeLogInfo_("[det] precursors n=" + String(precursors.size()) + " digest=" + String(detDigest(std::move(b))));
    }
    // [im-veto] Charge halving, corrected by ion mobility. See charge:im_charge_veto.
    if (getStringOption_("charge:im_charge_veto") == "true")
    {
      auto fit_line = [&](int z, double& a, double& b, double& sig, Size& n) -> bool {
        std::map<int, vector<double>> bins; // 25-Th bins of 1/K0, confident calls only
        for (const auto& pc : precursors)
        {
          if (pc.charge == z && pc.n_isotopes >= 3 && std::isfinite(pc.im) && pc.mono_mz > 0.0) { bins[(int)(pc.mono_mz / 25.0)].push_back(pc.im); }
        }
        vector<double> xs, ys;
        n = 0;
        for (auto& kv : bins)
        {
          auto& v = kv.second;
          if (v.size() < 30) { continue; }
          std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
          xs.push_back((kv.first + 0.5) * 25.0);
          ys.push_back(v[v.size() / 2]);
          n += v.size();
        }
        if (xs.size() < 4) { return false; }
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        const double m = (double)xs.size();
        for (Size i = 0; i < xs.size(); ++i)
        {
          sx += xs[i];
          sy += ys[i];
          sxx += xs[i] * xs[i];
          sxy += xs[i] * ys[i];
        }
        const double den = m * sxx - sx * sx;
        if (den <= 0.0) { return false; }
        b = (m * sxy - sx * sy) / den;
        a = (sy - b * sx) / m;
        vector<double> res;
        for (const auto& pc : precursors)
        {
          if (pc.charge == z && pc.n_isotopes >= 3 && std::isfinite(pc.im) && pc.mono_mz > 0.0) { res.push_back(fabs(pc.im - (a + b * pc.mono_mz))); }
        }
        if (res.empty()) { return false; }
        std::nth_element(res.begin(), res.begin() + res.size() / 2, res.end());
        sig = 1.4826 * res[res.size() / 2];
        return sig > 0.0;
      };
      double a2 = 0, b2 = 0, s2 = 0, a3 = 0, b3 = 0, s3 = 0;
      Size n2 = 0, n3 = 0;
      const bool ok2 = fit_line(2, a2, b2, s2, n2), ok3 = fit_line(3, a3, b3, s3, n3);
      const double kband = getDoubleOption_("charge:im_veto_band");
      Size n1 = 0, to2 = 0, to3 = 0;
      if (ok2)
      {
        for (auto& pc : precursors)
        {
          if (pc.charge != 1 || ! std::isfinite(pc.im)) { continue; }
          ++n1;
          const double r2 = fabs(pc.im - (a2 + b2 * pc.mono_mz)) / s2;
          const double r3 = ok3 ? fabs(pc.im - (a3 + b3 * pc.mono_mz)) / s3 : 1e9;
          if (r2 <= kband && r2 <= r3)
          {
            pc.charge = 2;
            ++to2;
          }
          // On the 3+ band: unset the charge; it becomes a guessed precursor, which is dropped.
          else if (ok3 && r3 <= kband)
          {
            pc.charge = 0;
            ++to3;
          }
        }
        writeLogInfo_("[im-veto] 2+ line 1/K0 = " + numStr(a2) + " + " + numStr(b2) + " * m/z (n=" + String(n2) + ", sigma " + numStr(s2) + ")"
                      + (ok3 ? "; 3+ line " + numStr(a3) + " + " + numStr(b3) + " * m/z (n=" + String(n3) + ", sigma " + numStr(s3) + ")"
                             : String("; no 3+ line"))
                      + ". Of " + String(n1) + " z=1 calls, " + String(to2) + " (" + numStr(100.0 * (double)to2 / (double)std::max<Size>(n1, 1))
                      + "%) sit on the 2+ band and are re-called 2+, " + String(to3) + " on the 3+ band are dropped, " + String(n1 - to2 - to3)
                      + " kept as 1+");
      }
      else
      {
        writeLogInfo_("[im-veto] too few confident 2+ calls to fit a mobility line; no z=1 re-calls");
      }
    }
    // [monocheck-v2] At the final charge: after the veto (which reads the walks' own monoisotopes, as without v2),
    // before the charge-less precursors are dropped (the veto edits in place, so the order still matches).
    if (cfg_mono_check >= 2)
    {
      monoPositionPass_(ms1_traces, ms1_store, precursors, mono_v2in, delta_rt, iso_im_tol, mass_ppm);
      vector<MonoV2In>().swap(mono_v2in);
    }
    // [emission] isotope support is required: drop the guessed (charge 0) singletons.
    {
      const Size before = precursors.size();
      precursors.erase(std::remove_if(precursors.begin(), precursors.end(), [](const PrecursorCand& pc) { return pc.charge == 0; }),
                       precursors.end());
      writeLogInfo_("[emission] dropped " + String(before - precursors.size()) + " guessed precursors, " + String(precursors.size())
                    + " isotope-supported remain");
    }
    // charge:min_charge filters the OUTPUT: precursors below it are assembled as usual (their isotope envelopes are what
    // the foreign-isotope removal of the other spectra looks for) and their spectra are dropped before the write
    const int zmin = getIntOption_("charge:min_charge");
    std::atomic<Size> n_below_zmin {0}; // written by the (possibly pipelined) tile writer
    // [ms1-funnel] diag:dump_ms1_tsv: MS1 traces and precursors, to attribute where a precursor is lost.
    const String ms1_dump = getStringOption_("diag:dump_ms1_tsv");
    if (! ms1_dump.empty())
    {
      ofstream ft((ms1_dump + ".traces.tsv").c_str());
      // [model] 12 significant digits: the default 6 quantises m/z 765 to ~1.3 ppm.
      ft.precision(12);
      ft << "mz\trt\tim\tintensity\tfwhm_s\tn_xic\trt_first\trt_last\n"; // [tile] span endpoints
      for (const auto& t : ms1_traces)
      {
        ft << t.mz << '\t' << rtOf(ms1_store, t) << '\t' << t.im << '\t' << intensityOf(ms1_store, t) << '\t' << xicFwhm(ms1_store, t) << '\t'
           << t.np() << '\t' << (t.span() ? rtAtSpan(ms1_store, t, 0) : 0.0) << '\t' << (t.span() ? rtAtSpan(ms1_store, t, t.span() - 1) : 0.0)
           << '\n';
      }
      ft.close();
      ofstream fp((ms1_dump + ".precursors.tsv").c_str());
      fp.precision(12);
      // [model] trace_idx makes the hypothesis -> seed-trace link explicit, so an offline replica
      // of inferPrecursors_ can be checked row-for-row rather than assumed to agree.
      fp << "mono_mz\tcharge\trt\tim\tguessed\ttrace_idx\n";
      for (const auto& pc : precursors)
      {
        fp << pc.mono_mz << '\t' << pc.charge << '\t' << pc.rt << '\t' << pc.im << '\t' << 0 << '\t' << pc.trace_idx << '\n';
      }
      fp.close();
      writeLogInfo_("MS1 funnel dump: " + String(ms1_traces.size()) + " traces, " + String(precursors.size()) + " precursors -> " + ms1_dump
                    + ".{traces,precursors}.tsv");
    }
    // [twin-rt] assembly:twin_rt_frames: the twin fold's frame tolerance, |dRT| <= (n + 0.5) cycles (0 s = one frame).
    // Auto (0): n = the median FWHM of the final precursors' seed traces (>= 3 points), in cycles, rounded, >= 1. The
    // median of a multiset does not depend on the container order, so neither does n.
    double twin_rt_tol = 0.0;
    String twin_rt_desc;
    // the fold joins frames only through v2's extension: without v2 no spectrum carries the mark, so the relation
    // is the same-frame one, and it is logged and stamped as such
    const bool twin_rt_inert = cfg_mono_check < 2;
    if (twin_im >= 0.0 && twin_rt_frames >= 0 && cycle > 0.0 && twin_rt_inert)
    {
      writeLogInfo_("[twins] assembly:twin_rt_frames: the same MS1 frame only (the cross-frame fold needs charge:mono_profile_check v2)");
    }
    if (twin_im >= 0.0 && twin_rt_frames >= 0 && cycle > 0.0 && ! twin_rt_inert)
    {
      int n = twin_rt_frames;
      String how = String(n) + " (set)";
      if (n == 0)
      {
        vector<double> fw;
        fw.reserve(precursors.size());
        for (const PrecursorCand& pc : precursors)
        {
          const double f = xicFwhm(ms1_store, ms1_traces[pc.trace_idx]);
          if (f > 0.0) { fw.push_back(f); }
        }
        double med = -1.0;
        if (! fw.empty())
        {
          const Size m = fw.size() / 2;
          std::nth_element(fw.begin(), fw.begin() + (long)m, fw.end());
          med = fw[m];
        }
        // rounded in double and capped before the int conversion (a near-zero cycle would overflow int)
        n = med > 0.0 ? (int)std::max(1.0, std::min(std::round(med / cycle), 1e6)) : 1;
        char b[200];
        std::snprintf(b, sizeof b, "%d (auto: median MS1 FWHM %.3f s over %zu seed traces = %.2f cycles)", n, med, fw.size(), med / cycle);
        how = String(b);
      }
      twin_rt_tol = (n + 0.5) * cycle;
      char b[96];
      std::snprintf(b, sizeof b, "|dRT| <= %.4f s", twin_rt_tol);
      twin_rt_desc = "MS1 frames up to " + how + " cycles apart (" + String(b) + ", cycle " + String(cycle) + " s"
                     + ", across frames precursor m/z within 0.05 ppm and one monoisotope placed by the v2 extension), one RT cell";
      cfg_twin_ext_mark = true;
      writeLogInfo_("[twins] assembly:twin_rt_frames: the fold joins " + twin_rt_desc);
    }
    // sort precursors by mono m/z so each window can binary-search its members instead of
    // rescanning all ~1M precursors per window (total order for determinism).
    // [outer] trace_idx is unique (one precursor per committed seed), so the order is total and sortExact sorts in
    // parallel; an equal-key pair (a non-finite key) would send it to std::sort itself.
    sortExact(precursors, [](const PrecursorCand& a, const PrecursorCand& b) {
      if (a.mono_mz != b.mono_mz) { return a.mono_mz < b.mono_mz; }
      if (a.rt != b.rt) { return a.rt < b.rt; }
      if (a.im != b.im) { return a.im < b.im; }
      return a.trace_idx < b.trace_idx;
    });
    vector<double> prec_mz(precursors.size());
    for (Size i = 0; i < precursors.size(); ++i)
    {
      prec_mz[i] = precursors[i].mono_mz;
    }
    auto prec_range = [&prec_mz](double lo, double hi) { // the precursors with mono m/z in [lo, hi]
      return std::pair<Size, Size> {(Size)(lower_bound(prec_mz.begin(), prec_mz.end(), lo) - prec_mz.begin()),
                                    (Size)(upper_bound(prec_mz.begin(), prec_mz.end(), hi) - prec_mz.begin())};
    };
    // [mem] Scoring reads only the profiles of traces a precursor points at: release the rest.
    // (ms1_map is already empty: detectTraces_ consumed it.)
    {
      vector<bool> needed(ms1_traces.size(), false);
      for (const auto& pc : precursors)
      {
        if (pc.trace_idx < ms1_traces.size()) { needed[pc.trace_idx] = true; }
      }
      // Compacted in place, in offset order (see compactUnreferenced()).
      const Size freed = compactUnreferenced(ms1_traces, ms1_store, needed);
      writeLogInfo_("Released XICs of " + String(freed) + " unreferenced MS1 traces." + rss());
#ifdef __GLIBC__
      if (trim_on) { trimLogged_("MS1 -> window loop"); }
#endif
    }
    writeLogInfo_(clockStamp() + " Detected " + String(ms1_traces.size()) + " MS1 traces -> " + String(precursors.size()) + " precursor hypotheses."
                  + rss());
    std::exception_ptr worker_err;

    //-------------------------------------------------------------
    // Window loop, tile by tile: the master admits windows, their work runs as tasks, results land in
    // index-addressed slots and each tile is canonically sorted, so output is thread-count-independent.
    // Exceptions are captured and rethrown serially.
    //-------------------------------------------------------------
    vector<pair<WinKey, PeakSlab*>> window_list;
    window_list.reserve(ms2_by_window.size());
    for (auto& kv : ms2_by_window)
    {
      window_list.emplace_back(kv.first, &kv.second);
    }
    const Size n_win = frozen.size(); // the run's windows (the resident list is empty under the tiled source)
    vector<WinKey> fkeys;
    fkeys.reserve(n_win);
    for (const auto& kv : frozen)
    {
      fkeys.push_back(kv.first);
    }
    // [heavy-first] Biggest windows (compact bytes) first, so the largest is not left alone in the tail; output is sorted per tile.
    auto heavy_first = [](const auto& a, const auto& b) { return a.second->bytes() > b.second->bytes(); };
    phaseAdd("PRECURSOR_INFER", _t_inf, _c_inf, false);
    win_phase.emplace("WINDOW_LOOP");
    std::stable_sort(window_list.begin(), window_list.end(), heavy_first);
    // MEMORY vs SPEED: every concurrently-processed window holds its own frames + traces + grid, so
    // peak RSS scales with the number of windows in flight, not with the window count. Capping
    // concurrency trades wall time for RAM. 0 = unlimited. [mem]
    const int max_conc = getIntOption_("perf:max_concurrent_windows");
    const int n_threads = omp_get_max_threads();
    // [bands] perf:trace_bands splits each window's m/z range into bands traced concurrently (the window
    // count is fixed by the method); 0 = auto, enough bands that windows x bands covers every thread.
    const int bands_opt = getIntOption_("perf:trace_bands");
    int n_bands = bands_opt;
    if (n_bands <= 0) // auto: enough bands that windows x bands covers every thread
    {
      n_bands = (int)std::ceil((double)n_threads / std::max<Size>(n_win, 1));
    }
    n_bands = std::max(1, n_bands);
    // [master/worker] n_conc only caps windows IN FLIGHT (perf:max_concurrent_windows, 0 = uncapped); the
    // free-RAM admission gate is the real bound.
    int n_conc = (int)n_win;
    if (max_conc > 0) { n_conc = std::min(n_conc, max_conc); }
    if (n_conc < 1) { n_conc = 1; }
    // outer loop over windows + inner loop over bands = two live levels
    if (n_bands > 1) { omp_set_max_active_levels(2); }
    writeLogInfo_("Parallelism: " + String(n_win) + " windows x " + String(n_bands) + " m/z bands = " + String((int)n_win * n_bands)
                  + " trace units; " + String(n_threads) + " worker threads over one task pool; up to " + String(n_conc)
                  + " windows in flight (memory-bounded)");
    bool integer_detector = integer_req;
    // A candidate tdf that could not be read at the preflight (e.g. still being copied) is tried once more.
    if (integer_detector && ! tofAxis().ok)
    {
      for (const std::string& c : tdf_cand)
      {
        if (loadTofAxis(c, why))
        {
          writeLogInfo_("flight-time axis from " + c + ": " + String(tofAxis().b_by_frame.size()) + " frame factors");
          break;
        }
      }
    }
#ifdef DIASPEXTRACT_WITH_MZPEAK
    // No tdf candidate calibrated the run: an mzPeak archive carries its own vendor/analysis.tdf.gz, and the
    // exact-m/z path has already read MzCalibration + Frames.T1 from it: the same calibration the .d path loads.
    if (integer_detector && ! tofAxis().ok)
    {
      if (const auto mzp = DIAspeXtract::lastMzPeakMeta(); mzp && mzp->exact)
      {
        const DIAspeXtract::MzPeakExactMz& ex = *mzp->exact;
        setTofAxis(ex.cal, ex.t1_by_frame);
        tofAxis().acq_lo = ex.acq_lo;
        tofAxis().acq_hi = ex.acq_hi;
        tofAxis().n_bins = ex.n_bins;
        tofAxis().bounds_why = ex.bounds_why;
        writeLogInfo_("flight-time axis from the mzPeak archive's embedded analysis.tdf: " + String(tofAxis().b_by_frame.size()) + " frame factors");
      }
    }
#endif
    // FAIL CLOSED: a frame without a parseable "frame=<Id>", or one the calibration table lacks, cannot be calibrated per
    // frame, and the integer detector must not calibrate it at reference temperature -- fall back to openms instead.
    if (integer_detector && tofAxis().ok && cstat.unmapped_frame > 0)
    {
      writeLogWarn_(
        "trace:detector integer needs a calibration for every frame, but " + String(cstat.unmapped_frame)
        + " frame(s) have no \"frame=<Id>\" in their nativeID (an mzML converted from a .d may lack it) or are missing from the calibration table. "
          "Falling back to trace:detector openms. Pass -trace:detector openms to silence this.");
      integer_detector = false;
    }
    if (integer_detector && ! tofAxis().ok)
    {
      // Fall back rather than refuse: the OpenMS detector needs no calibration, and spx:detector records which ran.
      writeLogWarn_("trace:detector integer needs the Bruker flight-time calibration (MzCalibration and Frames.T1 in "
                    "analysis.tdf), which this input lacks ("
                    + why
                    + "). Falling back to trace:detector openms, a different "
                      "algorithm, so the result is not comparable with an integer-detector run. Pass -trace:detector openms to "
                      "silence this.");
      integer_detector = false;
    }
    const double mem_frac = getDoubleOption_("perf:mem_fraction");
    std::atomic<Size> inflight {0};
    std::atomic<int> adm_hwm {0}; // high-water mark of admitted windows, for the log and the e2e
    // One malloc_trim when the last window has materialised returns the load's bookkeeping left on
    // glibc's free lists (one tile only).
    std::atomic<Size> n_materialised {0};
    // Projected footprint per window: ~44 B per compact peak sustained + 256 MiB.
    auto project = [&](const PeakSlab& sl) { return sl.peaks() * (Size)44 + (Size)256 * 1024 * 1024; };
    for (auto& kv : frozen)
    {
      FrozenWin& fw = kv.second;
      fw.b.resize(fw.frames());
      for (Size f = 0; f < fw.frames(); ++f)
      {
        fw.b[f] = tofAxis().factor(fw.frame_id[f]); // the value toTof() gives the slab
      }
      fw.nearest_local = nearestLocal(fw.frame_rt, ms1_store.frame_rt, delta_rt);
    }
    if (cell_pitch > 0.0 && ! integer_detector)
    {
      writeLogWarn_("tile:rt_sec is ignored: trace:detector openms traces whole windows, without retention-time cells.");
      cell_edges.clear();
      cell_boundaries.clear();
      cell_pitch = -1.0;
    }
    if (cell_pitch > 0.0 && bands_opt <= 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "perf:trace_bands 0 (auto) would make the output depend on the thread count, so it cannot be combined with "
                                        "tile:rt_sec > 0. Set a fixed perf:trace_bands.");
    }
    if (tile_source && integer_detector && bands_opt <= 0)
    {
      throw Exception::InvalidParameter(
        __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "perf:trace_bands 0 (auto) cannot be used when a .d is streamed (perf:stream_load true). Set a fixed perf:trace_bands.");
    }
    // [tile] The cells grouped into TILES. A tile owns the precursors whose MS1 frame time lies in
    // [lo - tile_drt, hi - tile_drt) -- the half-open ranges partition RT (lo = -inf for the first,
    // hi = +inf for the last) -- and traces its windows on the tile's frames only; a precursor's
    // candidates (|rtOf - rt| <= tile_drt) then lie in [lo - 2 tile_drt, hi): the tile's own
    // fragments plus those CARRIED from the previous tile (rtOf >= lo - 2 tile_drt, re-carried
    // transitively). One tile = the untiled loop, statement for statement. tile_drt is the WIDEST RT
    // window any precursor's gate may use: gate:faint_delta_rt under gate:adaptive (a faint precursor
    // 3-6 s from a cut line otherwise sees different candidates per tile grid).
    const double tile_drt = cfg_adaptive ? std::max(delta_rt, cfg_faint_drt) : delta_rt;
    struct Tile
    {
      double lo, hi;
      bool last;
    };
    vector<Tile> tiles;
    const int cells_per_tile = getIntOption_("tile:cells_per_tile");
    {
      const Size n_cells_run = cell_edges.size() + 1;
      const bool one_piece = (out_type != FileTypes::MZML); // mzPeak output is bulk-only: one tile
      const Size per = (cells_per_tile <= 0 || cell_pitch <= 0.0 || one_piece) ? n_cells_run : (Size)cells_per_tile;
      if (one_piece && cells_per_tile > 0 && cell_pitch > 0.0 && n_cells_run > 1)
      {
        writeLogWarn_("[tile] mzPeak output is written in one piece: the " + String(n_cells_run)
                      + " cells run as one tile and need the whole run's memory. Write mzML (-out_type mzML) to run in tiles.");
      }
      for (Size c0 = 0; c0 < n_cells_run; c0 += per)
      {
        const Size c1 = std::min(n_cells_run, c0 + per);
        tiles.push_back(Tile {c0 == 0 ? -std::numeric_limits<double>::infinity() : cell_edges[c0 - 1],
                              c1 == n_cells_run ? std::numeric_limits<double>::infinity() : cell_edges[c1 - 1], c1 == n_cells_run});
      }
      String tl;
      for (const Tile& t : tiles)
      {
        tl += (tl.empty() ? "" : " ") + String("[") + numStr(t.lo) + "," + numStr(t.hi) + ")";
      }
      writeLogInfo_("[tile] " + String(tiles.size()) + " tile(s) of up to " + String(per) + " cell(s): " + tl);
    }
    // [tile] Run-level per-window band edges from the acquisition, so every tile traces in the same bands. A resident
    // slab's bins are converted here, once, in place (toTof is idempotent).
    map<WinKey, pair<TofIdx, TofIdx>> wbands;
    if (integer_detector)
    {
      if (! (tofAxis().n_bins > 1))
      {
        throw Exception::InvalidParameter(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "trace:detector integer needs the acquisition m/z range and digitizer size (GlobalMetadata), which this input's calibration source lacks ("
            + tofAxis().bounds_why + "). Use -trace:detector openms.");
      }
      for (const WinKey& k : fkeys)
      {
        wbands[k] = {0xFFFFFFFFu, 0};
      }
      std::exception_ptr band_err; // [task-err] a throw escaping the region would abort without a message
#pragma omp parallel for schedule(dynamic, 1)
      for (long i = 0; i < (long)fkeys.size(); ++i)
      {
        try
        {
          const WinKey& k = fkeys[(Size)i];
          const auto rs = ms2_by_window.find(k); // a resident slab: converted here, once, in place
          if (rs != ms2_by_window.end()) { toTof(rs->second); }
          // tlo 0, thi the larger of the digitizer's last bin (DigitizerNumSamples - 1: no bin index can exceed it)
          // and the acquisition's upper m/z through the calibration at every frame of the frozen table
          // (n_bins > 1: checked above)
          const TofIdx last = (TofIdx)(tofAxis().n_bins - 1);
          TofIdx thi_mz = 0, thi_q = last;
          for (double b : frozen.at(k).b)
          {
            if (tofAxis().acq_hi > 0.0) { thi_mz = std::max(thi_mz, tofAxis().tofOf(tofAxis().acq_hi, b)); }
            // the last bin's m/z, through the store's 1e-5 quantum, back to a bin: on a coarse
            // calibration the round trip can land one bin higher
            thi_q = std::max(thi_q, tofAxis().tofOf(tofAxis().mzOf(last, b) + 1.0e-5, b));
          }
          const TofIdx thi_acq = std::max(thi_q, thi_mz);
          if (thi_acq >= 4294967294u)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "[tof] the acquisition's upper bin is not representable (the partition would wrap)", String(thi_acq));
          }
#pragma omp critical(winlog)
          writeLogInfo_("[tof] window " + winLabel(k) + " acquisition tlo=0 thi=" + String(thi_acq) + " (digitizer " + String(tofAxis().n_bins - 1)
                        + ", MzAcqRangeUpper " + String(thi_mz) + ")");
          wbands.at(k) = pair<TofIdx, TofIdx> {0, thi_acq};
        }
        catch (...)
        {
          taskErr(band_err);
        }
      }
      if (band_err) { std::rethrow_exception(band_err); }
    }
    // [tile] carried boundary fragments (records + spans), in for tile k = out of tile k - 1; and
    // the once-only ownership check per precursor (a window's m/z range is tile-independent)
    struct Carry
    {
      vector<Trace> rec;
      vector<float> inten;
    };
    // [pipe] One slot per window: `c` is carry_out(w, k - 1), and `ready` counts the window's tiles that have
    // published theirs, so tile k may take `c` once ready == k. Under the per-tile barrier that always holds;
    // with pipelined tiles it is the one dependency between a window's tiles -- tile k + 1 may slab, prepare,
    // detect, merge and split while tile k still runs; it waits only where it takes the carry (assembly).
    struct CarrySlot
    {
      Carry c;
      std::atomic<Size> ready {0};
    };
    // [pipe] A window's handle on its slot for one tile. take() waits until the window's previous tile has
    // published, then moves that carry out; publish() hands this tile's carry (put() before it, or none) to
    // the next tile. The destructor publishes on every exit of the window body, waiting first if the body left
    // before take(), so a slot is never published out of order and the next tile never waits forever.
    struct CarryGate
    {
      CarrySlot& s;
      const Size tk;
      bool took = false, done = false;
      void wait() const // [pipe] idle until the carry arrives: help the tile reader or writer meanwhile
      {
        while (s.ready.load(std::memory_order_acquire) < tk)
        {
          if (! teamAssist(4)) { std::this_thread::sleep_for(std::chrono::microseconds(200)); }
        }
      }
      Carry take()
      {
        wait();
        took = true;
        Carry c = std::move(s.c);
        s.c = Carry();
        return c;
      }
      void put(Carry&& c)
      { s.c = std::move(c); }
      void publish()
      {
        if (done) { return; }
        if (! took) { (void)take(); } // the body left early: drop the previous tile's carry, in order
        done = true;
        s.ready.store(tk + 1, std::memory_order_release);
      }
      ~CarryGate()
      { publish(); }
    };
    map<WinKey, CarrySlot> carry;
    map<WinKey, uint64_t> slab_chain; // [det] the per-window slab digest, chained across tiles
    map<WinKey, vector<uint8_t>> visited;
    for (const WinKey& k : fkeys)
    {
      carry[k];
      const auto [plo, phi] = prec_range(winBounds(k).first, winBounds(k).second);
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
    {
      out_exp.setMetaValue("spx:mz_calibration", DIAspeXtract::lastMzPeakCalibration()); // not BrukerTimsFile's state
    }
    else
#endif
    {
      out_exp.setMetaValue("spx:mz_calibration", BrukerTimsFile::lastMzCalibration());
    }
    out_exp.setMetaValue("spx:detector", String(integer_detector ? "integer" : "openms"));
    { // the test-only switches (docs/options.md): each changes the output or the run, so a file made under one says so
      String set;
      for (const char* v : {"DIASPEXTRACT_TILE_NO_CARRY", "DIASPEXTRACT_MS1_PRUNE_NO_WITNESS", "DIASPEXTRACT_MS1_PRUNE_NO_CHAIN", "DIASPEXTRACT_DET",
                            "DIASPEXTRACT_TASK_THROW"})
      {
        if (const char* e = std::getenv(v)) { set += (set.empty() ? "" : " ") + String(v) + "=" + e; }
      }
      if (! set.empty())
      {
        writeLogWarn_("Test-only variables are set; the output is not the tool's: " + set);
        out_exp.setMetaValue("spx:test_env", set);
      }
    }
    out_exp.setMetaValue("spx:tile_rt_sec", cell_pitch);
    out_exp.setMetaValue("spx:tile_cells_per_tile", cells_per_tile);
    out_exp.setMetaValue("spx:band_edges", integer_detector ? String("acquisition") : String("none")); // what was applied
    out_exp.setMetaValue("spx:tile_source", tile_source ? "stream" : "resident");
    out_exp.setMetaValue("spx:tiles", (int)tiles.size());
    out_exp.setMetaValue("spx:frame_table", frame_table_src); // the frozen tables' source (loader | mzpeak | slab)
    out_exp.setMetaValue("spx:pearson_G", "frames");
    out_exp.setMetaValue("spx:tile_boundaries", cell_boundaries.empty() ? String("none") : cell_boundaries);
    out_exp.setMetaValue("spx:corr_power", cfg_corr_power);
    out_exp.setMetaValue("spx:gate_adaptive",
                         cfg_adaptive_auto ? String(cfg_adaptive ? "auto:on" : "auto:off") : String(cfg_adaptive ? "faint" : "off"));
    out_exp.setMetaValue("spx:gate_faint", String(cfg_faint_frames) + " frames, " + String(cfg_faint_drt) + " s, power " + String(cfg_faint_cp));
    out_exp.setMetaValue("spx:gate_shrinkage", cfg_shrink ? String("fisher r0 ") + String(cfg_shrink_r0) + " sd " + String(cfg_shrink_sd) + " from "
                                                              + String(cfg_shrink_min_pts)
                                                          : String("off"));
    out_exp.setMetaValue("spx:gate_colocation", cfg_coloc ? String("ranges ") + String(cfg_coloc_overlap) : String("off"));
    out_exp.setMetaValue("spx:gate_short_traces", cfg_short_conc
                                                    ? String("concentric tol ") + String(cfg_short_center_tol) + " from " + String(cfg_short_min_pts)
                                                    : String("off"));
    out_exp.setMetaValue("spx:cap_rank", String("evidence"));
    out_exp.setMetaValue("spx:same_bin_merge",
                         String("im_slices ")
                           + (! integer_detector          ? String("n/a (openms detector)")
                              : cfg_im_slice_merge >= 0.0 ? String(cfg_im_slice_merge) + " (span " + String(2.0 * cfg_im_slice_merge) + ")"
                                                          : String("off"))
                           + ", assembly " + (cfg_merge_same_bin ? (integer_detector ? "on" : "on (exact m/z)") : "off"));
    out_exp.setMetaValue("spx:im_weight_sigma", cfg_im_weight_sigma);
    // [dnoise] what was applied: the port, "off", or skipped because the input is no .d
    out_exp.setMetaValue("spx:dnoise_ms1", dnoise.on ? String(String(DIAspeXtract::Dnoise::PORTED_FROM) + " MS1")
                                                     : String(dnoise_req ? "skipped (input is not a Bruker .d)" : "off"));
    out_exp.setMetaValue("spx:dnoise_ms1_params", dnoise.on ? dnoise.params() : String("none"));
    out_exp.setMetaValue("spx:dnoise_ms1_points", dnoise.on ? String("raw " + String(dnoise.total.raw) + " kept " + String(dnoise.total.kept)
                                                                     + " approximate " + String(dnoise.total.ambiguous))
                                                            : String("none"));
    out_exp.setMetaValue("spx:require_isotope_support", 1);
    if (! minlen_stamp.empty())
    {
      out_exp.setMetaValue("spx:ms1_min_length", minlen_stamp); // [minlen] auto only: a fixed value leaves the header unchanged
    }
    if (cfg_mono_check)
    {
      out_exp.setMetaValue("spx:monocheck", cfg_mono_check_stamp); // [monocheck] only when on
    }
    if (getDoubleOption_("charge:iso_im_tolerance") <= 0.0) // [isotol] auto only: a fixed value leaves no stamp
    {
      out_exp.setMetaValue("spx:iso_im_tolerance", iso_tol_stamp + ", first partner in m/z");
    }
    out_exp.setMetaValue("spx:twin_merge", twin_im < 0.0 ? String("off")
                                           : twin_rt_desc.empty()
                                             ? String("same frame, same charge, 20 ppm, 1/K0 <= " + String(twin_im) + ", fused 10 ppm")
                                             : String(twin_rt_desc + ", same charge, 20 ppm, 1/K0 <= " + String(twin_im)
                                                      + ", fused 10 ppm")); // [twin-rt] the same-frame text without the cross-frame fold
    if (cfg_foreign_w < 1.0)                                                // [foreign] stamped only when on (1 = off leaves no stamp)
    {
      out_exp.setMetaValue("spx:foreign",
                           String("weight " + String(cfg_foreign_w) + ", " + String(cfg_foreign_ppm) + " ppm, 1/K0 " + String(cfg_foreign_dim)
                                  + ", min_iso " + String(cfg_foreign_min_iso) + ", positions mono..M+3, own ion kept ("
                                  + (cfg_foreign_ppm > 30.0 ? String(cfg_foreign_ppm) : String("30")) + " ppm)")); // the own-ion tests' max(ppm, 30)
    }
    if (cfg_fuse_ppm > 0.0) // [fuse] stamped only when on (0 = off leaves no stamp)
    {
      String s = String("anchored ") + String(cfg_fuse_ppm) + " ppm";
      s += ", summed intensity, intensity-weighted m/z, cap key from the summed apex; after the same-bin merge, before min_fragments and "
           "max_fragments";
      out_exp.setMetaValue("spx:fragment_fuse", s);
    }
    std::unique_ptr<TileWriter> writer;
    // written as `<out>.part` and renamed on completion: a tile that throws must not leave a
    // structurally valid file holding the earlier tiles only
    const String out_part = out + ".part";
    if (out_type == FileTypes::MZML)
    {
      writer.reset(new TileWriter(out_part));
      writer->setExperimentalSettings(out_exp);
    }
    vector<MSSpectrum> all_out; // the bulk (mzPeak) path only
    Size n_out = 0;
    Size n_twins = 0;               // [twins] spectra folded into a twin, over every tile
    std::array<Size, 3> twin_st {}; // [twin-rt] mergeTwins's counters, over every tile
    writeLogInfo_("Processing " + String(n_win) + " isolation windows in " + String(tiles.size()) + " tile(s) from the "
                  + (tile_source ? "streaming" : "resident") + " source (OpenMP over windows, <=" + String(n_conc)
                  + " concurrent, admission bounded by " + String((int)(mem_frac * 100)) + "% of free RAM)..." + rss());
    ProgressLogger tile_progress; // one step per tile, set by the tile writer (serial, in tile order)
    tile_progress.setLogType(log_type_);
    tile_progress.startProgress(0, (SignedSize)tiles.size(), "extracting pseudo-spectra (tiles)");
    if (const Size cg = cgroupHeadroom(); cg != SIZE_MAX)
    {
      writeLogInfo_("[dyn-mem] cgroup memory limit: " + String((double)cg / (1ull << 30), false)
                    + " GiB of headroom, used as the free RAM where it is lower");
    }

    // [master/worker] The master walks the windows heaviest-first and admits them against memory BEFORE
    // spawning, so the byte budget and the window cap hold by construction; every unit of work inside a
    // window is an OpenMP task, so idle threads steal across windows. Results go to index-addressed slots
    // (win_out[wi], pslot[j], per[ci][b]), never appended in completion order.
    std::atomic<int> active {0};
    // [pipe] Pipelined tiles. With the per-tile barrier, every tile waits for its slowest window chain and the
    // read, sort, merge, write and trim between two tiles overlap nothing. Pipelined, ONE parallel region
    // spans the whole loop and its single thread drives it: it reads tile k + 1 and admits its windows while
    // tile k's still run, and finishes (sorts, merges, writes, trims) the tiles strictly in tile order as each
    // completes. A window's tiles depend on each other only through the carry (CarryGate); tiles own disjoint
    // precursors and are sorted and written per tile, so the output is the barrier's, bit for bit.
    // DIASPEXTRACT_PIPE_TILES=<n> (default 2, 1 = the per-tile barrier) bounds the tiles whose windows may run
    // at once: tile j is read only after every window of tile j - n has finished, which bounds the slabs
    // resident and the chain of carry waits; the window cap and the memory admission below still apply to
    // every window in flight, over all tiles. The reader's and writer's loops run as tasks of the window team
    // (teamFor); only the vendor reader's own decode loop (libOpenMS) opens a nested team.
    const int pipe_depth = []() {
      const char* e = std::getenv("DIASPEXTRACT_PIPE_TILES");
      const long v = e && *e ? std::atol(e) : 2L;
      return (int)std::max(1L, std::min(v, 16L));
    }();
    const bool pipe = pipe_depth > 1 && tiles.size() > 1 && n_threads > 1;
    std::atomic<bool> pipe_abort {false}; // a window or the driver failed: admit nothing more, skip unstarted windows
    writeLogInfo_("[pipe] "
                  + (pipe ? "tiles pipelined, up to " + String(pipe_depth) + " with windows in flight"
                          : String("off (")
                              + (pipe_depth <= 1     ? "DIASPEXTRACT_PIPE_TILES=1"
                                 : tiles.size() <= 1 ? "one tile"
                                                     : "one thread")
                              + "): a barrier after every tile"));
    // [pipe] One tile in flight: its slabs, its windows heaviest-first, their output slots, and how many of its
    // windows have not finished yet (the driver finishes the tile once that is 0).
    struct TileRun
    {
      Size tk = 0;
      map<WinKey, PeakSlab> tile_slabs;
      vector<pair<WinKey, PeakSlab*>> wl;
      vector<vector<MSSpectrum>> win_out;
      Size win_done = 0;
      std::atomic<Size> left {0};
      Size claim_next = 0, claim_end = 0; // [pipe] its windows' claims in `claims` (the driver's own bookkeeping)
    };
    // [pipe] A window admitted on the pipelined loop is CLAIMED by whoever starts it first: its task, or the driver while it
    // waits for that window's tile (run_claimed below). The claims outlive the tiles (a task whose window the driver
    // ran may still be queued when the tile is written); a deque keeps their addresses while the driver appends.
    struct WinClaim
    {
      TileRun* R;
      long wi;
      Size need;
      std::atomic<bool> taken {false};
      WinClaim(TileRun* r, long w, Size n): R(r), wi(w), need(n)
      {
      }
    };
    std::deque<WinClaim> claims;
    auto load_tile = [&](TileRun& tile_run) {
      const Size tk = tile_run.tk;
      const Tile& tile = tiles[tk];
      // [tile] the tile's slabs: one tile = the resident slabs themselves (moved out by the body);
      // several = a slice of each window's frames in [lo, hi) (the resident source keeps the whole
      // run: this arm is the exactness gate, not the memory configuration)
      map<WinKey, PeakSlab>& tile_slabs = tile_run.tile_slabs;
      vector<pair<WinKey, PeakSlab*>>& wl = tile_run.wl;
      if (tile_source)
      {
        // [stream] read the tile: a fresh consumer over the loader's RT range, lockstep-padded against
        // the frozen tables restricted to the tile, so a window's slab is exactly the table's frames
        // in [lo, hi); a window with none gets an empty slab (the body's empty-slab branch)
        Phase _ph("LOAD(tile)");
        const double _tr0 = phaseClock();
        PeakMap ms1_none;
        PickCompactConsumer tc(spicker, tile_slabs, ms1_none, cstat);
        tc.tablefn = [&]() {
          PickCompactConsumer::FrameTables t;
          for (const auto& kv : frozen)
          {
            PickCompactConsumer::FrameTab& tab = t[kv.first];
            for (Size f = 0; f < kv.second.frames(); ++f)
            {
              if (kv.second.frame_rt[f] >= tile.lo && kv.second.frame_rt[f] < tile.hi)
              {
                tab.rt.push_back(kv.second.frame_rt[f]);
                tab.id.push_back(kv.second.frame_id[f]);
              }
            }
          }
          return t;
        };
        BrukerTimsFile::Config c2 = cfg;
        c2.load_ms1 = false;
        c2.dia_ms2_rt_lo = tile.lo;
        c2.dia_ms2_rt_hi = tile.hi;
        // [pipe] On the pipelined driver the frames are decoded and sliced by the window team (teamFor: free threads and
        // windows at their break points help), not by a nested team that would put more threads than cores on the
        // machine (libomp spins then: Windows, 12 threads, 4/8 nested decode threads: 8:33/10:19 instead of 7:58 min).
        if (omp_in_parallel() && getIntOption_("perf:decode_threads") <= 0)
        {
          c2.dia_ms2_parallel_for = [](long n, const std::function<void(long, int)>& body) { teamFor(n, 1, body); };
          c2.dia_ms2_parallel_for_slots = omp_get_num_threads();
        }
        BrukerTimsFile().loadDIAStreaming(in, tc, c2); // a failure here is fatal: no fallback onto a half-read tile
        tc.finish();
        Size n_del = 0, n_pad = 0, n_absent = 0, n_peaks = 0;
        for (const auto& kv : tc.tables())
        {
          PeakSlab& sl = tile_slabs[kv.first]; // a window that delivered nothing still owns its table's frames
          if (sl.frames() == 0 && ! kv.second.rt.empty())
          {
            // every frame of the slice was skipped (an empty ion-mobility slice in each): pad them
            // all, exactly as the resident path does, so the frame sequence is the table's
            ++n_absent;
            for (Size f = 0; f < kv.second.rt.size(); ++f)
            {
              sl.append(kv.second.rt[f], kv.second.id[f], nullptr, nullptr, nullptr, 0, kv.second.rt.size());
            }
            n_pad += kv.second.rt.size();
          }
          else
          {
            n_del += kv.second.delivered;
            n_pad += kv.second.padded;
          }
          if (sl.rt != kv.second.rt || sl.frame_id != kv.second.id)
          {
            throw InputRefusal(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[3b] a tile's slab is not the frozen table's frames for the tile",
                               winLabel(kv.first));
          }
          n_peaks += sl.peaks();
        }
        for (auto& kv : tile_slabs)
        {
          indexRt(kv.second);
        }
        if (detOn() && integer_detector) // the OpenMS path needs the slabs in m/z quanta (materializeWindow)
        {
          // [det] the reader's exactness: the per-(window, tile) slab digest, chained over the tiles.
          // The tiles partition a window's frames in RT order and the digest is a sequential fold, so
          // the chain of a window's tiles EQUALS the whole-window digest a resident run prints -- the
          // check that the tiled reads are the resident peaks, frame for frame. The conversion to
          // flight-time bins has to happen first (it fills `b` and consumes `frame_id`); it is
          // idempotent, so the window body's own call becomes a no-op.
          for (auto& kv : tile_slabs)
          {
            toTof(kv.second);
            // the chain must START from the same FNV offset basis the whole-window digest starts
            // from, or it cannot equal it (a zero-initialised map entry would seed the fold with 0)
            uint64_t& h = slab_chain.try_emplace(kv.first, 1469598103934665603ULL).first->second;
            h = slabDigest(kv.second, h);
            writeLogInfo_("[det] window " + winLabel(kv.first) + " tile " + String(tk + 1) + " slab n=" + String(kv.second.peaks())
                          + " digest=" + String(slabDigest(kv.second)) + " chain=" + String(h));
          }
        }
        writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + " [" + numStr(tile.lo) + "," + numStr(tile.hi)
                      + ") read: " + String(tile_slabs.size()) + " windows, frames delivered " + String(n_del) + " padded " + String(n_pad)
                      + " absent-windows " + String(n_absent) + ", " + String(n_peaks) + " peaks, " + numStr(phaseClock() - _tr0) + " s" + rss());
      }
      else if (tiles.size() > 1)
      {
        Size tile_peaks = 0;
        for (const auto& w : window_list)
        {
          const PeakSlab& sl = *w.second;
          auto before = [](uint32_t i, double t) { return rtAxis()[i] < t; }; // rt_index ascends (indexRt)
          const auto rb = sl.rt_index.begin(), re = sl.rt_index.end();
          const Size f0 = std::lower_bound(rb, re, tile.lo, before) - rb;
          const Size f1 = std::lower_bound(rb + f0, re, tile.hi, before) - rb;
          PeakSlab& ts = tile_slabs[w.first];
          ts = sliceFrames(sl, f0, f1);
          tile_peaks += ts.peaks();
        }
        writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + " [" + String(tile.lo) + "," + String(tile.hi)
                      + "): " + String(tile_slabs.size()) + " windows, " + String(tile_peaks) + " peaks resident" + rss());
      }
      if (tile_source || tiles.size() > 1)
      {
        for (auto& kv : tile_slabs)
        {
          wl.emplace_back(kv.first, &kv.second);
        }
        std::stable_sort(wl.begin(), wl.end(), heavy_first);
      }
      else
      {
        wl = window_list;
      }
      tile_run.win_out.resize(wl.size());
      tile_run.left.store(wl.size());
    }; // load_tile
    // [pipe] the window body, defined below the admission that spawns it
    std::function<void(TileRun*, long, Size)> run_window;
    // [pipe] Set by the pipelined driver: runs one unstarted window of the tile it waits for (or helps the reader or
    // writer); true if it did. Called while admission is blocked, so the driver works instead of sleeping.
    std::function<bool()> admit_help;
    auto admit_tile = [&](TileRun* tile_run) {
      auto* rw = &run_window;
      tile_run->claim_next = tile_run->claim_end = claims.size();
      for (long wi = 0; wi < (long)tile_run->wl.size() && ! pipe_abort.load(); ++wi)
      {
        const Size need = project(*tile_run->wl[wi].second);
        for (;;)
        {
          const Size budget = (Size)(availableBytes() * mem_frac);
          // [dyn-mem] admit while fewer than n_conc windows are in flight and the projection fits mem_fraction
          // of MemAvailable (re-read at every admission); an empty pipeline always admits, however large the window
          if ((active.load() == 0) || (active.load() < n_conc && inflight.load() + need <= budget)) { break; }
          // A team of ONE (-threads 1, OMP_THREAD_LIMIT=1, OMP_DYNAMIC) has nobody else to run the windows admitted so
          // far, and libgomp's taskyield is a no-op: once the window cap or the memory budget held, this loop would spin
          // forever. In a team of one every queued window is a child of this thread, so taskwait runs them, which releases
          // their admission charges.
          if (omp_get_num_threads() < 2)
          {
#pragma omp taskwait
            continue;
          }
          // [pipe] Off the pipeline, yield: the driver must not be handed a window here (a no-op in libgomp, not in every runtime).
          // On it, the window cap is the window count, so the next tile's windows are admitted only as the current tile's
          // finish: the driver blocked here while ONE team mate ran windows (-threads 2: 1.25 busy cores). Run one instead.
          // clang-format off: cpplint wants braces on both sides of else; the OpenMS style would put the block on one line
          if (! pipe)
          {
#pragma omp taskyield
          }
          else if (admit_help && admit_help())
          {
            continue;
          }
          // clang-format on
          std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
        inflight += need;
        ++active;
        if (active.load() > adm_hwm.load()) adm_hwm.store(active.load());
        if (pipe)
        {
          WinClaim* c = &claims.emplace_back(tile_run, wi, need);
          tile_run->claim_end = claims.size();
#pragma omp task default(none) firstprivate(rw, c)
          if (! c->taken.exchange(true, std::memory_order_acq_rel)) (*rw)(c->R, c->wi, c->need);
        }
        else
        {
#pragma omp task default(none) firstprivate(rw, tile_run, wi, need)
          (*rw)(tile_run, wi, need);
        }
      }
    }; // admit_tile
    run_window = [&](TileRun* tile_run, const long wi, const Size need) {
      const Size tk = tile_run->tk;
      const Tile& tile = tiles[tk];
      vector<pair<WinKey, PeakSlab*>>& wl = tile_run->wl;
      vector<vector<MSSpectrum>>& win_out = tile_run->win_out;
      Size& win_done = tile_run->win_done;
      // A lambda, because `continue` is not a legal exit from a task's structured block and the
      // window body has one early exit (a window whose fragments all fail the finiteness filter).
      [&]() {
        try
        {
          const double win_lo = winBounds(wl[wi].first).first;
          const double win_hi = winBounds(wl[wi].first).second;
          // admitted by the master before this task existed; release both charges however we leave
          struct Release
          {
            std::atomic<Size>& b;
            std::atomic<int>& a;
            Size n;
            ~Release()
            {
              b -= n;
              --a;
            }
          } rel {inflight, active, need};
          // [pipe] this window's carry slot: its destructor publishes this tile's carry on every exit
          CarryGate gate {carry.at(wl[wi].first), tk};
          if (pipe_abort.load())
          {
            return; // [pipe] a window or the driver failed: the run is ending
          }
          // The completion line, printed on EVERY exit of the body (empty window, any assembly arm, exception).
          struct WinDone
          {
            std::function<void()> f;
            ~WinDone()
            { f(); }
          } win_done_guard {[&]() {
#pragma omp critical(winlog)
            writeLogInfo_("  window " + String(++win_done) + "/" + String(wl.size()) + " " + String(win_lo) + "-" + String(win_hi) + ": spectra "
                          + String(win_out[wi].size()) + clockStamp());
          }};
          // [tile] take ONLY this window's slab (a move: the store is the slab) and convert its
          // m/z quanta to flight-time bins in place. Only the OpenMS path needs the peaks converted
          // back into Peak1D doubles, which is the largest per-window allocation.
          PeakMap wmap;
          PeakSlab wslab;
          if (integer_detector)
          {
            wslab = std::move(*wl[wi].second);
            toTof(wslab);
            if (detOn())
            {
#pragma omp critical
              writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " slab n=" + String(wslab.peaks())
                            + " digest=" + String(slabDigest(wslab)));
            }
          }
          else
          {
            wmap = materializeWindow(*wl[wi].second);
            if (detOn()) // [det] compact store digest BEFORE tracing: splits "pick" from "trace+assemble"
            {
              // materializeWindow consumed the frames; digest the materialised map instead (same content)
              vector<uint64_t> b;
              for (const auto& sp : wmap)
              {
                for (const auto& pk : sp)
                {
                  pushBits(b, pk.getMZ());
                  float f = pk.getIntensity();
                  uint32_t w;
                  std::memcpy(&w, &f, 4);
                  b.push_back(w);
                }
              }
#pragma omp critical
              writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " compact n=" + String(b.size() / 2)
                            + " digest=" + String(detDigest(std::move(b))));
            }
          }
          if (++n_materialised == wl.size() && trim_on && tiles.size() == 1) trimLogged_("last window materialised -> store bookkeeping returned");
          wmap.sortSpectra();
          vector<Trace> frag_traces;
          TraceStore wst; // this window's frame table + arenas; outlives scoring
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
            wst.setFrames(fwin.rt_index, &fwin.b); // the FROZEN table, not the resident slab's
            // Banded like the OpenMS path. A band SEEDS only in its own flight-time core but may
            // extend anywhere, so no halo is needed and no trace is cut; the one approximation is
            // that two bands can take the same peak (each has its own `visited`), which is the same
            // class of error the OpenMS halo carries. Bands are
            // TASKS, so they fill the pool alongside every other window's work.
            // [tile] the window's RUN-LEVEL band edges (from the acquisition), the same in every tile
            const TofIdx tlo = wbands.at(wl[wi].first).first, thi = wbands.at(wl[wi].first).second;
            const int nb = (thi <= tlo) ? 1 : n_bands;
            {
              // The band edges are computed ONCE and shared by the seed segmentation and the detector,
              // so the two can never disagree about which band owns a seed.
              vector<TofIdx> bnd;
              if (nb > 1)
              {
                const double wdt = (double)(thi - tlo + 1) / nb;
                bnd.resize(nb + 1);
                for (int b = 0; b < nb; ++b)
                  bnd[b] = tlo + (TofIdx)(b * wdt);
                bnd[nb] = thi + 1;
              }
              else
              {
                bnd = {0, 0};
              }
              // [cell] Detection per fixed RT cell: the window's frames are cut at the run-level cell
              // edges (tile:rt_sec), each cell is traced on a slice of its own -- seeds, `visited`,
              // frame_live and the tof index are the cell's -- and every cell's band stores are absorbed
              // into the window store afterwards, in (cell, band) order, with the frame indices rebased
              // by the cell's first frame. One cell (tile:rt_sec=-1) is the whole window on the slab
              // itself: the same statements, no copy, byte-identical to the uncelled code.
              // the cut is SLAB-local (sliceFrames indexes the slab); the rebase below is a VALUE
              // lookup of each cell's first frame time in the frozen table
              vector<double> slab_rt(wslab.frames());
              for (Size f = 0; f < wslab.frames(); ++f)
                slab_rt[f] = rtAxis()[wslab.rt_index[f]];
              vector<Size> cut {0};
              for (double e : cell_edges)
              {
                const Size f = (Size)(lower_bound(slab_rt.begin(), slab_rt.end(), e) - slab_rt.begin());
                if (f > cut.back() && f < slab_rt.size()) cut.push_back(f);
              }
              cut.push_back(slab_rt.size());
              const Size n_cells = cut.size() - 1;
              const bool one_cell = (n_cells == 1);
              // Each band of each cell appends spans to a PRIVATE store; after the taskloops join, the
              // window store absorbs them in (cell, band) order and every offset is rebased.
              vector<vector<TraceStore>> bst(n_cells, vector<TraceStore>(nb));
              vector<vector<vector<Trace>>> per(n_cells, vector<vector<Trace>>(nb));
              for (Size ci = 0; ci < n_cells; ++ci)
              {
                const Size f0 = cut[ci], f1 = cut[ci + 1];
                PeakSlab cslab;
                if (! one_cell) cslab = sliceFrames(wslab, f0, f1);
                const PeakSlab& cs = one_cell ? wslab : cslab;
                if (! cs.tof.empty())
                {
                  // the window's band edges must enclose the cell's bins: band_of() and the tof index
                  // wrap for a bin below tlo. A run-level freeze must be a superset. In every mode, one
                  // cell included (the acquisition edges are metadata).
                  // [serial] the extremes as tasks: min and max are exact in any order
                  const Size np = cs.tof.size();
                  const long nr = parTasks(np, 1 << 18);
                  vector<TofIdx> clos(nr, 0xFFFFFFFFu), chis(nr, 0);
                  forTasks(nr, [&](long c) {
                    for (Size k = (Size)c * np / (Size)nr; k < (Size)(c + 1) * np / (Size)nr; ++k)
                    {
                      clos[c] = std::min(clos[c], cs.tof[k]);
                      chis[c] = std::max(chis[c], cs.tof[k]);
                    }
                  });
                  const TofIdx clo = *std::min_element(clos.begin(), clos.end()), chi = *std::max_element(chis.begin(), chis.end());
                  if (clo < tlo || chi > thi)
                  {
                    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "cell bins outside the window's band edges",
                                                          String(clo) + "-" + String(chi) + " vs " + String(tlo) + "-" + String(thi));
                  }
                }
                TracePrep tprep = prepareTracing(cs, ms2_noise, ms2_snr, bnd, tlo, thi); // once per cell
                if (nb == 1) per[ci][0] = detectTracesInteger(cs, tprep, bst[ci][0], mass_ppm, delta_im, ms2_noise, ms2_minlen);
                else
                {
                  vector<vector<Trace>>& pc = per[ci];
                  vector<TraceStore>& bc = bst[ci];
                  std::exception_ptr band_err; // [task-err] rethrown into the window body's catch
#pragma omp taskloop grainsize(1) default(shared)
                  for (int b = 0; b < nb; ++b)
                  {
                    try
                    {
                      teamAssist(1); // [pipe] a chunk of the pipelined tile reader or writer, if one is live
                      pc[b] = detectTracesInteger(cs, tprep, bc[b], mass_ppm, delta_im, ms2_noise, ms2_minlen, b);
                    }
                    catch (...)
                    {
                      taskErr(band_err);
                    }
                  }
                  if (band_err) std::rethrow_exception(band_err);
                }
              }
              // Each cell's first frame in the FROZEN table, by value (the slab-local cut is an index
              // into the resident slab, which may hold only some of the table's frames later); the
              // resident frames must be, frame for frame, the table's -- asserted per cell.
              vector<Size> f0z(n_cells);
              for (Size ci = 0; ci < n_cells; ++ci)
              {
                const double t0 = slab_rt[cut[ci]];
                f0z[ci] = (Size)(lower_bound(wst.frame_rt.begin(), wst.frame_rt.end(), t0) - wst.frame_rt.begin());
                for (Size j = 0; j < cut[ci + 1] - cut[ci]; ++j)
                {
                  if (f0z[ci] + j >= wst.frames() || wst.rt_index[f0z[ci] + j] != wslab.rt_index[cut[ci] + j])
                  {
                    throw OpenMS::Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                          "resident frames are not the frozen table's at cell " + String(ci) + " frame", String(j));
                  }
                }
              }
              // The slab is dead after the frame check above: release it BEFORE the merge, the window's peak allocation.
              vector<TofIdx>().swap(wslab.tof);
              vector<float>().swap(wslab.inten);
              vector<uint16_t>().swap(wslab.imq);
              vector<uint32_t>().swap(wslab.frame_off);
              vector<double>().swap(wslab.b);
              vector<uint32_t>().swap(wslab.rt_index);
              Size tot = 0, arena = 0;
              bool any_bins = false;
              for (Size ci = 0; ci < n_cells; ++ci)
              {
                for (int b = 0; b < nb; ++b)
                {
                  tot += per[ci][b].size();
                  arena += bst[ci][b].inten.size();
                  any_bins |= ! bst[ci][b].bins.empty();
                }
              }
              frag_traces.reserve(tot);
              wst.inten.reserve(arena);
              if (any_bins) wst.bins.reserve(arena); // as in splitIntegerTraces
              { // [serial] absorbed in (cell, band) order, rebased by the cell's first frame in the frozen table; copies as tasks
                vector<std::pair<TraceStore*, vector<Trace>*>> src;
                vector<uint32_t> fsh;
                for (Size ci = 0; ci < n_cells; ++ci)
                {
                  for (int b = 0; b < nb; ++b)
                  {
                    src.push_back({&bst[ci][b], &per[ci][b]});
                    fsh.push_back((uint32_t)f0z[ci]);
                  }
                }
                absorbAll(wst, frag_traces, src, fsh);
              }
            } // the band stores die here
            // Valley splitting, as on the OpenMS path; it rebuilds the arena, and the per-frame bins are dead after it.
            frag_traces = splitIntegerTraces(frag_traces, ms2_split, wst, n_bands * 4, splitter_scan_time);
            if (max_span > 0.0) // [serial] per trace, independent: as tasks over trace chunks
            {
              const Size n = frag_traces.size();
              const long nr = parTasks(n, 1 << 16);
              forTasks(nr, [&](long c) {
                for (Size i = (Size)c * n / (Size)nr; i < (Size)(c + 1) * n / (Size)nr; ++i)
                  trimToSpan(wst, frag_traces[i], max_span);
              });
            }
            vector<int16_t>().swap(wst.bins);
          }
          else
          {
            {
              vector<uint32_t> ri;
              ri.reserve(wmap.size()); // wmap is RT-sorted above
              for (const auto& sp : wmap)
                ri.push_back(rtIndex(sp.getRT()));
              wst.setFrames(ri, nullptr);
            }
            frag_traces = detectTraces_(wmap, wst, delta_im, ms2_noise, ms2_snr, ms2_minlen, ms2_split, n_bands, splitter_scan_time);
            if (max_span > 0.0)
            {
              for (auto& t : frag_traces)
                trimToSpan(wst, t, max_span);
            }
          }
          if (detOn())
          {
#pragma omp critical
            writeLogInfo_("[det] window " + String(win_lo) + "-" + String(win_hi) + " frag traces n=" + String(frag_traces.size())
                          + " digest=" + String(traceDigest(wst, frag_traces)));
          }
          // (compact frames already released inside materializeWindow)

          // [tile] the previous tile's boundary fragments join this window's list, spans copied
          // onto its arena; their frame indices are the frozen table's already (the same table in
          // every tile), and they enter BEFORE the content sort, so they take their canonical places.
          // [pipe] take() waits for this window's previous tile to publish them: the one point where a
          // window's tiles meet (everything above ran on this tile's slab alone).
          Size n_carried_in = 0;
          {
            const Carry ci = gate.take();
            if (! ci.rec.empty())
            {
              frag_traces.reserve(frag_traces.size() + ci.rec.size());
              for (const Trace& r : ci.rec)
              {
                if (wst.inten.size() + r.len >= (Size)std::numeric_limits<uint32_t>::max())
                  throw OpenMS::Exception::OutOfRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); // as absorbAll
                Trace t = r;
                t.off = (uint32_t)wst.inten.size();
                wst.inten.insert(wst.inten.end(), ci.inten.begin() + r.off, ci.inten.begin() + r.off + r.len);
                frag_traces.push_back(t);
              }
              n_carried_in = ci.rec.size();
            }
          }
          // [tile] the precursors this tile OWNS, decided with the SCORER's own arithmetic: a
          // fragment at frame time x is a candidate of rt iff fabs(x - rt) <= tile_drt (doubles), so
          // rt has no candidate at or beyond an edge T iff (T - rt) > tile_drt -- the same
          // subtraction, monotone in T. The first tile k with (T_{k+1} - rt) > tile_drt owns rt:
          // exactly one does (-inf never, +inf always). Marked here, before any early exit, so the
          // coverage check after the last tile sees every precursor exactly once. Only this
          // window's precursors (by mono m/z); NaN sorts to the end, excluded; charge 0 kept.
          const auto [plo, phi] = prec_range(win_lo, win_hi);
          vector<uint32_t> pidx;
          pidx.reserve(phi > plo ? phi - plo : 0);
          for (Size pi = plo; pi < phi; ++pi)
          {
            const PrecursorCand& pc = precursors[pi];
            if (! std::isfinite(pc.rt) || ! std::isfinite(pc.im)) continue;
            if (! (tk == 0 || ! (tile.lo - pc.rt > tile_drt))) continue; // an earlier tile owns it
            if (! (tile.last || (tile.hi - pc.rt > tile_drt))) continue; // a later tile owns it
            pidx.push_back((uint32_t)pi);
          }
          {
            vector<uint8_t>& vis = visited.at(wl[wi].first);
            for (uint32_t pi : pidx)
            {
              if (vis[pi - plo]++ != 0)
              {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[tile] a precursor is owned by two tiles",
                                              String(precursors[pi].rt, 6) + " s in window " + String(win_lo) + "-" + String(win_hi));
              }
            }
          }
          frag_traces.erase(
            remove_if(frag_traces.begin(), frag_traces.end(),
                      [&wst](const Trace& t) { return ! std::isfinite(t.im) || ! std::isfinite(t.mz) || ! std::isfinite(rtOf(wst, t)); }),
            frag_traces.end());
          if (frag_traces.empty()) return; // nothing traceable in this window

          // Content-only and TOTAL: traces equal in (im, mz, rt, apex intensity) fall through to the remaining
          // fields, so their order never depends on what else is resident; exact duplicates are interchangeable.
          // [serial] Sorted as tasks (parallelSort): it may order exact duplicates differently from std::sort, which
          // the comparator already allows -- the order of what else is resident moves them too.
          parallelSort(
            frag_traces.data(), frag_traces.size(),
            [&wst](const Trace& a, const Trace& b) {
              if (a.im != b.im) return a.im < b.im;
              if (a.mz != b.mz) return a.mz < b.mz;
              {
                const double ra = rtOf(wst, a), rb = rtOf(wst, b);
                if (ra != rb) return ra < rb;
              }
              {
                const double ia = intensityOf(wst, a), ib = intensityOf(wst, b);
                if (ia != ib) return ia < ib;
              }
              if (a.npts != b.npts) return a.npts < b.npts;
              if (a.len != b.len) return a.len < b.len;
              {
                const double fa = wst.frame_rt[a.frame0], fb = wst.frame_rt[b.frame0];
                if (fa != fb) return fa < fb;
              }
              // ...then the profile and the bin fields (two bands can share an apex with different flanks)
              for (Size k = 0; k < a.len; ++k)
              {
                const float xa = xv(wst, a, k), xb = xv(wst, b, k);
                if (xa != xb) return xa < xb;
              }
              if (a.tof != b.tof) return a.tof < b.tof;
              if (a.bframe != b.bframe) return a.bframe < b.bframe;
              if (a.rt_at != b.rt_at) return a.rt_at < b.rt_at;
              return false;
            },
            parTasks(frag_traces.size(), 1 << 15));
          // [tile] carry the boundary fragments to the next tile: every record (own or carried in)
          // whose RT is within 2 tile_drt of this tile's upper edge -- the candidates of a precursor
          // the next tile owns. Over ALL of frag_traces, so the carry is transitive by induction.
          // A precursor a later tile owns has (T - rt) <= tile_drt, its candidates x have
          // (rt - x) <= tile_drt: x >= T - 2 tile_drt up to two roundings of ~1e-12 s. The slack
          // over-carries by a hair (a carried fragment no precursor accepts changes nothing), never
          // under-carries.
          // [pipe] Taken and published HERE, after the content sort: scoring below only reads frag_traces
          // and the arena, so this is the carry the end of the body would produce, and the window's next
          // tile may assemble while this one scores.
          static const bool no_carry = std::getenv("DIASPEXTRACT_TILE_NO_CARRY") != nullptr; // e2e falsifier only
          if (! tile.last && ! no_carry)
          {
            Carry co;
            const double keep_from = tile.hi - 2.0 * tile_drt - 1e-6;
            for (const Trace& t : frag_traces)
            {
              if (rtOf(wst, t) >= keep_from)
              {
                Trace c = t;
                c.off = (uint32_t)co.inten.size();
                co.inten.insert(co.inten.end(), wst.inten.begin() + t.off, wst.inten.begin() + t.off + t.len);
                co.rec.push_back(c);
              }
            }
            if (n_carried_in || ! co.rec.empty())
            {
#pragma omp critical(winlog)
              writeLogInfo_("[tile] window " + String(win_lo) + "-" + String(win_hi) + " tile " + String(tk + 1) + ": carried in "
                            + String(n_carried_in) + ", carrying out " + String(co.rec.size()) + " fragments ("
                            + String((co.rec.size() * sizeof(Trace) + co.inten.size() * 4) / 1024) + " KB)");
            }
            gate.put(std::move(co));
          }
          gate.publish();
          vector<double> frag_im(frag_traces.size());
          {
            const Size n = frag_traces.size();
            const long nr = parTasks(n, 1 << 17); // [serial] as tasks
            forTasks(nr, [&](long c) {
              for (Size i = (Size)c * n / (Size)nr; i < (Size)(c + 1) * n / (Size)nr; ++i)
                frag_im[i] = frag_traces[i].im;
            });
          }
          FragRt frag_rt;
          frag_rt.build(frag_traces, wst, tile_drt, parTasks(frag_traces.size(), 1 << 17)); // [adaptive] buckets cover the wider window

          // Precompute the common correlation grid once per window; pdense is reused across the
          // window's precursors (this window runs on one thread). [perf]
          FragStats fg = buildFragStats(frag_traces, wst, fwin);
          vector<float> pdense(wst.frames(), 0.0f); // indexed by window-local frame
          // [foreign] the window's precursors by apex RT, owned by this tile or not: the same candidates in any tiling
          ForeignIdx_ fx;
          if (cfg_foreign_w < 1.0) fx.build(precursors, plo, phi);

          vector<MSSpectrum>& bucket = win_out[wi];
          {
            // Share-all (default, current best): each fragment goes to EVERY gated precursor.
            //
            // Precursors are the unit of parallelism (one window can hold a large share of them). Each result goes to
            // its own pre-sized slot, so task granularity cannot drop or reorder an emission; pdense is per-task scratch
            // and must stay firstprivate.
            vector<MSSpectrum> pslot(pidx.size()); // [tile] the owned precursors, in index order
            // [gomp-cliff] Bound the TASK COUNT: libgomp's GOMP_taskloop runs the whole loop inline in the
            // creating thread once team task_count + num_tasks > 64 x nthreads, silently serialising it
            // (grainsize does not bound the count). So <= n_threads worker tasks pull 32-precursor chunks from
            // a shared cursor; keep the team-wide total well under the threshold when adding task layers.
            const long n_prec_w = (long)pidx.size();
            const long n_task_w = std::max(1L, std::min((n_prec_w + 31) / 32, (long)n_threads));
            std::atomic<long> pcur {0};
            std::exception_ptr asm_err; // [task-err] rethrown into the window body's catch
#pragma omp taskloop num_tasks(n_task_w) firstprivate(pdense) default(shared)
            for (long w = 0; w < n_task_w; ++w)
            {
              try
              {
                taskFault("assemble");
                for (;;)
                {
                  teamAssist(1); // [pipe] a chunk of the pipelined tile reader or writer, if one is live
                  const long cs = pcur.fetch_add(32);
                  if (cs >= n_prec_w) break;
                  const long ce = std::min(cs + 32, n_prec_w);
                  for (long j = cs; j < ce; ++j)
                  {
                    const PrecursorCand& pc = precursors[pidx[(Size)j]];
                    MSSpectrum ms2;
                    assembleOne_(pc, win_lo, win_hi, frag_traces, wst, frag_im, frag_rt, ms1_traces, ms1_store, fg, gate_dim, delta_rt, min_corr,
                                 min_corr_pts, min_frags, max_frags, pdense, ms2, cfg_foreign_w < 1.0 ? &fx : nullptr);
                    if (! ms2.empty()) pslot[(Size)j] = std::move(ms2);
                  }
                }
              }
              catch (...)
              {
                taskErr(asm_err);
              }
            }
            if (asm_err) std::rethrow_exception(asm_err);
            for (auto& s : pslot)
              if (! s.empty()) bucket.push_back(std::move(s));
          }
        }
        catch (...)
        {
#pragma omp critical
          if (! worker_err) worker_err = std::current_exception();
          pipe_abort.store(true);
        }
      }(); // window body
      tile_run->left.fetch_sub(1, std::memory_order_acq_rel); // [pipe] last: the driver may finish the tile after this
    }; // run_window
    // [pipe] After every window of the tile has finished: the coverage check (last tile), the canonical sort,
    // the twin fold, the post-fusion cap and the write of the tile's block, then the between-tile trim.
    auto finish_tile = [&](TileRun& tile_run) {
      const Size tk = tile_run.tk;
      const Tile& tile = tiles[tk];
      vector<vector<MSSpectrum>>& win_out = tile_run.win_out;
      if (tile.last)
      {
        // [tile] coverage: every finite precursor of every window was owned by exactly one tile
        // (the per-tile check catches a second owner; this catches none)
        Size n_cov = 0;
        for (const WinKey& k : fkeys)
        {
          const auto [lo, hi] = winBounds(k);
          const auto [plo, phi] = prec_range(lo, hi);
          const vector<uint8_t>& vis = visited.at(k);
          for (Size pi = plo; pi < phi; ++pi)
          {
            if (! std::isfinite(precursors[pi].rt) || ! std::isfinite(precursors[pi].im)) { continue; }
            if (vis[pi - plo] != 1)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "[tile] a precursor was owned by no tile",
                                            String(precursors[pi].rt, 6) + " s in window " + String(lo) + "-" + String(hi));
            }
            ++n_cov;
          }
        }
        writeLogInfo_("[tile] coverage: " + String(n_cov) + " precursor-window pairs owned exactly once over " + String(tiles.size()) + " tile(s)");
        if (tile_source)
        {
          // the tile reads share `cstat` with pass 1; the pre-pass-2 report could not see them
          const Size lost2 = cstat.no_im_array + cstat.size_mismatch + cstat.bad_mz + cstat.bad_im;
          if (lost2)
          {
            writeLogWarn_("[compact] DROPPED " + String(lost2) + " peaks over the whole run (pass 1 + every tile): no_im_array="
                          + String(cstat.no_im_array) + " size_mismatch=" + String(cstat.size_mismatch) + " bad_mz=" + String(cstat.bad_mz)
                          + " bad_im=" + String(cstat.bad_im) + " (kept " + String(cstat.kept) + ")");
          }
          else
          {
            writeLogInfo_("[compact] no peak dropped over the whole run (kept " + String(cstat.kept) + ")");
          }
          if (detOn())
          {
            for (const auto& kv : slab_chain)
            {
              writeLogInfo_("[det] window " + winLabel(kv.first) + " slab chain=" + String(kv.second));
            }
          }
        }
        // MEMORY: every MS1 structure is dead once the last tile's window loop has joined -- scoring
        // was its only reader (assembleOne_ dereferences ms1_traces[pc.trace_idx]
        // and its span). Released BEFORE this tile's output is assembled and written. [mem]
        vector<Trace>().swap(ms1_traces);
        ms1_store = TraceStore();
        vector<PrecursorCand>().swap(precursors);
        vector<double>().swap(prec_mz);
        writeLogInfo_("[mem] MS1 traces, spans, precursors released after the window loop." + rss());
      }
      {
        // [tile] this tile's spectra in the canonical order: written as one block (mzML), or kept for
        // the bulk writer (mzPeak, one tile). Tiles are RT-disjoint, so the blocks concatenate to the
        // whole run's order.
        vector<MSSpectrum> tile_out;
        {
          Size n_tot = 0;
          for (const auto& bucket : win_out)
          {
            n_tot += bucket.size();
          }
          tile_out.reserve(n_tot);
        }
        for (auto& bucket : win_out)
        {
          for (auto& s : bucket)
          {
            const int z = s.getPrecursors().empty() ? 0 : s.getPrecursors()[0].getCharge();
            if (z > 0 && z < zmin)
            {
              ++n_below_zmin;
              continue;
            } // charge:min_charge
            tile_out.push_back(std::move(s));
          }
          bucket.clear();
          bucket.shrink_to_fit();
        }
        {
          Phase _ph("SORT(canonical)");
          canonicalSort(tile_out);
        }
        if (twin_im >= 0.0)
        {
          Phase _ph("MERGE(twins)");
          const Size before = tile_out.size(), gone = mergeTwins(tile_out, twin_im, twin_rt_tol, &cell_edges, tile_drt, &twin_st);
          n_twins += gone;
          writeLogInfo_("[twins] tile " + String(tk + 1) + "/" + String(tiles.size()) + ": " + String(gone) + " of " + String(before)
                        + " spectra folded into a twin");
        }
        if (writer)
        {
          for (auto& sp : tile_out)
          {
            sp.getDataProcessing().push_back(out_dp); // ONE entry, by pointer, for every spectrum
          }
          if (tiles.size() == 1)
          {
            writer->setExpectedSize(tile_out.size(), 0); // exact count: the bulk writer's bytes
          }
          {
            Phase _ph("WRITE(mzML)");
            writer->writeBlock(tile_out);
          }
          n_out += tile_out.size();
          writeLogInfo_("[tile] " + String(tk + 1) + "/" + String(tiles.size()) + ": " + String(tile_out.size()) + " spectra written" + clockStamp()
                        + rss());
        }
        else
        {
          all_out.reserve(all_out.size() + tile_out.size());
          for (auto& sp : tile_out)
          {
            all_out.push_back(std::move(sp));
          }
        }
      }
      tile_progress.setProgress((SignedSize)tk + 1);
#ifdef __GLIBC__
      // [trim] Between tiles, under perf:malloc_trim: nothing else returns the allocator's free pages here -- the in-loop
      // trim above runs only when there is one tile -- so a tiled loop's RSS would ratchet toward the arena pass 1 sized.
      // After tile k is written, before tile k+1 is read. Output-identical.
      // [pipe] Pipelined, the later tiles in flight keep their pages: the trim returns what tile k left behind.
      // [no-trim] Only under memory pressure, when the headroom (MemAvailable, cgroup) is below our own RSS: each trim took
      // 0.4-0.9 s at 64 threads and the next tile faulted the pages back in (64 threads: 2.4 s wall, a third of the page
      // faults; peak RSS +0.5-1 GB). Where memory is short, as on 16-32 GB machines, it still trims every tile.
      if (trim_on && tk + 1 < tiles.size() && availableBytes() < (Size)std::max(0L, rssMb()) * (Size)(1 << 20))
        trimLogged_("tile " + String(tk + 1) + "/" + String(tiles.size()) + " (memory pressure)");
#endif
    }; // finish_tile

    if (! pipe)
    {
      for (Size tk = 0; tk < tiles.size(); ++tk)
      {
        TileRun tile_run;
        tile_run.tk = tk;
        load_tile(tile_run);
#pragma omp parallel num_threads(n_threads)
#pragma omp single
        admit_tile(&tile_run);
        if (worker_err) { std::rethrow_exception(worker_err); }
        finish_tile(tile_run);
      }
    }
    else
    {
      // [pipe] The driver is thread 0 of ONE parallel region over all tiles: it reads tile fin + pipe_depth - 1
      // while tile fin's windows run, admits its windows as the cap and memory allow, and once tile fin's windows
      // are done hands the tile to the WRITER, a task that finishes the tiles strictly in tile order. So the read
      // of the next tile never waits behind a write.
      // The driver does not run windows: it reaches no taskgroup end or taskwait (teamFor, and parTasks keeps its
      // passes inline), so its n_threads - 1 team mates run the windows, the writer, and the reader's and writer's helper
      // tasks when they are free. One exception remains: past libgomp's task cliff (team task count > 64 x nthreads) a
      // spawned window or writer runs undeferred in the driver -- a stall, not a deadlock.
      // Thread 0 rather than `single`'s first arrival: the vendor reader's decode loop opens a nested team next
      // to the driver's place, so under OMP_PROC_BIND=close it shares the window team's cores, never others.
      // Deadlock-free: a window waits only on its own previous tile (CarryGate), which was admitted earlier and
      // is dequeued earlier (the team queue is FIFO), and a tile is read only after tile j - pipe_depth finished,
      // so that previous tile never waits itself; the writer of tile f waits only for the writer of tile f - 1,
      // spawned (and dequeued) before it.
      omp_set_max_active_levels(std::max(2, omp_get_max_active_levels()));
      vector<std::unique_ptr<TileRun>> runs(tiles.size());
      std::exception_ptr drv_err;
      auto drv_fail = [&]() {
#pragma omp critical(pipe_err)
        if (! drv_err) drv_err = std::current_exception();
        pipe_abort.store(true);
      };
      std::atomic<Size> written {0};      // tiles the writer has finished, in tile order
      auto write_tile = [&](const Size f) // the writer: outside the region, so a queued task never outlives it
      {
        while (written.load(std::memory_order_acquire) < f)
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        try
        {
          if (! pipe_abort.load()) { finish_tile(*runs[f]); }
        }
        catch (...)
        {
          drv_fail();
        }
        runs[f].reset();
        written.store(f + 1, std::memory_order_release); // on every exit: the next writer never waits forever
      };
#pragma omp parallel num_threads(n_threads)
      if (omp_get_thread_num() == 0)
      {
        // A team of one (OMP_DYNAMIC, a thread limit) has nobody to run the windows while the driver polls: then
        // one tile at a time, run by the driver's own taskwait and written by the driver -- the barrier's schedule.
        const bool solo = omp_get_num_threads() < 2;
        const Size depth = solo ? 1 : (Size)pipe_depth;
        Size next = 0; // the next tile to read
        auto issue = [&](const Size upto) {
          for (; next < tiles.size() && next < upto && ! pipe_abort.load(); ++next)
          {
            runs[next].reset(new TileRun);
            runs[next]->tk = next;
            load_tile(*runs[next]);
            admit_tile(runs[next].get());
          }
        };
        // [pipe] While it waits for tile R's windows the driver runs one that no thread has started yet, in a task of
        // its own (if(0): undeferred, on this thread) so that the window's taskwaits and taskgroups wait only for the
        // window's own tasks, never for the other windows and writers, which are children of the driver's implicit task.
        // Only windows of the tile it waits for: their carry came from the previous tile, whose windows are all done, so
        // the window never waits for another. Its passes run as tasks again (parTasks), off the driver flag.
        auto run_claimed = [&](TileRun& tile_run) -> bool {
          for (; tile_run.claim_next < tile_run.claim_end; ++tile_run.claim_next)
          {
            WinClaim* c = &claims[tile_run.claim_next];
            if (c->taken.exchange(true, std::memory_order_acq_rel)) { continue; }
            ++tile_run.claim_next;
            auto* rw = &run_window;
            onPipeDriver() = false;
#pragma omp task if (0) default(none) firstprivate(rw, c)
            (*rw)(c->R, c->wi, c->need);
            onPipeDriver() = true;
            return true;
          }
          return false;
        };
        Size help_fin = 0; // the tile the driver waits for: its windows' carry is complete, so it may run them
        if (! solo)
        {
          admit_help = [&]() { return teamAssist(8) || (help_fin < runs.size() && runs[help_fin] && run_claimed(*runs[help_fin])); };
        }
        onPipeDriver() = ! solo; // no taskloop on the driver (parTasks); solo runs its windows itself anyway
        try
        {
          for (Size fin = 0; fin < tiles.size() && ! pipe_abort.load(); ++fin)
          {
            help_fin = fin;
            // Backpressure on the writer: at most `depth` + 1 tiles whose windows are done wait for their write. Without it
            // the driver would keep reading and admitting while the writer lags (slow output storage), and finished spectra
            // and writer tasks spinning on their predecessor would pile up without a bound. It never binds while the writer
            // keeps pace. Deadlock-free: writers 0..fin-1 are spawned, and they wait only on each other, in spawn order.
            while (written.load(std::memory_order_acquire) + depth < fin && ! pipe_abort.load())
            {
              if (! teamAssist(8))
              {
                std::this_thread::sleep_for(std::chrono::milliseconds(1)); // help the writer meanwhile
              }
            }
            if (pipe_abort.load()) { break; }
            issue(fin + depth);
            if (! runs[fin])
            {
              break; // the reader stopped (a failure)
            }
            TileRun& tile_run = *runs[fin];
            if (solo)
            {
#pragma omp taskwait
            }
            while (tile_run.left.load(std::memory_order_acquire) != 0 && ! pipe_abort.load()) // help the reader or writer, else run a window
            {
              if (! teamAssist(8) && (solo || ! run_claimed(tile_run))) { std::this_thread::sleep_for(std::chrono::milliseconds(1)); }
            }
            if (pipe_abort.load()) { break; }
            if (solo) { write_tile(fin); }
            else
            {
              auto* wt = &write_tile;
#pragma omp task default(none) firstprivate(wt, fin)
              (*wt)(fin);
            }
          }
        }
        catch (...)
        {
          drv_fail();
        }
        admit_help = nullptr; // it refers to run_claimed, local to this block
        onPipeDriver() = false;
      } // the region's barrier: every admitted window and every writer has finished
      if (drv_err) { std::rethrow_exception(drv_err); }
      if (worker_err) { std::rethrow_exception(worker_err); }
    }

    win_phase.reset(); // [perf-instr] WINDOW_LOOP ends here: the tile loop's own phases were enclosed in it, the final write is not
    if (rtOffGrid().load() != 0)
    {
      writeLogWarn_("[soa-rt] " + String(rtOffGrid().load())
                    + " traces had a retention time between frames; "
                      "it was moved to the nearest frame time.");
    }
    else
    {
      writeLogInfo_("[soa-rt] every trace RT was exactly a frame time (0 snapped).");
    }
    writeLogInfo_("All windows done (admission high-water " + String(adm_hwm.load()) + " of " + String(n_conc) + " allowed concurrent windows)."
                  + clockStamp() + rss()); // [mem]
    tile_progress.endProgress();
#ifdef __GLIBC__
    if (trim_on) { trimLogged_("window loop -> write"); }
#endif

    if (writer)
    {
      {
        Phase _ph("WRITE(mzML)");
        writer->finish();
        writer.reset();
      } // the count patched, the index and footer written; throws on a failed write
      commitOutput_(out_part, out);
    }
    else
    {
      // the bulk path (mzPeak output, one tile): whole
      {
        Phase _ph("ASSEMBLE(PeakMap)");
        out_exp.setSpectra(std::move(all_out));
      }
      n_out = out_exp.size();
      for (Size i = 0; i < out_exp.size(); ++i)
      {
        out_exp[i].getDataProcessing().push_back(out_dp);
      }
      // to <out>.part like the mzML path (the unknown extension takes the one allowed type), then committed
      {
        Phase _ph("WRITE(mzPeak)");
        FileHandler().storeExperiment(out_part, out_exp, {out_type}, log_type_);
      }
      commitOutput_(out_part, out);
    }
    writeLogInfo_("Wrote " + String(n_out) + " pseudo-MS2 spectra to " + out + " (" + String(FileTypes::typeToName(out_type)) + ")");
    if (zmin > 1)
    {
      writeLogInfo_("charge:min_charge " + String(zmin) + ": " + String(n_below_zmin.load()) + " spectra of lower charge assembled and not written");
    }
    if (twin_im >= 0.0)
    {
      writeLogInfo_("[twins] " + String(n_twins) + " spectra folded into a twin (assembly:twin_im_tolerance " + String(twin_im) + ")");
    }
    if (twin_im >= 0.0 && twin_rt_tol > 0.0) // [twin-rt] the before/after of the near-frame relation
    {
      writeLogInfo_("[twins] assembly:twin_rt_frames: " + String(twin_st[0]) + " twin pairs at one MS1 frame, " + String(twin_st[1])
                    + " same-monoisotope pairs at nearby frames joined, " + String(twin_st[2]) + " left apart at an RT cell edge");
    }
    { // [merge] how much splitting the two merges removed
      auto& m = mergeCount();
      writeLogInfo_("[merge] im-slice merge (trace:im_slice_merge "
                    + (! integer_detector          ? String("n/a, openms detector")
                       : cfg_im_slice_merge >= 0.0 ? String(cfg_im_slice_merge)
                                                   : String("off"))
                    + "): " + String(m[0].load()) + " picked MS2 peaks merged away in " + String(m[1].load()) + " frames; assembly same-bin merge ("
                    + (cfg_merge_same_bin ? "on" : "off") + "): " + String(m[2].load()) + " fragments merged away");
    }
    if (cfg_foreign_w < 1.0)
    { // [foreign] what assembly:foreign_precursor_weight did
      auto& c = foreignCount();
      writeLogInfo_("[foreign] assembly:foreign_precursor_weight " + String(cfg_foreign_w) + " (" + String(cfg_foreign_ppm) + " ppm, 1/K0 "
                    + String(cfg_foreign_dim) + ", min_iso " + String(cfg_foreign_min_iso) + "): of " + String(c[6].load())
                    + " precursor-window assemblies with fragments, " + String(c[0].load()) + " had a foreign candidate (" + String(c[4].load())
                    + " candidates skipped as the precursor's own ion); " + String(c[1].load()) + " foreign envelopes matched; " + String(c[3].load())
                    + " fragments " + (cfg_foreign_w > 0.0 ? "attenuated" : "removed") + " in " + String(c[2].load()) + " assemblies, "
                    + String(c[5].load()) + " of which then fell below assembly:min_fragments");
    }
    if (cfg_fuse_ppm > 0.0)
    { // [fuse] what assembly:fragment_fuse_ppm did (per precursor-window assembly, before min_fragments and the cap)
      auto& c = fuseCount();
      const long long g = c[3].load();
      writeLogInfo_("[fuse] assembly:fragment_fuse_ppm " + String(cfg_fuse_ppm) + " (anchored): " + String(c[1].load()) + " fragments in, "
                    + String(c[2].load()) + " peaks out, in " + String(c[0].load()) + " precursor-window assemblies; " + String(g)
                    + " groups of >= 2 fused, mean size " + String(g ? (double)c[4].load() / (double)g : 0.0, false) + "; " + String(c[5].load())
                    + " assemblies fell below assembly:min_fragments");
    }
    reportPhases(phaseClock());
    return EXECUTION_OK;
  }
};

/// @endcond

#ifdef _WIN32
  #define environ _environ
#else
extern char** environ;
#endif

/**
  @brief Entry point: refuse old environment variable names, set up the environment and run the tool.

  @param[in] argc Argument count
  @param[in] argv Arguments
  @return The tool's exit code (UNKNOWN_ERROR for an exception the tool did not catch)
*/
int main(int argc, const char** argv)
{
  // The tool was SpeXtractor (variables SPEXTRACTOR_*) and then DIAspeXtractor (DIASPEXTRACTOR_*); it is DIAspeXtract
  // now, and every variable is DIASPEXTRACT_* (so are the three the patched libOpenMS reads). An old name would be
  // ignored without a word -- a chain's A/B would measure the default -- so any variable under either old prefix refuses
  // the run before anything else happens, naming the variable and its DIASPEXTRACT_ name. There is no dual read.
  static const struct
  {
    const char* prefix;
    Size len;
  } old_prefixes[] = {{"SPEXTRACTOR_", 12}, {"DIASPEXTRACTOR_", 15}};
  bool old_name = false;
  for (char** e = environ; e != nullptr && *e != nullptr; ++e)
  {
    for (const auto& op : old_prefixes)
    {
      if (std::strncmp(*e, op.prefix, op.len) != 0) { continue; }
      const int n = (int)std::strcspn(*e, "="); // the variable name
      const int m = n - (int)op.len;            // the part after the old prefix
      std::fprintf(stderr,
                   "Error: %.*s is set, but this tool reads only DIASPEXTRACT_* variables: "
                   "rename it to DIASPEXTRACT_%.*s or unset it. Nothing was run.\n",
                   n, *e, m, *e + op.len);
      old_name = true;
    }
    if (std::strncmp(*e, "SPX_NOIMS_", 10) == 0) // the non-IMS path's tuning variables, removed after 1.8.0: fixed values now
    {
      std::fprintf(stderr,
                   "Error: %.*s is set, but the non-IMS path no longer reads SPX_NOIMS_* variables (their values are "
                   "fixed): unset it. Nothing was run.\n",
                   (int)std::strcspn(*e, "="), *e);
      old_name = true;
    }
  }
  if (old_name) { return TOPPBase::ILLEGAL_PARAMETERS; }
  // TOPPBase's update check sends this tool's name and version to an OpenMS REST endpoint and, offline, prints Qt
  // QIODevice errors on the stderr the harness parses (build.yml asserts --help stays clean). Always off: the tool makes
  // no network request of its own (only an s3:// -in is read over the network; SECURITY.md), whatever the environment says.
#ifdef _WIN32
  _putenv_s("OPENMS_DISABLE_UPDATE_CHECK", "ON");
#else
  ::setenv("OPENMS_DISABLE_UPDATE_CHECK", "ON", 1);
#endif
  // [raw] Thermo .raw is read through .NET (openms-thermo-bridge). A release bundle carries its own runtime in dotnet/ next
  // to bin/, so nothing needs installing; DOTNET_ROOT, if the user set it, wins. Invariant globalization: the reader needs
  // no culture data, and without it .NET on Linux would need the system's ICU.
  {
    const std::filesystem::path bundled = std::filesystem::path(File::getExecutablePath()) / ".." / "dotnet";
    auto set_default = [](const char* k, const std::string& v) {
      const char* cur = std::getenv(k);
      if (cur && *cur) { return; }
#ifdef _WIN32
      _putenv_s(k, v.c_str());
#else
      ::setenv(k, v.c_str(), 1);
#endif
    };
    std::error_code ec;
    if (std::filesystem::is_directory(bundled, ec)) { set_default("DOTNET_ROOT", std::filesystem::canonical(bundled, ec).string()); }
    set_default("DOTNET_SYSTEM_GLOBALIZATION_INVARIANT", "1");
  }
  for (int i = 1; i < argc; ++i)
  {
    if (std::strncmp(argv[i], "-noims:", 7) == 0) { g_noims_cli += std::string(argv[i]) + " "; }
  }
  TOPPDIAspeXtract tool;
  // TOPPBase catches OpenMS exceptions only; one from a library (SQLite, opentims) would otherwise abort with no exit code
  try
  {
    return tool.main(argc, argv);
  }
  catch (const std::exception& e)
  {
    std::fprintf(stderr, "Error: %s\n", e.what());
    return TOPPBase::UNKNOWN_ERROR;
  }
}
