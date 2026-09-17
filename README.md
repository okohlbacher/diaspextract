# DIAspeXtract

Extract pseudo-DDA ("pseudo-MS/MS") spectra from Bruker timsTOF **diaPASEF** data.

DIAspeXtract was DIAspeXtractor (executable `diaspextractor`) in 1.2.0 and 1.2.1, and the tool released as
speXtract (executable `spextract`) up to v1.1.0; it was briefly SpeXtractor in development. speXtract and
SpeXtractor were taken by other tools; 1.3.0 only drops the "-or". DIAspeXtract 1.3.0 continues that version
line, and [CHANGELOG.md](CHANGELOG.md) maps the old names and defaults to the new ones.

DIAspeXtract turns a data-independent acquisition into a set of pseudo-tandem spectra that any
ordinary DDA search engine can read. Each emitted spectrum pairs one precursor hypothesis with the
fragment traces that co-elute with it in both retention time and ion mobility. The output is mzML
(or mzPeak, following the `-out` extension), so it feeds Sage, MSFragger, Comet or anything else that
reads mzML — no library, no spectral prediction, no enumeration of the search space in advance.

That last point is the reason the tool exists: because the search space is not enumerated before
the search runs, an **open or blind search** over the output can report variants and unexpected
modifications that a library-based DIA workflow cannot represent.

DIAspeXtract is BSD-3-Clause and is a standalone application. **OpenMS is a prerequisite, not a
host** — the tool links against an OpenMS installation but lives outside the OpenMS source tree.

## Requirements

| | |
|---|---|
| OpenMS | the development commit pinned in `patches/openms.lock` (bd4b895, 2026-04-24; it identifies itself as a 3.6.0 pre-release), patched by `scripts/apply_openms_patches.sh` and built `WITH_OPENTIMS` to read Bruker `.d`. No released OpenMS package works |
| Compiler | C++20 — GCC 13+ or Clang 16+ |
| OpenMP | required |
| CMake | 3.21+ (OpenMS's own requirement; the tool itself needs 3.16) |
| **Memory** | **~22 GB peak** on a 30-minute gradient (dataset D) and **~21 GB** on a 2-hour acquisition (TNBC 009) at the shipped defaults, 100 threads, with one 600-s retention-time cell per tile (25 GB on the 2-hour file with `-perf:malloc_trim false`); in one tile (`-tile:cells_per_tile 0`) the same files peak at ~37 and ~100 GB (the 30–60 min cohort spanned 68–109 GB at the 2026-09-04 defaults). It follows how dense the acquisition is more than how long it runs |
| Disk | measured on the release build at the shipped defaults: a 30-minute run (dataset D) writes **496,722** spectra as a **4.9 GB** mzML, a 2-hour run (TNBC 009) **3,010,923** as **27.8 GB** — about 9–10 kB per spectrum on both. That is about a quarter *fewer* spectra than the pre-1.2.0 default emitted (656k and ~3.9 M), but size from these figures, not from an older run's: per-spectrum size differs between acquisitions. mzPeak is about a third of the mzML (dataset D's 496,722 spectra are 1.76 GB) |

The memory figure is not a suggestion. A full-depth diaPASEF acquisition holds millions of mass traces
in flight. At the shipped defaults the window loop holds one 600-s tile at a time, so a 30-minute in-house
acquisition needs ~22 GB and a 2-hour one ~21 GB (25 GB with `-perf:malloc_trim false`); in one tile they
needed ~37 and ~100 GB, and low-load public files (PXD017703) ran in 10–25 GB under the earlier one-tile
defaults. Peak usage is reported at the end of every run, so you can
size a node from your own data. Memory figures throughout are peak RSS as `/usr/bin/time` reports
them (kB ÷ 10⁶).

## Install

```bash
git clone https://github.com/okohlbacher/diaspextract.git
cd diaspextract

# OpenMS first: the pinned development commit, patched, built. No released package works.
git clone https://github.com/OpenMS/OpenMS.git /path/to/OpenMS
git -C /path/to/OpenMS checkout $(sed -n 's/^OPENMS_BASE=//p' patches/openms.lock)
cmake -S /path/to/OpenMS -B /path/to/OpenMS/build -DCMAKE_BUILD_TYPE=Release \
      -DWITH_OPENTIMS=ON -DWITH_GUI=OFF -DHAS_XSERVER=OFF   # plus OpenMS's own dependencies

scripts/apply_openms_patches.sh /path/to/OpenMS    # see "OpenMS patches" below
cmake --build /path/to/OpenMS/build -j             # rebuild OpenMS: the patches change its headers and libOpenMS
cmake -B build -DOpenMS_DIR=/path/to/OpenMS/build
cmake --build build -j
export OPENMS_DATA_PATH=/path/to/OpenMS/share/OpenMS  # a source-built OpenMS has no installed share/; the binary dies on --help without it
```

The binary is `build/diaspextract`. `cmake --install build` puts it in `<install prefix>/bin` (set `-DCMAKE_INSTALL_PREFIX`).

### OpenMS patches

`scripts/apply_openms_patches.sh` installs `src/TdfMzCalibration.h` into the OpenMS tree and applies
four patches. **All four matter. The build stops without the first and the third; it does not
notice a missing second or fourth:**

- **Vendor m/z calibration.** Without it the Bruker reader converts flight time to m/z with a
  two-point linear-in-sqrt chord that is **−5 to −11 ppm biased** (m/z dependent) on every file
  measured, worth roughly **6–11% more closed-search Sage identifications** on the three files measured. Every emitted mzML
  records which calibration was used in the `spx:mz_calibration` userParam
  (`tdf_table_modeltype1` | `bruker_sdk` | `legacy_chord_APPROXIMATE`). The exact model covers
  `ModelType` 1 with the digitizer-temperature term, the `C4` mass offset and a quadratic term `C2` of
  either sign. A negative `C2`, as timsTOF Pro 2 acquisition software 2.0.53 writes it, is accepted since
  1.3.1, but only inside a checked domain (`src/TdfMzCalibration.h`). The golden tests pin each variant
  (`tests/calibration_golden.json`, `tests/test_calibration_cpp.cpp`). An unsupported calibration
  table **fails closed** rather than silently producing biased masses.
- **Lock-free elution-peak detection.** OpenMS guards a shared vector with a program-global critical
  section. Called from inside DIAspeXtract's parallel window loop, that one lock serialises the tool.
- **`MassTrace` move operations.** OpenMS declares a defaulted destructor and copy operations on
  `MassTrace`, which suppresses the implicit moves, so every `std::move` deep-copied its points. A
  `static_assert` stops the build without this one.
- **Parallel mzML serialisation.** Writing the spectrum list is CPU-bound, not I/O-bound — measured,
  it sustains ~400 MB/s against a filesystem that does 4.0 GB/s — and the stock writer does it on
  one thread. Encoding one spectrum is independent of every other, so this encodes chunks in
  parallel and appends them in index order, giving byte-identical output. Each thread needs its own
  handler and validator: the writer memoises CV-term validation in mutable state, and sharing one
  handler corrupts the heap. Measured on a 30-minute acquisition: **20.1 s → 5.9 s**.

Three environment variables reach the patched libOpenMS, none of them needed for a normal run:
`DIASPEXTRACT_ALLOW_CHORD_FALLBACK=1` (exactly `1`) opts into the legacy chord when the exact calibration
is unavailable, which otherwise fails closed; `DIASPEXTRACT_LOAD_BATCH=<n>` sets the frames per decode
batch of the parallel `.d` loader (default 256; the frames reach the detector in the same order at any
value); and `DIASPEXTRACT_SDK_PARALLEL` (any value) lets that loader decode in parallel while the Bruker
SDK converter is active, which it does serially by default because the SDK's thread-safety is not vouched
for.

Optionally, `OPENMS_BRUKER_SDK_PATH=/path/to/libtimsdata.so` uses Bruker's own library for the
conversion instead. It is an independent cross-check, not a requirement, and is not redistributed.

### Optional: `.mzpeak` input

Build the mzPeak C++ library once against the same Arrow/Parquet/libzip that OpenMS uses
(`scripts/build_mzpeak_lib.sh`), then configure with `-DMZPEAK_ROOT=<checkout>`. This also needs the
OpenMS mzPeak integration (the OpenMS/mzpeak fork), which the pinned commit does not carry; the plain recipe above builds without it,
reads `.d` and mzML, and writes mzML. `.d` remains the primary input and is faster.

## Run

```bash
diaspextract -in sample.d -out pseudo.mzML -threads 64
```

That is the complete command. Every default is the configuration the project benchmarks; a run that
needed extra flags to reproduce a published figure would be a bug. The end-to-end suite pins the
shipped detector and the charge floor (it passes a few extraction options of its own, so it is not a
test of the whole default set).

Then search the output like any DDA file:

```bash
sage sage.json -o results pseudo.mzML
```

### Options worth knowing

| option | default | what it does |
|---|---|---|
| `-threads` | all cores | worker threads; an explicit `-threads 1` on the command line is honoured (only the command line is inspected: `threads = 1` in an ini file is treated as unset). Parallelism is isolation windows × flight-time bands, and the window loop reaches ~78× (30-min) to ~89× (2-hour) thread occupancy at 100 threads; window concurrency is bounded by a free-RAM admission gate, not by core count |
| `-out_type` | follows the extension | `mzML` or `mzpeak`; in a build with mzPeak, `.mzpeak` when the extension says nothing (a plain build writes mzML only). Search engines read mzML, so write `.mzML` whenever a search comes next. mzPeak is written in one piece, so it runs as one tile (the whole run's memory, with a warning) |
| `-trace:detector` | `integer` | `integer` works on the instrument's flight-time bin; `openms` uses OpenMS `MassTraceDetection`. Different algorithms — see below |
| `-charge:min_charge` | 1 | lowest precursor charge to emit. Singly charged precursors are genuine on tryptic data (their mobility sits on the 1+ trend line) and are 14–20% of Sage peptides on a 2-hour acquisition, though only ~1.7% on a 30-minute one. Their stratum carries 3–9% entrapment FDR under a pooled 1% cut (file-dependent), so control FDR per charge downstream. `-charge:min_charge 2 -charge:im_charge_veto false` restores the old behaviour (the veto runs first, so `2` alone keeps its re-called 2+) |
| `-charge:im_charge_veto` | `true` | the 2+ and 3+ ion-mobility trend lines are fitted per run from the run's own confident precursors; a z=1 call on the 2+ line (11–26% of z=1 calls depending on the file: a halved 2+ envelope) is re-called 2+, one on the 3+ line is dropped. Genuine 1+ ions are untouched |
| `-charge:im_veto_band` | 2.5 | half-width of the veto's mobility band, in MAD-sigma of the fitted residuals |
| `-assembly:require_isotope_support` | `true` | drop precursor hypotheses with no isotope partner. `false` roughly doubles emission, was ~7× slower (measured 2026-09-04, before the scheduling fixes) and identifies fewer peptides |
| `-perf:malloc_trim` | `true` | return the allocator's free pages (glibc) at the two phase boundaries and between tiles. At the phase boundaries it trades **-30% peak memory for +8.5% wall** on a 2-hour acquisition in one tile (the loop re-faults the returned pages cold); between tiles it takes about 4 GB off TNBC 009 and 0.5 GB off dataset D at the defaults. Byte-identical output either way. Off when time, not memory, is the constraint |
| `-perf:stream_load` | `true` | read the `.d` frame by frame. `false` holds the whole run in memory (90 GB floor) **and changes the output**: measured, it is slightly *more* sensitive on a large file (+2.0% Sage peptides at unchanged entrapment FDR) for 1.75× the memory and ~1.7× the wall |
| `-perf:ms1_trace_bands` | 48 | band-parallel MS1 mass-trace detection (a halo partition; the band count moves the spectrum list slightly). Raised from 12 on 2026-09-07: MS1 tracing 2× faster; peptide counts within 0.01–2% on both engines and both benchmark acquisitions; 1 = off |
| `-dnoise:ms1` | `true` | **MS1 denoising, a bit-identical port of dnoise v0.1.0** (Garrett, Diedrich & Yates III, bioRxiv 2026.08.27.747603; MIT): per MS1 frame, an ion-mobility streak filter, a halo filter and the diaPASEF isolation-window gate, applied to the raw points before picking. It keeps 12.7–20.5% of the MS1 points. **Changes output**: the result is digest-identical to running the tool on dnoise's own filtered `.d`, where dataset D lost 4.3% (Sage) / 4.0% (MSFragger) of its peptides and TNBC 009 gained 1.5% / lost 0.5%. Peak memory at one cell per tile falls 13–21%. Bruker `.d` input only; other inputs pass through, and the header records it. `false` restores the previous output. The `dnoise:` tuning options mirror dnoise's own |
| `-perf:ms1_prune` | `true` | drop picked MS1 centroids at or below `trace:noise_threshold_int` at load, keeping a chain of witnesses so every MS1 band spectrum and the tracer's input stay exact; byte-identical output. Peak at one cell per tile, dnoise off: dataset D 31.99 → 25.90 GB, TNBC 009 69.11 → 31.97 GB |
| `-trace:max_span_sec` | 120 | trim a mass trace to this many seconds around its apex |
| `-trace:band_edges` | `acquisition` | the integer detector's per-window band edges come from the tdf's metadata (the digitizer's last bin) instead of the window's resident peaks, so a streaming reader can reproduce them. Measured neutral on both benchmark files (Sage +1/+1, MSFragger 0/−3, entrapment +0.00/+0.01 points); `slab` reproduces the pre‑2026‑09‑10 output |
| `-tile:cells_per_tile` | 1 | group the cells into tiles of this many cells and run the window loop tile by tile (each tile's spectra sorted and written before the next starts). Output-identical for any value (the cells are the unit; boundary fragments are carried between tiles); 0 = one tile. On a `.d` (the streaming source) **memory is one tile's**: at the defaults TNBC 009 peaks at 21.07 GB in 13 tiles with a trim between tiles against 99.70 GB in one (wall 9:57 against 9:42, same node), dataset D at 21.99 GB in 3 tiles against 37.44 GB (3:47 against 3:15). mzPeak output always runs as one tile |
| `-tile:rt_sec` | 600 | the integer detector traces each isolation window in fixed retention-time **cells** of this pitch, cut at the run's MS1 frame times and recorded in the mzML header (`spx:tile_boundaries`). The cut is part of the definition: any grouping of cells into tiles, resident or streamed, reproduces the same spectra bit for bit. Science price measured 2026-09-10 against whole-window tracing (`-1`): Sage +0.8% on both benchmark files, MSFragger +1.5% / flat, entrapment FDR +0.14 / +0.16 points (inside the release rule); peptides within 15 s of a cut line — 5% of the set on TNBC 009 — are lost at 9.3–9.4% against 6.0–7.9% elsewhere |

`diaspextract --helphelp` lists every option, most with a pointer to the measurement behind the default.

### Reading the output

Emitted spectra are MS2 with a synthetic precursor. Provenance is recorded as userParams on the run,
and the header is authoritative. A run at the defaults stamps sixteen: the detector and
calibration (`spx:detector`, `spx:mz_calibration`, `spx:require_isotope_support`, `spx:corr_power`,
`spx:pearson_G`, `spx:im_weight_sigma`), MS1 denoising (`spx:dnoise_ms1`, `spx:dnoise_ms1_params`,
`spx:dnoise_ms1_points`), the cell grid and tiling (`spx:tile_rt_sec`, `spx:tile_cells_per_tile`,
`spx:tiles`, `spx:tile_boundaries`, `spx:tile_source`), the band edges (`spx:band_edges`) and the
frozen frame table (`spx:frame_table`) — so a file can always be attributed to the configuration
that produced it.

## Two detectors

`-trace:detector` selects between genuinely different algorithms, not two implementations of one.

- **`integer`** (default) finds a peak's candidate traces by arithmetic on the instrument's own
  flight-time index, and never converts its compact store back to double m/z. It needs the vendor
  calibration and falls back to `openms` — loudly — without it.
- **`openms`** runs OpenMS `MassTraceDetection` on a materialised peak map.

They agree on about 85% of the union of identified peptides. Which one identifies more depends on
the search engine and on the file, so **the reason `integer` is the default is memory** — roughly
40% less — not a peptide gain. If you are chasing identifications on a particular dataset, it is
worth trying both.

## Performance

Measured at the shipped defaults (MS1 denoising and the MS1 prune on, one 600-s cell per tile with a trim between
tiles), 100 threads, on cluster nodes shared with other users, so the walls are indicative:

| | 30-minute gradient (dataset D) | 2-hour acquisition (PXD047793 run 009) |
|---|---|---|
| peak RSS | 22.0 GB (3 tiles) | 21.1 GB (13 tiles) |
| wall time | 3:47 | 9:57 |
| in one tile (`-tile:cells_per_tile 0`), same node | 37.4 GB, 3:15 | 99.7 GB, 9:42 |

Before MS1 denoising, the MS1 prune and the tiled default, a 30-minute gradient (33,553 frames) took 2:54 at 57 GB
and the 2-hour file 11:18 at 187 GB, emitting 0.85 M and 3.9 M pseudo-spectra.
Those walls were **11.6% and 17.7% faster than the build before them** respectively, from four
byte-identical changes in this release (see the CHANGELOG's "Measured"): the Savitzky-Golay
coefficients are cached, mzML spectra are serialised in parallel, the canonical sort sorts keys rather
than 750-byte spectrum objects, and the assembled PeakMap takes its vector whole. Every one reproduces
the pinned digests on both files, so identifications are unchanged.

On the 2-hour file, on the same node, the 2026-09-08 build was **0.61× the wall time and 0.55× the
peak memory** of the reference implementation, with 1.28× its Sage peptides at parity (0.996) with
MSFragger under per-charge FDR control (docs/REFERENCE-COMPARISON-2026-09-08.md); the changes above
take the wall ratio to roughly **0.54×** on the same measurement, and leave every peptide count
untouched. That comparison predates MS1 denoising and the tiled default, which move memory far
further down and cost peptides on one file and gain on the other -- see the CHANGELOG.

The allocator's retained free pages were the largest single memory item on a 2-hour acquisition;
`perf:malloc_trim` returns them at the phase boundaries (−30% peak for +8.5% wall, measured in one tile) and
between tiles (about 4 GB on the 2-hour file), and is on by default. Preloading `tcmalloc_minimal` measured −7% peak RSS (wall-neutral, byte-identical) before the
trim shipped; it has not been re-measured with it, so the two are not known to be additive.

Runtime is dominated by the window loop, which is parallel across isolation windows and, within a
window, across flight-time bands. The loader is 10–14%. Memory, not CPU, sets how many
windows can be in flight: an admission gate re-decides concurrency from free RAM at every window admission.

## Determinism

The spectrum data is **byte-identical across thread counts** for a fixed binary: on real data, dataset D
at the shipped defaults digests `e43672a0` at both 8 and 100 threads (2026-09-10), and at a 300-s cell
pitch the same spectra come out at 1, 2, 3 and 6 tiles at both thread counts; run 009 at 600 s matches
its one-tile run in 13 tiles. The end-to-end suite checks 1 vs 4 threads on synthetic input, which
exercises the fallback detector. The loader batch size is **output-neutral**: batches 256 and 4096
produce the identical digest once the picker's tail is flushed (the 2026-09-09 fix in this release),
so it is a resource knob, not an output knob. The pre-fix verdict — that 4096 changed the digest and
lost 7% of peptides — was the dropped tail, not the batch size; the 1024 arm has not been re-taken
since.
`bench/semantic_digest.py` hashes `<spectrumList>`..`</spectrumList>`, excluding the wall-clock
stamp, the recorded parameters and the trailing byte-offset index — none of which can match across
invocations.

This does **not** survive a code change: a last-ulp difference in one arithmetic path cascades to
~2% in peptide counts. Compare two builds by the accepted peptide **set** and the ppm median, never
by raw counts. The two statements are about different things and both have been measured.

## Tests

```bash
python3 test/test_diaspextract.py build/diaspextract
```

Twenty-five end-to-end checks against a synthetic acquisition (twenty-six in a build that writes mzPeak) cover:
- isotope ownership and the charge floor;
- parameter validation;
- thread- and tile-count invariance;
- mass-trace extension across missing cycles;
- that the shipped detector is the one that actually runs;
- that the MS1 prune changes no output, with negative controls that must change it;
- that a non-finite MS1 m/z is refused before the prune;
- that `dnoise:ms1` refuses a Bruker SDK calibration before the load;
- that a leftover `SPEXTRACTOR_*` or `DIASPEXTRACTOR_*` variable from before the renames refuses the run.

No fixtures, no framework, no network. The suite's only Bruker `.d` holds an `analysis.tdf` and no
frame data, so it cannot run the in-tool dnoise filter. That filter is covered by
`tests/test_dnoise_ms1.cpp` (gate boxes, filter rules and the per-point recovery), by
`tests/test_tdf_load.cpp` (the tdf reader on corrupt, incomplete and out-of-range tables; built only
where SQLite3 is found) and by the cluster digest gates.

`tests/` additionally holds the calibration, dnoise and tdf-reader unit tests and the entrapment-FDR
scoring used for the benchmark record.

## Layout

```
src/        the tool (one translation unit) plus four headers: the calibration model, the tdf reader,
            the MS1 denoiser and the mzPeak streaming loader
patches/    the four OpenMS patches and the pinned OpenMS commit (openms.lock)
scripts/    build helpers, OpenMS patch application, the shipped-option-token check
test/       the end-to-end suite
tests/      calibration and entrapment tests, golden calibration values
bench/      benchmark drivers and the semantic digest
docs/       measurements and design notes
```

## Documentation

- [docs/BASELINE.md](docs/BASELINE.md) — the decision record: every default, what was
  measured to justify it, and what was falsified
- [docs/MZ-AXIS-DESIGN.md](docs/MZ-AXIS-DESIGN.md) — the flight-time index and the integer detector
- [docs/charge-inference.md](docs/charge-inference.md) — charge assignment, its coupling to MS1
  splitting, and the levers that did not work
- [CHANGELOG.md](CHANGELOG.md)

## Citing

See [CITATION.cff](CITATION.cff).

## License

BSD-3-Clause. Derived from OpenMS (BSD-3-Clause) — see [NOTICE](NOTICE).
