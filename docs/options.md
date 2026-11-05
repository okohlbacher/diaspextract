# DIAspeXtract options

This page lists every option of DIAspeXtract with its default. The defaults are the recommended settings, so a
normal run needs only `-in` and `-out`:

```bash
diaspextract -in sample.d -out pseudo.mzML
```

`diaspextract --helphelp` prints the same list; `diaspextract --help` prints only the common options. Give a value
after the option name, for example `-gate:delta_rt 4`. Options with the values `true` and `false` need the word, for
example `-perf:malloc_trim false`. A flag takes no value. `-write_ini run.ini` writes every option to an ini file, and
`-ini run.ini` reads it back.

Options marked "set per run" print the value they chose in the log, and the header of the mzML output records it as
a `spx:` entry. Pass that value to fix it for another run.

How to install, run and read the output, and how to report a problem: [README.md](../README.md).

Input with ion mobility (diaPASEF) and input without it (Orbitrap, SWATH) take separate paths. All sections except
[Input and output](#input-and-output) and [Input without ion mobility](#input-without-ion-mobility) apply to diaPASEF
input only; input without ion mobility reads `-threads` and the options of its own section.

Sections: [Input and output](#input-and-output) · [Speed and memory](#speed-and-memory) ·
[Trace detection](#trace-detection) · [Precursor charge](#precursor-charge) · [Co-elution gate](#co-elution-gate) ·
[Assembly](#assembly) · [MS1 denoising](#ms1-denoising) · [Input without ion mobility](#input-without-ion-mobility) · [Diagnostics](#diagnostics) ·
[Environment variables](#environment-variables) · [Reproducing 1.5.0 spectra](#reproducing-150-spectra)

## Input and output

| option | default | what it does |
|---|---|---|
| `-in` | required | The DIA run: a Bruker `.d` directory (diaPASEF), an mzML file with or without ion mobility, or an mzPeak archive (mzPeak needs a build with mzPeak support). A `.d` is read frame by frame; an mzML is always read whole. Input without ion mobility at MS1 takes [its own path](#input-without-ion-mobility). |
| `-out` | required | The output spectra. The format follows the extension: `.mzML`, or `.mzpeak` in a build with mzPeak support. Search engines read mzML, so write `.mzML` if a search comes next. mzPeak output runs as one tile and needs the whole run's memory. |
| `-out_type` | (empty) | Sets the output format instead of the `-out` extension. Values: empty, `mzML`, `mzpeak`. Empty follows the extension; an unknown extension gives mzPeak (mzML in a build without mzPeak support). |

The standard OpenMS options work as usual: `-threads`, `-ini`, `-write_ini`, `-log`, `-no_progress` and `-debug`.
If `-threads` is not on the command line, the tool uses every physical core (hardware threads of one core, SMT, count
once; a taskset, cpuset or cgroup quota is respected); pass `-threads <n>` on a shared machine. A value above the
physical cores is lowered to them, with a warning: more threads only compete for the cores. A
`threads` value above 1 in an ini file still applies, but `threads 1` in an ini file counts as unset. To run on one
core, pass `-threads 1` on the command line.

## Speed and memory

These options trade run time against memory. They do not change the spectra, except the four marked
**Changes the output**.

| option | default | what it does |
|---|---|---|
| `-tile:rt_sec` | 600 | Width (seconds) of the fixed retention-time cells the integer detector works in. A trace that crosses a cell edge is cut there. 0 or negative: the whole run is one cell and one tile, at several times the memory. Ignored with `-trace:detector openms`. **Changes the output.** |
| `-tile:cells_per_tile` | 1 | Process the run in tiles of this many cells. Fewer cells per tile use less memory. 0: the whole run is one tile, at several times the memory. mzPeak output always runs as one tile. |
| `-perf:malloc_trim` | `true` | Return freed memory to the system between steps and tiles (Linux only). This lowers peak memory at some extra run time. Set `false` when time matters more than memory. Values: `true`, `false`. |
| `-perf:stream_load` | `true` | Read the input frame by frame (`.d`, or mzPeak in a build with mzPeak support; mzML is always read whole). `false` loads the whole run first, which needs much more memory and takes longer. Values: `true`, `false`. **Changes the output.** |
| `-perf:ms1_trace_bands` | 48 | Trace MS1 in this many m/z bands in parallel. 1 turns it off; 0 counts as 1. At least 0. **Changes the output** slightly. |
| `-perf:ms1_prune` | `true` | Drop MS1 peaks at or below `-trace:noise_threshold_int` while loading, which saves memory. `false` keeps every peak. Values: `true`, `false`. |
| `-perf:trace_bands` | 12 | Trace each isolation window's fragments in this many m/z bands in parallel. 1 turns it off. 0 picks the count automatically, which makes the output depend on the thread count; the integer detector refuses 0 when `-tile:rt_sec` is above 0 or a `.d` is streamed, as at the defaults. At least 0. **Changes the output.** |
| `-perf:mem_fraction` | 0.75 | Share of the free memory that the isolation windows may use, checked again before each window starts. Free memory is the kernel's MemAvailable, or the headroom under a cgroup limit (Slurm `--mem`, a container; page cache counts as free) when that is lower; on macOS, free plus inactive memory. One window always runs, however large. 0.05 to 0.95. |
| `-perf:max_concurrent_windows` | 0 | Most isolation windows processed at once. Fewer use less memory but take longer. 0: no cap (free memory still limits them). |

**Short of memory?** Keep `-tile:cells_per_tile 1`, lower `-perf:max_concurrent_windows` or `-perf:mem_fraction`, set
`DIASPEXTRACT_PIPE_TILES=1` (see [Environment variables](#environment-variables)), and write `.mzML`, not `.mzpeak`.
`-trace:detector openms` also needs much more memory. **Short of time?** Raise `-threads`, and set
`-perf:malloc_trim false` if memory is plentiful.

## Trace detection

The tool follows each ion over time (a trace), for precursors in the MS1 frames and for fragments in the MS2 frames.

| option | default | what it does |
|---|---|---|
| `-trace:mass_error_ppm` | 15 | m/z tolerance (ppm) of trace detection. At least 0. |
| `-trace:noise_threshold_int` | 100 | MS1 (precursor) noise level: peaks at or below it are not traced. At least 0. |
| `-trace:ms2_noise_threshold_int` | 10 | MS2 (fragment) noise level: peaks at or below it are not traced. At least 0. |
| `-trace:min_length_sec` | -1 | Shortest MS1 trace kept, in seconds from first to last point. Negative: set per run, `-trace:min_length_fwhm` × the run's median MS1 peak width, at least half an MS1 cycle; 3 s when the run has too few bright traces to measure the width. 0: no filter. 3 was the fixed value before 1.6.0. |
| `-trace:min_length_fwhm` | 0.4 | For the automatic `-trace:min_length_sec`: the shortest MS1 trace as this fraction of the run's median MS1 peak width. At least 0. |
| `-trace:ms1_chrom_peak_snr` | 3 | An MS1 trace's apex must exceed this multiple of `-trace:noise_threshold_int`. At least 0. |
| `-trace:ms2_chrom_peak_snr` | 1 | A fragment trace's apex must exceed this multiple of `-trace:ms2_noise_threshold_int`. Lower it to keep weaker fragments. At least 0. |
| `-trace:ms2_min_length_sec` | 0 | Shortest fragment trace kept (seconds). 0: no filter. |
| `-trace:detector` | `integer` | `integer` traces on the instrument's flight-time bins and uses much less memory. It needs the Bruker m/z calibration; without it (for example with a diaPASEF mzML) the tool falls back to `openms` and warns. `openms` uses OpenMS's mass-trace detection and works on whole isolation windows, without retention-time cells. The two find partly different peptides. Values: `integer`, `openms`. |
| `-trace:ms2_split_valleys` | 7 | Split fragment traces at valleys in their elution profile. The value is the expected peak width (seconds). 0: off. |
| `-trace:ms1_split_valleys` | 7 | Split precursor traces at valleys, so two peptides that elute apart at the same m/z stay separate. The value is the expected peak width (seconds). 0: off. |
| `-trace:max_span_sec` | 120 | Trim each trace to at most this many seconds around its apex. 0: no limit. At least 0. |
| `-trace:native_ms1_neighbors` | 0 | Add the n neighbouring MS1 frames on each side to every MS1 frame before peak picking. Bruker `.d` with `-perf:stream_load true` only, and it needs `-dnoise:ms1 false`. 0: off. |
| `-trace:im_slice_merge` | 0.015 | Integer detector: before tracing, merge a frame's fragment peaks in one flight-time bin whose 1/K0 values lie within this distance, so one ion does not become several traces. Negative: off. At most 1. |

## Precursor charge

The tool finds each precursor's isotope peaks, its charge and its monoisotopic peak. A precursor belongs to the
isolation windows whose bounds, rounded to 0.01 Th, contain its monoisotopic m/z; a window edge off that grid
(499.996) therefore counts as its rounded value (500.00).

| option | default | what it does |
|---|---|---|
| `-max_charge` | 5 | Highest precursor charge considered. At least 1. |
| `-charge:min_charge` | 1 | Lowest precursor charge written. Singly charged precursors are real 1+ ions but carry a higher false discovery rate, so control the FDR per charge after the search. 2 writes only 2+ and higher; the 2+ and higher spectra are the same as at 1 (lower charges are still assembled, for the foreign-isotope removal of the others). 1 to `-max_charge`. |
| `-charge:im_charge_veto` | `true` | Use ion mobility to correct charge calls. A 1+ call on the run's 2+ mobility trend becomes 2+; one on the 3+ trend is dropped. Genuine 1+ ions keep their charge. Values: `true`, `false`. |
| `-charge:im_veto_band` | 2.5 | For `-charge:im_charge_veto`: half-width of the 2+ and 3+ mobility bands, in robust standard deviations of the fit. |
| `-charge:mono_profile_check` | `v2` | Check each precursor's monoisotopic peak against the summed MS1 intensity at its isotope positions, and move it by one isotope when the expected isotope pattern disagrees. Values: `v2` (on), `v3` (opt-in: every candidate monoisotope position of the isotope walk is scored by its cosine with the averagine pattern times the share of that pattern observed; on four diaPASEF runs 0.5–2 % more peptides with Sage and 0.6–1.8 % with MSFragger, spectra change), `false` (off, as in 1.5.0). |
| `-charge:mono_position_min_ratio` | 0.2 | For `v2`: drop the lightest isotope when its intensity, relative to the expected pattern, is below this ratio; the monoisotopic peak then moves one isotope heavier. At least 0. |
| `-charge:mono_position_im_box` | 0 | For `v2`: 1/K0 half-width of the isotope intensity sums. 0: set per run (0.01 to 0.03). At least 0. |
| `-charge:mono_position_extend` | 0.6 | For `v2`: when nothing was dropped, move the monoisotopic peak one isotope lighter if that position's ratio is at least this. 0: off. At least 0. |
| `-charge:mono_position_iterate` | `true` | For `v2`: after a drop, test the new lightest isotope too, and repeat. `false` tests once. Values: `true`, `false`. |
| `-charge:iso_im_tolerance` | 0 | 1/K0 tolerance for matching a precursor's isotope peaks. 0: set per run (0.01 to 0.05). A value above 0 fixes it; 1.5.0 used 0.05. |

## Co-elution gate

The gate decides which fragments join a precursor's spectrum: a fragment must elute with the precursor, at the same
ion mobility.

| option | default | what it does |
|---|---|---|
| `-gate:delta_im` | 0.01 | Largest 1/K0 difference between a precursor and a fragment. It also sets the 1/K0 tolerance of peak picking and tracing, so it changes detection too. To widen only the fragment match, use `-gate:fragment_delta_im`. At least 0. |
| `-gate:delta_rt` | 4 | Largest difference (seconds) between the apex times of a precursor and a fragment. At least 0. 3 before 1.8.0. |
| `-gate:min_correlation` | 0.3 | Least correlation between the elution profiles of a precursor and a fragment. Raise it for cleaner but sparser spectra. -1 to 1. The correlation is a Pearson correlation over every frame of the isolation window in the run, frames without signal counting as zero; for two short elution peaks it is close to the cosine similarity of their profiles. |
| `-gate:min_correlation_points` | 3 | Fewest shared points a correlation needs. Where `-gate:shrinkage` applies, `-gate:shrink_min_points` is used instead. |

**Faint precursors.** A precursor with signal in only a few frames is faint. Its fragment profiles are too short for a
plain correlation, so three extra rules judge them: `-gate:shrinkage`, `-gate:colocation` and `-gate:short_traces`.
`-gate:adaptive` decides which precursors they apply to.

At the defaults the correlation test does not reject a fragment of a faint precursor: with at most 4 shared points
the shrinkage pulls every correlation, even -1, to at least 0.31, above `-gate:min_correlation`. Such a fragment is
kept when it elutes with the precursor (`-gate:faint_delta_rt`, `-gate:delta_im`) and passes `-gate:colocation` or,
when it has too few points for a correlation, `-gate:short_traces`. With `-gate:adaptive off` the same holds for
every precursor, up to 16 shared points for an uncorrelated fragment. A larger `-gate:shrink_prior_sd` (a weaker
pull), a smaller `-gate:shrink_prior_r` or a higher `-gate:min_correlation` makes the test bite; a smaller
`-gate:shrink_prior_sd` pulls harder toward the prior and lets even more fragments pass.

| option | default | what it does |
|---|---|---|
| `-gate:adaptive` | `faint` | `faint`: the extra rules apply to faint precursors only, which also use `-gate:faint_delta_rt` and `-gate:faint_corr_power`. `auto`: as `faint` on runs with many faint precursors (`-gate:auto_faint_fraction`), otherwise no extra rules. `off`: the rules apply to every precursor. Values: `off`, `faint`, `auto`. |
| `-gate:faint_max_frames` | 4 | A precursor with signal in at most this many frames counts as faint. At least 1. |
| `-gate:faint_delta_rt` | 6 | `-gate:delta_rt` for faint precursors (seconds). At least 0. |
| `-gate:faint_corr_power` | 1 | `-assembly:corr_power` for faint precursors. At least 0. |
| `-gate:auto_faint_fraction` | 0.45 | For `-gate:adaptive auto`: use the faint rules when at least this fraction of the run's precursors are faint. 0 to 1. |
| `-gate:shrinkage` | `fisher` | `fisher`: pull the correlation of a fragment with few shared points toward a typical value (`-gate:shrink_prior_r`) before the `-gate:min_correlation` test, instead of requiring `-gate:min_correlation_points` points. `off`: use the plain point count. Values: `off`, `fisher`. |
| `-gate:shrink_prior_r` | 0.6 | For `fisher`: the typical correlation. |
| `-gate:shrink_prior_sd` | 0.3 | For `fisher`: the spread of the typical value (Fisher z units). Smaller pulls harder. At least 0.01. |
| `-gate:shrink_min_points` | 1 | For `fisher`: fewest shared points to compute a correlation. At least 1. |
| `-gate:colocation` | `ranges` | `ranges`: keep a fragment only if its apex lies inside the precursor's retention-time span, the precursor's apex lies inside the fragment's (each with one frame of slack), and the two spans overlap by at least `-gate:colocation_overlap`. `off`: skip this test. Values: `off`, `ranges`. |
| `-gate:colocation_overlap` | 0.3 | For `ranges`: least overlap of the two spans, as a fraction of the shorter one. 0 to 1. |
| `-gate:short_traces` | `concentric` | Fragment traces too short for a correlation. `concentric`: keep one whose retention-time span lies inside the precursor's (one frame of slack) and is centred on it. `off`: drop them. Values: `off`, `concentric`. |
| `-gate:short_center_tol` | 0.5 | For `concentric`: how far the fragment's span centre may lie from the precursor's, as a fraction of the precursor's half-span (at least one cycle). At least 0. |
| `-gate:short_min_points` | 1 | For `concentric`: fewest points a short fragment trace needs. At least 1. |

**Experimental, off by default.**

| option | default | what it does |
|---|---|---|
| `-gate:fragment_delta_im` | -1 | Largest 1/K0 difference for the fragment match only; tracing keeps `-gate:delta_im`. Negative: use `-gate:delta_im`. |

## Assembly

Assembly turns each precursor and its gated fragments into one spectrum.

| option | default | what it does |
|---|---|---|
| `-assembly:min_fragments` | 3 | Write a spectrum only if it has at least this many fragments. At least 1. |
| `-assembly:max_fragments` | 500 | Keep at most this many fragments per spectrum, best first. The twin fold can exceed it: a folded spectrum holds the fragments of both twins. At least 1, and at least `-assembly:min_fragments`. |
| `-assembly:im_weight_sigma` | 0.005 | Weight each fragment's intensity by how close its 1/K0 is to the precursor's (a Gaussian of this width), before the fragment cap. 0: off, as in 1.5.0. |
| `-assembly:corr_power` | 2 | Weight each fragment's intensity by its elution correlation with the precursor, raised to this power. 0: off. |
| `-assembly:fragment_fuse_ppm` | 25 | Fuse fragments within this many ppm into one peak (intensities summed) before the fragment cap. One fragment ion often shows up as several peaks a few ppm apart. 0: off, as in 1.5.0. 0 to 50. |
| `-assembly:merge_same_bin` | `true` | Merge fragments at the same m/z (one flight-time bin) into one peak before the fragment cap. Leave it on: two peaks at one m/z break the fragment matching of OpenMS-based search engines. Values: `true`, `false`. |
| `-assembly:foreign_precursor_weight` | 0 | An isolation window also passes other precursors, and their unfragmented isotope peaks land among the fragments. Multiply such fragments by this weight: 0 removes them, 1 keeps them (as in 1.5.0). The precursor's own isotope peaks are never touched. 0 to 1. |
| `-assembly:foreign_precursor_ppm` | 10 | For `-assembly:foreign_precursor_weight`: m/z tolerance (ppm) between a foreign isotope position and a fragment. At least 0. |
| `-assembly:foreign_precursor_dim` | 0.015 | For `-assembly:foreign_precursor_weight`: largest 1/K0 difference between this precursor and a foreign one. At least 0. |
| `-assembly:foreign_precursor_min_iso` | 2 | For `-assembly:foreign_precursor_weight`: how many of a foreign precursor's four isotope positions (M to M+3) must carry a fragment. 1 to 4. |

**Twin fold.** Sometimes two spectra claim the same precursor. The twin fold merges them into one spectrum that holds
the fragments of both.

| option | default | what it does |
|---|---|---|
| `-assembly:twin_im_tolerance` | 0.005 | Fold two spectra with the same charge, precursor m/z within 20 ppm and 1/K0 within this tolerance. Negative: off. 0 folds only twins with identical 1/K0. At least -1. Twins are joined pair by pair, so one cluster can be a chain whose ends lie further apart than 20 ppm; on real runs such chains are one ion called at drifting m/z. |
| `-assembly:twin_rt_frames` | 0 | Also fold twins whose apexes lie up to n + 0.5 MS1 cycles apart and in the same retention-time cell (a pair across a cell edge stays apart), when their precursor m/z lie within 0.05 ppm (one flight-time bin) and the monoisotope check (`-charge:mono_profile_check v2` or `v3`) moved one of the two monoisotopic peaks one isotope lighter (see `-charge:mono_position_extend`). 0: set per run from the MS1 peak width. -1: the same MS1 frame only (as before 1.6.0). At least -1. |

## MS1 denoising

Removes noise from every MS1 frame before peak picking, with a port of dnoise v0.1.0 (Garrett, Diedrich & Yates III,
bioRxiv 2026.08.27.747603). It applies to Bruker `.d` input only; other inputs pass through unchanged. It also lowers
peak memory. Denoising needs the run's own m/z calibration table. The tool stops with a message when that table is
not usable, when `OPENMS_BRUKER_SDK_PATH` is set, or with `-trace:native_ms1_neighbors` above 0; add
`-dnoise:ms1 false` in those cases.

The port filters the raw points dnoise filters, recovered from what the `.d` reader delivers. Up to 100 ms
accumulation time every raw count is recovered exactly and the result is dnoise's, point for point. Above 100 ms the
reader's intensity correction maps several raw counts to one value; the port then takes the nearest count, so a
filter decision can differ from dnoise's. The log and the `spx:dnoise_ms1_points` header entry count these points
("approximate").

The filter settings are dnoise's defaults and are fixed.

| option | default | what it does |
|---|---|---|
| `-dnoise:ms1` | `true` | Denoise the MS1 frames. Values: `true`, `false`. |

## Input without ion mobility

Orbitrap or SWATH DIA in mzML, centroided, without ion mobility goes to a separate path (precursors: Biosaur2
isotope-pattern features of the MS1 scans; fragments: the mass traces of each isolation window that correlate with the
precursor's elution profile, weighted by r²⁴, at most 200 per spectrum; see the
[README](../README.md#input-without-ion-mobility)). The diaPASEF options do not apply to it. Profile spectra, ion
mobility at MS2 without it at MS1, and MS2 scans with several isolation windows (multiplexed DIA) are refused with
exit code 11.

| Option | Default | What it does |
|---|---|---|
| `-noims:min_corr` | 0.3 | Smallest Pearson correlation between a fragment's and the precursor's profile on the window's scan grid. 0 to 1. |
| `-noims:ppm` | 10 | m/z tolerance of the mass traces and isotope patterns (ppm). At least 0. |

## Diagnostics

| option | default | what it does |
|---|---|---|
| `-diag:dump_ms1_tsv` | (empty) | Write the MS1 traces and the precursors found in them to `<prefix>.traces.tsv` and `<prefix>.precursors.tsv`, to see where a precursor was lost. Empty: off. |
| `-diag:selftest_twins` | off (flag) | Run the built-in twin-fold self-test and exit. |
| `-diag:selftest_arena` | off (flag) | Run the built-in MS1 storage self-test and exit. |
| `-diag:spectrum_meta` | off (flag) | Write three userParams per spectrum: `spx_p_frames` (frames with precursor signal), `spx_faint` (1 when that is at most `-gate:faint_max_frames`), `spx_p_apex` (the precursor's apex intensity). The spectra do not change. |

## Environment variables

A normal run needs none of these.

| variable | default | what it does |
|---|---|---|
| `DIASPEXTRACT_PIPE_TILES` | 2 | How many tiles may have isolation windows running at once (1 to 16). 1 processes one tile at a time: lower peak memory, longer run. The output does not change. |
| `DIASPEXTRACT_LOAD_BATCH` | 256 | Frames the `.d` reader decodes per batch. The output does not change. |
| `OPENMS_BRUKER_SDK_PATH` | unset | Path to Bruker's `libtimsdata` library (not shipped with the tool). The `.d` reader then converts flight times to m/z with Bruker's library, as a cross-check. Needs `-dnoise:ms1 false`. |
| `DIASPEXTRACT_SDK_PARALLEL` | unset | Any value: decode frames in parallel even when Bruker's library does the m/z conversion. Without it, that conversion runs on one thread, because the library is not known to be thread-safe. |
| `DIASPEXTRACT_ALLOW_CHORD_FALLBACK` | unset | `1`: when the exact m/z calibration in `analysis.tdf` cannot be used, accept a less accurate approximation instead of stopping. Needs `-dnoise:ms1 false`. |
| `DIASPEXTRACT_MZPEAK_TDF` | unset | Path to the `analysis.tdf` (or, for an mzPeak archive, also `analysis.tdf.gz`) of the run, for an input that carries no calibration of its own: an mzPeak archive converted without the vendor files, or an mzML converted from a `.d` with `frame=` native IDs. The tool takes the exact m/z calibration from it, which the integer detector needs. |
| `DIASPEXTRACT_MZPEAK_EXACT` | unset | mzPeak builds: `0` always uses the archive's approximate m/z conversion instead of the exact calibration. The integer detector cannot run on it, so add `-trace:detector openms`. Unset, an archive whose exact calibration cannot be recovered stops the run, and the message names both ways out. |

`DIASPEXTRACT_LOAD_BATCH`, `OPENMS_BRUKER_SDK_PATH`, `DIASPEXTRACT_SDK_PARALLEL` and
`DIASPEXTRACT_ALLOW_CHORD_FALLBACK` are read by the patched OpenMS library, so they act only on `.d` input. The tool
refuses to start while a variable from an earlier release is set (`SPEXTRACTOR_*` or `DIASPEXTRACTOR_*`), and names
its new spelling.

**Test switches.** The end-to-end suite uses `DIASPEXTRACT_TILE_NO_CARRY`, `DIASPEXTRACT_MS1_PRUNE_NO_WITNESS`,
`DIASPEXTRACT_MS1_PRUNE_NO_CHAIN`, `DIASPEXTRACT_DET` and `DIASPEXTRACT_TASK_THROW` to break the tool on purpose. Each
changes the output or ends the run; a run with one set warns and records it as the `spx:test_env` header entry.
`DIASPEXTRACT_TEST_D` (a real diaPASEF `.d`) is read by the suite only.

## Reproducing 1.5.0 spectra

Version 1.6.0 changed seven defaults and 1.8.0 one more. Seven options set them back; the eighth needs no reset (see
below, and the v1.6.0 and v1.8.0 entries in [CHANGELOG.md](../CHANGELOG.md)):

```bash
diaspextract -in sample.d -out pseudo.mzML \
  -assembly:im_weight_sigma 0 -assembly:foreign_precursor_weight 1 -assembly:fragment_fuse_ppm 0 \
  -charge:mono_profile_check false -charge:iso_im_tolerance 0.05 -trace:min_length_sec 3 -gate:delta_rt 3
```

They restore 1.5.0's settings, not its exact bytes: 1.8.0's platform-independent math and tie-breaking change the
spectra slightly. Exact 1.5.0 output needs the 1.5.0 binary.

The eighth, `-assembly:twin_rt_frames`, needs no reset: with `-charge:mono_profile_check false` no twins fold across
MS1 frames.
