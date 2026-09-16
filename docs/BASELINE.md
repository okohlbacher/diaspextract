# dataset D is the decision file

> **Naming.** The tool was renamed DIAspeXtractor on 2026-09-11 and DIAspeXtract on 2026-09-15, and this record calls
> it DIAspeXtract in prose throughout. Commands, variables and paths quoted from earlier runs keep the names they ran
> under: SpeXtractor (speXtract before 2026-09-06, DIAspeXtractor from 2026-09-11 to 1.2.1), binary `spextractor`
> (earlier `spextract`, then `diaspextractor`), and `SPEXTRACTOR_*` / `DIASPEXTRACTOR_*` variables, now
> `DIASPEXTRACT_*`. CHANGELOG.md maps the old names to the new ones.
>
> **Not all of this is reproducible on the shipped binary.** Many measurements below were driven by
> `SPEXTRACTOR_*` probe environment variables (`FRAG_DUMP`, `WINDOW_FILE`, `PRECURSOR_LIST`,
> `MS2_RT_RANGE`, `PICK_SERIAL`, `PICK_STATIC`, `PICK_MZ_MODE`, `Z1_CHARGE_UNSET`, `TILE_TRIM`) and by
> options that 1.2.0 removed. The probes are gone, and the `SPEXTRACTOR_*` prefix itself is now
> **refused**: a leftover variable stops the run rather than being ignored. Those entries stand as the
> dated record of what was measured and decided, not as recipes to re-run.

**All progress/regression decisions are made on dataset D alone.** Other samples are for confirmation
once something has already won here — cross-sample comparison has produced three wrong
conclusions in this project, and single-sample iteration removes that class of error entirely.

## Reference numbers (frozen 2026-07-22)

| measure | value | how |
|---|---|---|
| **peptides @ peptide_q ≤ 0.01** | **10,072** | `sage_deiso.json`, current defaults |
| peptides, common-FDR procedure | 9,977 | DEPRECATED 2026-09-01 (procedure's primary outputs no longer exist) |
| the reference implementation, same procedure | 11,552 | DEPRECATED (doc-only; canonical ref = 11,517) — legacy ratio 86.4% not citable |
| MSFragger, same file | 10,470 | DEPRECATED (the reference implementation 13,014 doc-only; canonical = 13,932 form-level / 12,752 seq-level) |
| PSMs/peptide | 8.98 | the reference implementation 1.66 |
| within-cycle multiplicity | 2.85 | the defect |
| across-cycle spread | 2.66 | **correct** — MS1 FWHM is 2.61 cycles |
| MS1 FWHM | 3.61 s | flat with RT (1.04× across gradient) |

**Current defaults** (registered defaults, verified against the code 2026-09-03 --
`ms2_noise_threshold_int` and `ms2_split_valleys` had been benchmarked at 10 and 7.0 while the
registered defaults were still 100 and 0, so a default run could not reproduce any figure in this
file; the defaults were moved to the benchmarked values): `trace:ms1_split_valleys 7.0`,
`trace:ms2_split_valleys 7.0`, `trace:ms2_noise_threshold_int 10`, `trace:ms2_chrom_peak_snr 1.0`,
`trace:ms2_min_length_sec 0`, `trace:mz_estimator apex`, `charge:scoring count`,
`charge:min_charge 2`, `assembly:corr_power 2`, `perf:trace_bands 12`,
`assembly:require_isotope_support true`, `perf:stream_load true`, `trace:detector integer`.

**The 2026-09-03 defaults sweep only checked VALUE-BEARING options, not FLAGS, and two flags that
every benchmark passes were left defaulting OFF** -- found 2026-09-04 and fixed the same day. Both
are now `<true/false>` options defaulting to `true`, because `registerFlag_` cannot default to true.

`assembly:require_isotope_support` measured on dataset D for the first time in either direction:

| | spectra | wall | window-loop occupancy | Sage | MSFragger |
|---|---|---|---|---|---|
| **true** (every benchmark on record) | 655,776 | **5:15** | 65.7x | **12,482** | **12,337** |
| false (what the shipped default did) | 1,255,577 | **36:25** | **6.0x** | 12,212 | 12,148 |

**The shipped default was 6.9x slower AND worse on both engines.** The only prior measurement
(2026-07-26) put the flag's cost at −0.78% peptides; today it GAINS 2.2% Sage and 1.6% MSFragger.
The sign flipped somewhere between the charge gate, the apex estimator and the integer detector,
which is why it was re-measured rather than inherited. Note the occupancy collapse to 6.0x with the
gate off: 2.15 M precursor hypotheses instead of 1.10 M puts the window loop at 91.8% of the run
and it barely parallelises -- so the guessed singletons cost far more than their share of the work.

`perf:stream_load` is the frame-by-frame reader; with it off, `loadExperiment()` holds the whole run
at a MEASURED 90.3 GB floor at t=0, which alone exceeds the tool's entire current peak (80 GB).
Every figure in this file was taken with it on.

**And it is NOT output-neutral, contrary to the assumption every benchmark has rested on.** Measured
2026-09-04 on dataset D, both arms on one node:

| `perf:stream_load` | spectra | wall | peak RSS | digest |
|---|---|---|---|---|
| true (shipped, and what every figure here used) | 655,776 | 5:11 | 78.4 GB | -- |
| false | 656,371 | 10:30 | 121.8 GB | **DIFFERS** |

0.09% more spectra from the one-shot reader. The two paths pick peaks on different frame groupings,
so this is a second detector-like choice hiding in an option named for performance. **Which reader
is closer to the truth is UNMEASURED** -- neither arm has been through a search engine. What can be
said is that the shipped default is the arm every number in this file was taken with, and that the
alternative costs 2x the wall and 1.6x the memory. The help text now says CHANGES OUTPUT.

## Measurement methodology (2026-07-28) — variance model + estimator fix

* **BOTH engines are DETERMINISTIC: run-to-run sigma = 0.** Re-searching one mzML (dataset D corr_power=2)
  3x: Sage 10,802/10,802/10,802 AND MSFragger 11,560/11,560/11,560. So ALL benchmark peptide counts
  are EXACT; n=1 is justified for both. CONSEQUENCE: the MSFragger k-curve non-monotonicity
  (11,402/11,560/11,482/11,751 across k=1..4) is NOT run-to-run noise -- it is a REAL (non-smooth)
  response to the different inputs (now measured). The earlier "MSFragger
  is a noisy guardrail" framing was WRONG on the noise: MSFragger numbers are exact, just non-monotonic
  in k (its optimum k differs from Sage's). The ONLY stochastic uncertainty in the whole benchmark is
  the entrapment FDR estimate (small n), which the bootstrap CI below now quantifies.
* **Entrapment estimator was wrong TWO ways (fixed in `bench/entrapment.py`):**
  (1) it counted peptides SHARED between human and Arabidopsis as entrapment -> inflated the fraction
  ~2.5x (dataset D base entrap 222 -> correctly 86 foreign-ONLY). (2) it used the PROTEIN ratio (0.800); the
  correct PEPTIDE-hypothesis ratio (in-silico tryptic, Sage rules) is 0.6805, correction 1/r = 1.469.
* **Corrected absolute FDR at nominal 1% peptide_q is ~1.0-1.45%, NOT the ~2.3% previously reported** --
  the tool is WELL-CALIBRATED. With 95% bootstrap CIs (n_entrap ~ 60-90 -> noisy):
  dataset D base 1.34% [1.06-1.65] / corr2 1.20% [0.95-1.46]; dataset A base 1.01% / corr2 1.08%; dataset B base 1.45%
  [1.15-1.77] / corr2 1.18% [0.91-1.47]. **The arm FDR CIs OVERLAP** -> the earlier "FDR improved
  2.33->2.13%" was an estimator+noise artifact. Correct statement: **corr_power adds +9-10% target
  peptides at an entrapment FDR statistically INDISTINGUISHABLE from base** (more peptides, same FDR).
* **Consequence for future levers:** the +1-3% regime (im_weight increments, cap sweeps, rescoring
  features) needs this CI machinery + MSFragger replicate sigma to be falsifiable; do not decide small
  levers on point estimates.

## Decision rule

* **Progress** = Sage peptides ≥ 10,802 (standard bench config = tool defaults + `-assembly:require_isotope_support`; re-based 2026-09-01, was 10,072 pre-gate/pre-corr_power) and PSMs/peptide not worse (guardrail retained).
* **Regression** = peptides fall ≥ 1%.
* **Noise floor** = 0.06% same-node (from four runs of an effectively identical config:
  8,411 / 8,408 / 8,406 on one node). Differences below ~0.2% mean nothing.
* **Screen on Sage** (4.2× faster, signs reproduce); confirm winners on MSFragger before any
  default change. Sage is unreliable on effect SIZE — measured compression up to 3× on one arm.

## Confirmed WINS post-freeze (gate baseline: Sage 9,989 / MSFragger 11,033 / entrap 9,304 @ 2.33%)

* **`assembly:corr_power` + `im_weight_sigma` STACK — WIN on BOTH engines, entrapment-confirmed (2026-07-27).**
  Emitted fragment intensity *= co-elution^k. corr_power is EMIT-ONLY (excluded from the 500-cap key),
  so the corr_power arms preserve peak MEMBERSHIP (identical 924,255 spectra + same peak m/z set; pure
  intensity re-ranking). im_weight IS in the cap key, so with 73% of spectra at the 500-cap the STACK
  also REORDERS which fragments survive (its intended interferent-removal channel) -- same spectrum +
  peak COUNT, but not the same membership. So the stack's +310 over corr_power=2 has TWO channels
  (intensity reweight + cap membership), not one. The reference implementation on dataset D: Sage 11,517,
  MSFragger 13,932. MSFragger here is a NOISY guardrail (no replicate variance measured; screen on
  Sage+entrapment). The 2.33->2.13% entrap frac is WITHIN Poisson noise (n~85 entrapment) -- a guardrail
  (it passed: extra targets are not disproportionately entrapment), NOT an FDR win.
  | arm | Sage | MSFragger | entrap target @ actual-FDR |
  |---|---|---|---|
  | base (gate) | 9,989 | 11,033 | 9,304 @ 2.33% |
  | corr_power=1 | 10,622 | 11,402 | 9,894 @ 2.18% |
  | corr_power=2 | 10,802 | 11,560 | 10,235 @ 2.21% |
  | corr_power=3 | 10,901 | 11,482 | 10,288 @ 2.16% |
  | corr_power=4 | 10,875 | 11,751 | 10,262 @ 2.15% |
  | **corr_power=2 + im_weight 0.005** | **11,112 (96.5% of dt)** | **11,813 (84.8% of dt)** | **10,574 @ 2.13%** |
  corr_power alone PEAKS at k≈3 on Sage; the im_weight STACK is best on EVERY metric incl. the LOWEST
  entrap FDR (2.13%) -- the two are COMPLEMENTARY purity axes (RT co-elution vs IM), they ADD not erase.
  Sage 86.7% -> 96.5% of the reference implementation; MSFragger 79.2% -> 84.8%. Mechanism: a true fragment co-elutes with
  its precursor by construction (high c), so c^k down-weights mediocre-c INTERFERENTS. This OVERTURNED a
  UNANIMOUS prior prediction of -2..-5% MSFragger (premise "faint precursor -> low-c fragment"
  was wrong). MSFragger is NOISY here (non-monotonic) -- screen on Sage+entrapment.

  **GENERALISATION TEST PASSED for corr_power (frozen k=2/sigma=0.005 on untouched dataset A, dataset B; 2026-07-27):**
  Sage nominal 1%: dataset A base 8,822 -> k2 9,735 (+10.4%) -> stack 10,428; dataset B base 8,311 -> k2 9,032 (+8.7%)
  -> stack 9,159. So **corr_power=2 robustly generalises: +8-10% over base on ALL THREE files, both
  engines (MSFragger dataset A +5.6%, dataset B +7.6%), at nominal 1% AND <=2% empirical FDR.** The im_weight STACK
  increment is FILE-VARIABLE (+310 dataset D, +693 dataset A, +127 dataset B over k2; flat on MSFragger) -- a Sage-leaning
  bonus, NOT a robust default increment (dataset B +127 < the pre-registered +200 threshold). The <=1% empirical
  FDR numbers are erratic (k2 dips on dataset A, stack dips on dataset B) = estimator noise at ~50 entrapment n (the
  <=2% and nominal metrics are clean). **DECISION (pending review of the result): ship corr_power=2
  as the new default; im_weight_sigma=0.005 as a documented optional stack. Caveat: dataset A/dataset B are the same
  study series as dataset D, so this de-risks the dataset D family, not diaPASEF broadly -- a different-cohort file is
  the stronger (unavailable-today) test.**

## What is already falsified on this file or dataset B

Do not re-propose without new evidence: `competitive`, rank-pruning (`rp_max` 2/4/8),
`assembly:apportion`, `trace:frame_aggregation_n` (3, 5), `charge:ambiguity_margin` (0.5/1.0/2.0
— and it was *unreachable* in count mode, so it is untested rather than falsified),
`consolidate:delta_rt` (0.7 → −12.7%, 3.0 → −24.0%), **`merge:rt_window` 1.4 → −4.4%**,
mass recalibration, a monotonic ion-mobility charge prior, AlphaDIA-style log-sum co-elution
scoring (three implementations), **`gate:variance_support` (union-support Pearson — 2026-07-27:
DROP on unanimous analysis; it measures support co-LOCALISATION not co-elution SHAPE and reads
faint fragments' censored shoulders as anti-correlation, e.g. precursor-8/frag-5/overlap-3 scores
G-Pearson +0.47 but union −0.50. G-Pearson already ≈ robust cosine, so the change was cosmetic).**

## Spectrum COLLAPSE (same precursor across cycles) — FALSIFIED for quality (2026-07-28)
Tested on the OPEN-search entrapment results (simulation on existing PSMs, no re-search): group rank-1
PSMs by precursor coordinate ACROSS RT (charge + expmass 10 ppm + IM 0.01, single-link chained over a
10 s RT gap), keep one (or top-N) per group, RECOMPUTE the corrected entrapment FDR from scores.

| rule | spx kept | reduction | targets @1% FDR | dt targets |
|---|---|---|---|---|
| full            | 558,255 | 1.00x | **8,273** | 9,743 |
| apex (ID-agnostic) | 86,842 | 6.43x | 3,342 (-60%) | 4,491 (-54%) |
| best (ORACLE)   | 86,842 | 6.43x | 4,390 (-47%) | 5,089 |
| top3 / top5 / top10 (collapse + chimeric report_psms) | 132k/155k/190k | 4.2x/3.6x/2.9x | 5,448/5,767/6,603 | 6,068/6,672/7,422 |

**The multiple-testing-burden hypothesis is REFUTED.** Every arm's FDR CI sits at ~0.6-1.3%: collapse
removes true AND entrapment PSMs PROPORTIONALLY, so the FDR does not improve -- peptides are simply
lost. Emitting fewer spectra does NOT buy open-search sensitivity.
**It is INTRINSIC, not a DIAspeXtract defect:** the reference implementation collapses the same way (-24% at top10, -54% at
apex) and BOTH tools reduce to the SAME ~87k precursor features (spx 86,842 / dt 87,112) -- DIAspeXtract
just emits 1.33x more spectra per feature (6.43 vs 4.84). The across-cycle multiplicity is INDEPENDENT
EVIDENCE (each cycle = a different noise realisation + co-isolation mixture, another shot at a different
co-eluting peptide), not redundancy. Consistent with the long-standing "across-cycle spread 2.66 is
CORRECT (MS1 FWHM = 2.61 cycles)" finding.
**Residual option, untested:** TRUE merging (union of peaks -> better S/N) could in principle ADD IDs
rather than select among existing ones; this simulation only bounds the SELECT-among-existing case. But
it must beat a -20%..-60% deficit and it blends interferents across RT, so it is not promising.
**Practical use:** collapse is a SPEED/recall knob, not a quality lever -- top10 gives 2.94x fewer
spectra (much cheaper search) for -20% peptides, if total pipeline time ever matters more than recall.

## Monoisotope -1.003 Da open-search artefact: STOPPING guard FALSIFIED, scoring fix under test (2026-07-28)
The open-search re-baseline found our largest delta-mass artefact is a **-1.003 Da population** (top bin
7,616 PSMs; 8.20% isotope-shifted vs the reference implementation 3.46%): the `charge:scoring=count` walk defines the
monoisotope as the FURTHEST-LEFT peak it can reach, and `findPartner` matches m/z + RT + IM but **never
intensity**, so it latches onto noise or a co-eluting species one isotope below the true mono. Closed
search HIDES this (Sage enumerates `isotope_errors [-1,2]` and corrects it for free); open/blind search
turns it into a phantom -1 Da modification.

**Attempt 1 — `charge:mono_averagine_guard` (a STOPPING rule): FALSIFIED.** Require the leftward
candidate to reach the averagine minimum 1/lambda times the previous peak. At slack 0.5 on dataset D:
| metric | control | guard 0.5 |
|---|---|---|
| -1 Da bin | 7,616 | **2,734 (-64%, the target WORKED)** |
| +1 Da bin | 4,697 | **23,994 (+411%)** |
| total isotope-shift | 8.20% | **16.43% (WORSE)** |
| open targets @ FDR | 8,887 @1.29% | 8,939 @1.28% (flat, CIs overlap) |
| closed Sage | 10,802 | 11,697 (+8.3%) |
| emission | 924,255 | 1,005,353 (+8.8%) |
**A stopping rule can only trade over-reach for under-reach**, and worse: a walk that halts early leaves
the M+1 peak unconsumed by `used[]`, so it later seeds its OWN precursor one isotope too heavy -- that is
BOTH the +1 explosion and the +8.8% emission (one cause, two symptoms). Open search did not improve.
NOTE the closed-search 11,697 (> the reference implementation 11,517) is NOT a quality win to bank: peptides/1k spectra is
11.7 -> 11.6, i.e. the gain tracks the extra emission -- exactly the over-emission-rewards-closed-search
corruption. Do not ship the guard for it.

**Attempt 2 — `charge:mono_averagine_select` (a SCORING decision): under test.** Walk the full run, then
choose WHICH peak is the mono by averagine cosine of the envelope starting there. Corrects BOTH
directions, and every found peak is still consumed, so de-isotoping and emission stay unchanged (only the
reported precursor mass moves) -- which also removes the emission confound from the experiment.

### Attempt 2 result + the finding that STOPS this line of work (2026-07-29)
`charge:mono_averagine_select` raw cosine repeated the envelope scorer's documented short-vector trap
(`len` shrinks as the candidate mono moves right -> short vectors score high -> too-HEAVY mono). Adding
the same evidence weight the envelope scorer uses, (len-1)/3 capped at 1, largely fixed that:
| metric | control | guard 0.5 | msel raw | **msel + evidence weight** |
|---|---|---|---|---|
| emission | 924,255 | 1,005,353 | 924,716 | **924,552** |
| -1 bin | 7,616 | 2,734 | 2,314 | **4,015 (-47%)** |
| +1 bin | 4,697 | 23,994 | 21,591 | **12,188** |
| unmod | 31.2% | 24.3% | 26.4% | **32.0% (best)** |
| isotope-shift | **8.20%** | 16.43% | 15.24% | **9.48%** |
| open targets @FDR | 8,887 @1.29% | 8,939 @1.28% | 8,916 @1.33% | 8,932 @1.33% |
| closed Sage | 10,802 | 11,697 | 11,458 | 11,279 |
Pre-registered criterion was total isotope-shift BELOW 8.20% with both +-1 bins down: **NOT MET (9.48%)**.
Stopped rather than running a 4th variant. Both flags stay default-OFF, documented, for reproducibility.

**THE FINDING THAT MATTERS — the -1 artefact is NOT the open-search recall lever (the hypothesis was WRONG).**
Open-search targets are FLAT across all four arms (8,887 / 8,939 / 8,916 / 8,932; every CI overlaps)
even though the -1 population was cut 47-70%. Reason, obvious in hindsight: **a -1-shifted precursor is
still IDENTIFIED** -- it appears in the delta-mass histogram AS an identified PSM at -1.003 Da. Open
search does not lose the peptide, it mislabels it with a phantom modification. So the monoisotope
artefact is a **PTM-INTERPRETATION problem (real, and it matters for blind-PTM correctness, our actual
differentiator), NOT a recall problem.** It cannot explain the 8,887 vs 9,989 open-search gap; do not
spend more effort here expecting recall. The open-search recall gap needs a different explanation --
candidates not yet tested: our 924k vs the reference implementation's 700k emission competing in the huge open search
space, fragment quality on the faint tail, or precursor mass ACCURACY (ppm) rather than isotope offset.

## STEP 01 (2026-09-01): the metric re-denominated + denominators reconciled

### A. Delta-mass-BIN-level open-search scoring (bench/open_ptm_score.py v2, corrected estimator)
Unit = (peptide, delta-bin); classes scored separately: unmod (|d|<=8 mDa), nearzero (8-900 mDa =
mass-error tail, NOT a mod), artifact (iso lattice +-k*1.00335, k=1..9), known (curated Unimod list),
other (off-lattice, per-BIN FDR walks for top bins -- pooled class FDR let junk bins ride; the spx
-485..-499 Da window-edge wall at ~36% raw entrap collapses to ~30-80 @1% under per-bin walks).
| class @1% corrected FDR | DIAspeXtract | the reference implementation | ratio |
|---|---|---|---|
| unmodified | 4,885 | 7,015 | 70% |
| known-PTM (the raison-d'etre number) | 2,692 | 3,092 | **87%** |
| isotope-artifact (k<=9) | 4,073 | 3,611 | 1.13x (v1's "2.4x artifact" was a k<=3-grid artifact of its own: dt's big envelope ladder is ON-lattice at k>=4; ours skews to iso-1) |
| **nearzero mass-error (8-900 mDa)** | **7,099** | **2,814** | **2.5x -- NEW finding** |
| other (off-lattice) | 44,143 | 78,244 | 56% (interpret with care: correlated bins; contains e.g. -0.984 amidation-type masses off the iso lattice) |
**Findings:** (1) In the honest currency the reference implementation leads every legitimate class; the old
"89% targets@FDR" UNDERSTATED the unmod gap (70%) and roughly matched known-PTM (87%).
(2) NEW: spx's near-zero mass-error tail is 2.5x dt's, and spx unmod+nearzero = 11,984 vs dt 9,829 --
much of our "missing unmodified" population is IDENTIFIED but MASS-DEGRADED (precursor 8-900 mDa off),
i.e. an interpretation-quality defect (like the -1 Da class), pointing at PRECURSOR MASS ACCURACY.
NOTE: "mass recalibration" sits in the falsified list, but it was falsified ON THE CORRUPT CLOSED
METRIC -- exactly the class of falsification the 2026-09-01 review said may need re-audit. Flagged for
re-evaluation in THIS currency; not re-litigated yet. (3) The phantom-mod gap (artifact class) is
+13%, not 2.4x -- the earlier headline overstated it by using a too-short isotope grid.

### B. Denominators reconciled (the "no ratio is citable" blocker)
- **CANONICAL the reference implementation references from 2026-09-01 on: Sage 11,517** (dt_s30/results.sage.tsv,
  nominal peptide_q<=0.01 -- re-verified today; TSV archived) and **MSFragger 13,932**
  (msf_dt/fixed.tsv via score_msf_td target-decoy; archived).
- **11,552 / 13,014 are DEPRECATED**: products of the 2026-07-22 joint common-FDR procedure whose
  primary outputs (joint_results.json et al.) exist in no file -- doc-only numbers.
  Not to be used in ratios. The frozen table's 86.4%/80.5% ratios are legacy.
- **Decision-rule re-base:** the frozen threshold 10,072 predates the isotope-support gate AND
  corr_power (superseded twice). New reference for progress/regression: **Sage 10,802** = current
  tool defaults (corr_power=2) + the standard bench config (which passes -assembly:require_isotope_support
  explicitly; that flag remains OFF in the tool). Noise floor and +-1% rules unchanged.

### STEP 01 CLOSURE (v3 scorer, 2026-09-01) — SUPERSEDES the v2 table above
The v2 class table is NOT citable (review: Da-constant tolerance made every ratio a mass-accuracy
measurement biased against spx; per-bin @1% on zero-entrapment slivers vacuous; ontology gaps).
v3 (bench/open_ptm_score.py @ e9f5b64b): nearest-candidate assignment over {0, iso-lattice k<=15,
28 KNOWN masses}, ppm-scaled tol max(8 mDa, 12 ppm), conservative (e+1)/r tie-grouped walks with
bootstrap CIs, min-evidence rule (per-name @1% only when accepted e>=10), provenance headers.

**Citable v3 table (dataset D open search, @1% corrected class-FDR, CIs in ptm_score_v3.out):**
| class | spx | dt | reading |
|---|---|---|---|
| unmod | 6,907 | 7,890 | 87.5% at hypothesis level; PEPTIDE level ~parity (96-106% by definition: accepted-union 7,744 vs 8,079; best-hypothesis 6,547 vs 6,164) |
| knownPTM | 4,878 | 4,651 | "105%" is CONTAMINATED: spx's top-3 knownPTM bins are all iso-AMBIGUOUS (Didehydro/Amidation/Deamid = 34% of spx knownPTM units vs dt 14%) = our lattice spill relabeled. Ambiguity-excluded units: spx ~78%. Honest: BROADLY COMPARABLE SCALE, per-mod 57-131% (Ox 93%, Acetyl 131%, Phospho 57%), NO single mod FDR-certifiable (all per-name e<10) |
| isotope lattice (k<=15) | 7,737 | 8,826 | dt carries MORE lattice hypotheses (its +k up-ladder); spx signature = iso-1 (1,790 @1%, e=11, certified). NEITHER tool is clean |
| nearzero (mass-error) | 5,415 | 1,432 | 3.8x — THE spx-specific defect, mechanism nailed (below) |
| other | 40,160 | 72,747 | DESCRIPTIVE ONLY: dominated in BOTH tools by off-lattice integer-Da bins (dt's ladder to +21); not PTM discovery |

**The nearzero mechanism is PROVEN (three independent legs, all on existing data):**
(1) 70.1% of spx nearzero peptides ARE dt's unmod peptides; (2) 89.4% same-scan concordance with our
OWN closed search (same scan -> same peptide; 0.7% different -> chimeric explanation dead);
(3) deltas 93.8% negative, median -16.2 mDa, slope -18.3 ppm vs calcmass (dt control ~symmetric).
=> DIAspeXtract's reported precursor masses carry a SYSTEMATIC ~-18 ppm CALIBRATION BIAS introduced by
our pipeline (same raw data as dt). Fixable (per-run linear recalibration of reported mono m/z);
distinct from the old falsified "mass recalibration" lever (different target, corrupt metric).
Scheduled as its own arm AFTER step 02 (do not confound the pre-registered arm).

**MSFragger denominator reconciled:** 13,932/11,560 are PEPTIDOFORM-level (pair them);
sequence-level recount from archived TSVs: dt 12,752 vs spx 11,055 = 86.7%; legacy 13,014 was
sequence-level and is within 2% of 12,752 — unit mismatch closed, ratios now like-for-like.

**Step-01 bottom line:** the old "89% on open search" dissolves. At peptide level we are at ~parity
on unmodified; known-PTM scale is comparable (ambiguity-limited both ways); BOTH tools' open-search
delta dimension is dominated by non-PTM structure (lattice + integer-Da junk = 5-9x the legitimate
classes) — a tool-CLASS finding worth publishing; and DIAspeXtract's one distinctive open-search defect
is the -18 ppm mass bias (fixable), vs the reference implementation's up-ladder envelope artifacts. Step 02 (emission
arm) runs AS REGISTERED with a nearzero-stratified secondary readout, scored under this v3 metric.

## THE -8 ppm MASS BIAS: mechanism found + fixed (2026-09-01, follows the step-01 nearzero finding)
**Every DIAspeXtract-reported m/z to date carries a systematic, m/z-dependent -5..-11 ppm error**
(median -8.3 ppm on identified precursors; dt control -1.1 ppm on the same raw data).
**Mechanism (proven by direct comparison against Bruker's timsdata library via ctypes):** the loader
(BrukerTimsFile -> opentims++) converts TOF->m/z with `OpenSourceTof2MzConverter (linear-in-sqrt)` --
a TWO-POINT chord fit to just (MzAcqRangeLower=100, MzAcqRangeUpper=1700, DigitizerNumSamples): error
~0 at the anchors, bowing to -10.8 ppm mid-range, frame-independent -- reproducing every observed
feature (flat in RT and mass-in-ppm; the "IM trend" is the m/z-dependence in disguise; temperature was
a red herring, dC2=0 and dT1 tiny). The loader's SDK branch (OPENMS_BRUKER_SDK_PATH) upgraded ONLY the
ion-mobility converter -- the m/z converter was never swapped (and the env var was never set anyway).
**Fix: patches/openms-brukertims-sdk-mz.patch** -- 3 lines: the SDK branch now also installs
`BrukerTof2MzConverterFactory` (vendor per-frame calibration via tims_index_to_mz). Requires
OPENMS_BRUKER_SDK_PATH pointing at a local Bruker timsdata .so (if available locally;
used locally, never redistributed; without it the tool falls back to the open model as before).
The open-model inaccuracy is itself reportable upstream (opentims; also explains why mzdata disabled
its own m/z model port).
**Why mzPeak was NOT the answer (researched 2026-09-01):** spec
is draft 0.9 / prototype v0.1.0 / breaking layout change 2026-07-14 / no releases; the current
TDF->mzpeak converter stores RAW TOF indices and its own m/z path uses the same two-point model
(mzdata implemented the true vendor model 0.66.3 but disabled it 0.66.5 as "not consistently better");
the OpenMS in-tree MzPeakFile is unmerged/dormant and keyed to a superseded spelling. mzPeak re-poses
the calibration problem; it does not solve it.
**Impact on prior results:** [CORRECTED by review: the original "closed counts unaffected" claim
here was falsified by sdkcal itself, +12.3% -- candidate-window inclusion was never the only channel];
step-01 v3 open-search numbers stand AS MEASURED but the nearzero class (5,415 vs 1,432) was this bug.
Decisive arm `sdkcal` running: vendor-calibrated re-extraction -> expect ppm median -> ~-1, nearzero
to collapse toward dt's level, and possibly small closed gains (tighter effective tolerance).

### sdkcal CONFIRMED (2026-09-01): vendor m/z calibration is the new benchmark configuration
| metric | pre-fix (corr2) | **sdkcal (corr2 + vendor m/z)** | the reference implementation |
|---|---|---|---|
| ppm median (open IDs) | -8.29 [-12.5,-4.0] | **+1.74 [-2.4,+6.1]** | -1.10 [-3.8,+1.2] |
| closed Sage | 10,802 | **12,128 (+12.3%)** | 11,517 (scoped claim only -- see review verdict below) |
| closed Sage entrapment | 10,381 @1.20% | **11,435 @1.26% [1.00-1.51]** | -- (FDR unchanged -> gain is real) |
| closed MSFragger | 11,560 | 11,460 (-0.9%, FLAT) | 13,932 (82.3%) -> Sage-specific gain |
| emission | 924,255 | 922,132 (~unchanged) | 700,434 |
| open unmod @1% | 6,907 (87.5%) | **7,566 (95.9%)** | 7,890 |
| open nearzero @1% | 5,415 (3.8x dt) | **2,769 (1.9x)** | 1,432 |
| open knownPTM @1% | 4,878 | 6,058 (total; ambiguity-excluded units ~85% of dt) | 4,651 |
**Decision: OPENMS_BRUKER_SDK_PATH + patches/openms-brukertims-sdk-mz.patch = REQUIRED benchmark
config from now on; closed-Sage reference re-based 10,802 -> 12,128.** Caveats carried into review:
(1) the +12% is Sage-specific (MSFragger flat -- its gap is search/detection-side, consistent with
prior reviews); (2) EVERY pre-2026-09-01 number was measured at -8 ppm, INCLUDING the corr_power
validation and all falsifications -- corr_power=0 ablation under vendor-cal launched (sdkcal_cp0);
dataset A/dataset B holdout re-validation owed; falsified-lever re-audit only where a lever plausibly interacted
with mass accuracy; (3) knownPTM ambiguous fraction unchanged (34.7% vs dt 14%) but the bins are now
physically separable at +-4 ppm -- the 1,610-unit Amidation bin (1.4% raw entrap) deserves real
chemical scrutiny rather than automatic dismissal as spill.

### Calibration-step verdict (2026-09-01) — 12,128 is PROVISIONAL
Convergent verdict: the chord mechanism is the proven DOMINANT cause of the old raw-axis bias, but the
step over-claimed in four ways, now corrected:
1. **Attribution not isolated:** setting the env var flipped BOTH converters (m/z AND ion mobility) --
   the SDK IM swap engaged for the first time ever alongside the m/z patch. Part of the +12.3% could be
   IM-side; ALL prior IM-dependent results (incl. im_weight's validation) silently ran on the open IM
   model. IM-isolation arm required (m/z vendor + IM open).
2. **12,128 = provisional dataset D/vendor-SDK regression oracle, NOT a validated reference** until: cp0
   ablation (running), the IM-isolation arm, and dataset A/dataset B re-validation under vendor-cal.
3. **"Leads the reference implementation" is not citable.** Maximum defensible wording: "On dataset D, using Sage and
   the vendor-SDK calibration configuration, DIAspeXtract yielded 5.3% more closed-search peptides; under
   MSFragger it yielded 17.7% fewer." Per-1k-spectra efficiency: spx 13.2 vs dt 16.4 (dt +25%) -- the
   emission column stays in any public row. MSFragger flatness is UNINFORMATIVE (its default
   calibrate_mass auto-corrects input bias), not evidence against the fix.
4. **Residual defect remains:** sdkcal centers +1.74 ppm vs dt -1.10 (2.8 ppm apart, quartiles
   symmetric -- earlier "right-skewed" claim wrong) and nearzero is halved, not closed (1.9x dt).
   Candidates: our centroiding/mono reporting point. Paired same-feature spx-vs-dt reported-m/z
   comparison is the discriminator. The -18 ppm (step-01 regression slope on the truncated nearzero
   subset) vs -8.3 ppm (median, all IDs) are different statistics of the same defect -- reconciled here.
Hardening obligations recorded: benchmark must FAIL CLOSED on a missing/wrong SDK path (silent
fallback = known-wrong masses); opentims' BrukerTof2MzConverter DISCARDS tims_index_to_mz's return
code (silent incomplete-buffer risk -- fix when the patch is vendored properly); SDK identity pinned
in evidence/SHA256SUMS (sha a9708613..., PRERELEASE TDF-SDK 2.21.0 -- Bruker notes a tims_index_to_mz
implementation change in this line); patch base = OpenMS 3.6.0-pre-exported-20260717 tarball; runtime
headlines must be RE-MEASURED under the SDK path before publication. Determinism check (threads-8 vs
threads-100 semantic/byte hash) running. Two-tier calibration product conflict stands: license-blocked
users (the stated market) get -8 ppm masses from the open fallback -- an open true-polynomial port or
loud two-tier disclosure is a v0.2.0 gate.

### cp0 ablation result (2026-09-01): corr_power SURVIVES the calibration fix — provisionality condition 1/3 cleared
corr_power=0 + vendor-cal = 11,017 vs corr_power=2 + vendor-cal = 12,128. The 2x2 (closed Sage):
old-cal cp0 9,989 / cp2 10,802 (+8.1%); vendor-cal cp0 11,017 / cp2 12,128 (+10.1%). The levers are
independent and additive (mild positive synergy); corr_power's +8-10% validation was NOT a
calibration artifact. Note the cp0 arm also ran vendor IM, so the m/z-vs-IM attribution of the
calibration gain itself still awaits the IM-isolation arm (condition 2/3); dataset A/dataset B = condition 3/3.

### Table-model calibration: 3-FILE RESULT (open path, no vendor .so; 2026-09-01)
| file | old chord | **TDF-table model** | gain | vendor-SDK | ppm (table) |
|---|---|---|---|---|---|
| dataset D | 10,802 | **11,976** | +10.9% | 12,128 | +1.49 [-2.4,+5.4] |
| dataset A |  9,735 | **10,333** | +6.1%  | -- | +3.19 [-0.8,+7.1] |
| dataset B |  9,032 | **9,891**  | +9.5%  | -- | +2.29 [-1.7,+6.3] |
**The calibration gain GENERALISES: +6-11% closed Sage on all three files with NO vendor library.**
Provisionality conditions on the re-based reference: (1) cp0 ablation CLEARED; (2) IM-isolation
ANSWERED FOR FREE -- the table arm is vendor-exact m/z + OPEN (rational) IM, the sdkcal arm was
vendor both: dataset D 11,976 vs 12,128 => ~93% of the +12.3% is m/z-side, ~1.3% (152 peptides) is the
vendor-IM contribution; (3) dataset A/dataset B re-validated above. **New open-path reference: dataset D = 11,976**
(the 12,128 SDK number stays as the vendor-oracle upper bound, not the shipping config).
Residual +1.5..+3.2 ppm persists across all files and BOTH calibration paths -> consistent with the
centroiding/mono-reporting hypothesis (Option I), not the m/z axis. the reference implementation per-file ppm + dataset A/dataset B
head-to-head running to test whether part of it is instrument/per-file rather than ours.

### FIRST 3-FILE HEAD-TO-HEAD vs the reference implementation (Sage closed, corrected calibration, 2026-09-01)
the reference implementation had ONLY ever been measured on dataset D; dataset A/dataset B references created here under identical config.
| file | DIAspeXtract (table model) | the reference implementation | ratio | our ppm | dt ppm |
|---|---|---|---|---|---|
| dataset D | **11,976** | 11,517 | **104.0%** | +1.49 | -1.37 |
| dataset A | **10,333** | 10,242 | **100.9%** | +3.19 | +0.27 |
| dataset B | **9,891**  |  8,948 | **110.5%** | +2.29 | -0.60 |
**DIAspeXtract reaches 100.9-110.5% of the reference implementation on closed Sage across the three files** -- no vendor
library, open BSD path. STATISTICALLY (n=3 paired): mean 105.1%, 95% CI [93.0, 117.3] = consistent
with PARITY, not an advantage; the earlier phrasing "matches or exceeds on ALL THREE files" was
over-claiming on dataset A's +91. CORRECTION: the "ours 1.26% [1.00-1.51]" entrapment figure
quoted here came from the **sdkcal** arm, a different converter -- **the shipping table-model arm's
entrapment FDR has NOT been measured on any file**; that is an owed measurement, not a passed check.
CAVEATS THAT STAND: (a) emission is still ~1.32x dt, so per-spectrum efficiency still favours dt --
publish the emission column with any count row; (b) MSFragger arm pending (dt led 82% there; that
engine's gap is search/detection-side); (c) single cohort/instrument -- public-PXD replication is
still the external-claim gate; (d) our ppm residual (+1.5..+3.2) is 2-4 ppm ABOVE the reference implementation's
(-1.4..+0.3) on every file: dt's per-file spread (-1.37/+0.27/-0.60) is small and centred, so the
residual is OURS (centroiding/mono-reporting, Option I), NOT a per-file instrument property.

### PAIRED m/z bias vs the reference implementation (2026-09-02, 15:14) -- the residual is ours, in BOTH directions
Same (peptide, charge) identified by both tools at 1% (best PSM per tool), so selection bias is excluded:
| file | shared pairs | precursor ppm ours / dt / **paired ours−dt** | fragment ppm **paired ours−dt** | by charge (paired, precursor) | by m/z tercile |
|---|---|---|---|---|---|
| dataset A (tbl_s08 vs dt3_s08) | 10,091 | +5.17 / +2.13 / **+2.63** [p10 −2.6, p90 +10.3] | **−1.99** [−4.5, +0.5] | z2 +3.12, z3 +1.77, z4 +1.96 | low +1.87, mid +3.35, high +2.62 |
| dataset B (tbl_s23 vs dt3_s23) | 9,094 | +3.78 / +1.97 / **+1.48** [−3.4, +8.4] | **−2.08** [−4.4, +0.3] | z2 +1.96, z3 +0.53, z4 +0.29 | low +0.42, mid +1.88, high +1.87 |
| **dataset D** (d2_P1 = HEAD table arm vs dt_s30_closed, fresh the reference implementation Sage run = 11,517 peptides, reproducing the matrix reference) | 11,191 | +3.58 / +2.05 / **+1.12** [−3.7, +7.7] | **−2.41** [−4.7, +0.1] | z2 +1.59, z3 +0.31, z4 +0.54 | low +0.05, mid +1.62, high +1.49 |
Reading: on the same raw file and the same exact m/z axis, our REPORTED precursor m/z sits +1.5..+2.6 ppm
above the reference implementation's and our REPORTED fragment m/z sits ~2 ppm below -- two different centroid/reporting
conventions (precursor: MS1 trace m/z estimator, larger at z=2 and in the mid/high m/z terciles; fragment:
MS2 trace m/z estimator), not a calibration offset (which would move both the same way). Both feed Sage's
discriminant (precursor_ppm, average_ppm) and the 20 ppm fragment tolerance. Next: locate the two estimators
(`Trace.mz` for MS1 and MS2 traces: intensity-weighted mean over the trace vs apex) and A/B an apex /
top-k-weighted estimator behind the dataset D set gate, reporting ppm beside the count. dataset D pair running
Three-file pattern: the precursor bias is z=2-dominated (z3/z4 near zero on dataset D/dataset B) and absent in the lowest m/z tercile -- an isotope/mono-envelope effect (our reported mono for z=2 mid/high-m/z precursors sits a fraction of a ppm-scaled isotope offset high?) rather than a uniform centroid shift; the FRAGMENT bias is uniform −2.0..−2.4 ppm on all three files and is the one that touches the 20 ppm search tolerance and Sage's average_ppm feature (−3 ppm uniform shift cost −165 peptides on 09-01). Estimator A/B (apex / median vs weighted mean) queued behind the current gate.

### MEASURED (2026-09-02, 16:41): entrapment FDR of the SHIPPING table-model arm
`entrap_apply.py` (peptide-hypothesis ratio 0.6805, corrected estimator) on a Sage entrapment search of the HEAD
output (d2_P1, table-model calibration, parallel loader -- byte-identical to every run today):
| arm | targets @1% | entrapment | raw % | **FDR %** | 95% CI |
|---|---|---|---|---|---|
| **table model (shipping, HEAD)** | 11,489 | 107 | 0.92 | **1.37** | 1.10-1.64 |
| vendor-SDK calibration (09-01 arm) | 11,435 | 98 | 0.85 | 1.26 | 1.00-1.51 |
| **dataset A table model (shipping)** | 9,736 | 63 | 0.64 | **0.95** | 0.72-1.18 |
| dataset A the reference implementation | 9,953 | 84 | 0.84 | 1.24 | 0.99-1.54 |
| **dataset B table model (shipping)** | 9,462 | 85 | 0.89 | **1.32** | 1.04-1.59 |
| dataset B the reference implementation | 8,581 | 71 | 0.82 | 1.22 | 0.96-1.53 |
The owed measurement is now made: the shipping arm's true FDR at nominal 1% is ~1.4%, inside the
same band as the vendor-calibrated arm (CIs overlap), so the open-path calibration gain is not bought with
FDR. On all three files the shipping arm's true FDR at nominal 1% is 0.95-1.37%, inside the reference implementation's band on the two files where it exists (1.22-1.24%): the 3-file closed-search parity claim now carries its FDR control.

### `assembly:apportion` FALSIFIED CLEANLY (2026-09-02, 16:46) -- the first A/B with corr_power actually applied
Every earlier apportion/rp_max A/B compared corr_power=2 (share-all) against corr_power=0 (the variant) because those
branches bypassed the emit weights (fixed 15:05). Re-run with the weights applied: apportion=1.0 -> Sage 11,742
vs 12,217 share-all (**-3.9%**, only-share-all 1,234 / only-apportion 759). Same verdict as July, now clean:
cross-precursor intensity apportionment loses peptides. Consequence: BACKLOG item 1 (NNLS unmixing onto the MS1
basis) is cut -- its prerequisite ("apportion shows a signal") failed under the correct measurement. Side note: the
apportion path ran the window loop at 6.1x (3,658 s) -- an hour-long run -- irrelevant now that it is dead.

### Trace m/z estimator A/B (2026-09-02, 17:22): APEX m/z is +2.6% peptides (same binary, env switch)
`SPEXTRACTOR_MZ_ESTIMATOR` selects the reported m/z of every trace in `toTrace()`; default = OpenMS centroid
(intensity-weighted mean over the trace's members).
| estimator | spectra | Sage @1% | vs share-all default (12,217) | paired precursor ppm vs dt | paired fragment ppm vs dt |
|---|---|---|---|---|---|
| mean (default) | 922,902 | 12,217 | — | +1.12 | −2.41 |
| **apex** (max-intensity member's m/z) | 927,813 | **12,537 (+2.6%)** | only-default 757 / only-apex 1,077 | **+0.78** | −2.41 |
| median | 920,021 | 12,099 (−1.0%) | 756 / 638 | +1.18 | −2.36 |
Same binary -> the 0.06% same-binary floor applies; +2.6% is real on Sage. The precursor bias moves a third of
the way to the reference implementation; the FRAGMENT bias does not move at all -> the −2.4 ppm fragment offset is upstream of the
trace estimator (the IM-cluster pick's centroid, or the reference implementation reporting the apex RAW peak), a separate A/B.
Default change is gated on the standing rule: **MSFragger on d7_apex before apex becomes the default** (queued).

### APEX estimator passes the both-engines gate -> NEW DEFAULT (2026-09-02, 18:42)
| dataset D | Sage @1% | MSFragger @1% (raw-hyperscore walk) | precursor ppm paired vs dt |
|---|---|---|---|
| table model, mean estimator (until today) | 12,217 | 11,463 | +1.12 |
| **table model, APEX estimator** | **12,537 (+2.6%)** | **11,927 (+4.1%)** | **+0.78** |
| the reference implementation | 11,517 | 13,932 | — |
MSFragger ratio vs the reference implementation 82.3% -> 85.6%. Same binary (env switch), so the same-binary floor applies; both
engines move the same way; ppm beside the count improved. Default changed in code (`SPEXTRACTOR_MZ_ESTIMATOR=mean`
restores the OpenMS centroid). Owed next: dataset A/dataset B apex arms (Sage) and entrapment on the apex arm (queued).

### FALSIFIED (2026-09-02, 19:05): dropping the precursor's M+1..M+3 from fragment lists LOSES peptides
`SPEXTRACTOR_DROP_PREC_ISO=1` on the apex-default binary: isotope contamination of the emitted lists falls from
28.4% / 24.8% (M+1 / M+2 within 10 ppm) to 0.2% / 0.2%, and Sage drops **12,537 -> 12,249 (-2.3%**, only-apex
1,061 / only-dropiso 773; same binary, env switch). The unfragmented precursor isotopes in a pseudo-MS/MS spectrum are
useful to the engine (or their removal lets weaker peaks refill the 500-cap). Keep them; switch stays off.

### Pick-level m/z modes FALSIFIED (2026-09-02, 20:00) -- and the fragment offset is not in our estimator chain
`pickIMCluster:mz_mode` (OpenMS PeakPickerIM patch; default `weighted` = unchanged, verified byte-identical to the
apex-default run):
| pick m/z | spectra | Sage @1% | vs apex default (12,537) | paired precursor / fragment ppm vs dt |
|---|---|---|---|---|
| weighted (default) | 927,813 | 12,537 | identical | +0.78 / −2.41 |
| seed (most intense point) | 896,143 | 12,196 | −2.7% | +0.73 / −2.38 |
| top3 (weighted mean of 3) | 926,359 | 12,325 | −1.7% | +0.70 / −2.42 |
The −2.4 ppm fragment offset vs the reference implementation survives BOTH the trace estimator (apex/median/mean) and the pick centroid
(weighted/seed/top3) -- it is not produced by anything in our reporting chain after the TOF->m/z conversion, which is
the same exact model on both sides (2.5e-5 ppm). What remains: the reference implementation's own fragment m/z reporting/recalibration.
That is not a defect to chase on our side; the precursor axis is where our reporting differs (apex fixed a third).

### APEX estimator: 3-FILE RESULT + entrapment (2026-09-02, 20:57) -- the new open-path reference
| file | mean estimator (09-01 reference) | **apex (new default)** | gain | the reference implementation (Sage) | ratio | paired precursor ppm vs dt (mean → apex) | fragment |
|---|---|---|---|---|---|---|---|
| dataset D | 12,217 | **12,537** | +2.6% | 11,517 | **108.9%** | +1.12 → +0.78 | −2.41 |
| dataset A | 10,333 | **10,789** | +4.4% | 10,242 | **105.3%** | +2.63 → +1.47 | −1.95 |
| dataset B | 9,891 | **10,156** | +2.7% | 8,948 | **113.5%** | +1.48 → +1.17 | −2.13 |
Mean ratio 109.2% (n=3 paired; 09-01 was 105.1%). Entrapment of the dataset D apex arm: **1.28% [1.03-1.53]** at nominal 1%
(table/mean arm 1.37% [1.10-1.64]; the reference implementation 1.22-1.24% on dataset A/dataset B) -- the gain is not bought with FDR. MSFragger on
dataset D: 11,927 vs 11,463 (+4.1%; 85.6% of the reference implementation's 13,932). Runtime note: dataset A (1.23 M spectra, the largest file)
26:51 wall / 156 GB peak; dataset B 15:28 / 109 GB. The apex arm outputs are the new benchmark reference
(`d7_apex`, `d13_apex_s08`, `d13_apex_s23`); the mean-estimator numbers stay as the 09-01 row.

### Rescoring A/B, first two attempts (2026-09-02, 21:30)
FragPipe 24.0 headless (Basic-Search workflow, 25 ppm, our human_decoy FASTA) on the dataset D apex arm and on the reference implementation.
MSBooster (DIA-NN predictions) aborts on OUR pepXML both with and without spectrum prediction: "Prediction missing
in file for REM[15.9949]DQTM[15.9949]AANAQK|3" (a PSM whose peptide DIA-NN's predictor does not emit). the reference implementation's
arm completed: **18,670 peptides / 37,086 PSMs at 1% (Percolator + Philosopher sequential/picked)** vs 13,932
peptidoforms in the raw-hyperscore walk -- rescoring is worth +34% on the reference implementation's spectra, which is exactly why the
comparison must be made on both. Third attempt running: MSBooster off, Percolator on MSFragger features only, both
tools (like-for-like).

### RESCORING A/B, like-for-like (2026-09-02, 21:46): rescoring does NOT close the MSFragger gap
FragPipe 24.0 headless, Basic-Search (MSFragger 25 ppm, our human_decoy FASTA) -> Percolator -> Philosopher
(sequential, picked, 1% peptide), MSBooster OFF for both (it aborts on our pepXML, see above):
| dataset D arm | peptides @1% | PSMs @1% | PSMs / peptide | ratio spx/dt |
|---|---|---|---|---|
| DIAspeXtract apex | 13,211 | 113,685 | 8.6 | — |
| the reference implementation | 15,947 | 31,700 | 2.0 | **82.8%** |
| (raw-hyperscore walk, for reference) | 11,927 vs 13,932 | | | 85.6% |
| (the reference implementation + MSBooster RT+spectra) | 18,670 | 37,086 | 2.0 | ours unavailable |
Rescoring lifts both tools (+11% / +14%) and the ratio does not move (85.6 -> 82.8%). So the MSFragger deficit is
not a scoring-function artefact that a rescorer repairs: the 3,600 the reference implementation-only peptides have no competitive PSM
in our output at all. Together with the score comparison (identical scores on shared peptides) this points at CONTENT/coverage of the
faint tail -- emission competition or missing fragments -- not at search-side handling. Our 8.6 PSMs per peptide vs
the reference implementation's 2.0 is the over-emission signature: Percolator's peptide-level FDR gives redundancy nothing.
**Decision:** the next attribution arms are content-side: (c) the emission-controlled arm (pre-registered
step 02: quality-gate 924k -> ~700k + the reference implementation-precursor-matched sub-arm, both engines, full FDR curves) and (d) the
oracle-fragment arm. The learned charge/mono predictor stays conditional on (d). A MSBooster-tolerant rerun (var mods
capped identically for both tools) is owed for the +17% predictions bring the reference implementation.

### STEP 02, sub-arm 1 -- PRE-REGISTERED (2026-09-02, 22:00): precursor-quality gate by isotope-envelope depth
`SPEXTRACTOR_MIN_ISOTOPES=k` keeps precursors with >= k isotope peaks (k=2 is today's default via
require_isotope_support; k=3 and k=4 are the arms). Emission falls with k; fragment sharing is untouched.
Prediction if the MSFragger deficit is EMISSION COMPETITION: the spx/dt MSFragger ratio rises as emission falls
toward the reference implementation's 700k (peptides@1% hold or rise while spectra drop). Prediction if it is FAINT-TAIL CONTENT: the
ratio is flat or falls (the gate removes exactly the faint precursors). Gates: Sage set overlap vs the apex arm,
MSFragger raw walk vs msf_dt, emission count; both engines before any default change.

### STEP 02, sub-arm 2 -- PRE-REGISTERED (2026-09-02, 22:03): the reference implementation-precursor-MATCHED emission
`SPEXTRACTOR_PRECURSOR_LIST=<rt_sec mz z>` keeps only our precursors that match a listed the reference implementation dataset D precursor
(|dRT| <= 10 s, |dm/z| <= 10 ppm, same charge). The list is the reference implementation's own 700,434 pseudo-spectra (extracted from
its mzML; RT converted from minutes). This holds the PRECURSOR POPULATION fixed and varies only the spectra.
Prediction if the MSFragger deficit is EMISSION COMPETITION: with the same precursor set the spx/dt MSFragger ratio
moves toward parity (>= 95%). Prediction if it is PER-PRECURSOR CONTENT (fragment coverage of the faint tail): the
ratio stays near 85% even with matched precursors. Gates: emission count, Sage set vs the apex arm, MSFragger raw
walk vs msf_dt. Runs after sub-arm 1 (`bench_match.sh`, d15_matched).

### STEP 02 sub-arm 1 RESULT (2026-09-02, 22:34): EMISSION COMPETITION IS NOT THE MSFRAGGER MECHANISM
| dataset D arm | precursors kept | spectra | Sage @1% (vs apex 12,537) | MSFragger @1% (vs dt 13,932) | ratio | wall |
|---|---|---|---|---|---|---|
| apex (k=2, today's default) | 1,105,532 | 927,813 | 12,537 | 11,927 | 85.6% | 15:21 |
| **k=3** (>= 2 isotope partners) | 655,297 | **566,537 (−39%)** | 11,739 (−6.4%; only-apex 1,010 / only-k3 212) | **11,820 (−0.9%)** | 84.8% | 11:53 |
| k=4 | 399,837 | 353,689 (−62%) | 10,201 (−18.6%) | 10,538 (−11.6%) | 75.6% | 9:51 |
Pre-registered reading: cutting emission by 39% -- to BELOW the reference implementation's 700k -- moves the MSFragger ratio from 85.6%
to 84.8%. It does not rise. So the extra spectra are not what costs us on MSFragger (the emission-competition
hypothesis is refuted for the second time, now with a precursor-quality gate rather than merging), and the gate removes
real faint peptides on Sage (−6.4%). The deficit is per-precursor CONTENT of the faint tail, consistent with the
rescoring A/B (ratio unchanged after Percolator) and the score comparison (identical scores on shared peptides). Sub-arm 2 (matched
precursor population) is the direct test of that and is running. Side result: `SPEXTRACTOR_MIN_ISOTOPES=3` is a
legitimate speed/emission knob (−39% spectra, −23% wall, −0.9% MSFragger, −6.4% Sage) -- not a default.

### STEP 02 sub-arm 2 RESULT (2026-09-02, 22:49): SAME PRECURSORS, SAME 85% -- the deficit is per-precursor content
Reference: the reference implementation's 700,434 dataset D pseudo-spectra (RT/m/z/z). Of our 1,105,532 precursor hypotheses 445,008 match one
(10 s / 10 ppm / z); emitting only those gives 432,579 spectra.
| dataset D arm | spectra | Sage @1% | MSFragger @1% | MSFragger ratio vs dt |
|---|---|---|---|---|
| apex (all precursors) | 927,813 | 12,537 | 11,927 | 85.6% |
| **the reference implementation-precursor-matched** | 432,579 | 11,903 (103% of dt's 11,517) | **11,842** | **85.0%** |
| the reference implementation | 700,434 | 11,517 | 13,932 | — |
With the precursor POPULATION held to the reference implementation's, the MSFragger ratio does not move (85.0%). Together with sub-arm 1
(emission −39%: 84.8%), the rescoring A/B (Percolator: 82.8%) and the score comparison (identical scores on shared peptides), every
search-side and population-side explanation is now excluded. What is left is the CONTENT of the spectra we emit for the
precursors the reference implementation also emits: for ~3,600 peptides our spectrum of the same precursor does not reach MSFragger's
1% threshold while the reference implementation's does, and on Sage it does. Candidates, in test order: (1) the intensity reshaping
`corr_power=2` (a Sage win that MSFragger's rank-based hyperscore may dislike -- never measured on MSFragger; cp0 arm
launched); (2) fragment coverage of the faint tail (our ~55 fragments vs the reference implementation's ~500 per spectrum, top-150 cap in
MSFragger); (3) the oracle-fragment arm (true fragments at true charge/mono) to bound what content can give.
Note also: only 62% of the reference implementation's emissions have a matching precursor of ours -- the other 38% are either duplicates
on the reference implementation's side or precursors we never hypothesise (the MS1 funnel), a separate coverage question.

### CONTENT CANDIDATE 2 MEASURED (2026-09-03, 01:28): the faint tail is missing ~2-3 matched fragment ions per spectrum
MSFragger raw-walk TSVs (dataset D, best PSM per peptide+charge, target only; thresholds 15.2 ours / 14.0 dt):
| | apex arm (all precursors) | the reference implementation-precursor-matched arm |
|---|---|---|
| the reference implementation-identified / shared / **dt-only** | 14,974 / 11,037 / **3,937** | 14,974 / 11,367 / **3,607** |
| dt-only for which we HAVE a PSM (searched, below threshold) | **2,501 (64%)** | 1,567 (43%) |
| dt-only with NO PSM at all (never emitted or < min_frags) | 1,436 (36%) | 2,040 (57%) |
| dt-only, ours vs dt: hyperscore median | 12.5 vs 16.2 | 13.1 vs 16.2 |
| dt-only, ours vs dt: matched ions median (of 26 theoretical) | **7 vs 9** (deficit +2; p25 +1, p75 +4; same at every charge) | 6 vs 9 (deficit +3) |
| shared, ours vs dt: hyperscore / matched ions | 23.6 vs 22.9 / 12 vs 12 | 24.4 vs 22.4 / 12 vs 12 |
Reading: on the peptides we both find, our spectra are as good as the reference implementation's. On the ~3,900 we miss, two thirds
are searched and fall short by ~2 of ~9 fragment ions -- the faint tail is under-covered: the weaker fragments of
faint precursors are gated out (`gate:min_correlation` 0.3, `min_correlation_points`) or never traced
(`trace:ms2_noise_threshold_int` 10, `ms2_chrom_peak_snr` 1.0). The remaining third we never emit at all (the MS1
funnel: no hypothesis or < min_frags). This is the first mechanistic, per-peptide account of the MSFragger deficit.
**Pre-registered arms (running next, both engines):** (a) `gate:min_correlation 0.2`, (b) `trace:ms2_noise_threshold_int 5`,
(c) both. Prediction: the dt-only matched-ion deficit shrinks toward 0 and the MSFragger ratio rises above 90%;
falsifier: the deficit stays at +2 (then the missing ions are not in our traces at all -> MS2 aggregation / detector).

### CONTENT CANDIDATE 1 FALSIFIED (2026-09-03, 01:52): corr_power=2 helps MSFragger too
`-assembly:corr_power 0` on the apex arm: Sage 11,502 (−8.3% vs 12,537; only-cp2 1,640 / only-cp0 605) and
**MSFragger 11,429 (−4.2% vs 11,927)**. The correlation-power intensity reshaping is a win on both engines (first
MSFragger measurement of it; the July +8-10% was Sage-only). Intensity structure is not what MSFragger dislikes in our
faint-tail spectra. Remaining content candidate: fragment coverage (the +2 matched-ion deficit) -- arms running.

### COVERAGE ARMS FALSIFIED (2026-09-03, 03:10): the missing fragment ions are not in our traces at all
| dataset D arm | spectra | Sage vs apex 12,537 | MSFragger vs apex 11,927 | dt-only: ours vs dt matched ions |
|---|---|---|---|---|
| `gate:min_correlation 0.2` | 928,572 | 12,509 (145/117) | 11,935 | 7 vs 9, deficit +2 |
| `trace:ms2_noise_threshold_int 5` | 927,813 | **12,537, symdiff 0 (byte-identical)** | 11,927 | 7 vs 9, deficit +2 |
| both | 928,572 | 12,509 | 11,935 | 7 vs 9, deficit +2 |
Loosening the correlation gate adds nothing to the faint tail, and the intensity threshold 5 vs 10 changes NOTHING
(the SNR / min-length gates already decide). Per the pre-registered falsifier, the ~2 missing ions per faint spectrum are
not in our fragment TRACES: they are peaks that never become a trace (single-frame or sub-SNR) -- exactly the population
the reference implementation's ~500-peak spectra carry and our ~55-trace spectra do not. Next arm (content candidate 3): backfill each
pseudo-spectrum with the raw picked MS2 peaks of the precursor's apex frame inside its IM band (top-N not already
present), the DIA-Umpire/the reference implementation-style "peak-level" spectrum, both engines; falsifier: dt-only deficit stays +2.

### CONTENT CANDIDATE 3 -- PRE-REGISTERED (2026-09-03, 03:35): raw apex-frame peak backfill
`SPEXTRACTOR_BACKFILL_RAW=N` adds to each pseudo-spectrum the N most intense raw picked MS2 peaks of the frame nearest
the precursor's RT, inside its 1/K0 band (delta_im), not already present within 10 ppm -- untraced, single-frame
content (the population the coverage arms showed we never trace). Arms: N=50, N=150. Predictions: if the faint tail's
missing ions are these peaks, the dt-only matched-ion deficit (+2) shrinks and the MSFragger ratio rises (>= 90%);
Sage may fall if the added peaks are noise (the 150-peak search cap protects MSFragger more than Sage). Falsifier:
deficit stays +2 with N=150 -> the missing ions are not in the picked frame either (pick-level / MS2 aggregation).

### Backfill arm, first run INVALID (2026-09-03, 04:42) -- and the never-searched third attributed
N=50 and N=150 came back byte-identical to apex (Sage symdiff 0, MSFragger 11,927): the backfill never fired because
`detectTraces_` partitions the window's PeakMap into bands and clears it, so the map is empty by assembly time. Fixed
(a retained copy when the arm is on); rerun in flight. Not a result.
**MS1 funnel (funnel_dtonly.py on the apex arm's hypothesis dump, 1,105,532 hypotheses):** of the 1,436 the reference implementation-only
peptides for which we have NO PSM, 654 (46%) have a precursor hypothesis of ours within 10 s / 10 ppm / z (918 = 64%
within 30 s) -> hypothesised but no spectrum reached the search (assembly: < min_frags or no correlated fragments);
~500 (35%) have no hypothesis at all (MS1 trace / envelope). So the 3,937 dt-only peptides split ~64% searched-but-
under-covered / ~17% hypothesised-not-emitted / ~19% never hypothesised.

### ASSEMBLY-LOSS ARM -- PRE-REGISTERED (2026-09-03, 04:45): `assembly:min_fragments` 2
17% of the reference implementation-only peptides have a precursor hypothesis of ours but no spectrum reached the search. The gate that
drops a hypothesised precursor at assembly is `assembly:min_fragments` (a spectrum with fewer correlated fragments is
not emitted). Arm: min_fragments 2 (the default is already 3), both engines. Prediction: the dt-only "we have SOME PSM"
fraction rises from 64% and MSFragger gains if those spectra are identifiable; falsifier: the added spectra are
unidentifiable noise (Sage/MSFragger flat or down, emission up). Runs after the backfill rerun (`bench_minfrag.sh`).

### CONTENT CANDIDATE 3 CONFIRMED (2026-09-03, 05:39): raw apex-frame peak backfill lifts MSFragger to 90% -- at a Sage cost
| dataset D arm | Sage @1% (vs apex 12,537) | MSFragger @1% (vs apex 11,927) | ratio vs dt | shared with dt / dt-only | dt-only matched-ion deficit |
|---|---|---|---|---|---|
| apex | 12,537 | 11,927 | 85.6% | 11,037 / 3,937 | +2 |
| **backfill N=50** | 12,208 (−2.6%; 1,215 / 886) | **12,539 (+5.1%)** | **90.0%** | **12,188 / 2,786** | +3 (on the smaller remaining set) |
| backfill N=150 | 11,190 (−10.7%) | 12,042 (+1.0%) | 86.4% | 11,437 / 3,537 | +2 |
The untraced peaks of the precursor's apex frame ARE the missing content: adding the 50 most intense of them inside the
IM band converts 1,151 the reference implementation-only peptides into shared ones and moves the MSFragger ratio from 85.6% to 90.0% --
the first lever that moves it since the calibration fix. It costs Sage 2.6% (the added peaks are noise for the
spectra that were already rich: Sage's scoring is hurt by them, MSFragger's top-150 cap is not); N=150 is too much
for both. Not a default as measured -- the engines disagree. **Pre-registered next arm:** backfill only the faint
tail: `SPEXTRACTOR_BACKFILL_RAW=50` + `SPEXTRACTOR_BACKFILL_MAXFRAGS=K` (apply only to spectra with <= K assembled
fragments; K = 50 and 150 -- the first launch of this arm collided with a second chain on the same output paths at
06:18 and was stopped and relaunched cleanly at 06:19). Prediction: MSFragger keeps most of +5% while Sage returns to ~12,500; falsifier: Sage still
drops -> the noise cost is in the faint spectra themselves, and the selection must be per-peak (IM/RT co-elution of the
raw peak, i.e. MS2 aggregation across neighbouring frames).

### ASSEMBLY-LOSS ARM FALSIFIED (2026-09-03, 06:16): `min_fragments 2` changes nothing
928,378 spectra (+565), Sage and MSFragger identical to apex, dt-only PSM coverage unchanged (64%). The hypothesised-
but-not-emitted precursors (17% of the reference implementation-only) have NO correlated fragment trace at all, not too few -- the same
untraced-content problem as the searched tail. **Pre-registered follow-up:** `SPEXTRACTOR_BACKFILL_EMPTY=1` with
`SPEXTRACTOR_BACKFILL_RAW=50`: a hypothesised precursor with no correlated trace still gets a spectrum made of its apex
frame's IM-band peaks (min_fragments honoured). Prediction: the dt-only "no PSM at all" count (1,436) falls and
MSFragger rises further; falsifier: those spectra are unidentifiable (counts flat, emission up). Runs after the
faint-tail-only backfill arm.

### FAINT-TAIL-ONLY BACKFILL (2026-09-03, 07:11): removes the Sage cost, but also most of the MSFragger gain
`SPEXTRACTOR_BACKFILL_RAW=50` applied only to spectra with <= K assembled fragments:
| dataset D arm | Sage @1% (vs apex 12,537) | MSFragger @1% (vs 11,927) | dt-only deficit |
|---|---|---|---|
| unconditional N=50 (05:39) | 12,208 (−2.6%) | **12,539 (+5.1%)** | +3 on the remaining set |
| **K=150** | **12,724 (+1.5%; only-apex 11 / only-K150 198)** | 12,000 (+0.6%) | +2 |
| K=50 | 12,591 (+0.4%; 6 / 60) | 11,961 (+0.3%) | +2 |
Reading: the MSFragger gain of the unconditional arm comes from backfilling the RICH spectra (73% of spectra sit at
the 500-cap; share-all gives even faint precursors hundreds of fragments), which is exactly where Sage pays. Gating by
fragment count therefore does not reach the reference implementation-only spectra (deficit still +2). K=150 is a clean small win on
both engines (+1.5% / +0.6%, superset-like) -- a candidate default pending dataset A/dataset B + entrapment, not a mechanism.
**Pre-registered next arm:** gate the backfill by PRECURSOR intensity instead (`SPEXTRACTOR_BACKFILL_MAXQ=q`: only
precursors below the window's q-quantile of MS1 intensity), N=50, q = 0.5 and 0.25. Prediction: MSFragger keeps most
of +5% (the faint precursors are the dt-only ones) while Sage stays >= 12,500; falsifier: the gain tracks the rich
spectra regardless of precursor intensity (then MSFragger simply likes more peaks, and the two engines want
different spectra -- a per-engine output option, not a default).

### BACKFILL-EMPTY FALSIFIED (2026-09-03, 08:35): 53,410 extra spectra buy nothing on either engine
`SPEXTRACTOR_BACKFILL_EMPTY=1` (a hypothesised precursor with no correlated fragment trace still gets a spectrum built
from its apex frame's IM-band peaks) + `SPEXTRACTOR_BACKFILL_RAW=50`. Emission rises 927,813 -> **981,223 (+5.8%)**, so
the feature does what it says; the identifications do not follow.
| dataset D arm | spectra | Sage @1% (apex 12,537) | MSFragger @1% (apex 11,927; dt 13,932) | dt-only |
|---|---|---|---|---|
| empty + K=0 (backfill all) | 981,223 | 12,218 (−2.5%) | 12,600 (+5.6%, 90.4% of dt) | 3,309; some PSM 58% |
| empty + K=60 (faint tail) | 981,223 | 12,593 (+0.4%) | 11,982 (+0.5%) | 3,892; some PSM 64% |
| no-empty N=50 (reference) | 927,813 | 12,208 (−2.6%) | 12,539 (+5.1%) | some PSM 64% |
| no-empty K=150 (reference) | 927,813 | **12,724 (+1.5%)** | 12,000 (+0.6%) | +2 |
Reading: at equal backfill setting the empty arm reproduces the no-empty arm to within replicate noise on both engines
(12,600 vs 12,539 MSFragger; −2.5% vs −2.6% Sage), and the K-gated empty arm is *worse* than K=150 alone (12,593 vs
12,724 Sage). The 53k trace-free precursors are not identifiable: the pre-registered falsifier ("counts flat, emission
up") is met exactly. **Verdict: `SPEXTRACTOR_BACKFILL_EMPTY` stays off and is not a default candidate.** It also
confirms the assembly-loss finding from the other side -- the reference implementation-only precursors we never emit are not lost to
a fragment-count threshold, they have no correlated MS2 signal for us to find at all.

### PRECURSOR-INTENSITY-GATED BACKFILL (2026-09-03, 09:05): prediction held; q=0.5 is the best point on the frontier
`SPEXTRACTOR_BACKFILL_RAW=50` applied only to precursors below the window's q-quantile of MS1 intensity. Emission is
unchanged (927,813 spectra in every arm), so this is purely a content lever.
| dataset D arm | Sage @1% (apex 12,537) | MSFragger @1% (apex 11,927; dt 13,932) | dt-only | shared w/ dt |
|---|---|---|---|---|
| **q=0.5** | 12,442 (−0.8%; only-apex 460 / only-q 365) | **12,492 (+4.7%, 89.7% of dt)** | 2,940 | **12,034** |
| q=0.25 | 12,504 (−0.3%; 141 / 108) | 12,008 (+0.7%) | 3,913 | 11,061 |
| N=50 unconditional | 12,208 (−2.6%) | 12,539 (+5.1%) | — | — |
| K=150 fragment-gated | **12,724 (+1.5%)** | 12,000 (+0.6%) | — | — |
Reading: the pre-registered falsifier is NOT met. Backfilling only the faint HALF of precursors retains 92% of the
unconditional MSFragger gain (+4.7 of +5.1 points) while cutting the Sage cost by two thirds (−0.8% vs −2.6%), so the
gain does track precursor intensity rather than peak count alone. It is concentrated in the second quartile band:
q=0.25 (the faintest quarter only) collapses to +0.7%, i.e. the very faintest precursors' apex frames carry nothing
identifiable -- the same wall the backfill-empty arm hit. Sage degrades monotonically with backfill volume in all
arms. q=0.5 also gives the highest overlap with the reference implementation measured to date (12,034 shared, dt-only down to 2,940),
but the "no PSM at all" residue is 1,406, unmoved from baseline: this is rich-spectrum content, not the deficit
mechanism.
**Two candidate defaults now stand, and they are not the same trade:** K=150 is strictly better than baseline on both
engines but small (+1.5% / +0.6%); q=0.5 is much larger on MSFragger (+4.7%) at a −0.8% Sage cost. Per the both-engines
rule neither is adopted on one file. **Pre-registered confirmation arm:** both settings on dataset A and dataset B, both engines,
plus entrapment FDR. Adopt as default only if the dataset D sign holds on both files AND entrapment stays <= 1.4%
(shipping arm: dataset D 1.37%, dataset A 0.95%, dataset B 1.32%); a Sage regression on dataset A/dataset B larger tha dataset D's −0.8% kills q=0.5 and
leaves K=150 as the only candidate.

### ISOTOPE DUPLICATION: HYPOTHESIS (a) CONFIRMED, (b) REFUTED (2026-09-03, 09:35)
Full analysis in [ISOTOPE-DUPLICATION-2026-09-03.md](ISOTOPE-DUPLICATION-2026-09-03.md). 31.4% of
emitted spectra have a co-eluting same-charge partner 1-3 isotope steps below (decoy offset 3.8%,
excess 27.7%); removing them leaves 671,270 spectra against the reference implementation's 700,434, i.e. isotope-offset
duplication accounts for the entire emission excess. But they are NOT mergeable: content cosine
0.505 vs 0.446 for an arbitrary co-eluting neighbour (only 6.5% above 0.9), and every id-agnostic
collapse rule loses 5,011-5,173 of 14,944 peptides (oracle still loses 1,387). "Keep the lightest" is
backwards -- where only the heavy member is identified, Sage says the mass we reported was already
correct 70% of the time. Mechanism: `findPartner` skips `used[]` peaks, so a leftover one step above a
consumed run becomes its own monoisotope (heavy members carry fewer isotope partners, 3.53 vs 4.59).
Upstream ownership gate pre-registered; downstream merging stays falsified.

### BACKFILL CONFIRMATION ARM: BOTH CANDIDATES FALSIFIED (2026-09-03, 11:47) -- backfill stays OFF
Pre-registered arms `q=0.5` and `K=150` re-run on dataset A and dataset B, both engines, plus entrapment.
| file | arm | Sage @1% | vs apex | MSFragger @1% | vs apex |
|---|---|---|---|---|---|
| dataset D | apex | 12,537 | -- | 11,927 | -- |
| dataset D | q=0.5 | 12,442 | −0.8% | 12,492 | +4.7% |
| dataset D | K=150 | 12,724 | **+1.5%** | 12,000 | +0.6% |
| dataset A | apex | 10,789 | -- | 11,444 | -- |
| dataset A | q=0.5 | 10,516 | **−2.5%** | 11,620 | +1.5% |
| dataset A | K=150 | 10,618 | **−1.6%** | 11,483 | +0.3% |
| dataset B | apex | 10,156 | -- | 9,738 | -- |
| dataset B | q=0.5 | 10,134 | −0.2% | 10,178 | +4.5% |
| dataset B | K=150 | 10,166 | +0.1% | 9,751 | +0.1% |
Entrapment is flat everywhere (dataset A apex 0.91% / q05 0.91% / K150 0.97%; dataset B 1.13% / 1.13% / 1.14%),
so neither arm is bought with false positives -- but neither survives its own falsifier:
* **K=150 dies on the sign rule.** dataset D said +1.5% on Sage; dataset A says −1.6%. The "strictly better on both
  engines" claim was a one-file artefact.
* **q=0.5 dies on the magnitude rule.** The pre-registration killed it if any dataset A/dataset B Sage regression
  exceeded dataset D's −0.8%; dataset A is −2.5%. The MSFragger gain does replicate (+4.7 / +1.5 / +4.5), so the
  two engines genuinely want different spectra -- but that makes it a per-engine option at best, not a
  default, and nothing in the evidence says which engine to optimise.
**Verdict: raw apex-frame backfill stays OFF by default in every form tested** (unconditional,
fragment-gated, precursor-intensity-gated, empty-precursor). The MSFragger deficit is not closed by
adding peaks to spectra we already emit.

### NEW DEFAULT: `charge:min_charge` = 2 (2026-09-03, 15:05) -- confirmed on all three files, both engines
Singly-charged precursor hypotheses are ~30% of emission and ~1.7% of peptides. A z=1 pseudo-spectrum
is identified on 0.42% of its own spectra against 15.3% for z=2; tryptic peptides are essentially
never 1+ in ESI and `charge:scoring count` breaks ties toward the LOW charge, so the mis-assignments
land there. Dropping them removes their share of the multiple-testing burden too, so peptides go UP.
| file | spectra | Sage @1% | MSFragger @1% | entrapment |
|---|---|---|---|---|
| dataset D | 927,813 -> **656,254** (−29.3%) | 12,537 -> **12,642** (+0.8%) | 11,927 -> **12,073** (+1.2%) | 1.28% -> 1.38% |
| dataset A | 1,228,875 -> **863,319** (−29.7%) | 10,789 -> **11,084** (+2.7%) | 11,444 -> **12,061** (+5.4%) | 0.91% -> 0.99% |
| dataset B | 844,755 -> **586,069** (−30.6%) | 10,156 -> **10,294** (+1.4%) | 9,738 -> **9,822** (+0.9%) | 1.13% -> 1.28% |
dataset D wall time 15:21 -> 12:48 (−17%). Every pre-registered adoption condition is met: peptides rise on
BOTH engines on ALL THREE files, every entrapment estimate stays inside the corresponding apex 95% CI
and below the 1.4% bound. Against the reference implementation, MSFragger goes 85.6 -> 86.7% (dataset D), 83.7 -> 88.3% (dataset A),
91.8 -> 92.6% (dataset B); Sage on dataset D is 109.8% of the reference implementation. This is the first change in the emission line
that improves both engines while cutting emission AND runtime.
**Adopted as the default.** `charge:min_charge=1` restores the old behaviour; 3+ is catastrophic
(charge 2 carries 56.9% of peptides: dataset D Sage 5,088). The isotope collapse is NOT adopted in any form
-- on top of this gate it turns +0.8% into −3.8% (Sage) and +1.2% into −2.5% (MSFragger).


### WINDOW-LOOP PERFORMANCE LINE (2026-09-03, evening) -- one adopted, one held
| dataset D, 100 threads | wall | peak RSS | window-loop occupancy | Sage @1% |
|---|---|---|---|---|
| morning default (charge gate) | 12:53 | 103 GB | 64.8x | 12,642 |
| task pool alone | 15:50 | 152 GB | 86.2x | 12,642 (set identical) |
| **+ scoring gate reads a parallel array** | **7:17** | 164 GB | 80.7x | **12,642 (set identical)** |
| + integer detector (`trace:detector=integer`) | **6:53** | **105 GB** | 67.6x | 12,466 (−1.4%) |
**ADOPTED: the parallel-array gate and the task pool.** The gate rejected 99.4% of the fragments it
visited (867 billion visits for 5.5 billion survivors, measured) and read the field it needed out of
a 96 B record: ~83 TB of memory traffic per run. Set-identical output, replicated three times. The
task pool alone was a regression only because it amplified that traffic; with the gate fixed it is
neutral-to-positive, and the in-flight cap added on the earlier diagnosis is removed again.
**NOT ADOPTED (yet): the integer detector.** Faster and 36% lighter, but −1.4% Sage peptides against
the reference, with a symmetric ~8% set churn that is the residual of a faithful-not-identical
reimplementation. The only semantic difference between its 12,650 and 12,466 versions is whether a
frame with only sub-noise peaks counts as empty (the OpenMS rule) or as a miss; that A/B, MSFragger,
and entrapment are running. Adoption needs: peptides within the replicate spread on BOTH engines,
entrapment inside the charge-gate arm's interval, and dataset A replication of the runtime.


### WINDOW-LOOP PERFORMANCE, ROUND 2 (2026-09-04) -- three output-neutral changes, THREE files

Measured old-vs-new with both arms on the SAME host, three files concurrently on three machines,
all running one binary from a single shared install so the
arms cannot differ by build. **Every arm is digest-identical to its base**, which is the gate:
these changes are output-neutral by construction, so peptides are identical by definition and no
search was run.

| | dataset D (128 cores) | dataset A (224 cores) | dataset B (224 cores) |
|---|---|---|---|
| wall base -> perf (SINGLE PAIR -- see the retraction below) | 6:31 -> 6:07 | 8:36 -> 8:09 | 5:40 -> 5:22 |
| process peak RSS | 85.7 -> 83.3 GB (−2.8%) | 119.6 -> 114.2 GB (−4.5%) | 91.3 -> 88.6 GB (−2.9%) |
| **RSS at end of window loop** | 65.0 -> **48.3 GB** (−25.8%) | 85.6 -> **55.1 GB** (−35.6%) | 70.0 -> **49.8 GB** (−28.9%) |
| window loop wall | 209.5 -> 192.7 s (−8.0%) | 317.6 -> 303.5 s (−4.4%) | 175.6 -> 162.7 s (−7.3%) |
| window-loop occupancy | 67.6 -> 69.9x | 56.2 -> 61.3x | 56.9 -> 63.3x |
| system time | 7:12 -> 4:55 (−32%) | 11:34 -> 8:08 (−30%) | 8:17 -> 6:05 (−27%) |
| EPD MassTrace (ledger) | 7,502 -> 0.08 MB | 10,578 -> 0.08 MB | 9,192 -> 0.08 MB |
| sum of per-structure peaks | 36.2 -> 28.3 GB | 58.5 -> 46.9 GB | 41.3 -> 31.2 GB |
| digest vs base | **IDENTICAL** | **IDENTICAL** | **IDENTICAL** |

### RETRACTION (same day): the wall-clock claim does NOT survive replication

The single-pair table above was replicated as an interleaved **base/clean/base/clean** on dataset B,
on one 224-core node, with the load recorded before each run. Four runs, all four digest-identical:

| run | load before | wall | window loop | window-loop CPU | RSS at end of loop |
|---|---|---|---|---|---|
| base rep1 | 3.4 | 5:29 | 174.3 s | 9,850.8 | 69,999 MB |
| clean rep1 | 29.5 | 5:19 | 159.8 s | 10,269.2 | 48,215 MB |
| base rep2 | 41.8 | 5:25 | 168.3 s | 9,988.9 | 69,536 MB |
| clean rep2 | 44.7 | 5:36 | 173.8 s | 10,400.5 | 50,135 MB |

**Wall clock: base mean 5:27, clean mean 5:27.5. There is no wall-clock difference.** The within-arm
spread (base 5:25-5:29, clean 5:19-5:36) is larger than the −5% the single pairs appeared to show.
The window loop is 171.3 s vs 166.8 s (−2.6%) with overlapping ranges. **Window-loop CPU-seconds are
consistently ~4% HIGHER for the new code** (9,850/9,988 vs 10,269/10,400), which is the one timing
signal that reproduces, and it points the wrong way: per-trace `detectPeaks` calls and the gate's
extra build pass cost CPU that batching did not.

**What DOES replicate, tightly, is the memory.** RSS at the end of the window loop: 69,999 / 69,536
(base) vs 48,215 / 50,135 (clean) -- **−20.6 GB, −29.5%**, with the two reps of each arm within
2 GB of each other. That matches the −25.8% / −35.6% / −28.9% seen on the first pass and is
corroborated structurally by the EPD ledger lines falling 7.5-10.6 GB -> 0.08 MB, which is a change
in what exists, not a measurement.

**So the defensible claim is: identical output, ~30% less memory in the window loop, no measured
wall-clock change, and possibly ~4% more CPU.** The earlier "−5 to −6% wall on three files" was
three single pairs taken while another user's load drifted underneath them, and it is withdrawn.

#### THREE-PAIR REPLICATION, dataset D on the 224-core node (the settled numbers)

Six runs, interleaved base/clean x3, load recorded before each, **all five digest comparisons
identical**. Every metric below has NON-OVERLAPPING ranges between the arms, which is the bar the
single pairs failed:

| metric | base mean [range] | clean mean [range] | delta |
|---|---|---|---|
| wall | 321.3 s [320-323] | 312.7 s [310-317] | **−2.7%** |
| window-loop CPU-s | 9,228.9 [9,182-9,309] | 9,689.7 [9,624-9,802] | **+5.0%** |
| RSS at end of loop | 63,880 MB [63,715-64,137] | 46,922 MB [46,459-47,313] | **−26.5%** |
| process peak RSS | 81.3 GB [80.9-81.9] | 77.9 GB [76.4-78.9] | **−4.1%** |

**Final reading, both files together.** The memory result is large, tight and certain: ~27-30% off
the window loop on both files, ~4% off process peak, arms never overlapping. The wall-clock result
is a REAL BUT SMALL win on dataset D (−2.7%, half what the single pairs claimed) and ZERO on dataset B --
so "up to ~3% on one file, none on another", not "−5 to −6% on three files". The CPU cost is real
and reproduces on both files (+5.0% dataset D, +4.2% dataset B): per-trace `detectPeaks` calls and the gate's
extra build pass cost CPU that batching did not. **The trade is ~27% of the window loop's memory
for ~5% more CPU at roughly unchanged wall time.** For a tool whose concurrency is bounded by the
free-RAM admission gate, that is worth taking -- but it is a memory change and must be cited as one.

**THE MACHINES ARE SHARED.** The peer nodes were called "idle" on the strength of `nproc` and
`free` alone -- load average was never checked. They were not idle: the 128-core node was carrying
another user's three python jobs at load 42, the 224-core node three `jackhmmer` at load 21. A later unpaired run of the shipping binary came back at dataset D 6:28 / dataset A 8:40, i.e. back at
baseline, purely because it ran in a more contended window than the arms it was being compared to.
**Rule from here: quote a wall-clock delta only from arms run BACK TO BACK on one node, check
`uptime` before and after, and prefer the interleaved A/B/A/B below to a single pair.** The table
above satisfies the back-to-back condition; the unpaired 11:02 runs did not and are not quoted.

**Read the memory row that matters.** Process peak RSS barely moves (2-5 GB) because the peak is
set outside the window loop; the window loop itself now ENDS 26-36% lighter. That is the number
that governs how many windows can be in flight, not the process peak.

**CPU-seconds RISE slightly on two of three files** (dataset A 5:40:40 -> 5:51:55, dataset B 3:21:01 ->
3:23:16) while wall falls everywhere. That is not a contradiction and it is not free: occupancy
rose 2-6x, so the tool is using more of the machine for less elapsed time. Do not quote the CPU
column as an improvement.

The three changes:
1. **Arena reserves.** `TraceStore::absorb` appended without reserving into a destination that had
   just been swapped empty, so 12 band merges and 48 chunk merges regrew geometrically.
2. **Streamed valley splitting.** Was: build every trace of a chunk as an OpenMS `MassTrace`, then
   split the whole chunk, so the input payload and the entire split output were alive together.
   Now one trace is in flight at a time. `ElutionPeakDetection` is hoisted to chunk scope -- doing
   its four `setValue` calls per trace would have been ~4e4 `Param` round-trips per chunk and would
   have turned this into a slowdown.
3. **Quantised fragment RT gate.** The gate rejects 99.4% of what it visits and was reading an 8 B
   double to do it: 867e9 visits = 6.9 TB per dataset D run. A `uint16` bucket array makes the reject
   2 B; the exact test is unchanged and only runs on survivors. The score-gate counters are
   byte-identical between arms (862,781,140,969 visits, 5,467,207,106 survivors), which is the
   gate-level proof the candidate set did not move.

**Review verdict: OUTPUT-NEUTRAL, with two real defects in change 3 that the first safety proof had
missed.** That proof assumed finite intermediates; a denormal `delta_rt` gives `1/delta_rt = inf`, fragment
RTs spanning +/-DBL_MAX give an infinite span, and a NaN `pc.rt` reaches `(int)NaN` -- undefined
behaviour that on x86 rejects every fragment where the old test passed all of them. Fixed by making
`bucketOf` NaN-safe and degrading the whole field to disabled (`inv = 0` => every bucket 0 => every
`dq` 0 => the exact test alone decides), which costs no branch in the scan. Also fixed: the
`MEM_EPD_MT` guard was constructed after the loop and measured an already-destroyed payload, and
`absorb()` left `bins` short if a binless source followed a binful one.

**REJECTED on evidence, not opinion:**
* **Slab compaction to above-noise peaks** (was the largest proposed memory item). The existing
  `[mem] window ... seeds N of P peaks` line reads **100.0%** -- every peak is above noise, so
  there is nothing to drop and the slab's ~12 GB is irreducible. Killed for free, from a number the
  binary already printed.
* **Moving the seed sort into the band tasks.** Priced by its own supporting evidence at <1% of CPU,
  and it re-partitions seeds across bands -- the exact failure mode that cost 3.3% of peptides once.


## v0.3.0 FINAL BENCHMARK — six datasets, shipped defaults (2026-09-04)

Every dataset in the cohort, run with `spextractor -in <file>.d -out pseudo.mzML -threads 100` and
**nothing else**. If a figure below needed a flag, the defaults would be wrong; that is now a test.
Three nodes in parallel, both engines, `spx:detector=integer` and
`spx:require_isotope_support=1` recorded in every output.

| file | spectra | wall | peak RSS | window-loop occupancy | Sage @1% | MSFragger @1% |
|---|---|---|---|---|---|---|
| dataset A | 862,716 | 8:26 | 108.7 GB | 62.8x | 10,909 | 11,935 |
| dataset B | 585,503 | 5:37 | 84.6 GB | 61.3x | 10,272 | 9,691 |
| dataset C | 723,314 | 6:27 | 90.3 GB | 65.5x | 12,149 | 12,516 |
| dataset D | 655,776 | 5:21 | 79.0 GB | 60.7x | 12,482 | 12,337 |
| dataset E | 542,533 | 5:38 | 67.9 GB | 61.4x | 11,217 | 10,585 |
| dataset F | 597,267 | 5:59 | 73.8 GB | 66.7x | 11,362 | 11,049 |

**dataset C and dataset F had never been benchmarked before**; the cohort had only ever been exercised on
dataset A/dataset B/dataset D. Both behave like the rest, which is the first evidence that the defaults generalise
beyond the three files every decision in this file was made on.

Ranges across the cohort: wall 5:21-8:26, peak RSS 67.9-108.7 GB, occupancy 60.7-66.7x, emission
0.54-0.86 M spectra. Memory tracks acquisition size, and 108.7 GB on dataset A is what sets the
"80-125 GB" requirement now stated in the README.

**Read the two engines separately, as always.** Sage and MSFragger disagree on which files are
easy: dataset C is MSFragger's best (12,516) and Sage's second (12,149), while dataset B is the weakest on
both. Do not average them, and do not quote one as "the" peptide count.

Timings are NOT comparable between rows: the three nodes differ (128 vs 224 cores) and are shared
with other users, whose load is recorded in each run's own output file. The numbers that ARE
comparable across rows are spectra, peptides and peak RSS.


## THE READER QUESTION, SETTLED (2026-09-05)

`perf:stream_load` true/false disagree on ~0.09% of spectra and neither arm had ever been searched.
Both arms, both engines, two datasets, one node each, back to back -- plus entrapment on both.

| | dataset D | | dataset A | |
|---|---|---|---|---|
| | stream=true | stream=false | stream=true | stream=false |
| spectra | 655,776 | 656,371 | 862,716 | 865,119 |
| wall | 4:58 | 9:59 | 8:33 | 13:57 |
| peak RSS | 80.0 GB | 122.2 GB | 108.1 GB | 188.8 GB |
| Sage @1% | 12,482 | 12,485 | **10,909** | **11,132** |
| MSFragger @1% | 12,337 | 12,335 | 11,935 | 11,953 |
| peptide-set symmetric difference | 0.25% of union | | 4.24% of union | |
| **entrapment FDR** | **1.26% [1.01-1.51]** | **1.26% [1.02-1.49]** | **1.13% [0.89-1.40]** | **1.15% [0.93-1.40]** |

**Verdict: the one-shot reader is slightly BETTER, and its extra peptides are real.** On dataset D
the two are equivalent (3 peptides on Sage, 2 on MSFragger, 0.25% set churn -- noise). On dataset A,
the larger file, the one-shot reader finds **+223 Sage peptides (+2.0%)** and +18 MSFragger, and
**entrapment FDR is unchanged** (1.13% vs 1.15%, intervals essentially identical). The gain is not
manufactured false discovery; it is signal the streaming path loses.

Mechanism (inferred, not yet tested): the streaming reader peak-picks each frame on arrival, so it
has less context than a picker running over the whole loaded run. Dataset A is the largest
acquisition in the cohort, and it is where the loss appears -- consistent with something being lost
at frame-batch boundaries.

**The default stays `true`, but as an explicit TRADE, not a free choice**: streaming costs up to ~2%
of Sage peptides on a large file and buys 1.6-1.75x less memory (108 vs 189 GB on dataset A) and
~1.7x less wall time. On a machine that can hold it, `-perf:stream_load false` is the more sensitive
setting. That is now a documented user choice rather than an unexamined default.

## ENTRAPMENT ON THE SHIPPING CONFIGURATION, ALL SIX DATASETS (2026-09-05)

Every peptide count in this file rests on the engines' own 1% FDR; this is the check on it.
Sage entrapment search (human targets + Arabidopsis entrapment, peptide-hypothesis ratio 0.6805).

| dataset | target | entrapment hits | raw% | **FDR%** | 95% CI |
|---|---|---|---|---|---|
| A | 10,273 | 79 | 0.76% | **1.13%** | 0.89-1.40 |
| B | 9,859 | 86 | 0.86% | **1.28%** | 1.03-1.58 |
| C | 11,781 | 93 | 0.78% | **1.16%** | 0.92-1.42 |
| D | 11,927 | 102 | 0.85% | **1.26%** | 1.01-1.51 |
| E | 10,710 | 90 | 0.83% | **1.23%** | 1.00-1.48 |
| F | 10,869 | 97 | 0.88% | **1.31%** | 1.04-1.57 |

**True FDR is 1.13-1.31% at a nominal 1%** -- a mild, CONSISTENT under-estimate by the engine, with
every interval overlapping every other. Nothing anomalous, no dataset out of family, and the two
never-before-benchmarked datasets (C and F) sit inside the same band. The peptide counts in this
file are therefore honest to within roughly a quarter of a percentage point of FDR.


## PUBLIC BENCHMARK: PXD017703 (HeLa diaPASEF) RUNS END TO END (2026-09-05)

The standard public diaPASEF reference (Meier et al., Nat Methods 2020) is staged under
`benchmark/PXD017703/raw`: 9 Evosep acquisitions (60/100/200 SPD, 200 ng HeLa) and 18 "highsens"
acquisitions, 27 in all. PRIDE's advertised byte sizes were wrong for BOTH archives; `unzip -t`
CRC is the arbiter, size is only a hint, and a failed CRC resumes rather than unzips.

Making the tool run on it exposed three defects the in-house cohort could never show:

1. **"exactly one MzCalibration row" refused valid data.** These files carry two rows; every frame
   references Id 1. Three places (the OpenMS patch's converter, `loadTofAxis`, the mzPeak loader)
   now select the row `Frames.MzCalibration` references and refuse only if frames reference more
   than one. The OpenMS patch is regenerated and re-verified byte for byte against the pinned base.
2. **`C2 == 0` was rejected as "missing".** The 2020 timsTOF Pro firmware behind these files ships a
   purely linear-in-sqrt calibration with no quadratic term, stored as a real 0.0. The guard existed
   because a NULL C2 reads as 0.0 through the C API and would silently select the wrong law -- so
   the loader now converts NULL to NaN (refused with its own reason) and a stored 0.0 is accepted.
   The cancellation-free root already collapses to (t - C0)/b at C2 = 0; the golden tests pin it.
3. **The shipped detector never ran in a public build.** `loadTofAxis` sat behind
   `#if __has_include(<SQLiteCpp/SQLiteCpp.h>)`, false for an OpenMS INSTALLATION, so the integer
   detector silently fell back to `openms` and only one log line said so. Caught by the tenth
   end-to-end check on the stock-OpenMS CI. Both readers now use the sqlite3 C API through one
   helper; a missing SQLite is a compile error. Same class as the two default-flag findings.

First extraction, 100SPD 100 ng, shipped defaults: **136,818 pseudo-MS2 spectra, 57 s wall, 28 CPU-min,
16 GB peak**, `spx:mz_calibration=tdf_table_modeltype1`, `spx:detector=integer`.

Two operational lessons, recorded because each cost an hour: the shared install carries its own
`lib/libOpenMS.so`, so a converter fix in the OpenMS patch is invisible until the LIBRARY is
republished, not just the executable; and the reference implementation's JVM aborts with
`locale::facet::_S_create_c_locale name not valid` under a non-interactive environment -- run it
with `LC_ALL=C`.

**On reproducing the reference implementation's paper.** The paper (Nat Commun 2025, PMC11696033) contains NO HeLa
and no PXD017703. Its headline dataset is TNBC, PXD047793: 16 diaPASEF runs on timsTOF Pro,
FragPipe + the reference implementation direct DIA, **9,296 proteins per run on average, 10,341 total** (DIA-NN
quantification, 1% FDR at global protein, global precursor and run-specific precursor). Its
parameters for the reference implementation -- Delta Apex IM 0.01, Delta Apex RT 3, RF max 500, Corr threshold 0.3 -- are
exactly the ones this project's reference chain has always passed. So "numbers similar to the
paper" is a test on TNBC, being staged now (`benchmark/PXD047793`, ~306 GB; archive #001's
`analysis.tdf_bin` arrived with a bad CRC at full length and is re-fetched from scratch after the
rest). Their own converted mzMLs and FragPipe results (MassIVE MSV000094803) would allow a
like-for-like comparison; MassIVE was unreachable from this network on every route tried
(FTP timed out on both hosts, PROXI file listing 404) -- retry later or via the web UI.

### 2026-09-05, 15:30 -- calibration port is output-neutral; the C4 term is modelled; first public-data numbers

**Digest arm (back-to-back with the shipping build):** the calibration-row selection
(`Frames.MzCalibration`, not "exactly one row") plus the sqlite3 C-API port produce a
**byte-identical spectrum digest** on the in-house 30-min file: `ca609dc4…` for both
`final_S30/pseudo.mzML` and `digest_port/pseudo.mzML`, 655,776 spectra, 05:30 wall, 81.8 GB peak.
Output-neutral by construction and by measurement.

**C4 (18 of 27 public HeLa acquisitions refused as "C3/C4 != 0").** The 2019 timsTOF Pro set
(firmware 6.0.110) stores `C4 = -0.0905`, every frame referencing that row. No open implementation
models even C2 (timsrust and opentims are linear-in-sqrt), so the term was identified against the
vendor library as a local oracle: 6,000 TOF indices through `tims_index_to_mz`, C0/b/C2 held fixed,
residual regressed on one candidate correction at a time. Only `1/sqrt(m)` closes it (11.8 ns ->
0.004 ns; every other candidate stays at 5-12 ns), with coefficient `C4*b/2` -- the first-order image
of an **additive offset on the quadratic root, `m = u^2 - C4`**. Closed form vs vendor: +0.0723 ppm
constant on frame 1 without the temperature term, **0.0000 ppm with it** (the constant IS that frame's
`dC1*dT/1e6 = -7.2e-8`), 0.0000 ppm on the last frame. Dropping the term: **-53 .. -938 ppm**, worst
at low mass (0.0905 Da / 96 = 943 ppm, matches). `C4 == 0` is bit-identical to the old code
(`u*u - 0.0`), so nothing measured above changes. `C3` stays refused: never observed non-zero.
Pinned in both golden tests (10 vendor probes, 2 frames, ablation guard > 50 ppm). Build
pending; the 18 files are then extracted by the idempotent HeLa driver.

**Public HeLa (PXD017703), 2020 Evosep sets, shipped defaults, 100 threads, one node each, Sage 1%
peptides / MSFragger PSM rows (not peptides -- the driver counts rows; peptide/protein counts come with
the reference comparison):**

| acquisition | spectra | wall | RSS | Sage peptides | MSFragger PSM rows |
|---|---|---|---|---|---|
| 100SPD 100 ng S2-C1 (2731) | 136,818 | 0:42 | 16.2 GB | 12,861 | 115,396 |
| 100SPD 100 ng S2-C5 (2735) | 130,459 | 0:41 | 15.3 GB | 12,520 | 109,406 |
| 200SPD 200 ng S3-A1 (2737) | 75,462 | 0:30 | 10.8 GB | 6,154 | 63,161 |
| 200SPD 200 ng S3-A2 (2738) | 74,874 | 0:31 | 10.7 GB | 5,701 | 61,750 |
| 200SPD 200 ng S3-A4 (2740) | 72,238 | 0:27 | 10.2 GB | 6,364 | 60,366 |

2732 was lost to a SIGBUS at the minute the shared install was republished;
re-extracted by hand: 131,773 spectra, 0:43. The reference extractor ran on all 27 (exit 0,
0.4-4.9 GB mzML each, 1:02-1:25 wall / 34-57 GB on the Evosep files, 2:17-6:18 / 67-197 GB on the
2019 files) but the driver tested the output size in the same second the JVM exited and CephFS had
not surfaced it: every line says NO OUTPUT and no reference search ran. Driver is now idempotent
(skips extraction whose output exists, polls the size, skips completed searches); the rerun is queued
behind each node's chain and will fill the reference column.

### 2026-09-05, 15:45 -- the 2019 acquisitions extract: two defects, one per diaPASEF scheme

With C4 modelled, the 18 refused files split by scheme. **25pc** (4 window groups x 4 windows, each
group acquired four frames in a row, 16 distinct m/z): extracts on the C4 build as is -- 10 ng file
19: **220,291 spectra, 6:49 wall, 51 GB** (68,469 frames, a full 2-h gradient at 10 ng). **py3** (16
groups x 4 windows, **32 distinct m/z for 64 window rows**): every m/z window is acquired in TWO groups
with shifted, overlapping scan ranges (412.5: scans 649-918 in group 1, 770-918 in group 9) ~1.7 s
apart in the cycle. Keyed by isolation m/z alone the two slices landed in one window as all of
group A then all of group B; the retention-time order guard fired at frame index 4150 =
66,394 / 16 = frames per group -- the guard did its job. Merging by RT would be wrong: outside the
overlap a precursor sees every other frame empty. The window key is now (lo, hi, WindowGroup)
from the vendor native ID; ScanNumBegin was tried first and collapsed the eight highest windows,
whose two slices both start at scan 0 (56 windows, guard fired at window 48 -- predicted, observed).
For every other scheme each m/z sits in one group, so partition and order are unchanged; the dataset D
digest cannot move. e2e check 11 synthesises the scheme. **py3 10 ng file 22: 160,416 spectra,
1:44 wall, 26 GB, 64 windows.** Published 15:42 under the no-running-spextractor gate.

Operational: the HeLa driver is now idempotent and runs on three nodes (9, 11 and 7 files, the last
being the heavy 2019 acquisitions); one node's first pass ended without its completion marker because
the driver was edited in place while it ran, which is why the reruns are process-gated.

### 2026-09-05, 18:25 -- reference-side Sage was silently zero on every public file

Every reference line of the HeLa driver reported `Sage: 0` while MSFragger returned 250-350k PSM
rows. Cause: the reference extractor's mzML omits the standard "selected ion m/z" cvParam
(MS:1000744) and writes only the isolation-window target; Sage 0.14.6 then resolves the precursor
through `spectrumRef`, finds no MS1 spectrum in the file, and panics ("missing MS1 precursor for
frame=…"), leaving an empty peptide list that the driver counted as zero. The in-house head-to-head
had the same problem and patched it (a precursor-patching script that is not part of this
repository: one streaming pass injecting the value from the isolation window; 425,061 precursors on a
100 ng file). A separate reference-Sage pass with that patch ran all 27 files on one node in 30 min. Reference Sage peptides at
1%: 100SPD 100 ng 15,786 / 10,302 / 14,995; 200SPD 200 ng 7,130 / 4,882 / 7,323; 60SPD 200 ng
22,745 / 22,521 / 15,693; 25pc 100 ng 13,591 / 13,050 / 13,283; py3 100 ng 20,106 / 28,078 / 20,329.
The MSFragger column is PSM rows until the target-decoy scorer is run over the copied tables.

### 2026-09-05, 18:45 -- PUBLIC HeLa (PXD017703), all 27 acquisitions, both extractors, both engines

Full table in `docs/HELA-TABLE.md` (the table script is not part of this repository). Shipped
defaults, `-threads 100`, one node per chain; Sage peptides at 1% (rank 1, peptide q); MSFragger
peptides at 1% by peptide-level target-decoy on the raw tables (best hyperscore per peptide+mods,
decoy = all proteins `rev_`). The reference extractor's Sage column comes from the patched pass.

| scheme | n | Sage SpeX | Sage ref | ratio | MSFragger SpeX | MSFragger ref | ratio |
|---|---|---|---|---|---|---|---|
| 2019 25pc | 9 | 110,211 | 90,807 | **1.21** | 83,189 | 101,837 | 0.82 |
| 2019 py3 | 9 | 135,107 | 159,105 | 0.85 | 98,956 | 184,824 | **0.54** |
| 2020 100SPD | 3 | 38,213 | 41,473 | 0.92 | 34,964 | 51,143 | 0.68 |
| 2020 200SPD | 3 | 18,394 | 19,513 | 0.94 | 17,389 | 26,894 | 0.65 |
| 2020 60SPD | 3 | 62,678 | 61,551 | 1.02 | 58,345 | 75,034 | 0.78 |

What it says, in order of weight:
1. **The MSFragger deficit is real and larger on public data than in-house** (0.54-0.82 vs the
   in-house 0.85-0.90). It is worst on py3 and at 10 ng. Sage is at parity or better on three of
   five schemes; the engine disagreement is therefore NOT "natural variance" here -- MSFragger
   consistently gets ~1.4-1.9x more out of the reference spectra than out of ours, which points at
   spectrum CONTENT (fragment completeness per spectrum), the same conclusion as content candidate 2.
2. **Load dependence.** At 10 ng DIAspeXtract loses on both engines and both schemes (py3 10 ng: Sage
   4.8k vs 9.2k, MSFragger 2.7k vs 8.1k); at 50-100 ng on 25pc it wins with Sage by 21%. The
   isotope-support gate and the z>=2 floor, both wins in-house at 200 ng, are the obvious suspects
   at low load: on file 22 the gate dropped 2.09 M guessed precursors and kept 466 k.
3. **Emission** is 1.7-2.6x the reference's spectrum count at 50-100 ng 2019 (1.1-1.2 M vs the
   reference's ~425 k precursors), consistent with the over-generation finding; at 10 ng py3 it is
   160 k vs 63 k.
4. **Cost.** DIAspeXtract runs in 1/3 to 1/6 of the reference's RSS on every file (2019 100 ng:
   66-114 GB vs 168-197 GB; 10 ng: 27 GB vs 99-147 GB) and is faster on the 2020 files and py3;
   the 25pc 100 ng walls of 33-47 min are from a node running three jobs (uncontended: 13:38).
5. **Reference-side Sage replicates disagree by up to 35%** (2732: 10,394 vs 15,939/15,140; 2738:
   4,928 vs 7,194/7,391; 35: 28,359 vs 20,302/20,520) while its MSFragger replicates agree within
   5%. A deterministic re-search of 2732 is running; if it reproduces, the variance is in the
   reference extraction, not in Sage, and the reference totals carry that noise.

The reference implementation's paper reports no HeLa numbers; its headline dataset is TNBC (PXD047793), staging
now (7 of 16 by 18:20). That is the like-for-like paper comparison.

**Resolution of point 5 (19:00):** an isolated re-search of the 2732 reference file reproduces
Sage's own count exactly (10,302 target peptides both times; the table's 10,394 is the raw-table
recount, a definitional difference). Sage is deterministic; the replicate spread is in the
reference extractor's output for that acquisition, whose MSFragger count is nonetheless in line
with its siblings. The reference Sage totals therefore carry extraction noise the MSFragger totals
do not, which is one more reason the MSFragger column is the one to close.

## PUBLIC BENCHMARK 2: TNBC (PXD047793), the reference implementation's paper's headline data (2026-09-06, 07:10)

All 16 single-shot TNBC diaPASEF runs (72k frames, 130-min gradients, 630 GB staged) through both
extractors and both engines overnight on three nodes (16 h wall, ~30-38 min per file per extractor,
230-300 GB RSS on both sides at 100 threads). Reference-side Sage patched inline (selected ion m/z).
Driver-log numbers, Sage peptides at 1% (the MSFragger column is PSM rows until the peptide-level
table lands -- the peptide-level table script is running):

| | reference | DIAspeXtract |
|---|---|---|
| Sage peptides, 16 runs | 293,065 | **318,252** (ratio 1.086) |
| per run, mean | 18,317 | 19,891 |
| spectra per run | ~1.3-2.4 M PSM rows | 2.06-2.41 M pseudo-spectra |
| wall per run | 18-30 min | 26-38 min |
| RSS per run | 238-307 GB | 204-287 GB |

DIAspeXtract is ahead on Sage on 16 of 16 runs (smallest margin 001: 15,468 vs 14,633; largest
010: 20,352 vs 16,308). The paper's headline is proteins: FragPipe + the reference implementation quantified an
average of **9,296 proteins per run** (Spectronaut directDIA 8,997, DIA-NN library-free 9,520; all
with quantification/MBR, tryptic). Our comparable figure is Sage protein groups at 1% protein q per
run, pending in the table; identified-only, so expect it below a quantified-with-MBR count.

### TNBC FINAL (2026-09-06, 07:30): all 16 runs, peptide-level both engines, Sage protein groups

`docs/TNBC-TABLE.md`. Sage peptides **321,430 vs 295,781 (1.09)**, ahead on 16/16 runs; Sage protein
groups at 1% protein q ahead on 16/16 (per run 3.7-6.0k vs 3.5-5.8k); MSFragger peptides
**334,607 vs 392,897 (0.85)**, behind on 16/16 -- the in-house MSFragger ratio (0.85-0.90)
reproduces exactly on the paper's own data, so the engine split is a property of our spectra, not of
a dataset. Runtime on these 130-min files: DIAspeXtract 26-38 min vs 18-30 min (median 1.4x slower),
RSS 204-287 vs 238-307 GB (0.9x). The paper's 9,296 proteins/run is a quantified-with-MBR FragPipe
figure over the cohort and is not the same statistic as identified protein groups at 1% in one run.

### 2026-09-06: memory work, steps 1-4

**Step 1 (output-changing, the reason for a new baseline).** The MS1 arena compaction walked the
trace vector in CONTAINER order with a monotone write cursor, but the vector is sorted by m/z
(`ms1_traces` canonical sort) AFTER its spans were appended in detection order. Kept spans visited
after the cursor had passed their offset were overwritten before they were copied. Deterministic, so
no thread-count digest ever saw it. `compactUnreferenced()` now sorts the kept spans by offset,
validates bounds and disjointness before moving a byte, and memmoves down; `-diag:selftest_arena`
(e2e check 12) fails on the old algorithm. Measured effect on dataset D: **9,040 of 655,776 spectra
changed (1.4%), every one of them below precursor m/z 500** -- the lowest m/z band, where detection
order and m/z order disagree most, which is the mechanism's own prediction. New baselines: the dataset D
reference digest becomes `79f2a733…`, with a new TNBC 009 reference alongside it. The old dataset D
digest `ca609dc4…` reproduces exactly on the pre-fix binary, which is the control that the gate arms
are wired correctly.

**Steps 2-4 (output-identical, dataset D digest == `79f2a733…` on each).**

| step | what | measured on dataset D |
|---|---|---|
| 2 | instrumentation only | identical; see below |
| 3 | OpenMS `MassTrace` move operations (third patch) | identical; peak RSS 83.1 -> 82.7 GB, after-MS1 48.5 -> 48.1 GB |
| 4 | MS1 state released at the window-loop join; RT-axis reserve; per-frame release on the non-streaming path | building |

**What the instrumentation measures that the profile could only infer** (dataset D; TNBC is 2-3x
on every count): picked MS1 725,116,915 peaks in 1,343 frames = 13.8 GB as a PeakMap; compact store
11,975 MB for 1.256e9 peaks; per window at the band merge 12 arenas of ~330 MB held twice, ~6 M
traces. And the retention mechanism, as a number rather than an assertion: at the MS1 -> window-loop
boundary **33.4 GB of the 48.5 GB resident set is free-but-retained allocator memory** (arena 47.3
GB, mmapped 345 MB). Two thirds of RSS at that point is memory nobody is using.

**Not quotable from these runs:** wall-clock deltas. Two heavy jobs shared the node (the gate arms
ran alongside), so only digests and memory milestones are comparable. The step-3 wall effect is
measured on TNBC, where the MassTrace payloads are 8x larger.

#### COMPACTION-FIX GATE, dataset D (primary), 2026-09-06 14:10 — the fix PASSES

The pre-fix build against the fix (plus its selftest), one node, arms back to back, immutable arm
directories, binary/library/config/FASTA digests recorded. The pre-fix arm reproduced the historical
digest `ca609dc4…` exactly, which is the control that the arms are wired correctly.

| | pre-fix | fixed | ratio |
|---|---|---|---|
| spectra | 655,776 | 656,257 | 1.0007 |
| Sage peptides @1% | 12,605 | 12,594 | **0.9991** |
| MSFragger peptides @1% | 11,642 | 11,640 | **0.9998** |
| entrapment FDR | 1.26% [1.01-1.51] | 1.26% [1.03-1.51] | unchanged |

Set overlap: Sage 12,575 common, 30 only pre-fix, 19 only fixed (99.61% of the union); MSFragger 11,630
common, 12 and 10 (99.81%). Decision rule (>= 1% loss in either engine on D, or A contradicting D,
or entrapment more than 0.3 points above the base estimate): **none triggered**.

Read honestly: correcting 1.4% of the spectra costs 11 Sage peptides and 2 MSFragger peptides, and
gains 19 and 10 respectively. The corrupted XICs were producing spectra that scored about as well as
the correct ones -- the fix is justified by correctness, not by yield, and nothing downstream was
silently depending on the corruption. Dataset A is still running; its role is to contradict D or not.

#### COMPACTION-FIX GATE, dataset A: the trigger fired, and what the investigation found (14:15)

Dataset A moves far more than D: Sage 11,016 -> 11,218 (**+1.83%**), MSFragger 10,920 -> 10,689
(**-2.12%**), set overlap 95.5% and 89.6% of the union (D: 99.6% and 99.8%). The engines disagree in
sign, which is the decision rule's "A contradicts D" trigger, so the FDR walk was swept before going on:

| threshold | Sage pre-fix -> fixed | MSFragger pre-fix -> fixed |
|---|---|---|
| 0.1% | 8,262 -> 8,641 (**+4.6%**) | 9,768 -> 9,485 (-2.9%) |
| 0.5% | 10,081 -> 10,453 (**+3.7%**) | 11,210 -> 11,121 (-0.8%) |
| 1% | 11,016 -> 11,218 (**+1.8%**) | 11,935 -> 11,679 (-2.1%) |
| 2% | 11,915 -> 12,066 (**+1.3%**) | 12,571 -> 12,555 (-0.1%) |
| 5% | 13,671 -> 13,750 (**+0.6%**) | 13,915 -> 13,679 (-1.7%) |

**Sage's gain is monotone in the strictness of the threshold** -- largest where the evidence is
strongest -- which is what repairing corrupted precursor XICs should look like: better correlations
lift the best identifications most. **MSFragger's deficit is non-monotone** (-2.9, -0.8, -2.1, -0.1,
-1.7) with no gradient: scatter around roughly -1.5%, the signature of a ranking/threshold effect in
its hyperscore target-decoy walk rather than a systematic loss. Entrapment FDR is unchanged on both
arms (1.13% -> 1.12%, intervals overlapping).

**Verdict: the fix stands.** D is flat in both engines, A gains in Sage with a
gradient and scatters in MSFragger, entrapment is unchanged everywhere, and the change being measured
is a correctness fix. This is also the engine disagreement the project has recorded before (the sign
flips by file); it is not evidence against the fix. Dataset A being ~20x more affected than D is
worth a note of its own: whatever makes A's low-m/z MS1 traces denser also made it the file where the
corruption cost the most.

#### Steps 3 and 4 measured on dataset D (both output-identical, digest == `79f2a733…`)

| build | peak RSS | after-MS1 RSS | loop-end RSS |
|---|---|---|---|
| step 2 (instrumentation) | 83.11 GB | 48.55 | 48.98 |
| step 3 (MassTrace moves) | 82.68 GB | 48.09 | 48.91 |
| step 4 (MS1 state released at the loop join) | **81.73 GB** | 48.09 | 47.85 |

Cumulative -1.38 GB (-1.7%) on this file; the same steps are measured on TNBC separately, where the
mass-trace payloads are ~8x larger. Wall-clock is NOT quoted from these runs: the node was shared
with the gate arms for part of them.

One instructive detail from step 4's own log line: releasing the MS1 traces, their arena, the
precursors and the m/z index does **not** move RSS at that instant (47,847,848 kB before and after,
allocator free 40.0 GB) -- the bytes go to the per-thread arenas, not to the OS. What the release
buys is the peak: the sort, the optional merge and consolidate passes and the writer now reuse that
memory instead of growing on top of it, which is where the 0.9 GB comes from. That is the retention
mechanism seen from the other side, and it is why "free earlier" and "return to the OS" are two
different levers.

#### Step 8: the trace record 72 -> 56 bytes (output-identical, dataset D digest == `79f2a733…`)

| build | peak RSS on dataset D |
|---|---|
| step 2 (instrumentation) | 83.11 GB |
| step 4 (MS1 state released) | 81.73 GB |
| **step 8 (56-byte record)** | **78.99 GB** |

-2.74 GB from the record alone on this file, -4.12 GB (-5.0%) cumulative. The record is the largest
structure in the tool at scale (925 M parent + 1,217 M child traces on a 2-h acquisition = 143.6 GiB
at 72 bytes if simultaneous), so the same change is worth several times more there; that measurement
is running. The store is now an explicit argument to every span accessor, asserted in a debug build,
and the calibration factor is a store-local frame index -- kept explicit rather than derived from
frame0 + apex, because the unsplit integer path picks its calibration frame in `got` order while
makeSpan resolves equal-intensity apex ties in frame order, and the two disagree on exact ties.
A static_assert pins sizeof(Trace) == 56.

#### TNBC run 009: every step measured on the file class the work exists for

| build | peak RSS | vs step 2 | digest |
|---|---|---|---|
| step 2 (instrumentation) | 281.52 GB | -- | == TNBC reference |
| step 3 (MassTrace moves) | 279.23 GB | -2.29 | == TNBC reference |
| steps 4 + 8 (lifetimes, 56-byte record) | 270.26 GB | **-11.26** | == TNBC reference |
| + `malloc_trim` (experiment) | 206.66 GB | **-74.86** | == TNBC reference |
| steps 4 + 8 + trim as default | running | | |

Every arm is byte-identical to the step-1 baseline. The structural steps together are worth 11.3 GB
here; the allocator trim alone is worth 74.9 GB, which is the honest ranking and the reason it
shipped first. The record change measured -10.7 GB of that 11.3 rather than the 16-24 GB predicted
from the record count: not every window's records are resident at once, so the concurrency factor
applies to the saving as well as to the footprint.

#### Step 6: the admission gate now has bounds it can hold (output-identical)

Dataset D peak RSS: steps 1-4 and 8 give 75.3 GiB -> **plus the trim default and the fixed gate,
55.6 GiB, -26%**, digest == `79f2a733…`. The admission high-water reached 24 of 24 windows, i.e. neither the
window cap nor the recalibrated byte budget throttled anything: the default behaviour is unchanged
and the honest budget does not reject windows spuriously. That is the point -- the gate can now hold
a limit when one is asked for, without imposing one when it is not.

Three defects were in that one gate: `perf:max_concurrent_windows` was computed and logged but never
checked; the byte projection booked 160 B per peak against a sustained 44 B measured across 28
windows; and `availableBytes_()` added the process's own resident set back into its own budget.
Admission moved to the master, before the task is spawned, so both bounds hold by construction and no
worker sits in a spin loop holding a thread it cannot use.

**Verification note.** The suite reported 12/12 on a build that should have run 13 checks: the
deployment point executes its own copy of `test_spextractor.py` and that copy was stale, so the new
check silently did not exist. The deploy step now ships the suite with the source it tests.
A green result from a suite that quietly dropped a test is worse than a red one.

#### CORRECTION (2026-09-06, 19:55): the trim is NOT free -- it trades ~8.5% wall for ~30% memory

An earlier note here said the trim ran "not slower". That came from one unpaired run on a shared
node and was wrong. The interleaved A/B (baseline vs trim build, same node, same input, alternating)
and the build-by-build series agree:

| build | wall | user CPU | window-loop CPU | minor faults | peak |
|---|---|---|---|---|---|
| baseline (instrumented) | 31:16 | 40,671 s | 33,956 s | 363 M | 268.5 GB |
| + MassTrace moves | 31:36 | 43,984 | 37,203 | | 266.3 |
| + lifetimes + 56-byte record | 32:06 | 41,994 | 35,850 | | 257.7 |
| **+ the trim** | **33:59** | 46,662 | **41,092** | | **186.1** |
| interleaved: baseline r1 | 30:41 | 42,116 | 35,738 | 287 M | 264.9 |
| interleaved: trim build r1 | 33:18 | 47,553 | 40,988 | 362 M | 186.2 |

**Almost the whole memory win and almost the whole wall cost are the same change.** The trim returns
~99 GB at the MS1 boundary; the window loop then has to fault those pages back in (minor faults
+26%, 287 M -> 362 M) and touches them cold, so the loop's CPU rises ~15% even though concurrency is
unchanged (both arms admit all 28 windows, parallelism 37.2x vs 36.0x). System time barely moves
(+46 s), so this is cache and TLB warmth, not kernel overhead.

The second trim, at the loop -> write boundary, cannot affect the peak at all -- the peak happens
inside the loop -- so it is free but also pointless unless something after the loop needs the RAM.

The trade is therefore: **-82 GB of peak (-30.7%) for +2:37 of wall (+8.5%)** on a 2-hour
acquisition. Which side of that is right depends on whether memory or time is the binding
constraint, and that is a decision, not a measurement.

#### The 48-byte record, measured (output-identical on both files)

| build | peak RSS (GiB) |
|---|---|
| baseline (instrumented) | 268.5 |
| trim + admission gate + 56-byte record | 183.0-186.2 |
| **48-byte record + the split's reserve fixed** | **177.5** |

Cumulative **-91.0 GiB, -33.9%**, spectrum digest identical throughout.

The record change is confirmed directly rather than inferred: the per-window `[mem]` lines sum to
64,986 MB of trace records before and 55,699 MB after, a 9.1 GB reduction that is exactly 8/56 of
the total, and the peak fell 5.5 GiB of that -- the difference being the ~60% of windows resident at
the peak. Deriving the apex intensity from the arena costs nothing numerically: it is the same float
the field held, and the digest proves it on both files.

**A unit trap worth recording.** `/usr/bin/time -v` reports kilobytes; dividing by 1e6 gives a number
that looks like GB and is 7% smaller than the GiB the tool's own logs print. Comparing one against
the other made a real 5.5 GiB improvement look like no change at all. Every figure in this file is
GiB; when in doubt, compare kB against kB.

### Scheduling: the window loop was bounded by one window's serial chain (2026-09-07)

The window loop's wall time equalled the longest SINGLE window's stage chain, on both datasets --
prep, trace, split and score run sequentially inside a window, so the last window standing runs
its stages with the machine nearly idle:

| dataset | loop wall | longest window's own chain | that window |
|---|---:|---:|---|
| dataset D | 150.7 s | 147.4 s (98%) | 327-469, trace-bound (98.4 s) |
| TNBC 009 | 748.5 s | 741.9 s (99%) | 575-601, **score**-bound (581.0 s) |

On TNBC that score figure was not work, it was a scheduler cliff. libgomp's `GOMP_taskloop` runs
the whole loop inline in the creating thread once `task_count + num_tasks > 64 * nthreads`; at
`grainsize(32)` and 100 threads that is **204,800 precursors**, and exactly two of 28 windows sit
above it. Measured on the same libgomp with a standalone taskloop, 20 us per iteration:

| iterations | clause | tasks | parallelism |
|---:|---|---:|---:|
| 200,000 | `grainsize(32)` | 6,250 | 82.7x |
| 300,000 | `grainsize(32)` | 9,375 | **1.0x** |
| 300,000 | `num_tasks(400)` | 400 | 88.0x |

The two windows above the line scored 417.6 s and 581.0 s; every neighbouring window -- within 16%
on fragment count -- scored 12-28 s. An unclaused taskloop is safe (libgomp picks ~nthreads tasks,
verified to 10e6 iterations at 99.1x); only a data-dependent `grainsize` is dangerous.

#### Measured, byte-identical on both files

| build | dataset D total | TNBC 009 wall | TNBC window loop | TNBC PRECURSOR_INFER | TNBC peak RSS |
|---|---:|---:|---:|---:|---:|
| before | 312.1 s | 33:35 | 1150.0 s (31.2x) | 438.8 s (1.0x) | 181,777 MB |
| **+ bounded scoring pool + findPartner SoA** | **266.2 s** | **18:03** | **496.5 s (66.7x)** | **177.2 s** | 181,902 MB |

**-46% wall on TNBC, -15% on dataset D, at unchanged peak memory and an identical spectrum list on both.**
The scoring fix is a fixed pool of worker tasks pulling 32-precursor chunks from an atomic cursor,
so the task count no longer grows with the window size. `PRECURSOR_INFER` fell because
`findPartner`'s two gate fields (RT, IM) now sit in m/z order beside `sorted_mz` instead of being
chased through 48-byte records: same values, same order, same first match. Its three stages are now
logged -- on TNBC, m/z index 1.0 s, intensity sort 7.5 s, **greedy walk 156.3 s**, so the residue is
the greedy claim loop and not the sorts.

#### `perf:trace_bands` (dataset D) -- output-CHANGING, needs the peptide gate

| bands | window loop | total | CPU (loop) | spectra | digest |
|---:|---:|---:|---:|---:|---|
| 12 (default) | 136.0 s | 266.2 s | 9,236 s | 656,257 | reference |
| 24 | 117.0 s | 249.8 s | 9,201 s | 656,257 | differs |
| 48 | 108.7 s | 238.7 s | 9,311 s | 656,260 | differs |

CPU is flat within 1%, so the gain is purely tail parallelism, not less work -- which contradicts
the older note that banding reduces total trace work superlinearly (that measurement was on the
OpenMS detector path, not the integer one). A control re-run at 12 bands reproduced to within 2%
(4:41.7 vs 4:36.9) with an identical digest.

#### The two serial stages behind the remaining phases (2026-09-07)

Sub-stage timers went in before either change, because the phase milestones could not say where
the time went. On dataset D, `MS1_TRACE` = 75.9 s split as edges 1.0 | **distribute 39.0 (serial)** |
detect 28.2 | gather 0.2 | valley-split 2.4. That explains why raising `perf:ms1_trace_bands`
made the phase *worse* rather than better -- more bands is more of that same serial copy:

| `perf:ms1_trace_bands` | MS1_TRACE | parallelism | digest |
|---:|---:|---:|---|
| 12 (default) | 75.0 s | 3.9x | reference |
| 24 | 93.7 s | 3.5x | differs |
| 48 | 147.9 s | 2.5x | differs |

Chunking the distribution over spectra took it to 0.55 s. Segmenting the seed list per band took
per-window `prep` on TNBC from 60.1 s to 27.1 s at the worst window (904 -> 393 s summed).

#### Cumulative, byte-identical throughout (dataset D and TNBC 009 both match their references)

| build | change | dataset D total | TNBC wall | TNBC window loop | TNBC MS1_TRACE | TNBC INFER | TNBC peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|
| 0 | (before) | 312.1 s | 33:35 | 1150.0 s (31.2x) | 252.5 s (4.8x) | 438.8 s | 181,777 MB |
| 1 | bounded scoring pool, findPartner SoA | 266.2 s | 18:03 | 496.5 s (66.7x) | 247.3 s (4.8x) | 177.2 s | 181,902 MB |
| 2 | per-band seed segments | -- | -- | 450.4 s (78.1x) | 248.5 s (4.8x) | 154.0 s | -- |
| **3** | parallel band distribution | **237.5 s** | **15:01** | **467.5 s (69.5x)** | **106.2 s (11.2x)** | **155.7 s** | 187,928 MB |

**TNBC 33:35 -> 15:01 (-55%), dataset D 312.1 -> 237.5 s (-24%)**, spectrum list byte-identical on both
at every step. Peak memory is +3.4% (the chunked band sub-maps); `perf:malloc_trim` and the record
work still hold the -34% from the memory pass.

What is left, on TNBC: the window loop (467.5 s, 52%) is still one window's serial chain and that
chain is now 86% **trace**, capped at 12 bands; the greedy claim loop in precursor inference
(138.3 s of the phase's 155.7) is sequential by construction; the mzML write is 53.0 s at one
thread. Whole-run occupancy is 43x of 100 threads.

### `perf:trace_bands` 48: FAILS the gate, and stays off (2026-09-07)

The speed case was strong -- window loop -30% on both samples, occupancy 70.8x -> 94.2x (dataset D) and
67.6x -> 89.9x (dataset A), at 7% LOWER CPU -- so it went through the full gate protocol (both engines,
entrapment, datasets D and A, arms back to back on one node from one binary with only the flag
changed).

| sample | Sage 12 -> 48 | MSFragger 12 -> 48 | Sage set overlap | MSFragger set overlap | entrapment FDR |
|---|---|---|---:|---:|---|
| dataset D | 12,594 -> 12,607 (+0.10%) | 11,640 -> 11,643 (+0.03%) | 99.55% | 99.77% | 1.26% -> 1.29% |
| dataset A | 11,218 -> **11,027 (-1.70%)** | 10,689 -> **10,928 (+2.24%)** | 95.46% | 89.74% | 1.12% -> 0.99% |

**Rejected.** D is clean, but A both loses more than 1% in one engine and contradicts it in the
other -- the two engines move in opposite directions by 1.7% and 2.2%, on peptide sets that agree
only 95% and 90% between arms. That is the decision rule's "A contradicts D" case, and a pure
speed knob does not get to move 5-10% of the peptide set. `perf:trace_bands` stays at 12; the
numbers are recorded here so the trade is documented rather than rediscovered.

The speed it was buying has since been taken output-identically instead (below).

### Computed indices instead of binary searches (2026-09-07)

Two hot loops were binary-searching sorted arrays whose key is, or can be, an integer:

- `best()`, the innermost step of integer tracing, ran a `lower_bound` over its frame's ~25,000
  sorted flight-time bins on **every frame step of every trace extension**. The bin is an integer
  the instrument produced, so a per-frame direct-address row (~32 peaks per bucket, one global
  shift so boundaries align across frames) turns the search start into a subtract and a shift.
- `findPartner`, in the greedy claim loop, searched 22e6 doubles up to 50 times per seed. The key
  there is the top bits of the IEEE-754 pattern -- monotone in the value for positive doubles, so
  it is one load and one shift, and because the exponent leads, the bucket width is RELATIVE,
  which is the right shape for a ppm tolerance.

Both are exact: the bucket start is at or below the answer and the next bucket start is at or
above it, so the standard-library search runs on a sub-range that provably contains the same
position. Measured, byte-identical on both files:

| build | dataset D total | dataset D window loop | TNBC wall | TNBC window loop | TNBC INFER (greedy walk) |
|---|---:|---:|---:|---:|---:|
| before | 237.5 s | 143.9 s | 15:01 | 467.5 s (69.5x) | 155.7 s (138.3) |
| + m/z key | 234.2 s | -- | 14:32 | 459.8 s (70.4x) | 130.0 s (112.8) |
| **+ TOF index** | **205.2 s** | **116.7 s** (CPU -9%) | **14:11** | **412.9 s (74.3x)** | 148.7 s (128.8) |

`PRECURSOR_INFER` replicates to about +-15% on a shared node between runs of the same binary (130.0 vs
148.7 s for the same code), so read that column as a level, not a delta.

### Where the day ended

| | before | after | change |
|---|---:|---:|---|
| dataset D total | 312.1 s | 205.2 s | **-34%** |
| TNBC 009 wall | 33:38 / 26:45 (see below) | 14:11 | **-58% / -47%** |
| TNBC window loop | 1150.0 s (31.2x) | 412.9 s (74.3x) | -64% |
| TNBC MS1_TRACE | 252.5 s (4.8x) | 114.8 s (10.7x) | -55% |
| TNBC PRECURSOR_INFER | 438.8 s | 148.7 s | -66% |
| TNBC peak RSS | 181,777 MB | 190,849 MB | +5.0% |

**Read the TNBC delta with care.** Two runs of the SAME baseline binary with IDENTICAL digests
measured 33:38.00 and 26:44.58 -- 26% apart on node load alone, on a shared machine. The tables here
quote the slower one as the baseline, which flatters the result; against the faster one the 14:10.94
of the final build is
**-47%**. Both bracket the same conclusion, but no single-pair wall-clock claim on a shared node is worth
more than about +-15% unless the arms are interleaved.

Byte-identical spectrum list on both files at every step, 13/13 e2e throughout. The memory cost is
the chunked band sub-maps and the per-frame TOF index; it is the only regression and is reducible
by lowering the chunk count if RAM becomes binding again.

Still serial or thin, on TNBC: the greedy claim loop (~129 s of `PRECURSOR_INFER`, sequential by
construction -- an ordered-speculation scheme would make it exactly parallel but with unknown
replay cost), the mzML write (52 s, one thread), and the window loop's tail, which is one window's
prep -> trace -> split -> score chain with trace capped at 12 bands. Whole-run occupancy is 45x of
100 threads, up from 24x.

The **m/z key is a stand-in**: the natural key is the flight-time bin itself, which needs no
calibration at all, but `Trace::tof` is 0 on the OpenMS MS1 path (`toTrace` never sets it) and MS1
tracing still runs that path. Keying precursor inference on tof, and deferring the calibrated m/z
to the end, becomes available once MS1 goes through the integer detector.


### Ordered speculation on the greedy claim loop (2026-09-07)

Precursor inference was the last single-threaded phase. The loop is greedy -- each seed marks its
isotope partners used -- but it parallelises EXACTLY, because the only mutable state a seed reads is
the `used` bitmap and that bitmap is monotone. A batch is evaluated concurrently against a frozen
bitmap, recording the partner indices each seed's searches returned; the batch then commits in the
original order, re-evaluating only a seed whose recorded partners have since been claimed.

The dependency set is small and provable: `findPartner` reads `used[j]` only for candidates that
passed the RT/IM gate and returns the first unused one, so a call returning `j` depended on
`used[j]` being false, while a call returning -1 saw every candidate either fail the gate or already
be used -- both stable under a monotone bitmap.

| | greedy walk | phase | parallelism | re-evaluated |
|---|---:|---:|---:|---:|
| dataset D, before | 11.9 s | 15.7 s | 1.0x | -- |
| dataset D, after | **1.4 s** | **4.9 s** | **23.2x** | 2.5% |
| TNBC, before | 128.8 s | 148.7 s | 1.0x | -- |
| TNBC, after | **7.8 s** | **26.9 s** | **23.4x** | **0.7%** |

Byte-identical on dataset D AND on dataset A (the dataset A gate arm reproduced digest `61d57ca2...` exactly against
the earlier build), and on TNBC 009.

### TNBC 009 after the speculation change

| phase | wall | % | CPU | par |
|---|---:|---:|---:|---:|
| LOAD(stream) | 92.9 s | 13.1% | 5355.5 | 57.7x |
| MS1_TRACE | 117.9 s | 16.6% | 1182.4 | 10.0x |
| PRECURSOR_INFER | 26.9 s | 3.8% | 629.7 | 23.4x |
| WINDOW_LOOP | 380.7 s | 53.5% | 30372.0 | 79.8x |
| WRITE(mzML) | 52.4 s | 7.4% | 52.2 | 1.0x |
| **total** | **11:59.75** | | | peak RSS 188,593 MB |

The window loop is still one window's serial chain, and the `[tsub]` timers settle what that chain
is made of: for window 825-851, trace 304.6 s = slab 5.8 + **detect 286.8** + merge 5.3 + sort 10.8.
The band taskloop is 94% of it, so the serial merge and the record sort are NOT worth attacking --
band count is the only lever there.


### `perf:trace_bands` is settled: 24 fails dataset A the same way 48 did (2026-09-07)

Re-gated at the intermediate value on the current binary, arms back to back on one node:

| sample | bands | Sage | MSFragger | Sage overlap | MSFragger overlap | entrapment FDR |
|---|---|---|---|---:|---:|---|
| dataset A | 12 -> 24 | 11,218 -> **11,037 (-1.61%)** | 10,689 -> **10,929 (+2.25%)** | 95.56% | 89.66% | 1.12% -> 1.14% |
| dataset A | 12 -> 48 | 11,218 -> 11,027 (-1.70%) | 10,689 -> 10,928 (+2.24%) | 95.46% | 89.74% | 1.12% -> 0.99% |

The churn is the SAME SIZE at 24 and at 48, so it is not a dose effect of having more bands -- it is
a one-off change in which band owns a trace at a boundary. Any band count other than the one the
baseline was measured at moves ~4-10% of the peptide set and the two engines disagree about the
sign. **`perf:trace_bands` stays at 12**, and with it the window loop stays at its current floor for
output-identical work: the `[tsub]` timers show 94% of the critical window's trace stage is the band
taskloop itself, with the serial merge (5.3 s) and record sort (10.8 s) too small to matter.

That sensitivity is itself the finding: band boundaries change trace ownership because each band
carries its own `visited` set and a trace grown by one band can consume peaks in another band's core.
Making the partition genuinely band-independent would remove the churn AND turn band count into a
free speed knob.


### `perf:ms1_trace_bands`: the old sweep was measuring the serial distribute (2026-09-07)

The 12/24/48 sweep that concluded "more MS1 bands is worse" (75.0/93.7/147.9 s) was taken BEFORE the
band distribution was parallelised. That step was the thing growing with band count. Re-measured on
the current build, the conclusion inverts completely:

| `ms1_trace_bands` | MS1_TRACE | parallelism | detect | old sweep |
|---:|---:|---:|---:|---:|
| 12 (default) | 41.7 s | 9.8x | 31.4 s | 75.0 s |
| 24 | 27.1 s | 18.4x | 16.8 s | 93.7 s |
| 48 | **20.4 s** | **32.4x** | 8.4 s | 147.9 s |

And unlike `perf:trace_bands`, it is nearly output-neutral, because the MS1 path really does use the
halo + core-ownership scheme its comment claims:

| sample | Sage 12 -> 48 | MSFragger 12 -> 48 | Sage overlap | MSFragger overlap | entrapment |
|---|---|---|---:|---:|---|
| dataset D | 12,594 -> 12,593 (0.9999) | 11,640 -> 11,642 (1.0002) | **99.44%** | **99.37%** | 1.26% -> 1.24% |

Compare the same table for `trace_bands` on dataset D at 24: 99.87% / 99.71%, and on dataset A at 24: 95.56% /
89.66%. **Dataset A is intrinsically far more sensitive to MS2 band boundaries tha dataset D is** -- which is
why a knob that looks free on the primary dataset has to be gated on the second one.


### `ms1_trace_bands=48` adopted (2026-09-07)

Gate, arms back to back on one node from one binary with only the flag changed:

| sample | Sage 12 -> 48 | MSFragger 12 -> 48 | Sage overlap | MSFragger overlap | entrapment FDR |
|---|---|---|---:|---:|---|
| dataset D | 12,594 -> 12,593 (0.9999) | 11,640 -> 11,642 (1.0002) | 99.44% | 99.37% | 1.26% -> 1.24% |
| dataset A | 11,218 -> 11,211 (0.9994) | 10,689 -> **10,906 (1.0203)** | 99.24% | 89.33% | 1.12% -> **1.07%** |

Neutral-or-better on both engines and both datasets, with entrapment FDR improving on both. Adopted.

One honest caveat on the metric: MSFragger's accepted set on dataset A shows ~89-90% overlap between arms
in EVERY gate run today, including pairs where Sage overlaps 99.2-99.6% on the same two files. Its
peptide-level target-decoy walk is simply a much noisier statistic on that file than Sage's q-value
cut, so a 2% MSFragger move on dataset A should not be read as a 2% real gain.

#### dataset D with `ms1_trace_bands=48`

| phase | wall | % | par |
|---|---:|---:|---:|
| LOAD(stream) | 21.8 s | 10.7% | 74.3x |
| MS1_TRACE | **17.8 s** | 8.8% | **22.7x** |
| PRECURSOR_INFER | 5.6 s | 2.7% | 23.6x |
| WINDOW_LOOP | 129.6 s | 63.7% | 72.0x |
| WRITE(mzML) | 15.5 s | 7.6% | 1.0x |
| **total** | **3:27.48** | | |

Dataset D total 312.1 s before the day's work -> 192.5 s measured, **-38%**. The window loop is now 64% of the run
and is the floor: 94% of the critical window's trace stage is the band taskloop, and band count is
not adoptable.

### Trace-detector speculation: implemented, measured, removed

The technique that made precursor inference 23x is the right one in the wrong place. Precursor
inference runs ALONE on an otherwise idle machine; the band detector runs while 28 windows x 12
bands already saturate the pool, so speculative work there is pure waste -- a fixed width of 8 cost
**+31% window-loop CPU and made the wall 11% worse**. Gating the width on live occupancy removed the
penalty, but then it barely engages, because the window whose detect is the critical path is by
definition the last one: it only sees idle threads in its final ~20%.

Interleaved A/B, two rounds, one node, digests identical throughout:

| | round 1 | round 2 | mean |
|---|---:|---:|---:|
| control window loop | 133.0 s | 121.7 s | 127.4 s |
| window loop with adaptive speculation | 121.1 s | 121.1 s | 121.1 s |
| control critical-window detect | 89.7 s | 84.0 s | |
| speculation critical-window detect | 89.5 s | 94.1 s | |

The control's own spread (133.0 vs 121.7) is larger than the difference, and the stage it targets did
not move. Removed: ~200 lines of speculation machinery in the hottest loop is not worth an
unmeasurable win. The design is kept on record, because it becomes worthwhile the moment the loop
stops being tail-bound.


### New TNBC 009 reference digest (2026-09-07)

`ms1_trace_bands=48` changes the spectrum list, so the previous TNBC reference is retired.
**TNBC 009 reference digest = `da1b86fb5c38409f...`**. Measured on a 224-core node:

| phase | wall | % | CPU | par | (previous build, same file) |
|---|---:|---:|---:|---:|---:|
| LOAD(stream) | 94.5 s | 13.7% | 5688.6 | 60.2x | 92.9 s |
| MS1_TRACE | **59.3 s** | 8.6% | 1596.8 | **26.9x** | 117.9 s (10.0x) |
| PRECURSOR_INFER | 24.8 s | 3.6% | 769.9 | 31.0x | 26.9 s |
| WINDOW_LOOP | 412.0 s | 59.8% | 32968.6 | 80.0x | 380.7 s |
| WRITE(mzML) | 52.4 s | 7.6% | 52.2 | 1.0x | 52.4 s |
| **total** | **11:35.90** | | | | 11:59.75 |

peak RSS 190,852 MB. An independent run of the same build on another node gave **11:24.97 and 179.1 GB**
-- the two nodes differ by ~6% on RSS and ~2% on wall, which is the usual spread.

Where the time now is, on a 2-hour acquisition: the window loop (60%), and after it the mzML write
(7.6%, one thread) and the streaming load (13.7%, 60x). The window loop is tail-bound and its tail is
94% band taskloop, which band count cannot unlock without moving the peptide set.


### `SPEXTRACTOR_PICK_BATCH`: fast, and rejected (2026-09-07)

LOAD is 13.7% of a 2-hour run at ~60x of 100 threads: a read (serial) -> pick (parallel) pipeline
over batches of `kBatch` frames, where a batch not much larger than the thread count caps occupancy
and the reader between batches is a bubble. The knob is already an env var, so it was free to sweep
(dataset D):

| `PICK_BATCH` | LOAD | par | total | peak RSS | digest |
|---:|---:|---:|---:|---:|---|
| 256 (default) | 27.2 s / 19.8 s | 54.4x / 73.6x | 185.4 / 178.3 s | 55.5 GB | `9e5e5514` (both runs) |
| 1024 | 26.2 s | 56.4x | 183.8 s | 55.3 GB | `561446fb` |
| 4096 | **16.6 s** | **82.9x** | **165.3 s** | **51.6 GB** | `a3418b5e` |

Deterministic at a fixed batch size (256 reproduced its digest exactly), but the digest CHANGES with
batch size: the picker carries per-thread state, so which frames a thread sees changes its output --
the "~2% of intensities at the 1e-5..1e-4 level" the load comment already records for serial vs
parallel picking.

**Correction (2026-09-09):** the cause was not the picker -- the loader dropped its last
`n mod kBatch` MS2 spectra (two whole windows at 4096); see "tiling plan step 0", 0b below.

**Gated and REJECTED.** That perturbation is not benign:

| | Sage | MSFragger | Sage overlap | MSFragger overlap |
|---|---|---|---:|---:|
| dataset D, 256 -> 4096 | 12,593 -> **11,682 (-7.2%)** | 11,642 -> **10,930 (-6.1%)** | 91.7% | 81.3% |

A 7% peptide loss for a load-phase speedup is not a trade worth making, and in the gate's own arms
the wall was not even better (3:12.17 vs 3:13.32). `kBatch` stays at 256. The interesting part is
what it reveals: **the peak picker is materially sensitive to how frames are grouped across threads**,
which is a correctness smell worth its own investigation.


## Head to head with the reference implementation, measured tonight (2026-09-08)

Both tools on TNBC 009, on one 224-core node, within the same hour, both given 100 threads. The
recorded reference timing for this file was from a previous campaign, and comparing across time is
the error that produced one withdrawn claim -- so it was re-measured.

| | rep 1 | rep 2 | mean |
|---|---|---|---|
| reference implementation | 21:12.24 / 357.7 GB | 20:41.46 / 323.4 GB | 20:56.9 / 340.6 GB |
| **DIAspeXtract (min_charge 2)** | 11:24.97 / 179.1 GB | 10:36.92 / 178.5 GB | **11:00.9 / 178.8 GB** |

**1.90x faster at 0.53x the peak memory.** For context this file was 1.46x SLOWER than the reference
before this week's work (31:26.98 against 21:28.00).

One caveat on the memory ratio: the reference is a JVM run with `-Xmx400G`, so its peak RSS partly
reflects what it was allowed rather than what it needs -- the same flags gave 268.2 GB in the earlier
campaign. DIAspeXtract's 179 GB is a genuine requirement. The honest statement is "roughly half, and
not sensitive to a heap flag".


### `charge:min_charge` 1 + the ion-mobility charge veto -- dataset D (2026-09-08)

See docs/REFERENCE-COMPARISON-2026-09-08.md for why. Arms: the previous default (min_charge 2)
against min_charge 1 with the veto on, same file, same engines.

| | previous default | min_charge 1 + veto | |
|---|---:|---:|---|
| spectra | 656,257 | 927,712 | +41% |
| wall / peak RSS | 3:27.48 / 57.6 GB | **3:25.95 / 57.9 GB** | flat |
| veto | -- | of 407,478 z=1 calls: 107,040 -> 2+, 93,008 -> 3+, 207,430 kept | |
| Sage peptides (engine's global cut) | 12,593 | 12,536 (z=1: 202) | -0.45% |
| MSFragger peptides (global) | 11,642 | 11,581 (z=1: 163) | -0.52% |
| Sage, **per-charge** walk | 12,493 | **12,597** (z=1: 46) | **+0.8%** |
| MSFragger, **per-charge** walk | 12,360 | **12,700** (z=1: 28) | **+2.8%** |
| entrapment FDR, whole run | 1.24% | 1.38% | |
| entrapment by stratum (min_charge 1 + veto) | | z=1 **9.2%** (n=144, CI 3.9-16), z=2 1.05%, z=3 1.64% | |

Dataset D is the case the old default was tuned on, and it is where z=1 is rarest: 202 peptides against
5,000+ on a 2-hour TNBC file. Under the engines' own pooled cut the arm is -0.5% on both engines --
inside the gate's 1% rule -- and the loss is the pooled-threshold eviction documented in
docs/REFERENCE-COMPARISON-2026-09-08.md:
under a per-charge walk the same files give +0.8% and +2.8%. The z=1 stratum's 9% is real but small
(144 peptides) and is what per-charge control is for.

The veto's **3+ arm is the weak spot**: 93,008 re-calls to 3+ moved the z=3 count from 3,888 to
3,935 (global) and 3,817 to 3,832 (per-charge) -- essentially nothing -- so on this file they are
emission without identifications. The fitted 2+ and 3+ mobility lines are only ~1 sigma apart at low
m/z (at m/z 500: 0.867 vs 0.822 1/K0 against a 3+ sigma of 0.051), so a multiply charged z=1 in that
overlap is assigned 2 or 3 by whichever residual is smaller. The principled tie-breaker is the
intermediate-isotope test (a co-eluting trace at +ISO/2 says 2+, at +ISO/3 says 3+), which the
`SPEXTRACTOR_Z1_DIAG` machinery already computes; not yet wired in.


### New TNBC 009 reference digest (2026-09-08)

`charge:min_charge` 1 + the mobility veto change the spectrum list, so the previous digest is retired.
**TNBC 009 reference digest = `88f89abad8eca9ca...`**. Extraction facts:

| | previous (min_charge 2) | new (min_charge 1 + veto) |
|---|---:|---:|
| spectra | 2,343,822 | **4,035,234** (+72%) |
| veto | -- | 2,548,391 z=1 calls: 287,900 -> 2+, 178,226 -> 3+, 2,082,265 kept as 1+ |
| peak RSS | 190.9 GB | **190.3 GB** (flat: emission is not what bounds memory) |
| window loop | 412.0 s (80.0x) | 511.3 s (**90.6x**) |
| wall | 11:35.90 | 13:43.10 on another node, **contended** by a concurrent MSFragger -- not a head-to-head number; a clean rerun follows |

The 72% more spectra cost no memory and ~+11% wall when uncontended (the earlier z>=1 arm ran
in ~13:20 against 11:59 on the same node); the window loop's occupancy rose to 90.6x because
more precursors means more scoring work to fill the tail with.


### The veto's 3+ arm dropped -- dataset D (2026-09-08)

Same arms as the table above; the shipped default differs from it only in dropping z=1 calls that sit
on the 3+ mobility band instead of re-labelling them 3+.

| | previous default | veto, 3+ re-called | **shipped default** |
|---|---:|---:|---:|
| spectra | 656,257 | 927,712 | **854,817** |
| wall / peak RSS | 3:27.48 / 57.6 GB | 3:25.95 / 57.9 GB | **3:12.61 / 57.3 GB** |
| window loop (run log) | 129.6 s (72.0x) | -- | 115.3 s (78.5x) |
| Sage peptides, engine's pooled cut | 12,593 | 12,536 (-0.45%) | **12,609 (+0.13%)** |
| MSFragger peptides, pooled cut | 11,642 | 11,581 (-0.52%) | 11,564 (-0.67%) |
| Sage, per-charge walk | 12,493 | 12,597 (+0.8%) | **12,639 (+1.2%)** |
| MSFragger, per-charge walk | 12,360 | 12,700 (+2.8%) | **12,666 (+2.5%)** |
| entrapment, whole run | 1.24% | 1.38% | 1.33% |
| entrapment z=1 stratum | -- | 9.2% (n=144) | 8.8% (n=151) |

The 93,008 dropped calls were pure emission: peptide counts are unchanged or better on both walks,
73k fewer spectra than the 3+-re-calling arm, and the fastest dataset D wall measured. On dataset D the new default is
inside the gate's 1% rule under the engines' own pooled cut and positive on both engines under
per-charge control; the z=1 stratum here is tiny (~150 peptides) and dirty, which is what per-charge
FDR is for.


### The shipped default on TNBC 001 (2026-09-08)

Previous-default arm = min_charge 2; shipped default = min_charge 1 + mobility veto with the 3+ arm
dropped. Extraction on a shared, loaded node (not a timing number): 16:02, 170.1 GB, 3,409,189
spectra (previous default: 2,135,267). Veto: of 2,137,636 z=1 calls, 241,062 (11.2%) re-called 2+, 157,532 dropped, 1,739,042 kept.

| | previous default | **shipped default** | vs reference |
|---|---:|---:|---|
| Sage peptides, engine's pooled cut | 15,694 | **19,035 (+21.3%)** | 14,764 -> **1.289** |
| Sage proteins | 4,074 | 4,142 | 3,802 -> 1.089 |
| MSFragger peptides, pooled cut | 15,976 | **18,212 (+14.0%)** | 18,745 -> 0.972 |
| Sage, per-charge walk | 16,029 | **19,262 (+20.2%)** | |
| MSFragger, per-charge walk | 16,727 | **19,230 (+15.0%)** | |
| entrapment, whole run | 0.98% | 0.93% | |
| entrapment by stratum (shipped) | | z=1 2.91% (n=3,481), z=2 0.41%, z=3 0.57% | |
| per-charge walk entrapment (shipped) | | z=1 **0.87%**, z=2 1.14%, z=3 0.68%, whole 0.99% | |

Second file, same shape as 009: the z=1 stratum sits at ~2.9% under the pooled cut and at 0.87%
under per-charge control, with the gain intact on both engines under both walks.


### The shipped default on TNBC 009, and the new reference digest (2026-09-08)

**TNBC 009 reference digest = `057aad4beec52483...`** (the shipped v2.0.0 defaults: min_charge 1,
mobility veto, 3+-band calls dropped); it supersedes the previous default's digest. Extraction alone
on the node:

| | previous default | veto, 3+ re-called | **shipped default** |
|---|---:|---:|---:|
| spectra | 2,343,822 | 4,035,234 | **3,918,321** |
| wall, uncontended, same node | 11:24.97 / 10:36.92 | 13:07.40 / 12:36.77 | **12:42.21** |
| peak RSS | 179.1 / 178.5 GB | 177.4 / 180.0 GB | 187.9 GB |
| window loop | 412.0 s (80.0x) | 442.4 s (91.1x) | 443.4 s (89.2x) |
| Sage peptides, pooled | 24,213 | 27,808 | **27,855** |
| MSFragger peptides, pooled | 24,199 | 26,588 | **26,559** |
| Sage, per-charge walk | 24,484 | 27,908 | **27,831 (+13.7%)** |
| MSFragger, per-charge walk | 25,029 | 28,185 | **28,181 (+12.6%)** |
| entrapment, whole run | 1.01% | 1.00% | 1.03% |
| entrapment z=1 stratum (pooled / per-charge) | -- | 2.91% / 0.91% | 2.94% / 0.91% |

Emitting the singly charged population costs ~+15% wall on a 2-hour file (12:42 against 11:01 mean
for the previous default) and no memory beyond the usual +-10 GB run-to-run spread; the window loop's occupancy rose
to ~90x because there is more scoring work to fill its tail with.

### Head to head at the new default (2026-09-08, one node, both tools 100 threads)

| | reference implementation | DIAspeXtract | ratio |
|---|---|---|---|
| wall | 21:12.24 / 20:41.46 (mean 20:56.9) | **12:42.21** | **0.61x** |
| peak RSS | 357.7 / 323.4 GB (mean 340.6) | **187.9 GB** | **0.55x** |
| Sage peptides, per-charge 1% -- 009 | 21,743 | **27,831** | **1.28** |
| MSFragger peptides, per-charge 1% -- 009 | 28,294 | **28,181** | **1.00** |
| Sage, per-charge -- 001 | 15,418 | **19,262** | **1.25** |
| MSFragger, per-charge -- 001 | 19,652 | **19,230** | **0.98** |
| entrapment FDR, whole run (SpeX) | -- | 1.03% / 0.93% | at nominal 1% |


### v2.0.0 published (2026-09-08)

Published to the shared install, built from 8d79bbc (`src` sha256 `2a45e5e9...`, binary
`0da8cfd7...`); the release tag sits on the doc-only follow-up commit, whose source is byte-identical.
e2e 13/13. Dataset D on a lightly loaded node: **854,817 spectra, 3:17.41, 57.8 GB, window loop
120.6 s at 78.1x**, digest **`b499027f...`** -- identical to the release candidate built from the
same source with the old help strings, and the same count as the pre-release build. That digest is
the dataset D reference for v2.0.0. TNBC 009 with the published binary on a loaded node:
**3,918,321 spectra, 14:36.99, 189.9 GB, window loop 552.6 s at 91.9x**, digest **`057aad4b...`** --
the release reproduces the TNBC 009 reference digest exactly (the release candidate built from the
same source did too: 14:19.62 / 189.4 GB / identical). The uncontended wall for this file remains
12:42.21.

### Savitzky-Golay coefficient caching (2026-09-08, evening)

OpenMS rebuilt the smoothing filter for every mass trace. `ElutionPeakDetection::smoothData()`
constructs a `SavitzkyGolayFilter` and calls `setParameters()` per call, and `updateMembers_()`
recomputes the coefficients with one Eigen SVD per half-window each time. The coefficients depend
only on (frame_length, polynomial_order), which take a handful of values in a run, while
smoothData() is called once per trace -- hundreds of millions of times. A thread_local cache of the
configured filter leaves the coefficients, the filtering loop and every float conversion untouched.
Shipped in `patches/openms-epd-lockfree.patch`, which now carries both EPD changes.

Two smaller changes landed in the same build: the MS1 band distribution called `isSorted()` -- a
full traversal -- once per band per spectrum when the spectrum cannot change between bands, and the
PeakMap assembly moved spectra one at a time where the vector can be transferred whole.

**Interleaved A/B, TNBC 009, `data` alone at load 4, 100 threads, arms alternating:**

| | rep 1 | rep 2 | mean |
|---|---|---|---|
| before | 13:44.51 / 190.6 GB | 13:36.80 / 188.3 GB | 13:40.7 / 189.5 GB |
| **after** | 12:07.35 / 185.2 GB | 12:26.96 / 186.6 GB | **12:17.2 / 185.9 GB** |
| delta | | | **-83.5 s (-10.2%), -3.6 GB** |

All four runs reproduce the TNBC reference digest. On dataset D the phase CPU shows where it comes
from (arms at different node loads, so CPU seconds rather than wall):

| phase | before | after | |
|---|---:|---:|---|
| LOAD | 1769 | 1692 | -4.3% |
| MS1_TRACE | 411 | 228 | **-44.4%** |
| PRECURSOR_INFER | 179 | 114 | -36.3% |
| WINDOW_LOOP | 10132 | 9105 | **-10.1%** |
| ASSEMBLE | 0.8 | 0.0 | phase gone |
| total | 12513 | 11141 | **-11.0%** |

### Two detector hypotheses priced and dropped (2026-09-08, evening)

`perf_event_paranoid` is 4 on these nodes, so the detector was priced by ablation instead. The
window loop's trace stage is 93% `detectTracesInteger_` (dataset D: detect 1330 s of 1424 s summed
over windows), so both arms targeted its inner loop. Neither is worth taking:

- **`sdRobust`'s log-space form.** It computes `sqrt(wsum*sd^2 + w*(mz-mean)^2)` through three
  `log`s and two `exp`s (OpenMS's overflow-safe formulation). Replacing it with the direct algebraic
  form -- which would be output-changing, hence an ablation rather than a candidate -- measured
  **+0.9% window-loop CPU**, i.e. nothing. The transcendentals are not on the critical path.
- **The calibration call in `best()`.** Substituting a flight-time-bin distance for
  `mzOf()` measured **-1.8%**. That is the whole prize for a change the record already shows costs
  1.4% of peptides.

Together these say the detector is bound by memory traffic over its three parallel arrays and the
visited bitmap, not by arithmetic -- which is where any further work on it should aim.

### The detector's candidate scan: reordering its predicates LOSES (2026-09-08, night)

`best()` scans a frame's candidates in flight-time order and rejects on two predicates: intensity
above noise (`inten`, 4 B) and mobility inside the gate (`imq`, 2 B). Testing mobility first, and
replacing the per-step `tof[k] <= hi_t` condition with one `upper_bound` computed before the loop,
looks free -- both predicates are pure functions of k, so the candidate set, the visiting order and
the first-of-equals rule are all unchanged, and the spectrum list is byte-identical (verified).

It is a **10% regression**. Interleaved on an otherwise quiet node at matched load, dataset D:

| | rep 1 | rep 2 | mean wall | mean window-loop CPU |
|---|---|---|---|---|
| unchanged | 3:16.82 | 3:09.32 | **3:13.1** | 7,912 |
| reordered | 3:35.31 | 3:30.22 | 3:32.8 | 8,719 (+10.2%) |

The reason is the hardware prefetcher. The original loop reads `tof[k]` (the loop condition),
`inten[k]` and `imq[k]` in ascending k, so three sequential streams run and the L2 streamer keeps
all of them ahead of the loop. Testing the narrow array first and dropping the `tof` read breaks two
of those streams into dependent, scattered accesses: the ~5% of candidates that pass mobility then
fetch `inten` and `tof` at unpredictable offsets. Touching *fewer bytes* is not the same as touching
*fewer cache lines*, and here the sequential access pattern was worth more than the bytes saved.

Taken with the two ablations above, the detector's inner loop is bound by memory traffic that is
already being prefetched efficiently. Any further work there has to reduce the number of candidates
visited -- not the cost of visiting one.

### Parallel mzML serialisation, and the key-sort (2026-09-08, night)

Two changes to the serial tail, both byte-identical (dataset D reproduces the reference digest).

**The write.** It was the largest serial block left: 92.9 s of an 822 s run on the 2-hour file
(10.7%) at 1.0x. It is CPU-bound rather than I/O-bound, which is what made it worth attacking --
measured, the phase sustains ~400 MB/s while `dd` on the same filesystem does **4.0 GB/s**, a 10x
headroom. Serialising one spectrum is independent of every other; only the concatenation is ordered.
Spectra are now encoded in parallel into per-chunk buffers and appended in index order, in waves so
the buffered text stays bounded.

The obstacle, and why the first three attempts aborted: OpenMS's writer is not reentrant on one
handler. `writeSpectrum_`'s helpers are non-const, `MzMLHandler::cached_terms_` is a mutable memo of
CV-term validation written during serialisation, and the validator keeps its own state. A shared
handler corrupts the heap within seconds -- confirmed under gdb, aborting inside `writePrecursor_`'s
allocation, which is where the damage is *detected* rather than caused. Freezing `cached_terms_`
alone was not enough. Giving each thread its own handler and validator -- same members, same
`PeakFileOptions`, sharing only the read-only experiment -- is what works.

| dataset D, 100 threads | before | after |
|---|---:|---:|
| WRITE(mzML) wall | 20.1 s | **5.9 s** |
| WRITE parallelism | 1.0x | **51.4x** |

**The sort.** `std::sort` moves whatever it is handed, and an `MSSpectrum` object is ~750 bytes
against 32 for a key. Sorting a key array and permuting once gives the same permutation -- the
comparator returns exactly what the old one returned for every pair, the rare full tie still defers
to the peaks, and the keys start in the spectra's order, so introsort makes the same decisions.
SORT(canonical) 1.9 s -> 0.8 s on dataset D, and the key build parallelises where sorting objects
could not.

### The night's work, measured end to end (2026-09-08 -> 09, TNBC 009)

Interleaved A/B on `data`, arms alternating, 100 threads, against the released v2.0.0. **All four
runs reproduce the TNBC reference digest**, so the peptide sets are unchanged by construction and no
search was run.

| | rep 1 | rep 2 | mean |
|---|---|---|---|
| v2.0.0 | 13:43.65 / 189.2 GB | 13:42.13 / 188.6 GB | 822.9 s / 188.9 GB |
| **tonight** | 11:09.22 / 186.7 GB | 11:26.01 / 188.0 GB | **677.6 s / 187.4 GB** |
| delta | | | **-145.3 s, -17.7%** |

The improved arm ran at node load 59-65 against 17-19 for the baseline, so the figure is if anything
understated. Peak memory is unchanged within run-to-run spread.

Phase means, wall and CPU:

| phase | v2.0.0 wall | tonight wall | v2.0.0 CPU | tonight CPU |
|---|---:|---:|---:|---:|
| LOAD(stream) | 91.0 | 93.5 | 5309 | 5266 |
| MS1_TRACE | 58.6 | 52.9 | 1489 | **888 (-40%)** |
| PRECURSOR_INFER | 24.5 | 23.4 | 607 | 602 |
| WINDOW_LOOP | 491.5 | **432.0** | 44523 | **34491 (-22.5%)** |
| SORT(canonical) | 9.5 | **3.7** | 9.5 | 5.3 |
| ASSEMBLE(PeakMap) | 3.8 | **0.0** | 3.8 | 0.0 |
| WRITE(mzML) | 90.0 | **21.1** | 90 | **1037 at 49x** |
| total measured | 768.7 | **626.5** | | |

Two things to read from it. The window loop's CPU fell 22.5% -- more than the 10% measured on the
30-minute file, because longer gradients have longer traces and the coefficient cache pays per
trace. And its parallelism fell 90.0x -> 81.6x: with less work per window the tail is a larger share,
so scheduling now has room it did not have before. The write went from 1.0x to ~49x and is no longer
the serial floor; LOAD, at 93.5 s and 14% of the run, is now the second-largest phase.

### Parallelising the MS1 trace conversion: 17% slower (2026-09-09)

Converting 22 million OpenMS `MassTrace` objects into the tool's own records is a serial loop, and
per point it does a binary search over the frame table -- an obvious parallel target, and the
conversion is most of the gap between the MS1 phase's stage timers (33 s) and its wall (53 s).
Splitting it into a parallel half (coordinates and the (frame, intensity) point list, which only
READ the store) and a serial half (the arena append, which must stay ordered because the offsets are
the trace's identity) is exactly output-identical -- both files reproduce their digests.

It is slower. Interleaved on the 2-hour file, arms alternating:

| | rep 1 | rep 2 | mean wall | MS1_TRACE | parallelism |
|---|---|---|---|---|---|
| serial (shipped) | 11:02.29 | 10:47.44 | **654.9 s** | **50.5 s** | 18.8x |
| parallel halves | 11:14.13 | 12:00.08 | 697.1 s | 59.3 s | 15.5x |

The tell is that its CPU seconds *fell* (899 -> 879) while wall rose: the work really is cheaper, but
the loop cannot pay for its own bookkeeping. Per trace it is a handful of binary searches, and
against that the point-list buffers cost more than they save -- 22 million small vector allocations
move from one thread into all of them and contend in the allocator, and the parallel half writes
adjacent doubles in three coordinate arrays. A version that pre-sized one flat buffer and avoided
per-trace allocation was then written and measured too (2026-09-09): one chunked flat buffer, no
per-trace allocation at all, coordinates and slice-sort in the parallel half, arena append in the
serial half, again output-identical. It is **also no better** -- interleaved on dataset D, three
reps each: 168.2 s serial against 171.6 s, and the MS1 phase 14.8 s against 15.3 s.

So two independent parallelisations both lose to the serial loop, which is the answer: the per-trace
work is a few binary searches over a 4,797-entry table that stays in cache, and against that the
coordination costs more than the work. The serial loop also walks `mts` in order, which prefetches;
both parallel versions scatter that walk across threads. This is the same lesson the detector's
candidate scan taught tonight -- a sequential access pattern is worth more than the arithmetic it
appears to waste. The loop stays serial, with a comment recording both attempts.

### Second file: the same work on the 30-minute acquisition (2026-09-09)

Interleaved, three reps each, arms alternating; every run reproduces the dataset D reference digest.

| | rep 1 | rep 2 | rep 3 | mean | peak RSS |
|---|---|---|---|---|---|
| v2.0.0 | 3:17.31 | 3:16.67 | 3:15.34 | 196.4 s | 57.5 GB |
| **tonight** | 2:53.79 | 2:56.19 | 2:51.12 | **173.7 s** | 57.2 GB |
| delta | | | | **-22.7 s, -11.6%** | unchanged |

Less than the -17.7% on the 2-hour file, and the difference is explained rather than noise: the
Savitzky-Golay cache pays once per mass trace, and a 2-hour gradient carries 22.2 M MS1 traces
against 4.9 M here, with longer spans. The write parallelisation scales the same way -- 20.1 s of a
196 s run against 90 s of an 822 s one. Both files gain, the longer one gains more, and neither
changes a byte of output.

### Longest-first band dispatch: no effect (2026-09-09)

The window's 12 band tasks are created in band order, and with 24-28 windows in flight that is
288-336 tasks queueing on 100 threads -- so a heavy task created last should land at the end and
extend the makespan. Dispatching heaviest-first is free of output risk (each band writes a private
store and the merge absorbs them in band order regardless of completion), and the band's seed count
is known before any of them runs, so the work proxy costs nothing.

It changes nothing. Interleaved on dataset D, three reps each: 176.8 s against 176.4 s of wall --
inside the noise -- and the window loop's CPU *rose* 7,891 -> 8,239 for the sort. The pool is simply
deep enough: with ~300 tasks over 100 threads there is always work to steal, so creation order does
not determine the tail. Reverted.

**Three rejections in one night, and they agree.** Reordering the detector's candidate scan (-10%),
parallelising the MS1 conversion twice (no gain, once 17% worse), and now reordering band dispatch:
the window loop's remaining slack is not in task ordering, not in per-candidate arithmetic, and not
in the conversion. It is in the memory traffic of the detector's inner loop, which is already
prefetched about as well as a sequential scan can be. Further gains there need the algorithm to
visit fewer candidates -- which is output-changing and needs the full peptide gate.

## 2026-09-09: `Trace` 48 -> 40 bytes, the RT stored as a frame index (adopted, digest-identical)

`Trace::rt` was a cached double. Every producer reports an RT that is already one of the owning
store's frame times: the integer detector assigns `rtAxis()[sl.rt_index[fap]]` directly, and
ElutionPeakDetection -- which runs on BOTH MS1 and MS2 at the shipped defaults
(`trace:ms{1,2}_split_valleys = 7.0`) -- calls `updateSmoothedMaxRT()` on the unsplit trace
(ElutionPeakDetection.cpp:463) and on every split child (:541), which sets
`centroid_rt_ = trace_peaks_[max_idx].getRT()` (MassTrace.cpp:507), i.e. one of the points this tool
built from `frame_rt`. So the 8 bytes were a copy of `frame_rt[frame0 + k]`.

Store `k` instead, in the two bytes that were tail padding. `sizeof(Trace)` 48 -> 40, pinned by the
static_assert. `rtOf()` is one load from a ~38 KB per-window table; the 1.69e13-visit scoring scan
never reads it (it reads `FragRt::rt`, the SoA copy, now built from `rtOf()`).

**`rt_at` is SIGNED, and that is what makes it exact.** The first implementation clamped it into the
window `trimToSpan` kept, and the dataset D digest came back `69a6066d...` against the pinned `b499027f...`
with 855,493 spectra against 854,817. A split counter localised it in one run:

    [soa-rt] off-grid 23140 = not-a-frame 0 (of which strictly between two frames 0) + trimmed-out 23140

`not-a-frame` is ZERO -- every producer's RT really is a frame time -- and all 23,140 were my own
clamp. Making the offset signed leaves `frame0 + rt_at` (the absolute store frame, and therefore the
RT) invariant under trimming, because `frame0` gains `lo` and `rt_at` loses it. With that,
`[soa-rt] off-grid 0` and the digest is `b499027f...` exactly.

### Measured: interleaved A/B on `data`, arms alternating, released HEAD (`final1`) vs `soa1`

| | final1 | soa1 | |
|---|---|---|---|
| TNBC 009 (2 h), 2 reps -- peak RSS | 185.75 GB | **180.25 GB** | **-5.50 GB (-2.96%)** |
| TNBC 009 -- wall | 658.1 s | 675.5 s | +2.65% |
| dataset D (30 min), 3 reps -- peak RSS | 57.66 GB | **55.76 GB** | **-1.90 GB (-3.29%)** |
| dataset D -- wall | 171.4 s | 177.4 s | +3.49% |

All four 2-hour runs `SPECTRUM DATA IDENTICAL`; dataset D digest `b499027f...`; e2e 13/13.

**The wall difference is occupancy, not work.** CPU seconds on the 2-hour file, which are
load-robust where wall is not:

| phase | final1 CPU s | soa1 CPU s | |
|---|---|---|---|
| LOAD(stream) | 5,253.8 | 5,289.1 | +0.67% |
| MS1_TRACE | 957.8 | 925.1 | -3.42% |
| PRECURSOR_INFER | 596.5 | 606.6 | +1.70% |
| WINDOW_LOOP | 34,247.7 | 34,294.8 | **+0.14%** |
| WRITE(mzML) | 1,110.9 | 1,024.9 | -7.74% |
| **total** | **42,166.8** | **42,140.5** | **-0.06%** |

The window loop -- where all 2.14e9 records live -- costs the same CPU to within 0.14%, and total
CPU is flat at -0.06%. The wall gap is one outlier: `soa1_r2` ran the loop at a parallel factor of
77.2x against 81-84x for the other three. Note that run reported the LOWEST `uptime` 1-minute
average of the four (44.56), which is why that figure is not a usable covariate -- it is sampled at
the end of the run and reflects our own tail, not the competing load during it. Two reps cannot
settle wall; the CPU-second evidence says no work was added.

**Adopted on memory at CPU parity, not on wall.** -5.50 GB of a 185.75 GB peak for a
byte-identical output.

### What this does and does not establish

A fourth external review predicted -5.5 GiB by correcting the naive arithmetic: parents are freed
at the split (`vector<Trace>().swap(in)`), so the 2.14e9 records are NOT all simultaneously
resident; roughly 0.6 x children + 0.4 x parents ~ 49 GiB of records are at the peak, not 95.8. The
measured -5.50 GB confirms that correction and refutes the -16 GiB I first claimed. **Every further
per-record byte is therefore worth ~0.7 GB at the peak, not ~2 GB** -- use this coefficient, not the
record-count product, to price the remaining field removals (`im` -> side array, `bframe`, `mz`).

Three of four external reviewers were told this claim was load-bearing and asked
to check it, and all three repeated the author's original error that child RT is an off-grid
weighted mean. Only the fourth read `MassTrace.cpp`. The digest, not the review, was the arbiter -- and
it caught a second error (the clamp) that no reviewer saw either.

## 2026-09-09: d(peak RSS)/d(windows in flight) = 2.85 GB per window -- measured, no code change

The number every memory estimate in this project has been guessing at. `perf:max_concurrent_windows`
was already registered, resolved and enforced, and defaults to 0 = unlimited; the logs proved it had
never bound (`admission high-water 28 of 28`). Two runs of TNBC 009 on `data`, back to back on a
quiet node (load 6.3), build `soa1`, with a 250 ms `/proc/<pid>/status` VmRSS sidecar on each:

| | 28 windows | 11 windows | |
|---|---|---|---|
| peak RSS (`time -v`) | 182.20 GB | **133.76 GB** | **-48.44 GB (-26.6%)** |
| peak VmRSS (sidecar, 2566/3055 samples) | 181.88 GB | 133.80 GB | agrees to 0.2% |
| wall | 673.3 s | 799.1 s | +125.8 s (+18.7%) |
| WINDOW_LOOP wall | 431.3 s | 560.0 s | +128.7 s |
| WINDOW_LOOP CPU | 34,318.9 s | 33,762.2 s | **-1.6%** |
| WINDOW_LOOP parallel | 79.6x | 60.3x | |

Both `SPECTRUM DATA IDENTICAL`. Window-loop CPU is flat, so the entire wall cost is occupancy:
fewer windows in flight leave threads idle, they do not create work.

**The exchange rate, which is the usable result:**

    d(peak) / d(window in flight) = 2.85 GB per window
    d(wall) / d(window in flight) = 7.4 s per window
    => 0.385 GB of peak per second of wall, or 2.60 s of wall per GB saved

For comparison, `perf:malloc_trim` (on by default) trades ~30% of peak for ~8.5% of wall. The
concurrency knob is the same kind of trade at a similar rate, and it is a flag, not a change.

**What this re-prices.** Every structural saving should be read against 2.85 GB/window: the 40-byte
`Trace` (-5.50 GB) is worth about 1.9 windows of concurrency, i.e. ~14 s of wall it does NOT cost.
That is the argument for the structural work -- not that it is bigger than the knob, but that it is
free where the knob is not.

**What it does NOT establish.** The relation is measured only over 11-28 windows. Extrapolating the
2.85 GB/window slope to zero gives a 102 GB intercept against a measured ~46 GB loop-start floor, so
the fit is not linear down to one window: something in the loop is not per-window (the emitted
spectra accumulate regardless of concurrency, and allocator retention differs between the arms).
Do not use the slope outside the measured range, and do not read the intercept as a floor.


## 2026-09-09: Fix A -- return the drained compact store when the last window has materialised (adopted)

One `malloc_trim(0)` where the last window's `toSlab` returned (58c0981, arm `fixA` = `soa2` + this),
gated identical on D (`b499027f`, e2e 13/13) and measured interleaved on `data` against `soa2`,
TNBC 009, 100 threads, 250 ms VmRSS sidecar, `MemAvailable` 977-1079 GB at every start:

| | `soa2` r1 | `fixA` r1 | `soa2` r2 | `fixA` r2 | mean |
|---|---:|---:|---:|---:|---:|
| peak RSS (`time -v`) | 168.93 GB | 139.23 | 169.17 | 136.91 | **169.05 -> 138.07 GB (-31.0 GB, -18.3%)** |
| sidecar peak | 169.07 | 139.33 | 169.23 | 136.83 | agrees to 0.2% |
| wall | 10:54.7 | 11:15.8 | 11:39.9 | 10:49.4 | 677.3 -> 662.6 s (-2.2%, within the pair spread) |
| digest | identical | identical | identical | identical | |

The trim fires at the moment the last window's slab exists (t ~ 13 s into the loop) and returns
37.7 GB in 4.4 s (`101,532 -> 63,860 MB`), which is the arena corpse the lifetime audit priced at
25-43 GB. At 0.14 s per GB it is 20x cheaper than the concurrency knob's 2.60 s/GB. This is the
baseline every tiling number is read against (review correction C8): **138 GB, not 182 or 169.**
Superseded structurally by the direct slabs (2d-ii: no second representation to return), which
keep the `[trim]` site for the load's bookkeeping.

## 2026-09-09: tiling plan step 0 -- four zero-code measurements, one bug, one correction

The 2D-tiling plan starts with measurements that
need no architecture. All on dataset D unless stated; `soa2` = HEAD before Fix A.

### 0a. Allocator tunables on TNBC 009 (data, back to back, load 4.9, sidecar 250 ms)

| arm | wall | peak RSS (`time -v`) | sidecar peak | minor faults | LOAD | digest |
|---|---:|---:|---:|---:|---:|---|
| control (`soa2`) | 671.2 s | 169.70 GB | 169.89 GB | 251.6 M | 94.7 s | IDENTICAL |
| `MALLOC_MMAP_THRESHOLD_=65536 MALLOC_TRIM_THRESHOLD_=128M` | 846.3 s (+26%) | **154.15 GB (-15.5)** | 154.24 GB | 501.9 M | 261.7 s | IDENTICAL |

Production is glibc (`ldd`: nothing else linked, `LD_PRELOAD` empty). Serving every >= 64 kB block
from mmap does return memory (-15.5 GB, `mmap=12,513 MB` at the loop start against 1,869 MB), but
the compact store's ~126 kB chunks then each cost a page-fault storm: LOAD 94.7 -> 261.7 s and
2x the faults. So the corpse is real free-list, and the tunables are the wrong tool: Fix A (one
`malloc_trim(0)` when the last window has materialised, commit 58c0981) targets the same memory at
the moment it dies. Measured separately below.

### 0b. `SPEXTRACTOR_PICK_BATCH` 256 vs 4096 with `[det]` slab digests -- it was a dropped tail

The 2026-09-07 entry above ("fast, and rejected", -7.2% Sage peptides at 4096, "the picker carries
per-thread state") is **wrong about the cause**. Per-window peak counts at both batch sizes:

| | 256 | 4096 |
|---|---:|---:|
| windows in the compact store | 24 | **22** |
| compact-store peaks | 1,255,770,798 | 1,098,791,267 (-156,979,531) |
| windows with identical picked peak count | 20 of 22 present | |
| 684.29-708.32 / 1134.06-1400.62 | 1,237 frames each | **absent** |
| 662.33-685.29 / 1045.52-1135.06 | 1,342 frames each | **915** each |

The deficit closes exactly: 48,665,355 + 62,781,911 + 31,459,728 + 14,072,537 = 156,979,531. And the
frame arithmetic closes exactly too: 32,210 MS2 window spectra mod 256 = **210** = 2 x 105 (1,342 ->
1,237); mod 4096 = **3,538** = 2 x 1,342 + 2 x 427. `PickCompactConsumer::ensureMapsAreFilled_`
(the tail flush) is overridden but never invoked -- `FullSwathFileConsumer` reaches it only via
`retrieveSwathMaps`, which the streaming path never calls -- so **every run since streaming became
the default dropped its last `n mod kBatch` MS2 window spectra**. The `.d` reader emits window-group
major (the calibration patch's outer loop is over groups), so the loss lands entirely on the LAST
group: on D the last 105 cycles (~2.5 min) of two windows; at 4096, two whole windows plus 427
frames of two more -- 11% of the MS2 spectra, which is the 7.2% of peptides. The mzPeak reader
emits frame-major, so there the loss is the last ~9 cycles of every window instead (one more reason
the two readers' outputs differ). The differing slab digests of the 20 intact windows are the RT
axis: it is built from the loaded frames, so removing frames shifts every later `rt_index`.
Every pinned digest (D `b499027f`, TNBC 009 `057aad4b`) has 210 / ~170 spectra missing.

Fix: `finish()` called after both readers (commit 8219485). **Gate (node 1, build 01e455f):**

| | before (`soa2`) | after, batch 256 | after, batch 4096 |
|---|---:|---:|---:|
| windows / frames | 24; 21x1342, 1x1343, 2x1237 | 24; 22x1342, 2x1343 | identical to 256 |
| compact peaks | 1,255,770,798 | 1,258,053,975 | identical |
| spectra | 854,817 | **855,575 (+758)** | 855,575 |
| digest | `b499027f` | **`bc3c5ac9`** | **`bc3c5ac9` -- identical** |
| Sage peptides (rank 1, q <= 0.01) | 12,609 | **12,625 (+16)** | -- |
| wall / peak RSS | -- | 4:42 / 43.3 GB | 4:28 / 44.0 GB |

e2e 13/13. The batch size no longer changes a byte: the picker is per-frame pure and
`SPEXTRACTOR_PICK_BATCH` is a resource knob (4096: LOAD 48 -> 41 s on D, +0.7 GB peak). The
backlog item "the peak picker is sensitive to how frames are grouped" is closed as this bug, and
step 2d-i of the tiling plan (picker determinism under regrouping) with it: any batching is
equivalent. **New pinned digest for dataset D: `bc3c5ac9`**; TNBC 009 is re-pinned from the
probe's whole-run output below.

### 0c. `trace:max_span_sec` sweep -- shipped as CHANGES OUTPUT, never gated until now (node 2, `soa2`)

| `max_span_sec` | spectra | digest | Sage peptides (rank 1, q <= 0.01) | peak RSS | wall |
|---:|---:|---|---:|---:|---:|
| 0 (off) | 854,967 | `f435f4c2` | 12,619 | 51.60 GB | 2:56.6 |
| 60 | 854,030 | `443e448a` | 12,611 | 51.58 GB | 2:43.8 |
| **120 (default)** | 854,817 | `b499027f` | 12,609 | 50.70 GB | 2:45.0 |
| 240 | 855,143 | `94f899f7` | 12,619 | 51.32 GB | 2:52.9 |

Spread 10 peptides = 0.08%, below the set-noise of this gate: the trim is free in peptides, and
so is the halo it sizes. D's MS1 spans are longer than TNBC's (median 15.2 s, p90 59.6, p99 203.6,
max 1,243.9; 3.06% > 120 s, 0.66% > 240 s) -- the RT-local premise holds on both, with a heavier
tail here.

### 0d. One instrumented build (2afb302): G vs F, MS2 extents, per-window tof edges

- **`G == G_live == F` on every window** (`[gstat] F=1342 G=1342 G_live=1342`, all 24): every MS2
  frame of a window carries at least one fragment point. `buildFragStats`'s Pearson denominator can
  be injected from the frame table (`G_override = F_window`) with no output change -- step 2f is
  exact, and tiles can score without the whole window's fragments.
- **MS2 detection extents on the integer path** (never logged before; window 523-546, n=8.2 M):
  pre-trim span p50 = 5 s, p90 = 20, p99 = 55, max 345; **one-sided** (apex to the far end) p50 =
  5, p90 = 15, p99 = 45, max 330; fraction > 63 s = 0.27%, **> 123 s = 0.013%**, > 300 s = 2.4e-5.
  Post-trim: total span capped at 119 s, one-sided max 205 s (the asymmetric-trim case from review
  correction C1, real but 0.002%). So a 123 s halo mis-sees ~1 in 8,000 fragment traces at a tile
  edge and a 300 s halo ~1 in 40,000; the Sage gate prices what that costs.
- Window tof edges: `tlo` 0-3, `thi` 634,068-634,072 on every window -- the slab extremes ARE the
  acquisition range, so injecting `tlo/thi` from the method's m/z range changes nothing measurable.

### The parallel mzML writer wrote an empty index (v2.0.0 .. 2afb302)

Found by C4's "verify the shipped index" instruction: `patches/openms-mzml-parallel-write.patch`
concatenated the per-thread handlers' output but discarded their `spectra_offsets_`, so every
parallel-written file carried `<indexList count="0">` with one dummy entry -- still a well-formed
mzML, still read by Sage/MSFragger (which do not use the index), but not a valid indexedmzML, and
the e2e never checked it. Fixed (3c996c3) by rebasing each chunk's offsets by `os.tellp()` as it is
appended; validator (`step0/idxcheck.py`, reads every offset and checks it lands on
`<spectrum index="N"`) reports INDEX VALID on 854,817 offsets, spectrum digest unchanged
(`b499027f`), e2e 13/13. The validator is now part of the tail-fix gate chain.

## 2026-09-09: Fix A -- one `malloc_trim` when the last window has materialised: -31 GB on TNBC 009, wall-neutral

Step 0a showed the drained compact store is real free-list (the allocator tunables returned it,
at +26% wall). Fix A (58c0981) returns it at the one moment it dies: `toSlab` empties each window's
compact frames at window start, so when the LAST window has materialised the store is a corpse --
an atomic counter and one `malloc_trim(0)` under the existing `perf:malloc_trim` gate, ~6 lines.
Interleaved A/B on `data`, two reps each, 250 ms VmRSS sidecar, `MemAvailable` ~1 TB throughout:

| arm | rep | wall | peak RSS (`time -v`) | sidecar | digest |
|---|---|---:|---:|---:|---|
| `soa2` (HEAD before) | r1 | 10:54.7 | 168.93 GB | 169.07 | IDENTICAL |
| `soa2` | r2 | 11:39.9 | 169.17 GB | 169.23 | IDENTICAL |
| `fixA` | r1 | 11:15.8 | 139.23 GB | 139.33 | IDENTICAL |
| `fixA` | r2 | 10:49.4 | 136.91 GB | 136.83 | IDENTICAL |
| **mean** | | **11:17 -> 11:03 (-2%, noise)** | **169.05 -> 138.07 GB (-31.0 GB, -18.3%)** | | |

The trim itself: `101,532 -> 63,860 MB RSS in 4.4 s` at t~215 s (37.7 GB returned, 0.12 s/GB --
against the concurrency knob's 2.60 s/GB). Dataset D: 43.45 GB peak, digest `b499027f` (pre-tail
pin, identical), the trim returned 11.9 GB in 2.6 s. Wall is neutral within the two reps' spread
(the r1/r2 order of magnitude of node-load noise is 45 s). Published as arm `fixA`.

**What it re-prices.** The 2.85 GB/window exchange rate stands, but the baseline it applies to is
now 138 GB, not 169 (review correction C8: tiling is benchmarked against this). Step 2d-ii (the
loader builds the slabs directly, no second representation) makes the same saving structural --
no trim, no allocator dependence -- and is gated separately.

## 2026-09-09: tiling plan step 1 -- the RT-tile probe on TNBC 009 (no reader/writer code)

Review correction C7 restated step 1 as a probe: keep only the MS2 frames of one RT range plus a
halo, applied AFTER picking (batches unchanged); freeze the whole run's per-window tof edges and
Pearson support `G`; score only the precursors whose apex is in the core; MS1 stays whole-run
(pass 1 untiled). Arm `probe1` (e3a8917 = tail fix + Fix A + the env-gated probe), `data`, 100
threads, 250 ms VmRSS sidecar; 13 cores of 600 s at halo 123 s, the last core absorbing the
remainder; outputs concatenated by core (`bench/tile_concat.py`, self-check byte-exact on the whole
output) and searched once.

| run | loop peak RSS | loop start | loop wall | loop CPU | detect CPU | core precursors | spectra |
|---|---:|---:|---:|---:|---:|---:|---:|
| whole run (reference) | **138.81 GB** | 47.8 | 421.1 s | 34,811 s | 4,371 s | 5,445,436 | 3,919,816 |
| t0 0-600 s | 26.3 | 26.3 | 15.3 | 590 | 134 | 307,117 | 163,817 |
| t1 600-1200 | 23.3 | 9.5 | 37.7 | 1,561 | 322 | 501,647 | 325,387 |
| t2 1200-1800 | 26.4 | 9.8 | 36.8 | 1,892 | 392 | 510,645 | 357,867 |
| t3 1800-2400 | 26.2 | 9.2 | 46.9 | 2,051 | 432 | 499,855 | 358,782 |
| t4 2400-3000 | 28.8 | 12.1 | 46.2 | 2,415 | 530 | 502,970 | 371,347 |
| t5 3000-3600 | 28.0 | 16.9 | 48.5 | 2,507 | 566 | 489,708 | 373,525 |
| t6 3600-4200 | **29.8** | 10.7 | 51.5 | 2,682 | 588 | 484,063 | 378,118 |
| t7 4200-4800 | 28.2 | 10.2 | 50.4 | 2,419 | 533 | 449,372 | 353,260 |
| t8 4800-5400 | 28.7 | 10.6 | 50.1 | 2,502 | 548 | 418,722 | 327,660 |
| t9 5400-6000 | 28.4 | 12.3 | 45.8 | 2,323 | 492 | 391,518 | 301,048 |
| t10 6000-6600 | 27.6 | 17.5 | 40.0 | 2,027 | 415 | 345,945 | 250,137 |
| t11 6600-7200 | 26.5 | 11.2 | 34.2 | 1,763 | 338 | 300,919 | 202,190 |
| t12 7200-end | 21.3 | 10.3 | 25.3 | 1,263 | 234 | 242,955 | 156,658 |
| **13 tiles, sum** | max 29.8 (**-79%**) | | **528.7 s (+25.5%)** | **25,995 s (-25.3%)** | 5,523 s (+26.4%) | 5,445,436 (exact) | 3,919,796 (-20) |

Each tile holds 520 of the 4,797 MS2 frames per window (C+2H = 846 s / 1.626 s); the probe's
whole-process peak stays at 82 GB per tile because the load and the MS1 phase are still whole-run
(that is pass 1, step 5). Loop start ~10 GB = the persistent MS1 state + floor.

**What the probe settles.**
1. **Memory**: the window loop of a 600-s tile peaks at 21-30 GB against 138.8 GB whole -- the
   memory model's ~33 GB cell for C=600/H=123 (plus its 5-8 GB allowance) was conservative.
2. **CPU**: total loop CPU falls by a quarter under tiling despite the 41% frame halo. The detector
   itself is near-linear (+26% CPU for +41% frames: sub-linear); the saving is in everything after
   it -- merge, split, `buildFragStats` (per-fragment work over G frames), the candidate scans --
   which scale with the window's frames and fragments and are 9x smaller per tile: 30,441 -> 20,472 s
   (-33%). So the super-linearity was real, and it lives in the scoring path, not the detector.
3. **Wall**: +25.5% when the tiles run one after another, exactly the linear model's +23%: a tile's
   loop runs at ~49x against the whole run's 83x, because fewer windows are in flight. That is the
   occupancy the pipelined design (a tile's tail overlapped with the next tile's head) recovers;
   with the CPU already 25% lower, tiled wall should come in BELOW the whole run once two tiles
   overlap.
4. **Halo**: on the densest core (t6) and an early one (t1), H = 63 / 123 / 300 s:

   | | loop peak | loop CPU | loop wall | spectra |
   |---|---:|---:|---:|---:|
   | t6, H=63 | 26.7 GB | 2,352 s | 43.6 s | 378,110 |
   | t6, H=123 | 29.8 | 2,682 | 51.5 | 378,118 |
   | t6, H=300 | 35.7 | 3,441 | 71.4 | 378,128 |
   | t1, H=63 | 22.2 | 1,386 | 34.3 | 325,383 |
   | t1, H=123 | 23.3 | 1,561 | 37.7 | 325,387 |
   | t1, H=300 | 25.9 | 2,038 | 45.8 | 325,391 |

   The spectrum count moves by ~1e-5 across a 5x range of halo: the boundary is not where the
   spectra change. Memory and CPU price the halo at ~0.05 GB and ~4-6 CPU-s per second of halo on
   the densest tile; 123 s stays the default.
5. **Peptides (Sage, rank 1, q <= 0.01)**: tiled **27,931 vs whole 27,859 (+0.26%)**; set overlap
   26,083 common, 1,776 whole-only, 1,848 tiled-only (93.6%). For calibration, the max_span arms
   of step 0c -- a change to 0.1-0.6% of traces -- flip 30-112 of 12.6 k peptides (0.2-0.9%), and
   the project's 2026-09-02 thread-count comparison flipped 20% of the union at a 2% count
   difference: an FDR-thresholded set is a poor instrument for small perturbations. The per-spectrum comparison is the instrument, and it settled it:
**Per-spectrum attribution (`bench/mzml_specdiff.py`, whole vs tiled, 3,920,295 spectra aligned by
(RT, m/z, z)):** 0.29% identical, 61.0% differ in the m/z list, 38.7% in the peak count, 0.01% only
in one file -- and **uniformly in RT**: 99.7% of the spectra more than 123 s from any boundary
differ, the same as within 30 s of one. Dataset D, one tile [0, 900) at H = 123: 99.7% again,
including spectra 600-900 s from the only boundary, on a dataset whose longest post-split MS2 trace
spans 345 s. Per spectrum ~90% of the fragments are identical (same m/z, identical intensities) and
the rest are DIFFERENT TRACES -- so the split of traces, not their detection, depends on where the
window is cut. Cause, in the OpenMS source: `ElutionPeakDetection::detectElutionPeaks_` sets the
Savitzky-Golay window (and the local-extrema neighbourhood) from `chrom_fwhm /
mt.getAverageMS1CycleTime()`, and the tool hands EPD each parent's REAL points only, so a long
gappy parent (a background ion seen intermittently across the run) has an inflated cycle time, a
narrower window, and any sub-range of it -- what a tile sees -- is smoothed and split differently
from the whole. Parents split 1.25-1.30 ways on TNBC (20.8 M -> 26.0 M, 29.4 M -> 38.2 M on two
windows); children are short (TNBC, 1.22e9 traces: span p99 44 s, max 421 s; one-sided p99 35 s,
> 123 s 0.0034%), parents are not bounded at all. Fix (d85c2b4): EPD gains `scan_time` (0 = as
before) and the tool passes the MS level's median frame spacing (`trace:split_scan_time=frame`,
default; `trace` = the old behaviour), so the window is the same for every trace and the split is
local. External review (internal review record): the patch and plumbing are
correct and `trace` reproduces the old arithmetic; two defects fixed in 86444a5 (an upper median
of the gaps is 2x the cycle when most frames of a sparse window were skipped -> smallest positive
gap; `ceil(7 / 1.4)` sits on an integer and flips with the last bit -> backed off by 1e-6 on the
fixed path). Two findings that stand: (a) EPD's valley DECISIONS compare maxima across the whole
parent (ElutionPeakDetection.cpp:278-294, the retained-maximum logic), so a fixed window makes the
smoothing local but not the decision -- a counterexample splits the whole trace and not its
5-point-shorter crop with identical interior smoothing; full tile invariance therefore needs
BOUNDED PARENTS (a cap on the detector's extension around the seed), which is a design decision;
(b) a fixed point count is a fixed time span only on gap-free traces (SG convolves by index), so
the frame-spacing window over-smooths sparse traces and can merge two separable peaks that the
per-trace window kept apart -- the default is an empirical change to be judged on peptides, not a
neutral correction. Gates: the `trace` arm must reproduce `bc3c5ac9`; the `frame` default is
Sage-gated on D against 12,625; then the probe's tile t6 against the whole run under the new
default.
## 2026-09-09: tiling plan step 2 -- the output-identical refactors, each D-gated

Every arm below reproduces the tail-fix pin on dataset D (`bc3c5ac9`, `SPECTRUM DATA IDENTICAL`)
and passes the e2e suite; the reference is the `tail_S30` run's `pseudo.mzML`
(archived as `share/refs/S30_tail_bc3c5ac9.mzML`).

| step | commit | what changed | gate result |
|---|---|---|---|
| 2b/2c mzPeak `frame=` + flight-time axis from the archive's embedded tdf | fdb33b6 | `frameIdOf()` returned 0 on every mzPeak spectrum, so the integer detector never ran on that input; the exact-m/z path's calibration now feeds `setTofAxis()` | `.d` identical; mzPeak (dataset D v0.9.2 archive + sidecar tdf): **integer detector runs for the first time**, 855,561 spectra, Sage **12,537 vs 12,625 from the .d (-0.7%)**, set overlap 11,499 common / 1,126 d-only / 1,038 mz-only; the OpenMS detector on the same archive gives **12,756** (11,280 common with the .d) -- on this file the detectors' sign is the one the integer-detector memory records (openms +1.7%), and the archive's peaks differ from the .d's (the two detectors' outputs are not comparable peptide for peptide). LOAD 131 s / 63 GB peak against 45 s / 43 GB (the library decode). The July dataset A archive carries no per-window IM band and never streamed. |
| 2e `TileWriter` at one block | 3856893, 371653d | `PlainMSDataWritingConsumer` subclass: header once, blocks appended in canonical order, per-thread encoders, index offsets rebased, exact count when known, same-width placeholder otherwise; native IDs renewed as `spectrum=<rank>` (what the bulk writer emitted for our id-less spectra) | identical; **WRITE 19.6 -> 13.2 s (11.7x -> 29.4x)**; header differs only in completion time and the `out` name; INDEX VALID (855,575 offsets); e2e check 14 (count/id/index/indexList) added, 0f534bd |
| 2g content-only MS1 trace order | 7cd120f | `(mz, rt, im, intensity, npts, len, first frame time)`, equal-key neighbours logged | **0 equal-key neighbours** on D; identical |
| 2d-ii the loader builds the slabs directly | 98b34d5, eb04568, a77e089, 9c97b60 | `CompactFrame` store and `toSlab` deleted; frames appended to the window's `PeakSlab` as picked, m/z quanta -> bins in place at window start; one reusable per-batch scratch; 2.5x virtual over-reservation, no shrink; per-window append in parallel | identical, **24/24 per-window slab digests identical**; LOAD on D: 45 s (old) -> 80 (1.25x growth, per-frame temporaries) -> 65 (projected reserve + shrink) -> 59 (scratch) -> see #4; RSS at load end 31.1 -> 40.3 GB (11.7 GB of arena free-list appears during the load and is returned after the slab pass: 39.0 -> 28.3 GB); D peak 43.3 -> 46.2 GB (load-phase, below the loop peak on the 2-h file) |

Not done: 2a (the frame table comes with the tiled reader), 2f (exact by the metadata check --
frames per window group == non-empty frames on both datasets -- and lands with the reader).

2d-ii's remaining cost, measured on the same node in one chain (run #4, 9c97b60): LOAD 58.4 vs
46.3 s, window loop 188.9 vs 184.7 s wall and 12,001 vs 12,019 CPU-s (equal), total 291 vs 272 s.
The node has one NUMA domain, so placement is not it; the parallel per-window append changed
nothing (pick/flush 44.0 vs 45.0 s). What differs: the old store's per-frame vectors were
allocated -- and their pages first-touched -- by the 100 picking threads, while the direct slabs
are first-touched by at most 24 appending threads (one per window) or one; 12 GB of fresh pages
is ~3 M faults. A default-initialising allocator for the three slab arrays would let the copy run
100-wide (backlog). The 11.7-13 GB of arena free-list that appears during the load and goes with
the trim after the slab pass is not explained.

## 2026-09-09: the valley splitter's scan time -- gate on D and TNBC 009, and it does NOT make tiles invariant

`trace:split_scan_time` (d85c2b4, 86444a5, 3755121; EPD param `scan_time`): the smoothing window
of ElutionPeakDetection is `chrom_fwhm / scan time`, and the scan time was each parent's own
average cycle time over its real points -- inflated by gaps, so gappy parents got a narrower window
and any sub-range of one split differently. `frame` (new default) uses one run-wide value, the
smallest gap between consecutive MS1 frame times; `trace` is the old behaviour. Build 0796595
(arm `split1`), node 1 / node 3:

| | dataset D, `trace` | dataset D, `frame` | TNBC 009, `frame` (probe2 whole) |
|---|---:|---:|---:|
| digest | **`bc3c5ac9` (identical to the pin)** | differs | `dba82377` |
| MS1 traces | 4,905,492 | **3,092,075 (-37%)** | (precursors 4,433,660 vs 5,445,436, -19%) |
| spectra | 855,575 | **612,664 (-28%)** | **3,291,533 (-16%)** |
| Sage peptides (rank 1, q <= 0.01) | 12,625 | **13,388 (+6.0%)** | **28,845 (+3.5%)** vs 27,859 |
| set overlap with the old arm | | 11,564 common / 1,061 old-only / 1,824 new-only | 25,304 common |
| wall / peak | 4:35 / 46.9 GB | 4:02 / 46.6 GB | 11:16 / 130.2 GB |

So the fixed window is a large whole-run change, not a neutral correction (the reviews said so):
it splits far less -- most MS1 traces from MassTraceDetection have gaps, their per-trace windows
were narrow, and the run-wide window is wider -- and the fewer, longer traces identify MORE
peptides on both files. Whether the +6% / +3.5% is real signal or an FDR-calibration artefact of
a smaller search space is what the entrapment search and the second engine decide. **Entrapment
(Sage, human + Arabidopsis, dataset D, `bench/entrap_apply.py`, peptide-hypothesis ratio 0.6805):**

| arm | target peptides | entrapment hits | raw | FDR at nominal 1% | 95% CI |
|---|---:|---:|---:|---:|---|
| `trace` (old) | 11,722 | 108 | 0.91% | 1.35% | 1.11-1.62% |
| `frame` (new) | **12,496 (+6.6%)** | 99 | 0.79% | **1.16%** | 0.94-1.40% |

More peptides at a lower measured error rate: the gain is signal, not recalibration (the old arm's
1.35% matches the 1.33% recorded for D at the shipped defaults on 2026-09-05). **MSFragger 4.4.1
(second engine, `score_msf_td.py`, D): 12,236 -> 12,732 target peptides at 1% FDR (+4.1%).** Both
engines up, entrapment down: `frame` is the default on the development branch; the second
acquisition's MSFragger run and the entrapment on TNBC 009 are due before release.

**Tile invariance: no.** probe2 (arm `split1`, TNBC 009, core t6 3600-4200 s at H = 123 against
the whole run under the same default): 314,685 core spectra, 379 identical (0.1%), 210,822 differ
in the m/z list, 103,361 in the peak count -- the same 99.9% as with the old window, and uniform
in RT. As the external review predicted: the smoothing is now local, the DECISION is not (EPD's
retained-maximum walk), and detection itself is not (seeds, running centroid, `visited`). The
fragment-level number (the spectrum-level one saturates at one changed fragment): over the
314,685 core spectra present in both, **80.7 M fragments have the identical m/z in both, 46.4 M
exist only in the whole run and 46.2 M only in the tile -- 46.6% of the union agree exactly.**
That exact test also counts an apex-frame shift of the SAME fragment (the apex m/z is computed
from the apex frame's calibration) as one lost plus one gained; the 10-ppm-tolerant match and the
same metric on probe1's outputs (the old per-trace window) are recorded when they land.

## 2026-09-09: cause 2 vs cause 3 -- detection alone accounts for the larger half of the tile churn

Judge's open question 1 (internal review record): the whole run and the
core t6 probe (3600-4200 s, H = 123, frozen edges/G, MS1 whole-run) with `-trace:ms2_split_valleys
0` -- the valley splitter off, so the fragment traces ARE the detector's parents -- dumped per
window before trimming (`SPEXTRACTOR_FRAG_DUMP`, build 9981e97, node 1), symmetric difference
of the parents whose RT lies in the core, key (mz 1e-4, apex RT, im 1e-4, apex intensity):

| | parents in the core | whole-only | tile-only | common with a different extent |
|---|---:|---:|---:|---:|
| all 28 windows | **62.6 M common = 58.8% of the union** | 21.8 M | 22.0 M | 31,018 (0.05%) |
| 400-426 Th (sparsest) | 500,186 common | 25,325 (5%) | 26,264 | 16 |
| 1000-1026 Th | 2,611,394 | 960,254 (27%) | 965,755 | 1,582 |
| 750-951 Th (densest) | 2.8-3.4 M each | 27-28% each | | |

So with NO splitting involved, 41% of the parents a 600-s core sees differ from the whole run's --
detection (seeds in intensity order, the running centroid/sd, `visited` claims along connected
lanes, unbounded extension) is the larger half of the churn, and it scales with peak density; the
parents both runs share are the same objects (extents equal in 99.95%). A splitter-only fix (the
review's LPS option) therefore cannot deliver invariance, and the seed cap's residual would sit on
exactly these dense lanes. Ranking of the mechanisms for the decision: the cell-atomic grid makes
this churn part of the definition (parents are detected per cell in every mode, so it is the same
in the whole run and in any tiling); the halo designs fight it and lose.

## 2026-09-09: the cell-atomic grid, step 2 -- per-cell detection is exact on real data

Decision A (internal review record): the integer detector traces each
window per fixed RT cell (`tile:rt_sec`, cut at the run's MS1 frame times, stamped in the header),
so that any grouping of cells into tiles reproduces the whole run bit for bit. Commits 8779a0b,
fc88b82 (review fixes), ee49bda. Reviews (internal review record; two external reviews and a
four-lens refutation workflow: "sound-with-fixes", the fixes are in fc88b82).

**Gates (node 1, fc88b82):**

| gate | result |
|---|---|
| identity, `tile:rt_sec=-1`, `trace` window | `bc3c5ac9` -- **IDENTICAL** to the tail-fix pin |
| identity, `-1`, `frame` window (default) | **IDENTICAL** to the split arm's output |
| e2e | **15/15** -- check 15 holds 1 vs 4 threads byte-identical at a pitch that cuts the fixture in two |
| whole run at pitch 300 (6 cells) | 611,034 spectra (612,664 at -1); 5:04 / 47.3 GB |
| **cell exactness**, core cell 2 [602.0, 902.6) with cells 1-3 resident (probe, halo 0) vs the whole run at 300 | 152,360 spectra: **IDENTICAL**, 66.85 M fragments matched 100%, intensities 100% |
| cell exactness, core cell 3 [902.6, 1203.2) with cells 2-4 resident | 185,782 spectra: **IDENTICAL**, 84.32 M fragments 100% |

The same instrument gave 0.1% identical spectra and 51% matched fragments for the halo design
yesterday. The probe proves the MS2-detection half per cell (the review's M3): scoring support,
carry-over and pass 1 per cell are the tiled reader's contract (plan file, "tiled reader's
contract"), gated on the reader itself.

**The science price so far (arm `cell1` = 8779a0b, `frame` window):**

| | -1 (whole) | pitch 600 | |
|---|---:|---:|---|
| TNBC 009 spectra | 3,291,533 | 3,285,827 | -0.2% |
| TNBC 009 Sage peptides | 28,845 | **29,066** | +0.8%; set 27,058 common / 1,787 lost / 2,008 gained |
| within 15 s of a cut line (1,435 base peptides) | | lost 135 (9.4%) vs 6.0% elsewhere; gained 105 | the excess ~50 peptides = 0.17% of the set |
| TNBC 009 peak / wall | 130.2 GB / 11:16 | 124.5 GB / 11:37 | |
| **one 600-s cell (t6), no halo** (G3, sidecar) | | **loop peak 23.7 GB** (start 12.0), loop 37 s / 1,856 CPU-s, detect 446 s | vs 29.8 GB / 51 s / 2,682 / 588 with the 123-s halo |

The set churn (6%) is the FDR-set flip every perturbation shows; the cut lines add a small,
localised excess that the gains next to the lines nearly balance.

**Dataset D at -1 / 600 / 300 (arm `cell1`, node 1, Sage `sage_deiso`, entrapment
human+Arabidopsis, MSFragger 4.4.1 via `score_msf_td.py`):**

| pitch | cells | spectra | Sage peptides | entrapment FDR @ nominal 1% (95% CI) | MSFragger | digest |
|---:|---:|---:|---:|---|---:|---|
| -1 | 1 | 612,664 | 13,387 | 1.16% (0.95-1.40); 12,496 targets / 99 hits | 12,732 | `63d0ecfc` |
| **600** | 3 | 611,850 | **13,490 (+0.8%)** | **1.30% (1.05-1.55)**; 12,630 / 112 | **12,918 (+1.5%)** | **`f9e65631`** |
| 300 | 6 | 611,034 | 13,836 (+3.4%) | 1.38% (1.12-1.63); 12,918 / 121 | 12,921 (+1.5%) | `3f8d421d` |

Loss within 15 s of a cut line: 9.3% (600) / 9.6% (300) of the base peptides there, against
7.9% / 7.7% elsewhere -- the same ~2-point excess as on TNBC. More cuts give more Sage peptides
while the measured error rate drifts up (raw entrapment fraction 0.79 -> 0.88 -> 0.93%), so part
of the pitch-300 gain is calibration; at 600 both engines rise and the entrapment moves +0.14
points, inside its interval and the release rule (<= base + 0.3). **`tile:rt_sec` defaults to 600
from 2026-09-10** (the decision's pitch); new pins at the default: **D `f9e65631`**, **TNBC 009
`f0502f33`** (29,066 Sage). Gate on the default flip (58c6784, arm `grid1`): D at the default
digests `f9e65631` exactly (the review fixes are byte-neutral at 600 as at -1), `-1` reproduces
both window-mode pins, e2e 15/15.

**Second-acquisition release gate (arm `cell1`, TNBC 009, `data`, 2026-09-10 01:15):**

| | -1 (whole) | pitch 600 | rule |
|---|---:|---:|---|
| Sage peptides | 28,845 | 29,066 (+0.8%) | |
| entrapment FDR @ nominal 1% (95% CI) | 0.92% (0.79-1.07); 26,761 targets / 168 hits | 1.08% (0.93-1.23); 27,238 / 201 | <= base + 0.3 points: **+0.16, passes** |
| MSFragger 4.4.1 (`score_msf_td.py`) | 27,973 | 27,952 (-0.08%) | flat -- 21 peptides, inside the run-to-run spread |
| peak RSS / wall (600) | | 124.6 GB / 11:10 | |

So on the second acquisition the grid is a wash for MSFragger and +0.8% for Sage at +0.16
entrapment points; on D both engines rose. The default stands; nothing here argues for a
smaller pitch (see the pitch-300 calibration drift above).

**Fragment agreement under the old tiling, for the record** (probe1, old window, 13 tiles with a
123-s halo, `bench/mzml_specdiff.py` at HEAD): 51.6% of fragments matched within 10 ppm, 56.8%
of those with the same intensity -- the same picture as probe2 (51.0% / 51.8%); the churn was
never about the splitter window.


### Cell grid, step 1: frozen frame tables (2026-09-10)

The tiled reader's contract (internal review record, step 1) needs
every window's frame table -- the window store's frame sequence, the Pearson support `G`, and the
MS1 -> window `nearest_local` map -- to be a RUN-LEVEL constant, so a tile that holds some of the
frames scores exactly as the whole run does. Two commits:

- **1 (a), 0862aa8:** `FrozenWin` per window (frame times, vendor ids, RT-axis indices, per-frame
  calibration factors, `nearest_local` over ALL window frames, ties to the later frame); the
  window store's table is the frozen one, `buildFragStats` takes `G = F` and the table's map (the
  `used` policy and `G_override` are gone: G == G_live == F on 24/24 and 28/28 windows); the cell
  cut is slab-local with a value lookup of each cell's first frame time in the table (the
  ee49bda offset-rebase throw is deleted). Source of the table: the resident slab.
- **1 (b), 38907fb:** the table comes from the READER's metadata. The calibration patch exports
  `BrukerTimsFile::lastDIAFrameTable()` (per window group: frame id, `frame.time`, raw peak
  count, in the loop's own order) and `lastDIAWindows()`; the mzPeak sweep is cached as
  `MzPeakRunMeta` (with `number_of_peaks` when the archive records it); the consumer walks each
  window's delivered frames in LOCKSTEP with its table and appends every skipped table frame as
  an EMPTY frame (a raw frame with no peaks is not in the table, matching the loop's skip; an
  empty ion-mobility slice, or a probe-dropped frame, is padded), counted per window
  (`[frames] window W: delivered N, padded P`), the slab asserted equal to its table afterwards;
  header stamp `spx:frame_table=loader|mzpeak|slab` (mzML stays the slab: internal backlog).

**Gate 1 (a) (node 1, 0862aa8, chain `step1a`):** D at the default (600) `f9e65631` MATCH;
`-1` in both window modes IDENTICAL to their pins (`frame` = the split arm's output, `trace` =
`bc3c5ac9`); 300 `3f8d421d` MATCH; e2e 15/15 (check 15 now really runs the integer detector:
a773f6e). **The cell-exactness probes DIFFER at 1 (a), as they must:** the frozen table is the
resident slab's, and under the probe the slab holds three of six cells, so every window scores
on `F = G = 651` frames instead of 1,342 (`[gstat]`: 1,342/1,342 whole, 651/651 probe) -- the
Pearson mean and norm of every fragment move. 1 (a) therefore holds the whole-run pins but is
not yet tile-exact; that is exactly the gap 1 (b) closes (the probe's dropped frames are padded,
`F` is the run-level 1,342 again). The `SPEXTRACTOR_WINDOW_FILE` row's `G` is no longer read
(only `tlo/thi`); the probe's support now comes from the table.

**Gate 1 (b) (node 1, 38907fb, chain `step1b`, 01:17-02:13; the fix-up a91cfcf re-gated
below):**

| arm | `[frames]` | result |
|---|---|---|
| D at the default (600) | source=loader, 24 windows, 32,210 delivered, **0 padded** | `f9e65631` **MATCH**; 611,850 spectra, 3:57 / 47.1 GB |
| D `-1`, `frame` window | 0 padded | **IDENTICAL** to the split arm's output (612,664) |
| D `-1`, `trace` window | 0 padded | **IDENTICAL** to `bc3c5ac9` (855,575) |
| D at 300 | 0 padded | `3f8d421d` **MATCH** |
| cell exactness, core 2 (cells 1-3 resident) | 15,624 delivered, **16,586 padded** (the probe's dropped frames) | **IDENTICAL** to the whole run at 300 |
| cell exactness, core 3 (cells 2-4 resident) | 15,624 / 16,586 | **IDENTICAL** |
| mzPeak dataset D archive (v0.9.2 + sidecar tdf), new build vs the step-1a binary, default settings | source=mzpeak, 0 padded; the archive records `number_of_peaks`, so frames with none are skipped like the .d loop's | **IDENTICAL**; 611,445 spectra both (5:17 / 61.8 GB) |
| TNBC 009 at the default | source=loader, 28 windows, 134,312 delivered, **0 padded** | `f0502f33` **MATCH**; 3,285,827 spectra, 16:47 / 124.1 GB (a loaded node) |
| e2e | | **15/15** |

So the metadata table IS the delivered table on both acquisitions and on the archive (0 pads,
as the 2f metadata check predicted), today's output is byte-unchanged in every mode, and the
cell-exactness probe is restored: with the dropped frames padded, every window scores on the
run-level `F = 1,342` again and the core cells reproduce the whole run bit for bit. This is the
first build on which "only some of a window's frames are resident" gives the whole run's output
-- the property the tiled reader is built on. Review (two external reviewers; internal review record):
sound at P = 0; the structural throws had to be made fatal (the loader's fallback catch would
have appended a second full load onto the partial slabs) and the lockstep now matches (rt, frame
id) -- a91cfcf. **The fix-up re-gated (chain `step1b2`, 02:15-03:05): the same nine arms, the
same results -- `f9e65631` / IDENTICAL / IDENTICAL / `3f8d421d` / probes IDENTICAL / mzPeak
IDENTICAL / `f0502f33` / 15/15, 0 padded on 24/24 and 28/28.** Step 1 is closed; the dev build on
node 1 is a91cfcf's.


### Cell grid, step 2: the tile loop on the resident source (2026-09-10)

Commits 1e7dcb1 (the loop), 5570d25 and ad897d0 (the review fixes; internal review record),
c4997d1 (the gate tools ship with every deploy). `tile:cells_per_tile` groups the cells into
tiles; per tile the window loop's parallel region runs on the tile's frames (a slice of the
resident slabs; run-level band edges from the whole slab), owns the precursors decided with the
scorer's own arithmetic (the first tile k with (T_{k+1} - rt) > delta_rt), scores them with the
previous tile's boundary fragments carried over (rtOf >= T_k - 2 delta_rt - 1e-6 s, over own and
carried-in, so transitive), sorts the tile's spectra and writes them as one block through the
streaming writer (exact count with one tile; a same-width placeholder patched at the end
otherwise, which the digest tools normalise). 0 = one tile = today's loop statement for statement.

**Gate (node 1, ad897d0, chain `step2c`, 03:50-05:09; every digest with the current
`semantic_digest.py`):**

| arm | threads | wall | peak RSS | result |
|---|---:|---:|---:|---|
| D default (600), one tile | 100 | 3:19 | 47.1 GB | `f9e65631` **MATCH**, exact count |
| D `-1`, `frame` / `trace` window | 100 | 3:25 / 3:51 | 47.4 / 46.8 GB | **IDENTICAL** to both pins |
| D 300, one tile | 100 | 3:16 | 47.0 GB | `3f8d421d` **MATCH** |
| D 300 as **2 tiles** (3 cells each) | 8 / 100 | 11:59 / 3:46 | 38.8 / 47.6 GB | `3f8d421d` **MATCH** both; 24 windows carried |
| D 300 as **3 tiles** | 8 / 100 | 11:29 / 3:42 | 36.2 / 47.1 GB | `3f8d421d` **MATCH** both; 48 carry lines |
| D 300 as **6 tiles** (one cell each) | 8 / 100 | 10:52 / 3:45 | 36.0 / 47.2 GB | `3f8d421d` **MATCH** both; 120 carry lines |
| Sage on the 6-tile D-300 file | | | | **13,836 peptides** = the whole-run file's (the zero-padded count is accepted) |
| TNBC 009 at 600 as **13 tiles** | 100 | 11:59 | 125.0 GB | `f0502f33` **MATCH**; 3,285,827 spectra, count/index invariants OK |
| e2e | | | | **18/18** (16: 4 tiles on a 26-cycle fixture with the two precursors in different tiles, 1 and 4 threads, exact one-tile count, check 14's invariants on the tiled file; 17: `SPEXTRACTOR_TILE_NO_CARRY` changes the output and pitch 5 s < 2 delta_rt reproduces the one-tile run; 18: a window with no frames in a middle tile) |

Diagnostic pair before the chain (`diag2`): whole vs 2 tiles at 300, `bench/mzml_specdiff.py`:
611,034 of 611,034 spectra identical, 266.4 M fragments 100% matched with identical intensities,
0 differing at any distance to the boundary. The one `DIFFERS` in the first run of this chain
was the deployed digest tool (the 09-06 copy on shared storage, without the count normalisation; memory
`gate-tools-drift-at-deployment`); the deploy script ships the bench tools now.

Reading: the tile loop is exact on real data for any tiling of the cells and any thread count,
on both acquisitions, and a search engine reads the tiled file as the same file. The resident
source holds the whole run, so memory is not the claim here (the 8-thread arms are lower only
because fewer windows are in flight); the streaming readers (steps 3b/4) turn the tile into the
memory unit. Wall at 100 threads: +13-15% for 2-6 tiles on D (per-tile barrier tails and the
per-tile sort/write), the price the pipeline (step 5) is meant to recover.


### Cell grid, step 3a: the band partition from metadata (2026-09-10)

Commits 1641a00 (`trace:band_edges=slab|acquisition`), dde035e and 09947f5 (the review fixes;
internal review record). The integer detector traces each window in flight-time
bands, uniform in [tlo, thi]; a seed's band decides which `visited` set competes for it. Until
now tlo/thi were the window's resident peaks' extremes -- data a reader holding only a tile's
frames cannot see. `acquisition` takes tlo = 0 and thi = the larger of the digitizer's quantised
last bin (`DigitizerNumSamples - 1` through the store's 1e-5 Da quantum and back, per frame) and
the MzAcqRangeUpper conversion -- run-level metadata. Measured: on D the two coincide (634,072;
slab thi 634,068-634,072, tlo 0-3); on TNBC they do NOT: MzAcqRangeUpper converts to 399,199
while the digitizer has 399,967 bins and real peaks reach 399,965 -- the first version's TNBC
arm threw on the superset assert, exactly the hazard both reviewers named, and the digitizer
bound is what makes the option safe.

**Science gate, dataset D at the default grid (600), same build (chain `step3a` on 1641a00, whose
D edges equal the fix-up's; the fix-up re-checked in chain `step3afu`):**

| | `slab` (today) | `acquisition` | |
|---|---:|---:|---|
| digest | `f9e65631` (the pin) | **`8c1b047f`** | the re-pin |
| Sage peptides | 13,490 | 13,491 | set: 13,488 common, 2 slab-only, 3 acquisition-only |
| entrapment FDR @ nominal 1% | 1.30% (12,630 / 112) | 1.30% (12,629 / 112) | |
| MSFragger | 12,918 | 12,918 | |
| per-spectrum diff (`mzml_specdiff.py`, same precursors: 0 only-A/B) | | | 98.5% of 611,850 spectra identical; 0.99% with a different fragment count, 0.55% intensity-only; fragments 99.98% common |

The interior edges move by a bin or two (a uniform partition over [0, 634,072] instead of
[tlo, thi]), which re-assigns the seeds nearest each edge between two `visited` sets -- 1.5% of
spectra touched at the fragment level, five peptides in 13,490 at the search level, identical
error rate.

**Second acquisition, TNBC 009 at 600 (chain `step3afu` on the fix-up 09947f5; the slab control
re-run on the same build digests `f0502f33`, so the recorded slab searches are the base):**

| | `slab` (recorded, same spectra) | `acquisition` | rule |
|---|---:|---:|---|
| digest | `f0502f33` | **`bff54a1f`** | the re-pin |
| Sage peptides | 29,066 | 29,067 | |
| entrapment FDR @ nominal 1% | 1.08% (27,238 / 201) | 1.09% (27,237 / 202) | <= base + 0.3 points: **+0.01, passes** |
| MSFragger | 27,952 | 27,949 | -3 |
| wall / peak RSS (100 threads) | 13:12 / 125.0 GB | 13:14 / 124.9 GB | |

**Decision (2026-09-10 07:05): `trace:band_edges` defaults to `acquisition`.** The edges a
reader holding only a tile's frames can compute cost nothing measurable on either acquisition
(+1 / +1 Sage, 0 / -3 MSFragger, +0.00 / +0.01 entrapment points), and they remove the last
data-dependent input from the detector's partition. `slab` stays as the option that reproduces
the earlier pins. New pins at the default: **D `8c1b047f`**, **TNBC 009 `bff54a1f`**; the `-1`
and 300 pins are re-recorded on the flipped build below. Review: internal review record
(three rounds; the first version's TNBC arm threw on the superset assert because MzAcqRangeUpper
converts to 399,199 while the digitizer's last bin is 399,966 -- the reviewers' hazard, measured).

**Pins at the flipped defaults (chain `repin` on 7e71998, node 1, 07:06-07:40; the reference
files are kept on node 1 (`pins_acq/`) for the reader gates):**

| arm | digest |
|---|---|
| D at the default (600, acquisition edges) | **`8c1b047f`** (the value the science gate searched) |
| D `-1`, `frame` window | `5eb25eb0` |
| D `-1`, `trace` window | `20afff0b` |
| D at 300 | `8fde0944` |
| D at 600 with `-trace:band_edges slab` | `f9e65631` -- the pre-flip path, intact |
| TNBC 009 at the default | **`bff54a1f`** |

The e2e section of this chain (and of the two before it) came back EMPTY: with the flipped default
every integer-detector check refused, because the synthetic `analysis.tdf` fixture carried no
GlobalMetadata and `acquisition` needs it -- and the failure escaped `check()` as a traceback, so
a chain grepping for "FAIL" or "checks passed" saw nothing at all. Both halves are fixed (8e2f1c0):
the fixture writes the bounds a real tdf always has, and a suite that cannot run now prints
`0/0 checks passed; FAIL ...`. **A silent gate section is a failed gate, not a passed one.**


### Where the memory and the time actually go, measured per data structure (2026-09-10)

`SPEXTRACTOR_LEDGER=<file>` samples live bytes per data-structure category every 250 ms next to
VmRSS, the glibc arena figures and a phase marker; `bench/ledger_report.py` reads the trace.
Charges are atomic adds at allocation boundaries (never in a per-peak loop) and count TOUCHED
bytes, so the 2.5x slab over-reservation does not inflate them. Four arms on dataset D, node 2,
100 threads, 1 TB free, every arm digest-identical to its pin (chain `instr2`, binary
906e758+ledger):

| arm | wall | peak RSS | where the peak is |
|---|---:|---:|---|
| resident, 1 tile (the pre-3b shape) | 2:34 | 47.2 GB | the LOAD: MS1 map 13.8 + slabs 12.0 + raw batches 3.8 GB charged, ~17 GB unattributed |
| streaming, 1 tile | 2:27 | 41.1 GB | the window loop |
| streaming, 3 tiles | 2:57 | 31.8 GB | **the MS1 pass** |
| streaming, 6 tiles | 2:48 | 31.8 GB | **the MS1 pass**, t = 9 s, before any MS2 frame is read |

**The peak has moved off the window loop.** In the 6-tile arm the MS2 side never exceeds 23 GB of
RSS and 3.8 GB of live bytes per tile, while the MS1 pass peaks at 31.0 GB with 13.9 GB charged
(map 10.5, raw frame batches 3.4) -- 17.1 GB unattributed, of which the glibc free list is 5.6.
That is why 3 and 6 tiles have the SAME peak: the MS2 side is already under the MS1 floor, and
tiling it further changes nothing until pass 1 is tiled as well.

**The window loop's structures, 6-tile arm, per tile:** slab 3.5 GB, seed order + liveness + the
flight-time index 1.5, emitted spectra 1.1, window arena 0.9, trace records 0.8, band arenas 0.5,
scoring auxiliaries 0.3, the boundary carry 0.17. The slab is live for 34.5% of the run, the MS1
map for 14.7%; nothing else clears 5% of the peak for any measurable time.

**Allocator retention is the second finding.** Across the six tiles the charged bytes stay between
0.5 and 3.8 GB, but RSS climbs monotonically 3.3 -> 9.8 -> 16.2 -> 20.7 -> 22.8 -> 21.8 GB and only
falls at the final write. Freed tile memory is not returned; `perf:malloc_trim` fires at phase
boundaries, not per tile. On D this is hidden under the MS1 floor; on the 2-hour file the MS2 side
is ~4x larger and it would not be.

**The tiles are wildly unequal, and the grouping is free.** At pitch 300 the per-tile window loops
took 0.8, 12.6, 28.4, 25.8, 28.9 and 15.6 s, and at pitch 300 grouped in twos the two tiles held
3.30e8 and 9.28e8 picked peaks. Peak memory and wall are set by the densest tile, and because the
output is identical for ANY grouping of cells, cells can be grouped by cumulative peak count --
which the frame table's `num_peaks` gives before anything is read.

**Reads are small next to compute.** Per tile at 6 tiles: read 1.3-3.6 s, window loop 12.6-28.9 s,
write 0.5-2.0 s. A prefetch thread can hide the read; it cannot hide much else.


### What tiling actually costs, and why (2026-09-10)

The controlled pair (dataset D, pitch 600, same digest `8c1b047f`): one tile runs the window loop
in 89.5 s at 67.6x parallelism; three tiles run it in 121.7 s at 43.1x -- **+32.2 s of occupancy
loss while the CPU total FALLS 13%** (5,990 -> 5,207 CPU-s). Tiling does not make the work bigger;
it makes the machine idler.

The cause is a straggler, measured from the tool's own per-window stage lines: in tile 2 the
slowest window takes 61.2 s against a 22.5 s median, and the tile cannot close until it finishes.
Summed over the three tiles the straggler gap is **71.7 s**. The culprit is the acquisition's own
geometry -- one isolation window is **141.3 Da wide** against a ~26 Da typical, and carries 71,294
of the tile's spectra with 51.0 s of tracing.

Two consequences. First, a prefetch thread is not the lever: all the I/O in a 3-tile run is 14.1 s
and the barrier costs 32.2. Second, the barrier is an artefact of the loop structure rather than a
data dependency -- **the carry is keyed by WINDOW** (`carry_out.at(wl[wi].first)`), so window w of
tile k+1 depends only on window w of tile k: 24 independent chains, not a barrier. A per-window
ready flag over a global (window, tile) work list removes the edge for zero redundant work, which
is what makes the redundant-boundary-cell idea (56% extra peaks re-read) unnecessary on one
machine.


### Cell grid, step 3b: the .d streaming source, gated (2026-09-10)

Commits 933bba5 (the two-pass source) and dbcfa9c (the review fixes; internal review record).
Pass 1 reads MS1 only and the loader exports the whole DIA frame table anyway; pass 2 calls the
patched loader once per tile over `Config::dia_ms2_rt_lo/hi`, into a fresh consumer lockstep-padded
against the tile's slice of the frozen tables. `SPEXTRACTOR_TILE_RESIDENT=1` keeps the resident
source for a same-build control.

**Gate (node 1, 54a7721, chain `step3b`, 07:40-09:23): 14 arms, every one SPECTRUM DATA
IDENTICAL to its pin, none differing.** G0 one tile on D at 600 / 300 / -1 `trace` and on TNBC at
600; G1 D at 300 as 2, 3 and 6 tiles at 8 and 100 threads; the resident control on the same build;
TNBC at 600 as 13 tiles. Every arm reported its frames: 32,210 window frames delivered on D and
134,312 on TNBC, 0 padded, 0 absent, and `[tile] coverage` accounted every precursor-window pair
exactly once (710,356 on D at one tile; 4,206,205 on TNBC over 13 tiles).

**What it costs and what it saves:**

| | resident | streaming, 1 tile | streaming, tiled |
|---|---:|---:|---:|
| D at 600, wall / peak | 2:34 / 47.2 GB | 2:27 / 41.1 GB | 2:57 / 31.8 GB (3 tiles) |
| **TNBC 009 at 600, wall / peak** | 13:12 / **124.1 GB** | 13:49 / **109.9 GB** | 11:33 / **70.9 GB** (13 tiles) |

On the two-hour file the streaming source alone takes 14 GB off the peak (the MS1 map and the MS2
slabs no longer coexist) and tiling takes another 39 GB: **124.1 -> 70.9 GB, -43%, at -12% wall**.
The remaining 70.9 GB is the MS1 pass, which pass 2 never touches.

**Every failure the chain reported was in its own instrumentation, not the tool:** eleven arms
failed a stamp assertion whose grep did not allow for OpenMS's `type="xsd:string"` attribute
between `name` and `value`; the chained slab digest was seeded from a zero-initialised map entry
instead of the FNV offset basis it must share with the whole-window digest; and the e2e suite
could not run at all (the synthetic tdf fixture predated the band-edge default's metadata
requirement, and an early refusal fired on mzML inputs that have no tdf and are supposed to fall
back quietly). All three are fixed (2b1f6e3, 8e2f1c0, 906e758) and re-checked by the `verify`
chain on the rebuilt binary. **A gate whose failures are all in the gate is not a passing gate:
it is a gate that has to be re-run, which is what the verify chain is for.**

**Re-run and CLOSED (chain `verify`, node 1, f07da66, 09:24-09:36): `VERIFY_RESULT PASS`.**
All three settle:

- the stamps read `spx:tile_source = stream`, `band_edges = acquisition`, `frame_table = loader`,
  `tiles = 3`, `tile_rt_sec = 600.0` -- the grep, not the tool, had been wrong (it did not allow
  for OpenMS's `type="xsd:string"` between `name` and `value`, the second time that has cost a
  gate a false failure);
- the chained per-tile slab digests equal the resident whole-window digests over all 24 windows:
  **the six tiled reads deliver exactly the peaks one whole-window read delivers, frame for
  frame**, and the two files are SPECTRUM DATA IDENTICAL;
- e2e **18/18**.

The 3-tile streaming arm on the rebuilt binary reproduces the pin at 3:43 and **32.1 GB**, with
frames delivered 10,416 / 10,416 / 11,378, 0 padded and 0 absent in every tile.

### Tiling has hit its floor on D: the MS1 pass (2026-09-10, chain `instr3`, 4258d16)

Five arms, every one reproducing its pinned digest, on node 2 at 100 threads with the repaired
ledger. Two controlled pairs, one per cell pitch, so grouping is the only variable inside a pair.

| arm | wall | peak RSS |
|---|---:|---:|
| resident, 1 tile | 2:29.83 | 46.77 GB |
| stream, 1 tile, pitch 600 | 2:30.71 | 41.71 GB |
| stream, 3 tiles, pitch 600 | 2:55.30 | **31.76 GB** |
| stream, 1 tile, pitch 300 | 2:36.90 | 38.02 GB |
| stream, 6 tiles, pitch 300 | 2:55.49 | **31.76 GB** |

**The 3-tile and 6-tile arms peak within 1 MB of each other** (31,757,788 vs 31,758,716 kB) and
finish within 0.2 s, at two different pitches and two different digests, because both are pinned by
the same event: the MS1 pass at t≈9 s. Inside the window loop the six-tile arm is much smaller
(slab 3,548 MiB against 5,779 at three tiles) and none of that reaches the process peak. **More
tiles is now free and worth nothing on this file; the next memory gain must come out of the MS1
pass.** Grouping is worth −23.9% peak for +16.3% wall at pitch 600, and −16.5% for +18.5% at 300.

The cell pitch moves the ONE-TILE peak (41.71 -> 38.02 GB, −8.8%, because a smaller cell is a
smaller detection working set) and stops mattering once tiled. So `tile:rt_sec` is a science knob,
not a memory knob, wherever `tile:cells_per_tile` is set.

**The resident and streaming peaks are made of different things.** The resident arm's 45,884 MiB at
t=28.0 s falls to 28,052 two seconds later with the charged total UNCHANGED: 17.3 GiB of it was one
load's worth of allocator free list, not live data. The streaming arm's 30,974 MiB at t=8.8 s has
only 5,495 MiB of glibc free, so **~10.7 GiB there is live and uncharged**. It is mostly
`PeakPickerIM::pickIMCluster`'s own per-spectrum temporaries -- 88.125 B per RAW peak across six
simultaneously live arrays, times up to 100 concurrent picks -- plus the loader's 64-frame MS1 decode
stage (`patch:340`). (A first version of this paragraph blamed a 256-frame loader decode batch; that
batch is the MS2 block and never runs in pass 1, and the correction is in data2.md.) A trim can
only help the resident shape. On the streaming shape the picker flush batch is the lever:
`SPEXTRACTOR_PICK_BATCH=64` takes 3.0 GB off for +2.3% wall, and `SPEXTRACTOR_LOAD_BATCH` is inert on
this peak (data3-batch.md).

### TNBC 009 with the repaired ledger: on the long file the floor is MS1 tracing (2026-09-10, chain `tnbcled`)

| arm | wall | peak | process peak |
|---|---:|---:|---|
| 1 tile | 10:00.62 | 109.12 GB | 106,699 MiB, window loop |
| 13 tiles | 9:12.61 | 70.24 GB | 68,666 MiB at 74.7 s, MS1 partition phase |

Both reproduce `bff54a1f`. **At 13 tiles the process peak is the MS1 partition phase**, 2-4 GiB above
the load peak (64,186 MiB) -- the reverse of dataset D, where the partition phase sits 8.4 GiB below
the load at the default batch. So the staging levers that pay on D (picker batch 64, narrower picker
scratch) are worth ~nothing on TNBC's tiled peak until MS1 tracing comes down. At the peak instant
the RSS is 29,451 MiB of charged partition, **28,605 MiB of glibc free list** (the drained MS1 map)
and 10,610 MiB live but uncharged (MassTraceDetection's internal band copies): a trim inside MS1
tracing is the lead candidate for the long file, with a measured ceiling of 28.6 GiB.

Tiling cut the scorer's scan 12.016e12 -> 1.120e12 fragment visits (−90.7%) with identical survivors
(13,857,222,667), loop CPU −33.7% (−9.7% on D: the untiled scan grows with run length squared). The
13-tile loop's RSS climbs 15.5 -> 51.7 GiB tile over tile while charging ~6 GiB, because nothing
trims between tiles. Details: internal review record.

### The composition on D: MS1-only picker batch + in-loop trim, −15.3% peak (2026-09-10, chain `compose`, 531dad4)

| arm | wall | peak |
|---|---:|---:|
| control (two, bracket 0.07 GB / 3.0 s) | 2:57.05 / 3:00.05 | 31.83 / 31.91 GB |
| `SPEXTRACTOR_PICK_BATCH_MS1=128` | 2:59.03 | 30.55 GB |
| `SPEXTRACTOR_PICK_BATCH_MS1=64` | 3:14.72 | 29.19 GB |
| **`PICK_BATCH_MS1=64` + `SPEXTRACTOR_LOAD_TRIM=16`** | **3:13.09** | **26.99 GB** |
| `PICK_BATCH_MS1=64` + `LOAD_TRIM=4` | 3:32.77 | 26.23 GB |

All six reproduce `8c1b047f`. The MS1-only knob removes the shared knob's tile-read cost (16.6 s
against 23.2 s) and leaves an MS1-pass cost of +1.9 s at 128 and +4.5 s at 64. A trim every 16th
flush breaks the soft window-loop cap the verifiers identified: **26.99 GB, −15.3%, for +14.5 s
(+8.1%)**, all of it attributable to the two knobs. The peak moves back to the end of the load
(13.3 GiB MS1 map + 10.1 GiB free list, ~2.5 GiB live uncharged), which is D's next floor. Trimming
every 4th flush is dominated. **Untested on TNBC 009**, whose tiled floor is the MS1 partition phase,
where neither knob acts. Details: data3-batch.md section 8.

### TNBC 009: the pass-1 levers take 6.3 GB off, smaller tiles take nothing (2026-09-10, chain `tnbc2`, d7cc14d)

| 13 tiles unless stated | wall | peak |
|---|---:|---:|
| control, first / last | 10:16 / 9:24 | 68.94 / 68.02 GB |
| MS1 trim after band distribution | 9:39 | 65.94 GB |
| MS1 trim + MS1-only picker batch 64 + trim every 16th flush | 10:20 | **62.22 GB** |
| pitch 300 (25 tiles; first recorded as 26) | 8:44 | 68.72 GB |
| pitch 300 + the same levers | 9:52 | 62.78 GB |

Pitch-600 arms reproduce `bff54a1f`; the pitch-300 pair agree (`725298dd`). The levers take **−6.3 GB
(−9.1%)** against the mean control, output identical, wall within the controls' 52 s spread; the peak moves
from the MS1 partition phase to the end of the load, where the MS1 map (44.2 GiB at the instant) binds.
**Smaller tiles do not lower the peak** at either lever setting, because the peak is in pass 1. Pitch 300
halves the scorer's fragment visits (−47%) and trims loop CPU ~10%, but it changes the output and is a
science knob. Details: internal review record.

### dnoise v0.1.0 as a preprocessor: every MS1/MS2 combination, both files (2026-09-10, chains `dnoise_d`, `dnoise2_d`, `dnoise_tnbc`, `dnoise2_tnbc`, d7cc14d, `-tile:cells_per_tile 1`)

| configuration | D Sage | D MSFragger | D peak | TNBC Sage | TNBC MSFragger | TNBC peak |
|---|---:|---:|---:|---:|---:|---:|
| control | 13,491 | 12,918 | 31.84 GB | 29,067 | 27,949 | 68.82 GB (node 2), 68.41 (node 3) |
| no-op, every stage off | identical | identical | 32.17 (12.91) | identical | identical | 69.73 (27.53) |
| MS1 only, dnoise's default | −4.3% | −4.0% | 21.41 (12.90) | +1.5% | −0.5% | 25.08 (25.24) |
| MS2 only | −5.4% | −13.0% | 32.13 (10.51) | −10.0% | −20.9% | 67.80 (26.97) |
| MS1 + MS2 | −7.2% | −15.3% | 8.13 (10.29) | −7.8% | −15.3% | 19.87 (26.32) |

Peak is DIAspeXtract's tiled peak, with dnoise's own run in brackets; the pipeline peak is the larger of the
two. TNBC's MS1-only and MS1 + MS2 arms ran on node 3, the rest on node 2, each against a control on its own
node. The no-op reproduces both pins, so dnoise's rewrite changes nothing DIAspeXtract reads. Entrapment
never rises above the control's 95% interval: D's arms stay inside it (at most +0.16 points), and TNBC's
all fall just below it (0.92-0.93% against 1.09%).

**MS1 filtering unlocks the memory.** Pass 1 sets the tiled peak, so MS2 filtering alone changes nothing.
With MS1 filtered, dnoise's own run sets the pipeline peak on TNBC; on D it does so once both are filtered.
**MS2 filtering buys the time and costs MSFragger** 13-21% of its peptides on both files. On diaPASEF,
dnoise's MS2 mode is its MS/MS streak filter inside each window. Sage's lost peptides are the faint,
low-scoring ones: 23-56% of the control's faintest intensity quintile, at most 2.5% of the brightest. Only
TNBC's MS1-only arm gains net peptides, and only on Sage. None of this was measured at the shipped single
tile, where the window loop sets the peak.

As a preprocessor, MS1 only fails the S gate on D. The user decided to run it inside the tool by default anyway
(next section). Details: internal review record.

### dnoise MS1 inside the tool, on by default (2026-09-10, stage-2 build fd30bf0e, arm `dnoise_s2a`)

`src/DnoiseMs1.h` ports dnoise v0.1.0's default MS1 path exactly. The tool recovers each MS1 point's flight-time
index, scan index and raw intensity from the loader's values and filters before picking (`dnoise:ms1`, default
true). The port reproduces dnoise's keep masks with 0 differing points on 28 D and 24 TNBC 009 frames. The
whole-file gates reproduce DIAspeXtract run on dnoise's own filtered file:

| file | arm (one cell per tile, 100 threads) | digest | kept MS1 points | peak | wall |
|---|---|---|---:|---:|---:|
| D | `-dnoise:ms1 false` | `8c1b047f` (old pin) | all | 26.46 GB | 4:23 |
| D | defaults | `e43672a0` (= dnoise's file) | 118,053,254 | 22.41 GB | 3:44 |
| D | defaults, 8 threads | `e43672a0` | 118,053,254 | 13.29 GB | 12:07 |
| TNBC 009 | `-dnoise:ms1 false` | `bff54a1f` (old pin) | all | 34.67 GB | 9:47 |
| TNBC 009 | defaults | `3aafd0b5` (= dnoise's file) | 649,615,267 | 25.41 GB | 8:18 |

- **New reference digests at the defaults:** D `e43672a0`, TNBC 009 `3aafd0b5`.
- **Searches carry over** from dnoise's filtered files:
  - D: Sage 12,917 (−4.3%), MSFragger 12,404 (−4.0%);
  - TNBC 009: Sage 29,510 (+1.5%), MSFragger 27,808 (−0.5%).
- **Cost and speed:** the filter costs 355 CPU-s on D and 1,121 on TNBC 009 at 100 threads. The MS1 load still gets
  faster (D 30.9 → 16.5 s wall), because the picker sees an eighth of the points.
- **Final-build numbers** (arm `dnoise_s2b`, one cell per tile, every arm at its pin;
  internal review record):
  - Defaults: D 22.48 GB at 3:41, TNBC 009 25.17 GB at 10:46.
  - With `-dnoise:ms1 false`: D 25.90 GB, TNBC 009 31.97 GB.
  - TNBC 009's MS1 load gets faster on node 3 (55.7 → 37.0 s) but not on node 2 (58.5 → 68.4 s), so its wall
    change depends on the node.

### MS1 prune, output-identical (2026-09-10)

`perf:ms1_prune` (default true) drops picked MS1 centroids at or below `trace:noise_threshold_int` at load. A
witness chain and a pre-prune edge sample keep band membership and the tracer's input exact
(internal review record).
- **Identity:** digests are identical with the prune off and on for D (`8c1b047f`, with all 245 `[det]`/`[ms1-edges]`
  lines identical), for dnoise-filtered D (`e43672a0`) and for TNBC 009 (`bff54a1f`).
- **Peak** at one cell per tile, dnoise off (final build, back to back): D 31.99 → 25.90 GB (−19%), TNBC 009
  69.11 → 31.97 GB (−54%). With dnoise on, the prune takes only 0.14 GB more off D and 1.0 GB more off TNBC 009.

### Why more tiles still leave a large footprint (2026-09-10, instr_tnbc2)

Pass 1 processes the whole run's MS1 and is never tiled. The tile loop's RSS follows the glibc arena that pass 1
leaves, because freed chunks are re-dirtied tile after tile and nothing trims between tiles.

Peak by tile count, at identical output:

| input | peaks |
|---|---|
| TNBC 009, original | 107.2 / 78.1 / 69.1 / 68.8 / 69.5 GB at 1 / 2 / 4 / 7 / 13 tiles |
| TNBC 009, dnoise-filtered | 97.8 / 48.7 / 24.9 / 20.4 GB at 1 / 4 / 13 / 25 tiles |
| D, dnoise-filtered | 36.9 / 29.3 / 21.4 / 14.5 GB at 1 / 2 / 3 / 6 tiles |

- **Trimming off:** with `perf:malloc_trim false` the loop runs flat at 60-64 GB instead of climbing from 15 to 49 GB.
- **Eager trim threshold** (`MALLOC_TRIM_THRESHOLD_=67108864 MALLOC_TOP_PAD_=0`, zero code): it cuts 14-22% of the
  peak at identical output for 3-6% wall (TNBC 68.7 → 59.0 GB, filtered TNBC 25.2 → 21.4, filtered D 14.5 → 11.3).
- **`malloc_trim` between tiles:** it takes only 0.8 GB off filtered D, where the growth happens inside a tile.
- **On the final build, at the defaults** (arm `dnoise_s2b`, identical output in every arm):
  - TNBC 009 on `data`: 99.70 GB and 9:42 at one tile, against 21.07 GB and 9:57 at 13 tiles with a trim between
    tiles (20.73 GB with the eager threshold, which nearly triples the MS1 load);
  - D on node 1: 37.44 GB and 3:15 at one tile, against 22.48 GB and 3:41 at 3 tiles, and 21.99 GB and 3:47 at
    3 tiles with a trim between tiles (18.14 GB and 4:06 with the eager threshold).
  - Decided 2026-09-11: `tile:cells_per_tile` 1 and the trim between tiles are the defaults (TNBC 009 21.07 GB / 9:57,
    D 21.99 GB / 3:47); the eager threshold is not adopted (CHANGELOG).

Details: internal review record.
