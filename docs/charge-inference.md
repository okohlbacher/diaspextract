# READ THIS FIRST (2026-09-05): the confusion table below points the OTHER WAY

> **Note (2026-09-11, release 1.2.0).** Two options this document argues about --
> `assembly:open_search_safe` and `charge:ambiguity_margin` -- were **removed in 1.2.0**; the shipped
> binary refuses either on the command line or in an ini file. The sections below are kept as the
> falsification record that removed them, and the help strings and comments they quote are gone with
> the options. Everything said about them is history, not the current tool.

Re-verified against the code 2026-09-05, every claim below checked by hand:
**most of this document argues against the wrong failure.**

**1. `4->2` means WE said 4 and the truth was 2 -- an OVERCALL, not an undercall.**
The benchmark collation keys the counter `zc[(z, best["Precursor.Charge"])]` with `z` read off OUR
mzML, and prints it `f"{a}->{t}"`. So the key is (ours, truth). The delta-mass a search
engine sees is ours minus truth = **+1198 Da** at m/z 600, not the -1198 that was quoted in the
`assembly:open_search_safe` help string (removed in 1.2.0) and repeated below. The repository does not agree with
itself on this.

**2. That table belongs to a retired algorithm whose tie-break is the OPPOSITE of the shipped one.**
The table is the `envelope` arm. Envelope breaks score ties toward the HIGHEST z
(`if (a.z != b.z) return a.z > b.z;` in the envelope comparator). The shipping `count` arm breaks
ties toward the LOWEST z (strict `>` over an ascending z loop from `best_n = 1`).
**The shipped arm's confusion cells have never been measured -- they are `--` in the table below.**
So the structural bias of the code that actually ships is toward UNDERCALL, and the fixes discussed
in this file are aimed at an overcall produced by a tie-break the shipped code does not contain.

**3. The envelope 49.7% measures a bug, so "averagine cosine cannot break z vs 2z" is unproven.**
`xicCorr` returns a sentinel `-2.0` when fewer than three frames overlap or a variance is
zero. That sentinel is then averaged in and clamped:
`mean_corr = cn ? std::max(0.0, csum / cn) : 0.0`. One isotope with <3 shared frames therefore
zeroes `cos_sim * mean_corr * evidence` outright, every poisoned hypothesis ties at exactly 0.0, and
the tie-break hands the seed the HIGHEST charge that found two partners. That is a complete
mechanical account of the 15,734 `4->2` cells. It fires often: `trace:min_length_sec` is 3.0 s
against a measured MS1 FWHM of 3.61 s, so MS1 traces routinely sit on the `m < 3` cliff. The
reference tool reaches 82.6% using exactly the criteria this file concluded were insufficient.

**4. The one unambiguous code defect: the count walk always starts at z=1.**
The walk is `for (int z = 1; z <= max_charge; ++z)` regardless of `charge:min_charge`. A seed won by
z=1 marks its seed and every partner `used[]`, and only afterwards are z<zmin precursors erased.
**The peaks are spent on a precursor that is then deleted, and no z>=2 hypothesis is ever formed for
that seed.** This is an ownership defect and it is one line. It is dormant at the shipped default
(`charge:min_charge=1` emits the z=1 call rather than erasing it) and live whenever `min_charge` is
raised above 1.

**5. `assembly:open_search_safe` was a strict no-op, and the 2026-09-04 default made it worse.**
It set `pc.charge = 0` only where `pc.charge == 0` already, and marked exactly that set
`guessed`; the emission path then erased every `guessed` precursor whenever
`assembly:require_isotope_support == "true"` -- the default since 2026-09-04. The precursors the
flag existed to protect were deleted before it could act, and its comment still read "OFF by
default". **This is what removed it in 1.2.0.**

**6. Charge errors cost CLOSED-search identifications today, in the published numbers.**
`msfragger.params` has `override_charge = 0`, so MSFragger TRUSTS our annotation. Measured
2026-09-02 (the charge-corruption test in `docs/BENCHMARK-MATRIX-2026-09-01.md`): forcing the engine to re-derive charge
instead (`override_charge=1, z=1..5`) moved DIAspeXtract 10,750 -> **11,238 (+4.5%)** against the
reference tool's 13,004 -> 13,276 (+2.1%), narrowing the gap 2,254 -> 2,038. So roughly 10% of the
MSFragger deficit is our charge labels, and it is recoverable. The counter-claim that stood in the
`assembly:open_search_safe` help string -- "Closed search barely notices (charge unset cost only
1.3%: 8,123->8,019)" -- had no provenance anywhere in the repository.

**7. Every published charge number describes a tool that no longer exists.** 74.6% and the arm table
predate `charge:min_charge=2`, `assembly:require_isotope_support=true`, the apex m/z estimator, the
integer detector and `perf:stream_load` -- the last of which is not output-neutral. The measured
population is disjoint from the shipped one, so 74.6% cannot be differenced against anything. It is
nonetheless still quoted as a live adoption gate in `bench/README.md`.

**8. Charge is the only metric in the harness with no chance floor.** The shifted-query control in
the benchmark collation records `matched_decoy` but not `charge_pairs`. Correcting this can only
help: on a chance match our z is independent of truth, so chance agreement is
`sum_z p_ours(z) p_truth(z) <= max_z p_truth(z)` -- chance matches drag the score toward or below
the constant baseline, so the honest margin over "always answer z=2" is **wider** than +5.0 points.

**THE ADJUDICATION EXPERIMENT** that would settle the direction question empirically: pair our
spectra with the reference tool's on co-identified precursors and decide from the spectra which is
right. Note its central trap -- both engines TRUST the declared charge, so the adjudicating search
must re-derive it or the test merely echoes our own annotation back at us.

**WHAT TO DO, in order.** (a) Fix the z=1 ownership defect -- start the walk at `min_charge` -- and
gate it on both engines plus entrapment, since it changes what is emitted. (b) Re-measure the
confusion table on the SHIPPED arm before designing anything else; the direction of the error is
currently unknown. (c) Fix the `xicCorr` sentinel, then re-test envelope scoring, which has never
been fairly evaluated. (d) Decide whether to ship `override_charge=1` in the benchmark config, or
treat the +4.5% as the size of the prize for fixing charge properly. Do NOT build a learned charge
model against the 74.6% gate: the number is stale and the gate's 0.17% noise floor is itself retired.

---

# Charge inference: what the 15.9% is, and what will not fix it

The MS1 funnel analysis localises **15.9%** of precursor loss to wrong charge assignment and
**10.5%** to failed monoisotope identification — against only 2.2% undetected. Truth-set
charge agreement is **47.8% vs the reference implementation's 82.6%**, dominant confusion `4→2` (15,734 cases).

## Why 4→2 happens

Isotope spacing is 1/z. At z=4 that is 0.25 Da; at z=2, 0.50 Da. A spurious peak halfway
between two real isotopes is indistinguishable from a genuine z=4 partner **on spacing
alone** — and the higher-charge hypothesis always has more chances to find one, so the
error is systematically biased toward overcalling charge.

Cost in open search: ΔM = (z′−z)(m/z − 1.00728). Measured at m/z 600:

| ours | true | phantom ΔM |
|---|---|---|
| z=4 | z=2 | **−1198.0 Da** |
| z=3 | z=2 | −599.0 Da |
| z=1 | z=2 | +599.0 Da |

## RETRACTED FALSIFICATION: ion mobility was retired on a wrong null

**This section previously read "FALSIFIED: ion mobility cannot break the tie". That conclusion
rested on a statistical error and is withdrawn.**

The nearest-neighbour test asks whether a precursor's nearest neighbour in (m/z, IM) **shares**
its charge. That is a coincidence rate between two draws, so the null is `sum(p_z^2)`. The
script compared it against `max(p_z)` instead — the null for *predicting* one label:

| | value |
|---|---|
| null used (`max p_z`) | 69.6% |
| **correct null (`sum p_z^2`)** | **55.8%** |
| observed | 75.0% |
| margin as published | +5.4 points |
| **margin corrected** | **+19.2 points** |

So the (m/z, IM) plane carries **substantially more** charge information than reported, and
"do not build an IM-based charge prior" was wrong.

What still stands, because it is a direct observation rather than a null comparison:

Hypothesis was that diaPASEF's mobility dimension constrains charge — the acquisition
tiles themselves track the charge-2 line, so charge ought to be separable in (m/z, 1/K0).
Tested against DIA-NN dataset A ground truth (43,499 precursors with m/z, IM **and** true charge):

| m/z 600–700 | z=2 | z=3 | z=4 |
|---|---|---|---|
| 1/K0 | 0.993 ± 0.031 | 0.937 ± 0.044 | **0.966 ± 0.035** |

**z=4 lies between z=2 and z=3 — mobility is not monotonic in charge.** At fixed m/z a
higher-charge peptide is proportionally heavier (mass = z × m/z), so its larger CCS
cancels the extra charge, non-monotonically.

Nearest-neighbour test (m/z ±2 Da, closest IM): **75.0% share charge vs a 69.6% chance
baseline — +5.4 points.** The plane carries almost no charge information beyond the base
rate. Do not build an IM-based charge prior.

This also retires the broader "use `ScanNumBegin`/`ScanNumEnd` to gate charge" idea. The
tile IM range is still worth using as a *transmission* filter (a precursor outside the
tile's mobility band could not have produced its fragments), but not as charge evidence.

## What is left

Isotope-envelope **intensity** matching: for a genuine z=4 precursor of mass ~2400 Da the
averagine model predicts specific M/M+1/M+2 abundance ratios. A spurious +0.25 Da peak
will not match them. Spacing is degenerate between z and 2z; **abundance is not**.

## CORRECTION (2026-07-21): envelope scoring is default-ON, not default-OFF

An earlier version of this file claimed `charge:scoring=envelope` was an untested lever
waiting to be enabled. It **was** the default until 2026-07-21; the code now ships
`registerStringOption_("charge:scoring", "<mode>", "count", ...)`. A benchmark arm passing
`-charge:scoring envelope`
returned byte-identical output to baseline (6509 spectrum_q / 5813 peptide_q), which is
what exposed it.

This matters for how every number here is read:

* The **15.9% wrong-charge rate and 47.8% charge agreement are what envelope scoring
  already produces.** They are not a "before" measurement.
* Averagine-shape cosine x isotope-XIC co-elution, with gapped envelopes allowed, is
  therefore **not sufficient** to break the z vs 2z degeneracy. The abundance argument
  above is correct in principle but the current implementation of it is not working.
* The informative experiment is the opposite direction: `-charge:scoring count` (the
  legacy partner-count heuristic, ties favouring LOW charge, break-on-first-miss) to
  establish whether envelope is better or worse than the alternative. Until that runs,
  we do not know that envelope scoring is helping at all.

The dominant confusion being `4->2` (overcalling charge) is consistent with a scorer that
rewards finding MORE envelope members: a z=4 hypothesis has twice as many candidate slots
as z=2 in the same m/z span, so gapped-envelope tolerance may actively favour it. Worth
checking whether the evidence weighting penalises hypothesis complexity enough.


---

# Moved here from the README (2026-09-04)

The README became a tool description; these measurements belong with the topic. Verbatim,
except that the comparison tool's name is written as `the reference tool`.

## Charge inference: the two defaults are coupled, not independent

Measured against a DIA-NN dataset B reference (`charge agreement` = our charge vs truth on matched
spectra; always answering z=2 would score **69.6%**):

| arm | charge agreement | `4->2` errors |
|---|---|---|
| ms1split + ambiguity_margin | **41.7%** | 34,554 |
| ms1split alone | 45.8% | 19,995 |
| baseline (`envelope`) | 49.7% | 10,199 |
| `count` | 71.8% | -- |
| **`ms1split` + `count`** (default) | **74.6%** | -- |

**MS1 splitting must not ship without `count`** — an *empirical* coupling, not an architectural
one. The two touch separate code paths (`trace:ms1_split_valleys` in `detectTraces_`,
`charge:scoring` in `inferPrecursors_`); nothing in the code enforces the pairing, so the
combination is a benchmark finding that a future edit could silently break. On its own MS1
splitting drops charge agreement below baseline and
more than triples the `4->2` confusion, which fabricates a ~-1198 Da phantom modification at
m/z 600 -- invisible to closed search (the engine enumerates charge) and fatal to open search.
Shipping `trace:ms1_split_valleys` without `charge:scoring count` would have been actively
harmful for this tool's actual purpose.

### `charge:ambiguity_margin` -- falsified, four tests (option removed in 1.2.0)

| margin | peptide_q | vs `split_count` 8,411 |
|---|---|---|
| 0.5 | 8,408 | -0.04% |
| 1.0 | 8,406 | -0.06% |
| 2.0 | 8,397 | -0.2% |

2.0 is the first value that can genuinely fire against an integer partner count (it admits any
hypothesis within 2 partners of the best) and it still changes nothing. Hedging charge does not
help; deciding it well does. Left at 0.
