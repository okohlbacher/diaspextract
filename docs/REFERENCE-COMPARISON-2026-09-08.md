# DIAspeXtractor vs the reference implementation, measured 2026-09-07/08

The first three sections (runtime, quality, the MSFragger deficit) were measured on the **previous
default** (`charge:min_charge=2`, `perf:ms1_trace_bands=48`); the Decision section on the two builds
that introduced the mobility veto, the second of which is the **shipped v2.0.0 default**
(`charge:min_charge=1`, `charge:im_charge_veto` on). Every figure is against the reference implementation's own
output, using one scoring definition applied to BOTH sides:
Sage = unique peptides at rank 1 with `peptide_q <= 0.01`; MSFragger = the peptide-level target-decoy
walk. Entrapment FDR uses the human-target + Arabidopsis-entrapment search with the peptide-hypothesis
ratio 0.6805.

## Runtime and memory (TNBC 009, one 224-core node, 100 threads both tools, same hour)

| | rep 1 | rep 2 | mean |
|---|---|---|---|
| reference implementation | 21:12.24 / 357.7 GB | 20:41.46 / 323.4 GB | 20:56.9 / 340.6 GB |
| **DIAspeXtractor, previous default** (min_charge 2) | 11:24.97 / 179.1 GB | 10:36.92 / 178.5 GB | **11:00.9 / 178.8 GB** |
| ratio, previous default | | | **0.53x wall, 0.53x RSS** |
| **DIAspeXtractor, shipped default** (one run) | 12:42.21 / 187.9 GB | | 12:42.21 / 187.9 GB |
| ratio, shipped default | | | **0.61x wall, 0.55x RSS** |

Re-measured rather than taken from the earlier campaign: a cross-time comparison is what produced
one withdrawn claim. Before this week's work the same file ran **31:26.98** against the reference's
21:28.00, i.e. 1.46x SLOWER; at the previous default it was 1.90x faster, and at the shipped default
(z=1 emitted, +67% spectra) it is 1.65x faster (20:56.9 / 12:42.21).

Caveat on the memory ratio: the reference is a JVM run with `-Xmx400G`, so its peak RSS partly
reflects what it was allowed rather than what it needs (the same flags gave 268.2 GB in the earlier
campaign). DIAspeXtractor's 179 GB (previous default) / 188 GB (shipped default) is a genuine
requirement and is not sensitive
to a heap flag.

## Quality, two files

The previous default (min_charge 2):

| | TNBC 009 | | | TNBC 001 | | |
|---|---:|---:|---|---:|---:|---|
| | SpeX | ref | ratio | SpeX | ref | ratio |
| Sage peptides | 24,213 | 18,920 | **1.280** | 15,694 | 14,764 | **1.063** |
| Sage proteins | 5,813 | 4,944 | **1.176** | 4,074 | 3,802 | **1.072** |
| MSFragger peptides | 24,199 | 27,427 | **0.882** | 15,976 | 18,745 | **0.852** |
| entrapment FDR (SpeX) | 1.01% (0.86-1.20) | | | 0.98% (0.80-1.19) | | |

The shipped default (min_charge 1 + mobility veto), each engine's own pooled cut:

| | TNBC 009 | | | TNBC 001 | | |
|---|---:|---:|---|---:|---:|---|
| | SpeX | ref | ratio | SpeX | ref | ratio |
| Sage peptides | 27,855 | 18,920 | **1.472** | 19,035 | 14,764 | **1.289** |
| Sage proteins | not on record | 4,944 | | 4,142 | 3,802 | **1.089** |
| MSFragger peptides | 26,559 | 27,427 | **0.968** | 18,212 | 18,745 | **0.972** |
| entrapment FDR (SpeX) | 1.03% | | | 0.93% | | |

No shipped-default Sage protein count for 009 is on record. The previous-default 009 protein figure
of 5,813 above appears only in this file -- it is in neither docs/BASELINE.md nor docs/TNBC-TABLE.md
-- so treat it as unverified.

The direction is consistent across both files and the magnitude is not: DIAspeXtractor finds more with
Sage and fewer with MSFragger, and at the previous default the Sage advantage was 1.28 on one file
and 1.06 on the other. At the shipped default the pooled Sage ratios are 1.47 (009) and 1.29 (001),
per-charge 1.28 and 1.25 (Decision
section below). Entrapment FDR is controlled at the nominal 1% in both, so the Sage advantage is not
FDR inflation.

## Where the MSFragger deficit comes from

Characterising the 6,749 MSFragger peptides the reference finds on 009 and we do not (previous
default, min_charge 2):

| set | n | median hyperscore | z=1 | z=2 | z>=3 |
|---|---:|---:|---:|---:|---:|
| only reference | 6,749 | 18.19 | **21.1%** | 55.4% | 23.5% |
| common | 21,641 | 26.30 | ~0% | 70.0% | 29.9% |
| only DIAspeXtractor | 3,428 | 18.56 | ~0% | 70.1% | 27.8% |

Two things. First, both exclusive sets are **marginal** -- median hyperscore ~18 against ~26 for the
peptides both tools find -- so a large part of the "gap" is churn at the acceptance boundary rather
than signal we never emitted. Second, **21.1% of what we miss is charge 1**, against essentially none
of what we share: a class `charge:min_charge=2` never emits. The 21.1% z=1 component of this set is
what the shipped default now emits; the residual MSFragger deficit at the shipped default (0.968
pooled, 1.00 per-charge) has not been re-characterised.

Confirmed by measurement, on 009 (both arms on the previous-default binary; the z>=1 arm is
min_charge 1 without the veto, so its figures are not the shipped ones, which are in the Decision
section):

| | z>=2 (the previous default) | z>=1 (same binary, min_charge 1, no veto) |
|---|---:|---:|
| Sage peptides | 24,245 | **27,685 (+14.2%)** |
| MSFragger peptides | 24,121 | **26,332 (+9.2%)** |
| Sage proteins | 5,802 | 5,641 (-2.8%) |
| entrapment FDR | 1.03% | **1.01%** |
| spectra | 2,343,822 | 4,035,058 (+72%) |
| wall | 11:59 | ~13:20 (+11%) |

Repeated on a second file, and it is stronger there:

| | TNBC 009 z>=2 -> z>=1 (previous default, no veto) | TNBC 001 z>=2 -> z>=1 (previous default, no veto) |
|---|---|---|
| Sage peptides | 24,245 -> 27,685 (**+14.2%**) | 15,694 -> 18,854 (**+20.1%**) |
| Sage ratio vs ref | 1.280 -> **1.463** | 1.063 -> **1.277** |
| MSFragger peptides | 24,121 -> 26,332 (**+9.2%**) | 15,976 -> 17,068 (**+6.8%**) |
| MSFragger ratio vs ref | 0.882 -> **0.960** | 0.852 -> **0.911** |
| Sage proteins | 5,802 -> 5,641 (-2.8%) | 4,074 -> 4,118 (+1.1%) |
| entrapment FDR | 1.03% -> **1.01%** | 0.98% -> **0.90%** |
| spectra | +72% | +64% |
| wall | +11% | +20% |

So on this dataset class the z>=2 default costs 14-20% of Sage peptides and most of the MSFragger
deficit, with entrapment FDR unchanged or better on both files. That directly contradicts the IH2
measurement that established the default (where dropping z=1 GAINED peptides on both engines), so
`charge:min_charge` is dataset-dependent. Resolved below: the default changed to 1. Cost at the
shipped default on 009: +67% spectra (2,343,822 -> 3,918,321), +15% wall (12:42 against the 11:01
previous-default mean), no memory beyond run-to-run spread. The 009 protein-count drop (-2.8% in the
previous-default arm above) was NOT re-measured at the shipped default -- no 009 protein figure is on
record for the shipped build; on 001 the protein
count rose 4,074 -> 4,142.


## What the charge-1 spectra are (resolved 2026-09-08)

The +14-20% Sage gain from emitting z=1 raised a question: had we not already fixed a
charge-assignment defect, and was this the same class of error? Direct analysis of the extracted
spectra settled it. There are **two separate things**
inside the z=1 stratum, and the whole-run entrapment figure hid one of them.

### 1. The identifications are genuine singly charged peptides

| evidence | result |
|---|---|
| ion mobility vs the 2+ trend line (1/K0 = 0.5026 + 0.000729 m/z, fit on 15,001 z=2 IDs) | z=1 IDs sit **+0.29 Vs/cm2 above it** in every 100-Th bin; 99.4% above, 0.6% on the 2+ band. A halved 2+ ion would carry the 2+ mobility. |
| the reference implementation's own 1+ IDs | identical offset (+0.283 to +0.296), 99.9% above the 2+ line -- the same population |
| sequence chemistry | 92% are 7-9 residues, ~920 Da, 86% end in K, only ~12% carry any internal K/R/H (69% for z=2): single-basic-site tryptic peptides, the class that flies as 1+ |
| same sequence seen at both charges | 570 (009) / 422 (001) peptides: co-elute at dRT = 0, with the 1+ ion 0.65 1/K0 above its own 2+ in 100% of pairs |
| match quality | MSFragger matches 75% of theoretical b/y for z=1 (55% for z=2, 33% for decoys); extracted DTAs show mono m/z within ~5 ppm and b/y ladders at 1+ |
| the discriminating experiment | same spectra emitted with the charge UNSET so the engines search 2-4: **only 319 of 5,947 (5.4%) come back**, all at 2+ -- the dual-charge minority -- and the arm falls below the previous default (23,938 vs 24,245) |
| why the engines report 1+ at all | `precursor_charge 2-4` / `override_charge 0` is the fallback when the file carries no charge; both tools write charge 1 on these spectra and both engines honour it |

Neither tool has 1+ identifications below m/z 700, because a 1+ there is under the engines' mass
floor (500 Da, 7 residues).

### 2. Charge halving is real -- in the calls, not the identifications

The count walk calls a 2+ envelope z=1 whenever its M+1 (+0.5017) is missing or already claimed and
its M+2 (+1.0034) is present: z=1 then wins outright (M0, M2, M4 = 3 partners against z=2's 1).
Measured on 009: **~13% of the 2.55 M z=1 calls sit on the 2+ mobility band** (p2 = 51% there against
a 2.3% chance rate on the 1+ band), 60% of them below m/z 750 where a 1+ cannot be scored at all --
dead emission, and a lost 2+. Dropping z=1 (the old default) did not undo it: the walk had already
consumed the 2+ envelope's peaks, so the old default lost those 2+ ions too. **Fix:**
`charge:im_charge_veto` fits the 2+/3+ mobility lines per run from the run's own confident calls
(>= 3 isotopes, 25-Th bin medians, MAD residual scale) and, at the shipped default, re-calls a z=1
whose mobility lies within 2.5 sigma of the 2+ line to 2+ and drops a z=1 on the 3+ line (its charge
is unset and `require_isotope_support` removes it; the first veto build re-called both, see the
Decision section).
Genuine 1+ ions are 7-8 sigma away and untouched.

### 3. The z=1 stratum runs ~2.8% entrapment FDR at a nominal 1% cut -- and the whole-run figure hid it

| stratum | 009, z>=1 arm (previous default, no veto) | 001, z>=1 arm (previous default, no veto) |
|---|---|---|
| z=1 | 5,145 peptides, **2.80%** (CI 2.25-3.30) | 3,465, **2.84%** (2.24-3.49) |
| z=2 | 0.43% | 0.39% |
| z=3 | 0.91% | 0.52% |
| whole run | 1.01% | 0.90% |

At the shipped default the record is: 009 z=1 2.94% pooled / 0.91% per-charge; 001 z=1 2.91% (n=3,481)
pooled / 0.87% per-charge, z=2 0.41%, z=3 0.57%, whole run 0.93%.

This is not an extraction artefact: decoys reproduce it (3.3-4.0% of z=1 targets clear the global
threshold against 0.3% of z=2; z=1 is in fact the best-calibrated stratum, entrapment/decoy ratio
1.3 against 2.2 for z=2), zero z=1 peptides sit on the 2+ band, and the false z=1 hits are themselves
on the 1+ mobility line -- genuine 1+ ions matched to the wrong sequence, not halved 2+ ions.

It is also NOT "short peptides are intrinsically hard": at matched discriminant score *and* matched
length, z=1 8-10-mers are 3-5x dirtier than z=2 8-10-mers (009 score octiles: 11.0/2.9/4.1/2.8%
against 2.2/1.7/1.3/0.0%), and within z=1 at equal score 7-mers are no worse than 8-10-mers. The
mechanism is a **lower true-match prior in the z=1 stratum under one pooled threshold**: 63% of
emitted z=1 spectra produce no rank-1 row at all (32% for z>=2) and 7.5% of the rest reach
spectrum q <= 0.01 (15.9% for z>=2). Any subgroup with a lower prior sits above 1% under a global
cut. Adding the z=1 search space also refits Sage's LDA (discriminant scores change for 100% of the
shared z>=2 PSMs) and *lowered* the z>=2 strata (0.85 -> 0.43%), which is exactly why the pooled
1.01% looked clean. The reference implementation's own 1+ stratum shows the same thing (decoy FDR
1.4% on 009, 7.6% on 001), so it is an engine-side calibration property of the class.

**Remedy, measured (shipped default on 001):** a per-charge FDR walk puts z=1 at 0.87%, z=2 at 1.14%, z=3 at
0.68%, whole run 0.99%, while keeping the gain -- 19,262 Sage peptides per-charge against 16,029 for
the previous default under the same walk (**+20.2%**). The SpeX-only z=1 peptides carry protein
co-support (another accepted z>=2 peptide on the same protein) in 92.9% of cases, against 39% for
accepted decoys.
**A whole-run entrapment figure must never again be quoted without its charge strata.**

One caveat on the 001 MSFragger numbers (previous-default arms): MSFragger's parameter optimisation chose
`intensity_transform=1` for the default arm and `0` for the z>=1 arm, so its +6.8% on 001 mixes the
z=1 effect with an engine-parameter flip triggered by the changed spectrum population; 009 is clean.
Whether the same flip occurred between the previous-default and shipped-default 001 arms
(15,976 -> 18,212) is not recorded.

### Decision

`charge:min_charge` default 2 -> 1, `charge:im_charge_veto` on. The z=1 peptides are real and the
reference implementation reports them too; the halved 2+ calls are corrected by physics the
instrument already measured; the residual stratum FDR is a search-engine calibration matter and is
now documented rather than hidden.

One refinement after the first measurement: the veto's 3+ arm was re-calling 7% of z=1
calls to 3+, but only ~1.8% of z=1 calls are genuinely on the 3+ band (the 2+ and 3+ mobility lines
are ~1 sigma apart at low m/z), and the re-calls identified nothing on either file (z=3 peptides
5,313 -> 5,322 on 009, 3,888 -> 3,935 on IH2) while adding 178k / 93k re-called precursors
(117k / 73k spectra). A z=1 call on
the 3+ band is now dropped rather than re-labelled. The IH2 pair is on record in
docs/BASELINE.md; the 009 pair exists only in the source comment beside the veto's 3+ branch in
`src/spextractor.cpp` and in no docs/BASELINE.md entry.

The like-for-like comparison, both tools under the same per-charge FDR walk, on 009 (the shipped
build):

| | DIAspeXtractor | reference | ratio |
|---|---:|---:|---|
| Sage peptides, per-charge 1% | 27,831 | 21,743 | **1.28** |
| MSFragger peptides, per-charge 1% | 28,181 | 28,294 | **1.00** |
| Sage, each engine's own pooled cut | 27,855 | 18,920 | 1.47 |
| MSFragger, pooled cut | 26,559 | 27,427 | 0.97 |

Per-charge control moves the reference more than us on Sage (its z>=2 strata sat below 1% under its
pooled cut), which is why the honest Sage ratio is 1.28 rather than 1.47. On MSFragger the two tools
are now at parity.

The same on the second file, at the shipped default (the 3+ arm dropped), against the reference on 001:

| | DIAspeXtractor | reference | ratio |
|---|---:|---:|---|
| Sage peptides, per-charge 1% | 19,262 | 15,418 | **1.25** |
| MSFragger peptides, per-charge 1% | 19,230 | 19,652 | **0.98** |
| Sage, pooled cut | 19,035 | 14,764 | 1.29 |
| MSFragger, pooled cut | 18,212 | 18,745 | 0.97 |

Two files, two engines, two walks: +25-28% on Sage and within 2.2% of parity on MSFragger, with the
per-charge z=1 stratum at 0.87-0.91% and whole-run entrapment at 0.93-1.03%.
