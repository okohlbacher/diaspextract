# The emission excess is duplicate claims on one monoisotope, and folding them together is nearly free

**Question (2026-09-20).** DIAspeXtract emits about 1.45x the spectra of the reference implementation for about the same
peptides. How much of that emission is the same precursor emitted twice, which spectra can go without losing
identifications, and what does removing them cost?

Line numbers below refer to `src/diaspextract.cpp` at v1.3.2, the source this study examined (`git show
v1.3.2:src/diaspextract.cpp`); 1.4.0 added the fold described in section 7.

**Verdict.** The excess is real, it is far larger than the isotope duplication the project has been chasing, and it
is removable. **45.9% of inferred precursors have a competitor claiming the same monoisotope.** Deleting the
duplicates costs 1.3-1.9% of peptides; *merging* them -- keeping one precursor and unioning the fragments --
removes 23.8% of the emission for **0.21%** of Sage peptides and 0.18% of MSFragger peptides, or 31.6% for 0.53%
and 0.40%. At the wider setting the file drops to 719,631 spectra, **below the reference implementation's 725,731
on the same run**, while its 50,533 Sage peptides are 98.4% of the reference implementation's 51,347 there (both
extractors searched with the same configuration in one benchmark run). The isotope-offset rules the project already rejected are
confirmed rejected against a proper control for the first time: every one of them is at or below a size- and
charge-matched random cut.

**The fix existed and was deleted unmeasured.** `assembly:dedup_precursors` / `dedupPrecursors_` (in the
source before 1.2.0, called with 15 ppm / 3.0 s / 0.01) is exactly this de-duplication, placed immediately after precursor inference, with a better tiebreak than
anything measured offline here: the higher MS1 trace apex intensity wins. It was pre-registered as a probable loss,
never run, and removed on 2026-09-08 in the pass that dropped 24 options
"with no measurement on record". Its help text carries a redundancy observation, not a falsification.

All measurements on PXD029836 run 1418 (120 min, 32 windows) with the released v1.3.1 binary, one extraction kept and
re-used. An **arm** is one filtered or merged copy of that extraction, searched with Sage and MSFragger in closed search
(the configuration of `bench/sage_closed.json`, ±25 ppm, fully tryptic) and counted the same way throughout: Sage,
distinct peptides at rank 1 with peptide q <= 0.01; MSFragger, the best hyperscore per peptide and modification, the
deepest target-decoy prefix at 1 %. **Exploratory: the arms were not timed, so no wall-time claim comes from them.** Scripts: `bench/iso_frontier.py` (price list), `bench/iso_rule.py` and
`bench/neigh_rule.py` (drop lists), `bench/mzml_drop.py` (exact spectrum filter), `bench/twin_merge.py` (the merge
rewriter), `bench/twin_origin.py` (origin classification), `bench/twin_sim.py` (content similarity),
`bench/funnel_dup.py` (the mechanism at the source), `bench/rho.py` (redundancy vs coverage).

## 1. The mechanism, counted at the source

`-diag:dump_ms1_tsv` gives every inferred precursor with its seed trace, before windows and before assembly:
**4,514,709 MS1 traces -> 1,145,141 precursors, and 1,145,141 distinct seed traces** -- one precursor per seed,
exactly. Grouping precursors by (mono m/z within 20 ppm, same charge, |dRT| <= 1.5 s, |d1/K0| <= 0.01):

| | precursors | share |
|---|---|---|
| no competitor | 619,387 | 54.1% |
| 1 competitor | 234,533 | 20.5% |
| 2 competitors | 128,910 | 11.3% |
| 3 competitors | 75,375 | 6.6% |
| 4 competitors | 41,635 | 3.6% |
| 5+ competitors | 45,301 | 4.0% |

Of the 577,208 resulting pairs, classified by what the two seed traces were:

| seed relation | pairs | share |
|---|---|---|
| one to three isotope steps apart, both walked down to the same mono | 370,078 | 64.2% |
| two distinct traces at the same m/z, RT and mobility | 204,469 | 35.4% |
| other | 2,661 | 0.5% |

**Why nothing stops it.** `prec.setMZ(pc.mono_mz)` (src/diaspextract.cpp:3038) emits the *inferred* monoisotope,
and the walk tries the seed plus up to three isotopes below it as the mono candidate (:2723). `findPartner` (:2620)
claims partners at nonzero offsets only and skips `used[j]` -- so `used[]` owns **traces**, and nothing owns a
**monoisotope position**. A seed whose best hypothesis resolves to an already-claimed mono still emits. This is why
the duplication kept measuring as "isotope" at ~22% in the emitted m/z: the larger, cleaner signal sits at offset
zero, where nobody looked.

**A second, independent source.** Adjacent diaPASEF windows on this run overlap by exactly 1.00 m/z, and precursor
selection at :4554 takes `[win_lo, win_hi]` closed at both ends with no cross-window ownership guard -- where the
tile path immediately below it has an exact ownership predicate and *throws* if two tiles claim one precursor
(:4565-4568). Classifying emitted twin pairs by isolation window and by the stamped RT cell edges
(`spx:tile_boundaries`, which drift up to 26 s from multiples of the nominal 600 s pitch): **85.8%** same window
and away from any cell edge, **13.2% cross-window**, **1.0%** cell cuts.
The cell grid is close to free; the window seam is not.

**Every linked pair shares its MS1 frame exactly.** Of the pairs the relation links within |dRT| <= 1.5 s -- 563,634
at d1/K0 <= 0.01, 354,462 at 0.005 -- every one has a text-identical RT, and none sits one frame apart
(`bench/twin_rt_exact.py`). So "same MS1 frame" is the same relation on this file, and it is the one 1.4.0 ships: it
can never link two RT cells, which keeps the fold invariant to tile grouping by construction.

## 2. The price list, and what it does and does not say

`bench/iso_frontier.py` prices 14 removal rules by peptides LOST at the fixed baseline threshold -- a peptide is
lost only when every spectrum identifying it is removed -- against a size- and charge-matched random cut.
Same-position twins are the only family that beats random; every isotope-offset rule is at or below it, and the
census-ceiling rule (k=1..3 heavy) is 0.52x, actively worse than random.

| rule | removed | % spectra | LOST | % peptides | random LOST | selectivity |
|---|---|---|---|---|---|---|
| same position, d1/K0 <= 0.005, dRT <= 1.5 | 236,720 | 22.5% | 1,529 | 2.66% | 3,317 | 2.17x |
| same position, d1/K0 <= 0.01, dRT <= 1.5 | 314,857 | 29.9% | 2,096 | 3.65% | 4,681 | 2.23x |
| same position, d1/K0 <= 0.01, dRT <= 3 | 480,911 | 45.7% | 4,150 | 7.22% | 8,386 | 2.02x |
| isotope M+1 heavy, d1/K0 <= 0.005 | 81,218 | 7.7% | 970 | 1.69% | 1,032 | 1.06x |
| isotope M+1 heavy, d1/K0 <= 0.02 | 150,408 | 14.3% | 2,354 | 4.10% | 2,003 | 0.85x |
| isotope k=1..3 heavy, d1/K0 <= 0.02 | 232,462 | 22.1% | 6,450 | 11.22% | 3,334 | 0.52x |

**What this table cannot say.** `lost()` counts a peptide only when *every* identifying spectrum is removed, so it
is invariant to *which* copy a rule deletes: the whole 2.23x is a statement that the rule concentrates on
well-covered peptides and under-samples thin ones, not that the copies were redundant. It also uses a different
peptide universe (spectrum-level q, 57,484) from the arms (peptide-level q, 50,800). Redundancy is established in
sections 3 and 4, not here.

## 3. Redundancy or coverage: rho

Of the identified spectra a rule removes, the fraction whose **retained** twin carries the **same** peptide
(`bench/rho.py`). The project's anchor for the both-identified chance level is 96% (isotope-linked pairs,
docs/ISOTOPE-DUPLICATION-2026-09-03.md:56-59).

| | d1/K0 <= 0.01 | d1/K0 <= 0.005 |
|---|---|---|
| identified spectra removed | 165,785 (52.7%) | 125,228 (52.9%) |
| retained twin, **same** peptide | 92.3% | 94.5% |
| retained twin identified, different peptide | 0.6% | 0.5% |
| retained twin exists, unidentified | 5.5% | 5.1% |
| **orphan** -- no retained twin at all | 1.6% | 0.0% (7 spectra) |
| conditional on both identified | **99.4%** | **99.5%** |
| singly-identified peptides losing their only spectrum | 13.5% | 9.9% |
| overall identified-spectrum removal rate | 43.2% | 32.6% |

The removed spectra are duplicate claims: where the retained twin is identified at all, it is the same peptide
99.4% of the time. The coverage effect is real too -- thin peptides are under-sampled ~3x -- but it is not the
explanation.

**The greedy pairwise rule does not keep a survivor per cluster.** On a chain it can remove both members of a pair
(1.6% of identified removals here, and its drop list depends on input row order). The union-find clustering in
`bench/twin_merge.py` keeps exactly one member per connected component by construction, which is one more reason
the merge form is the one to carry forward.

## 4. The arms

Nine spectrum-filter arms over the one kept extraction, both engines, counted as above. The baseline is the
unfiltered file; `identity` is the same re-serialisation every filtered arm gets (`index=` renumbered, the
`indexedmzML` wrapper dropped) with **nothing** removed, so file form is not a confound.

| arm | what it removes | spectra | Sage | vs base | MSFragger | vs base | lost per 1,000 removed |
|---|---|---|---|---|---|---|---|
| baseline | -- | 1,052,672 | 50,800 | -- | 45,447 | -- | -- |
| identity | nothing (re-serialised) | 1,052,672 | 50,800 | **0.00%** | 45,447 | **0.00%** | -- |
| merge, d1/K0 <= 0.005 | twin clusters, fragments unioned | 802,099 | 50,691 | -0.21% | 45,367 | -0.18% | **0.4** |
| merge, d1/K0 <= 0.01 | twin clusters, fragments unioned | 719,631 | 50,533 | -0.53% | 45,264 | -0.40% | **0.8** |
| drop, d1/K0 <= 0.005 | weaker twin dropped | 815,952 | 50,155 | -1.27% | 45,023 | -0.93% | 2.7 |
| drop, same clusters as merge | survivor keeps its own peaks | 719,631 | 49,834 | -1.90% | 44,676 | -1.70% | 2.9 |
| drop, d1/K0 <= 0.01 | weaker twin dropped | 737,815 | 49,863 | -1.84% | 44,739 | -1.56% | 3.0 |
| isotope M+1 heavy | the heavy member | 971,454 | 50,215 | -1.15% | 44,905 | -1.19% | 7.2 |
| random, size + charge matched | -- | 737,815 | 47,807 | -5.89% | 43,100 | -5.16% | 9.5 |
| off-lattice neighbour | co-eluting, co-mobile, not on the lattice | 895,869 | 48,654 | -4.22% | 43,562 | -4.15% | **13.7** |

Three things fall out, all on both engines:

1. **What the union is worth.** The merge arm and the drop arm on the same clusters remove the *identical* 333,041
   spectra. Unioning the fragments rather than discarding them is worth **+699 Sage and +588 MSFragger peptides**
   -- 72% and 76% of the drop's cost. The two copies are not copies: their peak counts differ in 69% of pairs
   (median 30 peaks), content cosine is 0.824 median against 0.695 for an off-lattice co-eluting neighbour, and
   pairs where only one member identifies sit at 0.651 -- complementary fragments, not a duplicate list.
2. **It is not "well-covered peptides survive".** The off-lattice neighbour arm removes spectra that are just as
   well covered but are not a second claim on the same monoisotope, and it costs **13.7 peptides per 1,000
   removed, worse than random (9.5)**. Having a co-eluting neighbour makes a spectrum more valuable; only the
   same-monoisotope relation marks it as redundant. (That population is also half the size: 156,803 candidates
   against the twin rule's 314,857.)
3. **The isotope rules stay rejected.** 7.2 peptides per 1,000 removed, against 2.7-3.0 for twin drops and 0.4-0.8
   for twin merges.

## 5. The error rate does not move

Entrapment (Sage against the target database concatenated with a foreign proteome, 20,416 target / 16,343 foreign
proteins, peptide-hypothesis ratio 0.6805,
scored with the shipped `bench/entrapment.py` estimator) on the baseline and on the arm that would ship:

| arm | accepted | entrapment hits | raw | corrected FDR | 95% CI |
|---|---|---|---|---|---|
| baseline | 48,313 | 319 | 0.66% | **0.98%** | [0.87, 1.09] |
| merge, d1/K0 <= 0.01 | 48,208 | 327 | 0.68% | **1.00%** | [0.89, 1.11] |

**+0.03 points**, against the standing release rule of <= +0.30. The nominal 1% is well calibrated on both, so the
merge arm's peptides are real and the small loss is a real loss -- not the threshold sliding on a thinned decoy
population, which is the failure mode `bench/fdrcheck.py` caught on the min_correlation sweep.

## 6. What this does NOT establish

- **The arms price a filter, not the fix.** Inference is a greedy claim loop over a monotone `used[]` bitmap
  (:2589-2596, :2626-2627, :2820-2821), so suppressing a precursor upstream frees its claimed traces and changes
  which *other* precursors form. Deleting a finished spectrum is exact for the arm and is not a model of
  `dedupPrecursors_` restored. The upstream change has to be measured on its own.
- **No wall-time or memory claim.** The arms were not timed.
- **One file.** 1418 only; the emission ratio and the window overlap are file properties.
- **The merge is a fragment union, not a re-assembly.** The principled version folds the duplicate seeds' XICs
  together *before* fragment correlation, which the offline union only approximates.

## 7. In the tool: 1.4.0rc1 (2026-09-21)

The fold shipped as `-assembly:twin_im_tolerance` (default 0.005), run per tile on the canonically sorted block with a
same-MS1-frame RT clause -- the relation measured here, since every linked pair shared its frame. On run 1418 it
equals `bench/twin_merge.py --im=0.005` for 802,099 of 802,099 spectra and every fragment; fold off equals the v1.3.1
spectra; entrapment FDR is 0.98 % [0.87, 1.09] with the fold at 0.005 against 0.98 % without (+0.01 points unrounded). On two
further runs (dataset D, see docs/BASELINE.md, and run 009 of PXD047793) the fold off reproduces the 1.3 output
digests exactly, and the fold on gives the same digest at every tile grouping and thread count tried.

## 8. Proposed on 2026-09-20 (superseded by what 1.4.0 shipped)

Done in 1.4.0 as the offline union moved into the tool (section 7), rather than as the precursor-level
restoration first proposed here; the window seam is folded by the same relation.

1. **Restore `dedupPrecursors_`** at its original site (after inference, before the isotope-support gate), and
   extend `Precursor_` to carry the merged members' trace indices so `assembleOne_` can correlate fragments against
   the summed XIC instead of one seed's. The drop form is worth -1.8% of peptides for -30% of spectra; the merge
   form is worth -0.5%. Placement relative to the IM charge veto matters -- the historical code ran *before* it --
   and the post-tile coverage assertion (:4695-4711) has to learn the same predicate or the run aborts.
2. **Give the window loop the ownership rule the tile loop already has** (:4554 vs :4565-4568): with a 1.00 m/z
   seam between adjacent windows, a precursor in the overlap is assembled twice. 13.2% of twin pairs.
3. Both are deterministic by construction if the grouping runs over the canonically sorted `precursors` vector
   (:3908) rather than inside the speculative loop.
