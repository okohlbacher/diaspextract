# Search-engine policy

**Develop on Sage. Refresh the baseline with both engines before any release or headline claim.**

## Measured basis (joint benchmark, 6 cells, one common FDR procedure)

| | value |
|---|---|
| peptide correlation | r = 0.972 |
| **rank correlation** | **Spearman ρ = 1.000** |
| PSM correlation | r = 0.9994 |
| ratio (sX/dT) correlation | r = 0.978 |
| **MSFragger offset on the ratio** | **−6.4 points**, sd 1.22 |
| **runtime** | Sage 133 s vs MSFragger 559 s → **4.2× faster** |

Sage never ordered two configurations differently from MSFragger.

## The rule

* **Iteration / screening → Sage alone.** 4.2× throughput, perfect rank agreement. Reject bad
  ideas cheaply.
* **Every status refresh → run BOTH**, and refresh the recorded
  baseline from that pair — the whole status table, not just the changed arm.
* **Never quote a Sage-only ratio as a headline.** MSFragger puts DIAspeXtract/the reference implementation **6.4
  points lower** on all three samples. Sage systematically flatters us.

## Why the rank agreement is not sufficient on its own

**Every paired point in that set spans a 14–61% difference.** Development changes are 1–5%, and
there is no paired data at that scale. The offset itself drifts **1.22 points** between samples —
comparable to the effect size a real change produces.

There is also a mechanism-specific risk. Sage's reported score is **discriminant-rescored**;
MSFragger's hyperscore is **raw**. A change producing many *marginal* PSMs can be exploited by
Sage's rescoring and vanish under MSFragger. `trace:ms1_split_valleys` — the current default — is
exactly that kind of change: it quadruples emitted groups, and most of the extra identifications
are marginal.

## How

`bench/joint_bench.py` scores both engines under **one** common target-decoy procedure against
**one** shared target+decoy FASTA (Sage's `generate_decoys` off), so the engine is the only
variable.

Documented asymmetry, not corrected: **no rescoring is applied to either side**, which favours
Sage on absolute counts. Adding Percolator to MSFragger alone would reintroduce the confound this
design removes; adding it to both is the correct fix and is not done.

Operational note: **MSFragger cannot read DIAspeXtract mzML** without
`bench/fix_mzml_cvparams.py` — OpenMS writes valueless cvParams, which is legal mzML but crashes
batmass-io.
