# Benchmarks

## `sage_closed.json` — the closed metric (what every README number came from)

Byte-for-byte the configuration behind every `peptide_q` figure in the repo.

What it structurally cannot see: `precursor_tol` is ±25 ppm. A charge error fabricates
ΔM = (z′−z)(m/z − 1.00728) — hundreds of Da — so those PSMs are never reported. **Closed search
cannot score the failure mode this tool exists to avoid.**

## `open_bench.py` — entrapment-controlled open search

Built because the closed benchmark rewards the opposite of the goal: we emit 4.6 PSMs/peptide
against the reference implementation's 1.46, so redundancy is free score under `peptide_q`. That one fact predicts
in advance every result observed — each redundancy-reducing lever losing, and `ms1split` winning
while simultaneously raising co-isolation 52% and dropping charge below baseline. *Six
falsifications sharing one biased instrument are one falsification.*

Three measurements, ordered by how much they assume:

1. **Delta-mass geometry** — no ground truth, no decoy, no entrapment. Charge errors are
   *linear in m/z* and form **rays** through (1.00728, 0) with integer slopes; real
   modifications are constant shifts and form **horizontal bands**. Separable by construction.
   Validated on synthetic data: 300/300 slope −2 rays recovered, 125× over the chance floor,
   pure-noise control at 0.0000 vs a 0.0123 floor. **Reports its own chance floor.**
2. **Entrapment FDR** — 16,343 *A. thaliana* proteins vs 20,416 human (ratio 0.800). Independent
   of the decoy model, which degrades as open search grows the candidate space ~10³×
   (Kong 2017; Chick 2015; Wen et al., *Nat Methods* 2025).
3. **Modified peptides at ≤1% entrapment FDR** — the primary metric.

### Two windows, deliberately

`sage_open.json` is −150/+500 Da (primary). `sage_open_wide.json` is −3000/+1000 Da and exists
**only** as the positive control: a 4→2 error at m/z 500–1400 lands at −1000 to −2800 Da and is
**not reported at all** inside ±500. A benchmark run only at ±500 would be blind to the very
failure it was built to detect.

### Pre-registered falsification — written before the first run

Closed ratio `split_count/base` on dataset B = 8,411/5,817 = **1.446**.

| modified-peptide ratio | conclusion |
|---|---|
| **≥ 1.446** | "closed benchmark rewards shotgunning" is **falsified**; defaults vindicated for both purposes; the six falsified levers stand; stop reducing redundancy |
| **1.0 – 1.446** | real, but overstated by the closed metric |
| **< 1.0** | shipped defaults are **actively harmful** for the stated purpose |

Replicate spread is 0.17% (8,411 / 8,408 / 8,406 / 8,397), so differences under ~0.2% are noise.

**Positive control, non-negotiable.** The `envelope` arm has 10,199 known 4→2 confusions. They
must appear as slope −2 ray mass in the wide window and be markedly weaker in `count`. **If the
control does not fire, no other number from this script may be quoted.**

### Run

```bash
python3 open_bench.py fasta --target human.fasta --foreign arabidopsis.fasta --out human_entrap.fasta
./run_open.sh
python3 open_bench.py analyse RESULTS.tsv --label arm --ratio 0.800
```

All three arms are extracted on **one node with one binary** — comparing across nodes or builds
is the error class this benchmark exists to eliminate.


## `joint_bench.py` — the success metric across the three datasets

**3 samples × {the reference implementation, DIAspeXtract} × {Sage, MSFragger}**, reporting the two numbers the
project is scored on:

* **PSMs** at ≤1% PSM-level FDR
* **unique peptides** at ≤1% peptide-level FDR

### The design decision that matters

Both engines are scored by **one identical target-decoy procedure** implemented in
`joint_bench.py`, *not* by each engine's own FDR machinery (Sage's internal linear discriminant;
MSFragger via PeptideProphet/Philosopher or Percolator). Those pipelines differ in rescoring, in
PSM→peptide aggregation, and in decoy handling — comparing their outputs would confound *"which
engine identifies more"* with *"whose FDR is more permissive"*.

For the same reason both engines search **one pre-built target+decoy FASTA** with a shared
`rev_` prefix, and Sage's internal decoy generation is turned **off**, so the two see
byte-identical search spaces.

Peptide-level FDR is "picked": best-scoring PSM per distinct sequence, then target-decoy over
that reduced list. That axis is the point — DIAspeXtract emits ~4.6 PSMs/peptide against
the reference implementation's 1.46, so PSM counts reward redundancy and peptide counts do not. **Both are
reported so the gap is visible rather than hidden.**

Validated before use: null control (targets and decoys from one distribution) collapses to 1
accepted; peptide-level correctly reduces 5 PSMs over 2 sequences to 2.

### Known asymmetry — read before quoting an engine comparison

**No rescoring is applied.** Sage's reported score is already discriminant-rescored; MSFragger's
hyperscore is not, which favours Sage on absolute counts. Adding Percolator to MSFragger alone
would reintroduce the confound this design removes; adding it to *both* is the correct fix and
is not done here.

**The tool comparison is unaffected** — the reference implementation vs DIAspeXtract is made *within* an engine, and
that is the comparison this harness exists for. The engine axis is a secondary readout.

### Run

```bash
python3 joint_bench.py decoydb --target human.fasta --out human_decoy.fasta
./run_joint.sh
python3 joint_bench.py collate joint_results.json
```
