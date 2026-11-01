# bench/

Small standard-library tools for checking DIAspeXtract output. Each script's docstring has its usage.

- `semantic_digest.py`: sha256 of the `<spectrumList>`; the output-identity gate (the header is ignored).
- `mzml_header_diff.py`: diffs everything outside the spectrum list, completion time masked; the header check the digest cannot do.
- `mzml_specdiff.py`: per-spectrum diff of two outputs, to locate a digest difference; spectra are paired by RT,
  precursor m/z, charge, 1/K0 and isolation window.
- `entrapment.py`: entrapment FDR estimator (foreign-proteome false positives, bootstrap CI). The peptide-space ratio
  follows the Sage config's enzyme and bounds, digested as Sage does (N-terminal cleavage with `c_terminal: false`,
  semi-enzymatic, no cleavage with `cleave_at: "$"`, non-specific with `cleave_at: ""`, refused when too large to
  enumerate); results are reused only for the same inputs.
- `fix_mzml_cvparams.py`: adds `value=""` to valueless cvParams so MSFragger can read OpenMS mzML, and rebuilds the
  index for the moved bytes.
- `tile_concat.py`: concatenates RT-tiled outputs, keeping each tile's core precursors.
- `slice_d.py`: cuts a retention-time range out of a Bruker `.d` into a small, self-contained `.d` (the end-to-end
  suite's real input; quick A/B screens).
- `sage_closed.json`: a closed-search Sage configuration (set `fasta` before use).
