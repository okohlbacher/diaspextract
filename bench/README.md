# bench/

Small standard-library tools for checking DIAspeXtract output. Each script's docstring has its usage.

- `semantic_digest.py`: sha256 of the `<spectrumList>`; the output-identity gate (the header is ignored).
- `mzml_header_diff.py`: diffs everything outside the spectrum list, completion time masked; the header check the digest cannot do.
- `mzml_specdiff.py`: per-spectrum diff of two outputs, to locate a digest difference.
- `entrapment.py`: entrapment FDR estimator (foreign-proteome false positives, bootstrap CI).
- `fix_mzml_cvparams.py`: adds `value=""` to valueless cvParams so MSFragger can read OpenMS mzML.
- `tile_concat.py`: concatenates RT-tiled outputs, keeping each tile's core precursors.
- `sage_closed.json`: a closed-search Sage configuration (set `fasta` before use).
