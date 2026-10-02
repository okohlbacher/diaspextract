# Contributing to DIAspeXtract

Thanks for your interest. This page covers what is specific to this project.

## Before you start

**You need a patched OpenMS.** DIAspeXtract does not build against a released OpenMS package. It needs
OpenMS at the commit pinned in `patches/openms.lock`, patched by `scripts/apply_openms_patches.sh`.
[README.md](README.md#install) lists the steps, and `.github/workflows/build.yml` is the reference recipe.

**You need a large machine for real data.** The Requirements table in [README.md](README.md#requirements) gives
the memory that real runs need, and [docs/options.md](docs/options.md#speed-and-memory) lists the settings that
lower it. The tests use a small synthetic input and run on any machine.

## Changing a default

Judge a change on identified peptides, not on spectrum counts. You may screen ideas with Sage alone, because it
is fast. Before you change a default or quote a gain, show all three:

- more identified peptides with both Sage and MSFragger;
- entrapment FDR inside the previous interval (`bench/entrapment.py`);
- the same result on a second public run.

If a change should leave the output alone, prove it with `bench/semantic_digest.py`: the digests must match.
Equal peptide counts do not prove equal output, because a last-bit arithmetic difference can move them.
[bench/README.md](bench/README.md) describes these scripts, including the one MSFragger needs to read the output.

## Data

Use only public data sets from ProteomeXchange (PXD accessions) for benchmarks, tests and development, so that
anyone can reproduce every number. Do not add, stage or cite unpublished data.

## Running the tests

```bash
python3 test/test_diaspextract.py /path/to/diaspextract     # 29 end-to-end checks (30 where the build writes mzPeak), synthetic input
cmake -B build -DDIASPEXTRACT_TESTS_ONLY=ON && cmake --build build && ctest --test-dir build
```

The second command needs no OpenMS. CI runs it on Linux, macOS and Windows for every pull request. The
end-to-end checks need OpenMS, so CI runs them only in the full build (`build.yml`: weekly, on release tags and on
demand). Run them yourself before you open a pull request.

Run the end-to-end checks from a shell with no `DIASPEXTRACT_*` variable set. The suite passes your
environment to the tool, so a stray variable changes every check.

## Pull requests

- Say what you measured, on which data, and what result would have proven you wrong. "Should be faster" is not
  a measurement. Neither, usually, is the wall-time difference of two runs on a shared machine.
- Keep the defaults reproducible: a run with only `-in`, `-out` and `-threads` must give the benchmarked
  configuration.
- If you add, change or remove an option, update its help text, [docs/options.md](docs/options.md) and any
  mention in README.md. CI fails when README.md or this page names an option the tool does not have. If the
  change alters the spectra, add a CHANGELOG entry.
- Do not commit acquisition or specimen identifiers. Keep sample paths and names in untracked local files.

## Reporting a bug

Include the DIAspeXtract version, the OpenMS version, the full command line, and the `spx:*` userParams from the
header of the output mzML. They record the settings that actually ran, which usually explains the problem.
Report security issues as [SECURITY.md](SECURITY.md) describes, not in a public issue.
