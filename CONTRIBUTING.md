# Contributing to DIAspeXtract

Thanks for your interest. A few things are specific to this project and will save you time.

## Before you start

**DIAspeXtract cannot be built against a released OpenMS package.** It includes a header that
`scripts/apply_openms_patches.sh` installs *into* the OpenMS tree, calls an OpenMS method the same
patch *adds*, uses APIs that postdate the newest OpenMS release, and refuses to compile against a
`MassTrace` without move operations (the third patch). You need an OpenMS source tree at the commit
pinned in `patches/openms.lock` (develop, 2026-04-24), patched by that script — the patches are cut
against exactly that commit. `.github/workflows/build.yml` does exactly this and is the reference
recipe.

**It needs a large machine.** At the shipped defaults (100 threads, one 600-s cell per tile) peak memory is about 22 GB on the 30-minute dataset D and about 21 GB on the 2-hour TNBC 009 (25 GB with `-perf:malloc_trim false`); in one tile (`-tile:cells_per_tile 0`) about 37 and 100 GB (docs/BASELINE.md). The test suite runs on
a synthetic input and needs none of that.

## The one rule that matters

**A change to extraction is judged on identified peptides, never on spectrum counts, and never on
one search engine alone.** This project has repeatedly found changes that improve one engine and
harm the other, and changes whose sign flips between acquisitions. The standing bar for altering a
default is:

- both search engines,
- entrapment FDR inside the previous interval,
- and confirmation on a second acquisition.

If a change is *supposed* to be output-neutral, prove it with a spectrum-list digest
(`bench/semantic_digest.py`) rather than with peptide counts — counts move by ~2% from a last-ulp
arithmetic difference, so equal counts are not equality.

## Running the tests

```bash
python3 test/test_diaspextract.py /path/to/diaspextract     # 25 end-to-end checks (26 where the build writes mzPeak), synthetic input
cmake -B build -DDIASPEXTRACT_TESTS_ONLY=ON && cmake --build build && ctest --test-dir build
```

The second form needs no OpenMS and is what CI runs on four platforms.

Which CI legs run depends on the repository. The two Linux legs of `tests.yml` run on every push and
pull request everywhere. The macOS and Windows legs, and the weekly OpenMS-from-source build of
`build.yml`, run on the public repository (`okohlbacher/diaspextract`) and on a manual
`workflow_dispatch`; the private development repository runs Linux only on push, because its Actions
minutes are paid and the release gate is the public repository's CI on the exported code, green on all
legs before a release is tagged. To run the full matrix on the development repository, dispatch
`tests.yml` by hand.

Run the end-to-end suite from a shell with no `DIASPEXTRACT_*` variable set: it passes the caller's
environment through to the tool (check 22 strips only the two old prefixes it tests), so a stray probe
variable would configure every check. CI runs it under `env -i`, and the cluster publish script starts
from a clean environment.

## Pull requests

- Say what you measured, on what, and what would have falsified it. "Should be faster" is not a
  measurement; a wall-clock delta from two runs on a shared machine usually isn't either.
- Keep defaults reproducible: a run with `-in`, `-out` and `-threads` must reproduce the benchmarked
  configuration. There is a test for this, because it has been broken before.
- Don't commit acquisition or specimen identifiers; keep sample paths and names in untracked local files.

## Reporting a bug

Include the DIAspeXtract version, the OpenMS version, the full command line, and the `spx:*` userParams
from the output mzML — they record which detector and which calibration actually ran, which is
usually the answer.
