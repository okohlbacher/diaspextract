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
python3 bench/slice_d.py /path/to/run.d /tmp/slice.d 3000 3060   # once: one minute of a public run as a real .d
DIASPEXTRACT_TEST_D=/tmp/slice.d python3 test/test_diaspextract.py /path/to/diaspextract
cmake -B build -DDIASPEXTRACT_TESTS_ONLY=ON && cmake --build build && ctest --test-dir build
```

The end-to-end suite has 32 checks (33 where the build writes mzPeak). All but one run on small synthetic inputs; the
last runs a real diaPASEF `.d` through the default path, and is skipped without `DIASPEXTRACT_TEST_D`. CI gives it the
5.6-minute HeLa test run that OpenMS publishes. The `ctest` line needs no OpenMS. CI runs it on Linux, macOS and Windows for every pull request.

The end-to-end checks need OpenMS. `build.yml` builds the patched OpenMS and the tool on Linux (x64, arm64), macOS
(arm64, x64) and Windows, runs `ctest` and the end-to-end suite on each, packages the binary with its libraries and
runs the package with nothing from the build. It runs on the public repository for every push to main and every pull
request, weekly, on release tags and on demand. Its first run after a change to the OpenMS pin, the patches or a
dependency rebuilds OpenMS, one to three hours per platform; later runs take the cached build. Run the checks yourself
before you open a pull request.

`DIASPEXTRACT_TEST_D_DIGEST` pins that run's spectra digest: when set, the last check compares the output with it
(the check prints the digest when it is unset). CI takes it from the repository variable of the same name, set from
the first green build's log; all five platforms must give the same digest. A change that alters the spectra updates
it together with its CHANGELOG entry.

## Releasing

Set `project(VERSION)` in `CMakeLists.txt`, add a `## vX.Y.Z` section to `CHANGELOG.md`, and push the tag `vX.Y.Z`.
`build.yml` refuses a tag that does not match the version. It signs and notarizes the macOS binaries (it fails without
the `MACOS_*` secrets), signs the Windows executable through SignPath when `SIGNPATH_API_TOKEN` is set, attaches the
seven packages and `SHA256SUMS` to the release, creating it from the CHANGELOG section if it does not exist, and asks
the Homebrew tap to follow. Before a release that changes the packaging, run the workflow by hand with
`release_dry_run`: it takes the same path and returns the packages as workflow artifacts.

Run the end-to-end checks from a shell with no `DIASPEXTRACT_*` variable set other than `DIASPEXTRACT_TEST_D` and
`DIASPEXTRACT_TEST_D_DIGEST`. The
suite passes your environment to the tool, so a stray variable changes every check.

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
