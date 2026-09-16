## What this changes

<!-- One or two sentences. -->

## Evidence

<!-- Delete whichever does not apply. -->

**Output-neutral change** — spectrum-list digest is identical:

```
bench/semantic_digest.py before.mzML after.mzML
```

**Change that moves identifications** — both engines, and a second acquisition:

| | before | after |
|---|---|---|
| engine 1, dataset 1 | | |
| engine 2, dataset 1 | | |
| engine 1, dataset 2 | | |
| engine 2, dataset 2 | | |
| entrapment FDR | | |

## Checks

- [ ] `python3 test/test_diaspextract.py <binary>` passes
- [ ] A run with only `-in`, `-out` and `-threads` still reproduces the benchmarked configuration
- [ ] No acquisition or specimen identifiers added to tracked files
