---
name: Bug report
about: Something went wrong during extraction
labels: bug
---

**Command line** (the complete one):

```
diaspextract -in ... -out ... -threads ...
```

**Versions**: DIAspeXtract / OpenMS / OS / compiler

**Provenance from the output mzML** — these record what actually ran, and are usually the answer:

```
spx:detector, spx:require_isotope_support, spx:mz_calibration
```

**What happened**, and what you expected instead. Include the tail of the run log; it names the
calibration used and the detector selected, and warns when it falls back.

**Machine**: cores and RAM. Peak memory is 68–110 GB on real acquisitions, and an out-of-memory
kill can look like an unrelated failure.
