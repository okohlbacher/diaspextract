# A global discrete m/z axis: design

**Origin.** User proposal, 2026-09-03: a TOF instrument measures integer flight times, so m/z is
discrete at the source. Define one global m/z axis and a spectrum becomes a vector of intensities.
This follows the global RT axis (already landed) and applies the same idea to the other dimension.

## 1. What is already discrete, and what is thrown away

The compact peak store is *already* on a uniform integer m/z grid:

    constexpr double MZ_Q = 1e5;                  // m/z quantum = 1e-5 Da
    f.mzq.push_back((uint32_t)llround(mz * MZ_Q));

10 bytes per peak: `uint32` m/z index, `float` intensity, `uint16` mobility. 0.007 ppm at m/z 1400
against a 20 ppm trace tolerance, so the quantisation is 2,800x finer than anything downstream.

And then `materializeWindow()` throws it away:

    p.setMZ((double)f.mzq[i] / MZ_Q);

because the next stage is OpenMS `MassTraceDetection`, which takes a `PeakMap` of `double` m/z. So
every peak is converted back to a double, re-sorted, and matched by a parts-per-million tolerance
with a binary search per comparison -- on data that arrived as integers. That conversion is the
boundary this design removes.

## 2. The axis is the flight-time index

**Decision (user, 2026-09-03): the native flight-time index, not a uniform m/z grid.** The tool only
runs on ion-mobility TOF data, so the native axis is always available and there is no reason to
approximate it. The earlier draft of this document chose a uniform 1e-5 Da grid to keep one code
path for hypothetical non-TOF input; that was the wrong trade.

| | uniform m/z grid | **native flight-time bin** |
|---|---|---|
| index | `round(mz * 1e5)` | the instrument's TOF bin |
| size | ~1.7e8 over m/z 100-1700 | **~5e5** |
| lossless? | no -- a second discretisation on top of the instrument's own | **yes: the bin IS the measurement** |
| same ion across frames | different m/z as the digitizer temperature drifts | **same bin** |
| tolerance | a ppm window in a derived quantity | **a handful of integer bins** |

Three reasons the native axis is better in kind, not merely smaller:

1. **It is lossless.** No quantum has to be justified, because there is no quantisation. The
   compact store's 1e-5 Da grid was a coarse approximation of a *derived* quantity.
2. **The same ion gives the same bin in every frame.** A calibrated m/z wobbles with the digitizer
   temperature -- that per-frame factor is exactly what `TdfMzCalibration::frameFactor` models, and
   it is why recovering the exact calibration was worth 12% of identifications. In m/z space part
   of the ppm tolerance is spent absorbing that drift. In bin space it is not spent at all.
3. **An index of 5e5 is a usable subscript.** 1.7e8 is not.

`diaspextract::TdfMzCalibration` already provides both directions, `tofToMz(tof, b)` and
`mzToTof(mz, b)`, plus `frameFactor(T1)`. So the axis is:

```cpp
using TofIdx = uint32_t;
struct TofAxis {
  diaspextract::TdfMzCalibration cal;
  vector<double> b_by_frame;     // per vendor frame Id
  double factor(size_t frame_id) const;
  double mzOf(TofIdx tof, double b) const;      // report time only
  TofIdx tofOf(double mz, double b) const;      // load time only
};
```

**The extraction never needs a calibrated m/z at all.** Grouping, gating and correlation all work on
bins; the calibration is applied once per trace, at report time, with the apex frame's own factor --
which is also the most accurate place to apply it.

**Fail closed.** Without the vendor calibration there is no flight-time axis, and `trace:detector=
integer` refuses to run rather than inventing one. Inventing one would silently reintroduce the
approximation the detector exists to avoid.

**The loader does not have to change.** The compact store's 1e-5 Da quantum is ~400x finer than one
TOF bin at m/z 600, so inverting the calibration on a stored value recovers the original bin
exactly. The conversion happens once, when a window's slab is built.

## 3. The data structures

```cpp
/// The m/z axis is global, uniform and integer: index i means m/z = i / MZ_Q. Stored nowhere,
/// because it is an affine function of the index -- unlike the RT axis, which has to be a table.
using TofIdx = uint32_t;                      // the instrument's flight-time bin
// mzOf / tofOf live on TofAxis (section 2): they need the frame's calibration factor.
```

A frame is a row of `PeakSlab`'s frame table. A **peak** never becomes a `Peak1D` again:

```cpp
/// One window's peaks, frame-major, in index space. Replaces the materialised PeakMap.
struct PeakSlab
{
  vector<uint32_t> frame_off;   ///< size F+1; frame f owns [frame_off[f], frame_off[f+1])
  vector<TofIdx>   tof;         ///< flight-time bin, ascending within each frame
  vector<double>   b;           ///< size F: the frame's calibration factor
  vector<float>    inten;
  vector<uint16_t> imq;         ///< the same uint16 mobility quantum as the compact store
  vector<uint32_t> rt_index;    ///< size F: each frame's index on the GLOBAL RT axis
};
```

Per-peak cost is unchanged at 10 bytes, but the ~20 B/peak `Peak1D` + `FloatDataArray`
materialisation disappears, along with its per-`MSSpectrum` overhead. On dataset A that materialisation is
the largest single per-window allocation.

A **trace** keeps the shape it has after the RT-axis change, with m/z now an index too:

```cpp
struct Trace
{
  double   mz;          ///< reported m/z: the apex BIN calibrated with the apex frame's factor
  uint32_t rt;          ///< apex frame, on the global RT axis
  uint16_t imq;
  float    intensity;
  vector<uint32_t> xi;  ///< frame indices, ascending  (already landed)
  vector<float>    xv;
};
```

`double mz`, `double rt`, `double im` per trace become 4+4+2 bytes. At 6.4 M MS1 traces and ~8.5 M
fragment traces per window this is worth having, but it is not the point: the point is what it does
to the loops.

## 4. What the integer axis buys in the loops

**Trace detection becomes bucketing, not searching.** The current detector accepts a peak if it is
within `mass_error_ppm` of the trace's running centroid, which costs a binary search per candidate
and makes the cost superlinear in peaks per band (measured: 4 bands 18:25 wall / 163 GB against 12
bands 12:53 / 103 GB on dataset D). On an integer axis, peaks of one ion across frames differ by at most
`ppmSpan(mz, tol)` index units, so:

  * bucket peaks by `mz >> SHIFT` with `1 << SHIFT` chosen just above the largest tolerance span;
  * a trace's candidates are in its own bucket and the two neighbours -- O(1), no search;
  * the running centroid stays a `double` internally, so the acceptance decision is unchanged.

**Fragment-to-precursor m/z tests become integer comparisons.** The scorer excludes the precursor
peak (`fabs(f.mz - pc.mono_mz) < 0.01`) and optionally its isotopes; assembly dedups at 10 ppm.
These become integer differences against a precomputed span.

**Banding can be made exact, or provably not.** The band edges become index cut points and the halo
becomes an exact index count. Whether that makes banding exact depends on whether a trace's total
index span is boundable -- the drifting-centroid question. If it is
not, the honest outcome is a bound enforced by construction (cap a trace's span in index units), and
banding becomes exact *by definition of the detector* rather than by hope.

## 5. Refactoring run, in order, each step compiling and testable

1. **Axis primitives + `Trace` in index space.** Mechanical; `mzOf`/`idxOf` at every boundary.
   Output must be byte-identical: the values are the same doubles, just reconstructed.
2. **`PeakSlab` and a slab-based window materialiser**, used only by a new detector; the OpenMS
   path stays live behind `trace:detector=openms|integer` so the two can be A/B'd on the same build.
3. **The integer-axis detector**: bucket, extend, split. This is the step that changes results, so
   it is gated and compared on identified peptides with both engines, never on a digest.
4. **Scorer and assembly in index space**, removing the remaining double comparisons.
5. **Banding on index cut points**, with the span bound made explicit.

Steps 1, 2, 4 are output-neutral and can be verified by the semantic digest. Step 3 is a new
algorithm and must clear the both-engines rule. Step 5 depends on whether the span bound of step 4
holds.

## 6. How it is to be measured

Against the current head, on dataset D and dataset A, 100 threads, same node:

| metric | why |
|---|---|
| wall time | the headline; the reference implementation is 9:46 on dataset A against our 22:39 |
| peak RSS | 152 GB is the number that forces the memory admission gate |
| window-loop parallel speedup | occupancy; 67.9x of 100 before the task pool |
| stage breakdown | whether the gain is where the design predicts (trace, 47%) |
| Sage + MSFragger peptides | steps 1/2/4 must not move them at all; step 3 must not lose them |
| entrapment FDR | any peptide gain must not be bought with false positives |

**Falsifier for the whole design:** if steps 1, 2 and 4 land byte-identical but neither runtime nor
peak memory improves by more than the ~2% replicate spread, then the double conversions were never
the cost and the remaining time is genuinely in the detector's arithmetic -- in which case step 3 is
the only thing worth doing, and it should be judged on its own.
