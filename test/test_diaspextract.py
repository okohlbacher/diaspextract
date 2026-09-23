#!/usr/bin/env python3
"""End-to-end checks on a synthetic input. No framework, no fixtures: one file, plain asserts.

Verification on this project has been whole-file cluster benchmarks, which catch regressions in
aggregate but cannot say WHY a number moved, and cannot run before a commit. These are the smallest
checks that fail if the things most likely to break do break:

  1. isotope ownership   one z=2 envelope (M, M+1, M+2) must yield ONE precursor at the mono, not
                         three, and the M+1/M+2 peaks must not reseed precursors of their own
  2. min_charge          a supported z=1 envelope is emitted at the default (1) and dropped at 2, while a
                         z=2 envelope is unaffected either way
  3. parameter bounds    a negative count must be REJECTED, not wrapped through Size into an
                         empty-but-successful run
  4. determinism         the same input at 1 and 4 threads must give byte-identical spectrum data
  5. trace extension     one missing cycle must leave ONE trace and its gapped fragment in the
                         spectrum; seven missing cycles must give TWO traces. Every version of the
                         integer detector that lost 92% of the peptides passed checks 1-4.
  ... 14. writer           count, id="spectrum=K", index="K" and the <indexList> offsets agree
  ... 19. MS1 prune        perf:ms1_prune leaves the spectra byte-identical; dropping its witnesses does not
     19b. prune chain      a band kept alive only by interior chain links: dropping just the links changes it
     19c. stored order     MS1 spectra stored in descending m/z give the ascending fixture's result, prune off or on
     19d. non-finite m/z   a NaN or infinite MS1 m/z is refused right after the pick, before the prune, prune off or on
      20. dnoise:ms1       on by default, it filters a Bruker .d only: an mzML is skipped, logged, stamped, unchanged
      21. dnoise SDK       a .d holding only its analysis.tdf: with a Bruker SDK path set, dnoise:ms1 refuses before the
                           load; without one the run gets past the dnoise setup and fails in the loader
      22. old variables    a set SPEXTRACTOR_* (SpeXtractor) or DIASPEXTRACTOR_* (DIAspeXtractor) variable refuses the run
                           before anything, naming it and its DIASPEXTRACT_ name; DIASPEXTRACT_DET runs and prints its
                           [det] lines

These exercise the default (OpenMS) detector only. `trace:detector=integer` refuses to run
without the vendor flight-time calibration, which a synthetic mzML cannot carry, so it is judged on
real data in the benchmark harness and not here.

An mzML is never streamed: every check here takes the RESIDENT load (the stream attempt fails on it
before any frame is delivered), so the MS1 prune inside the streaming consumer (flushMS1_, its
cross-batch sample offsets) is covered only by the cluster gates on real .d input.

Usage: test_diaspextract.py /path/to/diaspextract [workdir]
Exit status is 0 only if every check passes.
"""
import base64, hashlib, json, os, re, struct, subprocess, sys, tempfile

ISO = 1.0033548
PROTON = 1.007276


def b64f(vals, double=True):
    fmt = "<%d%s" % (len(vals), "d" if double else "f")
    return base64.b64encode(struct.pack(fmt, *vals)).decode()


def spectrum(idx, ms_level, rt, peaks, im, prec=None, group=0):
    """One mzML spectrum. `peaks` is [(mz, intensity)], `im` one 1/K0 per peak."""
    mz = b64f([p[0] for p in peaks])
    inten = b64f([p[1] for p in peaks], double=False)
    imarr = b64f(im)
    pre = ""
    if prec:
        target, lo, hi = prec
        pre = f"""<precursorList count="1"><precursor><isolationWindow>
<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{target}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{target-lo}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{hi-target}" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
</isolationWindow><activation><cvParam cvRef="MS" accession="MS:1000044" name="dissociation method"/></activation></precursor></precursorList>"""
    # "frame=<N>" is the vendor key frameIdOf() parses; without it every frame is unmappable and
    # the integer detector falls back, which is how all 8 checks silently ran on the OTHER detector.
    return f"""<spectrum id="frame={idx + 1} windowGroup={group} scan=0" index="{idx}" defaultArrayLength="{len(peaks)}">
<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum"/>
<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{ms_level}"/>
<cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum"/>
<scanList count="1"><cvParam cvRef="MS" accession="MS:1000795" name="no combination"/>
<scan><cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{rt}" unitAccession="UO:0000010" unitName="second" unitCvRef="UO"/></scan></scanList>
{pre}
<binaryDataArrayList count="3">
<binaryDataArray encodedLength="{len(mz)}"><cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{mz}</binary></binaryDataArray>
<binaryDataArray encodedLength="{len(inten)}"><cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000521" name="32-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{inten}</binary></binaryDataArray>
<binaryDataArray arrayLength="{len(im)}" encodedLength="{len(imarr)}"><cvParam cvRef="MS" accession="MS:1002816" name="mean inverse reduced ion mobility array" unitAccession="MS:1002814" unitName="volt-second per square centimeter" unitCvRef="MS"/><cvParam cvRef="MS" accession="MS:1000523" name="64-bit float"/><cvParam cvRef="MS" accession="MS:1000576" name="no compression"/><binary>{imarr}</binary></binaryDataArray>
</binaryDataArrayList></spectrum>"""


def synth(path, n_cycles=14, cycle_s=1.4, gap_cycles=(), groups=None, z1_shift=0, drop_ms2=None, ms1_extra=None,
          ms1_descending=False):
    """Two precursors that co-elute in one isolation window: a z=2 envelope and a z=1 envelope.

    Both are given a 3-peak isotope envelope so `require_isotope_support` keeps them, matching
    fragments that co-elute with the precursor, and a distinct ion mobility.

    `gap_cycles`: cycles in which the z=2 precursor AND its first fragment are ABSENT. This is
    the input that separates a correct trace-extension rule from a wrong one: OpenMS tolerates up
    to five consecutive missing frames before it ends a trace, so one gap must still give ONE
    trace, and seven consecutive gaps must give TWO. Every version of the integer detector that
    lost 92% of peptides passed the other checks in this file; none would have passed this.

    `groups`: two WindowGroup ids; the MS2 frames then alternate between them by cycle,
    which is how one diaPASEF scheme (PXD017703 "py3") acquires the SAME m/z window twice per
    cycle with shifted ion-mobility slices. They must become two windows, not one.

    `z1_shift`: cycles by which the z=1 precursor's elution apex is moved (negative = earlier),
    so the two precursors can be owned by different tiles. `drop_ms2`: (lo, hi) in seconds; MS2
    frames with lo <= rt < hi are omitted (a window with no frames in one tile).

    `ms1_extra`: a function of the cycle returning extra MS1 peaks [(mz, intensity, 1/K0)].
    `ms1_descending`: store every MS1 spectrum in DESCENDING m/z order, the IM array in lockstep.
    """
    win = (600.0, 590.0, 620.0)                    # target, lo, hi
    z2_mono, z2_im = 601.3, 0.90                   # a doubly-charged precursor
    z1_mono, z1_im = 610.6, 1.20                   # a singly-charged one
    z2_frags = [(233.11, 0.7), (348.19, 1.0), (461.27, 0.55)]
    z1_frags = [(288.14, 0.8), (401.22, 0.9)]
    specs, idx = [], 0
    for c in range(n_cycles):
        rt = 1.0 + c * cycle_s
        shape = 1.0 - abs(c - n_cycles / 2.0) / (n_cycles / 2.0 + 1.0)   # a triangular elution
        amp = max(shape, 0.05)
        shape1 = 1.0 - abs(c - (n_cycles / 2.0 + z1_shift)) / (n_cycles / 2.0 + 1.0)
        amp1 = max(shape1, 0.05)
        ms1, ms1_im = [], []
        gap = c in gap_cycles
        for k, rel in enumerate((1.0, 0.55, 0.20)):                      # averagine-ish envelope
            if not gap:
                ms1.append((z2_mono + k * ISO / 2.0, 6.0e5 * amp * rel)); ms1_im.append(z2_im)
            ms1.append((z1_mono + k * ISO / 1.0, 4.0e5 * amp1 * rel)); ms1_im.append(z1_im)
        for mz, it, im in (ms1_extra(c) if ms1_extra else ()):
            ms1.append((mz, it)); ms1_im.append(im)
        order = sorted(range(len(ms1)), key=lambda i: ms1[i][0])
        if ms1_descending: order.reverse()
        specs.append(spectrum(idx, 1, rt, [ms1[i] for i in order], [ms1_im[i] for i in order])); idx += 1
        if drop_ms2 and drop_ms2[0] <= rt < drop_ms2[1]: continue
        ms2, ms2_im = [], []
        for fi, (mz, rel) in enumerate(z2_frags):
            if gap and fi == 0: continue                                   # first fragment gapped too
            ms2.append((mz, 2.0e5 * amp * rel)); ms2_im.append(z2_im)
        for mz, rel in z1_frags: ms2.append((mz, 1.5e5 * amp1 * rel)); ms2_im.append(z1_im)
        order = sorted(range(len(ms2)), key=lambda i: ms2[i][0])
        specs.append(spectrum(idx, 2, rt, [ms2[i] for i in order], [ms2_im[i] for i in order],
                              prec=win, group=groups[c % 2] if groups else 0)); idx += 1
    body = "\n".join(specs)
    open(path, "w").write(f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<indexedmzML xmlns="http://psi.hupo.org/ms/mzml"><mzML version="1.1.0" id="diaspextract_test">
<cvList count="1"><cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/></cvList>
<fileDescription><fileContent><cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum"/></fileContent></fileDescription>
<softwareList count="1"><software id="so_test" version="0"><cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="diaspextract-test"/></software></softwareList>
<instrumentConfigurationList count="1"><instrumentConfiguration id="ic_test"><cvParam cvRef="MS" accession="MS:1000031" name="instrument model"/></instrumentConfiguration></instrumentConfigurationList>
<dataProcessingList count="1"><dataProcessing id="dp_test"><processingMethod order="0" softwareRef="so_test"><cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML"/></processingMethod></dataProcessing></dataProcessingList>
<run id="run_test" defaultInstrumentConfigurationRef="ic_test">
<spectrumList count="{len(specs)}" defaultDataProcessingRef="dp_test">
{body}
</spectrumList></run></mzML></indexedmzML>
""")
    return z2_mono, z1_mono


def mzcal_table(db):
    """Write dataset D's MzCalibration row (tests/calibration_golden.json) into db; return the golden entry."""
    g = json.load(open(os.path.join(os.path.dirname(__file__), "..", "tests", "calibration_golden.json")))[0]
    db.execute("CREATE TABLE MzCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, DigitizerTimebase REAL,"
               " DigitizerDelay REAL, C0 REAL, C1 REAL, C2 REAL, T1 REAL, dC1 REAL, dC2 REAL, C3 REAL, C4 REAL)")
    db.execute("INSERT INTO MzCalibration VALUES (1,?,?,?,?,?,?,?,?,0,0,0)",
               (g["model_type"], g["timebase"], g["delay"], g["C0"], g["C1"], g["C2"], g["T1_ref"], g["dC1"]))
    return g


def synth_tdf(path, n_frames):
    """A minimal analysis.tdf: the real dataset D MzCalibration row plus a Frames table.

    The calibration constants are the measured dataset D ones from tests/calibration_golden.json, not
    invented numbers -- an invented row would either fail isSupported() or silently define a
    different mass scale. T1 is given a small per-frame spread so the per-frame factor is actually
    exercised rather than collapsing to the reference.
    """
    import sqlite3
    if os.path.exists(path):
        os.remove(path)
    db = sqlite3.connect(path)
    g = mzcal_table(db)
    # Real tdfs reference the calibration row per frame (Frames.MzCalibration); mirror that, since
    # the loader selects the row the frames reference rather than "the only row".
    db.execute("CREATE TABLE Frames (Id INTEGER, T1 REAL, MzCalibration INTEGER)")
    db.executemany("INSERT INTO Frames VALUES (?,?,1)",
                   [(i + 1, g["T1_ref"] + 0.03 * (i % 3)) for i in range(n_frames)])
    # GlobalMetadata: the run-level bounds of the flight-time axis. A real tdf always has them and
    # trace:band_edges=acquisition (the default since 2026-09-10) REFUSES a run without them, so a
    # fixture that omits them cannot exercise the integer detector at all. Values are TEXT, as the
    # vendor writes them; the digitizer sample count is dataset D's.
    db.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    db.executemany("INSERT INTO GlobalMetadata VALUES (?,?)",
                   [("MzAcqRangeLower", "99.990834"), ("MzAcqRangeUpper", "1700.000000"),
                    ("DigitizerNumSamples", "634073")])
    db.commit(); db.close()
    return path


def synth_fake_d(d):
    """A Bruker .d holding only its analysis.tdf: no analysis.tdf_bin, so the loader cannot read a single frame.

    The metadata carries everything the tool decides before the load: dataset D's MzCalibration row and the
    GlobalMetadata bounds (as in synth_tdf), and what dnoise:ms1's setup reads -- dataset D's TimsCalibration row,
    ten Frames (ids 1 and 6 MS1) of 944 scans with integer NumScans/NumPeaks and AccumulationTime 99.958 ms, and 24
    diaPASEF window rows. The shapes are those of the dnoise port's second-review fixtures.
    """
    import sqlite3
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "analysis.tdf")
    if os.path.exists(path):
        os.remove(path)
    db = sqlite3.connect(path)
    db.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    db.executemany("INSERT INTO GlobalMetadata VALUES (?,?)",
                   [("AcquisitionSoftware", "timsTOF"), ("DigitizerNumSamples", "634073"),
                    ("MzAcqRangeLower", "99.990834"), ("MzAcqRangeUpper", "1700.000000"),
                    ("OneOverK0AcqRangeLower", "0.600000"), ("OneOverK0AcqRangeUpper", "1.400000")])
    g = mzcal_table(db)
    db.execute("CREATE TABLE TimsCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER NOT NULL, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9)")
    db.execute("INSERT INTO TimsCalibration VALUES (1, 2, ?,?,?,?,?,?,?,?,?,?)",
               (1, 943, 234.09826168388614, 95.59372590131042, 33.9622641509434, 1, -0.0031464178402676644,
                167.9496150068565, 16.646316600032645, 2241.865411900982))
    db.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, MsMsType INTEGER, NumScans INTEGER, NumPeaks INTEGER,"
               " AccumulationTime REAL, TimsCalibration INTEGER, MzCalibration INTEGER, T1 REAL)")
    db.executemany("INSERT INTO Frames VALUES (?,?,?,944,1000,99.958,1,1,?)",
                   [(i, 0.1 * i, 0 if i % 5 == 1 else 9, g["T1_ref"] + 0.03 * (i % 3)) for i in range(1, 11)])
    db.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER,"
               " IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL)")
    db.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?,?)",
                   [(1 + k % 12, 34 + (k * 7) % 400, 500 + (k * 13) % 400, 400.0 + 0.5 * k, 25.0, 30.0) for k in range(24)])
    db.commit(); db.close()
    return d


def run(binary, inp, out, extra=(), threads=1, expect_fail=False, env=None):
    cmd = [binary, "-in", inp, "-out", out, "-threads", str(threads),
           "-assembly:min_fragments", "2", "-assembly:require_isotope_support", "true",
           "-trace:ms2_min_length_sec", "0", "-trace:ms1_split_valleys", "0",
           "-trace:ms2_split_valleys", "0"] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if expect_fail:
        assert r.returncode != 0, "expected a non-zero exit, got success:\n" + r.stdout[-800:]
        return r
    assert r.returncode == 0, "run failed:\n" + (r.stdout + r.stderr)[-2000:]
    return r


def precursors(mzml):
    """(mz, charge) of every emitted spectrum."""
    out, mz = [], None
    for line in open(mzml, errors="replace"):
        if "selected ion m/z" in line: mz = float(re.search(r'value="([^"]+)"', line).group(1))
        elif "charge state" in line and mz is not None:
            out.append((mz, int(re.search(r'value="([^"]+)"', line).group(1)))); mz = None
    return out


def peaks_of(mzml, prec_mz, tol=0.02):
    """m/z list of the first spectrum whose selected ion is within tol of prec_mz."""
    want, arr, prec = False, None, 8
    for line in open(mzml, errors="replace"):
        if "<spectrum " in line: want = False
        elif "selected ion m/z" in line:
            v = float(re.search(r'value="([^"]+)"', line).group(1)); want = abs(v - prec_mz) < tol
        elif want and 'name="m/z array"' in line: arr = "mz"
        elif want and 'name="intensity array"' in line: arr = "it"
        elif want and 'name="64-bit float"' in line: prec = 8
        elif want and 'name="32-bit float"' in line: prec = 4
        elif want and arr == "mz" and "<binary>" in line:
            b = line[line.find("<binary>") + 8: line.find("</binary>")]
            raw = base64.b64decode(b)
            return list(struct.unpack(f"<{len(raw)//prec}{'d' if prec == 8 else 'f'}", raw))
    return []


def ms1_traces_near(tsv, mz, tol=0.02):
    """(n_xic, ...) of every MS1 trace in a diag:dump_ms1_tsv dump near an m/z."""
    out = []
    with open(tsv) as f:
        hdr = next(f).rstrip("\n").split("\t")
        for line in f:
            q = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if abs(float(q["mz"]) - mz) < tol: out.append(int(q["n_xic"]))
    return out


def digest(mzml):
    """Hash the spectrum data only -- the header carries a wall-clock stamp and parameters."""
    h, on, done = hashlib.sha256(), False, False
    for line in open(mzml, "rb"):
        if b"<spectrumList" in line:
            on = True
            line = re.sub(rb'count="0+(\d)', rb'count="\1', line, count=1)   # a tiled run's zero-padded count
        if on: h.update(line)
        if b"</spectrumList>" in line: done = True; break
    assert on and done, f"{mzml}: no complete <spectrumList> -- the empty hash would compare equal to another empty hash"
    return h.hexdigest()


def stamp(path, key):
    """An output's spx:<key> header stamp, or None."""
    m = re.search(r'spx:%s"[^>]*value="([^"]*)"' % key, open(path, encoding="utf-8", errors="replace").read(200000))
    return m.group(1) if m else None


def main():
    if len(sys.argv) < 2: sys.exit(__doc__)
    binary = sys.argv[1]
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="diaspextract_test_")
    os.makedirs(work, exist_ok=True)
    inp = os.path.join(work, "synth.mzML")
    z2_mono, z1_mono = synth(inp)
    print(f"synthetic input: {inp}  (z=2 mono {z2_mono}, z=1 mono {z1_mono})")
    fails, ran = [], []

    def check(name, fn):
        ran.append(name)
        try:
            fn(); print(f"  PASS  {name}")
        except AssertionError as e:
            print(f"  FAIL  {name}: {e}"); fails.append(name)

    def tdf_env():   # the synthetic analysis.tdf, handed to the tool through the sidecar variable
        return dict(os.environ, DIASPEXTRACT_MZPEAK_TDF=synth_tdf(os.path.join(work, "analysis.tdf"), n_frames=64))

    # 1 + 2: the default drops the z=1 envelope and keeps one precursor per envelope
    out_def = os.path.join(work, "default.mzML")
    run(binary, inp, out_def)
    pd = precursors(out_def)

    def c1():
        z2 = [p for p in pd if p[1] == 2]
        assert z2, f"no z=2 precursor emitted; got {pd}"
        assert all(abs(m - z2_mono) < 0.02 for m, _ in z2), \
            f"a z=2 precursor is not at the monoisotope {z2_mono}: {z2} -- an M+1/M+2 peak reseeded"
    check("isotope ownership: one z=2 precursor, at the mono", c1)

    # [2026-09-08] The default is min_charge 1: singly charged calls are emitted, because on a 2-hour
    # tryptic acquisition they are 14-20% of the Sage peptides and physically 1+ by ion mobility.
    # The synthetic fixture is far too small for the mobility veto to fit a trend line, so the veto
    # is inert here and the z=1 envelope must come through untouched.
    def c2():
        assert [p for p in pd if p[1] == 1], \
            f"charge:min_charge=1 is the default but no z=1 precursor was emitted: {pd}"
    check("min_charge default 1 emits the z=1 envelope", c2)

    out_z2 = os.path.join(work, "z2only.mzML")
    run(binary, inp, out_z2, extra=["-charge:min_charge", "2"])
    p2 = precursors(out_z2)

    def c3():
        assert not [p for p in p2 if p[1] == 1], f"charge:min_charge=2 must drop the z=1 envelope: {p2}"
        assert len([p for p in p2 if p[1] == 2]) == len([p for p in pd if p[1] == 2]), \
            "raising min_charge changed the z=2 result"
    check("min_charge 2 drops it and leaves z=2 alone", c3)

    # 3: a negative count must be rejected, not wrapped through Size
    def c4():
        run(binary, inp, os.path.join(work, "bad.mzML"),
            extra=["-assembly:min_fragments", "-1"], expect_fail=True)
    check("a negative min_fragments is rejected", c4)

    def c5():
        run(binary, inp, os.path.join(work, "bad2.mzML"),
            extra=["-gate:delta_rt", "-5"], expect_fail=True)
    check("a negative gate:delta_rt is rejected", c5)

    # 4: determinism across thread counts for a fixed binary
    out_t4 = os.path.join(work, "t4.mzML")
    run(binary, inp, out_t4, threads=4)

    def c6():
        a, b = digest(out_def), digest(out_t4)
        assert a == b, f"spectrum data differs between 1 and 4 threads:\n  {a}\n  {b}"
    check("byte-identical at 1 vs 4 threads", c6)

    # 5: trace extension across gaps -- the check every broken integer detector would have failed
    # 20 cycles, so that after a seven-cycle gap BOTH halves are long enough to pass the MS1
    # minimum trace length (a 3-frame remnant is correctly dropped by OpenMS, which is not the
    # behaviour under test here).
    for label, gaps, want_traces in (("one missing cycle -> ONE trace", (10,), 1),
                                     ("seven missing cycles -> TWO traces", tuple(range(7, 14)), 2)):
        inp_g = os.path.join(work, f"gap{len(gaps)}.mzML")
        synth(inp_g, n_cycles=20, gap_cycles=gaps)
        out_g = os.path.join(work, f"gap{len(gaps)}.out.mzML")
        dump = os.path.join(work, f"gap{len(gaps)}")
        run(binary, inp_g, out_g, extra=["-diag:dump_ms1_tsv", dump])

        def cg(label=label, want=want_traces, dump=dump, out_g=out_g, gaps=gaps):
            tr = ms1_traces_near(dump + ".traces.tsv", z2_mono)
            assert len(tr) == want, f"expected {want} MS1 trace(s) at {z2_mono} with gaps {gaps}, got {len(tr)} (n_xic {tr})"
            if want == 1:
                assert tr[0] >= 20 - len(gaps) - 1, f"the surviving trace is too short: n_xic {tr[0]}"
                # the gapped FRAGMENT must still reach the z=2 spectrum: a rule that ends a trace at
                # its first gap drops it, which is how 92% of the peptides went missing
                pk = peaks_of(out_g, z2_mono)
                assert any(abs(m - 233.11) < 0.02 for m in pk), f"gapped fragment 233.11 missing from the z=2 spectrum: {pk}"
        check(f"trace extension: {label}", cg)

    # 9 + 10: the SHIPPED detector, and the calibration failing closed.
    # Until 2026-09-04 every check above ran on spx:detector=openms -- the fallback -- because the
    # synthetic nativeIDs carried no "frame=" key and there was no tdf, so the integer detector
    # (the default since that morning) was never once exercised by the suite.
    def c9():
        out_i = os.path.join(work, "integer.mzML")
        env = tdf_env()
        run(binary, inp, out_i, env=env)
        got = stamp(out_i, "detector")
        assert got == "integer", (
            f"the SHIPPED default detector did not run: spx:detector={got}. With a parseable "
            f"frame= nativeID and a valid tdf there is nothing left to fall back for.")
    check("the shipped detector (integer) actually runs", c9)

    def c10():
        # No tdf: the integer detector cannot calibrate, and must SAY so rather than invent one.
        out_f = os.path.join(work, "fallback.mzML")
        r = run(binary, inp, out_f, extra=("-trace:detector", "integer"))
        assert stamp(out_f, "detector") == "openms", "no calibration available, yet integer still ran"
        assert "alling back" in (r.stdout + r.stderr), \
            "fell back to the OpenMS detector without warning -- a silent detector switch"
    check("no calibration -> falls back to openms, loudly", c10)

    # 11: one isolation m/z acquired in two window groups (ion-mobility slices) is two windows.
    def c11():
        n_win = lambda r: int(re.search(r"Split into MS1 \+ (\d+) windows", r.stdout).group(1))
        inp_s = os.path.join(work, "shared.mzML")
        synth(inp_s, groups=(1, 9))
        r1 = run(binary, inp, os.path.join(work, "one_slice.mzML"))
        r2 = run(binary, inp_s, os.path.join(work, "two_slices.mzML"))
        assert n_win(r1) == 1, f"the plain input should be ONE window, got {n_win(r1)}"
        assert n_win(r2) == 2, (
            f"two ion-mobility slices of one m/z window were routed to {n_win(r2)} window(s); keyed by "
            f"m/z alone they collapse into one whose frames are not in retention-time order")
    check("one m/z window in two IM slices -> two windows", c11)

    # 12: the MS1 arena compaction walks kept spans in OFFSET order. The vector is sorted by m/z after
    #     its spans were appended in detection order, and the previous container-order walk overwrote
    #     spans it had not copied yet (~1% of precursor XICs, deterministic). The selftest's fixture
    #     presents traces in an order that defeats that walk.
    def c12():
        r = run(binary, inp, os.path.join(work, "selftest_arena.mzML"), extra=("-diag:selftest_arena",))
        assert "[arena] selftest passed" in (r.stdout + r.stderr), \
            "arena selftest did not report success:\n" + (r.stdout + r.stderr)[-1500:]
    check("MS1 arena compaction survives container order != offset order", c12)

    # 12b: the twin fold (assembly:twin_im_tolerance). The default fixture has one z=2 and one z=1 precursor, so no
    #      twin: the fold must not touch it. With groups=(1, 9) the same m/z window is acquired in two window groups,
    #      so each precursor is assembled twice -- twins on RT, m/z, charge and 1/K0 -- and the fold must fire; the
    #      logged count must be exactly the spectra removed, and the header must name the setting either way.
    def c12b():
        r = run(binary, inp, os.path.join(work, "selftest_twins.mzML"), extra=("-diag:selftest_twins",))
        assert "[twins] selftest passed" in (r.stdout + r.stderr), \
            "twin selftest did not report success:\n" + (r.stdout + r.stderr)[-1500:]
        n = lambda f: int(re.search(rb'<spectrumList count="0*(\d+)"', open(f, "rb").read()).group(1))
        folded = lambda r: int(re.search(r"\[twins\] (\d+) spectra folded into a twin", r.stdout).group(1))
        stamp = lambda f: re.search(rb'name="spx:twin_merge"[^>]*value="([^"]*)"', open(f, "rb").read(40000))
        on, off = os.path.join(work, "twins_on.mzML"), os.path.join(work, "twins_off.mzML")
        r_on = run(binary, inp, on)
        run(binary, inp, off, extra=("-assembly:twin_im_tolerance", "-1"))
        assert folded(r_on) == 0 and digest(on) == digest(off), "the fold changed a fixture that has no twin"
        assert stamp(on) and stamp(on).group(1).startswith(b"same frame") and stamp(off) and stamp(off).group(1) == b"off", \
            "the header does not record the twin setting"
        inp_g = os.path.join(work, "twins_groups.mzML")
        synth(inp_g, groups=(1, 9))
        g_on, g_off = os.path.join(work, "twins_g_on.mzML"), os.path.join(work, "twins_g_off.mzML")
        rg = run(binary, inp_g, g_on)
        run(binary, inp_g, g_off, extra=("-assembly:twin_im_tolerance", "-1"))
        assert folded(rg) > 0, "the fold did not fire on a precursor assembled in two window groups"
        assert n(g_off) - n(g_on) == folded(rg), f"off {n(g_off)} - on {n(g_on)} != the logged fold count {folded(rg)}"
    check("twin fold: inert without twins, fires on a precursor assembled twice, count reconciles", c12b)


    # 13: the concurrency cap is enforced. It used to be computed, logged and never checked: the
    #     only bound was a byte budget that over-booked every window, so nothing ever throttled.
    def c13():
        hw = lambda r: int(re.search(r"admission high-water (\d+) of (\d+)", r.stdout).group(1))
        cap = lambda r: int(re.search(r"admission high-water (\d+) of (\d+)", r.stdout).group(2))
        r1 = run(binary, inp, os.path.join(work, "cap1.mzML"), extra=("-perf:max_concurrent_windows", "1"), threads=4)
        assert cap(r1) == 1 and hw(r1) <= 1, f"cap 1 admitted {hw(r1)} windows at once"
        rN = run(binary, inp, os.path.join(work, "capN.mzML"), threads=4)
        assert digest(os.path.join(work, "cap1.mzML")) == digest(os.path.join(work, "capN.mzML")), \
            "the concurrency cap changed the output; results must be written to index-addressed slots"
    check("the concurrency cap bounds admitted windows and does not change output", c13)

    # 14: the writer's invariants. The parallel writer shipped an EMPTY <indexList> for a release
    #     (count="0", one dummy entry) and nothing here looked; the streaming writer (TileWriter) now
    #     owns the count, the ids and the index, so all three are checked on the default output:
    #     declared count == spectra, id="spectrum=K" and index="K" are the rank from 0, and every
    #     index offset lands on its own <spectrum.
    def invariants(path, label):
        data = open(path, "rb").read()
        n_decl = int(re.search(rb'<spectrumList count="(\d+)"', data).group(1))
        specs = [(m.start(), m.group(1), m.group(2)) for m in re.finditer(rb'<spectrum id="([^"]*)" index="(\d+)"', data)]
        assert n_decl == len(specs) > 0, f"{label}: spectrumList count={n_decl} but {len(specs)} <spectrum> elements"
        for k, (_, sid, idx) in enumerate(specs):
            assert int(idx) == k and sid == b"spectrum=%d" % k, f"{label} spectrum {k}: id={sid!r} index={idx!r}"
        offs = re.findall(rb'<offset idRef="([^"]*)">(\d+)</offset>', data)
        assert len(offs) == len(specs), f"{label}: indexList has {len(offs)} offsets for {len(specs)} spectra"
        for k, (ref, off) in enumerate(offs):
            assert ref == b"spectrum=%d" % k and int(off) == specs[k][0], f"{label} index entry {k}: {ref!r} at {off!r}, spectrum at {specs[k][0]}"
        # OpenMS records the offset one byte early (at the newline before the tag); accept that.
        ilo = int(re.search(rb'<indexListOffset>(\d+)</indexListOffset>', data).group(1))
        assert data[ilo:ilo + 16].lstrip().startswith(b"<indexList"), f"{label}: indexListOffset {ilo} does not point at <indexList"
        return len(specs), re.search(rb'<spectrumList count="([^"]*)"', data).group(1)

    def c14():
        invariants(out_def, "default output")
    check("count, ids, index and indexList agree", c14)

    # 15: the cell grid (tile:rt_sec) is thread-invariant: each cell is traced on its own slice with
    #     its own seeds/visited, band stores are absorbed in a fixed order, and the band count is
    #     explicit -- so 1 and 4 threads must give byte-identical spectra at a pitch that cuts the
    #     fixture into two cells (14 cycles of 1.4 s; pitch 10 s = 7 frames -> 2 cells). Both runs keep
    #     the two cells in ONE tile (-tile:cells_per_tile 0): the only 1-thread run of a window body with
    #     two or more cells. At the default of one cell per tile it silently became a tiled run.
    def c15():
        # the cells exist only on the integer detector: without the tdf the tool falls back to
        # openms and clears the grid AFTER logging it (the first version of this check passed vacuously)
        env = tdf_env()
        r1 = run(binary, inp, os.path.join(work, "cell1.mzML"), extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "0"), threads=1, env=env)
        r4 = run(binary, inp, os.path.join(work, "cell4.mzML"), extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "0"), threads=4, env=env)
        m = re.search(r"\[cell\] pitch [0-9.]+ s = (\d+) MS1 frames .*?: (\d+) cells", r1.stdout + r1.stderr)
        assert m and int(m.group(2)) >= 2, "the fixture was not cut into >= 2 cells:\n" + (r1.stdout + r1.stderr)[-800:]
        mt = re.search(r"\[tile\] (\d+) tile\(s\) of up to (\d+) cell", r1.stdout + r1.stderr)
        assert mt and int(mt.group(1)) == 1 and int(mt.group(2)) >= 2, "check 15 must trace >= 2 cells in ONE tile:\n" + (r1.stdout + r1.stderr)[-800:]
        assert stamp(os.path.join(work, "cell1.mzML"), "detector") == "integer", "cells need the integer detector; it did not run"
        assert stamp(os.path.join(work, "cell1.mzML"), "tile_boundaries") not in (None, "", "none"), "the header records no cut lines"
        assert digest(os.path.join(work, "cell1.mzML")) == digest(os.path.join(work, "cell4.mzML")), \
            "tile:rt_sec output differs between 1 and 4 threads"
    check("the cell grid is thread-invariant", c15)

    # 16-18: the TILE loop. A longer fixture (26 cycles of 1.4 s) with the z=1 precursor's apex
    #     moved 6 cycles earlier: pitch 10 s (7 frames) cuts it at 10.8 / 20.6 / 30.4 s, one cell per
    #     tile gives four tiles; z=1 (rt 10.8) is owned by tile 2 and z=2 (rt 19.2) by tile 3, so two
    #     tiles emit and two write nothing. z=2's fragments crest just BEFORE the 20.6-s cut, so its
    #     spectrum in tile 3 depends on fragments CARRIED from tile 2.
    def tiled_fixture(name, **kw):
        path = os.path.join(work, name)
        synth(path, n_cycles=26, z1_shift=-6, **kw)
        return path

    # 15b: the fold keeps thread-count and tile-grouping invariance with real clusters in play (the checks further
    #      down run on fixtures without twins).
    def c15b():
        env = tdf_env()
        inp_g = tiled_fixture("tiles_twins.mzML", groups=(1, 9))
        outs = []
        for name, extra, thr in (("tw_one_t1", ("-tile:cells_per_tile", "0"), 1), ("tw_one_t4", ("-tile:cells_per_tile", "0"), 4),
                                 ("tw_cells_t4", ("-tile:cells_per_tile", "1"), 4)):
            o = os.path.join(work, name + ".mzML")
            r = run(binary, inp_g, o, extra=("-tile:rt_sec", "10") + extra, threads=thr, env=env)
            assert int(re.search(r"\[twins\] (\d+) spectra folded into a twin", r.stdout).group(1)) > 0, f"{name}: the fold did not fire"
            outs.append(o)
        assert digest(outs[0]) == digest(outs[1]) == digest(outs[2]), "the folded output depends on the thread count or the tile grouping"
    check("twin fold: 1 vs 4 threads and one tile vs one cell per tile give the same spectra", c15b)

    # 16: tile-count x thread invariance, two non-empty blocks, the count/index invariants on the
    #     tiled file, the exact (unpadded) count on the one-tile file
    def c16():
        env = tdf_env()
        inp_t = tiled_fixture("tiles.mzML")
        one = os.path.join(work, "tiles_one.mzML")
        # every one-tile reference names 0: the default is one cell per tile since 2026-09-11 (checks 15-18)
        run(binary, inp_t, one, extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "0"), threads=4, env=env)
        n_one, cnt_one = invariants(one, "one-tile file")
        assert cnt_one == str(n_one).encode(), f"one-tile count is not exact: {cnt_one!r}"
        outs = []
        for thr in (1, 4):
            o = os.path.join(work, f"tiles{thr}.mzML")
            r = run(binary, inp_t, o, extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "1"), threads=thr, env=env)
            log = r.stdout + r.stderr
            m = re.search(r"\[tile\] (\d+) tile\(s\) of up to 1 cell", log)
            assert m and int(m.group(1)) >= 3, "the fixture did not run as >= 3 tiles:\n" + log[-800:]
            written = [int(x) for x in re.findall(r"\[tile\] \d+/\d+: (\d+) spectra written", log)]
            assert sum(1 for w in written if w > 0) >= 2, f"fewer than two tiles wrote spectra: {written}"
            assert "coverage:" in log, "no coverage line"
            outs.append(o)
        hdr = open(outs[0], encoding="utf-8", errors="replace").read(200000)
        mt = re.search(r'spx:tiles"[^>]*value="(\d+)"', hdr)
        assert mt and int(mt.group(1)) >= 3, "the header does not record the tile count"
        d1 = digest(one)
        for o in outs:
            assert digest(o) == d1, f"tiled output {os.path.basename(o)} differs from the one-tile run"
        n_t, cnt_t = invariants(outs[1], "tiled file")
        assert n_t == n_one
    check("the tile loop is tile-count- and thread-invariant", c16)

    # 17: the carry is load-bearing: the same tiled run WITHOUT the carry (a test-only switch)
    #     must differ from the one-tile run, and the normal run must report carried fragments
    def c17():
        env = tdf_env()
        inp_t = os.path.join(work, "tiles.mzML"); one = os.path.join(work, "tiles_one.mzML")
        assert os.path.exists(inp_t) and os.path.exists(one), "check 16 must run first"
        r = run(binary, inp_t, os.path.join(work, "tiles_c.mzML"), extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "1"), threads=2, env=env)
        carried = [int(x) for x in re.findall(r"carried in (\d+)", r.stdout + r.stderr)]
        assert carried and max(carried) > 0, "no tile reported carried-in fragments:\n" + (r.stdout + r.stderr)[-800:]
        env_nc = dict(env, DIASPEXTRACT_TILE_NO_CARRY="1")
        nc = os.path.join(work, "tiles_nocarry.mzML")
        run(binary, inp_t, nc, extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "1"), threads=2, env=env_nc)
        assert digest(nc) != digest(one), "switching the carry off changed nothing: the fixture does not exercise it"
        # the TRANSITIVE carry: pitch 5 s (4 frames) is below 2 x gate:delta_rt (6 s), so a
        # fragment must be re-carried across more than one tile to reach its precursor
        one5 = os.path.join(work, "tiles5_one.mzML"); t5 = os.path.join(work, "tiles5_t.mzML")
        run(binary, inp_t, one5, extra=("-tile:rt_sec", "5", "-tile:cells_per_tile", "0"), threads=2, env=env)
        r5 = run(binary, inp_t, t5, extra=("-tile:rt_sec", "5", "-tile:cells_per_tile", "1"), threads=2, env=env)
        m = re.search(r"\[tile\] (\d+) tile\(s\) of up to 1 cell", r5.stdout + r5.stderr)
        assert m and int(m.group(1)) >= 5, "pitch 5 did not give >= 5 tiles"
        assert digest(t5) == digest(one5), "the transitive carry (pitch < 2 delta_rt) does not reproduce the one-tile run"
    check("the boundary carry is load-bearing", c17)

    # 18: a window with NO frames in a middle tile (its MS2 frames in [20.6, 30.4) dropped -- the
    #     range is given with a margin, 1.0 + 14 * 1.4 is 20.599999999999998): the
    #     empty-slab branch runs, the tile scores with carried fragments only, and the tiled run
    #     still equals the one-tile run of the same file
    def c18():
        env = tdf_env()
        inp_e = tiled_fixture("tiles_empty.mzML", drop_ms2=(20.0, 30.0))   # the frames at 20.6 (= 20.599999999999998) .. 29.0
        one = os.path.join(work, "tiles_empty_one.mzML")
        run(binary, inp_e, one, extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "0"), threads=2, env=env)
        t = os.path.join(work, "tiles_empty_t.mzML")
        r = run(binary, inp_e, t, extra=("-tile:rt_sec", "10", "-tile:cells_per_tile", "1"), threads=2, env=env)
        assert "no frames in the tile" in (r.stdout + r.stderr), "the empty-slab branch did not run:\n" + (r.stdout + r.stderr)[-800:]
        assert digest(t) == digest(one), "the tiled run of a window with an empty tile differs from the one-tile run"
        invariants(t, "empty-tile file")
    check("a window with no frames in a tile", c18)

    # 18b: mzPeak output is written in one piece, so a run whose cells make several tiles at the default
    #      tile:cells_per_tile runs them as ONE tile, with a warning; it used to be refused after the load and the
    #      MS1 tracing. Only a build that writes mzPeak can run this; a plain build skips it and says so.
    hlp = subprocess.run([binary, "--help"], capture_output=True, text=True, stdin=subprocess.DEVNULL)
    if re.search(r"valid formats: [^)]*'mzpeak'", re.sub(r"\s+", " ", hlp.stdout + hlp.stderr)):
        def c18b():
            inp_m = tiled_fixture("tiles_mzp.mzML")
            o = os.path.join(work, "tiles.mzpeak")
            r = run(binary, inp_m, o, extra=("-tile:rt_sec", "10"), threads=2, env=tdf_env())
            log = r.stdout + r.stderr
            mt = re.search(r"\[tile\] (\d+) tile\(s\) of up to (\d+) cell", log)
            assert mt and int(mt.group(1)) == 1 and int(mt.group(2)) >= 2, "mzPeak output did not run its cells as one tile:\n" + log[-800:]
            assert "mzPeak output is written in one piece" in log, "the cells ran as one tile without the warning:\n" + log[-800:]
            assert os.path.exists(o) and re.search(r"Wrote [1-9]\d* pseudo-MS2 spectra", log), "no mzPeak spectra written:\n" + log[-800:]
        check("mzPeak output runs its cells as one tile instead of refusing", c18b)
    else:
        print("  SKIP  mzPeak output runs its cells as one tile: this build writes no mzPeak (--help lists no mzpeak format)")

    # 19: the MS1 prune (perf:ms1_prune, default true) drops picked MS1 peaks at or below
    #     trace:noise_threshold_int at load and must leave the output byte-identical. The fixture makes it
    #     load-bearing: an MS1 trace at m/z 400.25 with hits in cycles 6, 10 (apex) and 14 only, while
    #     cycles 7-9 and 11-13 carry nothing but sub-threshold peaks inside its band's m/z range (400.2 at
    #     50 and at exactly 100.0; 400.24-400.38 at 50, of which the chain drops 400.24). A band spectrum
    #     whose peaks are all at or below the threshold is still a SCAN to MassTraceDetection: 3 hits in 9
    #     scans fail min_sample_rate 0.5. Without those spectra it is 3 of 3 and a trace -- so the
    #     no-witness control (DIASPEXTRACT_MS1_PRUNE_NO_WITNESS=1, survivors only) must change the MS1
    #     traces or the band spectra, and the witness chain must change nothing. Each cycle's lowest peak
    #     is its stride-97 m/z sample and places the band edges: 450 in cycles 0-5 / 15-20 (an edge above
    #     the trace, and misses that end its extension), 300 in 6/10/14, 200 in 8. The hit in cycle 14 is
    #     float32 nextafter(100): a survivor. 1 band is exact trivially; 2 and 48 are the real test.
    #     The no-witness control drops the anchors too, and its change is pinned exactly: one more MS1 trace (400.25,
    #     3 of 3) and one more precursor, both [det] digests and the band spectra differ. Like every check here this
    #     runs the RESIDENT load (an mzML is never streamed): the prune inside the streaming consumer (flushMS1_) is
    #     covered only by the cluster gates on real .d input.
    PRUNE_RE = r"\[ms1-prune\] picked (\d+) survivors (\d+) witnesses (\d+) at_noise (\d+) nan (\d+) bad_mz (\d+): kept (\d+)"
    PRUNE_KEYS = {"traces": r"\[det\] MS1 traces n=\S+ digest=\S+", "sorted": r"\[det\] MS1 traces sorted by content: [^\n]*",
                  "prec": r"\[det\] precursors n=\S+ digest=\S+", "bands": r"\[det\] MS1 band spectra \S*", "edges": r"\[ms1-edges\] [^\n]*"}
    nx100 = struct.unpack("<f", struct.pack("<I", struct.unpack("<I", struct.pack("<f", 100.0))[0] + 1))[0]
    def prune_extra(c):
        if c <= 5 or c >= 15: return [(450.0, 150.0, 0.80)]
        if c in (6, 10, 14): return [(300.0, 150.0, 0.80), (400.25, {6: 150.0, 10: 500.0, 14: nx100}[c], 0.95)]
        if c == 8: return [(200.0, 150.0, 0.80), (400.2, 50.0, 0.95)]
        if c in (7, 9): return [(400.2, 50.0, 0.95)]
        return [(400.2, 100.0, 0.95), (400.24, 50.0, 0.95), (400.28, 50.0, 0.95), (400.38, 50.0, 0.95)]
    def prune_arm(inp_p, tag, bands, prune, env_extra=None):
        """One DIASPEXTRACT_DET run: (output, {key: line or None}, every [det] and [ms1-edges] line, log)."""
        out = os.path.join(work, f"{tag}_b{bands}_{prune}.mzML")
        env = dict(os.environ, DIASPEXTRACT_DET="1")
        for v in ("DIASPEXTRACT_MS1_PRUNE_NO_WITNESS", "DIASPEXTRACT_MS1_PRUNE_NO_CHAIN"): env.pop(v, None)
        env.update(env_extra or {})
        r = run(binary, inp_p, out, extra=("-perf:ms1_trace_bands", str(bands), "-perf:ms1_prune", prune, "-diag:dump_ms1_tsv", out), env=env)
        log = r.stdout + r.stderr
        got = {k: (m.group(0) if m else None) for k, m in ((k, re.search(p, log)) for k, p in PRUNE_KEYS.items())}
        assert got["traces"] and got["prec"], f"{tag} bands {bands} prune {prune}: no [det] MS1/precursor line under DIASPEXTRACT_DET:\n" + log[-800:]
        assert (got["bands"] is not None) == (bands > 1) and (got["edges"] is not None) == (bands > 1), \
            f"{tag} bands {bands} prune {prune}: [det] MS1 band spectra / [ms1-edges] expected iff bands > 1: {got}"
        return out, got, re.findall(r"\[(?:det|ms1-edges)\][^\n]*", log), log
    prune_ref = {}   # (bands, prune) -> (digest, every [det]/[ms1-edges] line) of check 19's fixture, for 19c
    def c19():
        inp_p = os.path.join(work, "prune.mzML")
        synth(inp_p, n_cycles=21, ms1_extra=prune_extra)
        n_of = lambda line: int(re.search(r"n=(\d+)", line).group(1))
        for bands in (1, 2, 48):
            off, g_off, det_off, log_off = prune_arm(inp_p, "prune", bands, "false")
            on, g_on, det_on, log_on = prune_arm(inp_p, "prune", bands, "true")
            assert "[ms1-prune]" not in log_off, "perf:ms1_prune false still pruned"
            m = re.search(PRUNE_RE, log_on)
            assert m, "the prune did not log its counters:\n" + log_on[-800:]
            picked, surv, wit, at_noise, _, bad_mz, kept = (int(x) for x in m.groups())
            assert kept == surv + wit and picked > kept and wit > 0 and at_noise == 3 and bad_mz == 0, \
                f"the fixture does not exercise the prune (picked {picked}, survivors {surv}, witnesses {wit}, at_noise {at_noise}, bad_mz {bad_mz}, kept {kept})"
            assert digest(on) == digest(off), f"bands {bands}: perf:ms1_prune true changed the spectra"
            for k in PRUNE_KEYS:
                assert g_on[k] == g_off[k], f"bands {bands}: {k} differs with the prune on:\n  off {g_off[k]}\n  on  {g_on[k]}"
            assert ms1_traces_near(off + ".traces.tsv", 400.25) == ms1_traces_near(on + ".traces.tsv", 400.25) == [], \
                f"bands {bands}: with every witness kept the 400.25 trace must fail min_sample_rate (3 hits in 9 scans)"
            prune_ref[(bands, "false")], prune_ref[(bands, "true")] = (digest(off), det_off), (digest(on), det_on)
            if bands > 1:
                nw, g_nw, _, log_nw = prune_arm(inp_p, "prune_nowit", bands, "true", {"DIASPEXTRACT_MS1_PRUNE_NO_WITNESS": "1"})
                assert "DIASPEXTRACT_MS1_PRUNE_NO_WITNESS: survivors only" in log_nw, "the no-witness control did not announce itself"
                near_nw = ms1_traces_near(nw + ".traces.tsv", 400.25)
                assert n_of(g_nw["traces"]) == n_of(g_off["traces"]) + 1 and n_of(g_nw["prec"]) == n_of(g_off["prec"]) + 1 and near_nw == [3], \
                    f"bands {bands}: dropping the witnesses must add exactly the 400.25 trace (3 of 3) and one precursor:\n" \
                    f"  off   {g_off['traces']} / {g_off['prec']}\n  nowit {g_nw['traces']} / {g_nw['prec']}, near 400.25 {near_nw}"
                assert g_nw["traces"] != g_off["traces"] and g_nw["prec"] != g_off["prec"] and g_nw["bands"] != g_off["bands"], \
                    f"bands {bands}: the no-witness control must change both [det] digests and the band spectra"
    check("the MS1 prune is output-identical, and its witnesses are load-bearing", c19)

    # 19b: the chain's INTERIOR links are load-bearing on their own (check 19 cannot show it: its anchors and gap
    #      endpoints keep every band alive). Every frame has < 97 peaks, so its stride-97 sample is its lowest peak:
    #      399.9 in cycles 0-5, 200.0 in 6-14, 400.6 in 15-20. At 48 bands the trace at 400.25 (hits in cycles 6, 10 =
    #      apex, 14) then owns the core band [399.9, 400.6), range [399.78003, 400.72018). Cycles 7-9 / 11-13 hold,
    #      inside that range, nothing but the links of a sub-threshold chain 399.70..400.80 at 0.05 Th (hop at 400 =
    #      0.12 Th): its first peak 200.0, its gap far side 399.70 and its last peak (612.61) all lie OUTSIDE the range.
    #      With every link kept those six frames stay SCANS (3 hits in 9: no trace); with only the links dropped
    #      (DIASPEXTRACT_MS1_PRUNE_NO_CHAIN=1, anchors kept) they vanish from 26 bands and the trace is 3 of 3. At 1 and 2
    #      bands the outer ranges always hold a frame's first or last peak, so links can never matter there.
    def c19b():
        chain = [round(399.70 + 0.05 * k, 2) for k in range(23)]
        def extra(c):
            if c <= 5: return [(399.9, 150.0, 0.80)]
            if c >= 15: return [(400.6, 150.0, 0.80)]
            if c in (6, 10, 14): return [(200.0, 50.0, 0.80), (400.25, {6: 150.0, 10: 500.0, 14: 150.0}[c], 0.95)]
            return [(200.0, 50.0, 0.80)] + [(mz, 50.0, 0.95) for mz in chain]
        inp_c = os.path.join(work, "prune_chain.mzML")
        synth(inp_c, n_cycles=21, ms1_extra=extra)
        near = lambda out: ms1_traces_near(out + ".traces.tsv", 400.25)
        n_of = lambda line: int(re.search(r"n=(\d+)", line).group(1))
        for bands in (1, 2, 48):
            off, g_off, _, _ = prune_arm(inp_c, "chain", bands, "false")
            on, g_on, _, log_on = prune_arm(inp_c, "chain", bands, "true")
            assert digest(on) == digest(off), f"bands {bands}: perf:ms1_prune true changed the spectra"
            for k in PRUNE_KEYS:
                assert g_on[k] == g_off[k], f"bands {bands}: {k} differs with the prune on:\n  off {g_off[k]}\n  on  {g_on[k]}"
            assert near(off) == near(on) == [], f"bands {bands}: with every link kept the 400.25 trace must be invalid: {near(off)} {near(on)}"
        # the control, at 48 bands (the loop's last arms)
        nc, g_nc, _, log_nc = prune_arm(inp_c, "chain_nochain", 48, "true", {"DIASPEXTRACT_MS1_PRUNE_NO_CHAIN": "1"})
        m_on, m_nc = re.search(PRUNE_RE, log_on), re.search(PRUNE_RE, log_nc)
        assert m_on and m_nc and "DIASPEXTRACT_MS1_PRUNE_NO_CHAIN: interior chain links dropped" in log_nc, \
            "no [ms1-prune] counters, or the no-chain control did not announce itself:\n" + log_nc[-800:]
        surv_on, wit_on, surv_nc, wit_nc = int(m_on.group(2)), int(m_on.group(3)), int(m_nc.group(2)), int(m_nc.group(3))
        assert surv_nc == surv_on and wit_nc == 9 + 6 and wit_on > wit_nc, \
            f"the control must drop the links ONLY: survivors {surv_on} -> {surv_nc}, witnesses {wit_on} -> {wit_nc} (9 first peaks + 6 gap far sides stay)"
        assert g_nc["edges"] == g_off["edges"], "the control moved the edges; it must change band membership only"
        bands_of = lambda line: [int(x) for x in line.split()[-1].split(",")]
        deficit = [a - b for a, b in zip(bands_of(g_off["bands"]), bands_of(g_nc["bands"]))]
        assert set(deficit) == {0, 6} and deficit.count(6) == 26, f"expected the 6 chain frames to leave exactly 26 bands: {deficit}"
        assert n_of(g_nc["traces"]) == n_of(g_off["traces"]) + 1 and near(nc) == [3], \
            f"without the links the 400.25 trace must become valid (3 of 3): {g_off['traces']} -> {g_nc['traces']}, near {near(nc)}"
    check("the MS1 prune's interior chain links are load-bearing on their own", c19b)

    # 19c: stored order. Check 19's fixture with every MS1 spectrum stored in DESCENDING m/z (the IM array in lockstep)
    #      must give exactly the ascending fixture's result at 2 and 48 bands, prune off and on: the digest and every
    #      [det] and [ms1-edges] line. An as-stored stride-97 sample of it would give other band edges, so this pins
    #      the whole pipeline -- the mzML loader's sort, sortSpectra before MS1 tracing, the prune's own sort -- not
    #      frame()'s sort alone.
    def c19c():
        assert prune_ref, "check 19 must run first"
        inp_d = os.path.join(work, "prune_desc.mzML")
        synth(inp_d, n_cycles=21, ms1_extra=prune_extra, ms1_descending=True)
        first = next(l for l in open(inp_d) if 'name="m/z array"' in l)   # spectrum 0 is an MS1 frame
        raw = base64.b64decode(first[first.find("<binary>") + 8: first.find("</binary>")])
        mzs = list(struct.unpack(f"<{len(raw) // 8}d", raw))
        assert len(mzs) > 1 and mzs == sorted(mzs, reverse=True), f"the fixture's MS1 spectrum is not stored descending: {mzs}"
        for bands in (2, 48):
            for prune in ("false", "true"):
                out, _, det, _ = prune_arm(inp_d, "prune_desc", bands, prune)
                want_digest, want_det = prune_ref[(bands, prune)]
                assert digest(out) == want_digest, f"bands {bands} prune {prune}: the descending fixture's spectra differ from the ascending one's"
                assert det and det == want_det, \
                    f"bands {bands} prune {prune}: [det]/[ms1-edges] lines differ from the ascending fixture:\n  asc  {want_det}\n  desc {det}"
    check("MS1 spectra stored in descending m/z give the ascending result, prune off and on", c19c)

    # 19d: a non-finite MS1 m/z is refused after the pick and BEFORE the prune, prune off and on alike. synth packs m/z as
    #      float64 (19c decodes '<d'), so NaN and inf survive the mzML round trip. The [ms1-prune] line is logged after the
    #      resident prune, so its absence shows the refusal happens in the pick loop. The stream-fallback passthrough runs
    #      only on streamed .d/mzPeak input and stays covered by the cluster gates, as flushMS1_ already is.
    def c19d():
        for tag, bad in (("nan", float("nan")), ("inf", float("inf"))):
            inp_n = os.path.join(work, f"prune_{tag}.mzML")
            synth(inp_n, n_cycles=21, ms1_extra=lambda c: prune_extra(c) + ([(bad, 5000.0, 0.95)] if c == 10 else []))
            for bands in (1, 48):
                for prune in ("false", "true"):
                    r = run(binary, inp_n, os.path.join(work, f"{tag}_b{bands}_{prune}.mzML"),
                            extra=("-perf:ms1_trace_bands", str(bands), "-perf:ms1_prune", prune), expect_fail=True)
                    log = r.stdout + r.stderr
                    assert "non-finite m/z" in log, f"{tag} bands {bands} prune {prune}: not refused for its m/z:\n" + log[-800:]
                    assert "[ms1-prune]" not in log, f"{tag} bands {bands} prune {prune}: refused only after the load summary, not at the pick"
    check("a non-finite MS1 m/z is refused before the prune, prune off and on", c19d)

    # 20: dnoise:ms1 (default true) denoises the raw MS1 frames of a Bruker .d, which this suite cannot build. What it
    #     can pin is the rest of the contract: an mzML is not filtered but skipped, logged and stamped, and
    #     -dnoise:ms1 false stamps "off" -- both with the spectra of the default run.
    def c20():
        on, off = os.path.join(work, "dnoise_on.mzML"), os.path.join(work, "dnoise_off.mzML")
        r_on, r_off = run(binary, inp, on), run(binary, inp, off, extra=("-dnoise:ms1", "false"))
        log_on, log_off = r_on.stdout + r_on.stderr, r_off.stdout + r_off.stderr
        assert "[dnoise] MS1 denoising skipped: the input is not a Bruker .d" in log_on, "the default run did not log the skip:\n" + log_on[-800:]
        assert stamp(on, "dnoise_ms1") == "skipped (input is not a Bruker .d)", f"spx:dnoise_ms1 = {stamp(on, 'dnoise_ms1')!r} on an mzML at the default"
        assert "[dnoise]" not in log_off, "-dnoise:ms1 false still logged a [dnoise] line"
        assert stamp(off, "dnoise_ms1") == "off", f"spx:dnoise_ms1 = {stamp(off, 'dnoise_ms1')!r} with -dnoise:ms1 false"
        assert digest(on) == digest(off) == digest(out_def), "the spectra depend on dnoise:ms1 on an input it does not filter"
    check("dnoise:ms1 skips an mzML (logged and stamped), 'false' stamps off, the spectra are unchanged", c20)

    # 21: dnoise:ms1 inverts the loader's m/z through the tdf's MzCalibration table model, but with a Bruker SDK path set
    #     the loader calibrates through the SDK, so the setup refuses before the load. The input is a .d holding only its
    #     analysis.tdf (synth_fake_d: no analysis.tdf_bin, the loader cannot read a frame). With OPENMS_BRUKER_SDK_PATH the
    #     run must be refused naming dnoise:ms1 and the SDK, with no dnoise setup line and no loader output; without it the
    #     same input must get past the dnoise setup and fail in the loader, which shows the refusal comes first.
    def c21():
        fake = synth_fake_d(os.path.join(work, "fake.d"))
        env = dict(os.environ)
        for v in ("OPENMS_BRUKER_SDK_PATH", "DIASPEXTRACT_MZPEAK_TDF", "DIASPEXTRACT_TILE_RESIDENT"): env.pop(v, None)
        setup_line = "MS1 path on the raw frames, before the pick"
        r = run(binary, fake, os.path.join(work, "fake_sdk.mzML"), expect_fail=True, env=dict(env, OPENMS_BRUKER_SDK_PATH="/nonexistent"))
        log = r.stdout + r.stderr
        assert "dnoise:ms1" in log and "Bruker SDK" in log, "OPENMS_BRUKER_SDK_PATH set: not refused naming dnoise:ms1 and the Bruker SDK:\n" + log[-800:]
        assert setup_line not in log and "[stream]" not in log, "OPENMS_BRUKER_SDK_PATH set: refused only after the dnoise setup or the load:\n" + log[-800:]
        r = run(binary, fake, os.path.join(work, "fake_table.mzML"), expect_fail=True, env=env)
        log = r.stdout + r.stderr
        assert setup_line in log, "no Bruker SDK path: the fake .d did not get past the dnoise setup:\n" + log[-1500:]
        assert "[stream] streaming load failed" in log, "no Bruker SDK path: the run did not fail in the loader (no analysis.tdf_bin):\n" + log[-800:]
    check("dnoise:ms1 refuses a Bruker SDK calibration before the load, and gets past its setup without one", c21)

    # 22: the renames' fail-closed guard. Every variable is DIASPEXTRACT_* since the tool became DIAspeXtract; before that
    #     it was SPEXTRACTOR_* (SpeXtractor) and DIASPEXTRACTOR_* (DIAspeXtractor), and an old name would be ignored
    #     without a word. A run with one set under EITHER old prefix -- and with both -- must refuse before anything else,
    #     naming each variable and its DIASPEXTRACT_ name; the same switch under the new name must run and print its
    #     [det] lines.
    def c22():
        env = {k: v for k, v in os.environ.items() if not k.startswith(("SPEXTRACTOR_", "DIASPEXTRACTOR_"))}
        for tag, old_vars in (("spextractor", {"SPEXTRACTOR_DET": "1"}),
                              ("diaspextractor", {"DIASPEXTRACTOR_DET": "1"}),
                              ("both", {"SPEXTRACTOR_LEDGER": "x", "DIASPEXTRACTOR_TILE_RESIDENT": "1"})):
            old = os.path.join(work, f"old_env_{tag}.mzML")
            r = run(binary, inp, old, expect_fail=True, env=dict(env, **old_vars))
            log = r.stdout + r.stderr
            for v in old_vars:
                new_name = "DIASPEXTRACT_" + v.split("_", 1)[1]
                assert re.search(r"(?<![A-Z])" + v + r"\b", log) and re.search(r"(?<![A-Z])" + new_name + r"\b", log), \
                    f"{v}=1: not refused naming the variable and its {new_name} name:\n" + log[-800:]
            assert "[det]" not in log and not os.path.exists(old), f"{tag}: refused only after the run had started:\n" + log[-800:]
        r = run(binary, inp, os.path.join(work, "new_env.mzML"), env=dict(env, DIASPEXTRACT_DET="1"))
        assert "[det] " in r.stdout + r.stderr, "DIASPEXTRACT_DET=1 ran but printed no [det] line:\n" + (r.stdout + r.stderr)[-800:]
    check("a leftover SPEXTRACTOR_* or DIASPEXTRACTOR_* variable refuses the run; its DIASPEXTRACT_ name works", c22)

    print(f"\n{len(ran) - len(fails)}/{len(ran)} checks passed" + (f"; FAILED: {', '.join(fails)}" if fails else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    # A suite that CANNOT RUN must not look like a suite that ran: a failure outside a check() --
    # a fixture that no longer matches the tool's requirements, a binary that refuses an option --
    # used to escape as a traceback, and a gate chain grepping for "checks passed" or "FAIL" saw
    # nothing at all (2026-09-10: the band-edge default started requiring tdf metadata the
    # synthetic fixture did not write, and three chains reported an empty e2e section).
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except BaseException as e:
        import traceback
        traceback.print_exc()
        print(f"\n0/0 checks passed; FAIL the suite could not run: {type(e).__name__}: {e}")
        sys.exit(2)
