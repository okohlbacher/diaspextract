#!/usr/bin/env python3
"""Delete listed spectra from a pseudo-MS/MS mzML; everything else byte-identical.

    mzml_drop.py SRC.mzML DST.mzML DROP.txt        (DROP.txt: one original spectrum index per line)

The exact form of "what if the extractor had not emitted these precursors": extraction is per-precursor, so removing a
precursor upstream and deleting its spectrum downstream give the same file. `id="spectrum=N"` is KEPT, so a search
engine's scan reference still names the original spectrum and results join back to the unfiltered run; `index=` is
renumbered, `<spectrumList count=` is set, and the indexedmzML wrapper is dropped (its offsets would be wrong, and the
search engines do not read it). Self-check: an empty DROP reproduces SRC's semantic digest.
"""
import re, sys

IDX = re.compile(rb'(<spectrum [^>]*?index=")(\d+)(")')
COUNT = re.compile(rb'(<spectrumList count=")(\d+)(")')

def main(src, dst, drop_path):
    drop = {int(x) for x in open(drop_path).read().split()}
    total = 0
    with open(src, "rb") as f:                       # pass 1: the declared count must be exact
        for line in f:
            if line.lstrip().startswith(b"<spectrum "):
                total += int(IDX.search(line).group(2)) not in drop
    n_out = n_drop = 0
    with open(src, "rb") as f, open(dst, "wb") as out:
        keep, in_spec = True, False
        for line in f:
            s = line.lstrip()
            if s.startswith(b"<spectrum "):
                in_spec = True
                keep = int(IDX.search(line).group(2)) not in drop
                if keep:
                    line = IDX.sub(lambda m: m.group(1) + str(n_out).encode() + m.group(3), line, count=1); n_out += 1
                else: n_drop += 1
            elif not in_spec:
                if s.startswith(b"<indexedmzML"): continue
                if s.startswith(b"<indexList") : break       # past </mzML>: the index wrapper's tail
                if b"<spectrumList" in line: line = COUNT.sub(lambda m: m.group(1) + str(total).encode() + m.group(3), line)
            if keep or not in_spec: out.write(line)
            if s.startswith(b"</spectrum>"): in_spec, keep = False, True
    assert n_out == total and n_drop == len(drop), (n_out, total, n_drop, len(drop))
    print(f"kept {n_out:,}, dropped {n_drop:,} -> {dst}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) != 4: raise SystemExit(__doc__)
    main(*sys.argv[1:])
