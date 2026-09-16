#!/usr/bin/env python3
"""Concatenate RT-tiled DIAspeXtract outputs into one mzML, keeping each tile's CORE precursors.

    tile_concat.py OUT.mzML TILE.mzML:LO:HI [TILE.mzML:LO:HI ...]

A spectrum belongs to the tile whose half-open core [LO, HI) contains its scan start time
(MS:1000016, seconds -- our writer's unit). Tiles must be given in RT order; the output is then in
the canonical order because RT is the primary sort key and tiles are disjoint. The header comes
from the first file; `index=` is renumbered from 0 and `<spectrumList count=` is set to the total;
the indexedmzML wrapper is dropped (its offsets would be wrong, and search engines do not read it).
Two passes (count, then write) so the declared count is exact.
Self-check: `tile_concat.py out whole.mzML:0:1e12` reproduces whole.mzML's semantic digest exactly.
"""
import re, sys

RT = re.compile(rb'accession="MS:1000016"[^>]*value="([^"]+)"')
# OpenMS writes <spectrum id="spectrum=N" index="N" ...>: both carry the rank, both are renumbered
IDX = re.compile(rb'(<spectrum [^>]*?)(index=")(\d+)(")')
SID = re.compile(rb'(<spectrum [^>]*?id="spectrum=)(\d+)(")')
COUNT = re.compile(rb'(<spectrumList count=")(\d+)(")')

def blocks(path, lo, hi, on_head=None, on_tail=None):
    """Yield the <spectrum ...>...</spectrum> blocks of `path` whose scan start time is in [lo, hi)."""
    with open(path, "rb") as f:
        block = None
        for line in f:
            s = line.lstrip()
            if block is None:
                if s.startswith(b"<spectrum "): block = [line]
                elif s.startswith(b"</spectrumList>"):
                    if on_tail: on_tail([line] + [l for l in f][:2])   # </spectrumList>, </run>, </mzML>
                    return
                elif on_head: on_head(line)
                continue
            block.append(line)
            if s.startswith(b"</spectrum>"):
                data = b"".join(block); block = None
                m = RT.search(data)
                if not m: raise SystemExit(f"{path}: spectrum without MS:1000016")
                if lo <= float(m.group(1)) < hi: yield data

def main(out_path, specs):
    tiles = []
    for sp in specs:
        path, lo, hi = sp.rsplit(":", 2)
        tiles.append((path, float(lo), float(hi)))
    # pass 1: count, so the declared count is exact (no placeholder; the self-check is byte-exact)
    counts = [sum(1 for _ in blocks(path, lo, hi)) for path, lo, hi in tiles]
    total = sum(counts)
    head, tail = [], []
    n_out = 0
    with open(out_path, "wb") as out:
        for ti, (path, lo, hi) in enumerate(tiles):
            def on_head(line, ti=ti):
                if ti or line.lstrip().startswith(b"<indexedmzML"): return
                out.write(COUNT.sub(lambda m: m.group(1) + str(total).encode() + m.group(3), line) if b"<spectrumList" in line else line)
            def on_tail(lines): tail[:] = lines
            for data in blocks(path, lo, hi, on_head, on_tail):
                data = IDX.sub(lambda mm: mm.group(1) + mm.group(2) + str(n_out).encode() + mm.group(4), data, count=1)
                data = SID.sub(lambda mm: mm.group(1) + str(n_out).encode() + mm.group(3), data, count=1)
                out.write(data)
                n_out += 1
            print(f"tile {ti} {path} [{lo},{hi}): kept {counts[ti]}", file=sys.stderr)
        out.write(b"".join(tail))
    assert n_out == total, (n_out, total)
    print(f"wrote {n_out} spectra to {out_path}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv) < 3: raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2:])
