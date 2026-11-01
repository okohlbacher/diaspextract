#!/usr/bin/env python3
"""Cut a retention-time range out of a Bruker timsTOF .d into a smaller, self-contained .d.

    slice_d.py IN.d OUT.d RT_FROM RT_TO      (seconds, Frames.Time)

Keeps the frames with RT_FROM <= Time < RT_TO, copies their blocks into a compacted analysis.tdf_bin and rewrites
Frames.TimsId to the new offsets; rows of other tables that name a dropped frame are deleted. Frame ids are kept. The
result reads like the original run cut short: a fast real input for the end-to-end checks and for A/B screens.
"""
import os, shutil, sqlite3, struct, sys


def main(src, dst, t0, t1):
    if os.path.exists(dst):
        sys.exit(f"{dst} exists")
    os.makedirs(dst)
    try:
        cut(src, dst, t0, t1)
    except BaseException:          # a refusal (sys.exit) or an error: no partial .d is left behind
        shutil.rmtree(dst, ignore_errors=True)
        raise


def cut(src, dst, t0, t1):
    shutil.copyfile(os.path.join(src, "analysis.tdf"), os.path.join(dst, "analysis.tdf"))
    db = sqlite3.connect(os.path.join(dst, "analysis.tdf"))
    keep = db.execute("SELECT Id, TimsId FROM Frames WHERE Time >= ? AND Time < ? ORDER BY TimsId", (t0, t1)).fetchall()
    if not keep:
        sys.exit("no frame in that range")
    if db.execute("SELECT Value FROM GlobalMetadata WHERE Key = 'TimsCompressionType'").fetchone() != ("2",):
        sys.exit("only TimsCompressionType 2 is supported (the block length is the first 4 bytes)")
    ids = {i for i, _ in keep}
    new_off = {}
    with open(os.path.join(src, "analysis.tdf_bin"), "rb") as fi, open(os.path.join(dst, "analysis.tdf_bin"), "wb") as fo:
        for fid, off in keep:
            fi.seek(off)
            (n,) = struct.unpack("<I", fi.read(4))
            fi.seek(off)
            block = fi.read(n)
            if len(block) != n:
                sys.exit(f"frame {fid}: block cut short in the input")
            new_off[fid] = fo.tell()
            fo.write(block)
    db.execute("CREATE TEMP TABLE keep(Id INTEGER PRIMARY KEY)")
    db.executemany("INSERT INTO keep VALUES (?)", [(i,) for i in ids])
    for (table,) in db.execute("SELECT name FROM sqlite_master WHERE type = 'table'").fetchall():
        cols = [c[1] for c in db.execute(f"PRAGMA table_info({table})")]
        col = "Id" if table == "Frames" else "Frame" if "Frame" in cols else None
        if col:
            db.execute(f"DELETE FROM {table} WHERE {col} NOT IN (SELECT Id FROM keep)")
    db.executemany("UPDATE Frames SET TimsId = ? WHERE Id = ?", [(o, i) for i, o in new_off.items()])
    db.commit()
    db.execute("VACUUM")
    db.close()
    print(f"{dst}: {len(keep)} frames, {os.path.getsize(os.path.join(dst, 'analysis.tdf_bin'))} bytes of frame data")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
