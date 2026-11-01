#!/usr/bin/env python3
"""Add value="" to every valueless <cvParam> so MSFragger's batmass-io can read an OpenMS mzML.

    fix_mzml_cvparams.py IN.mzML OUT.mzML

OpenMS writes boolean-type cvParams with NO value attribute (legal mzML 1.1: e.g.
<cvParam ... name="centroid spectrum"/>). batmass-io (MSFragger's reader) calls String.trim() on
the missing value, NPEs on the spectrum, and silently reports "Scans = 0" -- the file looks empty
and MSFragger searches nothing. This is the reader's bug, not the file's, but the fix lives here:
give every cvParam a value.

Streaming, line-oriented: the 14.8 GB file never lands in memory. Only cvParam TAGS are touched;
the base64 binary arrays (the bulk) pass through untouched. A cvParam that already has value= is
left exactly as-is, so RT/mz/intensity values are never altered. The added bytes move every spectrum, so an
indexedmzML gets its <indexList> and <indexListOffset> rebuilt from the offsets written here (an indexed reader
seeks by them), and its <fileChecksum> is recomputed: the SHA-1 of the output up to and including "<fileChecksum>".
"""
import hashlib, re, sys

SRC, DST = sys.argv[1], sys.argv[2]
TAG = re.compile(rb"<cvParam\b[^>]*?/?>")                              # a whole <cvParam ...> or <cvParam .../> tag
ELEM = re.compile(rb'<(spectrum|chromatogram) [^>]*?\bid="([^"]*)"')   # an indexed element and its id


def fix(tag):
    if b"value=" in tag:
        return tag                                  # already has a value -> untouched
    # insert value="" just before the closing (/> or >), preserving self-closing form
    if tag.endswith(b"/>"):
        return tag[:-2].rstrip() + b' value=""/>'
    return tag[:-1].rstrip() + b' value=""' + b'>'


n_lines = n_fixed = 0
offsets = {b"spectrum": [], b"chromatogram": []}   # (id, byte offset of its '<' in DST)
sha1 = hashlib.sha1()


class Hashed:
    """The output file; every byte written also goes into the checksum."""
    def __init__(self, f): self.f = f
    def tell(self): return self.f.tell()
    def write(self, b): sha1.update(b); self.f.write(b)
    def writelines(self, bs):
        for b in bs: self.write(b)


with open(SRC, "rb") as fi, open(DST, "wb") as raw:
    fo = Hashed(raw)
    for line in fi:
        n_lines += 1
        if line.lstrip().startswith(b"<indexList"):  # the old index points into SRC: write the new one
            at = fo.tell()
            fo.write(b'<indexList count="%d">\n' % sum(1 for v in offsets.values() if v))
            for kind, ids in offsets.items():
                if ids:
                    fo.write(b'\t<index name="%s">\n' % kind)
                    fo.writelines(b'\t\t<offset idRef="%s">%d</offset>\n' % io for io in ids)
                    fo.write(b"\t</index>\n")
            fo.write(b"</indexList>\n")
            for rest in fi:                          # skip the old index; keep what follows it
                if rest.lstrip().startswith(b"<indexListOffset>"):
                    fo.write(b"<indexListOffset>%d</indexListOffset>\n" % at)
                elif b"<fileChecksum>" in rest:      # the checksum covers everything up to and including this tag
                    head = rest[:rest.index(b"<fileChecksum>") + len(b"<fileChecksum>")]
                    fo.write(head)
                    raw.write(sha1.hexdigest().encode() + b"</fileChecksum>\n")
                elif not (rest.lstrip().startswith((b"<offset", b"<index ", b"<indexList")) or
                          rest.strip() in (b"</index>", b"</indexList>")):
                    fo.write(rest)
            break
        if b"<cvParam" in line and b"value=" not in line:
            line, k = TAG.subn(lambda m: fix(m.group(0)), line)
            n_fixed += k
        m = ELEM.search(line)
        if m:
            offsets[m.group(1)].append((m.group(2), fo.tell() + m.start()))
        fo.write(line)
        if n_lines % 20_000_000 == 0:
            sys.stderr.write("  ...%d lines, %d cvParams fixed\n" % (n_lines, n_fixed))
            sys.stderr.flush()
sys.stderr.write("done: %d lines, %d valueless cvParams patched -> %s\n" % (n_lines, n_fixed, DST))
