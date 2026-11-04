# Security policy

## Reporting a vulnerability

Please report security issues privately via GitHub's "Report a vulnerability" button on the
Security tab, rather than opening a public issue. If that button is not available, open an issue
titled "security contact request" with no details in it, and we will arrange a private channel.

## Scope

DIAspeXtract is a batch command-line tool. It reads instrument data and writes mzML or mzPeak; it runs no server.
It makes a network request only when `-in` is an `s3://` URI: it then reads that one object through Arrow's S3
client, with the credentials of the usual AWS environment, into a temporary file, and logs the URI without any
credentials it carries. The realistic security surface is therefore **untrusted input files**: a malformed or hostile
Bruker `.d`, Thermo `.raw`, `.mzML` or `.mzpeak` reaching the parser. A `.raw` file is read by Thermo's
RawFileReader, a .NET library hosted in the tool's process (on macOS the tool carries the `allow-jit` entitlement for
it; library validation stays on). Reports of crashes, out-of-bounds reads or unbounded allocation triggered by an input file
are in scope and welcome.

Before the reader runs, the tool checks a `.d` itself: that it is diaPASEF, that every frame lies inside
`analysis.tdf_bin`, that no frame claims more points (`Frames.NumPeaks`) than its block can decode to (256 bytes per
byte of the block, 8 bytes per point) or than its scans have flight-time bins, and that the calibration tables hold
numbers. mzPeak members are read into buffers that grow with the data actually read, not with the size the archive
claims. Most parsing is still done by OpenMS and by the vendor libraries. If the flaw is there, we will help route it
upstream.

## Known issues

- **Allocation from `NumPeaks`.** The checks above bound it by the frame's block in the file, not by memory: a hostile
  `.d` whose one frame occupies a 134 MB block can still make the reader allocate arrays for 2^32 points (about 96 GB),
  and fail with an out-of-memory error.
- **The opentims frame decoder** (vendored by OpenMS, `opentims++/opentims.cpp`, `decompress_frame`) trusts the per-scan
  peak counts inside a decompressed frame: it writes `scan_ids`, `tofs` and `intensities` for as many points as the
  scan headers add up to and reads the decompressed buffer by those offsets, without checking either against
  `Frames.NumPeaks` or the buffer's length. A crafted frame can therefore write and read out of bounds. Real
  instrument files are consistent; the check belongs upstream (opentims), and we have not patched the vendored copy.
  Do not run DIAspeXtract on `.d` files from untrusted sources outside a sandbox.

## Data handling

DIAspeXtract reads local files and writes local files. It does not phone home, collect telemetry, or
transmit anything; it switches OpenMS's update check off whatever the environment says. If you find otherwise, that is a security bug — report it.
