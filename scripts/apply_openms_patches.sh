#!/bin/bash
# Apply DIAspeXtract's OpenMS patches to an OpenMS source tree and install the shared calibration header.
#
# WITHOUT the calibration, MassTrace and TOPP-external patches the tool does not build; without the EPD patch it builds and runs, but the Bruker .d reader falls back to a two-point
# linear-in-sqrt TOF->m/z chord that is -5..-11 ppm biased on the files we measured -- which costs
# ~6-11% of closed-search peptide identifications (docs/BENCHMARK-MATRIX-2026-09-01.md). The build is
# silent about it; only the run log says which calibration was used, and every emitted mzML carries
# the answer in the spx:mz_calibration userParam.
#
# Usage: scripts/apply_openms_patches.sh /path/to/OpenMS
#
# The patches are cut against the exact public OpenMS commit named in patches/openms.lock. Check
# that commit out first; against anything else they may fail, or worse, apply with fuzz.
set -euo pipefail
cd "$(dirname "$0")/.."
TREE=${1:?usage: apply_openms_patches.sh /path/to/OpenMS-source-tree}
[ -f "$TREE/src/openms/source/FORMAT/BrukerTimsFile.cpp" ] || {
  echo "not an OpenMS source tree (no src/openms/source/FORMAT/BrukerTimsFile.cpp): $TREE" >&2; exit 1; }

# the calibration model is header-only and lives in THIS repo; the patch includes it
install -m 0644 src/TdfMzCalibration.h "$TREE/src/openms/include/OpenMS/FORMAT/TdfMzCalibration.h"
echo "installed src/TdfMzCalibration.h -> $TREE/src/openms/include/OpenMS/FORMAT/"

PIN=$(sed -n 's/^OPENMS_BASE=//p' patches/openms.lock)
HEAD_SHA=$(git -C "$TREE" rev-parse HEAD 2>/dev/null || echo unknown)
if [ -n "$PIN" ] && [ "$HEAD_SHA" != unknown ] && [ "$HEAD_SHA" != "$PIN" ]; then
  echo "WARNING: $TREE is at $HEAD_SHA, not the pinned $PIN. The patches may apply with fuzz." >&2
fi

for p in patches/*.patch; do
  if patch -p1 -N --dry-run -d "$TREE" < "$p" >/dev/null 2>&1; then
    patch -p1 -N -d "$TREE" < "$p"; echo "applied  $p"
  elif patch -p1 -R --dry-run -d "$TREE" < "$p" >/dev/null 2>&1; then
    echo "already applied, skipping  $p"
  else
    echo "FAILED to apply $p -- the tree may be a different OpenMS version" >&2; exit 1
  fi
done
echo "done. Rebuild OpenMS (the patches change its headers and libOpenMS), then build diaspextract and verify the run log says:"
echo "  TOF m/z calibration: TDF MzCalibration table model (exact, license-free)"
