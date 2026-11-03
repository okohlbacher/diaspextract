#!/bin/bash
# Copy the non-system library closure of <bundle>/bin/diaspextract into <bundle>/lib and point every reference at
# @executable_path/../lib/, so the bundle runs from wherever it is unpacked. Used by .github/workflows/build.yml.
#
#   .github/scripts/bundle-macos.sh <bundle> <search dir>...
#
# dylibbundler does the dylibs. It does NOT follow frameworks (it leaves "@rpath/QtCore.framework/..." in place), and
# Homebrew's Qt is frameworks, so they are handled here: each framework's binary is copied flat into lib/, the
# references are rewritten, and dylibbundler is run again on the copies for their own dylib dependencies.
# Every Mach-O is then left with one LC_RPATH (dylibbundler rewrites each old rpath to the same value, and dyld refuses
# duplicates) and an ad-hoc signature (install_name_tool invalidates signatures; arm64 runs nothing unsigned). The
# Developer ID signature, when there is one, replaces it later. Finally every reference is checked.
set -euo pipefail
D=${1:?usage: bundle-macos.sh <bundle> <search dir>...}; shift
EXE="$D/bin/diaspextract"; L="$D/lib"; P='@executable_path/../lib/'
SEARCH=("$@")
S=(); for s in "${SEARCH[@]}"; do S+=(-s "$s"); done

# </dev/null: dylibbundler prompts for a library it cannot find; a prompt must fail the job, not hang it.
dylibbundler -od -b -x "$EXE" -d "$L" -p "$P" "${S[@]}" < /dev/null

machos() { find "$D/bin" "$L" -type f | while IFS= read -r f; do if file -b "$f" | grep -q Mach-O; then echo "$f"; fi; done; }
# dependency lines are tab-indented; the others are headers (one per architecture in a universal binary)
deps() { otool -L "$1" | awk '/^\t/{print $1}' | sort -u; }
resolve() {   # the file a framework reference names, through the search dirs for @rpath and friends
  local ref=$1 rel s
  case "$ref" in
    @*) rel=${ref#*/}
        for s in "${SEARCH[@]}"; do [ -f "$s/$rel" ] && { echo "$s/$rel"; return 0; }; done
        return 1 ;;
    *)  [ -f "$ref" ] && echo "$ref" ;;
  esac
}

while :; do
  new=()
  while IFS= read -r f; do
    while IFS= read -r ref; do
      case "$ref" in *.framework/*) ;; *) continue ;; esac
      case "$ref" in /System/*) continue ;; esac
      name=$(basename "$ref")
      if [ ! -f "$L/$name" ]; then
        src=$(resolve "$ref") || { echo "::error::cannot find framework $ref (needed by $f)" >&2; exit 1; }
        cp -L "$src" "$L/$name"; chmod u+w "$L/$name"
        install_name_tool -id "$P$name" "$L/$name" 2> /dev/null
        new+=("$L/$name")
        echo "framework: $ref -> lib/$name"
      fi
      install_name_tool -change "$ref" "$P$name" "$f" 2> /dev/null
    done < <(deps "$f")
  done < <(machos)
  [ ${#new[@]} -gt 0 ] || break
  X=(); for n in "${new[@]}"; do X+=(-x "$n"); done
  dylibbundler -of -b "${X[@]}" -d "$L" -p "$P" "${S[@]}" < /dev/null
done

while IFS= read -r f; do
  while [ "$(otool -l "$f" | awk '/LC_RPATH/{g=1;next} g&&/^ *path /{print $2; g=0}' | grep -cx "$P")" -gt 1 ]; do
    install_name_tool -delete_rpath "$P" "$f"
  done
  codesign --force --sign - "$f" 2> /dev/null
done < <(machos)

# Every dependency is either the system's or a file in lib/ referenced through $P.
bad=0
while IFS= read -r f; do
  self=$(otool -D "$f" | grep -v ':$' | sort -u || true)
  while IFS= read -r ref; do
    [ "$ref" = "$self" ] && continue
    case "$ref" in
      /usr/lib/*|/System/*) ;;
      "$P"*) [ -f "$L/${ref#"$P"}" ] || { echo "::error::$f needs $ref, which is not in lib/" >&2; bad=1; } ;;
      *) echo "::error::$f still references $ref" >&2; bad=1 ;;
    esac
  done < <(deps "$f")
  while IFS= read -r rp; do
    [ "$rp" = "$P" ] || { echo "::error::$f has LC_RPATH $rp" >&2; bad=1; }
  done < <(otool -l "$f" | awk '/LC_RPATH/{g=1;next} g&&/^ *path /{print $2; g=0}')
done < <(machos)
exit $bad
