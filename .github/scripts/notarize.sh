# Sourced by build.yml (macOS legs of "build" and "desktop"). Needs MACOS_APPLE_ID, MACOS_TEAM_ID, MACOS_NOTARY_PASSWORD.

# notarytool takes the app-specific password only as an argument: it is in that process's argv on this throwaway runner
# for the duration of the call, and masked in the log.
notarize() {
  local out id status
  out=$(xcrun notarytool submit "$1" --apple-id "$MACOS_APPLE_ID" --team-id "$MACOS_TEAM_ID" \
          --password "$MACOS_NOTARY_PASSWORD" --wait --timeout 80m --output-format json) || true
  id=$(printf '%s' "$out" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("id",""))' 2>/dev/null || true)
  status=$(printf '%s' "$out" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("status",""))' 2>/dev/null || true)
  echo "notarization of $(basename "$1"): ${status:-no answer} (submission ${id:-none})"
  if [ "$status" != Accepted ]; then
    # a rejection names only the submission; the log names the offending file
    [ -z "$id" ] || xcrun notarytool log "$id" --apple-id "$MACOS_APPLE_ID" --team-id "$MACOS_TEAM_ID" --password "$MACOS_NOTARY_PASSWORD"
    echo "::error::notarization of $(basename "$1") was not accepted"; return 1
  fi
}

# stapler reaches Apple's ticket server (CloudKit), which the macos-15-intel runner loses again and again (Error 68 /
# -1005, three release runs in a row). The item is notarized at this point; the staple only lets it open offline, as
# Gatekeeper fetches the ticket online otherwise. So retry, then ship it unstapled with a warning.
staple() {
  local i
  for i in 1 2 3 4; do
    xcrun stapler staple "$1" && xcrun stapler validate "$1" && return 0
    [ "$i" = 4 ] || sleep $((i * 30))
  done
  echo "::warning title=not stapled::$1 is notarized, but Apple's ticket server was unreachable; Gatekeeper checks it online on first open"
  return 1
}
