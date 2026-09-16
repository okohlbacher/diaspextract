#!/usr/bin/env python3
"""Every `-<group>:<option>` token in the shipped instructions is an option the tool registers.

An option removed from the tool but left inside a shipped example is not a documentation wart: the
tool refuses an unknown option and aborts before it does any work ("Unknown option(s) '[-a:b]' given.
Aborting!"), so the example is dead on arrival. The 1.2.0 audit found exactly that --
bench/run_open.sh.example still passed `-assembly:open_search_safe`, removed in the same release --
after a release commit had edited that very file for something else.

Scope is deliberately the files that tell a reader what to RUN: the bench examples, README.md and
CONTRIBUTING.md. CHANGELOG.md and docs/ are the historical record and name removed options on
purpose; flagging those would train everyone to ignore this check.

    scripts/check_option_tokens.py [--ini ini.xml] [file ...]

With --ini (a `-write_ini` dump from the built binary) the option inventory is the binary's own;
without it, the inventory is parsed from the register*_ calls in src/diaspextract.cpp, so the check
runs in a tree with no build. Exit 1 if any token is not a registered option.
"""
import argparse
import glob
import pathlib
import re
import sys
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_FILES = sorted(glob.glob(str(ROOT / "bench" / "*.example"))) + [
    str(ROOT / "README.md"),
    str(ROOT / "CONTRIBUTING.md"),
]
TOKEN = re.compile(r"-([a-z][a-z_0-9]*):([a-z][a-z_0-9]*)")


def from_ini(path):
    """group:option for every ITEM under a NODE, from a -write_ini dump."""
    known = set()

    def walk(node, prefix):
        for child in node:
            name = child.get("name")
            if child.tag in ("ITEM", "ITEMLIST") and prefix:
                known.add(prefix + name)
            elif child.tag == "NODE":
                # the two outer nodes are the tool name and the instance number, not option groups
                walk(child, name + ":" if prefix is not None else None)

    root = ET.parse(path).getroot()
    for tool in root:                       # <NODE name="DIAspeXtract">
        for inst in tool:                   # <NODE name="1">
            walk(inst, None)
            for child in inst:
                if child.tag == "NODE":
                    walk(child, child.get("name") + ":")
    return known


def from_source(path):
    """group:option for every registered option in the tool's source."""
    text = pathlib.Path(path).read_text(encoding="utf-8")
    return {m for m in re.findall(r'register(?:\w+)?_\(\s*"([a-z][a-z_0-9]*:[a-z][a-z_0-9]*)"', text)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ini", help="a -write_ini dump from the built binary")
    ap.add_argument("files", nargs="*", default=None)
    args = ap.parse_args()

    if args.ini:
        known, where = from_ini(args.ini), args.ini
    else:
        src = ROOT / "src" / "diaspextract.cpp"
        known, where = from_source(src), str(src)
    if not known:
        print(f"no options found in {where} -- the check would pass vacuously", file=sys.stderr)
        return 2

    files = args.files or DEFAULT_FILES
    bad = 0
    for f in files:
        for n, line in enumerate(pathlib.Path(f).read_text(encoding="utf-8").splitlines(), 1):
            for group, opt in TOKEN.findall(line):
                tok = f"{group}:{opt}"
                if tok not in known:
                    print(f"{f}:{n}: -{tok} is not a registered option ({where})")
                    bad += 1
    print(f"checked {len(files)} file(s) against {len(known)} registered options: "
          f"{bad or 'no'} unknown token(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
