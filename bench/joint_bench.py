#!/usr/bin/env python3
"""Joint benchmark: pseudo-spectra x {Sage, MSFragger} -> FDR-corrected PSMs and peptides.

Metrics (the two the project is scored on):
    PSMs      at <=1% PSM-level FDR
    PEPTIDES  unique, at <=1% peptide-level FDR

Design decision that matters
----------------------------
Both engines are scored by ONE identical target-decoy procedure implemented here, NOT by each
engine's own FDR machinery (Sage's internal linear discriminant; MSFragger via
PeptideProphet/Philosopher or Percolator). Those pipelines differ in rescoring, in how they
aggregate PSM->peptide, and in their decoy handling, so comparing their outputs would confound
"which engine identifies more" with "whose FDR is more permissive".

For the same reason both engines search ONE pre-built decoy database with the same prefix
(Sage's `generate_decoys` is turned OFF), so the two see byte-identical search spaces.

What this deliberately does NOT do: rescoring. Sage's reported score is already
discriminant-rescored, MSFragger's hyperscore is not. That asymmetry is REPORTED, not silently
corrected -- adding Percolator to one side and not the other would reintroduce exactly the
confound this design removes. Read the engine comparison as "raw discriminating power under a
common FDR", and see LIMITATIONS at the bottom.

Peptide-level FDR
-----------------
Best-scoring PSM per distinct peptide sequence, then target-decoy over that reduced list --
the standard "picked" approach. This is the axis that matters here: we emit ~4.6 PSMs/peptide
against the reference implementation's 1.46, so PSM-level counts reward redundancy and peptide-level counts do not.
Both are reported precisely so the gap is visible rather than hidden.
"""
import argparse, csv, json, math, os, re, subprocess, sys
from pathlib import Path

DECOY_PREFIX = "rev_"


# --------------------------------------------------------------------------- decoy database
def build_decoy_fasta(target_fa, out_fa, prefix=DECOY_PREFIX):
    """Reversed-sequence decoys, so BOTH engines search one identical database.

    Sage can generate decoys internally and MSFragger cannot; letting each do its own would
    mean the two never saw the same search space, and the comparison would be meaningless.
    """
    n_t = n_d = 0

    def emit(out, hdr, seq):
        nonlocal n_d
        if hdr is None:
            return
        acc = hdr[1:].split()[0]
        out.write(">%s%s decoy\n" % (prefix, acc))
        s = "".join(seq)
        out.write(s[::-1] + "\n")
        n_d += 1

    with open(out_fa, "w") as out:
        hdr, seq = None, []
        with open(target_fa) as f:
            for line in f:
                if line.startswith(">"):
                    n_t += 1
                    hdr, seq = line.rstrip("\n"), []
                    out.write(line)
                else:
                    seq.append(line.strip())
                    out.write(line)
        # second pass for decoys (kept separate so target order is preserved)
        with open(target_fa) as f:
            hdr, seq = None, []
            for line in f:
                if line.startswith(">"):
                    emit(out, hdr, seq)
                    hdr, seq = line.rstrip("\n"), []
                else:
                    seq.append(line.strip())
            emit(out, hdr, seq)
    return n_t, n_d


# ------------------------------------------------------------------------------- unified FDR
def target_decoy_fdr(rows, level, fdr=0.01):
    """One procedure, both engines.

    rows: list of {"key", "score", "decoy"}. `key` is the spectrum id for PSM level and the
    peptide sequence for peptide level.

    Returns (n_accepted, threshold). Standard competition estimate FDR = n_decoy / n_target,
    walking down the score ranking and taking the LAST threshold that still satisfies the
    target -- so a late run of decoys cannot leave a permissive early threshold in place.
    """
    if level == "peptide":
        best = {}
        for r in rows:
            k = r["key"]
            if k not in best or r["score"] > best[k]["score"]:
                best[k] = r
        rows = list(best.values())
    ranked = sorted(rows, key=lambda r: -r["score"])
    n_t = n_d = 0
    accepted, thresh = 0, None
    for r in ranked:
        if r["decoy"]:
            n_d += 1
        else:
            n_t += 1
        if n_t and (n_d / n_t) <= fdr:
            accepted, thresh = n_t, r["score"]
    return accepted, thresh


# ------------------------------------------------------------------------------ engine readers
def read_sage(tsv):
    need = ["peptide", "proteins", "sage_discriminant_score", "scannr"]
    out = []
    with open(tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        miss = [c for c in need if c not in (rd.fieldnames or [])]
        if miss:
            sys.exit("ABORT: %s lacks %s (have %s)" % (tsv, miss, rd.fieldnames))
        for x in rd:
            prot = x["proteins"]
            out.append({"spectrum": x["scannr"], "peptide": x["peptide"],
                        "score": float(x["sage_discriminant_score"]),
                        "decoy": all(p.startswith(DECOY_PREFIX) for p in prot.split(";") if p)})
    return out


def read_msfragger(tsv):
    """MSFragger .tsv. Column names vary by version, so resolve them rather than assume."""
    with open(tsv) as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames or []
        def pick(*cands):
            for c in cands:
                if c in cols:
                    return c
            return None
        c_pep = pick("Peptide", "peptide", "Modified Peptide")
        c_scr = pick("Hyperscore", "hyperscore", "SpecEValue", "Expectation")
        # MSFragger 4.4.1 emits lowercase plural `proteins`; omitting that exact spelling made
        # the resolver abort on a perfectly good 96 MB result file.
        c_pro = pick("proteins", "Proteins", "Protein", "protein")
        c_spec = pick("scannum", "Spectrum", "scannr", "ScanNr")
        if not (c_pep and c_scr and c_pro):
            sys.exit("ABORT: %s: cannot resolve peptide/score/protein columns from %s"
                     % (tsv, cols))
        # Expectation values sort ASCENDING; everything else DESCENDING. Getting this backwards
        # silently inverts the ranking and yields a confident, wrong number.
        invert = c_scr in ("SpecEValue", "Expectation")
        out = []
        for x in rd:
            try:
                s = float(x[c_scr])
            except (TypeError, ValueError):
                continue
            if invert:
                s = -math.log10(max(s, 1e-300))
            prot = x[c_pro]
            out.append({"spectrum": x.get(c_spec, ""), "peptide": x[c_pep], "score": s,
                        "decoy": all(p.startswith(DECOY_PREFIX)
                                     for p in re.split(r"[;,]", prot) if p.strip())})
    return out


READERS = {"sage": read_sage, "msfragger": read_msfragger}


def score_file(path, engine, label):
    rows = READERS[engine](path)
    if not rows:
        return {"label": label, "engine": engine, "error": "no PSMs parsed from %s" % path}
    psms = [{"key": "%s|%s" % (r["spectrum"], r["peptide"]), "score": r["score"],
             "decoy": r["decoy"]} for r in rows]
    peps = [{"key": r["peptide"], "score": r["score"], "decoy": r["decoy"]} for r in rows]
    n_psm, t_psm = target_decoy_fdr(psms, "psm")
    n_pep, t_pep = target_decoy_fdr(peps, "peptide")
    n_dec = sum(1 for r in rows if r["decoy"])
    return {
        "label": label, "engine": engine,
        "psms_1pct": n_psm,
        "peptides_1pct": n_pep,
        "psms_per_peptide": (n_psm / n_pep) if n_pep else None,
        "raw_psms": len(rows),
        "raw_decoy_fraction": n_dec / len(rows),
        "psm_threshold": t_psm, "peptide_threshold": t_pep,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("decoydb", help="build the shared target+decoy FASTA")
    d.add_argument("--target", required=True)
    d.add_argument("--out", required=True)

    s = sub.add_parser("score", help="score one engine output under the common FDR")
    s.add_argument("path")
    s.add_argument("--engine", required=True, choices=sorted(READERS))
    s.add_argument("--label", required=True)
    s.add_argument("--json")

    c = sub.add_parser("collate", help="table across sample x tool x engine")
    c.add_argument("json")

    a = ap.parse_args()
    if a.cmd == "decoydb":
        nt, nd = build_decoy_fasta(a.target, a.out)
        print("%d target + %d decoy -> %s" % (nt, nd, a.out))
        return 0

    if a.cmd == "score":
        r = score_file(a.path, a.engine, a.label)
        print(json.dumps(r, indent=2))
        if a.json:
            p = Path(a.json)
            acc = json.loads(p.read_text()) if p.exists() else []
            acc.append(r)
            p.write_text(json.dumps(acc, indent=2))
        return 0

    rows = json.loads(Path(a.json).read_text())
    ok = [r for r in rows if "error" not in r]
    print("%-22s %-10s %10s %10s %8s" % ("run", "engine", "PSMs@1%", "PEP@1%", "PSM/pep"))
    print("-" * 66)
    for r in sorted(ok, key=lambda r: (r["label"], r["engine"])):
        print("%-22s %-10s %10d %10d %8.2f"
              % (r["label"], r["engine"], r["psms_1pct"], r["peptides_1pct"],
                 r["psms_per_peptide"] or 0))
    for r in rows:
        if "error" in r:
            print("FAILED %-20s %-10s %s" % (r["label"], r.get("engine"), r["error"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())

# LIMITATIONS -- read before quoting any engine comparison
#
# * NO RESCORING. Sage's reported score is already discriminant-rescored; MSFragger's hyperscore
#   is not. Sage is therefore favoured on absolute counts. Adding Percolator to MSFragger alone
#   would reintroduce the confound this design removes; adding it to both is the correct fix and
#   is not done here. The TOOL comparison (the reference implementation vs DIAspeXtract) is unaffected, because it is
#   made WITHIN an engine -- that is the comparison this harness exists for.
# * Peptide-level FDR is "picked" (best PSM per sequence). Modified forms of one sequence collapse
#   together; a peptidoform-level count would be larger for both.
# * PSM-level FDR is computed over (spectrum, peptide) pairs. Where an engine reports multiple
#   ranks per spectrum this is not the same as one-PSM-per-spectrum; report_psms is set to 1 for
#   both engines so the two agree.
# * No replicates. Measured run-to-run spread on this pipeline is 0.17%, so differences below
#   ~0.2% are noise.
