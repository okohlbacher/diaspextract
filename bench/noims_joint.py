#!/usr/bin/env python3
"""Joint benchmark of non-IMS extraction arms: per arm and sample, extract once, then search the output four times --
Sage and MSFragger on the human target-decoy database, and both again on human + entrapment proteome -- and record
spectra, wall time, peptides at 1 % FDR and the entrapment FDR of each engine.

  noims_joint.py run   --arms arms.json --sample PXD021394 --raw X.mzML --out DIR [--threads 32]
  noims_joint.py table DIR          # markdown table of every record in DIR/records
  noims_joint.py rescore DIR        # recount MSFragger in older records from their kept TSVs
  noims_joint.py selftest

Sage is counted with the frozen harness's rules (harness/corpus.py sage_peptides, entrap_strata), the entrapment
estimator is bench/entrapment.py (peptide-hypothesis ratio for the search's digestion rule, bootstrap CI). MSFragger is
counted by msf_count (best hyperscore per peptide+mods, target-decoy competition at 1 %, classes from the protein
accession) with the same estimator: est. false targets = n_entrap / ratio. `rescore DIR` recounts older records.

arms.json: {"name": {"cmd": "bash template with {in} {out} {threads} {dir}", "env": "/abs/env.sh" (optional)}}.
--profile tims: .d input (copied as a directory) and the timsTOF search configs (tools/sage_deiso.json, msfragger.params).
"""
import argparse, csv, gzip, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time
from pathlib import Path

B = Path("/pfs/10/project/bw26g016/diaspx-bench")
HARNESS = B / "results/orbi171/harness/992f314129fb3a0f/corpus.py"
TOOLS = {"sage": B / "tools/sage/sage", "sage_config": B / "results/baseline-kit/sage_orbi.json",
         "msfragger_jar": B / "tools/MSFragger-4.4.1.jar", "msfragger_params": B / "results/baseline-kit/msfragger_orbi.params",
         "fix": B / "tools/fix_mzml_cvparams.py", "entrap_fasta": B / "db/human_entrap.fasta"}
# timsTOF (.d input, IMS path): the harness's standard search configs of the corpus runs
TIMS = {"sage_config": B / "tools/sage_deiso.json", "msfragger_params": B / "tools/msfragger.params"}
JAVA = "/opt/bwhpc/common/devel/java_jdk/21.0.2/bin/java"
TAG = "ENTRAP_"


def load(path, name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, str(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sha(p):
    with open(p, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def timed(cmd, log, cwd, env=None):
    """bash command under GNU time; returns wall s and peak RSS GB; a failure stops the arm"""
    tl = Path(cwd) / (Path(log).stem + ".time")
    with open(log, "w") as f:
        r = subprocess.run(["/usr/bin/time", "-v", "-o", str(tl), "bash", "-c", cmd], cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT)
    t = tl.read_text()
    if r.returncode:
        raise RuntimeError("%s failed (%d), see %s" % (cmd[:80], r.returncode, log))
    wall = re.search(r"Elapsed \(wall clock\) time.*: (?:(\d+):)?(\d+):([\d.]+)", t)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    return {"wall_s": round(int(wall.group(1) or 0) * 3600 + int(wall.group(2)) * 60 + float(wall.group(3)), 2),
            "rss_gb": round(int(rss.group(1)) / 1048576, 2)}


def spectra(mz):
    with open(mz, "rb") as f:
        for line in f:
            if b"<spectrumList" in line:
                return int(re.search(rb'count="(\d+)"', line).group(1))
    raise RuntimeError("no spectrumList in %s" % mz)


def isotope_share(tsv):
    """Sage: of the rank-1 target PSMs at peptide q <= 1 %, the share whose precursor needed a nonzero isotope error
    (a wrong monoisotope the search had to correct)"""
    n = nz = 0
    with open(tsv, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["rank"] == "1" and r["label"] == "1" and float(r["peptide_q"]) <= 0.01:
                n += 1
                nz += float(r["isotope_error"]) != 0.0
    return {"psms_1pct": n, "isotope_error_nonzero_pct": round(100.0 * nz / max(n, 1), 2)}


def msf_count(tsv, ratio=None, fdr=0.01):
    """MSFragger peptides at `fdr`: best hyperscore per (peptide, modifications), target-decoy competition with the cut
    taken after complete equal-score groups. A PSM is classified by its protein's accession (the first token of the
    `proteins` header: descriptions contain commas, so splitting the header would turn decoys into targets): `rev_`
    decoy, `ENTRAP_` entrapment, else target; an equal-score tie between classes counts as decoy. With `ratio`, also the
    entrapment FDR estimate of bench/entrapment.py (foreign peptides / ratio / accepted non-foreign targets)."""
    best = {}
    with (gzip.open(tsv, "rt", newline="") if str(tsv).endswith(".gz") else open(tsv, newline="")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            prot = (row.get("proteins") or row.get("protein") or "").strip()
            if not prot:
                raise RuntimeError("MSFragger TSV row without a protein (columns: %s)" % ",".join(row.keys()))
            acc = prot.split()[0]
            cls = "D" if acc.startswith("rev_") else "E" if acc.startswith(TAG) else "T"
            key = (row["peptide"], row.get("modification_info", "") or row.get("assigned_modifications", ""))
            hs = float(row["hyperscore"])
            if key not in best or hs > best[key][0] or (hs == best[key][0] and cls == "D"):
                best[key] = (hs, cls)
    groups = {}
    for hs, cls in best.values():
        groups.setdefault(hs, []).append(cls)
    d = t = e = 0
    n_t = n_e = 0
    for hs in sorted(groups, reverse=True):
        g = groups[hs]
        d += g.count("D"); t += g.count("T"); e += g.count("E")
        if t + e and d / (t + e) <= fdr:
            n_t, n_e = t, e
    out = {"peptides": n_t + n_e, "keys": len(best), "decoy_keys": sum(v[1] == "D" for v in best.values()),
           "parser": "accession"}
    if ratio is not None:
        if n_t <= 0:   # recorded, not raised: the arm's other searches stay usable
            return dict(out, error="0 target peptides at 1 %: a failed search", fdr_pct=float("nan"))
        out.update(target=n_t, entrap=n_e, fdr_pct=round(100.0 * n_e / ratio / n_t, 3))
    return out


def entrap_db(out: Path):
    """human+entrapment FASTA with reversed decoys for MSFragger (Sage makes its own), in the job's own work dir: MSFragger
    writes its peptide index beside the FASTA, and jobs of several samples would race on a shared one"""
    db = out / "db" / "human_entrap_decoy.fasta"
    if db.is_file():
        return db
    db.parent.mkdir(parents=True, exist_ok=True)
    recs, h, s = [], None, []
    for line in open(TOOLS["entrap_fasta"]):
        if line.startswith(">"):
            if h: recs.append((h, "".join(s)))
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h: recs.append((h, "".join(s)))
    tmp = db.with_suffix(".tmp%d" % os.getpid())   # jobs of several samples may build it at the same time
    with open(tmp, "w") as f:
        for h, s in recs:
            f.write(">%s\n%s\n" % (h, s))
        for h, s in recs:
            f.write(">rev_%s\n%s\n" % (h, s[::-1]))
    os.replace(tmp, db)
    return db


def run(a):
    out = Path(a.out).resolve()
    (out / "records").mkdir(parents=True, exist_ok=True)
    c, est = load(HARNESS, "corpus_frozen"), load(Path(__file__).with_name("entrapment.py"), "entrapment_est")
    arms = json.load(open(a.arms))
    work = Path(tempfile.mkdtemp(prefix="nj-", dir=os.environ.get("TMPDIR", "/tmp")))
    if a.profile == "tims":
        TOOLS.update(TIMS)
    staged = work / Path(a.raw.rstrip("/")).name
    (shutil.copytree if os.path.isdir(a.raw) else shutil.copy)(a.raw, staged)
    # entrapment configs: the standard ones with only the database swapped
    sc = json.load(open(TOOLS["sage_config"]))
    sc["database"]["fasta"] = str(TOOLS["entrap_fasta"])
    se = work / "sage_entrap.json"; se.write_text(json.dumps(sc, indent=1))
    mdb = entrap_db(work)
    me = work / "msfragger_entrap.params"
    me.write_text(re.sub(r"(?m)^database_name = .*$", "database_name = %s" % mdb, open(TOOLS["msfragger_params"]).read()))
    # the peptide-hypothesis ratio for Sage's digestion rule; MSFragger's entrapment FDR uses it too, so its rule must match
    rule = est.digest_rule(sc)
    mp = dict(re.findall(r"(?m)^(\w+)\s*=\s*([^#\n]*?)\s*(?:#.*)?$", open(TOOLS["msfragger_params"]).read()))
    mrule = {"mc": int(mp.get("allowed_missed_cleavage_1", 2)), "lo": int(mp.get("digest_min_length", 7)),
             "hi": int(mp.get("digest_max_length", 50)), "cleave": mp.get("search_enzyme_cut_1", "KR"),
             "restrict": mp.get("search_enzyme_nocut_1", "P"), "semi": int(mp.get("num_enzyme_termini", 2)) < 2,
             "c_terminal": mp.get("search_enzyme_sense_1", "C") == "C"}
    mrule["mlo"], mrule["mhi"] = (float(x) for x in mp.get("digest_mass_range", "500.0 5000.0").split())
    if any(rule[k] != v for k, v in mrule.items()):
        raise SystemExit("MSFragger and Sage digest differently (%s vs %s): one entrapment ratio does not fit both" % (mrule, rule))
    ratio = est.peptide_ratio(str(TOOLS["entrap_fasta"]), TAG, rule)
    ident = {k: sha(v) for k, v in TOOLS.items()}
    ident.update(harness=sha(HARNESS), entrapment_py=sha(Path(__file__).with_name("entrapment.py")), msfragger_entrap_db=sha(mdb),
                 joint=sha(__file__), raw=sha(staged) if staged.is_file() else "%s+%s" % (sha(staged / "analysis.tdf"),
                                                                                         sha(staged / "analysis.tdf_bin")), threads=a.threads)
    for name, arm in arms.items():
        recp = out / "records" / ("%s.%s.json" % (a.sample, name))
        env = dict(os.environ)
        if arm.get("env"):
            r = subprocess.run(["bash", "-c", 'source "$1" && env -0', "_", arm["env"]], capture_output=True, check=True)
            env = dict(kv.split("=", 1) for kv in r.stdout.decode().split("\0") if "=" in kv)
        # the arm's identity: its command, its env file, the extractor binary it resolves to (on the arm's PATH) and every
        # effective SPX_* setting; without a resolvable SPX_BIN the record is never reused
        exe = shutil.which(env["SPX_BIN"], path=env.get("PATH")) if env.get("SPX_BIN") else None
        arm_id = {"cmd": arm["cmd"], "env": sha(arm["env"]) if arm.get("env") else None, "exe": sha(exe) if exe else None,
                  "spx_env": sorted((k, v) for k, v in env.items() if k.startswith("SPX_"))}
        if recp.is_file() and exe:   # reuse a finished record only for the same inputs, tools and arm
            old = json.load(open(recp))
            if old.get("status") == "ok" and (old.get("ident"), old.get("arm_id")) == (ident, json.loads(json.dumps(arm_id))):
                continue
        d = work / name; shutil.rmtree(d, ignore_errors=True); d.mkdir()
        rec = {"sample": a.sample, "arm": name, "cmd": arm["cmd"], "env": arm.get("env"), "threads": a.threads,
               "node": os.uname().nodename, "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "ident": ident, "arm_id": arm_id,
               "entrap_ratio": ratio}
        try:
            if arm.get("env"):
                rec["arm_provenance"] = Path(arm["env"]).with_name("provenance.txt").read_text() if Path(arm["env"]).with_name("provenance.txt").is_file() else None
            mz = d / "pseudo.mzML"
            cmd = arm["cmd"].format(**{"in": staged, "out": mz, "threads": a.threads, "dir": d})
            rec["extract"] = timed(cmd, d / "extract.log", d, env)
            rec["spectra"] = spectra(mz)
            for eng in ("sage", "sage_entrap"):
                sd = d / eng; sd.mkdir()
                cfg = TOOLS["sage_config"] if eng == "sage" else se
                rec[eng] = timed("%s %s -o %s %s %s" % (TOOLS["sage"], cfg, sd, c.SAGE_TELEMETRY, mz), sd / "sage.log", sd)
                tsv = sd / "results.sage.tsv"
                if eng == "sage":
                    rec[eng].update(c.sage_peptides(tsv), **isotope_share(tsv))
                else:
                    rec[eng]["entrapment"] = c.entrap_strata(tsv, est, str(TOOLS["entrap_fasta"]))
            fixed = d / "fixed.mzML"
            subprocess.run(["python3", str(TOOLS["fix"]), str(mz), str(fixed)], check=True)
            for eng, params in (("msfragger", TOOLS["msfragger_params"]), ("msfragger_entrap", me)):
                md = d / eng; md.mkdir()
                os.link(fixed, md / "in.mzML")
                rec[eng] = timed("%s -Xmx200G -jar %s %s ./in.mzML" % (JAVA, TOOLS["msfragger_jar"], params), md / "msf.log", md)
                tsv = md / "in.tsv"
                rec[eng].update(msf_count(tsv) if eng == "msfragger" else msf_count(tsv, ratio))
            rec["status"] = "ok"
        except Exception as e:                                     # the arm fails, the next one still runs
            rec.update(status="failed", error=str(e)[-2000:])
        rec["ended"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        recp.write_text(json.dumps(rec, indent=1))
        keep = out / "tsv" / ("%s.%s" % (a.sample, name))
        keep.mkdir(parents=True, exist_ok=True)
        for eng, fn in (("sage", "results.sage.tsv"), ("sage_entrap", "results.sage.tsv"), ("msfragger", "in.tsv"), ("msfragger_entrap", "in.tsv")):
            if (d / eng / fn).is_file():
                with open(d / eng / fn, "rb") as src, gzip.open(keep / (eng + ".tsv.gz"), "wb") as dst:
                    shutil.copyfileobj(src, dst)
        shutil.rmtree(d, ignore_errors=True)
        print("%s %s %s" % (a.sample, name, rec["status"]), flush=True)
    shutil.rmtree(work, ignore_errors=True)


def table(a):
    rows = [json.load(open(p)) for p in sorted(Path(a.out, "records").glob("*.json"))]
    print("| sample | arm | spectra | wall s | Sage pep | Sage iso-err % | Sage entrap FDR % (95 % CI) | MSFragger pep | MSFragger entrap FDR % |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        if r["status"] != "ok":
            print("| %s | %s | failed: %s |||||" % (r["sample"], r["arm"], r["error"][:60])); continue
        e = r["sage_entrap"]["entrapment"]["all"]
        print("| %s | %s | %s | %s | %s | %s | %.2f (%.2f–%.2f) | %s | %.2f |" % (
            r["sample"], r["arm"], format(r["spectra"], ","), r["extract"]["wall_s"], format(r["sage"]["peptides"], ","),
            r["sage"].get("isotope_error_nonzero_pct", "-"),
            e["fdr_pct"], e["ci_pct"][0], e["ci_pct"][1], format(r["msfragger"]["peptides"], ","), r["msfragger_entrap"]["fdr_pct"]))


def rescore(a):
    """recount MSFragger peptides and entrapment FDR of every record in DIR from its kept TSVs with msf_count (records
    written before the accession rule split protein headers at commas); the old values are kept under msf_comma_split"""
    for p in sorted(Path(a.out, "records").glob("*.json")):
        r = json.load(open(p))
        tdir = Path(a.out, "tsv", "%s.%s" % (r["sample"], r["arm"]))
        if r.get("status") != "ok" or not (tdir / "msfragger.tsv.gz").is_file() or r["msfragger"].get("parser") == "accession":
            continue
        r["msf_comma_split"] = {"msfragger": r["msfragger"].get("peptides"), "msfragger_entrap": r["msfragger_entrap"].get("fdr_pct")}
        r["msfragger"].update(msf_count(tdir / "msfragger.tsv.gz"))
        r["msfragger_entrap"].update(msf_count(tdir / "msfragger_entrap.tsv.gz", r["entrap_ratio"]))
        p.write_text(json.dumps(r, indent=1))
        print("%s %s msfragger %s -> %s, entrap FDR %s -> %s" % (r["sample"], r["arm"], r["msf_comma_split"]["msfragger"],
              r["msfragger"]["peptides"], r["msf_comma_split"]["msfragger_entrap"], r["msfragger_entrap"]["fdr_pct"]))


def selftest(_=None):
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.tsv"
        rows = [("AAA", 50, "sp|H1|A_HUMAN Kinase, mitochondrial"), ("CCC", 40, "ENTRAP_sp|X|B_ARATH Synthase, chloroplastic"),
                ("DDD", 30, "ENTRAP_sp|Y|C_ARATH Protein"), ("EEE", 20, "rev_sp|H3|D_HUMAN Phosphatase, placental type"),
                ("GGG", 20, "sp|H5|E_HUMAN Protein"), ("FFF", 10, "sp|H4|F_HUMAN Protein")]
        p.write_text("peptide\thyperscore\tproteins\n" + "".join("%s\t%s\t%s\n" % r for r in rows))
        # 3 non-decoys above the first decoy; the decoy's equal-score group (with GGG) gives 1/4 > 1 %, so the cut stays above
        m = msf_count(p, 0.5)
        assert (m["target"], m["entrap"], m["peptides"], m["decoy_keys"]) == (1, 2, 3, 1), m   # a comma does not hide a decoy
        assert m["fdr_pct"] == 400.0, m              # 2 / 0.5 / 1
    print("selftest ok")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run"); r.add_argument("--arms", required=True); r.add_argument("--sample", required=True)
    r.add_argument("--raw", required=True); r.add_argument("--out", required=True); r.add_argument("--threads", type=int, default=32)
    r.add_argument("--profile", choices=["orbi", "tims"], default="orbi")
    t = sp.add_parser("table"); t.add_argument("out")
    rs = sp.add_parser("rescore"); rs.add_argument("out")
    sp.add_parser("selftest")
    a = ap.parse_args()
    {"run": run, "table": table, "rescore": rescore, "selftest": selftest}[a.cmd](a)
