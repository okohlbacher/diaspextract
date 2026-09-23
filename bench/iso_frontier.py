#!/usr/bin/env python3
"""Price list of spectrum-removal rules on an emitted pseudo.mzML (exact: erasing a precursor upstream == deleting
its spectrum, because extraction is per-precursor independent). For each rule: spectra removed, peptides LOST at the
FIXED baseline threshold (a peptide is lost when every spectrum identifying it is removed), and the same for a
size- and charge-matched RANDOM control, so selectivity is visible. No re-search: the threshold shift is a separate
question answered by physically re-searching the filtered file.
Usage: iso_frontier.py <iso.tsv> <results.sage.tsv> [--ppm=20] [--seed=1]
"""
import sys, csv, bisect, random, collections
ISO = 1.0033548
iso_p, sage_p = sys.argv[1], sys.argv[2]
kw = dict(a.lstrip("-").split("=") for a in sys.argv[3:])
PPM, SEED = float(kw.get("ppm", 20)), int(kw.get("seed", 1))

S = []                                   # (index, rt, mz, z, im, npeaks, niso)
with open(iso_p) as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        S.append((int(p[0]), float(p[1]), float(p[2]), int(p[3] or 0), float(p[4] or 0), int(p[5] or 0), int(p[7] or 0)))
n = len(S)
pep_of = {}
with open(sage_p) as f:
    for x in csv.DictReader(f, delimiter="\t"):
        if x["label"] != "1" or float(x["spectrum_q"]) > 0.01: continue
        pep_of[int(x["scannr"].split("=")[-1])] = x["peptide"]
specs_of = collections.defaultdict(set)
for i, p in pep_of.items(): specs_of[p].add(i)
NPEP = len(specs_of)

byz = collections.defaultdict(list)
for s in S: byz[s[3]].append(s)
for z in byz: byz[z].sort(key=lambda s: s[2])
mzs = {z: [s[2] for s in v] for z, v in byz.items()}

def partners(s, k, rt_tol, im_tol):
    """spectra k isotope steps BELOW s (k=0: same position), same charge, co-eluting, co-mobile"""
    idx, rt, mz, z, im = s[0], s[1], s[2], s[3], s[4]
    if z <= 0: return []
    t = mz - k * ISO / z; tol = t * PPM * 1e-6
    arr, keys = byz[z], mzs[z]
    out = []
    for j in range(bisect.bisect_left(keys, t - tol), bisect.bisect_right(keys, t + tol)):
        o = arr[j]
        if o[0] != idx and abs(o[1] - rt) <= rt_tol and abs(o[4] - im) <= im_tol: out.append(o)
    return out

def lost(removed):
    return sum(1 for p, sp in specs_of.items() if sp <= removed)

def control(removed, reps=3):
    """random removal matched in size AND charge composition; mean over reps"""
    comp = collections.Counter(z_of[i] for i in removed)
    rng, tot = random.Random(SEED), 0
    for _ in range(reps):
        R = set()
        for z, c in comp.items(): R.update(rng.sample(ids_by_z[z], c))
        tot += lost(R)
    return tot / reps

z_of = {s[0]: s[3] for s in S}
ids_by_z = collections.defaultdict(list)
for s in S: ids_by_z[s[3]].append(s[0])
weaker = lambda a, b: (a[6], a[5], -a[0]) < (b[6], b[5], -b[0])     # fewer isotopes, then fewer peaks, then later index

rules = collections.OrderedDict()
def add(name, R): rules[name] = set(R)

# k = 0: same lattice position. Remove the weaker twin only (never both).
for nm, rt_tol, im_tol in (("k0 twin  dIM<=0.005 dRT<=1.5", 1.5, 0.005), ("k0 twin  dIM<=0.01  dRT<=1.5", 1.5, 0.01),
                           ("k0 twin  dIM<=0.01  dRT<=3", 3.0, 0.01), ("k0 twin  dIM<=0.02  dRT<=3", 3.0, 0.02),
                           ("k0 twin  dIM<=0.05  dRT<=3", 3.0, 0.05)):
    R = set()
    for s in S:
        for o in partners(s, 0, rt_tol, im_tol):
            if weaker(s, o) and o[0] not in R: R.add(s[0]); break
    add(nm, R)
# k = 1: remove the HEAVY member, the LIGHT member, by mobility distance
for nm, im_tol in (("k1 heavy dIM<=0.005", 0.005), ("k1 heavy dIM<=0.01", 0.01), ("k1 heavy dIM<=0.02", 0.02), ("k1 heavy dIM<=0.05", 0.05)):
    add(nm, [s[0] for s in S if partners(s, 1, 3.0, im_tol)])
for nm, im_tol in (("k1 LIGHT dIM<=0.01", 0.01), ("k1 LIGHT dIM<=0.02", 0.02)):
    R = set()
    for s in S:
        for o in partners(s, 1, 3.0, im_tol): R.add(o[0])
    add(nm, R)
add("k1..3 heavy dIM<=0.02 (census ceiling)", [s[0] for s in S if any(partners(s, k, 3.0, 0.02) for k in (1, 2, 3))])
# weak-envelope strata, no lattice at all
add("niso == 2 (minimum envelope), all", [s[0] for s in S if s[6] == 2])
add("niso == 2 AND has k1 partner below", [s[0] for s in S if s[6] == 2 and partners(s, 1, 3.0, 0.02)])

print(f"{n:,} spectra; {len(pep_of):,} identified; {NPEP:,} peptides at the fixed 1% spectrum-FDR threshold\n")
print(f"{'rule':<42}{'removed':>9}{'% spec':>8}{'dark%':>7}{'LOST':>7}{'% pep':>7}{'random':>8}{'% pep':>7}{'select.':>9}")
for nm, R in rules.items():
    if not R: print(f"{nm:<42}{0:>9}"); continue
    L, C = lost(R), control(R)
    dark = 100.0 * sum(1 for i in R if i not in pep_of) / len(R)
    sel = (C / L) if L else float("inf")
    print(f"{nm:<42}{len(R):>9,}{100*len(R)/n:>7.1f}%{dark:>6.1f}%{L:>7,}{100*L/NPEP:>6.2f}%{C:>8,.0f}{100*C/NPEP:>6.2f}%{sel:>8.2f}x")
print("\nselect. = random LOST / rule LOST: > 1 means the rule is SAFER than removing the same number of random spectra of the same charges.")
