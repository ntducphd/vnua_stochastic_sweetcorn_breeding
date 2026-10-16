# -*- coding: utf-8 -*-
"""Compare rerun replicates with the deposited per-cycle records, field for field.

A rerun is made from the saved founder object with the deposited code, for example
    Rscript sim/run_application.R --seed 4001 --out out/rerun_4001 --founders out/founders \n        --scenarios baseline,resource_matched --pipeline-layout
Each raw_w<seed>.csv found in the directory is compared with that seed's records in
out/application/pooled_records.csv (baseline and visual sorting from run 1, the matched scenarios
from run 2), and each raw_m<seed>.csv with out/matched_controls/pooled_records.csv, for
the scenarios the rerun contains. Writes <dir>/comparison.txt, whose last line is
ALL RECORDS IDENTICAL or RECORDS DIFFER.

Usage: python scripts/compare_rerun.py [dir]          (default out/rerun_4001)
"""
import csv
import glob
import io
import os
import re
import sys

DIR = sys.argv[1] if len(sys.argv) > 1 else "out/rerun_4001"
FIELDS = ("meanG", "gain", "genicVar", "genicPct", "totalVarG", "accuracy", "donorPct", "pctSweet",
          "nAdmitted", "nCand", "nClearBar", "nClearTrue", "accDonor", "biasDonor", "admitTrueGap",
          "gapCycle", "nF2screened", "nSweetUsable", "nProgeny", "plotsGrown")
STORED = {"w": "out/application/pooled_records.csv", "m": "out/matched_controls/pooled_records.csv"}


def read(path):
    with io.open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


stored = {k: read(v) for k, v in STORED.items()}
key = lambda r: (r["scenario"], r["strategy"], r["cycle"])
lines, oks = [], []
for path in sorted(glob.glob(os.path.join(DIR, "raw_[wm]*.csv"))):
    kind, seed = re.match(r"raw_([wm])(\d+)\.csv", os.path.basename(path)).groups()
    new = read(path)
    sc = {r["scenario"] for r in new}
    old = [r for r in stored[kind] if r["seed"] == seed and r["scenario"] in sc]
    O = {key(r): r for r in old}
    same = sum(key(r) in O and all(r[f] == O[key(r)][f] for f in FIELDS) for r in new)
    ok = len(new) == len(old) and same == len(new)
    oks.append(ok)
    lines.append("%-30s seed %s: rerun %4d records, stored %4d, identical field for field %4d"
                 % ("application" if kind == "w" else "count-matched controls", seed,
                    len(new), len(old), same))
if not oks:
    sys.exit("compare_rerun: no raw_w*.csv or raw_m*.csv in %s" % DIR)
lines.append("ALL RECORDS IDENTICAL" if all(oks) else "RECORDS DIFFER")
text = "\n".join(lines) + "\n"
io.open(os.path.join(DIR, "comparison.txt"), "w", encoding="utf-8", newline="\n").write(text)
print(text, end="")
