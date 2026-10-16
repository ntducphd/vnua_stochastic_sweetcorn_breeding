# -*- coding: utf-8 -*-
"""Write out/application/pooled_records.csv, the table every result of Section 3.4 is computed from.

out/application/raw_w<seed>.csv holds, for one replicate, the per-cycle records of the five
programmes in the four scenarios (sim/run_application.R). The pooled table stacks the ten
replicates in seed order. Replicate identity is the seed; rep is written as seed - 4000.

Usage: python R/01_pool_application.py [outfile]
"""
import csv
import glob
import io
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/application/pooled_records.csv"
SEED_BASE = 4000
EXPECT = 4 * 5 * 20                 # scenarios x arms x cycles, per replicate

fields, rows = None, []
paths = sorted(glob.glob("out/application/raw_w4*.csv"))
for path in paths:
    with io.open(path, encoding="utf-8", newline="") as fh:
        rd = csv.DictReader(fh)
        if fields is None:
            fields = rd.fieldnames
        elif rd.fieldnames != fields:
            sys.exit("column mismatch in %s" % path)
        block = list(rd)
    if len(block) != EXPECT:
        sys.exit("%s: %d rows, expected %d" % (path, len(block), EXPECT))
    for r in block:
        r["rep"] = str(int(r["seed"]) - SEED_BASE)
    rows.extend(block)

with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fields)
    w.writeheader()
    w.writerows(rows)
print("wrote %s: %d rows from %d replicates" % (OUT, len(rows), len(paths)))
