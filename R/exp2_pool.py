# -*- coding: utf-8 -*-
"""Join the four Experiment 2 shards into the table that R/exp2_summary.py reads.

sim/run_exp2.R was run as four shards of ten replicates, with base seeds 500000, 510000,
520000 and 530000, each writing out/exp2/shard<k>.csv. This script stacks them and adds the
replicate identifier uid = <seedBase>_<rep>, which pairs an arm with the closed programme of
the same replicate.

Usage: python R/exp2_pool.py [out/exp2/pooled.csv]
"""
import csv
import glob
import io
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/exp2/pooled.csv"
FIELDS = ["design", "seedBase", "rep", "arm", "gain", "genicPct", "accLate"]

rows = []
for path in sorted(glob.glob(OUT.rsplit("/", 1)[0] + "/shard[0-9].csv")):
    with io.open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append([r[k] for k in FIELDS] + ["%s_%s" % (r["seedBase"], r["rep"])])
if len(rows) != 320:
    sys.exit("expected 320 rows (4 shards x 10 replicates x 2 designs x 4 arms), found %d"
             % len(rows))

with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, lineterminator="\r\n")   # CRLF, as in checksums/outputs.sha256
    w.writerow(FIELDS + ["uid"])
    w.writerows(rows)
print("wrote %s with %d rows" % (OUT, len(rows)))
