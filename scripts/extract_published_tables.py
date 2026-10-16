# -*- coding: utf-8 -*-
"""Transcribe Tables 2 and 3 of Peixoto et al. (2024, G3 14: jkae128) into data/.

The article is open access under CC BY 4.0. The two tables give, for the final generation of
parents after 20 years, the mean over 50 replicates of the genetic mean and of the average genetic
variance, with a spread in parentheses, by heritability, programme (Conv, GSTC, GSF1),
genotype-by-environment interaction (absent, present) and ratio of old to new parents. Table 2
is the programme with 50 crosses per generation and Table 3 the one with 200. The rows are rebuilt
from the coordinates of the words on page 6 of the published PDF and checked against the
page image.

Columns are written in the order printed in each table: Table 2 prints 1:3, 1:1, 3:1, 0:1 and
Table 3 prints 3:1, 1:1, 1:3, 0:1. `ratio_printed` keeps the printed header. `ratio` reads both
tables in the order of Table 3, the order the article's text implies for Table 2: the text ranks
1:3 first under Conv and GSTC, which Table 2 matches only with its first and third headers
exchanged. (The authors' scripts keep 13, 25, 37 and 50 new parents of 50 for the four ratios and
do not build the tables.) The footnote of Table 2 calls the bracketed numbers the standard error of the
estimates; that of Table 3 calls them the standard deviation and the standard error of the
estimates.

Usage: python scripts/extract_published_tables.py <pdf> [out.csv]   (needs PyMuPDF)
"""
import csv
import io
import re
import sys

import pymupdf as fitz

if len(sys.argv) < 2:
    sys.exit("usage: python scripts/extract_published_tables.py <PDF of Peixoto et al. 2024, doi:10.1093/g3journal/jkae128> [out.csv]")
PDF = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "data/peixoto2024_tables2_3.csv"

page = fitz.open(PDF)[5]
words = page.get_text("words")


def rows_between(y_top, y_bot):
    ws = sorted((w for w in words if y_top <= w[1] < y_bot), key=lambda w: (round(w[1]), w[0]))
    rows, cur, cy = [], [], None
    for w in ws:
        if cy is None or abs(w[1] - cy) > 2.5:
            if cur:
                rows.append([x[4] for x in sorted(cur, key=lambda x: x[0])])
            cur, cy = [w], w[1]
        else:
            cur.append(w)
    if cur:
        rows.append([x[4] for x in sorted(cur, key=lambda x: x[0])])
    return rows


caps = sorted(w[1] for w in words if w[4] == "Table")[:2]
out = []
for table, crosses, y0, y1 in ((2, 50, caps[0], caps[1] - 60), (3, 200, caps[1], caps[1] + 200)):
    h2 = param = None
    ratios = None
    for row in rows_between(y0 + 10, y1):
        text = " ".join(row)
        if row[:2] == ["Parental", "parameters"]:
            ratios = row[3:]
            assert len(ratios) == 8, row
            continue
        if text.startswith("High heritability"):
            h2 = "high"
            continue
        if text.startswith("Low heritability"):
            h2 = "low"
            continue
        if text.startswith("Genetic mean"):
            param, row = "genetic_mean", row[2:]
        elif text.startswith("Average genetic"):
            param, row = "genetic_variance", row[2:]
        elif text.startswith("variance "):
            row = row[1:]
        if not row or row[0] not in ("Conv", "GSTC", "GSF1"):
            continue
        prog, cells = row[0], row[1:17]
        vals = [float(c) for c in cells[0::2]]
        spr = [float(re.sub(r"[()]", "", c)) for c in cells[1::2]]
        assert len(vals) == 8 and len(spr) == 8, row
        order = ["3:1", "1:1", "1:3", "0:1"]           # Table 3's order, used for both tables
        for j in range(8):
            out.append(dict(table=table, crosses=crosses, heritability=h2, parameter=param,
                            program=prog, ge="absent" if j < 4 else "present",
                            ratio_printed=ratios[j], ratio=order[j % 4], mean=vals[j],
                            spread=spr[j]))

assert len(out) == 2 * 2 * 2 * 3 * 8, len(out)
with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
print("wrote %s: %d cells (2 tables x 2 heritabilities x 2 parameters x 3 programmes x 8 columns)"
      % (OUT, len(out)))
