# -*- coding: utf-8 -*-
"""Bootstrap intervals for the pairing factors and the relative variability of Table 4.

Replicates are resampled with replacement, keeping every arm and scenario of a replicate
together. 4,000 resamples, seed 20261013, 95% percentile intervals.

  pairing factor   E = [Var(Y_A) + Var(Y_C)] / Var(Y_A - Y_C), equation (17), for cumulative
                   gain: Experiment 2 (40 paired replicates, out/exp2/pooled.csv) and the
                   application (10 replicates, out/matched_controls/per_run_summary.csv)
  Table 4          median over the seven distinct contrasts of sd/|mean| of the paired
                   contrast, for the endpoint, the late mean and the late slope
                   (out/application/per_run_summary.csv), and the slope-to-endpoint ratio of medians

Usage: python R/05_bootstrap.py [out/matched_controls/bootstrap_intervals.csv]
"""
import csv
import io
import sys

import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/matched_controls/bootstrap_intervals.csv"
B, SEED = 4000, 20261013
rng = np.random.default_rng(SEED)


def read(path):
    with io.open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def factor(a, c):
    return (np.var(a, ddof=1) + np.var(c, ddof=1)) / np.var(a - c, ddof=1)


out = []


def boot(label, n, stat, point):
    idx = rng.integers(0, n, size=(B, n))
    vals = np.array([stat(i) for i in idx])
    vals = vals[np.isfinite(vals)]
    lo, hi = np.percentile(vals, [2.5, 97.5])
    out.append(dict(quantity=label, estimate=point, lo=lo, hi=hi, resamples=len(vals)))


# ---- Experiment 2: paired-only factor on cumulative gain -------------------------------
ex2 = [r for r in read("out/exp2/pooled.csv") if r["design"] == "paired"]
uids = sorted({r["uid"] for r in ex2})
g = {(r["uid"], r["arm"]): float(r["gain"]) for r in ex2}
C = np.array([g[(u, "CLOSED")] for u in uids])
for arm in ("DIRECT", "BRIDGE", "MAS"):
    A = np.array([g[(u, arm)] for u in uids])
    boot("exp2 gain factor %s" % arm, len(uids), lambda i, A=A: factor(A[i], C[i]),
         factor(A, C))

# ---- application: paired-only factor on twenty-cycle gain ------------------------------
pm = read("out/matched_controls/per_run_summary.csv")
reps = sorted({int(r["rep"]) for r in pm})
e = {(r["scenario"], r["strategy"], int(r["rep"])): float(r["gainEnd"]) for r in pm}
for sc in ("baseline", "resource_matched"):
    Cc = np.array([e[(sc, "CLOSED", k)] for k in reps])
    for arm in ("DIRECT", "BRIDGE", "MAS"):
        A = np.array([e[(sc, arm, k)] for k in reps])
        boot("application gain factor %s %s" % (sc, arm), len(reps),
             lambda i, A=A, Cc=Cc: factor(A[i], Cc[i]), factor(A, Cc))

# ---- Table 4: median relative sd of the paired contrast over the seven contrasts --------
per = read("out/application/per_run_summary.csv")
EST = ("estEndpoint", "estLateMean", "estSlope")
CONTRASTS = [("baseline", a) for a in ("ELITE2", "DIRECT", "BRIDGE", "MAS")] + \
            [("resource_matched", a) for a in ("DIRECT", "BRIDGE", "MAS")]
v = {(r["scenario"], r["strategy"], int(r["rep"]), f): float(r[f]) for r in per for f in EST}
reps = sorted({int(r["rep"]) for r in per})
D = {f: np.array([[v[(sc, a, k, f)] - v[(sc, "CLOSED", k, f)] for k in reps]
                  for sc, a in CONTRASTS]) for f in EST}      # contrasts x replicates


def median_cv(f, i):
    d = D[f][:, i]
    return np.median(np.std(d, axis=1, ddof=1) / np.abs(d.mean(axis=1)))


allrep = np.arange(len(reps))
for f in EST:
    boot("table4 median relative sd %s" % f, len(reps), lambda i, f=f: median_cv(f, i),
         median_cv(f, allrep))
boot("table4 ratio of medians slope/endpoint", len(reps),
     lambda i: median_cv("estSlope", i) / median_cv("estEndpoint", i),
     median_cv("estSlope", allrep) / median_cv("estEndpoint", allrep))

with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, ["quantity", "estimate", "lo", "hi", "resamples"],
                       lineterminator="\n")
    w.writeheader()
    for r in out:
        w.writerow({k: ("%.6f" % x if isinstance(x, float) else x) for k, x in r.items()})
for r in out:
    print("%-48s %8.3f  (%7.3f to %7.3f)" % (r["quantity"], r["estimate"], r["lo"], r["hi"]))
print("wrote %s (B = %d, seed %d)" % (OUT, B, SEED))
