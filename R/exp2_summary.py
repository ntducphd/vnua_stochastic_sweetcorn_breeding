# -*- coding: utf-8 -*-
"""Summarise experiment 2: what pairing on a common founder population removes.

The variance-reduction factor is a ratio of two sample variances on n-1 degrees of freedom
and is given with its F interval.

Usage: python R/exp2_summary.py [out/exp2/pooled.csv]
"""
import csv
import io
import statistics as st
import sys

from scipy import stats

SRC = sys.argv[1] if len(sys.argv) > 1 else "out/exp2/pooled.csv"
OUT = SRC.rsplit("/", 1)[0] + "/variance_reduction_ci.csv"
ARMS = ["DIRECT", "BRIDGE", "MAS"]
OUTCOMES = [("gain", "cumulative gain"),
            ("genicPct", "genic variance"),
            ("accLate", "late prediction accuracy")]

rows = list(csv.DictReader(io.open(SRC, encoding="utf-8")))


def series(design, arm, field):
    return {r["uid"]: float(r[field]) for r in rows
            if r["design"] == design and r["arm"] == arm
            and r[field] not in ("", "NA", "NaN")}


out = []
print("%-26s %-7s %5s %10s %10s %8s %7s %16s %s"
      % ("outcome", "arm", "n", "sd unpair", "sd paired", "cor", "VR", "95% CI on VR",
         "verdict"))
for field, label in OUTCOMES:
    for arm in ARMS:
        pA, pC = series("paired", arm, field), series("paired", "CLOSED", field)
        uA, uC = series("unpaired", arm, field), series("unpaired", "CLOSED", field)
        kp = sorted(set(pA) & set(pC))
        ku = sorted(set(uA) & set(uC))
        if len(kp) < 5 or len(ku) < 5:
            continue
        dp = [pA[k] - pC[k] for k in kp]
        du = [uA[k] - uC[k] for k in ku]
        sdP, sdU = st.stdev(dp), st.stdev(du)
        rho = st.correlation([pA[k] for k in kp], [pC[k] for k in kp])
        vr = (sdU / sdP) ** 2
        n = min(len(kp), len(ku))
        lo_f = stats.f.ppf(0.025, n - 1, n - 1)
        hi_f = stats.f.ppf(0.975, n - 1, n - 1)
        cl, ch = vr / hi_f, vr / lo_f
        verdict = ("pairing helps" if cl > 1 else
                   "pairing hurts" if ch < 1 else "not resolved")
        print("%-26s %-7s %5d %10.4f %10.4f %8.3f %7.2f %7.2f to %-6.2f %s"
              % (label, arm, n, sdU, sdP, rho, vr, cl, ch, verdict))
        out.append(dict(outcome=field, label=label, arm=arm, n=n,
                        sdUnpaired=sdU, sdPaired=sdP, corPaired=rho,
                        varianceReduction=vr, ciLo=cl, ciHi=ch,
                        nUnpairedForSamePrecision=vr * n, verdict=verdict))

with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, list(out[0].keys()))
    w.writeheader()
    w.writerows(out)
print("\nwrote %s" % OUT)

