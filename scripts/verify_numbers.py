# -*- coding: utf-8 -*-
"""Recompute every number the manuscript quotes, straight from the output files and scripts.

Sections follow the paper. Numbers taken from the literature (dates, the a1-sh2 distances,
the survey of Riggs et al.) are not checked here.

Usage: python scripts/verify_numbers.py
"""
import csv
import glob
import io
import math
import os
import re
import statistics as st
from collections import Counter

from scipy import stats

FAIL = []


def ok(label, got, want, tol=0.02):
    good = abs(got - want) <= max(tol, abs(want) * tol)
    print("  %-62s paper %-10s data %-10s %s"
          % (label, ("%g" % want), ("%.4g" % got), "OK" if good else "*** MISMATCH ***"))
    if not good:
        FAIL.append(label)


def exact(label, got, want, dp=2):
    """The paper prints `want` to `dp` decimals; `got` must round to it."""
    good = abs(got - want) <= 0.5 * 10 ** -dp + 1e-9
    print("  %-62s paper %-10s data %-10s %s"
          % (label, ("%.*f" % (dp, want)), ("%.6g" % got), "OK" if good else "*** MISMATCH ***"))
    if not good:
        FAIL.append(label)


def holds(label, cond):
    print("  %-84s %s" % (label, "OK" if cond else "*** DOES NOT HOLD ***"))
    if not cond:
        FAIL.append(label)


def read(path):
    return list(csv.DictReader(io.open(path, encoding="utf-8")))


def tcrit(p, df):
    return stats.t.ppf(p, df)


def tost(m, s, n, delta):
    se = s / math.sqrt(n)
    return max(stats.t.sf((m + delta) / se, n - 1), stats.t.cdf((m - delta) / se, n - 1))


def ci(x):
    n = len(x)
    m = st.mean(x)
    h = tcrit(0.975, n - 1) * st.stdev(x) / math.sqrt(n)
    return m, m - h, m + h


rows = read("out/application/pooled_records.csv")
per = read("out/application/per_run_summary.csv")
eq = read("out/application/power.csv")
cost = read("out/application/locus_cost.csv")
sup = read("out/application/supply.csv")
found = read("out/application/founder_summary.csv")
scale = read("out/application/founder_scale.csv")
t4 = {(r["scenario"], r["arm"]): r for r in read("out/application/prediction_bias_tha.csv")}
est = read("out/application/estimator_comparison.csv")
vr = read("out/exp2/variance_reduction_ci.csv")
ex2 = read("out/exp2/pooled.csv")
pr = read("out/matched_controls/per_run_summary.csv")
mc = {(r["scenario"], r["arm"], r["contrast"].split(":")[0]): r for r in
      read("out/matched_controls/matched_contrasts.csv")}
t3 = {r["estimator"]: r for r in read("out/matched_controls/estimator_distinct.csv")}
R05 = io.open("R/02_analyse_application.R", encoding="utf-8").read()

# The engine's parameter sets (sim/R/config.R), read from the engine itself
import json
import subprocess
import tempfile
with tempfile.TemporaryDirectory() as tmp_:
    cfg_json = os.path.join(tmp_, "configs.json").replace("\\", "/")
    subprocess.run(["Rscript", "-e",
                    'source("sim/R/load.R"); load_sim("sim"); jsonlite::write_json(list('
                    'application = config_application(), exp1 = config_exp1(), exp2 = config_exp2()),'
                    ' "%s", auto_unbox = TRUE, digits = NA)' % cfg_json],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    CFG = json.load(io.open(cfg_json, encoding="utf-8"))
A, E1, E2 = CFG["application"], CFG["exp1"], CFG["exp2"]

# ------------------------------------------------------------------------------------------
print("Abstract, Highlights and Section 2.1: model and parameters (Table 1)")
g, fo, pg, ig, rm = A["genome"], A["founders"], A["programme"], A["introgression"], A["resourceMatching"]
exact("QTL in the yield trait, |Q|", g["nChr"] * g["nQtlPerChr"], 1000, 0)
exact("markers in the panel, |M|", g["nChr"] * g["nSnpPerChr"], 1000, 0)
exact("chromosomes x segregating sites (application)", g["nChr"] * 1000 + g["segSites"], 10250, 0)
exact("generations since the pools split", g["splitGenerations"], 100, 0)
exact("founders simulated per pool", fo["nSimulatedPerPool"], 180, 0)
exact("founders kept per pool", fo["nKeptPerPool"], 80, 0)
exact("elite founder mean, t/ha", fo["eliteMean"], 14.0, 1)
exact("elite founder additive variance, t^2 ha^-2", fo["eliteVarA"], 1.0, 1)
exact("heritability in the founder population", fo["h2"], 0.30, 2)
holds("every run recorded h2 = 0.3", {r["h2"] for r in rows} == {"0.3"})
h2c = [float(r["h2_cycle1"]) for r in scale]
exact("heritability, first-cycle elite nursery, lowest", min(h2c), 0.12)
exact("heritability, first-cycle elite nursery, highest", max(h2c), 0.30)
exact("replicates whose first-cycle heritability is below 0.30", sum(x < 0.30 for x in h2c), 9, 0)
exact("sugar additive effect a = (31.89 - 20.25)/2", (31.89 - 20.25) / 2, 5.82)
exact("sugar additive effect in the configuration", A["sugar"]["additive"], 5.82)
exact("sugar residual SD in the configuration", A["sugar"]["residualSd"], 0.46)
exact("burn-in cycles", fo["burnInCycles"], 8, 0)
exact("parents n_p", pg["nParents"], 40, 0)
exact("crosses n_c", pg["nCrosses"], 40, 0)
exact("progeny per cross n_o", pg["nProgeny"], 15, 0)
exact("training lines per cycle", pg["nTrain"], 250, 0)
exact("training archive, records", pg["archiveSize"], 1000, 0)
exact("cycles", pg["nCycles"], 20, 0)
exact("years per cycle", pg["yearsPerCycle"], 3, 0)
holds("year = 3 x cycle in every record", all(int(r["year"]) == 3 * int(r["cycle"]) for r in rows))
exact("donor crosses per cycle", ig["nDonorCrosses"], 6, 0)
exact("F2 per donor cross", ig["f2PerCross"], 150, 0)
exact("BRIDGE families in the second step", ig["bridgeFamilies"], 20, 0)
exact("BRIDGE second-step progeny per family", ig["bridgeProgeny"], 45, 0)
exact("MAS seedlings screened, multiple of the F2", ig["masScreenFactor"], 4, 0)
exact("admission cap kappa", pg["cap"], 2, 0)
budget = pg["nCrosses"] * pg["nProgeny"]
exact("control nursery, plots", budget, 600, 0)
exact("resource matching, introgression plots", rm["introgressionPlots"], 200, 0)
exact("resource matching, progeny per elite cross", rm["nurseryProgeny"], 10, 0)
exact("resource matching, DIRECT F2 per cross", rm["directF2PerCross"], 33, 0)
exact("resource matching, BRIDGE F2 per cross in step 1", rm["bridgeF2PerCross"], 16, 0)
exact("resource matching, MAS seedlings per cross", rm["masSeedlingsPerCross"], 132, 0)
exact("resource matching, BRIDGE second-step progeny per family", rm["bridgeProgeny"], 5, 0)
holds("matched inflow arms ran with 10 progeny per cross from cycle 2",
      {r["nProgeny"] for r in rows if r["scenario"] == "resource_matched" and r["cycle"] != "1"
       and r["strategy"] in ("DIRECT", "BRIDGE", "MAS")} == {"10"})
holds("cumulative response is zero at cycle 1 in every run",
      all(abs(float(r["gain"])) < 1e-9 for r in rows if r["cycle"] == "1"))
exact("Bonferroni critical value, family of seven", tcrit(1 - 0.025 / 7, 9), 3.462, 3)
exact("Bonferroni critical value, family of six", tcrit(1 - 0.025 / 6, 9), 3.364, 3)
exact("nominal t critical value, 9 df", tcrit(0.975, 9), 2.262, 3)
holds("bootstrap: 4,000 resamples, seed 20260925",
      "B = 4000" in R05 and "set.seed(20260925)" in R05)
# Experiment 1 configuration (Table 1, Section 2.3)
G1 = E1["exp1"]["generator"]
holds("Exp 1: founder seed 20260925, programme seed 7",
      E1["seeds"]["founderSeed"] == 20260925 and E1["seeds"]["programmeSeed"] == 7)
holds("Exp 1: generators called five times for 40 individuals, 3 chromosomes, 150 sites",
      (G1["calls"], G1["nInd"], G1["nChr"], G1["segSites"]) == (5, 40, 3, 150))
holds("Exp 1: programme of 4 cycles after 3 burn-in, 120 founders per pool, 60 kept",
      E1["programme"]["nCycles"] == 4 and E1["founders"]["burnInCycles"] == 3 and
      E1["founders"]["nSimulatedPerPool"] == 120 and E1["founders"]["nKeptPerPool"] == 60)
# Experiment 2 configuration
g2, f2, p2, i2 = E2["genome"], E2["founders"], E2["programme"], E2["introgression"]
exact("Exp 2: chromosomes", g2["nChr"], 4, 0)
exact("Exp 2: segregating sites per chromosome", g2["segSites"], 120, 0)
exact("Exp 2: QTL per chromosome", g2["nQtlPerChr"], 60, 0)
exact("Exp 2: markers per chromosome", g2["nSnpPerChr"], 40, 0)
exact("Exp 2: founders simulated per pool", f2["nSimulatedPerPool"], 150, 0)
exact("Exp 2: founders kept per pool", f2["nKeptPerPool"], 40, 0)
exact("Exp 2: cycles", p2["nCycles"], 10, 0)
exact("Exp 2: burn-in cycles", f2["burnInCycles"], 5, 0)
exact("Exp 2: crosses", p2["nCrosses"], 20, 0)
exact("Exp 2: progeny per cross", p2["nProgeny"], 10, 0)
exact("Exp 2: candidates per cycle", p2["nCrosses"] * p2["nProgeny"], 200, 0)
exact("Exp 2: parents", p2["nParents"], 20, 0)
exact("Exp 2: training lines per cycle", p2["nTrain"], 150, 0)
exact("Exp 2: donor crosses", i2["nDonorCrosses"], 4, 0)
exact("Exp 2: F2 per donor cross", i2["f2PerCross"], 60, 0)
exact("Exp 2: BRIDGE families", i2["bridgeFamilies"], 10, 0)
exact("Exp 2: BRIDGE second-step progeny per family", i2["bridgeProgeny"], 24, 0)
holds("Exp 2: shard base seeds 500000, 510000, 520000, 530000, ten replicates each",
      Counter(r["seedBase"] for r in ex2) == Counter({b: 80 for b in
                                                     ("500000", "510000", "520000", "530000")})
      and {r["rep"] for r in ex2} == {str(i) for i in range(1, 11)})
exact("Exp 2: founder populations drawn (40 paired + 4 x 40 unpaired)", 40 + 4 * 40, 200, 0)
holds("Exp 2: founder populations not saved (no founder files under out/exp2)",
      not glob.glob("out/exp2/*.rds"))
# environment
holds("environment: R 4.5.1, AlphaSimR 2.1.0 recorded with every run",
      {r["RVersion"] for r in rows} == {"4.5.1"} and
      {r["AlphaSimR"] for r in rows} == {"2.1.0"} and {r["ompThreads"] for r in rows} == {"1"})
sess = io.open("environment/sessionInfo.txt", encoding="utf-8").read()
holds("environment: rrBLUP 4.6.3 on x86_64-w64-mingw32 recorded in environment/sessionInfo.txt",
      "rrBLUP" in sess and "4.6.3" in sess and "x86_64-w64-mingw32" in sess)
sizes = [os.path.getsize(p) / 1e6 for p in glob.glob("out/founders/founder_4*.rds")]
holds("ten founder objects of about 0.5 MB (%.2f to %.2f MB)" % (min(sizes), max(sizes)),
      len(sizes) == 10 and 0.4 < min(sizes) and max(sizes) < 0.65)
holds("Experiment 1 programme-layer founder object deposited (out/exp1/founder_20260925.rds)",
      os.path.exists("out/exp1/founder_20260925.rds"))

# ------------------------------------------------------------------------------------------
print("\nSection 2.6, the replicates")
holds("ten replicates, seeds 4001 to 4010",
      [int(r["seed"]) for r in found] == list(range(4001, 4011)))
holds("no run failed: 10 replicates x 4 scenarios x 5 arms x 20 cycles = 4,000 records",
      len(rows) == 4000 and Counter((r["scenario"], r["strategy"]) for r in rows) ==
      Counter({(s, a): 200 for s in ("baseline", "resource_matched", "visual_sort",
                                     "visual_matched")
               for a in ("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS")}))

# ------------------------------------------------------------------------------------------
import hashlib
import subprocess
import sys
import tempfile
with tempfile.TemporaryDirectory() as tmp:
    rebuilt = os.path.join(tmp, "pooled_records.csv")
    subprocess.run([sys.executable, "R/01_pool_application.py", rebuilt], check=True,
                   stdout=subprocess.DEVNULL)
    md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
    holds("R/01_pool_application.py rebuilds out/application/pooled_records.csv byte for byte "
          "from the ten raw_w<seed>.csv files",
          md5(rebuilt) == md5("out/application/pooled_records.csv"))

print("\nSection 2.2, the programmes and scenarios")
adm = lambda sc, a: st.mean(float(r["admitted"]) for r in pr
                            if r["scenario"] == sc and r["strategy"] == a)
for sc, arm, want in (("baseline", "ELITE2", 40), ("baseline", "DIRECT", 3.3),
                      ("baseline", "BRIDGE", 36.8), ("baseline", "MAS", 5.0),
                      ("resource_matched", "ELITE2", 40), ("resource_matched", "DIRECT", 2.1),
                      ("resource_matched", "BRIDGE", 21.5), ("resource_matched", "MAS", 3.4)):
    exact("%s %s lines admitted" % (sc, arm), adm(sc, arm), want, 1)
    if arm != "ELITE2":
        exact("%s ELITE2M_%s lines admitted" % (sc, arm), adm(sc, "ELITE2M_" + arm), want, 1)
tests_ = io.open("sim/tests/test_results.txt", encoding="utf-8").read()
holds("engine tests: ELITE2M with all-zero and all-twos schedules reproduces CLOSED and ELITE2 exactly",
      "[PASS] ELITE2M identities" in tests_ and " 0 failed" in tests_)
holds("engine tests: a run from a saved founder object is identical in a new R session",
      "[PASS] a programme run from a saved founder object is identical in a new R session" in tests_)
mrows_ = read("out/matched_controls/pooled_records.csv")
sched_ = lambda rs, sc, a, s: [r["nAdmitted"] for r in sorted(
    (r for r in rs if r["scenario"] == sc and r["strategy"] == a and r["seed"] == s),
    key=lambda r: int(r["cycle"]))]
holds("every ELITE2M arm admitted, cycle by cycle, what its inflow arm admitted (60 schedules)",
      all(sched_(mrows_, sc, "ELITE2M_" + a, s) == sched_(rows, sc, a, s) and
          len(sched_(rows, sc, a, s)) == 20
          for sc in ("baseline", "resource_matched") for a in ("DIRECT", "BRIDGE", "MAS")
          for s in [str(x) for x in range(4001, 4011)]))
c = {(r["scenario"], r["arm"]): r for r in cost}
pc = lambda sc, a: float(c[(sc, a)]["plotsPerCyc"])
excess = [pc("resource_matched", a) / 600 - 1 for a in ("DIRECT", "BRIDGE", "MAS")]
exact("matched arms, extra plots over twenty cycles, lowest (%)", 100 * min(excess), 0.7, 1)
exact("matched arms, extra plots over twenty cycles, highest (%)", 100 * max(excess), 1.3, 1)
later = [st.mean(float(r["plotsGrown"]) for r in rows if r["scenario"] == "resource_matched"
                 and r["strategy"] == a and r["cycle"] != "1") for a in ("DIRECT", "BRIDGE", "MAS")]
holds("from cycle 2 on, matched DIRECT and MAS grow about 598 plots a cycle and BRIDGE about 594 (%s)"
      % ", ".join("%.1f" % x for x in later),
      abs(later[0] - 598) < 1 and abs(later[1] - 594) < 1 and abs(later[2] - 598) < 1)
exact("visual_matched DIRECT plots per cycle", pc("visual_matched", "DIRECT"), 459, 0)
exact("visual_matched BRIDGE plots per cycle", pc("visual_matched", "BRIDGE"), 532, 0)
exact("visual_matched MAS plots per cycle", pc("visual_matched", "MAS"), 608, 0)

# ------------------------------------------------------------------------------------------
print("\nSection 3.1 and Table 2, reproducibility by layer")
e1 = read("out/exp1/within_session.csv")
lay = lambda name: [r["signature"] for r in e1 if r["layer"].startswith(name)]
exact("runMacs calls", len(lay("runMacs")), 5, 0)
exact("runMacs distinct genomes", len(set(lay("runMacs"))), 5, 0)
exact("quickHaplo calls", len(lay("quickHaplo")), 5, 0)
exact("quickHaplo distinct genomes", len(set(lay("quickHaplo"))), 1, 0)
for arm in ("CLOSED", "DIRECT", "BRIDGE", "MAS"):
    sig = lay("programme " + arm)
    holds("run_programme() %s: 6 calls, 1 distinct result to ten decimals" % arm,
          len(sig) == 6 and len(set(sig)) == 1 and len(sig[0].split("|")[0].split(".")[1]) == 10)
holds("genome signatures cover 3 x 150 = 450 sites",
      all(s.endswith("|450") for s in lay("runMacs") + lay("quickHaplo")))
# Experiment 1 across sessions (sim/run_exp1.R): one call of each
# founder generator in each of three separate R sessions
e1b = read("out/exp1/across_sessions.csv")
lab = lambda name: [r["signature"] for r in e1b if r["layer"].startswith(name)]
exact("runMacs calls in separate sessions", len(lab("runMacs")), 3, 0)
exact("runMacs distinct genomes in separate sessions", len(set(lab("runMacs"))), 3, 0)
exact("quickHaplo calls in separate sessions", len(lab("quickHaplo")), 3, 0)
holds("quickHaplo returned one genome in all eight calls, in and across sessions",
      len(set(lab("quickHaplo") + lay("quickHaplo"))) == 1)
holds("three separate sessions (distinct process ids)", len({r["pid"] for r in e1b}) == 3)
dos = [float(s.split("|")[0]) for s in lay("runMacs") + lab("runMacs")]
exact("runMacs mean dosage over all eight calls, lowest", min(dos), 0.629, 3)
exact("runMacs mean dosage over all eight calls, highest", max(dos), 0.701, 3)
holds("spread about a tenth of the value (%.3f)" % ((max(dos) - min(dos)) / st.mean(dos)),
      0.08 <= (max(dos) - min(dos)) / st.mean(dos) <= 0.13)


def traj(sc, arm):
    z = sorted((r for r in rows if r["scenario"] == sc and r["strategy"] == arm),
               key=lambda r: (int(r["rep"]), int(r["cycle"])))
    return [(float(r["gain"]), float(r["genicPct"])) for r in z]


worst, pairs = 0.0, 0
for a, b in (("baseline", "visual_sort"), ("resource_matched", "visual_matched")):
    for arm in ("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS"):
        x, y = traj(a, arm), traj(b, arm)
        worst = max([worst] + [max(abs(p[0] - q[0]), abs(p[1] - q[1])) for p, q in zip(x, y)])
        if a == "baseline":
            pairs += len(x) // 20
exact("baseline/visual-sorting replicate-by-arm pairs", pairs, 50, 0)
holds("Fig. 4 and Section 2.5: resource matching leaves CLOSED and ELITE2 identical to baseline",
      all(traj("baseline", a) == traj("resource_matched", a) for a in ("CLOSED", "ELITE2")))
holds("their gain and genic-variance trajectories are identical, also under matched resources (largest %.2g)"
      % worst, worst == 0.0)

# ------------------------------------------------------------------------------------------
print("\nSection 3.2 and Table 3, pairing")
v = {(r["outcome"], r["arm"]): r for r in vr}
# last element: replicates the unpaired design needs for the paired precision where the F
# interval excludes one, None where it includes one (Section 3.2: gain of DIRECT and BRIDGE only)
T3 = {("gain", "DIRECT"): (1.235, 0.893, 0.504, 1.91, 1.01, 3.62, 77),
      ("gain", "BRIDGE"): (1.414, 0.948, 0.552, 2.23, 1.18, 4.21, 89),
      ("gain", "MAS"): (1.274, 0.935, 0.531, 1.86, 0.98, 3.51, None),
      ("genicPct", "DIRECT"): (14.51, 12.20, 0.205, 1.42, 0.75, 2.68, None),
      ("genicPct", "BRIDGE"): (11.62, 11.11, 0.309, 1.09, 0.58, 2.07, None),
      ("genicPct", "MAS"): (12.02, 12.78, 0.125, 0.88, 0.47, 1.67, None),
      ("accLate", "DIRECT"): (0.119, 0.098, 0.251, 1.49, 0.79, 2.82, None),
      ("accLate", "BRIDGE"): (0.127, 0.116, 0.011, 1.20, 0.63, 2.27, None),
      ("accLate", "MAS"): (0.115, 0.115, -0.157, 1.01, 0.53, 1.91, None)}
for (o, a), (su, sp, co, f, lo, hi, ne) in T3.items():
    r = v[(o, a)]
    dpsd = 2 if o == "genicPct" else 3
    exact("Table 3 %s %s sd unpaired" % (o, a), float(r["sdUnpaired"]), su, dpsd)
    exact("Table 3 %s %s sd paired" % (o, a), float(r["sdPaired"]), sp, dpsd)
    exact("Table 3 %s %s correlation" % (o, a), float(r["corPaired"]), co, 3)
    exact("Table 3 %s %s variance reduction" % (o, a), float(r["varianceReduction"]), f)
    exact("Table 3 %s %s F interval, low" % (o, a), float(r["ciLo"]), lo)
    exact("Table 3 %s %s F interval, high" % (o, a), float(r["ciHi"]), hi)
    holds("Table 3 %s %s F interval excludes one only for the gain of DIRECT and BRIDGE" % (o, a),
          (float(r["ciLo"]) > 1) == (ne is not None))



def series(design, arm, field):
    return {r["uid"]: float(r[field]) for r in ex2
            if r["design"] == design and r["arm"] == arm and r[field] not in ("", "NA", "NaN")}


recomputed = 0.0
PO = {}
for (o, a) in T3:
    pA, pC = series("paired", a, o), series("paired", "CLOSED", o)
    uA, uC = series("unpaired", a, o), series("unpaired", "CLOSED", o)
    kp, ku = sorted(set(pA) & set(pC)), sorted(set(uA) & set(uC))
    dp = [pA[k] - pC[k] for k in kp]
    du = [uA[k] - uC[k] for k in ku]
    f = (st.stdev(du) / st.stdev(dp)) ** 2
    n = min(len(kp), len(ku))
    recomputed = max(recomputed, abs(f - float(v[(o, a)]["varianceReduction"])),
                     abs(f / stats.f.ppf(0.975, n - 1, n - 1) - float(v[(o, a)]["ciLo"])),
                     abs(f / stats.f.ppf(0.025, n - 1, n - 1) - float(v[(o, a)]["ciHi"])))
    A, C = [pA[k] for k in kp], [pC[k] for k in kp]
    PO[(o, a)] = (st.variance(A) + st.variance(C)) / st.variance(dp)
holds("Table 3 regenerates from out/exp2/pooled.csv within 4e-15 (largest %.2g)" % recomputed,
      recomputed <= 4e-15)
gainf_ = [PO[("gain", a)] for a in ("DIRECT", "BRIDGE", "MAS")] + \
         [float(v[("gain", a)]["varianceReduction"]) for a in ("DIRECT", "BRIDGE", "MAS")]
otherf_ = [PO[(o, a)] for o in ("genicPct", "accLate") for a in ("DIRECT", "BRIDGE", "MAS")] + \
          [float(v[(o, a)]["varianceReduction"]) for o in ("genicPct", "accLate")
           for a in ("DIRECT", "BRIDGE", "MAS")]
holds("Highlights, Table 7, Conclusions: gain factors about twofold (%.2f to %.2f), paired variance about half"
      % (min(gainf_), max(gainf_)), 1.8 < min(gainf_) and max(gainf_) < 2.3)
holds("Table 7: every genic-variance and accuracy factor is smaller than every gain factor",
      max(otherf_) < min(gainf_))
shards = []
for path in sorted(glob.glob("out/exp2/shard[0-9].csv")):
    shards += [(r["design"], r["seedBase"], r["rep"], r["arm"], r["gain"], r["genicPct"],
                r["accLate"], "%s_%s" % (r["seedBase"], r["rep"])) for r in read(path)]
holds("out/exp2/pooled.csv is the four shards joined with uid = seedBase_rep",
      shards == [(r["design"], r["seedBase"], r["rep"], r["arm"], r["gain"], r["genicPct"],
                  r["accLate"], r["uid"]) for r in ex2])
for (o, a), want in ((("gain", "DIRECT"), 1.95), (("gain", "BRIDGE"), 2.23),
                     (("gain", "MAS"), 2.13), (("genicPct", "DIRECT"), 1.25),
                     (("genicPct", "BRIDGE"), 1.44), (("genicPct", "MAS"), 1.14),
                     (("accLate", "DIRECT"), 1.33), (("accLate", "BRIDGE"), 1.01),
                     (("accLate", "MAS"), 0.86)):
    exact("paired-only factor, equation (17), %s %s" % (o, a), PO[(o, a)], want)
b9 = stats.f.ppf(1 - 0.025 / 9, 39, 39)
for a, want in (("DIRECT", 0.77), ("BRIDGE", 0.90), ("MAS", 0.75)):
    exact("Bonferroni over nine, F lower limit, gain %s" % a,
          float(v[("gain", a)]["varianceReduction"]) / b9, want)
exact("F acceptance region for a ratio of one, low (39, 39 df)", stats.f.ppf(0.025, 39, 39), 0.53)
exact("F acceptance region for a ratio of one, high", stats.f.ppf(0.975, 39, 39), 1.89)
holds("one of the six genic-variance and accuracy factors is below one, genic variance of MAS",
      [(o, a) for o in ("genicPct", "accLate") for a in ("DIRECT", "BRIDGE", "MAS")
       if float(v[(o, a)]["varianceReduction"]) < 1] == [("genicPct", "MAS")])


def fisher(rho, n=40):
    z, h = math.atanh(rho), 1.96 / math.sqrt(n - 3)
    return math.tanh(z - h), math.tanh(z + h)


cg = [float(v[("gain", a)]["corPaired"]) for a in ("DIRECT", "BRIDGE", "MAS")]
exact("gain correlations, lowest", min(cg), 0.50)
exact("gain correlations, highest", max(cg), 0.55)
exact("gain correlations, lowest Fisher limit", min(fisher(x)[0] for x in cg), 0.23)
exact("gain correlations, highest Fisher limit", max(fisher(x)[1] for x in cg), 0.74)
cgv = [float(v[("genicPct", a)]["corPaired"]) for a in ("DIRECT", "BRIDGE", "MAS")]
exact("genic-variance correlations, lowest (MAS)", min(cgv), 0.12)
exact("genic-variance correlations, highest (BRIDGE)", max(cgv), 0.31)
holds("genic-variance correlations: lowest MAS, highest BRIDGE",
      min(cgv) == cgv[2] and max(cgv) == cgv[1])
cal = [float(v[("accLate", a)]["corPaired"]) for a in ("DIRECT", "BRIDGE", "MAS")]
exact("late-accuracy correlations, lowest (MAS)", min(cal), -0.16)
exact("late-accuracy correlations, highest (DIRECT)", max(cal), 0.25)
holds("late-accuracy correlations: lowest MAS, highest DIRECT",
      min(cal) == cal[2] and max(cal) == cal[0])
lo_, hi_ = fisher(float(v[("genicPct", "MAS")]["corPaired"]))
exact("MAS genic-variance correlation, Fisher low", lo_, -0.19)
exact("MAS genic-variance correlation, Fisher high", hi_, 0.42)


def ends(sc, arm):
    z = sorted((r for r in pr if r["scenario"] == sc and r["strategy"] == arm),
               key=lambda r: int(r["rep"]))
    return [float(r["gainEnd"]) for r in z]


APP = {}
APPR = {}
for sc, arm, rho, fac in (("baseline", "DIRECT", -0.26, 0.79), ("baseline", "BRIDGE", 0.27, 1.37),
                          ("baseline", "MAS", 0.60, 2.36),
                          ("resource_matched", "DIRECT", 0.41, 1.64),
                          ("resource_matched", "BRIDGE", 0.67, 3.01),
                          ("resource_matched", "MAS", 0.45, 1.80)):
    A, C = ends(sc, arm), ends(sc, "CLOSED")
    d = [x - y for x, y in zip(A, C)]
    APPR[(sc, arm)] = st.correlation(A, C)
    exact("application correlation with CLOSED, %s %s" % (sc, arm), APPR[(sc, arm)], rho)
    APP[(sc, arm)] = (st.variance(A) + st.variance(C)) / st.variance(d)
    exact("application factor, equation (17), %s %s" % (sc, arm), APP[(sc, arm)], fac)
exact("application factors, smallest (0.79)", min(APP.values()), 0.79)
exact("application factors, largest (3.01)", max(APP.values()), 3.01)
holds("the only application factor below one is baseline DIRECT, whose correlation is negative",
      [k for k in APP if APP[k] < 1] == [("baseline", "DIRECT")] and APPR[("baseline", "DIRECT")] < 0)
exact("Section 4.4: application correlations, lowest", min(APPR.values()), -0.26)
exact("Section 4.4: application correlations, highest", max(APPR.values()), 0.67)

# ------------------------------------------------------------------------------------------
# out/matched_controls/bootstrap_intervals.csv (R/05_bootstrap.py, seed 20261013): the application
# factors of Fig. 1d and the Experiment 2 paired-only factors. Table 3 uses the intervals of
# R/06_inference.py (seed 20261014, keys po_gain_*).
bs = {r["quantity"]: r for r in read("out/matched_controls/bootstrap_intervals.csv")}
holds("bootstrap file: 4,000 resamples, seed 20261013, kept in R/05_bootstrap.py",
      all(r["resamples"] == "4000" for r in bs.values()) and
      "B, SEED = 4000, 20261013" in io.open("R/05_bootstrap.py", encoding="utf-8").read())
for a, lo, hi in (("DIRECT", 1.17, 3.17), ("BRIDGE", 1.18, 4.32), ("MAS", 1.34, 3.88)):
    r = bs["exp2 gain factor %s" % a]
    holds("bootstrap point estimate equals equation (17), Exp 2 %s" % a,
          abs(float(r["estimate"]) - PO[("gain", a)]) < 1e-5)
    exact("Exp 2 paired-only factor %s, bootstrap low" % a, float(r["lo"]), lo)
    exact("Exp 2 paired-only factor %s, bootstrap high" % a, float(r["hi"]), hi)
inc = [float(bs["application gain factor %s %s" % (s_, a)]["lo"]) < 1 <
       float(bs["application gain factor %s %s" % (s_, a)]["hi"])
       for s_ in ("baseline", "resource_matched") for a in ("DIRECT", "BRIDGE", "MAS")]
exact("application factors whose bootstrap interval includes one", sum(inc), 4, 0)
holds("every application bootstrap interval includes the Exp 2 paired-only factor of the same arm",
      all(float(bs["application gain factor %s %s" % (s_, a)]["lo"]) < PO[("gain", a)] <
          float(bs["application gain factor %s %s" % (s_, a)]["hi"])
          for s_ in ("baseline", "resource_matched") for a in ("DIRECT", "BRIDGE", "MAS")))
for s_, a in APP:
    holds("bootstrap point estimate equals equation (17), application %s %s" % (s_, a),
          abs(float(bs["application gain factor %s %s" % (s_, a)]["estimate"]) - APP[(s_, a)]) < 1e-5)

print("\nSection 3.3 and Table 4, three response statistics")
for e_, est_, lo, hi in (("estEndpoint", 4.09, 1.31, 11.27), ("estLateMean", 2.83, 1.28, 8.86),
                         ("estSlope", 3.93, 1.03, 10.48)):
    r = bs["table4 median relative sd %s" % e_]
    holds("bootstrap point estimate equals Table 4 median, %s" % e_,
          abs(float(r["estimate"]) - float(t3[e_]["medianRelSd"])) < 1e-5)
    exact("Table 4 %s median, bootstrap low" % e_, float(r["lo"]), lo)
    exact("Table 4 %s median, bootstrap high" % e_, float(r["hi"]), hi)
r = bs["table4 ratio of medians slope/endpoint"]
exact("slope-to-endpoint ratio of medians", float(r["estimate"]), 0.96, 2)
exact("ratio of medians, bootstrap low", float(r["lo"]), 0.16)
exact("ratio of medians, bootstrap high", float(r["hi"]), 4.88)
for e, med, lo, hi, least, most in (("estEndpoint", 4.09, 1.81, 550.13, 1, 3),
                                    ("estLateMean", 2.83, 1.34, 27.15, 3, 1),
                                    ("estSlope", 3.93, 1.52, 23.15, 3, 3)):
    r = t3[e]
    exact("%s contrasts" % e, float(r["nContrasts"]), 7, 0)
    exact("%s median relative sd" % e, float(r["medianRelSd"]), med)
    exact("%s smallest relative sd" % e, float(r["minRelSd"]), lo)
    exact("%s largest relative sd" % e, float(r["maxRelSd"]), hi)
    exact("%s smallest in" % e, float(r["leastVariable"]), least, 0)
    exact("%s largest in" % e, float(r["mostVariable"]), most, 0)
dist = [r for r in est if r["scenario"] == "baseline" or
        (r["scenario"] == "resource_matched" and r["arm"] != "ELITE2")]
holds("seven distinct contrasts out of sixteen nominal ones",
      len(dist) == 21 and len(est) == 48)
cvt = {}
for r in dist:
    cvt.setdefault((r["scenario"], r["arm"]), {})[r["estimator"]] = r
rs = [x["estSlope"]["cv"] for x in cvt.values()]
ratio_s = [float(x["estSlope"]["cv"]) / float(x["estEndpoint"]["cv"]) for x in cvt.values()]
ratio_l = [float(x["estLateMean"]["cv"]) / float(x["estEndpoint"]["cv"]) for x in cvt.values()]
exact("slope/endpoint relative-sd ratio by contrast, median", st.median(ratio_s), 1.07)
exact("slope/endpoint relative-sd ratio by contrast, lowest", min(ratio_s), 0.00)
exact("slope/endpoint relative-sd ratio by contrast, highest", max(ratio_s), 2.18)
exact("late-mean/endpoint relative-sd ratio by contrast, median", st.median(ratio_l), 0.74)
exact("late-mean/endpoint relative-sd ratio by contrast, lowest", min(ratio_l), 0.02)
exact("late-mean/endpoint relative-sd ratio by contrast, highest", max(ratio_l), 1.26)
exact("contrasts where the slope's relative sd exceeds the endpoint's",
      sum(x > 1 for x in ratio_s), 4, 0)
sq = [x * x for x in ratio_s]
exact("squared relative-sd ratio, slope over endpoint, lowest", min(sq), 0.0, 1)
exact("squared relative-sd ratio, slope over endpoint, highest", max(sq), 5, 0)
exact("squared relative-sd ratio, slope over endpoint, median", st.median(sq), 1.1, 1)
# random-walk benchmark: zero at cycle 1, 19 independent increments to cycle 20
k = [15, 16, 17, 18, 19]
w = [-2, -1, 0, 1, 2]
rel_end = math.sqrt(19) / 19
rel_late = math.sqrt(sum(min(a, b) for a in k for b in k) / 25) / (sum(k) / 5)
rel_slope = math.sqrt(sum(w[i] * w[j] * min(k[i], k[j]) for i in range(5)
                          for j in range(5)) / 100)
exact("benchmark: late mean relative to the endpoint", rel_late / rel_end, 1.03)
exact("benchmark: slope relative to the endpoint", rel_slope / rel_end, 2.22)
exact("benchmark: replicates for the slope, times the endpoint's", (rel_slope / rel_end) ** 2, 4.9, 1)
sl = {k_: float(x["estSlope"]["meanDiff"]) for k_, x in cvt.items()}
en = {k_: float(x["estEndpoint"]["meanDiff"]) for k_, x in cvt.items()}
exact("mean slope contrasts, largest (t/ha/yr)", max(sl.values()), 0.021, 3)
exact("mean slope contrasts, smallest (t/ha/yr)", min(sl.values()), 0.002, 3)
holds("Section 3.3: mean late slope positive and mean endpoint negative on all seven contrasts",
      all(x > 0 for x in sl.values()) and all(x < 0 for x in en.values()))
worst_k = max(cvt, key=lambda k_: float(cvt[k_]["estSlope"]["cv"]))
holds("the worst slope value comes from the slope contrast smallest in magnitude (baseline MAS)",
      worst_k == min(sl, key=lambda k_: abs(sl[k_])) == ("baseline", "MAS"))
holds("on resource-matched MAS, where the slope contrast is largest, the slope is least variable",
      max(sl, key=lambda k_: abs(sl[k_])) == ("resource_matched", "MAS") and
      float(cvt[("resource_matched", "MAS")]["estSlope"]["cv"]) ==
      min(float(cvt[("resource_matched", "MAS")][e]["cv"]) for e in ("estEndpoint", "estLateMean",
                                                                     "estSlope")))

# ------------------------------------------------------------------------------------------
print("\nSection 3.4.1, the simulated base")
chrs = Counter(int(r["locusChr"]) for r in found)
exact("chromosomes carrying the endosperm locus", len(chrs), 7, 0)
exact("replicates with the locus on chromosome 10", chrs[10], 0, 0)
exact("replicates with the locus on chromosome 1", chrs[1], 1, 0)
holds("locusChr in the founder objects agrees with the summary",
      [int(r["locusChr"]) for r in scale] == [int(r["locusChr"]) for r in found])
f_ = lambda k_: [float(r[k_]) for r in found]
for k_, m, lo, hi in (("baseMean", 20.41, 18.66, 21.43), ("donorMean", 14.34, 13.04, 16.99),
                      ("realisedGap", 6.07, 1.67, 8.10)):
    exact("%s, mean" % k_, st.mean(f_(k_)), m)
    exact("%s, lowest" % k_, min(f_(k_)), lo)
    exact("%s, highest" % k_, max(f_(k_)), hi)
sda = [float(r["sdA_tha"]) for r in scale]
exact("post-burn-in additive SD, lowest t/ha", min(sda), 0.69)
exact("post-burn-in additive SD, highest t/ha", max(sda), 1.07)
nsd = st.mean(float(a["gap"]) / float(a["sdA_tha"]) for a in scale)
holds("gap of on average about seven additive SDs (mean %.2f)" % nsd, 6.5 <= nsd < 7.5)
gend = [float(r["gapEnd"]) for r in per if r["scenario"] == "baseline" and r["strategy"] == "CLOSED"]
exact("closed programme above the donor pool at cycle 20, t/ha", st.mean(gend), 14.54)
cul = f_("cullDonor")
exact("donor founders culled, lowest %", 100 * min(cul), 0, 0)
exact("donor founders culled, highest % (up to 21%)", 100 * max(cul), 21, 0)
nr = [int(r["nRedraw"]) for r in found]
exact("locus redraws, fewest", min(nr), 0, 0)
exact("locus redraws, most", max(nr), 17, 0)
exact("locus redraws, median", st.median(nr), 2, 0)
exact("draws in which the locus fell on a QTL or marker", sum(nr), 44, 0)
exact("draws made", sum(nr) + len(nr), 54, 0)
exact("QTL plus marker sites per chromosome", g["nQtlPerChr"] + g["nSnpPerChr"], 200, 0)
holds("44 of 54 exceeds the 200/250 a random site would give", 44 / 54 > 200 / 250)

print("\nSection 3.4.1, the recessive constraint")
for arm, scr_w, sw_w, frac, pval in (("DIRECT", 180000, 45128, 0.251, 0.49),
                                     ("MAS", 720000, 179699, 0.250, 0.41)):
    z = [r for r in rows if r["scenario"] == "baseline" and r["strategy"] == arm]
    scr = sum(int(float(r["nF2screened"])) for r in z)
    sw = sum(int(float(r["nSweetUsable"])) for r in z)
    exact("%s F2 screened" % arm, scr, scr_w, 0)
    exact("%s sweet homozygotes" % arm, sw, sw_w, 0)
    exact("%s recovery fraction" % arm, sw / scr, frac, 3)
    exact("%s exact binomial p against 1/4" % arm, stats.binomtest(sw, scr, 0.25).pvalue, pval,
          3 if pval < 0.1 else 2)
holds("neither p is below 0.025", all(stats.binomtest(s_, n_, 0.25).pvalue > 0.025
                                      for s_, n_ in ((45128, 180000), (179699, 720000))))
z = [r for r in rows if r["scenario"] == "baseline" and r["strategy"] == "BRIDGE"]
exact("BRIDGE recovery fraction", sum(float(r["nSweetUsable"]) for r in z) /
      sum(float(r["nF2screened"]) for r in z), 0.6263, 4)
exact("BRIDGE expected fraction (900 x 0.25 + 900 x 1)/1800", (900 * 0.25 + 900) / 1800, 0.625, 3)
holds("every elite pool homozygous sweet in every cycle of every run",
      all(abs(float(r["pctSweet"]) - 100) < 1e-9 for r in rows))

print("\nSection 3.4.2, cost under the two discard assumptions")
for sc, arm, want in (("baseline", "DIRECT", 1500), ("baseline", "BRIDGE", 2400),
                      ("visual_sort", "DIRECT", 826), ("visual_sort", "BRIDGE", 1727),
                      ("baseline", "MAS", 1498), ("visual_sort", "MAS", 1498)):
    exact("%s %s plots per cycle" % (sc, arm), pc(sc, arm), want, 0)
exact("introgression plots, DIRECT field discard", pc("baseline", "DIRECT") - 600, 900, 0)
exact("introgression plots, BRIDGE field discard", pc("baseline", "BRIDGE") - 600, 1800, 0)
exact("introgression plots, DIRECT visual sorting", pc("visual_sort", "DIRECT") - 600, 226, 0)
exact("introgression plots, BRIDGE visual sorting", pc("visual_sort", "BRIDGE") - 600, 1127, 0)
exact("DIRECT introgression nursery cut, fold", 900 / (pc("visual_sort", "DIRECT") - 600), 4.0, 1)
exact("BRIDGE introgression nursery cut, %", 100 * (1 - (pc("visual_sort", "BRIDGE") - 600) / 1800),
      37, 0)
for sc, arm, est_, lo, hi in (("baseline", "DIRECT", 9091, 7143, 12000),
                              ("visual_sort", "DIRECT", 5004, 3934, 6876),
                              ("baseline", "BRIDGE", 1304, 1263, 1356),
                              ("visual_sort", "BRIDGE", 939, 909, 976),
                              ("visual_sort", "MAS", 5994, 5259, 7132),
                              ("baseline", "MAS", 5994, 5259, 6976)):
    r = c[(sc, arm)]
    exact("%s %s plots per admitted line" % (sc, arm), float(r["plotsPerLine"]), est_, 0)
    exact("%s %s bootstrap low" % (sc, arm), float(r["plotsPerLineLo"]), lo, 0)
    exact("%s %s bootstrap high" % (sc, arm), float(r["plotsPerLineHi"]), hi, 0)
ipl = lambda sc, a: (pc(sc, a) - 600) * 20 / float(c[(sc, a)]["admitted"])
for sc, arm, want in (("baseline", "DIRECT", 5455), ("visual_sort", "DIRECT", 1368),
                      ("baseline", "BRIDGE", 978), ("visual_sort", "BRIDGE", 613)):
    exact("%s %s introgression plots per admitted line" % (sc, arm), ipl(sc, arm), want, 0)
exact("DIRECT replicates admitting no line, baseline", float(c[("baseline", "DIRECT")]["zeroAdmitReps"]), 0, 0)
exact("MAS plots over DIRECT plots under visual sorting", pc("visual_sort", "MAS") /
      pc("visual_sort", "DIRECT"), 1.8, 1)
holds("MAS costs most per admitted line under visual sorting; intermediate under field discard",
      float(c[("visual_sort", "MAS")]["plotsPerLine"]) > max(
          float(c[("visual_sort", a)]["plotsPerLine"]) for a in ("DIRECT", "BRIDGE")) and
      float(c[("baseline", "BRIDGE")]["plotsPerLine"]) < float(c[("baseline", "MAS")]["plotsPerLine"])
      < float(c[("baseline", "DIRECT")]["plotsPerLine"]))
holds("its interval overlaps that of direct crossing under visual sorting",
      float(c[("visual_sort", "MAS")]["plotsPerLineLo"]) < float(c[("visual_sort", "DIRECT")]["plotsPerLineHi"]))
A, B = ends("baseline", "MAS"), ends("baseline", "DIRECT")
m, lo, hi = ci([x - y for x, y in zip(A, B)])
exact("MAS minus DIRECT gain, baseline", m, 0.55)
exact("MAS minus DIRECT gain, low", lo, -0.60)
exact("MAS minus DIRECT gain, high", hi, 1.69)

print("\nSection 3.4.3 and Table 5, the predictor inside the candidate pool")
sd = {(r["scenario"], r["arm"]): r for r in sup}
for arm, cand, cb, ct, acc in (("DIRECT", 226, 107, 139, 0.302), ("MAS", 898, 408, 495, 0.334),
                               ("BRIDGE", 900, 12373, 5843, 0.284), ("ELITE2", 600, 3213, 27504, None)):
    r = sd[("baseline", arm)]
    exact("%s candidates per cycle" % arm, float(r["cand"]), cand, 0)
    exact("%s clearing on GEBV, total over 200 replicate-cycles" % arm, 200 * float(r["clearBar"]), cb, 0)
    exact("%s clearing on true value, total" % arm, 200 * float(r["clearTrue"]), ct, 0)
    if acc is not None:
        exact("%s accuracy within the pool" % arm, float(r["accDonor"]), acc, 3)
for arm, b, blo, bhi, gp, glo, ghi, ng in (("DIRECT", 4.36, 3.61, 5.12, -3.88, -5.27, -2.48, 10),
                                           ("MAS", 4.40, 3.71, 5.10, -4.34, -4.93, -3.75, 10),
                                           ("BRIDGE", 1.92, 1.64, 2.20, -1.41, -1.61, -1.21, 10),
                                           ("ELITE2", None, None, None, -0.09, -0.14, -0.04, 10)):
    r = t4[("baseline", arm)]
    if b is not None:
        exact("%s GEBV minus true value, t/ha" % arm, float(r["bias"]), b)
        exact("%s bias, low" % arm, float(r["biasLo"]), blo)
        exact("%s bias, high" % arm, float(r["biasHi"]), bhi)
    exact("%s admitted minus true threshold, t/ha" % arm, float(r["gap"]), gp)
    exact("%s gap, low" % arm, float(r["gapLo"]), glo)
    exact("%s gap, high" % arm, float(r["gapHi"]), ghi)
    exact("%s replicates averaged for the gap" % arm, float(r["gapN"]), ng, 0)
# the t/ha values follow from the internal ones and each replicate's scale
sc_ = {i + 1: float(r["scale"]) for i, r in enumerate(scale)}
for arm in ("DIRECT", "MAS", "BRIDGE"):
    z = [r for r in per if r["scenario"] == "baseline" and r["strategy"] == arm]
    conv = st.mean(float(r["biasDonor"]) * sc_[int(r["rep"])] for r in z)
    holds("%s bias recomputed from per-run values and scales (%.3f)" % (arm, conv),
          abs(conv - float(t4[("baseline", arm)]["bias"])) < 1e-9)
sca = [float(r["scale"]) for r in scale]
holds("one internal unit is 0.89 to 1.82 t/ha across replicates (%.3f to %.3f)" % (min(sca), max(sca)),
      round(min(sca), 2) == 0.89 and round(max(sca), 2) == 1.82)
exact("half-donor bias, Discussion (4.4)", float(t4[("baseline", "DIRECT")]["bias"]), 4.4, 1)
exact("half-donor bias, MAS (4.4)", float(t4[("baseline", "MAS")]["bias"]), 4.4, 1)
exact("quarter-donor bias (1.9)", float(t4[("baseline", "BRIDGE")]["bias"]), 1.9, 1)
exact("bias cut by halving the donor fraction, %",
      100 * (1 - float(t4[("baseline", "BRIDGE")]["bias"]) / float(t4[("baseline", "DIRECT")]["bias"])), 56, 0)
exact("admitted lines below the parents, DIRECT (3.88)", -float(t4[("baseline", "DIRECT")]["gap"]), 3.88)
exact("admitted lines below the parents, MAS (4.34)", -float(t4[("baseline", "MAS")]["gap"]), 4.34)
exact("admitted lines below the parents, BRIDGE (1.41)", -float(t4[("baseline", "BRIDGE")]["gap"]), 1.41)
gaps_ = [-float(t4[("baseline", a)]["gap"]) for a in ("DIRECT", "MAS", "BRIDGE")]
exact("Section 4.6: admitted lines short of the true threshold, least (1.4)", min(gaps_), 1.4, 1)
exact("Section 4.6: admitted lines short of the true threshold, most (4.3)", max(gaps_), 4.3, 1)
exact("elite lines ELITE2 admits below the parents (0.09)", -float(t4[("baseline", "ELITE2")]["gap"]), 0.09)
rat = {a: float(sd[("baseline", a)]["clearBar"]) / float(sd[("baseline", a)]["clearTrue"])
       for a in ("DIRECT", "MAS", "BRIDGE")}
exact("DIRECT ratio, GEBV over true value", rat["DIRECT"], 0.8, 1)
exact("MAS ratio", rat["MAS"], 0.8, 1)
exact("BRIDGE ratio", rat["BRIDGE"], 2.1, 1)
holds("fewer DIRECT and MAS candidates clear the threshold on GEBV than on true value",
      rat["DIRECT"] < 1 and rat["MAS"] < 1 and float(t4[("baseline", "DIRECT")]["bias"]) > 0
      and float(t4[("baseline", "MAS")]["bias"]) > 0)
exact("ELITE2 clearing on GEBV per cycle", float(sd[("baseline", "ELITE2")]["clearBar"]), 16.1, 1)
exact("ELITE2 clearing on true value per cycle", float(sd[("baseline", "ELITE2")]["clearTrue"]), 137.5, 1)
accE = st.mean(float(r["accuracy"]) for r in rows if r["scenario"] == "baseline" and r["strategy"] == "ELITE2")
exact("elite accuracy, ELITE2 baseline (about 0.41)", accE, 0.41)
p_ = 40 / 600
zc = stats.norm.ppf(1 - p_)
i_ = stats.norm.pdf(zc) / p_
exact("normal model, lines above the GEBV threshold", 600 * stats.norm.sf(i_), 16, 0)
exact("normal model, lines above the true threshold", 600 * stats.norm.sf(i_ * accE), 128, 0)
z = [r for r in rows if r["scenario"] == "baseline" and r["strategy"] == "DIRECT"]
holds("DIRECT admitted in 20 cycles, in all ten replicates",
      sum(r["nAdmitted"] != "0" for r in z) == 20 and
      len({r["rep"] for r in z if r["nAdmitted"] != "0"}) == 10)

print("\nSection 3.4.4, Table 6 and Figure 5, the control arms")
e = {(r["scenario"], r["arm"]): float(r["diff"]) for r in eq}
for sc, arm, want in (("baseline", "ELITE2", -0.00), ("baseline", "DIRECT", -0.60),
                      ("baseline", "BRIDGE", -0.34), ("baseline", "MAS", -0.05)):
    exact("%s %s gain difference" % (sc, arm), e[(sc, arm)], want)
m, lo, hi = ci([x - y for x, y in zip(ends("baseline", "ELITE2"), ends("baseline", "CLOSED"))])
exact("ELITE2 minus CLOSED, low", lo, -0.69)
exact("ELITE2 minus CLOSED, high", hi, 0.69)
exact("Section 4.2: ELITE2 minus CLOSED, half-width of the 95% interval (0.69)", hi - m, 0.69)
ad = [int(float(r["admitted"])) for r in pr if r["scenario"] == "baseline" and r["strategy"] == "DIRECT"]
exact("DIRECT admissions per replicate, fewest", min(ad), 1, 0)
exact("DIRECT admissions per replicate, most", max(ad), 6, 0)
T6 = {("baseline", "DIRECT"): ((-0.60, -1.87, 0.67), (0.09, -0.68, 0.86), (-0.70, -1.71, 0.32), (-0.60, -1.40, 0.20)),
      ("baseline", "BRIDGE"): ((-0.34, -1.26, 0.57), (-0.23, -1.06, 0.60), (-0.12, -0.78, 0.54), (-0.34, -0.95, 0.26)),
      ("baseline", "MAS"): ((-0.05, -0.90, 0.79), (0.01, -0.75, 0.78), (-0.07, -1.02, 0.88), (-0.05, -0.58, 0.48)),
      ("resource_matched", "DIRECT"): ((-0.15, -0.89, 0.59), (-0.14, -0.80, 0.51), (-0.01, -0.62, 0.61), (-0.15, -0.70, 0.41)),
      ("resource_matched", "BRIDGE"): ((-0.50, -1.15, 0.15), (-0.25, -1.14, 0.64), (-0.25, -0.99, 0.49), (-0.50, -1.18, 0.18)),
      ("resource_matched", "MAS"): ((-0.26, -1.04, 0.51), (-0.27, -0.91, 0.37), (0.01, -0.52, 0.54), (-0.26, -0.86, 0.33))}
for (sc, arm), cells in T6.items():
    for key, (d, lo, hi) in zip(("total", "recycling", "donor", "unmatched donor"), cells):
        r = mc[(sc, arm, key)]
        exact("Table 6 %s %s %s" % (sc, arm, key), float(r["diff"]), d)
        exact("Table 6 %s %s %s, low" % (sc, arm, key), float(r["lo"]), lo)
        exact("Table 6 %s %s %s, high" % (sc, arm, key), float(r["hi"]), hi)
donor = [mc[k_] for k_ in mc if k_[2] == "donor"]
recyc = [mc[k_] for k_ in mc if k_[2] == "recycling"]
unm = [mc[k_] for k_ in mc if k_[2] == "unmatched donor"]
holds("every donor-contrast interval spans zero", all(float(r["lo"]) < 0 < float(r["hi"]) for r in donor))
holds("every unmatched donor interval spans zero", all(float(r["lo"]) < 0 < float(r["hi"]) for r in unm))
exact("smallest donor contrast", min(float(r["diff"]) for r in donor), -0.70)
exact("largest donor contrast", max(float(r["diff"]) for r in donor), 0.01)
ub = [float(mc[("baseline", a, "unmatched donor")]["diff"]) for a in ("DIRECT", "BRIDGE", "MAS")]
um = [float(mc[("resource_matched", a, "unmatched donor")]["diff"]) for a in ("DIRECT", "BRIDGE", "MAS")]
exact("unmatched donor reading, baseline, smallest", min(ub), -0.60)
exact("unmatched donor reading, baseline, largest", max(ub), -0.05)
exact("unmatched donor reading, matched, smallest", min(um), -0.50)
exact("unmatched donor reading, matched, largest", max(um), -0.15)
exact("unmatched donor reading over both scenarios, smallest", min(ub + um), -0.60)
exact("unmatched donor reading over both scenarios, largest", max(ub + um), -0.05)
holds("no recycling interval excludes zero, at baseline or under matched resources",
      [(r["scenario"], r["arm"]) for r in recyc if not float(r["lo"]) < 0 < float(r["hi"])] == [])
holds("no recycling or donor interval excludes zero after Bonferroni over six",
      all(float(r["bonLo"]) < 0 < float(r["bonHi"]) for r in recyc + donor))
for sc, arm, want, lo, hi in (("baseline", "DIRECT", 0.10, -0.52, 0.72),
                              ("baseline", "MAS", 0.02, -0.59, 0.62),
                              ("baseline", "BRIDGE", -0.23, -0.53, 0.08)):
    m, l_, h_ = ci([x - y for x, y in zip(ends(sc, "ELITE2M_" + arm), ends(sc, "ELITE2"))])
    exact("ELITE2M minus ELITE2, %s %s" % (sc, arm), m, want)
    exact("ELITE2M minus ELITE2, %s %s, low" % (sc, arm), l_, lo)
    exact("ELITE2M minus ELITE2, %s %s, high" % (sc, arm), h_, hi)
    holds("ELITE2M minus ELITE2 at baseline, %s, spans zero" % arm, l_ < 0 < h_)
mm = []
for arm in ("DIRECT", "BRIDGE", "MAS"):
    m, l_, h_ = ci([x - y for x, y in zip(ends("resource_matched", "ELITE2M_" + arm),
                                          ends("resource_matched", "ELITE2"))])
    mm.append(m)
    holds("ELITE2M minus ELITE2 under matched resources, %s, spans zero" % arm, l_ < 0 < h_)
exact("ELITE2M minus ELITE2 under matched resources, smallest", min(mm), -0.27)
exact("ELITE2M minus ELITE2 under matched resources, largest", max(mm), -0.14)
holds("ELITE2 keeps 600 plots, each matched ELITE2M 400 from cycle 2 on (15 and 10 progeny per cross)",
      {r["nProgeny"] for r in rows if r["scenario"] == "resource_matched" and r["strategy"] == "ELITE2"} == {"15"}
      and {r["nProgeny"] for r in read("out/matched_controls/pooled_records.csv")
           if r["scenario"] == "resource_matched" and r["cycle"] != "1"} == {"10"})
# replicate 4007, matched resources: DIRECT admits nothing, and its ELITE2M replays an all-zero
# schedule on the same 400-plot nursery, so the two differ only through their random draws
z7 = {r["strategy"]: r for r in pr if r["scenario"] == "resource_matched" and r["rep"] == "7"}
holds("seed 4007 is replicate 7 in the pooled records",
      {r["seed"] for r in rows if r["rep"] == "7"} == {"4007"})
exact("seed 4007, DIRECT lines admitted under matched resources", float(z7["DIRECT"]["admitted"]), 0, 0)
exact("seed 4007, its ELITE2M lines admitted", float(z7["ELITE2M_DIRECT"]["admitted"]), 0, 0)
exact("seed 4007, DIRECT minus its ELITE2M at cycle 20, t/ha",
      float(z7["DIRECT"]["gainEnd"]) - float(z7["ELITE2M_DIRECT"]["gainEnd"]), 0.60)

print("\nSections 2.6 and 3.4.5, power, equivalence and multiplicity")
mde = [float(r["mde80"]) for r in eq]
exact("smallest MDE against CLOSED", min(mde), 0.90)
exact("largest MDE against CLOSED", max(mde), 1.76)
G = st.mean(float(r["estEndpoint"]) for r in per if r["scenario"] == "baseline" and r["strategy"] == "CLOSED")
exact("closed programme's twenty-cycle gain, t/ha", G, 8.46)
exact("smallest MDE as % of that gain", 100 * min(mde) / G, 11, 0)
exact("largest MDE as % of that gain", 100 * max(mde) / G, 21, 0)
dmde = [(tcrit(0.975, 9) + tcrit(0.80, 9)) * float(r["sd"]) / math.sqrt(10) for r in donor]
exact("smallest MDE, donor contrasts", min(dmde), 0.74)
exact("largest MDE, donor contrasts", max(dmde), 1.41)
exact("Section 4.4: smallest MDE over both families (donor)", min(mde + dmde), 0.74)
exact("Section 4.4: largest MDE over both families (against CLOSED)", max(mde + dmde), 1.76)
# (Equivalence bounds are computed in R/06_inference.py.)
seven = [r for r in eq if r["scenario"] == "baseline" or
         (r["scenario"] == "resource_matched" and r["arm"] != "ELITE2")]
exact("distinct contrasts in the family", len(seven), 7, 0)
sig = [r for r in seven if r["significant"] == "TRUE"]
holds("no contrast against CLOSED nominally different from zero", sig == [])
tt7 = {(r["scenario"], r["arm"]): abs(float(r["diff"])) / (float(r["sd"]) / math.sqrt(10)) for r in seven}
exact("largest |t|, contrasts against CLOSED", max(tt7.values()), 1.75)
holds("the largest |t| against CLOSED is BRIDGE under matched resources",
      max(tt7, key=tt7.get) == ("resource_matched", "BRIDGE"))
exact("largest |t|, donor contrasts", max(abs(float(r["diff"])) / (float(r["sd"]) / math.sqrt(10))
                                          for r in donor), 1.55)
rb = mc[("baseline", "BRIDGE", "recycling")]
exact("t, recycling at BRIDGE's count", float(rb["diff"]) / (float(rb["sd"]) / math.sqrt(10)), -0.62)
bm = [(tcrit(1 - 0.025 / 7, 9) + tcrit(0.80, 9)) * float(r["sd"]) / math.sqrt(10) for r in seven]
exact("MDE at the Bonferroni level, smallest", min(bm), 1.24)
exact("MDE at the Bonferroni level, largest", max(bm), 2.44)
exact("Bonferroni critical value in estimator_distinct.csv",
      float(t3["estEndpoint"]["bonferroniT"]), 3.462, 3)
# the Bonferroni columns that R/02_analyse_application.R writes to out/application/
pcv = read("out/application/paired_contrasts.csv")
exact("arm-by-scenario contrasts in out/application/paired_contrasts.csv", len(pcv), 16, 0)
exact("distinct contrasts among them", len({round(float(r["d_estEndpoint"]), 9) for r in pcv}), 7, 0)
tb16 = [(float(r["d_estEndpoint"]) - float(r["bonLo_estEndpoint"])) / (float(r["sd_estEndpoint"]) / math.sqrt(10))
        for r in pcv]
holds("every row of that file uses one Bonferroni t (%.4f to %.4f)" % (min(tb16), max(tb16)),
      max(tb16) - min(tb16) < 1e-9)
exact("Bonferroni t in that file (family of 16)", tb16[0], 3.997, 3)
exact("Bonferroni t for a family of 16, recomputed", tcrit(1 - 0.025 / 16, 9), 3.997, 3)

print("\nSection 4, Discussion, and Section 5, Conclusions")
inflow = [mc[(s_, a, "total")] for s_ in ("baseline", "resource_matched") for a in ("DIRECT", "BRIDGE", "MAS")]
exact("inflow contrasts excluding zero, nominal", sum(float(r["hi"]) < 0 for r in inflow), 0, 0)
holds("none after correction", all(float(r["bonHi"]) > 0 for r in inflow))
exact("recycling contrasts excluding zero, nominal", sum(float(r["hi"]) < 0 for r in recyc), 0, 0)
exact("first-order recycling cost per cycle, t/ha",
      2 / 42 * -float(t4[("baseline", "ELITE2")]["gap"]), 0.004, 3)
exact("first-order recycling cost over nineteen rounds, t/ha",
      19 * 2 / 42 * -float(t4[("baseline", "ELITE2")]["gap"]), 0.08)
_, e2lo, e2hi = ci([x - y for x, y in zip(ends("baseline", "ELITE2"), ends("baseline", "CLOSED"))])
holds("the first-order expected loss lies well inside the interval of ELITE2 minus CLOSED",
      e2lo < -19 * 2 / 42 * -float(t4[("baseline", "ELITE2")]["gap"]) < e2hi)
# Section 4.5: the horizon simulated here is that of Allier et al. (60 years)
exact("horizon simulated, years (cycles x years per cycle)", pg["nCycles"] * pg["yearsPerCycle"], 60, 0)
# Table 1: archive, heritability, cycle length and cap
holds("archive 1,000, h2 0.30, 3-year cycles and cap 2 in the engine configuration",
      (pg["archiveSize"], fo["h2"], pg["yearsPerCycle"], pg["cap"]) == (1000, 0.30, 3, 2))
exact("visual sorting cut DIRECT total plots, %", 100 * (1 - pc("visual_sort", "DIRECT") / 1500), 45, 0)
exact("visual sorting cut BRIDGE total plots, %", 100 * (1 - pc("visual_sort", "BRIDGE") / 2400), 28, 0)
zz = st.NormalDist().inv_cdf(0.975) + st.NormalDist().inv_cdf(0.80)
need = [math.ceil((zz * float(r["sd"]) / 0.3) ** 2) for r in donor]
exact("fewest replicates for 0.3 t/ha", min(need), 48, 0)
exact("most replicates for 0.3 t/ha", max(need), 175, 0)
holds("0.3 t/ha is about two-thirds of one cycle of closed gain (%.2f)" % (0.3 / (G / 19)),
      0.6 <= 0.3 / (G / 19) <= 0.7)

# Section 4.3: contrasts each statistic resolves at ten replicates (|t| = sqrt(10) / relative sd)
est7 = [r for r in est if r["scenario"] == "baseline" or
        (r["scenario"] == "resource_matched" and r["arm"] != "ELITE2")]
tB7 = tcrit(1 - 0.025 / 7, 9)
for e_, want in (("estSlope", 0), ("estLateMean", 1), ("estEndpoint", 0)):
    tt = [math.sqrt(10) / float(r["cv"]) for r in est7 if r["estimator"] == e_]
    holds("%s: seven distinct contrasts, |t| from cv equals |mean|/se" % e_,
          len(tt) == 7 and all(abs(math.sqrt(10) / float(r["cv"]) -
                                   abs(float(r["meanDiff"])) / (float(r["sdDiff"]) / math.sqrt(10))) < 1e-9
                               for r in est7 if r["estimator"] == e_))
    exact("%s contrasts different from zero, nominal" % e_, sum(t > tcrit(0.975, 9) for t in tt), want, 0)
    exact("%s contrasts different from zero after Bonferroni over seven" % e_, sum(t > tB7 for t in tt), 0, 0)
z_ = [r for r in est7 if r["estimator"] == "estLateMean" and math.sqrt(10) / float(r["cv"]) > tcrit(0.975, 9)]
holds("the one nominally resolved late-mean contrast is BRIDGE under matched resources",
      [(r["scenario"], r["arm"]) for r in z_] == [("resource_matched", "BRIDGE")])

print("\nSection 4.7, limitations")
exact("QTL spacing on a 2 M chromosome, cM", 200 / g["nQtlPerChr"], 2, 0)
holds("2 cM is more than twenty times 0.09 cM", 2 / 0.09 > 20)
exact("bridging cycles in which the cap binds, baseline %",
      100 * st.mean(r["nAdmitted"] == "2" for r in rows if r["scenario"] == "baseline" and r["strategy"] == "BRIDGE"), 89, 0)
exact("bridging cycles in which the cap binds, matched %",
      100 * st.mean(r["nAdmitted"] == "2" for r in rows if r["scenario"] == "resource_matched" and r["strategy"] == "BRIDGE"), 43, 0)
exact("sh2 wx1 double recessive in a selfed double heterozygote", 1 / 16, 0.0625, 4)

print("\nNumbers inserted from out/numbers.json, and the replicate-4001 rerun")
import json
with tempfile.TemporaryDirectory() as tmp:
    saved = json.load(io.open("out/numbers.json", encoding="utf-8"))
    subprocess.run([sys.executable, "R/06_inference.py"], check=True, stdout=subprocess.DEVNULL)
    again = json.load(io.open("out/numbers.json", encoding="utf-8"))
holds("R/06_inference.py regenerates out/numbers.json exactly (%d numbers)" % len(again),
      saved == again)
NJ = again
for k_, want in (("vr_gain_DIRECT", 1.91), ("vr_gain_BRIDGE", 2.23), ("vr_gain_MAS", 1.86),
                 ("po_gain_DIRECT", 1.95), ("po_gain_BRIDGE", 2.23), ("po_gain_MAS", 2.13),
                 ("rho_gain_DIRECT", 0.504), ("closed_gain20", 8.46), ("mde_min_closed7", 0.901),
                 ("mde_max_closed7", 1.76), ("need_min", 48), ("need_max", 175)):
    exact("numbers.json %s agrees with the independent recomputation" % k_, NJ[k_]["value"], want,
          3 if want < 1 else 2)
# claims printed around the inserted numbers
holds("Section 3.2: all three paired-only gain intervals exclude one",
      all(NJ["po_gain_%s_lo" % a]["value"] > 1 for a in ("DIRECT", "BRIDGE", "MAS")))
holds("Section 3.2: of the genic-variance and accuracy paired-only intervals, only BRIDGE genic excludes one",
      [(o, a) for o in ("genic", "acc") for a in ("DIRECT", "BRIDGE", "MAS")
       if not NJ["po_%s_%s_lo" % (o, a)]["value"] < 1 < NJ["po_%s_%s_hi" % (o, a)]["value"]]
      == [("genic", "BRIDGE")])
holds("Section 3.2: F interval excludes one for gain of DIRECT and BRIDGE, bootstrap interval for BRIDGE only",
      [a for a in ("DIRECT", "BRIDGE", "MAS") if NJ["vr_gain_%s_Flo" % a]["value"] > 1] == ["DIRECT", "BRIDGE"]
      and [a for a in ("DIRECT", "BRIDGE", "MAS") if NJ["vr_gain_%s_blo" % a]["value"] > 1] == ["BRIDGE"])
holds("Section 3.2: Brown-Forsythe rejects for BRIDGE and MAS, not DIRECT",
      [a for a in ("DIRECT", "BRIDGE", "MAS") if NJ["vr_gain_%s_bfp" % a]["value"] < 0.05] == ["BRIDGE", "MAS"])
holds("Section 3.2: after Bonferroni over nine only BRIDGE keeps a bootstrap lower limit above one",
      [a for a in ("DIRECT", "BRIDGE", "MAS") if NJ["vr_gain_%s_b9lo" % a]["value"] > 1] == ["BRIDGE"])
holds("Section 3.2: every genic-variance and accuracy unpaired-to-paired interval spans one",
      all(NJ["vr_%s_%s_Flo" % (o, a)]["value"] < 1 < NJ["vr_%s_%s_Fhi" % (o, a)]["value"] and
          NJ["vr_%s_%s_blo" % (o, a)]["value"] < 1 < NJ["vr_%s_%s_bhi" % (o, a)]["value"]
          for o in ("genic", "acc") for a in ("DIRECT", "BRIDGE", "MAS")))
holds("Section 3.3: slope factor below the constant-increment value, interval includes it and excludes one",
      1 < NJ["ratio_slope_gm_lo"]["value"] < NJ["ratio_slope_gm"]["value"] < NJ["rw_slope"]["value"]
      < NJ["ratio_slope_gm_hi"]["value"])
holds("Section 3.3: observed-increment prediction inside the interval of the geometric mean",
      NJ["ratio_slope_gm_lo"]["value"] < NJ["rw_het"]["value"] < NJ["ratio_slope_gm_hi"]["value"])
exact("Conclusions: slope factor, about 3.7", NJ["ratio_slope_gm"]["value"], 3.7, 1)
exact("Conclusions: Proposition 2 with observed increment variances, 4.1", NJ["rw_het"]["value"], 4.1, 1)
exact("Conclusions: Proposition 2 with constant increments, 4.9", NJ["rw_slope"]["value"], 4.9, 1)
eqb_ = [NJ["eqb_donor6_%s_%s" % (s_, a)]["value"] for s_ in ("b", "m") for a in ("DIRECT", "BRIDGE", "MAS")]
eqbc_ = [NJ["eqbc_donor6_%s_%s" % (s_, a)]["value"] for s_ in ("b", "m") for a in ("DIRECT", "BRIDGE", "MAS")]
holds("Table 7: all corrected and five of six nominal donor equivalence bounds exceed one cycle of gain",
      all(x > NJ["one_cycle_gain"]["value"] for x in eqbc_) and
      sum(x > NJ["one_cycle_gain"]["value"] for x in eqb_) == 5)
holds("Section 3.4.5: both joint bounds exceed twice the closed programme's gain per cycle",
      min(max(eqb_), max(eqbc_)) > 2 * NJ["one_cycle_gain"]["value"])
holds("Proposition 2 constant: 6K(w^2+1)/(5w(w^2-1)) = 4.94 for K = 19, w = 5",
      abs(NJ["rw_slope"]["value"] - 6 * 19 * 26 / (5 * 5 * 24)) < 1e-12)
holds("Proposition 2 constant equals the relative-sd value squared (%.4f)" % (2.22 ** 2),
      abs(NJ["rw_slope"]["value"] - NJ["rw_relsd_slope"]["value"] ** 2) < 1e-9)


print("\nSection 3.5: the published comparison (Peixoto et al. 2024, Tables 2 and 3), recomputed here")
pub_ = read("data/peixoto2024_tables2_3.csv")
holds("192 cells transcribed (2 tables x 2 heritabilities x 2 parameters x 3 pipelines x 8)",
      len(pub_) == 192)
# cells read off the page image, with the printed header and the order of Table 3
for t_, h_, par_, prog_, ge_, printed_, ratio_, m_, s_ in (
        (2, "high", "genetic_mean", "Conv", "absent", "1:3", "3:1", 127.8, 2.8),
        (2, "high", "genetic_mean", "Conv", "absent", "3:1", "1:3", 131.9, 3.1),
        (2, "low", "genetic_variance", "GSF1", "present", "3:1", "1:3", 54.0, 23.9),
        (3, "high", "genetic_mean", "GSTC", "absent", "1:3", "1:3", 152.4, 4.4),
        (3, "low", "genetic_variance", "Conv", "present", "3:1", "3:1", 32.7, 21.0),
        (3, "high", "genetic_mean", "Conv", "absent", "1:1", "1:1", 139.2, 3.0),
        (2, "low", "genetic_mean", "GSTC", "absent", "1:1", "1:1", 114.9, 5.7)):
    c_ = [r for r in pub_ if (int(r["table"]), r["heritability"], r["parameter"], r["program"],
                              r["ge"], r["ratio_printed"]) == (t_, h_, par_, prog_, ge_, printed_)]
    holds("Table %d %s %s %s %s printed %s (read %s) = %.1f (%.1f)"
          % (t_, h_, par_, prog_, ge_, printed_, ratio_, m_, s_),
          len(c_) == 1 and c_[0]["ratio"] == ratio_ and float(c_[0]["mean"]) == m_
          and float(c_[0]["spread"]) == s_)
holds("Table 3 is read as printed and Table 2 with 1:3 and 3:1 exchanged",
      all(r["ratio"] == r["ratio_printed"] for r in pub_ if r["table"] == "3") and
      all(r["ratio"] == {"1:3": "3:1", "3:1": "1:3"}.get(r["ratio_printed"], r["ratio_printed"])
          for r in pub_ if r["table"] == "2"))
vc_ = [(float(r["mean"]), float(r["spread"])) for r in pub_ if r["parameter"] == "genetic_variance"]
exact("variance cells whose spread exceeds the mean (cannot be a standard error)",
      sum(s > m for m, s in vc_), 2, 0)
exact("variance cells whose spread equals the mean", sum(s == m for m, s in vc_), 1, 0)
Z975_, Z80_ = 1.959964, 0.841621
ZB_ = st.NormalDist().inv_cdf(1 - 0.025 / 24)
exact("Bonferroni critical z over 24", ZB_, 3.08)
G_ = {}
for r in pub_:
    G_.setdefault((r["crosses"], r["heritability"], r["parameter"], r["program"], r["ge"]), {})[
        r["ratio"]] = (float(r["mean"]), float(r["spread"]))


def z_(a, b, rho=0.0):
    return (a[0] - b[0]) / math.sqrt((a[1] ** 2 + b[1] ** 2 - 2 * rho * a[1] * b[1]) / 50)


best2 = []
for k_, c_ in G_.items():
    if k_[2] != "genetic_mean":
        continue
    (r1, a1), (r2, a2) = sorted(c_.items(), key=lambda kv: -kv[1][0])[:2]
    best2.append((k_, r1, r2, a1, a2))
holds("24 settings on the genetic mean", len(best2) == 24)
d_ = [a1[0] - a2[0] for k_, r1, r2, a1, a2 in best2]
exact("median best minus second best, units", st.median(d_), 2.2, 1)
exact("... as % of the best mean", st.median([100 * (a1[0] - a2[0]) / a1[0] for k_, r1, r2, a1, a2 in best2]), 1.5, 1)
exact("... in founder additive SDs (sqrt 30)", st.median(d_) / math.sqrt(30), 0.40)
zz_ = [z_(a1, a2) for k_, r1, r2, a1, a2 in best2]
exact("best against second best significant, nominal, of 24", sum(z > Z975_ for z in zz_), 15, 0)
exact("... after Bonferroni over 24", sum(z > ZB_ for z in zz_), 6, 0)
mde_ = [(Z975_ + Z80_) * math.sqrt((a1[1] ** 2 + a2[1] ** 2) / 50) for k_, r1, r2, a1, a2 in best2]
exact("smallest detectable difference at 50 replicates, minimum", min(mde_), 1.8, 1)
exact("smallest detectable difference at 50 replicates, maximum", max(mde_), 3.4, 1)
exact("observed difference below the detectable difference, of 24",
      sum(d < m for d, m in zip(d_, mde_)), 17, 0)
n1_ = sorted(math.ceil((Z975_ + Z80_) ** 2 * (a1[1] ** 2 + a2[1] ** 2) / (0.01 * a1[0]) ** 2)
             for k_, r1, r2, a1, a2 in best2)
exact("replicates for a 1% difference, minimum", n1_[0], 80, 0)
exact("replicates for a 1% difference, maximum", n1_[-1], 437, 0)
exact("replicates for a 1% difference, median", st.median(n1_), 202, 0)
exact("with a correlation of 0.7, significant after correction, of 24",
      sum(z_(a1, a2, 0.7) > ZB_ for k_, r1, r2, a1, a2 in best2), 17, 0)
rec_ = [b for b in best2 if b[0][3] in ("Conv", "GSTC")]
holds("16 settings of Conv and GSTC", len(rec_) == 16)
exact("1:3 ranks first, of 16", sum(r1 == "1:3" for k_, r1, r2, a1, a2 in rec_), 15, 0)
exact("... with 1:1 second", sum(r1 == "1:3" and r2 == "1:1" for k_, r1, r2, a1, a2 in rec_), 15, 0)
exact("median lead of the best over the second, Conv and GSTC, units",
      st.median([a1[0] - a2[0] for k_, r1, r2, a1, a2 in rec_]), 1.5, 1)
exact("... as % of the mean", st.median([100 * (a1[0] - a2[0]) / a1[0] for k_, r1, r2, a1, a2 in rec_]), 1.3, 1)
exact("... significant, nominal, of 16", sum(z_(a1, a2) > Z975_ for k_, r1, r2, a1, a2 in rec_), 8, 0)
exact("... after correction, of 16", sum(z_(a1, a2) > ZB_ for k_, r1, r2, a1, a2 in rec_), 1, 0)
v01_ = [(G_[k_]["1:3"], G_[k_]["0:1"]) for k_, r1, r2, a1, a2 in rec_]
exact("1:3 above 0:1, of 16", sum(a[0] > b[0] for a, b in v01_), 16, 0)
exact("median lead of 1:3 over 0:1, units", st.median([a[0] - b[0] for a, b in v01_]), 4.7, 1)
exact("... in founder additive SDs", st.median([a[0] - b[0] for a, b in v01_]) / math.sqrt(30), 0.86)
exact("... significant after correction, of 16", sum(z_(a, b) > ZB_ for a, b in v01_), 15, 0)
exact("... significant at the nominal level, of 16", sum(z_(a, b) > Z975_ for a, b in v01_), 16, 0)
exact("Bonferroni-significant rankings that are 0:1 leading under GSF1, of 6",
      sum(z_(a1, a2) > ZB_ and k_[3] == "GSF1" and r1 == "0:1" for k_, r1, r2, a1, a2 in best2), 5, 0)
v1311_ = [(G_[k_]["1:3"], G_[k_]["1:1"]) for k_, r1, r2, a1, a2 in rec_]
exact("1:3 above 1:1, of 16", sum(a[0] > b[0] for a, b in v1311_), 15, 0)
exact("1:3 against 1:1 significant after correction with a correlation of 0.7, of 16",
      sum(z_(a, b, 0.7) > ZB_ for a, b in v1311_), 10, 0)
exact("Bonferroni over 16 changes no count: 1:3 against 0:1",
      sum(z_(a, b) > st.NormalDist().inv_cdf(1 - 0.025 / 16) for a, b in v01_), 15, 0)
exact("Bonferroni over 16 changes no count: best against second best",
      sum(z_(a1, a2) > st.NormalDist().inv_cdf(1 - 0.025 / 16) for k_, r1, r2, a1, a2 in rec_), 1, 0)
gv_ = []
for k_, c_ in G_.items():
    if k_[2] == "genetic_variance":
        lo_ = min((v for r_, v in c_.items() if r_ != "0:1"), key=lambda v: v[0])
        gv_.append((lo_, c_["0:1"]))
exact("0:1 lowest genetic variance, of 24", sum(b[0] < a[0] for a, b in gv_), 24, 0)
exact("... drop significant, nominal", sum(z_(a, b) > Z975_ for a, b in gv_), 23, 0)
exact("... drop significant after correction", sum(z_(a, b) > ZB_ for a, b in gv_), 22, 0)
for k_, want in (("pub_res", 15), ("pub_res_bonf", 6), ("pub_below_mde", 17), ("pub_rec_first13", 15),
                 ("pub_rec_res_bonf", 1), ("pub_rec_vs01_res_bonf", 15), ("pub_var_res_bonf", 22),
                 ("pub_var_welch_res_bonf", 21), ("pub_rho07_res_bonf", 17), ("pub_spread_gt_mean", 2)):
    exact("numbers.json %s agrees with the recomputation" % k_, NJ[k_]["value"], want, 0)
holds("Section 3.1: out/rerun_4001/comparison.txt exists", os.path.exists("out/rerun_4001/comparison.txt"))
if os.path.exists("out/rerun_4001/comparison.txt"):
    cmpx = io.open("out/rerun_4001/comparison.txt", encoding="utf-8").read()
    holds("replicate 4001 rerun from the deposited code reproduces every stored record",
          "ALL RECORDS IDENTICAL" in cmpx)
    nre_ = [int(x) for x in re.findall(r"identical field for field\s+(\d+)", cmpx)]
    exact("replicate 4001 rerun: records compared and identical (200 + 120)", sum(nre_), 320, 0)
    rw_ = read("out/rerun_4001/raw_w4001.csv")
    holds("the rerun covers the baseline and resource-matched scenarios",
          {r["scenario"] for r in rw_} == {"baseline", "resource_matched"})

print("\nSection 3.4.4, genic variance at cycle 20 (Figs 3b and 4b)")
ge = [st.mean(float(r["genicEnd"]) for r in per if r["scenario"] == sc and r["strategy"] == a)
      for sc in ("baseline", "resource_matched") for a in ("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS")]
exact("genic variance at cycle 20, lowest arm mean (% of base)", min(ge), 12, 0)
exact("genic variance at cycle 20, highest arm mean (% of base)", max(ge), 20, 0)

print("\n%s" % ("ALL NUMBERS REPRODUCE" if not FAIL else
                "FAILED: %d value(s) do not reproduce: %s" % (len(FAIL), FAIL)))
sys.exit(1 if FAIL else 0)
