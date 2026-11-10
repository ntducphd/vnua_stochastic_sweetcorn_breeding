# -*- coding: utf-8 -*-
"""Statistical inference. Writes the tables behind Sections 3.2 to 3.5, the data of Figures 1
and 2, and out/numbers.json, which holds every number the paper quotes from these analyses
(key -> {"value": number, "text": the number as the paper prints it}).

  A  Experiment 2 (out/exp2/pooled.csv)
       variance-reduction factor of pairing with its F interval and a two-sample bootstrap
       interval, Brown-Forsythe test, normality and kurtosis of the contrasts, leave-one-out
       range, Fisher-z intervals for the correlations, and the lower bound on the share of
       variance carried by the shared inputs that Proposition 2.2 gives.
  B  Response statistics (out/application/per_run_summary.csv, out/application/pooled_records.csv)
       Proposition 2.3 on the drift scale: Var(slope per cycle) / Var(endpoint / 19) and
       Var(late mean / 17) / Var(endpoint / 19) for the seven distinct contrasts, with
       bootstrap intervals and Pitman-Morgan tests; increment variances by cycle; the
       heteroscedastic and random-drift predictions; stability of the "largest in six of
       seven" count.
  C  Application contrasts (out/matched_controls/per_run_summary.csv and out/application/)
       exact sign-flip max-T adjustment for each family, Bonferroni as a sensitivity
       analysis, equivalence bounds (the smallest symmetric margin that two one-sided tests
       accept), minimum detectable effects, replicates needed with the chi-square
       uncertainty of the standard deviation, and the remaining quoted quantities.
  D  Published comparison (data/peixoto2024_tables2_3.csv), Section 3.5
       z and Welch tests of the best against the second-best ratio, rounding limits,
       minimum detectable differences and replicates needed.

Bootstrap seed 20261014 (numpy default_rng).

Usage: python R/06_inference.py            (run from the project root)
"""
import csv
import io
import itertools
import json
import math
import os
from collections import OrderedDict

import numpy as np
from scipy import stats

OUTDIR = "out/inference"
os.makedirs(OUTDIR, exist_ok=True)
SEED = 20261014
rng = np.random.default_rng(SEED)
N = OrderedDict()          # key -> {"value": float, "text": str}


def put(key, value, fmt):
    """Record a quoted number with the format in which the manuscript prints it."""
    if isinstance(fmt, str):
        text = fmt % value
    else:
        text = fmt(value)
    if text.startswith("-") and not text.strip("-0."):      # a value that rounds to zero prints as 0.00
        text = text[1:]
    text = text.replace("-", "−") if text.startswith("-") else text
    N[key] = {"value": float(value), "text": text}


def read(path):
    with io.open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write(path, rows):
    keys = list(rows[0].keys())
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, keys, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.6g" % v if isinstance(v, (float, np.floating)) else v)
                        for k, v in r.items()})


ARMS = ("DIRECT", "BRIDGE", "MAS")
OUTCOMES = (("gain", "gain"), ("genicPct", "genic"), ("accLate", "acc"))

# =========================================================================================
# A  Experiment 2
# =========================================================================================
ex2 = read("out/exp2/pooled.csv")


def series(design, arm, field):
    z = sorted((r for r in ex2 if r["design"] == design and r["arm"] == arm),
               key=lambda r: (r["seedBase"], int(r["rep"])))
    return np.array([float(r[field]) for r in z])


F975 = stats.f.ppf(0.975, 39, 39)
F025 = stats.f.ppf(0.025, 39, 39)
B9 = stats.f.ppf(1 - 0.025 / 9, 39, 39)
put("F_accept_lo", F025, "%.2f")
put("F_accept_hi", F975, "%.2f")
B_EX2 = 20000
rowsA = []
for field, tag in OUTCOMES:
    for arm in ARMS:
        pA, pC = series("paired", arm, field), series("paired", "CLOSED", field)
        uA, uC = series("unpaired", arm, field), series("unpaired", "CLOSED", field)
        dp, du = pA - pC, uA - uC
        vr = du.var(ddof=1) / dp.var(ddof=1)
        bi = rng.integers(0, 40, size=(B_EX2, 40))
        bj = rng.integers(0, 40, size=(B_EX2, 40))
        bvr = du[bj].var(axis=1, ddof=1) / dp[bi].var(axis=1, ddof=1)
        blo, bhi = np.percentile(bvr, [2.5, 97.5])
        b9lo = np.percentile(bvr, 100 * 0.025 / 9)
        loo = [np.delete(du, k).var(ddof=1) / dp.var(ddof=1) for k in range(len(du))]
        rho = np.corrcoef(pA, pC)[0, 1]
        zf, hf = math.atanh(rho), 1.96 / math.sqrt(40 - 3)
        po = (pA.var(ddof=1) + pC.var(ddof=1)) / dp.var(ddof=1)
        q = pA.std(ddof=1) / pC.std(ddof=1)
        ii = rng.integers(0, 40, size=(4000, 40))           # paired replicates resampled
        bpo = (pA[ii].var(axis=1, ddof=1) + pC[ii].var(axis=1, ddof=1)) / \
            (pA[ii] - pC[ii]).var(axis=1, ddof=1)
        po_lo, po_hi = np.percentile(bpo, [2.5, 97.5])
        row = OrderedDict(outcome=field, arm=arm, vr=vr, F_lo=vr / F975, F_hi=vr / F025,
                          boot_lo=blo, boot_hi=bhi, F_bonf9_lo=vr / B9, boot_bonf9_lo=b9lo,
                          bf_p=stats.levene(du, dp, center="median").pvalue,
                          shapiro_p_unpaired=stats.shapiro(du).pvalue,
                          shapiro_p_paired=stats.shapiro(dp).pvalue,
                          exkurt_unpaired=stats.kurtosis(du), exkurt_paired=stats.kurtosis(dp),
                          loo_min=min(loo), loo_max=max(loo), rho=rho,
                          rho_lo=math.tanh(zf - hf), rho_hi=math.tanh(zf + hf),
                          paired_only=po, paired_only_lo=po_lo, paired_only_hi=po_hi,
                          sd_ratio_q=q,
                          one_over_1mrho=1 / (1 - rho),
                          exact_from_rho=(1 + q * q) / (1 + q * q - 2 * q * rho))
        rowsA.append(row)
        key = "%s_%s" % (tag, arm)
        put("vr_%s" % key, vr, "%.2f")
        put("vr_%s_Flo" % key, vr / F975, "%.2f")
        put("vr_%s_Fhi" % key, vr / F025, "%.2f")
        put("vr_%s_blo" % key, blo, "%.2f")
        put("vr_%s_bhi" % key, bhi, "%.2f")
        put("vr_%s_b9lo" % key, b9lo, "%.2f")
        put("vr_%s_F9lo" % key, vr / B9, "%.2f")
        put("vr_%s_bfp" % key, row["bf_p"], "%.3f")
        put("vr_%s_loomin" % key, min(loo), "%.2f")
        put("vr_%s_loomax" % key, max(loo), "%.2f")
        put("rho_%s" % key, rho, "%.2f")
        put("rho_%s_lo" % key, row["rho_lo"], "%.2f")
        put("rho_%s_hi" % key, row["rho_hi"], "%.2f")
        put("po_%s" % key, po, "%.2f")
        put("po_%s_lo" % key, po_lo, "%.2f")
        put("po_%s_hi" % key, po_hi, "%.2f")
        put("sw_unpaired_%s" % key, row["shapiro_p_unpaired"], "%.3f")
        put("kurt_unpaired_%s" % key, row["exkurt_unpaired"], "%.2f")
write(os.path.join(OUTDIR, "exp2_robust.csv"), rowsA)
g = [r for r in rowsA if r["outcome"] == "gain"]
put("vr_gain_min", min(r["vr"] for r in g), "%.1f")
put("vr_gain_max", max(r["vr"] for r in g), "%.1f")
put("po_gain_min", min(r["paired_only"] for r in g), "%.1f")
put("po_gain_max", max(r["paired_only"] for r in g), "%.1f")
put("vr_gain_blo_min", min(r["boot_lo"] for r in g), "%.2f")
put("vr_gain_blo_max", max(r["boot_lo"] for r in g), "%.2f")
put("vr_gain_b9lo_min", min(r["boot_bonf9_lo"] for r in g), "%.2f")
put("rho_gain_min", min(r["rho"] for r in g), "%.2f")
put("rho_gain_max", max(r["rho"] for r in g), "%.2f")
put("rho_gain_lo_min", min(r["rho_lo"] for r in g), "%.2f")
put("rho_gain_hi_max", max(r["rho_hi"] for r in g), "%.2f")
# Proposition 2.2: rho <= sqrt(kappa_A kappa_C), so the geometric-mean share of variance carried
# by the shared inputs is at least rho; the lower Fisher limit bounds it with 95% confidence.
put("kappa_lb_min", min(r["rho"] for r in g), "%.2f")
put("kappa_lb_max", max(r["rho"] for r in g), "%.2f")
put("kappa_lb_conf_min", min(r["rho_lo"] for r in g), "%.2f")
# Marginal variance (sum of the two arms' variances) of the unpaired and paired samples,
# DIRECT, gain.
uD, uC = series("unpaired", "DIRECT", "gain"), series("unpaired", "CLOSED", "gain")
pD, pC = series("paired", "DIRECT", "gain"), series("paired", "CLOSED", "gain")
put("margvar_unpaired_DIRECT", uD.var(ddof=1) + uC.var(ddof=1), "%.2f")
put("margvar_paired_DIRECT", pD.var(ddof=1) + pC.var(ddof=1), "%.2f")
# Pooled-marginal estimator: Var(Y_A) + Var(Y_C) from all 80 runs of each arm.
for arm in ARMS:
    A = np.r_[series("paired", arm, "gain"), series("unpaired", arm, "gain")]
    C = np.r_[series("paired", "CLOSED", "gain"), series("unpaired", "CLOSED", "gain")]
    dp = series("paired", arm, "gain") - series("paired", "CLOSED", "gain")
    put("pooledmarg_gain_%s" % arm, (A.var(ddof=1) + C.var(ddof=1)) / dp.var(ddof=1), "%.2f")

# =========================================================================================
# B  Proposition 2.3: response statistics as estimators of the drift
# =========================================================================================
K, W, YRS = 19, 5, 3                      # increments to cycle 20, window, years per cycle


def c_weights(K, w):
    """Coefficients of the OLS slope per step on the increments of a random walk observed at
    steps 0..K, window = the last w points (proof of Proposition 2.3, Section 2.5)."""
    k = np.arange(K - w + 1, K + 1)
    kb, sxx = k.mean(), ((k - k.mean()) ** 2).sum()
    return np.array([sum(kk - kb for kk in k if kk >= j) / sxx for j in range(1, K + 1)])


def late_weights(K, w):
    k = np.arange(K - w + 1, K + 1)
    return np.array([np.mean(k >= j) for j in range(1, K + 1)])


cw = (c_weights(K, W) ** 2).sum()                       # 6(w^2+1)/(5w(w^2-1))
assert abs(cw - 6 * (W * W + 1) / (5 * W * (W * W - 1))) < 1e-12
lw = (late_weights(K, W) ** 2).sum()                    # (K-w+1) + (w-1)(2w-1)/(6w)
assert abs(lw - ((K - W + 1) + (W - 1) * (2 * W - 1) / (6 * W))) < 1e-12
late_mean_steps = np.mean(np.arange(K - W + 1, K + 1))  # 17
rw_slope = K * cw                                        # drift scale: slope / (endpoint/K)
rw_late = (lw / late_mean_steps ** 2) / (1 / K)          # late/17 vs endpoint/19
rw_window = K / (W - 1)                                  # endpoint-of-window rate
put("rw_cw", cw, "%.2f")
put("rw_slope", rw_slope, "%.2f")
put("rw_late", rw_late, "%.2f")
put("rw_window", rw_window, "%.2f")
put("rw_ols_ineff", cw * (W - 1), "%.2f")
# relative SD (sd/mean) of the late mean and the slope under constant drift
put("rw_relsd_late", math.sqrt(lw) / late_mean_steps / (math.sqrt(K) / K), "%.2f")
put("rw_relsd_slope", math.sqrt(cw) / (math.sqrt(K) / K), "%.2f")

per = read("out/application/per_run_summary.csv")
raw = read("out/application/pooled_records.csv")
SEVEN = [("baseline", a) for a in ("ELITE2", "DIRECT", "BRIDGE", "MAS")] + \
        [("resource_matched", a) for a in ARMS]
REPS = list(range(1, 11))
pv = {(r["scenario"], r["strategy"], int(r["rep"]), f): float(r[f]) for r in per
      for f in ("estEndpoint", "estLateMean", "estSlope")}


def contrast(sc, arm, f):
    return np.array([pv[(sc, arm, k, f)] - pv[(sc, "CLOSED", k, f)] for k in REPS])


def traj(sc, arm):
    out = {}
    for r in raw:
        if r["scenario"] == sc and r["strategy"] == arm:
            out.setdefault(int(r["rep"]), {})[int(r["cycle"])] = float(r["gain"])
    return np.array([[out[k][c] for c in range(1, 21)] for k in sorted(out)])


def pitman_morgan(x, y, ratio0=1.0):
    """Test Var(x)/Var(y) = ratio0 for paired samples (Pitman 1939; Morgan 1939)."""
    x = x / math.sqrt(ratio0)
    r = np.corrcoef(x + y, x - y)[0, 1]
    n = len(x)
    t = r * math.sqrt((n - 2) / (1 - r * r))
    return 2 * stats.t.sf(abs(t), n - 2)


B_APP = 4000
rowsB, INC, VD = [], [], []
for sc, arm in SEVEN:
    end = contrast(sc, arm, "estEndpoint") / K                 # drift per cycle
    late = contrast(sc, arm, "estLateMean") / late_mean_steps
    slope = contrast(sc, arm, "estSlope") * YRS                # per year -> per cycle
    idx = rng.integers(0, 10, size=(B_APP, 10))
    rs = slope.var(ddof=1) / end.var(ddof=1)
    rl = late.var(ddof=1) / end.var(ddof=1)
    brs = slope[idx].var(axis=1, ddof=1) / end[idx].var(axis=1, ddof=1)
    brl = late[idx].var(axis=1, ddof=1) / end[idx].var(axis=1, ddof=1)
    D = traj(sc, arm) - traj(sc, "CLOSED")                     # 10 x 20, D[:,0] = 0
    inc = np.diff(D, axis=1)
    ic = inc - inc.mean(axis=0)
    r1 = np.sum(ic[:, 1:] * ic[:, :-1]) / np.sum(ic * ic)
    INC.append(inc.var(axis=0, ddof=1))
    VD.append(D.var(axis=0, ddof=1))
    rowsB.append(OrderedDict(
        scenario=sc, arm=arm, ratio_slope=rs, ratio_slope_lo=np.percentile(brs, 2.5),
        ratio_slope_hi=np.percentile(brs, 97.5),
        pm_p_vs_rw=pitman_morgan(slope, end, rw_slope), pm_p_vs_1=pitman_morgan(slope, end, 1.0),
        ratio_late=rl, ratio_late_lo=np.percentile(brl, 2.5), ratio_late_hi=np.percentile(brl, 97.5),
        mean_slope_per_cycle=slope.mean(), mean_end_over_K=end.mean(),
        lag1_r_increments=r1, var_D20_over_D10=VD[-1][-1] / VD[-1][9]))
INC, VD = np.array(INC), np.array(VD)
s2 = INC.mean(axis=0)                                          # increments 1..19
c = c_weights(K, W)
rw_het = K * K * (c * c * s2).sum() / s2.sum()
put("rw_het", rw_het, "%.2f")
lwt = late_weights(K, W)
put("rw_het_late", (K * K / late_mean_steps ** 2) * (lwt * lwt * s2).sum() / s2.sum(), "%.2f")
put("inc_var_early", INC[:, :9].mean(), "%.3f")
put("inc_var_late", INC[:, 14:].mean(), "%.3f")
tt = np.arange(20).astype(float)
y = VD.mean(axis=0)
lin = np.linalg.lstsq(tt[:, None], y, rcond=None)[0]
quad = np.linalg.lstsq(np.c_[tt, tt ** 2], y, rcond=None)[0]
vdrift = max(quad[1], 0.0)
rw_rd = (vdrift + cw * quad[0]) / (vdrift + quad[0] / K)
put("rw_randomdrift", rw_rd, "%.2f")
rat = np.array([r["ratio_slope"] for r in rowsB])
ratl = np.array([r["ratio_late"] for r in rowsB])
put("ratio_slope_min", rat.min(), "%.2f")
put("ratio_slope_max", rat.max(), "%.2f")
put("ratio_slope_median", np.median(rat), "%.2f")
put("ratio_late_min", ratl.min(), "%.2f")
put("ratio_late_max", ratl.max(), "%.2f")
put("ratio_late_median", np.median(ratl), "%.2f")
put("pm_reject_rw", sum(r["pm_p_vs_rw"] < 0.05 for r in rowsB),
    lambda v: ["no", "one", "two", "three", "four", "five", "six", "seven"][int(v)])
put("pm_reject_1", sum(r["pm_p_vs_1"] < 0.05 for r in rowsB),
    lambda v: ["no", "one", "two", "three", "four", "five", "six", "seven"][int(v)])
put("lag1_min", min(r["lag1_r_increments"] for r in rowsB), "%.2f")
put("lag1_max", max(r["lag1_r_increments"] for r in rowsB), "%.2f")
# geometric mean of the seven drift-scale ratios, bootstrap over replicates
SL = np.array([contrast(sc, a, "estSlope") * YRS for sc, a in SEVEN])
EN = np.array([contrast(sc, a, "estEndpoint") / K for sc, a in SEVEN])
gm = math.exp(np.mean(np.log(SL.var(axis=1, ddof=1) / EN.var(axis=1, ddof=1))))
bgm = []
for b in range(B_APP):
    i = rng.integers(0, 10, 10)
    bgm.append(math.exp(np.mean(np.log(SL[:, i].var(axis=1, ddof=1) / EN[:, i].var(axis=1, ddof=1)))))
put("ratio_slope_gm", gm, "%.2f")
put("ratio_slope_gm_lo", np.percentile(bgm, 2.5), "%.2f")
put("ratio_slope_gm_hi", np.percentile(bgm, 97.5), "%.2f")
# mean check: constant drift implies E[slope per cycle] = E[endpoint]/19
mr = [(r["mean_slope_per_cycle"], r["mean_end_over_K"]) for r in rowsB]
put("meanratio_median", np.median([a / b for a, b in mr]), "%.2f")
put("n_slope_smaller", sum(abs(a) < abs(b) for a, b in mr),
    lambda v: ["no", "one", "two", "three", "four", "five", "six", "seven"][int(v)])
# stability of "the slope has the largest relative sd in >= 6 of 7 contrasts" (Table 4 metric)
REL = {f: np.array([contrast(sc, a, f) for sc, a in SEVEN]) for f in ("estEndpoint", "estLateMean", "estSlope")}
cnt = []
for b in range(B_APP):
    i = rng.integers(0, 10, 10)
    M = np.vstack([REL[f][:, i].std(axis=1, ddof=1) / np.abs(REL[f][:, i].mean(axis=1))
                   for f in ("estEndpoint", "estLateMean", "estSlope")])
    cnt.append(int((M.argmax(axis=0) == 2).sum()))
cnt = np.array(cnt)
put("six_of_seven_boot_pct", 100 * np.mean(cnt >= 6), "%.0f")
put("four_of_seven_boot_pct", 100 * np.mean(cnt >= 4), "%.0f")
write(os.path.join(OUTDIR, "rw_contrasts.csv"), rowsB)
write(os.path.join(OUTDIR, "rw_by_cycle_contrast.csv"),
      [OrderedDict(scenario=sc, arm=a, t=int(t), varD=VD[i, t],
                   inc_var=(INC[i, t - 1] if t >= 1 else float("nan")))
       for i, (sc, a) in enumerate(SEVEN) for t in range(20)])
write(os.path.join(OUTDIR, "rw_by_cycle.csv"),
      [OrderedDict(t=int(t), varD_mean=y[t], varD_min=VD[:, t].min(), varD_max=VD[:, t].max(),
                   fit_linear=lin[0] * t, fit_quadratic=quad[0] * t + quad[1] * t * t,
                   inc_var_mean=(s2[t - 1] if t >= 1 else float("nan")))
       for t in range(20)])

# =========================================================================================
# C  Application contrasts: families, max-T, equivalence bounds, power
# =========================================================================================
pm = read("out/matched_controls/per_run_summary.csv")
ge = {(r["scenario"], r["strategy"], int(r["rep"])): float(r["gainEnd"]) for r in pm}
adm = {(r["scenario"], r["strategy"], int(r["rep"])): float(r["admitted"]) for r in pm}


def ends(sc, arm):
    return np.array([ge[(sc, arm, k)] for k in REPS])


def tstat(D):
    return D.mean(axis=1) / (D.std(axis=1, ddof=1) / math.sqrt(D.shape[1]))


FLIPS = np.array(list(itertools.product((1, -1), repeat=10)))   # 1,024 joint sign flips
t975, t95, t80 = stats.t.ppf(0.975, 9), stats.t.ppf(0.95, 9), stats.t.ppf(0.80, 9)


def family(name, items):
    """items: list of (label, scenario, arm, vector of 10 paired differences)."""
    D = np.array([d for _, _, _, d in items])
    tobs = tstat(D)
    mx = np.array([np.max(np.abs(tstat(D * s))) for s in FLIPS])
    crit = np.quantile(mx, 0.95)
    m = len(items)
    tB = stats.t.ppf(1 - 0.025 / m, 9)
    tBeq = stats.t.ppf(1 - 0.05 / m, 9)
    rows = []
    for (lab, sc, arm, d), t in zip(items, tobs):
        mean, sd = d.mean(), d.std(ddof=1)
        se = sd / math.sqrt(10)
        praw = 2 * stats.t.sf(abs(t), 9)
        rows.append(OrderedDict(
            family=name, contrast=lab, scenario=sc, arm=arm, n=10, mean=mean, sd=sd,
            lo=mean - t975 * se, hi=mean + t975 * se, t=t, p_raw=praw,
            p_bonferroni=min(1.0, m * praw), p_maxT=np.mean(mx >= abs(t) - 1e-12),
            maxT_crit=crit, bonf_crit=tB,
            eqbound_nominal=abs(mean) + t95 * se, eqbound_corrected=abs(mean) + tBeq * se,
            mde80_nominal=(t975 + t80) * se, mde80_maxT=(crit + t80) * se,
            mde80_bonferroni=(tB + t80) * se))
    return rows


C7 = [("arm-closed", sc, a, ends(sc, a) - ends(sc, "CLOSED")) for sc, a in SEVEN]
DON = [("arm-elite2m", sc, a, ends(sc, a) - ends(sc, "ELITE2M_" + a))
       for sc in ("baseline", "resource_matched") for a in ARMS]
REC = [("elite2m-closed", sc, a, ends(sc, "ELITE2M_" + a) - ends(sc, "CLOSED"))
       for sc in ("baseline", "resource_matched") for a in ARMS]
UNM = [("arm-elite2", sc, a, ends(sc, a) - ends(sc, "ELITE2"))
       for sc in ("baseline", "resource_matched") for a in ARMS]
E2E = [("elite2m-elite2", sc, a, ends(sc, "ELITE2M_" + a) - ends(sc, "ELITE2"))
       for sc in ("baseline", "resource_matched") for a in ARMS]
rowsC = family("closed7", C7) + family("donor6", DON) + family("recycling6", REC) + \
        family("unmatched6", UNM) + family("elite2m_vs_elite2_6", E2E)
write(os.path.join(OUTDIR, "contrasts_inference.csv"), rowsC)
F = {(r["family"], r["scenario"], r["arm"]): r for r in rowsC}
ab = {"baseline": "b", "resource_matched": "m"}
for r in rowsC:
    k = "%s_%s_%s" % (r["family"], ab[r["scenario"]], r["arm"])
    put("mean_" + k, r["mean"], "%.2f")
    put("lo_" + k, r["lo"], "%.2f")
    put("hi_" + k, r["hi"], "%.2f")
    put("t_" + k, r["t"], "%.2f")
    put("pmaxT_" + k, r["p_maxT"], "%.3f")
    put("pbonf_" + k, r["p_bonferroni"], "%.3f")
    put("eqb_" + k, r["eqbound_nominal"], "%.2f")
    put("eqbc_" + k, r["eqbound_corrected"], "%.2f")
    put("mde_" + k, r["mde80_nominal"], "%.2f")
for fam in ("closed7", "donor6", "recycling6"):
    z = [r for r in rowsC if r["family"] == fam]
    put("maxTcrit_" + fam, z[0]["maxT_crit"], "%.3f")
    put("bonfcrit_" + fam, z[0]["bonf_crit"], "%.3f")
    put("nsig_maxT_" + fam, sum(r["p_maxT"] < 0.05 for r in z), "%d")
    put("nsig_bonf_" + fam, sum(r["p_bonferroni"] < 0.05 for r in z), "%d")
    put("nsig_raw_" + fam, sum(r["p_raw"] < 0.05 for r in z), "%d")
    put("eqb_min_" + fam, min(r["eqbound_nominal"] for r in z), "%.2f")
    put("eqb_max_" + fam, max(r["eqbound_nominal"] for r in z), "%.2f")
    put("eqbc_min_" + fam, min(r["eqbound_corrected"] for r in z), "%.2f")
    put("eqbc_max_" + fam, max(r["eqbound_corrected"] for r in z), "%.2f")
    put("mde_min_" + fam, min(r["mde80_nominal"] for r in z), "%.2f")
    put("mde_max_" + fam, max(r["mde80_nominal"] for r in z), "%.2f")
    put("mdeT_min_" + fam, min(r["mde80_maxT"] for r in z), "%.2f")
    put("mdeT_max_" + fam, max(r["mde80_maxT"] for r in z), "%.2f")
for fam in ("unmatched6", "elite2m_vs_elite2_6"):
    z = [r for r in rowsC if r["family"] == fam]
    put("mean_min_" + fam, min(r["mean"] for r in z), "%.2f")
    put("mean_max_" + fam, max(r["mean"] for r in z), "%.2f")
don = [r for r in rowsC if r["family"] == "donor6"]
put("donor_mean_min", min(r["mean"] for r in don), "%.2f")
put("donor_mean_max", max(r["mean"] for r in don), "%.2f")
put("donor_maxabs_t", max(abs(r["t"]) for r in don), "%.2f")
G = np.mean(ends("baseline", "CLOSED"))
put("closed_gain20", G, "%.2f")
put("one_cycle_gain", G / K, "%.2f")
put("mde_pct_min_closed7", 100 * N["mde_min_closed7"]["value"] / G, "%.0f")
put("mde_pct_max_closed7", 100 * N["mde_max_closed7"]["value"] / G, "%.0f")
# replicates to detect 0.3 t/ha (80% power, two-sided 0.05, normal approximation), with the
# 95% chi-square interval of each donor contrast's variance on 9 df
zz = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)
need = [(zz * r["sd"] / 0.3) ** 2 for r in don]
lo_f, hi_f = 9 / stats.chi2.ppf(0.975, 9), 9 / stats.chi2.ppf(0.025, 9)
put("need_min", math.ceil(min(need)), "%d")
put("need_max", math.ceil(max(need)), "%d")
put("need_min_lo", math.ceil(min(need) * lo_f), "%d")
put("need_min_hi", math.ceil(min(need) * hi_f), "%d")
put("need_max_lo", math.ceil(max(need) * lo_f), "%d")
put("need_max_hi", math.ceil(max(need) * hi_f), "%d")
# MAS minus DIRECT at baseline
d = ends("baseline", "MAS") - ends("baseline", "DIRECT")
se = d.std(ddof=1) / math.sqrt(10)
put("mas_minus_direct", d.mean(), "%.2f")
put("mas_minus_direct_lo", d.mean() - t975 * se, "%.2f")
put("mas_minus_direct_hi", d.mean() + t975 * se, "%.2f")
# replicate 4004 (rep 4): DIRECT admitted no line yet ended below CLOSED; its count-matched
# control replays CLOSED exactly because the schedule is all zeros
put("rep4004_direct_admitted", adm[("baseline", "DIRECT", 4)], "%d")
put("rep4004_direct_minus_closed", ge[("baseline", "DIRECT", 4)] - ge[("baseline", "CLOSED", 4)], "%.2f")
put("rep4004_e2m_minus_closed", abs(ge[("baseline", "ELITE2M_DIRECT", 4)] - ge[("baseline", "CLOSED", 4)]), "%.1g")
# first-order recycling cost of ELITE2 (two lines 0.08 t/ha below the 40 parents they join)
t4 = {(r["scenario"], r["arm"]): r for r in read("out/application/prediction_bias_tha.csv")}
gap = -float(t4[("baseline", "ELITE2")]["gap"])
put("firstorder_cycle", 2 / 42 * gap, "%.3f")
put("firstorder_19", 19 * 2 / 42 * gap, "%.2f")
put("elite2_vs_firstorder", abs(F[("closed7", "baseline", "ELITE2")]["mean"]) / (19 * 2 / 42 * gap), "%.0f")
# cap binding in BRIDGE and Mendelian recovery
for sc, k in (("baseline", "b"), ("resource_matched", "m")):
    z = [r for r in raw if r["scenario"] == sc and r["strategy"] == "BRIDGE"]
    put("capbind_%s" % k, 100 * np.mean([r["nAdmitted"] == "2" for r in z]), "%.0f")
for arm in ("DIRECT", "MAS"):
    z = [r for r in raw if r["scenario"] == "baseline" and r["strategy"] == arm]
    scr = sum(int(float(r["nF2screened"])) for r in z)
    sw = sum(int(float(r["nSweetUsable"])) for r in z)
    put("binom_p_%s" % arm, stats.binomtest(sw, scr, 0.25).pvalue, "%.3f")
# normal model for the elite pool (Remark 1): 40 of 600 selected, unbiased predictor
accE = np.mean([float(r["accuracy"]) for r in raw if r["scenario"] == "baseline" and r["strategy"] == "ELITE2"])
p_ = 40 / 600
i_ = stats.norm.pdf(stats.norm.ppf(1 - p_)) / p_
put("elite_acc", accE, "%.2f")
put("normal_above_gebv", 600 * stats.norm.sf(i_), "%.0f")
put("normal_above_true", 600 * stats.norm.sf(i_ * accE), "%.0f")

# =========================================================================================
# D  A published comparison: what 50 replicates resolve (Section 3.5)
# =========================================================================================
# Peixoto et al. (2024, G3 14: jkae128), Tables 2 and 3, transcribed into data/ by
# scripts/extract_published_tables.py; the column `ratio` reads both tables in the order of
# Table 3 (Section 2.7). The spreads are read as standard deviations over the 50 replicates
# and the replicates of two scenarios as independent (Sections 2.7 and 4.4).
# Founder additive variance 30 (their GlobalParameters.R, VarG = 30): effects are also given in
# founder additive standard deviations.
pub = read("data/peixoto2024_tables2_3.csv")
NREP_PUB = 50
SD_FOUNDER = math.sqrt(30.0)
cells = {}
for r in pub:
    k = (int(r["crosses"]), r["heritability"], r["parameter"], r["program"], r["ge"])
    cells.setdefault(k, {})[r["ratio"]] = (float(r["mean"]), float(r["spread"]))
z975, z80 = stats.norm.ppf(0.975), stats.norm.ppf(0.80)
zbonf = stats.norm.ppf(1 - 0.025 / 24)                   # 24 settings per parameter


def ztest(a, b, rho=0.0):
    """z for the difference of two scenario means over 50 replicates, spreads as SDs."""
    (m1, s1), (m2, s2) = a, b
    return (m1 - m2) / math.sqrt((s1 ** 2 + s2 ** 2 - 2 * rho * s1 * s2) / NREP_PUB)


def zrange(a, b):
    """Smallest and largest z within the rounding of the published values (to 0.1)."""
    (m1, s1), (m2, s2) = a, b
    lo = max(m1 - m2 - 0.1, 0) / math.sqrt(((s1 + 0.05) ** 2 + (s2 + 0.05) ** 2) / NREP_PUB)
    hi = (m1 - m2 + 0.1) / math.sqrt((max(s1 - 0.05, 0.01) ** 2 + max(s2 - 0.05, 0.01) ** 2) / NREP_PUB)
    return lo, hi


def welch_p(a, b):
    (m1, s1), (m2, s2) = a, b
    v1, v2 = s1 ** 2 / NREP_PUB, s2 ** 2 / NREP_PUB
    df = (v1 + v2) ** 2 / (v1 ** 2 / (NREP_PUB - 1) + v2 ** 2 / (NREP_PUB - 1))
    return 2 * stats.t.sf(abs(m1 - m2) / math.sqrt(v1 + v2), df)


pubrows = []
for k in sorted(cells):
    c = cells[k]
    if k[2] == "genetic_mean":                 # best ratio against the second best
        (r1, a), (r2, b) = sorted(c.items(), key=lambda x: -x[1][0])[:2]
    else:                                      # replacing every parent (0:1) against the lowest other
        b_ = min(((r_, v) for r_, v in c.items() if r_ != "0:1"), key=lambda x: x[1][0])
        (r1, a), (r2, b) = b_, ("0:1", c["0:1"])
    d = a[0] - b[0]
    s2 = a[1] ** 2 + b[1] ** 2
    zlo, zhi = zrange(a, b)
    pubrows.append(dict(
        crosses=k[0], heritability=k[1], parameter=k[2], program=k[3], ge=k[4],
        first=r1, second=r2, mean_first=a[0], mean_second=b[0], diff=d,
        diff_pct=100 * d / a[0], diff_founder_sd=d / SD_FOUNDER,
        z=ztest(a, b), z_rounding_lo=zlo, z_rounding_hi=zhi, welch_p=welch_p(a, b),
        z_rho07=ztest(a, b, 0.7), mde50=(z975 + z80) * math.sqrt(s2 / NREP_PUB),
        n_1pct=math.ceil((z975 + z80) ** 2 * s2 / (0.01 * a[0]) ** 2),
        z_if_se=d / math.sqrt(s2)))
write(OUTDIR + "/published_reanalysis.csv", pubrows)
gm = [r for r in pubrows if r["parameter"] == "genetic_mean"]
gv = [r for r in pubrows if r["parameter"] == "genetic_variance"]
assert len(gm) == 24 and len(gv) == 24
var_cells = [v for k, c in cells.items() if k[2] == "genetic_variance" for v in c.values()]
put("pub_n", len(gm), "%d")
put("pub_bonf_z", zbonf, "%.3f")
put("pub_founder_sd", SD_FOUNDER, "%.1f")
put("pub_spread_gt_mean", sum(s > m for m, s in var_cells), "%d")
put("pub_spread_eq_mean", sum(s == m for m, s in var_cells), "%d")
# all 24 settings: best against second best
put("pub_d_median", np.median([r["diff"] for r in gm]), "%.1f")
put("pub_dpct_median", np.median([r["diff_pct"] for r in gm]), "%.1f")
put("pub_d_median_sd", np.median([r["diff_founder_sd"] for r in gm]), "%.2f")
put("pub_res", sum(r["z"] > z975 for r in gm), "%d")
put("pub_res_bonf", sum(r["z"] > zbonf for r in gm), "%d")
put("pub_res_round_min", sum(r["z_rounding_lo"] > z975 for r in gm), "%d")
put("pub_res_round_max", sum(r["z_rounding_hi"] > z975 for r in gm), "%d")
put("pub_resb_round_min", sum(r["z_rounding_lo"] > zbonf for r in gm), "%d")
put("pub_resb_round_max", sum(r["z_rounding_hi"] > zbonf for r in gm), "%d")
put("pub_welch_res", sum(r["welch_p"] < 0.05 for r in gm), "%d")
put("pub_welch_res_bonf", sum(r["welch_p"] < 0.05 / 24 for r in gm), "%d")
put("pub_mde_min", min(r["mde50"] for r in gm), "%.1f")
put("pub_mde_max", max(r["mde50"] for r in gm), "%.1f")
put("pub_mde_sd_min", min(r["mde50"] for r in gm) / SD_FOUNDER, "%.2f")
put("pub_mde_sd_max", max(r["mde50"] for r in gm) / SD_FOUNDER, "%.2f")
put("pub_below_mde", sum(r["diff"] < r["mde50"] for r in gm), "%d")
put("pub_n1pct_min", min(r["n_1pct"] for r in gm), "%d")
put("pub_n1pct_max", max(r["n_1pct"] for r in gm), "%d")
put("pub_n1pct_median", np.median([r["n_1pct"] for r in gm]), "%.0f")
put("pub_rho07_res_bonf", sum(r["z_rho07"] > zbonf for r in gm), "%d")
# the 16 settings of Conv and GSTC, behind the recommendation to keep a quarter of the old parents
rec = [r for r in gm if r["program"] in ("Conv", "GSTC")]
assert len(rec) == 16
put("pub_rec_n", len(rec), "%d")
put("pub_rec_first13", sum(r["first"] == "1:3" for r in rec), "%d")
put("pub_rec_second11", sum(r["first"] == "1:3" and r["second"] == "1:1" for r in rec), "%d")
put("pub_rec_d_median", np.median([r["diff"] for r in rec]), "%.1f")
put("pub_rec_dpct_median", np.median([r["diff_pct"] for r in rec]), "%.1f")
put("pub_rec_d_median_sd", np.median([r["diff_founder_sd"] for r in rec]), "%.2f")
put("pub_rec_res", sum(r["z"] > z975 for r in rec), "%d")
put("pub_rec_res_bonf", sum(r["z"] > zbonf for r in rec), "%d")
vs01 = []
for r in rec:
    c = cells[(r["crosses"], r["heritability"], "genetic_mean", r["program"], r["ge"])]
    vs01.append((c["1:3"][0] - c["0:1"][0], ztest(c["1:3"], c["0:1"])))
put("pub_rec_vs01_res", sum(z > z975 for d, z in vs01), "%d")
put("pub_res_bonf_gsf1", sum(r["z"] > zbonf and r["program"] == "GSF1" and r["first"] == "0:1" for r in gm), "%d")
v1311 = [(cells[(r["crosses"], r["heritability"], "genetic_mean", r["program"], r["ge"])]["1:3"],
          cells[(r["crosses"], r["heritability"], "genetic_mean", r["program"], r["ge"])]["1:1"]) for r in rec]
put("pub_rec_13v11_pos", sum(a[0] > b[0] for a, b in v1311), "%d")
put("pub_rho07_rec13v11_res_bonf", sum(ztest(a, b, 0.7) > zbonf for a, b in v1311), "%d")
put("pub_rec_vs01_lead", sum(d > 0 for d, z in vs01), "%d")
put("pub_rec_vs01_d_median", np.median([d for d, z in vs01]), "%.1f")
put("pub_rec_vs01_d_median_sd", np.median([d for d, z in vs01]) / SD_FOUNDER, "%.2f")
put("pub_rec_vs01_res_bonf", sum(z > zbonf for d, z in vs01), "%d")
# genetic variance: replacing every parent against the lowest other ratio
put("pub_var_lowest", sum(r["second"] == "0:1" and r["diff"] > 0 for r in gv), "%d")
put("pub_var_res", sum(r["z"] > z975 for r in gv), "%d")
put("pub_var_res_bonf", sum(r["z"] > zbonf for r in gv), "%d")
put("pub_var_welch_res_bonf", sum(r["welch_p"] < 0.05 / 24 for r in gv), "%d")

with io.open("out/numbers.json", "w", encoding="utf-8") as fh:
    json.dump(N, fh, indent=1, ensure_ascii=False)
print("wrote %s/exp2_robust.csv, rw_contrasts.csv, rw_by_cycle.csv, contrasts_inference.csv"
      % OUTDIR)
print("wrote out/numbers.json (%d numbers; seed %d)" % (len(N), SEED))
for k in ("vr_gain_DIRECT", "vr_gain_DIRECT_blo", "vr_gain_DIRECT_bhi", "po_gain_min", "po_gain_max",
          "rw_slope", "rw_het", "rw_randomdrift", "ratio_slope_median", "ratio_slope_gm",
          "six_of_seven_boot_pct", "maxTcrit_closed7", "nsig_maxT_closed7",
          "eqb_min_donor6", "eqb_max_donor6", "need_min", "need_max"):
    print("  %-28s %s" % (k, N[k]["text"]))
