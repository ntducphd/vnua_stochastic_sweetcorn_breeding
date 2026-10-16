################################################################################
## 02_analyse_application.R -- analysis of the application.
##
## Writes the founder summary, per-run summaries, paired contrasts against CLOSED with
## Bonferroni intervals, minimum detectable effects at 80% power, the comparison of three
## late-horizon response statistics, the cost per admitted line (ratio estimator with a
## bootstrap interval) and the supply of admissible lines with prediction accuracy and
## bias inside the donor-derived pool.
##
## Usage: Rscript R/02_analyse_application.R [infile] [primaryOutcome]
## The default input is written by R/01_pool_application.py.  The multiplicity family printed
## here counts all 16 arm-by-scenario contrasts; nine of them repeat one of the seven
## distinct ones exactly; the paper uses the seven in R/04_analyse_matched_controls.R (Table 4) and
## R/06_inference.py (adjusted p-values, family closed7).
################################################################################

args   <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "out/application/pooled_records.csv"
PRIM   <- if (length(args) >= 2) args[2] else "gain"

d <- read.csv(infile, stringsAsFactors = FALSE)
dir.create("out/application", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE)

ARMS <- intersect(c("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS"), unique(d$strategy))
d$strategy <- factor(d$strategy, levels = ARMS)
SCEN <- unique(d$scenario)
maxC <- max(d$cycle)
nRep <- length(unique(d$rep))
cat(sprintf("replicates %d | cycles %d | arms %s | scenarios %s\n",
            nRep, maxC, paste(ARMS, collapse = ","), paste(SCEN, collapse = ",")))

## -----------------------------------------------------------------------------
## Founder base, including how often the endosperm locus had to be redrawn to keep
## it out of the QTL and marker sets.
## -----------------------------------------------------------------------------
fs <- d[!duplicated(d$rep), c("rep", "seed", "locusChr", "maxFreqDiff", "cullDonor",
                              "baseMean", "donorMean", "realisedGap", "nRedraw")]
cat("\n--- founder base ---\n"); print(format(fs, digits = 4), row.names = FALSE)
write.csv(fs, "out/application/founder_summary.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Per-run summary.  Three estimators of the late-horizon rate are kept side by side.
## -----------------------------------------------------------------------------
late <- function(g) g$cycle > (maxC - 5)
perRun <- do.call(rbind, lapply(split(d, list(d$rep, d$scenario, d$strategy), drop = TRUE),
  function(g) {
    g <- g[order(g$cycle), ]; L <- late(g)
    data.frame(
      rep = g$rep[1], scenario = g$scenario[1], strategy = g$strategy[1],
      ## three response statistics: two levels and a rate
      estEndpoint = g$gain[nrow(g)],
      estLateMean = mean(g$gain[L]),
      estSlope    = unname(coef(lm(gain ~ year, data = g[L, ]))[2]),
      genicEnd    = g$genicPct[nrow(g)],
      genicLate   = mean(g$genicPct[L]),
      accEnd      = g$accuracy[nrow(g)],
      accLate     = mean(g$accuracy[L], na.rm = TRUE),
      donorEnd    = g$donorPct[nrow(g)],
      ## uncensored supply and the predictor inside the donor pool
      candPerCyc     = mean(g$nCand),
      clearBarPerCyc = mean(g$nClearBar),
      clearTruePerCyc= mean(g$nClearTrue),
      accDonor       = mean(g$accDonor, na.rm = TRUE),
      biasDonor      = mean(g$biasDonor, na.rm = TRUE),
      admitTrueGap   = mean(g$admitTrueGap, na.rm = TRUE),
      admitted       = sum(g$nAdmitted),
      plotsPerCyc    = mean(g$plotsGrown),
      f2screened     = mean(g$nF2screened),
      sweetUsable    = mean(g$nSweetUsable),
      gapStart = g$gapCycle[1], gapEnd = g$gapCycle[nrow(g)],
      visualSort = g$visualSort[1])
  }))
write.csv(perRun, "out/application/per_run_summary.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Paired contrasts against CLOSED, with an explicit multiplicity family.
## The family is every arm x scenario contrast on the PRIMARY outcome.
## -----------------------------------------------------------------------------
pick <- function(sc, a, f) {
  z <- perRun[perRun$scenario == sc & perRun$strategy == a, ]
  z[[f]][order(z$rep)]
}
contrasts <- setdiff(ARMS, "CLOSED")
mFam <- length(contrasts) * length(SCEN)
tNom <- qt(0.975, nRep - 1)
tBon <- qt(1 - 0.025 / mFam, nRep - 1)
cat(sprintf("\nmultiplicity family: %d contrasts (%d arms x %d scenarios) on '%s'; t* nominal %.3f, Bonferroni %.3f\n",
            mFam, length(contrasts), length(SCEN), PRIM, tNom, tBon))

pairTab <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(contrasts, function(a) {
  out <- list(scenario = sc, arm = a, n = nRep)
  for (f in c("estEndpoint", "estLateMean", "estSlope", "genicEnd", "genicLate",
              "accLate", "donorEnd")) {
    cv <- pick(sc, "CLOSED", f); av <- pick(sc, a, f)
    dd <- av - cv
    dd <- dd[is.finite(dd)]
    if (length(dd) < 2) { out[[paste0("d_", f)]] <- NA; next }
    se <- sd(dd) / sqrt(length(dd))
    out[[paste0("d_", f)]]      <- mean(dd)
    out[[paste0("lo_", f)]]     <- mean(dd) - tNom * se
    out[[paste0("hi_", f)]]     <- mean(dd) + tNom * se
    out[[paste0("bonLo_", f)]]  <- mean(dd) - tBon * se
    out[[paste0("bonHi_", f)]]  <- mean(dd) + tBon * se
    out[[paste0("sd_", f)]]     <- sd(dd)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}))))
write.csv(pairTab, "out/application/paired_contrasts.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Minimum detectable effect of each contrast on the primary outcome, and whether it
## differs from zero at the nominal level.
## -----------------------------------------------------------------------------
primField <- switch(PRIM, gain = "estEndpoint", slope = "estSlope", genic = "genicEnd",
                    "estEndpoint")
cat(sprintf("\n--- primary outcome '%s': power ---\n", PRIM))
cat(sprintf("%-17s %-7s %8s %8s %9s %9s\n", "scenario", "arm", "diff", "sd", "MDE80", "sig?"))
pw <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(contrasts, function(a) {
  dd <- pick(sc, a, primField) - pick(sc, "CLOSED", primField)
  dd <- dd[is.finite(dd)]
  n <- length(dd); if (n < 3) return(NULL)
  m <- mean(dd); s <- sd(dd); se <- s / sqrt(n)
  ## smallest true effect this design detects with 80% power at alpha = 0.05
  mde <- (qt(0.975, n - 1) + qt(0.80, n - 1)) * se
  sig <- abs(m) > qt(0.975, n - 1) * se
  cat(sprintf("%-17s %-7s %8.3f %8.3f %9.3f %9s\n", sc, a, m, s, mde, ifelse(sig, "yes", "no")))
  data.frame(scenario = sc, arm = a, diff = m, sd = s, mde80 = mde, significant = sig,
             stringsAsFactors = FALSE)
}))))
write.csv(pw, "out/application/power.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Estimator comparison for the late-horizon rate.
## The three statistics measure different quantities (two levels and a rate), so they
## are compared on RELATIVE sampling variability: the standard deviation of the paired
## contrast divided by the absolute mean contrast, which sets the replicates needed.
## -----------------------------------------------------------------------------
cat("\n--- estimators of the late-horizon response, relative sampling variability ---\n")
cat(sprintf("%-17s %-7s %-12s %10s %10s %10s\n", "scenario", "arm", "estimator",
            "mean diff", "sd diff", "|sd/mean|"))
estTab <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(contrasts, function(a) {
  do.call(rbind, lapply(c("estEndpoint", "estLateMean", "estSlope"), function(f) {
    dd <- pick(sc, a, f) - pick(sc, "CLOSED", f); dd <- dd[is.finite(dd)]
    if (length(dd) < 2) return(NULL)
    cv <- sd(dd) / abs(mean(dd))
    cat(sprintf("%-17s %-7s %-12s %10.4f %10.4f %10.2f\n", sc, a, f,
                mean(dd), sd(dd), cv))
    data.frame(scenario = sc, arm = a, estimator = f, meanDiff = mean(dd),
               sdDiff = sd(dd), cv = cv, stringsAsFactors = FALSE)
  }))
}))))
write.csv(estTab, "out/application/estimator_comparison.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## Cost of the recessive constraint per arm and scenario: plots per admitted line with a
## bootstrap interval, and replicates that admitted nothing.
## -----------------------------------------------------------------------------
bootRatio <- function(num, den, B = 4000) {
  keep <- is.finite(num) & is.finite(den)
  num <- num[keep]; den <- den[keep]
  if (!length(num) || sum(den) == 0) return(c(NA, NA, NA))
  est <- sum(num) / sum(den)
  bs <- replicate(B, { i <- sample.int(length(num), replace = TRUE)
                       if (sum(den[i]) == 0) NA else sum(num[i]) / sum(den[i]) })
  c(est, quantile(bs, 0.025, na.rm = TRUE), quantile(bs, 0.975, na.rm = TRUE))
}
set.seed(20260925)   # bootstrap seed
cat("\n--- cost of the recessive constraint (ratio estimator, bootstrap 95% CI) ---\n")
cat(sprintf("%-17s %-7s %9s %9s %8s %9s %7s %9s %s\n", "scenario", "arm", "F2/cyc",
            "sweet/cyc", "usable", "plots/cyc", "admit", "plots/line", "zero-admit reps"))
costTab <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(ARMS, function(a) {
  z <- perRun[perRun$scenario == sc & perRun$strategy == a, ]
  if (!nrow(z)) return(NULL)
  zeroes <- sum(z$admitted == 0)
  r <- bootRatio(z$plotsPerCyc * maxC, z$admitted)
  cat(sprintf("%-17s %-7s %9.0f %9.1f %8.3f %9.0f %7.1f %9s %d of %d\n",
              sc, a, mean(z$f2screened), mean(z$sweetUsable),
              ifelse(mean(z$f2screened) > 0, mean(z$sweetUsable) / mean(z$f2screened), NA),
              mean(z$plotsPerCyc), mean(z$admitted),
              ifelse(is.na(r[1]), "-", sprintf("%.0f", r[1])), zeroes, nrow(z)))
  data.frame(scenario = sc, arm = a, f2screened = mean(z$f2screened),
             sweetUsable = mean(z$sweetUsable), plotsPerCyc = mean(z$plotsPerCyc),
             admitted = mean(z$admitted), plotsPerLine = r[1],
             plotsPerLineLo = r[2], plotsPerLineHi = r[3],
             zeroAdmitReps = zeroes, nReps = nrow(z), stringsAsFactors = FALSE)
}))))
write.csv(costTab, "out/application/locus_cost.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## The supply of admissible lines, uncensored, and the predictor inside the donor pool.
## -----------------------------------------------------------------------------
cat("\n--- supply of admissible lines, and the predictor where the decision is made ---\n")
cat(sprintf("%-17s %-7s %9s %11s %12s %9s %9s %11s\n", "scenario", "arm", "cand/cyc",
            "clearBar/cyc", "clearTrue/cyc", "accDonor", "biasDonor", "admit-trueBar"))
supTab <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(setdiff(ARMS, "CLOSED"),
  function(a) {
    z <- perRun[perRun$scenario == sc & perRun$strategy == a, ]
    if (!nrow(z)) return(NULL)
    cat(sprintf("%-17s %-7s %9.0f %11.1f %12.1f %9.3f %9.3f %11.3f\n", sc, a,
                mean(z$candPerCyc), mean(z$clearBarPerCyc), mean(z$clearTruePerCyc),
                mean(z$accDonor, na.rm = TRUE), mean(z$biasDonor, na.rm = TRUE),
                mean(z$admitTrueGap, na.rm = TRUE)))
    data.frame(scenario = sc, arm = a, cand = mean(z$candPerCyc),
               clearBar = mean(z$clearBarPerCyc), clearTrue = mean(z$clearTruePerCyc),
               accDonor = mean(z$accDonor, na.rm = TRUE),
               biasDonor = mean(z$biasDonor, na.rm = TRUE),
               admitTrueGap = mean(z$admitTrueGap, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))))
write.csv(supTab, "out/application/supply.csv", row.names = FALSE)


cat("\nwrote out/application/*.csv\n")
cat(R.version.string, "| AlphaSimR", as.character(packageVersion("AlphaSimR")),
    "| rrBLUP", as.character(packageVersion("rrBLUP")), "|", R.version$platform, "\n")
