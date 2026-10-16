################################################################################
## 04_analyse_matched_controls.R -- decomposition of each inflow arm's contrast against the
## closed programme, using the count-matched controls (ELITE2M arms) of sim/run_application.R.
##
##   arm - CLOSED              total effect (as in 02_analyse_application.R)
##   ELITE2M_arm - CLOSED      cost of recycling as many extra lines as the arm admitted
##   arm - ELITE2M_arm         value of donor-derived lines over elite ones, at equal
##                             admission counts and equal elite nursery
##   arm - ELITE2, ELITE2 - CLOSED   the decomposition against the unmatched control
##
## Paired over the ten replicates (common founders and seed). Nominal 95% t intervals
## and Bonferroni intervals over the six donor-value contrasts (3 arms x 2 scenarios).
## The visual-sorting scenarios are checked to be genetically identical to baseline and
## resource_matched, so the decomposition is reported for those two.
##
## Usage: Rscript R/04_analyse_matched_controls.R
################################################################################

app <- read.csv("out/application/pooled_records.csv", stringsAsFactors = FALSE)
mf <- Sys.glob("out/matched_controls/raw_m4*.csv")   # the per-seed files, not the pooled output below
if (!length(mf)) stop("no out/matched_controls/raw_m4*.csv; run sim/run_application.R first")
m  <- do.call(rbind, lapply(mf, read.csv, stringsAsFactors = FALSE))
write.csv(m, "out/matched_controls/pooled_records.csv", row.names = FALSE)
nRep <- length(unique(m$rep)); maxC <- max(app$cycle)
cat(sprintf("matched controls: %d rows, %d replicates, arms %s\n", nrow(m), nRep,
            paste(sort(unique(m$strategy)), collapse = ",")))

## ---- visual sorting changes accounting only --------------------------------
chk <- function(a, b) {
  out <- c()
  for (arm in c("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS")) {
    x <- app[app$scenario == a & app$strategy == arm, ]; x <- x[order(x$rep, x$cycle), ]
    y <- app[app$scenario == b & app$strategy == arm, ]; y <- y[order(y$rep, y$cycle), ]
    out[arm] <- max(abs(x$gain - y$gain), abs(x$genicPct - y$genicPct), na.rm = TRUE)
  }
  out
}
idA <- chk("baseline", "visual_sort"); idB <- chk("resource_matched", "visual_matched")
cat("max |difference| in gain and genic variance, baseline vs visual_sort:\n"); print(signif(idA, 3))
cat("resource_matched vs visual_matched:\n"); print(signif(idB, 3))
stopifnot(all(idA < 1e-9), all(idB < 1e-9))

## ---- the replayed schedules match the pooled records -------------------------
## Each ELITE2M arm replays the admissions of its inflow arm, cycle by cycle (sim/run_application.R).
for (sc in c("baseline", "resource_matched")) for (arm in c("DIRECT", "BRIDGE", "MAS")) {
  a <- app[app$scenario == sc & app$strategy == arm, ]; a <- a[order(a$seed, a$cycle), ]
  b <- m[m$scenario == sc & m$strategy == paste0("ELITE2M_", arm), ]; b <- b[order(b$seed, b$cycle), ]
  stopifnot(identical(as.integer(a$nAdmitted), as.integer(b$nAdmitted)), identical(a$seed, b$seed))
}
cat("schedules: every ELITE2M arm admitted exactly what its inflow arm admitted, seed by seed and cycle by cycle\n")

## ---- the seven distinct contrasts of the application -----------------------
## Visual sorting changes plot counts only, and resource matching leaves ELITE2 (which grows
## no introgression nursery) unchanged, so nine of the sixteen arm-by-scenario contrasts of
## 02_analyse_application.R repeat one of the seven below exactly. Table 4 of the paper and the multiplicity family
## are computed over the seven.
est <- read.csv("out/application/estimator_comparison.csv", stringsAsFactors = FALSE)
dup <- function(a, b) { x <- est[est$scenario == a, ]; y <- est[est$scenario == b, ]
  x <- x[order(x$arm, x$estimator), ]; y <- y[order(y$arm, y$estimator), ]; max(abs(x$cv - y$cv)) }
e2 <- est[est$arm == "ELITE2", ]
stopifnot(dup("baseline", "visual_sort") < 1e-12, dup("resource_matched", "visual_matched") < 1e-12,
          max(abs(e2$cv[e2$scenario == "baseline"][order(e2$estimator[e2$scenario == "baseline"])] -
                  e2$cv[e2$scenario == "resource_matched"][order(e2$estimator[e2$scenario == "resource_matched"])])) < 1e-12)
dist <- rbind(est[est$scenario == "baseline", ],
              est[est$scenario == "resource_matched" & est$arm != "ELITE2", ])
wide <- reshape(dist[, c("scenario", "arm", "estimator", "cv")], idvar = c("scenario", "arm"),
                timevar = "estimator", direction = "wide")
names(wide) <- sub("^cv\\.", "", names(wide))
E <- c("estEndpoint", "estLateMean", "estSlope")
least <- table(factor(E[apply(wide[, E], 1, which.min)], levels = E))
most  <- table(factor(E[apply(wide[, E], 1, which.max)], levels = E))
t3 <- data.frame(estimator = E, nContrasts = nrow(wide),
                 medianRelSd = sapply(E, function(x) median(wide[[x]])),
                 minRelSd = sapply(E, function(x) min(wide[[x]])),
                 maxRelSd = sapply(E, function(x) max(wide[[x]])),
                 leastVariable = as.integer(least), mostVariable = as.integer(most))
t3$ratioToEndpoint <- t3$medianRelSd / t3$medianRelSd[1]
t3$bonferroniT <- qt(1 - 0.025 / nrow(wide), length(unique(app$rep)) - 1)
write.csv(t3, "out/matched_controls/estimator_distinct.csv", row.names = FALSE)
cat("\nTable 4 over the seven distinct contrasts:\n"); print(format(t3, digits = 3), row.names = FALSE)

## ---- per-run endpoints -----------------------------------------------------
all <- rbind(app[app$scenario %in% c("baseline", "resource_matched"), names(m)], m)
ep <- do.call(rbind, lapply(split(all, list(all$rep, all$scenario, all$strategy), drop = TRUE),
  function(g) {
    g <- g[order(g$cycle), ]
    data.frame(rep = g$rep[1], scenario = g$scenario[1], strategy = g$strategy[1],
               gainEnd = g$gain[nrow(g)], genicEnd = g$genicPct[nrow(g)],
               admitted = sum(g$nAdmitted), stringsAsFactors = FALSE)
  }))
write.csv(ep, "out/matched_controls/per_run_summary.csv", row.names = FALSE)

pick <- function(sc, a) { z <- ep[ep$scenario == sc & ep$strategy == a, ]; z[order(z$rep), ] }
tNom <- qt(0.975, nRep - 1)
tBon <- qt(1 - 0.025 / 6, nRep - 1)
row <- function(sc, arm, label, a, b) {
  x <- pick(sc, a); y <- pick(sc, b)
  stopifnot(identical(x$rep, y$rep))
  d <- x$gainEnd - y$gainEnd
  se <- sd(d) / sqrt(length(d))
  data.frame(scenario = sc, arm = arm, contrast = label, n = length(d),
             diff = mean(d), sd = sd(d), lo = mean(d) - tNom * se, hi = mean(d) + tNom * se,
             bonLo = mean(d) - tBon * se, bonHi = mean(d) + tBon * se,
             admittedA = mean(x$admitted), admittedB = mean(y$admitted),
             stringsAsFactors = FALSE)
}
tab <- do.call(rbind, lapply(c("baseline", "resource_matched"), function(sc)
  do.call(rbind, lapply(c("DIRECT", "BRIDGE", "MAS"), function(arm) {
    ctl <- paste0("ELITE2M_", arm)
    rbind(row(sc, arm, "total: arm - CLOSED", arm, "CLOSED"),
          row(sc, arm, "recycling: ELITE2M - CLOSED", ctl, "CLOSED"),
          row(sc, arm, "donor: arm - ELITE2M", arm, ctl),
          row(sc, arm, "unmatched donor: arm - ELITE2", arm, "ELITE2"),
          row(sc, arm, "unmatched recycling: ELITE2 - CLOSED", "ELITE2", "CLOSED"))
  }))))
write.csv(tab, "out/matched_controls/matched_contrasts.csv", row.names = FALSE)
op <- options(width = 160)
print(format(tab[, c("scenario", "arm", "contrast", "diff", "lo", "hi", "bonLo", "bonHi",
                     "admittedA", "admittedB")], digits = 3), row.names = FALSE)
options(op)
cat("\nwrote out/matched_controls/matched_contrasts.csv, per_run_summary.csv, pooled_records.csv\n")
cat(R.version.string, "| AlphaSimR", as.character(packageVersion("AlphaSimR")), "\n")
