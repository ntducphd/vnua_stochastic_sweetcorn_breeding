################################################################################
## 07_figures_1_2.R -- Figures 1 and 2 of the paper.
##
## Figure 1: what pairing removes (Experiment 2 and the application).
##   a-c  cumulative gain of each arm against the closed programme, 40 paired replicates
##        (filled) and unpaired ones (open), one panel per inflow arm.
##   d    variance-reduction factor of pairing on a log scale: estimate with its F interval
##        and two-sample bootstrap interval, the paired-only estimate of equation (17), and
##        the application's paired-only factors with bootstrap intervals.
## Figure 2: Proposition 2.3 tested on the seven distinct contrasts of the application.
##   a    variance across replicates of the contrast D_t against the number of increments,
##        with the random-walk (linear) and random-drift (quadratic) fits
##   b    variance of the increments by cycle
##   c    drift-scale variance ratio of the late slope to the endpoint, with bootstrap
##        intervals and the predictions of Proposition 2.3
##   d    mean late slope per cycle against mean endpoint / 19 (constant drift: identity)
##
## Style: R/fig_style.R.
##
## Reads   out/exp2/pooled.csv, out/inference/exp2_robust.csv, out/matched_controls/bootstrap_intervals.csv,
##         out/inference/rw_by_cycle_contrast.csv, out/inference/rw_by_cycle.csv, out/inference/rw_contrasts.csv,
##         out/numbers.json (all written by R/06_inference.py or earlier scripts)
## Writes  figures/R02_fig1_pairing.{png,pdf}, figures/R02_fig2_randomwalk.{png,pdf}
## Usage:  Rscript R/07_figures_1_2.R
################################################################################

ex  <- read.csv("out/exp2/pooled.csv", stringsAsFactors = FALSE)
rb  <- read.csv("out/inference/exp2_robust.csv", stringsAsFactors = FALSE)
bs  <- read.csv("out/matched_controls/bootstrap_intervals.csv", stringsAsFactors = FALSE)
rc  <- read.csv("out/inference/rw_by_cycle_contrast.csv", stringsAsFactors = FALSE)
rbc <- read.csv("out/inference/rw_by_cycle.csv", stringsAsFactors = FALSE)
rw  <- read.csv("out/inference/rw_contrasts.csv", stringsAsFactors = FALSE)
N   <- jsonlite::fromJSON("out/numbers.json")
num <- function(k) N[[k]]$value
dir.create("figures", showWarnings = FALSE)
source("R/fig_style.R")         # FIG_W, FONT, PT, COL, letter(), both()
ARMS3 <- c("DIRECT", "BRIDGE", "MAS")
tha   <- function(s) bquote(.(s) ~ "(t" ~ ha^-1 * ")")

## ---- Figure 1 ------------------------------------------------------------------
gainOf <- function(design, arm) {
  z <- ex[ex$design == design & ex$arm == arm, ]
  z <- z[order(z$seedBase, z$rep), ]
  z$gain
}
scatter <- function(arm, l) {
  pA <- gainOf("paired", arm); pC <- gainOf("paired", "CLOSED")
  uA <- gainOf("unpaired", arm); uC <- gainOf("unpaired", "CLOSED")
  r  <- range(ex$gain[ex$arm %in% c("CLOSED", ARMS3)])         # one scale for a-c
  plot(NA, xlim = r, ylim = r, axes = FALSE, xlab = "", ylab = "", asp = 1)
  abline(0, 1, lwd = AX_LWD, col = "grey60")
  points(uC, uA, pch = 1, col = "grey45", cex = 0.8)
  points(pC, pA, pch = 16, col = adjustcolor(COL[[arm]], 0.85), cex = 0.8)
  axis(1, lwd = AX_LWD, tcl = -0.3); axis(2, lwd = AX_LWD, tcl = -0.3, las = 1)
  box(bty = "l", lwd = AX_LWD)
  title(xlab = expression(paste("CLOSED cumulative gain (t ", ha^-1, ")")), line = 2)
  title(ylab = bquote(.(arm) ~ "cumulative gain (t" ~ ha^-1 * ")"), line = 2.4)
  q <- rb[rb$outcome == "gain" & rb$arm == arm, ]
  legend("topleft", bty = "n", cex = 1, pch = c(16, 1),
         col = c(COL[[arm]], "grey45"),
         legend = as.expression(list(
           bquote("paired, " * italic(r) == .(sprintf("%.2f", q$rho)) ~
                  "(" * .(sprintf("%.2f", q$rho_lo)) ~ "to" ~ .(sprintf("%.2f", q$rho_hi)) * ")"),
           "unpaired")))
  letter(l)
}
forest <- function(l) {
  OUT <- c(gain = "Gain", genicPct = "Genic variance", accLate = "Late accuracy")
  e2 <- do.call(rbind, lapply(names(OUT), function(o) rb[rb$outcome == o, ]))
  e2$lab <- paste(OUT[e2$outcome], e2$arm)
  app <- do.call(rbind, lapply(c("baseline", "resource_matched"), function(sc)
    do.call(rbind, lapply(ARMS3, function(a) {
      z <- bs[bs$quantity == sprintf("application gain factor %s %s", sc, a), ]
      data.frame(lab = sprintf("Gain %s, %s", a, if (sc == "baseline") "base" else "matched"),
                 arm = a, est = z$estimate, lo = z$lo, hi = z$hi)
    }))))
  n2 <- nrow(e2); na <- nrow(app)
  y2 <- rev(seq_len(n2)) + na + 1 + c(rep(2, 3), rep(1, 3), rep(0, 3))
  ya <- rev(seq_len(na))
  xl <- c(0.15, 30)
  plot(NA, xlim = log10(xl), ylim = c(0.5, max(y2) + 0.8), axes = FALSE, xlab = "", ylab = "")
  rect(log10(num("F_accept_lo")), min(y2) - 0.6, log10(num("F_accept_hi")), max(y2) + 0.8,
       col = adjustcolor("grey80", 0.5), border = NA)
  abline(v = 0, lwd = AX_LWD)
  tk <- c(0.2, 0.5, 1, 2, 5, 10, 20)
  axis(1, at = log10(tk), labels = tk, lwd = AX_LWD, tcl = -0.3)
  title(xlab = "Variance-reduction factor of pairing (log scale)", line = 2)
  segments(log10(e2$F_lo), y2, log10(e2$F_hi), y2, lwd = 0.7, col = "grey30")
  segments(log10(e2$boot_lo), y2 + 0.22, log10(e2$boot_hi), y2 + 0.22, lwd = 2.2,
           col = adjustcolor(COL[e2$arm], 0.55), lend = 1)
  points(log10(e2$vr), y2, pch = 16, col = COL[e2$arm], cex = 0.9)
  points(log10(e2$paired_only), y2, pch = 5, col = "black", cex = 0.8)
  segments(log10(app$lo), ya, log10(app$hi), ya, lwd = 2.2,
           col = adjustcolor(COL[app$arm], 0.55), lend = 1)
  points(log10(app$est), ya, pch = 5, col = "black", cex = 0.8)
  axis(2, at = c(y2, ya), labels = c(e2$lab, app$lab), las = 1, lwd = 0, cex.axis = 1)
  mtext("Experiment 2, n = 40", side = 4, at = mean(y2), line = 0.4, cex = 1, las = 0)
  mtext("Application, n = 10", side = 4, at = mean(ya), line = 0.4, cex = 1, las = 0)
  legend(log10(3.6), y2[4] + 0.6, bty = "n", cex = 1,     # in the empty space beside rows 4-9
         legend = c("unpaired / paired variance ratio", "F interval", "bootstrap interval",
                    "paired runs only, equation (17)"),
         pch = c(16, NA, NA, 5), lty = c(NA, 1, 1, NA), lwd = c(NA, 0.7, 2.2, NA),
         col = c("grey20", "grey30", adjustcolor("grey30", 0.55), "black"))
  letter(l)
}
both("figures/R02_fig1_pairing", FIG_W, 150, function() {
  layout(matrix(c(1, 2, 3, 4, 4, 4), 2, byrow = TRUE), heights = c(1, 1.25))
  par(mar = c(3.2, 4.2, 1.6, 0.8), mgp = c(2, 0.55, 0), las = 1, cex = 1)
  scatter("DIRECT", "a"); scatter("BRIDGE", "b"); scatter("MAS", "c")
  par(mar = c(3.2, 10.5, 1.6, 7))
  forest("d")
})

## ---- Figure 2 ------------------------------------------------------------------
LTYsc <- c(baseline = 1, resource_matched = 2)
PCHsc <- c(baseline = 16, resource_matched = 1)
key   <- unique(rc[, c("scenario", "arm")])
lab7  <- sprintf("%s%s", key$arm, ifelse(key$scenario == "baseline", "", " (m)"))
panelVar <- function(l) {
  plot(NA, xlim = c(0, 19), ylim = c(0, max(rc$varD) * 1.02), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(nrow(key))) {
    z <- rc[rc$scenario == key$scenario[i] & rc$arm == key$arm[i], ]
    lines(z$t, z$varD, col = adjustcolor(COL[[key$arm[i]]], 0.8), lty = LTYsc[[key$scenario[i]]],
          lwd = 0.8)
  }
  lines(rbc$t, rbc$varD_mean, lwd = 1.6)
  lines(rbc$t, rbc$fit_linear, lwd = 1, lty = 3, col = "grey20")
  lines(rbc$t, rbc$fit_quadratic, lwd = 1, lty = 5, col = "grey20")
  axis(1, at = c(0, 4, 9, 14, 19), lwd = AX_LWD, tcl = -0.3); axis(2, lwd = AX_LWD, tcl = -0.3, las = 1)
  box(bty = "l", lwd = AX_LWD)
  title(xlab = "Increments accumulated (cycle − 1)", line = 2)
  title(ylab = expression(paste("Var(", italic(D)[italic(k)], ") across replicates (t"^2, " ha"^-2, ")")), line = 2.6)
  legend("topleft", bty = "n", cex = 1, lty = c(1, 3, 5), lwd = c(1.6, 1, 1),
         legend = expression("mean of the seven contrasts",
                             paste("random walk, linear in ", italic(k)),
                             paste("varying drift, ", italic(ak) + italic(bk)^2)))
  letter(l)
}
panelInc <- function(l) {
  z <- rc[rc$t >= 1, ]
  plot(NA, xlim = c(2, 20), ylim = c(0, max(z$inc_var) * 1.02), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_len(nrow(key))) {
    q <- z[z$scenario == key$scenario[i] & z$arm == key$arm[i], ]
    lines(q$t + 1, q$inc_var, col = adjustcolor(COL[[key$arm[i]]], 0.6),
          lty = LTYsc[[key$scenario[i]]], lwd = 0.7)
  }
  m <- rbc[rbc$t >= 1, ]
  points(m$t + 1, m$inc_var_mean, pch = 16, cex = 0.8)
  segments(2, num("inc_var_early"), 10, num("inc_var_early"), lwd = 1.6)
  segments(16, num("inc_var_late"), 20, num("inc_var_late"), lwd = 1.6)
  axis(1, at = c(2, 5, 10, 15, 20), lwd = AX_LWD, tcl = -0.3); axis(2, lwd = AX_LWD, tcl = -0.3, las = 1)
  box(bty = "l", lwd = AX_LWD)
  title(xlab = "Cycle", line = 2)
  title(ylab = expression(paste("Variance of the increment (t"^2, " ha"^-2, ")")), line = 2.6)
  letter(l)
}
panelRatio <- function(l) {
  x <- seq_len(nrow(rw)); xg <- nrow(rw) + 1.2
  yl <- c(0.15, 60)
  plot(NA, xlim = c(0.5, xg + 0.5), ylim = log10(yl), axes = FALSE, xlab = "", ylab = "")
  abline(h = log10(num("rw_slope")), lwd = 1, col = "grey20")
  abline(h = log10(num("rw_het")), lwd = 1, lty = 2, col = "grey20")
  abline(h = 0, lwd = AX_LWD, lty = 3)
  segments(x, log10(rw$ratio_slope_lo), x, log10(rw$ratio_slope_hi),
           col = COL[rw$arm], lwd = 1.2)
  points(x, log10(rw$ratio_slope), pch = PCHsc[rw$scenario], col = COL[rw$arm], cex = 1)
  tk <- c(0.2, 0.5, 1, 2, 5, 10, 20)
  axis(2, at = log10(tk), labels = tk, lwd = AX_LWD, tcl = -0.3, las = 1)
  segments(xg, log10(num("ratio_slope_gm_lo")), xg, log10(num("ratio_slope_gm_hi")), lwd = 1.2)
  points(xg, log10(num("ratio_slope_gm")), pch = 18, cex = 1.3)
  axis(1, at = c(x, xg), labels = c(lab7, "Geo. mean"), lwd = AX_LWD, tcl = -0.3, las = 2, cex.axis = 1)
  box(bty = "l", lwd = AX_LWD)
  title(ylab = "Slope / endpoint variance ratio (log scale)", line = 2.6)
  legend("topright", bty = "n", cex = 1, lty = c(1, 2, 3), lwd = c(1, 1, 0.6),
         legend = c(sprintf("constant increments, %.2f", num("rw_slope")),
                    sprintf("observed increment variances, %.2f", num("rw_het")),
                    "equal precision, 1"))
  letter(l)
}
panelMeans <- function(l) {
  r <- range(c(rw$mean_slope_per_cycle, rw$mean_end_over_K, 0))
  plot(NA, xlim = r, ylim = r, axes = FALSE, xlab = "", ylab = "", asp = 1)
  abline(0, 1, lwd = AX_LWD, col = "grey50"); abline(h = 0, v = 0, lwd = 0.4, lty = 3)
  ## legend on a white ground over the reference lines, drawn before the axes
  PCHarm <- c(ELITE2 = 23, DIRECT = 21, BRIDGE = 22, MAS = 24)
  legend("bottomright", bty = "o", bg = "white", box.col = NA, cex = 1,
         pch = c(21, 21, PCHarm), pt.bg = c("grey30", "white", COL[c("ELITE2", ARMS3)]),
         col = c("grey30", "grey30", COL[c("ELITE2", ARMS3)]),
         legend = c("baseline", "resource-matched", "ELITE2", ARMS3))
  points(rw$mean_end_over_K, rw$mean_slope_per_cycle, pch = PCHarm[rw$arm],
         col = COL[rw$arm], bg = ifelse(rw$scenario == "baseline", COL[rw$arm], "white"),
         cex = 1.1)
  axisM(1, lwd = AX_LWD, tcl = -0.3); axisM(2, lwd = AX_LWD, tcl = -0.3, las = 1)
  box(bty = "l", lwd = AX_LWD)
  title(xlab = expression(paste("Mean endpoint contrast / 19 (t ", ha^-1, " per cycle)")), line = 2)
  title(ylab = expression(paste("Mean late slope (t ", ha^-1, " per cycle)")), line = 2.6)
  letter(l)
}
both("figures/R02_fig2_randomwalk", FIG_W, 135, function() {
  par(mfrow = c(2, 2), mar = c(3.4, 4.2, 1.6, 0.8), mgp = c(2, 0.55, 0), las = 1)
  par(cex = 1)
  panelVar("a"); panelInc("b")
  par(mar = c(5.6, 4.2, 1.6, 0.8))
  panelRatio("c")
  par(mar = c(3.4, 4.2, 1.6, 0.8))
  panelMeans("d")
})

cat("wrote figures/R02_fig1_pairing.{png,pdf}, figures/R02_fig2_randomwalk.{png,pdf}\n")
cat(R.version.string, "\n")
