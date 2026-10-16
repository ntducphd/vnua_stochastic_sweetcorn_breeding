################################################################################
## 08_figures_3_5.R -- Figures 3 to 5 of the paper, drawn from the pooled output files.
##
## Figures 3 and 4: gain, genic variance, supply of admissible candidates and prediction
## accuracy per cycle, baseline and resource-matched scenarios, under field discard (the
## visual-sorting scenarios are identical in every plotted quantity).
## Figure 5: the decomposition of each inflow arm's contrast with the count-matched
## controls (Table 6), with nominal 95% t intervals.
##
## Style: R/fig_style.R; each arm also has its own line type.
##
## Reads   out/application/pooled_records.csv, out/matched_controls/matched_contrasts.csv
## Writes  figures/R10_fig3_baseline.{png,pdf}, figures/R10_fig4_resource_matched.{png,pdf},
##         figures/R10_fig5_decomposition.{png,pdf}
## Usage:  Rscript R/08_figures_3_5.R
################################################################################

d  <- read.csv("out/application/pooled_records.csv", stringsAsFactors = FALSE)
mc <- read.csv("out/matched_controls/matched_contrasts.csv", stringsAsFactors = FALSE)
dir.create("figures", showWarnings = FALSE)
source("R/fig_style.R")         # FIG_W, FONT, PT, COL, letter(), both()

ARMS <- c("CLOSED", "ELITE2", "DIRECT", "BRIDGE", "MAS")
LTY  <- c(CLOSED = 1, ELITE2 = 2, DIRECT = 3, BRIDGE = 4, MAS = 5)
LWD  <- c(CLOSED = 1.3, ELITE2 = 1.1, DIRECT = 1.1, BRIDGE = 1.1, MAS = 1.1)
tha  <- expression(paste("t ", ha^-1))

summ <- function(z, col) {
  a <- aggregate(z[[col]], by = list(cycle = z$cycle, strategy = z$strategy),
                 FUN = function(x) c(m = mean(x, na.rm = TRUE),
                                     lo = unname(quantile(x, 0.1, na.rm = TRUE)),
                                     hi = unname(quantile(x, 0.9, na.rm = TRUE))))
  data.frame(cycle = a$cycle, strategy = a$strategy,
             m = a$x[, 1], lo = a$x[, 2], hi = a$x[, 3])
}

panel <- function(z, col, ylab, l, arms = ARMS, logy = FALSE, legend = FALSE, ylim = NULL) {
  s <- summ(z[z$strategy %in% arms, ], col)
  tr <- if (logy) function(x) log10(1 + x) else identity
  drawn <- c(s$m, if ("CLOSED" %in% arms) unlist(s[s$strategy == "CLOSED", c("lo", "hi")]))
  yl <- if (is.null(ylim)) range(tr(drawn), na.rm = TRUE) else tr(ylim)
  plot(NA, xlim = c(1, max(z$cycle)), ylim = yl, xlab = "Cycle", ylab = "",
       axes = FALSE, bty = "n")
  axis(1, at = c(1, 5, 10, 15, 20), lwd = AX_LWD, tcl = -0.3)
  if (logy) {
    tk <- c(0, 1, 3, 10, 30, 100, 300)
    tk <- tk[tr(tk) <= yl[2] * 1.02]
    if (max(tk) == 10 && tr(20) <= yl[2] * 1.02) tk <- c(tk, 20)  # top label when 30 is off scale
    axis(2, at = tr(tk), labels = tk, las = 1, lwd = AX_LWD, tcl = -0.3)
  } else axis(2, las = 1, lwd = AX_LWD, tcl = -0.3)
  box(bty = "l", lwd = AX_LWD)
  title(ylab = ylab, line = 2.6)
  if ("CLOSED" %in% arms) {                    # spread of the control only, to keep it legible
    b <- s[s$strategy == "CLOSED", ]
    polygon(c(b$cycle, rev(b$cycle)), tr(c(b$lo, rev(b$hi))),
            col = adjustcolor("#000000", alpha.f = 0.12), border = NA)
  }
  for (a in arms) {
    q <- s[s$strategy == a, ]
    q <- q[order(q$cycle), ]
    lines(q$cycle, tr(q$m), col = COL[[a]], lty = LTY[[a]], lwd = LWD[[a]])
  }
  if (legend) legend("topleft", legend = arms, col = COL[arms], lty = LTY[arms],
                     lwd = LWD[arms], bty = "n", seg.len = 2.6, cex = 1)
  letter(l)
}

figure <- function(sc, stem) {
  z <- d[d$scenario == sc, ]
  both(stem, FIG_W, 120, function() {
    par(mfrow = c(2, 2), mar = c(3.2, 4.6, 1.6, 0.8), mgp = c(2, 0.55, 0), las = 1)
    par(cex = 1)
    panel(z, "gain", expression(paste("Cumulative gain (t ", ha^-1, ")")), "a", legend = TRUE,
          ylim = c(0, 10))
    panel(z, "genicPct", "Genic variance (% of base)", "b", ylim = c(0, 100))
    panel(z, "nClearBar", expression(paste("Candidates clearing, log(1 + ", italic(x), ") axis")), "c",
          arms = setdiff(ARMS, "CLOSED"), logy = TRUE)
    panel(z, "accuracy", "Prediction accuracy, elite candidates", "d", ylim = c(0, 0.7))
  })
}
figure("baseline", "figures/R10_fig3_baseline")
figure("resource_matched", "figures/R10_fig4_resource_matched")

## ---- Figure 5: the decomposition with count-matched controls ----------------
KIND <- c("total: arm - CLOSED", "recycling: ELITE2M - CLOSED", "donor: arm - ELITE2M",
          "unmatched donor: arm - ELITE2")
LAB  <- c("arm − CLOSED", "ELITE2M − CLOSED", "arm − ELITE2M",
          "arm − ELITE2")
PCH  <- c(16, 15, 17, 1)
fig5 <- function(sc, l) {
  rows <- do.call(rbind, lapply(c("DIRECT", "BRIDGE", "MAS"), function(a)
    do.call(rbind, lapply(seq_along(KIND), function(k) {
      r <- mc[mc$scenario == sc & mc$arm == a & mc$contrast == KIND[k], ]
      data.frame(arm = a, k = k, diff = r$diff, lo = r$lo, hi = r$hi)
    }))))
  rows$y <- rev(seq_len(nrow(rows))) + rep(c(2, 1, 0), each = 4)  # gap between arms
  sel <- mc[mc$contrast %in% KIND & mc$arm %in% c("DIRECT", "BRIDGE", "MAS"), ]
  xl <- range(pretty(c(sel$lo, sel$hi)))              # one scale for both panels
  plot(NA, xlim = xl, ylim = c(0.5, max(rows$y) + 0.5), axes = FALSE,
       xlab = expression(paste("Difference in cumulative gain at cycle 20 (t ", ha^-1, ")")),
       ylab = "")
  abline(v = 0, lwd = AX_LWD)
  segments(rows$lo, rows$y, rows$hi, rows$y, lwd = 1, col = COL[rows$arm])
  points(rows$diff, rows$y, pch = PCH[rows$k], col = COL[rows$arm], cex = 0.9)
  axis(1, at = seq(-2, 1.5, 0.5), labels = FALSE, lwd = AX_LWD, tcl = -0.2)
  axisM(1, at = -2:1, lwd = 0, lwd.ticks = AX_LWD, tcl = -0.3)
  mtext(c(baseline = "Baseline", resource_matched = "Resource-matched")[[sc]], side = 3,
        line = 0.35, cex = 1)
  axis(2, at = rows$y, labels = LAB[rows$k], las = 1, lwd = 0, cex.axis = 1)
  for (a in c("DIRECT", "BRIDGE", "MAS")) {
    yy <- rows$y[rows$arm == a]
    mtext(a, side = 4, at = mean(yy), las = 1, line = 0.3, col = COL[[a]], font = 2,
          cex = 1)
  }
  letter(l)
}
both("figures/R10_fig5_decomposition", FIG_W, 92, function() {
  par(mfrow = c(1, 2), mar = c(3.4, 8.8, 1.6, 4.2), mgp = c(2, 0.55, 0))
  par(cex = 1)
  fig5("baseline", "a")
  fig5("resource_matched", "b")
})

cat("wrote figures/R10_fig3_baseline, R10_fig4_resource_matched and",
    "R10_fig5_decomposition (.png and .pdf)\n")
cat(R.version.string, "\n")
