################################################################################
## 03_founder_scale.R -- the yield scale of each deposited founder object, and
## the prediction-bias columns of Table 5 of the paper on that scale.
##
## The yield trait is simulated on AlphaSimR's internal genetic scale and rescaled to
## t/ha per replicate (makeFounders: toYield(g) = g * scale + shift). Gains are written
## in t/ha by runProgram, but the prediction bias inside the candidate pool (biasDonor)
## and the admitted lines' distance from the true bar (admitTrueGap) are recorded on the
## internal scale. Their t/ha values are the internal values times that replicate's
## scale, so the conversion is made replicate by replicate here.
##
## Reads   out/founders/founder_<seed>.rds  (the deposited founder objects)
##         out/application/per_run_summary.csv  (written by R/02_analyse_application.R)
## Writes  out/application/founder_scale.csv    one row per replicate
##         out/application/prediction_bias_tha.csv  Table 5 columns in t/ha, mean and 95% t
##                                                  interval over replicates
##
## Usage: Rscript R/03_founder_scale.R
################################################################################

suppressMessages(library(AlphaSimR))

files <- sort(Sys.glob("out/founders/founder_4*.rds"))
if (!length(files)) stop("no founder objects under out/founders/")

sc <- do.call(rbind, lapply(files, function(path) {
  f <- readRDS(path)
  data.frame(seed      = f$seed,
             scale     = 1 / f$scale$sdA0,          # t/ha per internal unit
             sdA_tha   = f$base$sdA,                # additive SD of the cycle-1 nursery, t/ha
             h2_cycle1 = f$base$h2Cycle1,           # heritability in that nursery
             varE      = f$scale$varE,              # residual variance, internal scale
             locusChr  = f$locus$chr,
             nRedraw   = f$nRedraw,
             baseMean  = f$base$meanG,
             donorMean = f$base$donorMean,
             gap       = f$base$realisedGap,
             AlphaSimR = f$software$AlphaSimR,
             R         = f$software$RString)
}))
write.csv(sc, "out/application/founder_scale.csv", row.names = FALSE)

pr <- read.csv("out/application/per_run_summary.csv", stringsAsFactors = FALSE)
pr$scale <- sc$scale[match(pr$rep, seq_len(nrow(sc)))]   # rep r ran on the r-th seed
stopifnot(!anyNA(pr$scale))

ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x); m <- mean(x)
  h <- if (n > 1) qt(0.975, n - 1) * sd(x) / sqrt(n) else NA_real_
  c(n = n, mean = m, lo = m - h, hi = m + h)
}
tab <- do.call(rbind, lapply(c("baseline", "resource_matched"), function(s)
  do.call(rbind, lapply(c("DIRECT", "MAS", "BRIDGE", "ELITE2"), function(a) {
    z <- pr[pr$scenario == s & pr$strategy == a, ]
    z <- z[order(z$rep), ]
    b <- ci(z$biasDonor * z$scale)
    g <- ci(z$admitTrueGap * z$scale)
    data.frame(scenario = s, arm = a,
               biasN = b[["n"]], bias = b[["mean"]], biasLo = b[["lo"]], biasHi = b[["hi"]],
               gapN = g[["n"]], gap = g[["mean"]], gapLo = g[["lo"]], gapHi = g[["hi"]],
               clearBar = mean(z$clearBarPerCyc), clearTrue = mean(z$clearTruePerCyc),
               stringsAsFactors = FALSE)
  }))))
write.csv(tab, "out/application/prediction_bias_tha.csv", row.names = FALSE)

print(format(sc[, c("seed", "scale", "sdA_tha", "h2_cycle1", "locusChr", "nRedraw")], digits = 4),
      row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)
cat("wrote out/application/founder_scale.csv and out/application/prediction_bias_tha.csv\n")
cat(R.version.string, "| AlphaSimR", as.character(packageVersion("AlphaSimR")), "\n")
