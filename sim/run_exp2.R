#!/usr/bin/env Rscript
## =============================================================================
## run_exp2.R -- Experiment 2: variance removed by pairing on common founders
##
## One shard of the reduced configuration (Table 1, Experiment 2 column):
##   paired    replicate r: one founder population and one programme seed, base + r,
##             shared by CLOSED, DIRECT, BRIDGE and MAS
##   unpaired  replicate r, arm k (k = 1..4 for CLOSED, DIRECT, BRIDGE, MAS): its own
##             founder population and programme seed, base + 100000 k + r
## Founder populations are drawn with runMacs and not saved (Section 2.4).
## Writes <out>/shard<k>.csv with columns design, seedBase, rep, arm, gain (cumulative
## gain at the last cycle, t/ha), genicPct (genic variance at the last cycle, % of
## post-burn-in), accLate (mean accuracy over the last five cycles) and seed, the layout
## R/exp2_pool.py reads. Rows are appended one (design, replicate) block at a time; a
## rerun skips the blocks already in the file.
##
## Usage:
##   Rscript sim/run_exp2.R --shard 0 --out out_sim/exp2 [--reps 10] [--config exp2]
##       [--designs paired,unpaired]
## Shards 0-3 use base seeds 500000, 510000, 520000, 530000.
## =============================================================================

Sys.setenv(OMP_NUM_THREADS = "1")
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  simDir <- if (length(f)) dirname(normalizePath(f[1])) else "sim"
  source(file.path(simDir, "R", "load.R"))
  load_sim(simDir)
})

opt <- parse_cli(commandArgs(trailingOnly = TRUE))
if (is.null(opt$shard) || is.null(opt$out)) {
  cat("usage: Rscript sim/run_exp2.R --shard <0-3> --out <dir> [--reps 10] [--designs paired,unpaired]\n")
  quit(status = 2)
}
cfg <- get_config(if (is.null(opt$config)) "exp2" else opt$config)
shard <- as.integer(opt$shard)
base <- cfg$seeds$shardBases[shard + 1L]
if (is.na(base)) stop("shard must be between 0 and ", length(cfg$seeds$shardBases) - 1L)
reps <- if (is.null(opt$reps)) cfg$seeds$repsPerShard else as.integer(opt$reps)
designs <- if (is.null(opt$designs)) c("paired", "unpaired") else split_list(opt$designs)
arms <- cfg$exp2$arms
outFile <- file.path(opt$out, sprintf("shard%d.csv", shard))
timing <- file.path(opt$out, sprintf("shard%d_timing.tsv", shard))

## Outcomes of one run (Table 3): gain and genic variance at the last cycle, accuracy
## averaged over the late window.
summarise_run <- function(rec, design, r, arm, seed) {
  rec <- rec[order(rec$cycle), ]
  last <- nrow(rec)
  late <- seq.int(last - cfg$exp2$lateWindow + 1L, last)
  data.frame(design = design, seedBase = base, rep = r, arm = arm,
             gain = rec$gain[last], genicPct = rec$genicPct[last],
             accLate = mean(rec$accuracy[late], na.rm = TRUE), seed = seed,
             stringsAsFactors = FALSE)
}

done <- if (file.exists(outFile)) {
  d <- read.csv(outFile, stringsAsFactors = FALSE)
  unique(paste(d$design, d$rep))
} else character(0)

for (r in seq_len(reps)) for (design in designs) {
  if (paste(design, r) %in% done) {
    log_msg("exp2 shard %d: %s replicate %d already written", shard, design, r); next
  }
  t0 <- Sys.time()
  if (design == "paired") {
    seed <- base + r
    founder <- build_founders(cfg, seed)
    block <- do.call(rbind, lapply(arms, function(a)
      summarise_run(run_programme(founder, cfg, a, "baseline", programme_seed(cfg, seed), rep = r),
                    design, r, a, seed)))
  } else {
    block <- do.call(rbind, lapply(seq_along(arms), function(k) {
      seed <- base + cfg$seeds$unpairedStride * k + r
      founder <- build_founders(cfg, seed)
      summarise_run(run_programme(founder, cfg, arms[k], "baseline", programme_seed(cfg, seed), rep = r),
                    design, r, arms[k], seed)
    }))
  }
  append_csv(block, outFile)
  append_timing(timing, shard = shard, design = design, rep = r,
                seconds = round(seconds_since(t0), 1))
  log_msg("exp2 shard %d: %s replicate %d written (%.0f s)", shard, design, r, seconds_since(t0))
}
log_msg("exp2 shard %d: complete (%s)", shard, outFile)
