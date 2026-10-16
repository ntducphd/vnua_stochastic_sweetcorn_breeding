#!/usr/bin/env Rscript
## =============================================================================
## run_application.R -- one replicate of the donor-introgression application
##
## Runs CLOSED, ELITE2, DIRECT, BRIDGE and MAS in every scenario, then the count-matched
## controls ELITE2M_<arm>, all from the replicate's founder object (read from
## --founders if saved there, otherwise built and saved with a SHA-256 checksum) and the
## replicate's programme seed. Writes <out>/raw_<seed>.csv (all arms and scenarios) and
## appends per-run timings to <out>/timing_<seed>.tsv. A replicate whose CSV exists is
## skipped, so the driver can be rerun over a list of seeds until all are done.
##
## Usage:
##   Rscript sim/run_application.R --seed 4001 --out out_sim/application \
##       --founders out_sim/founders [--config application]
##       [--scenarios baseline,resource_matched,visual_sort,visual_matched]
##       [--arms CLOSED,ELITE2,DIRECT,BRIDGE,MAS]
##       [--elite2m-scenarios baseline,resource_matched | all | none]
##       [--cycles N] [--pipeline-layout] [--force]
##
##   --pipeline-layout  also write raw_w<seed>.csv and raw_m<seed>.csv, the layout read
##                      by R/01_pool_application.py and R/04_analyse_matched_controls.R
##   --cycles N         shorten the run (timing estimates only; not for analysis)
##   --force            rerun even if the replicate's CSV exists
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
if (is.null(opt$seed) || is.null(opt$out)) {
  cat("usage: Rscript sim/run_application.R --seed <seed> --out <dir> [--founders <dir>] [options]\n")
  quit(status = 2)
}
seed <- as.integer(opt$seed)
cfg <- get_config(if (is.null(opt$config)) "application" else opt$config)
if (!is.null(opt$cycles)) cfg <- modify_config(cfg, list(programme = list(nCycles = as.integer(opt$cycles))))
scenarios <- if (is.null(opt$scenarios)) cfg$scenarios else split_list(opt$scenarios)
arms <- if (is.null(opt$arms)) MAIN_ARMS else split_list(opt$arms)
e2m <- if (is.null(opt[["elite2m-scenarios"]])) {
  intersect(c("baseline", "resource_matched"), scenarios)
} else if (identical(opt[["elite2m-scenarios"]], "all")) {
  scenarios
} else if (identical(opt[["elite2m-scenarios"]], "none")) {
  character(0)
} else split_list(opt[["elite2m-scenarios"]])
founderDir <- if (is.null(opt$founders)) file.path(opt$out, "founders") else opt$founders

outFile <- replicate_path(opt$out, seed)
timing <- file.path(opt$out, sprintf("timing_%d.tsv", seed))
if (file.exists(outFile) && !isTRUE(opt$force)) {
  log_msg("seed %d: %s exists, skipped", seed, outFile)
  quit(status = 0)
}
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

log_msg("seed %d: configuration %s, %d cycles, scenarios %s, arms %s, ELITE2M in %s",
        seed, cfg$name, cfg$programme$nCycles, paste(scenarios, collapse = ","),
        paste(arms, collapse = ","), if (length(e2m)) paste(e2m, collapse = ",") else "none")
t0 <- Sys.time()
founder <- load_or_build_founders(cfg, seed, founderDir)
append_timing(timing, seed = seed, scenario = "-", arm = paste0("founders (", attr(founder, "source"), ")"),
              seconds = round(seconds_since(t0), 2))

records <- run_replicate(cfg, seed, founder = founder, scenarios = scenarios, arms = arms,
                         elite2mScenarios = e2m,
                         onRun = function(sc, a, s) append_timing(timing, seed = seed, scenario = sc,
                                                                  arm = a, seconds = round(s, 2)))
write_records(records, outFile)
if (isTRUE(opt[["pipeline-layout"]])) split_for_pipeline(records, opt$out, seed)
total <- seconds_since(t0)
append_timing(timing, seed = seed, scenario = "-", arm = "total", seconds = round(total, 2))
log_msg("seed %d: wrote %s (%d rows) in %.1f min", seed, outFile, nrow(records), total / 60)
