#!/usr/bin/env Rscript
## =============================================================================
## run_exp1.R -- Experiment 1: reproducibility under one fixed seed, layer by layer
##
## Layers (Section 2.3):
##   1 runMacs     founder haplotypes from MaCS (compiled code), generator seed
##   2 quickHaplo  founder haplotypes from R's generator (positive control), same seed
##   3 programme   CLOSED, DIRECT, BRIDGE and MAS, each run programmeCalls times on one
##                 founder object (Table 1, Experiment 1 column) with programme seed 7
##   4 programme   the same once more on that object after saveRDS/readRDS
## A genome is summarised as "mean dosage|dosage sum of individual 1|number of sites", a
## programme run as "final cumulative gain|final genic variance" (t/ha, 10 decimals).
##
## Usage:
##   Rscript sim/run_exp1.R within   --out out_sim/exp1 [--config exp1]
##   Rscript sim/run_exp1.R across   --label 1 --file out_sim/exp1/across_sessions.csv [--config exp1]
##   Rscript sim/run_exp1.R sessions --out out_sim/exp1 [--n 3] [--config exp1]
##   Rscript sim/run_exp1.R all      --out out_sim/exp1 [--config exp1]
## Outputs: within_session.csv (layer, rep, signature, note), across_sessions.csv (session,
## pid, layer, signature, AlphaSimR, R), table2.csv (layer, calls, distinct, reproducible),
## and the founder object founder_<seed>.rds with its checksum.
## =============================================================================

Sys.setenv(OMP_NUM_THREADS = "1")
SIM_DIR <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  if (length(f)) dirname(normalizePath(f[1])) else "sim"
})
source(file.path(SIM_DIR, "R", "load.R"))
load_sim(SIM_DIR)

opt <- parse_cli(commandArgs(trailingOnly = TRUE), positional = "mode")
mode <- if (is.null(opt$mode)) "all" else opt$mode
cfg <- get_config(if (is.null(opt$config)) "exp1" else opt$config)
G <- cfg$exp1$generator

genome_signature <- function(mapPop) {
  M <- pullSegSiteGeno(mapPop)
  sprintf("%.10f|%d|%d", mean(M), as.integer(sum(M[1, ])), ncol(M))
}
programme_signature <- function(rec) {
  last <- rec[rec$cycle == max(rec$cycle), ]
  sprintf("%.10f|%.10f", last$gain, last$genicVar)
}
call_runMacs <- function() {
  seed_rng(cfg$exp1$generatorSeed)
  runMacs(nInd = G$nInd, nChr = G$nChr, segSites = G$segSites, species = cfg$genome$species,
          split = cfg$genome$splitGenerations, nThreads = 1L)
}
call_quickHaplo <- function() {
  seed_rng(cfg$exp1$generatorSeed)
  quickHaplo(nInd = G$nInd, nChr = G$nChr, segSites = G$segSites)
}
GEN_LAYERS <- c(runMacs = "runMacs (MaCS coalescent, compiled)",
                quickHaplo = "quickHaplo (R generator)")

## ---- layers 1-4 within one session --------------------------------------------
run_within <- function(out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  add <- function(layer, rep, signature, note = "")
    rows[[length(rows) + 1]] <<- data.frame(layer = layer, rep = rep, signature = signature,
                                            note = note, stringsAsFactors = FALSE)
  for (k in seq_len(G$calls)) add(GEN_LAYERS[["runMacs"]], k, genome_signature(call_runMacs()))
  for (k in seq_len(G$calls)) add(GEN_LAYERS[["quickHaplo"]], k, genome_signature(call_quickHaplo()))
  log_msg("exp1: founder layers done (%d calls each)", G$calls)

  fseed <- cfg$seeds$founderSeed; pseed <- cfg$seeds$programmeSeed
  founder <- build_founders(cfg, fseed)
  path <- founder_path(out, fseed)
  save_founders(founder, path)
  reloaded <- read_founders(path, cfg, fseed)
  for (a in cfg$exp1$arms) {
    layer <- paste("programme", a)
    for (k in seq_len(cfg$exp1$programmeCalls))
      add(layer, k, programme_signature(run_programme(founder, cfg, a, "baseline", pseed)),
          "founder object in memory")
    add(layer, cfg$exp1$programmeCalls + 1L,
        programme_signature(run_programme(reloaded, cfg, a, "baseline", pseed)),
        "founder object after saveRDS/readRDS")
    log_msg("exp1: programme %s done", a)
  }
  d <- do.call(rbind, rows)
  write_csv_atomic(d, file.path(out, "within_session.csv"))
  d
}

## ---- one call of each generator, appended to a CSV (one R session per call) -----------
run_across <- function(label, file) {
  d <- data.frame(session = label, pid = Sys.getpid(),
                  layer = unname(GEN_LAYERS),
                  signature = c(genome_signature(call_runMacs()), genome_signature(call_quickHaplo())),
                  AlphaSimR = as.character(utils::packageVersion("AlphaSimR")),
                  R = R.version.string, stringsAsFactors = FALSE)
  append_csv(d, file)
  d
}

## ---- n separate Rscript sessions ---------------------------------------------------
run_sessions <- function(out, n) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(out, "across_sessions.csv")
  if (file.exists(file)) file.remove(file)
  rscript <- file.path(R.home("bin"), "Rscript")
  for (s in seq_len(n)) {
    status <- system2(rscript, c(shQuote(file.path(SIM_DIR, "run_exp1.R")), "across",
                                 "--label", s, "--file", shQuote(file), "--config", cfg$name))
    if (status != 0) stop("session ", s, " failed")
  }
  read.csv(file, stringsAsFactors = FALSE)
}

## ---- Table 2: distinct results per layer ------------------------------------------
summarise_layers <- function(within, across, out) {
  tab <- function(d, label) do.call(rbind, lapply(split(d, d$layer), function(z)
    data.frame(layer = paste0(z$layer[1], label), calls = nrow(z),
               distinct = length(unique(z$signature)),
               reproducible = length(unique(z$signature)) == 1, stringsAsFactors = FALSE)))
  t2 <- rbind(if (!is.null(within)) tab(within, ", one session"),
              if (!is.null(across)) tab(across, ", separate sessions"))
  write_csv_atomic(t2, file.path(out, "table2.csv"))
  print(t2, row.names = FALSE)
  t2
}

if (mode == "within") {
  w <- run_within(opt$out); summarise_layers(w, NULL, opt$out)
} else if (mode == "across") {
  print(run_across(opt$label, opt$file), row.names = FALSE)
} else if (mode == "sessions") {
  s <- run_sessions(opt$out, if (is.null(opt$n)) G$sessions else as.integer(opt$n))
  summarise_layers(NULL, s, opt$out)
} else if (mode == "all") {
  w <- run_within(opt$out)
  s <- run_sessions(opt$out, if (is.null(opt$n)) G$sessions else as.integer(opt$n))
  summarise_layers(w, s, opt$out)
} else {
  stop("mode must be within, across, sessions or all")
}
