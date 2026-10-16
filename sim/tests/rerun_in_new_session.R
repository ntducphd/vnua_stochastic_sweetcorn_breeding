## =============================================================================
## rerun_in_new_session.R -- helper for run_tests.R
##
## Called with Rscript (a fresh R session) to rerun a replicate from a saved founder
## object, so that the test can compare the records with those of the parent session.
##   Rscript rerun_in_new_session.R <simDir> <config.rds> <founderDir> <seed> <out.csv>
## =============================================================================
Sys.setenv(OMP_NUM_THREADS = "1")
args <- commandArgs(trailingOnly = TRUE)
simDir <- args[1]
source(file.path(simDir, "R", "load.R"))
load_sim(simDir)
cfg <- readRDS(args[2])
seed <- as.integer(args[4])
founder <- read_founders(founder_path(args[3], seed), cfg, seed)
records <- run_replicate(cfg, seed, founder = founder, verbose = FALSE)
write_records(records, args[5])
cat("rerun written:", args[5], nrow(records), "rows, pid", Sys.getpid(), "\n")
