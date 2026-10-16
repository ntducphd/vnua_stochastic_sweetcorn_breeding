## =============================================================================
## load.R -- load the engine: one thread, packages, then the modules in dependency order
##
##   source("sim/R/load.R"); load_sim("sim")
## =============================================================================

SIM_MODULES <- c("utils.R", "config.R", "genetics.R", "prediction.R", "founders.R",
                 "outcomes.R", "arms.R", "cycle.R", "programme.R", "io.R")

load_sim <- function(simDir, envir = globalenv()) {
  if ("AlphaSimR" %in% loadedNamespaces() && !identical(Sys.getenv("OMP_NUM_THREADS"), "1"))
    warning("AlphaSimR was loaded before OMP_NUM_THREADS was set to 1; ",
            "the engine still passes nThreads = 1 to every AlphaSimR call")
  Sys.setenv(OMP_NUM_THREADS = "1")      # before AlphaSimR is loaded
  suppressPackageStartupMessages({
    library(AlphaSimR)
    library(rrBLUP)
  })
  for (m in SIM_MODULES) sys.source(file.path(simDir, "R", m), envir = envir)
  invisible(simDir)
}

## Directory of the script being run by Rscript (for the drivers and tests).
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[startsWith(a, "--file=")])
  if (!length(f)) return(normalizePath("."))
  dirname(normalizePath(f[1]))
}
