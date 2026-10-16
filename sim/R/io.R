## =============================================================================
## io.R -- writing and reading the engine's outputs
##
## Every file is written to "<path>.partial" and renamed when complete, so an interrupted
## run never leaves a file that a resumed run would mistake for a finished one.
## =============================================================================

write_csv_atomic <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".partial")
  utils::write.csv(df, tmp, row.names = FALSE)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("could not move ", tmp, " to ", path)
  invisible(path)
}

## Per-cycle records of one replicate, columns in RECORD_COLUMNS order.
write_records <- function(records, path) {
  missing <- setdiff(RECORD_COLUMNS, names(records))
  if (length(missing)) stop("records lack columns: ", paste(missing, collapse = ", "))
  write_csv_atomic(records[, RECORD_COLUMNS], path)
}

read_records <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)

## The replicate file written by run_application.R.
replicate_path <- function(dir, seed) file.path(dir, sprintf("raw_%d.csv", seed))

## The layout the existing analysis scripts read: raw_w<seed>.csv (CLOSED, ELITE2, DIRECT,
## BRIDGE, MAS; read by R/01_pool_application.py) and raw_m<seed>.csv (ELITE2M arms; read
## by R/04_analyse_matched_controls.R). Both keep RECORD_COLUMNS.
split_for_pipeline <- function(records, dir, seed) {
  m <- startsWith(records$strategy, MATCHED_PREFIX)
  paths <- c(w = file.path(dir, sprintf("raw_w%d.csv", seed)),
             m = file.path(dir, sprintf("raw_m%d.csv", seed)))
  write_records(records[!m, , drop = FALSE], paths[["w"]])
  if (any(m)) write_records(records[m, , drop = FALSE], paths[["m"]])
  invisible(paths)
}

## Append rows to a CSV (header only when the file is new); used for the Experiment 1
## across-session log and the Experiment 2 shards, which grow one block at a time.
append_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  new <- !file.exists(path)
  utils::write.table(df, path, sep = ",", row.names = FALSE, col.names = new, append = !new,
                     qmethod = "double")
  invisible(path)
}

## Tab-separated timing log of a driver.
append_timing <- function(path, ...) {
  row <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., stringsAsFactors = FALSE)
  new <- !file.exists(path)
  utils::write.table(row, path, sep = "\t", row.names = FALSE, col.names = new, append = !new,
                     quote = FALSE)
}
