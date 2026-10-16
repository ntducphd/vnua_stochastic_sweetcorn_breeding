## =============================================================================
## utils.R -- threads, the random number generator, logging and small helpers
## =============================================================================

## One thread everywhere (Section 2.8). OMP_NUM_THREADS must be set before AlphaSimR is
## loaded; load.R does that, and every SimParam created by the engine gets nThreads = 1.
use_one_thread <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1")
  invisible(TRUE)
}

## Seed R's generator with an explicit kind, so that a run does not depend on the
## RNGkind() defaults of the session that happens to run it.
seed_rng <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  invisible(seed)
}

## Copy of the generator state (used by the trace that tests where two arms' streams part).
rng_state <- function() get(".Random.seed", envir = globalenv())

## Versions written into every record and founder object.
software_versions <- function() {
  list(R = paste(R.version$major, R.version$minor, sep = "."),
       RString = R.version.string,
       AlphaSimR = as.character(utils::packageVersion("AlphaSimR")),
       rrBLUP = as.character(utils::packageVersion("rrBLUP")),
       platform = R.version$platform,
       ompThreads = Sys.getenv("OMP_NUM_THREADS", unset = NA_character_))
}

## Time-stamped progress message on stderr (stdout stays free for results).
log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
}

## Seconds elapsed since t0 (a proc.time() or Sys.time() value).
seconds_since <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

## Indices of the n largest values, largest first; ties keep the original order.
top_n <- function(x, n) {
  if (n <= 0 || !length(x)) return(integer(0))
  order(x, decreasing = TRUE, method = "radix")[seq_len(min(n, length(x)))]
}

## Population variance (divisor N), the convention of AlphaSimR's varA()/varG(), which also
## defines the founder variance of 1 that the trait is scaled to.
pop_var <- function(x) {
  if (!length(x)) return(NA_real_)
  m <- mean(x)
  mean((x - m)^2)
}

## Pearson correlation that returns NA (silently) when it is undefined.
safe_cor <- function(x, y) {
  if (length(x) < 2 || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y)
}

## Minimal "--key value" / "--flag" command-line parser for the drivers.
parse_cli <- function(args, positional = character(0)) {
  out <- list(); pos <- character(0); k <- 1
  while (k <= length(args)) {
    a <- args[k]
    if (startsWith(a, "--")) {
      key <- substring(a, 3)
      if (grepl("=", key, fixed = TRUE)) {
        out[[sub("=.*", "", key)]] <- sub("^[^=]*=", "", key)
      } else if (k < length(args) && !startsWith(args[k + 1], "--")) {
        out[[key]] <- args[k + 1]; k <- k + 1
      } else {
        out[[key]] <- TRUE
      }
    } else {
      pos <- c(pos, a)
    }
    k <- k + 1
  }
  for (j in seq_along(positional)) if (j <= length(pos)) out[[positional[j]]] <- pos[j]
  out$.positional <- pos
  out
}

## Comma-separated list option -> character vector (NULL stays NULL).
split_list <- function(x) if (is.null(x) || isTRUE(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]
