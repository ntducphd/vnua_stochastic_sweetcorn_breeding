## =============================================================================
## programme.R -- running an arm, and a whole replicate of the application
## =============================================================================

## Run one arm in one scenario for cfg$programme$nCycles cycles from a founder object.
##   founder        object from build_founders()/load_or_build_founders(); not modified
##   programmeSeed  set immediately before cycle 1, after the founders exist
##   schedule       admissions per cycle for an ELITE2M arm
##   trace          keep per-cycle internals (training sample, parents, generator states)
## Returns the per-cycle records (data.frame, RECORD_COLUMNS) with attribute "trace".
run_programme <- function(founder, cfg, armName, scenarioName, programmeSeed,
                          schedule = NULL, rep = NA_integer_, trace = FALSE) {
  scenario <- SCENARIOS[[scenarioName]]
  if (is.null(scenario)) stop("unknown scenario: ", scenarioName)
  if (scenario$resourceMatch && !isTRUE(cfg$resourceMatching$enabled))
    stop("resource matching is not enabled in configuration ", cfg$name)
  arm <- build_arm(armName, cfg, scenario, schedule)

  ## Each run works on its own copy of the SimParam (an R6 reference object), so runs
  ## never see each other's ID counter or any other state change.
  sp <- founder$SP$clone(deep = TRUE)
  sp$nThreads <- 1L
  env <- list(cfg = cfg, scenario = scenario, sp = sp, founder = founder)
  state <- list(nursery = founder$nursery, nurseryProgeny = founder$nurseryProgeny,
                archive = new_archive(cfg$programme$archiveSize))

  seed_rng(programmeSeed)
  rows <- vector("list", cfg$programme$nCycles)
  traces <- if (trace) vector("list", cfg$programme$nCycles) else NULL
  for (t in seq_len(cfg$programme$nCycles)) {
    step <- run_cycle(state, t, arm, env, trace = trace)
    state <- step$state
    rows[[t]] <- step$record
    if (trace) traces[[t]] <- step$trace
  }
  out <- assemble_records(rows, founder, cfg, armName, scenario, founder$seed, rep)
  if (trace) attr(out, "trace") <- traces
  out
}

## The programme seed of a replicate (Section 2.6: one founder object and one programme
## seed per replicate, shared by every arm and scenario).
programme_seed <- function(cfg, seed) seed + cfg$seeds$programmeOffset

## Admission schedule of an inflow arm, as replayed by its ELITE2M control.
admission_schedule <- function(records) as.integer(records$nAdmitted[order(records$cycle)])

## One replicate of the application: every main arm in every requested scenario, then the
## count-matched controls, which replay the admissions their inflow arm made in the same
## call (same founders, same programme seed, same scenario).
##   elite2mScenarios  scenarios in which ELITE2M arms are run (default: baseline and
##                     resource_matched; under visual sorting they would repeat these)
##   onRun             optional callback(scenario, arm, seconds) for logging
run_replicate <- function(cfg, seed, founderDir = NULL, scenarios = cfg$scenarios,
                          arms = MAIN_ARMS,
                          elite2mScenarios = intersect(c("baseline", "resource_matched"), scenarios),
                          founder = NULL, onRun = NULL, verbose = TRUE) {
  stopifnot(all(arms %in% MAIN_ARMS), all(elite2mScenarios %in% scenarios))
  if (is.null(founder)) founder <- load_or_build_founders(cfg, seed, founderDir, verbose = verbose)
  pseed <- programme_seed(cfg, seed)
  rep <- as.integer(seed - cfg$seeds$repOffset)
  runs <- list()
  timed <- function(sc, a, expr) {
    t0 <- Sys.time()
    res <- force(expr)
    s <- seconds_since(t0)
    if (verbose) log_msg("seed %d  %-16s %-15s %6.1f s", seed, sc, a, s)
    if (!is.null(onRun)) onRun(sc, a, s)
    res
  }
  for (sc in scenarios) for (a in arms)
    runs[[paste(sc, a)]] <- timed(sc, a, run_programme(founder, cfg, a, sc, pseed, rep = rep))
  for (sc in elite2mScenarios) for (a in intersect(INFLOW_ARMS, arms)) {
    sched <- admission_schedule(runs[[paste(sc, a)]])
    name <- paste0(MATCHED_PREFIX, a)
    res <- timed(sc, name, run_programme(founder, cfg, name, sc, pseed, schedule = sched, rep = rep))
    if (!identical(admission_schedule(res), sched))
      stop(sprintf("%s did not replay the schedule of %s (seed %d, %s)", name, a, seed, sc))
    runs[[paste(sc, name)]] <- res
  }
  out <- do.call(rbind, unname(runs))
  rownames(out) <- NULL
  out
}
