## =============================================================================
## config.R -- parameter sets and their validation
##
## A configuration is a plain nested list with seven blocks:
##   genome           chromosomes, segregating sites, QTL and markers per chromosome
##   founders         coalescent draw, cull at the endosperm locus, burn-in, t/ha scale
##   programme        cycle operator of Algorithm 1 (parents, crosses, training archive, cap)
##   introgression    donor crosses, F2 sizes, the BRIDGE second step, the MAS screen
##   resourceMatching the 200-plot introgression budget and the smaller elite nursery
##   seeds            replicate seeds and how rep / programme seed are derived from them
##   sugar            the sugar trait of the endosperm locus (recorded, never simulated)
## Three sets follow Table 1 of the paper: config_application(), config_exp1() and
## config_exp2(). modify_config() derives variants (used by the tests) and re-validates.
## Nothing in this file draws a random number.
## =============================================================================

## The four scenarios of Section 2.2: two discard assumptions x two resource settings.
SCENARIOS <- list(
  baseline         = list(name = "baseline",         visualSort = FALSE, resourceMatch = FALSE),
  resource_matched = list(name = "resource_matched", visualSort = FALSE, resourceMatch = TRUE),
  visual_sort      = list(name = "visual_sort",      visualSort = TRUE,  resourceMatch = FALSE),
  visual_matched   = list(name = "visual_matched",   visualSort = TRUE,  resourceMatch = TRUE)
)

INFLOW_ARMS <- c("DIRECT", "BRIDGE", "MAS")
MAIN_ARMS   <- c("CLOSED", "ELITE2", INFLOW_ARMS)
MATCHED_PREFIX <- "ELITE2M_"

## -----------------------------------------------------------------------------
## Table 1, application column
## -----------------------------------------------------------------------------
config_application <- function() {
  cfg <- list(
    name = "application",
    genome = list(
      nChr = 10L, segSites = 250L,          # 10 chromosomes x 250 segregating sites
      nQtlPerChr = 100L, nSnpPerChr = 100L, # disjoint QTL and marker sets
      species = "MAIZE",                    # MaCS maize history built into runMacs
      splitGenerations = 100L               # the two pools split 100 generations ago
    ),
    founders = list(
      nSimulatedPerPool = 180L,             # founders simulated per pool ...
      nKeptPerPool = 80L,                   # ... and the first survivors kept per pool
      h2 = 0.30,                            # in the combined founder population, before the cull
      burnInCycles = 8L,                    # phenotypic selection in the elite pool only
      eliteMean = 14.0,                     # t/ha of the kept elite founders
      eliteVarA = 1.0,                      # (t/ha)^2 of the kept elite founders
      maxRedraws = 10000L,                  # safety limit on QTL/marker redraws
      maxCoalescentAttempts = 10L           # safety limit on coalescent draws (cull failures)
    ),
    programme = list(
      nCycles = 20L, yearsPerCycle = 3L,
      nParents = 40L,                       # n_p
      nCrosses = 40L, nProgeny = 15L,       # n_c x n_o = 600-plot elite nursery
      nTrain = 250L, archiveSize = 1000L,   # training lines per cycle; latest records kept
      cap = 2L                              # kappa, admissions per cycle
    ),
    introgression = list(
      nDonorCrosses = 6L, f2PerCross = 150L,
      bridgeFamilies = 20L, bridgeProgeny = 45L,
      masScreenFactor = 4L                  # MAS screens 4 x the F2 by marker
    ),
    resourceMatching = list(
      enabled = TRUE,
      introgressionPlots = 200L,            # of the control's 600 plots
      nurseryProgeny = 10L,                 # elite nursery 40 x 10
      directF2PerCross = 33L,
      bridgeF2PerCross = 16L,               # BRIDGE step 1; step 2 gets the other half
      bridgeProgeny = 5L,                   # 20 families x 5 = 100 plots (see SPEC_NOTES)
      masSeedlingsPerCross = 132L
    ),
    scenarios = names(SCENARIOS),
    seeds = list(
      replicates = 4001:4010,
      repOffset = 4000L,                    # rep = seed - repOffset
      programmeOffset = 0L                  # programme seed = replicate seed + offset
    ),
    sugar = list(additive = 5.82, dominance = -5.82, residualSd = 0.46)
  )
  validate_config(cfg)
}

## -----------------------------------------------------------------------------
## Table 1, Experiment 1 column (programme layer) plus the founder-generator layers
## -----------------------------------------------------------------------------
config_exp1 <- function() {
  cfg <- config_application()
  cfg$name <- "exp1"
  cfg$founders$nSimulatedPerPool <- 120L
  cfg$founders$nKeptPerPool <- 60L
  cfg$founders$burnInCycles <- 3L
  cfg$programme$nCycles <- 4L
  cfg$resourceMatching$enabled <- FALSE      # "not run"
  cfg$scenarios <- "baseline"
  cfg$seeds <- list(replicates = 20260925L, repOffset = 0L, programmeOffset = 0L,
                    founderSeed = 20260925L, programmeSeed = 7L)
  cfg$exp1 <- list(
    arms = c("CLOSED", "DIRECT", "BRIDGE", "MAS"),
    programmeCalls = 5L,                     # in-memory calls; one more after saveRDS/readRDS
    generatorSeed = 20260925L,               # layers 1 and 2
    generator = list(nInd = 40L, nChr = 3L, segSites = 150L, calls = 5L, sessions = 3L)
  )
  validate_config(cfg)
}

## -----------------------------------------------------------------------------
## Table 1, Experiment 2 column (reduced configuration)
## -----------------------------------------------------------------------------
config_exp2 <- function() {
  cfg <- config_application()
  cfg$name <- "exp2"
  cfg$genome$nChr <- 4L
  cfg$genome$segSites <- 120L
  cfg$genome$nQtlPerChr <- 60L
  cfg$genome$nSnpPerChr <- 40L
  cfg$founders$nSimulatedPerPool <- 150L
  cfg$founders$nKeptPerPool <- 40L
  cfg$founders$burnInCycles <- 5L
  cfg$programme$nCycles <- 10L
  cfg$programme$nParents <- 20L
  cfg$programme$nCrosses <- 20L
  cfg$programme$nProgeny <- 10L
  cfg$programme$nTrain <- 150L
  cfg$introgression <- list(nDonorCrosses = 4L, f2PerCross = 60L, bridgeFamilies = 10L,
                            bridgeProgeny = 24L, masScreenFactor = 4L)
  cfg$resourceMatching$enabled <- FALSE      # "not run"
  cfg$scenarios <- "baseline"
  cfg$seeds <- list(replicates = integer(0), repOffset = 0L, programmeOffset = 0L,
                    shardBases = c(500000L, 510000L, 520000L, 530000L),
                    repsPerShard = 10L, unpairedStride = 100000L)
  cfg$exp2 <- list(arms = c("CLOSED", "DIRECT", "BRIDGE", "MAS"),
                   lateWindow = 5L)          # accLate = mean accuracy over the last 5 cycles
  validate_config(cfg)
}

## Look a configuration up by name (used by the command-line drivers).
get_config <- function(name) {
  switch(name,
         application = config_application(),
         exp1 = config_exp1(),
         exp2 = config_exp2(),
         test = config_test(),
         stop(sprintf("unknown configuration '%s' (application, exp1, exp2, test)", name),
              call. = FALSE))
}

## A small configuration for smoke tests of the drivers (minutes, not hours). It keeps the
## structure of the application (all four scenarios, resource matching) at toy sizes.
config_test <- function() {
  cfg <- config_application()
  cfg$name <- "test"
  cfg$genome <- list(nChr = 2L, segSites = 60L, nQtlPerChr = 15L, nSnpPerChr = 15L,
                     species = "MAIZE", splitGenerations = 100L)
  cfg$founders$nSimulatedPerPool <- 50L
  cfg$founders$nKeptPerPool <- 20L
  cfg$founders$burnInCycles <- 2L
  cfg$programme <- list(nCycles = 4L, yearsPerCycle = 3L, nParents = 10L, nCrosses = 10L,
                        nProgeny = 6L, nTrain = 30L, archiveSize = 100L, cap = 2L)
  cfg$introgression <- list(nDonorCrosses = 3L, f2PerCross = 40L, bridgeFamilies = 8L,
                            bridgeProgeny = 10L, masScreenFactor = 4L)
  cfg$resourceMatching <- c(list(enabled = TRUE, introgressionPlots = 60L, nurseryProgeny = 4L),
                            derive_resource_matching(60L, 3L, 8L, 4L))
  cfg$seeds <- list(replicates = 1:3, repOffset = 0L, programmeOffset = 0L,
                    shardBases = c(500000L, 510000L, 520000L, 530000L),
                    repsPerShard = 1L, unpairedStride = 100000L,
                    founderSeed = 11L, programmeSeed = 7L)
  cfg$exp1 <- list(arms = c("CLOSED", "DIRECT"), programmeCalls = 2L, generatorSeed = 20260925L,
                   generator = list(nInd = 20L, nChr = 1L, segSites = 40L, calls = 2L,
                                    sessions = 2L))
  cfg$exp2 <- list(arms = c("CLOSED", "DIRECT", "BRIDGE", "MAS"), lateWindow = 2L)
  validate_config(cfg)
}

## Recursively replace entries of a configuration and validate the result.
## Example: modify_config(cfg, list(programme = list(nCycles = 3L)))
modify_config <- function(cfg, changes) {
  validate_config(utils::modifyList(cfg, changes))
}

## -----------------------------------------------------------------------------
## Resource matching (Section 2.2): the introgression budget fixes the F2 counts.
## DIRECT spends it on F2; BRIDGE splits it equally between its two steps (step 2 keeps
## the bridge families and shrinks their progeny); MAS screens masFactor x DIRECT's F2.
## -----------------------------------------------------------------------------
derive_resource_matching <- function(budget, nDonorCrosses, bridgeFamilies, masFactor) {
  half <- budget %/% 2L
  list(directF2PerCross = as.integer(budget %/% nDonorCrosses),
       bridgeF2PerCross = as.integer(half %/% nDonorCrosses),
       bridgeProgeny = as.integer(half %/% bridgeFamilies),
       masSeedlingsPerCross = as.integer(masFactor * (budget %/% nDonorCrosses)))
}

## -----------------------------------------------------------------------------
## Validation: types, ranges and the relations between parameters the engine relies on.
## Returns the configuration (invisibly) or stops with every problem found.
## -----------------------------------------------------------------------------
validate_config <- function(cfg) {
  problems <- character(0)
  need <- function(ok, msg) if (!isTRUE(ok)) problems <<- c(problems, msg)
  is_count <- function(x, min = 1) is.numeric(x) && length(x) == 1 && !is.na(x) &&
    x == round(x) && x >= min

  for (b in c("genome", "founders", "programme", "introgression", "resourceMatching", "seeds"))
    need(is.list(cfg[[b]]), sprintf("block '%s' is missing", b))
  if (length(problems)) stop(paste(c("invalid configuration:", problems), collapse = "\n  "),
                             call. = FALSE)
  g <- cfg$genome; f <- cfg$founders; p <- cfg$programme
  i <- cfg$introgression; r <- cfg$resourceMatching

  ## genome: the endosperm locus must be able to fall outside both sets
  for (k in c("nChr", "segSites", "nQtlPerChr", "nSnpPerChr")) need(is_count(g[[k]]), paste0("genome$", k))
  need(g$nQtlPerChr + g$nSnpPerChr < g$segSites,
       "genome: QTL + markers per chromosome must leave at least one free site for the locus")
  need(is.character(g$species) && length(g$species) == 1, "genome$species")
  need(is_count(g$splitGenerations), "genome$splitGenerations")

  ## founders
  for (k in c("nSimulatedPerPool", "nKeptPerPool", "burnInCycles", "maxRedraws",
              "maxCoalescentAttempts")) need(is_count(f[[k]]), paste0("founders$", k))
  need(f$nKeptPerPool <= f$nSimulatedPerPool, "founders: kept per pool exceeds simulated per pool")
  need(f$nKeptPerPool >= 2, "founders: at least two founders per pool must be kept")
  need(is.numeric(f$h2) && f$h2 > 0 && f$h2 < 1, "founders$h2 must lie in (0, 1)")
  need(is.numeric(f$eliteMean) && length(f$eliteMean) == 1, "founders$eliteMean")
  need(is.numeric(f$eliteVarA) && f$eliteVarA > 0, "founders$eliteVarA must be positive")

  ## programme
  for (k in c("nCycles", "yearsPerCycle", "nParents", "nCrosses", "nProgeny", "nTrain",
              "archiveSize")) need(is_count(p[[k]]), paste0("programme$", k))
  need(is_count(p$cap, 0), "programme$cap")
  need(p$nParents >= 2, "programme: at least two parents are needed for random mating")
  nursery <- p$nCrosses * p$nProgeny
  need(p$nCrosses <= choose(p$nParents, 2), "programme: more crosses than parent pairs")
  need(p$nTrain <= nursery, "programme: training lines exceed the elite nursery")
  need(p$archiveSize >= p$nTrain, "programme: archive smaller than one training sample")
  need(p$nParents + p$cap <= nursery,
       "programme: the nursery cannot supply the parents plus the elite controls' admissions")

  ## introgression
  for (k in c("nDonorCrosses", "f2PerCross", "bridgeFamilies", "bridgeProgeny", "masScreenFactor"))
    need(is_count(i[[k]]), paste0("introgression$", k))

  ## resource matching (only checked when the scenario can be run)
  need(is.logical(r$enabled) && length(r$enabled) == 1, "resourceMatching$enabled")
  if (isTRUE(r$enabled)) {
    for (k in c("introgressionPlots", "nurseryProgeny", "directF2PerCross", "bridgeF2PerCross",
                "bridgeProgeny", "masSeedlingsPerCross"))
      need(is_count(r[[k]]), paste0("resourceMatching$", k))
    d <- derive_resource_matching(r$introgressionPlots, i$nDonorCrosses, i$bridgeFamilies,
                                  i$masScreenFactor)
    for (k in names(d))
      need(identical(as.integer(r[[k]]), d[[k]]),
           sprintf("resourceMatching$%s = %s, but the %d-plot budget implies %d",
                   k, format(r[[k]]), r$introgressionPlots, d[[k]]))
    need(p$nTrain <= p$nCrosses * r$nurseryProgeny,
         "resourceMatching: training lines exceed the smaller elite nursery")
    need(p$nParents + p$cap <= p$nCrosses * r$nurseryProgeny,
         "resourceMatching: the smaller nursery cannot supply parents plus admissions")
  }

  ## scenarios
  need(is.character(cfg$scenarios) && length(cfg$scenarios) >= 1 &&
         all(cfg$scenarios %in% names(SCENARIOS)),
       sprintf("scenarios must be among %s", paste(names(SCENARIOS), collapse = ", ")))
  if (!isTRUE(r$enabled))
    need(!any(vapply(SCENARIOS[cfg$scenarios], `[[`, logical(1), "resourceMatch")),
         "resource-matched scenarios requested but resourceMatching$enabled is FALSE")

  ## seeds
  need(is_count(cfg$seeds$repOffset, 0) || identical(cfg$seeds$repOffset, 0L), "seeds$repOffset")
  need(is.numeric(cfg$seeds$programmeOffset), "seeds$programmeOffset")

  if (length(problems))
    stop(paste(c("invalid configuration:", problems), collapse = "\n  "), call. = FALSE)
  invisible(cfg)
}

## The part of a configuration that defines a founder object. A saved founder file is
## only reused under a configuration whose founder part is identical.
founder_config <- function(cfg) {
  list(genome = cfg$genome,
       founders = cfg$founders[c("nSimulatedPerPool", "nKeptPerPool", "h2", "burnInCycles",
                                 "eliteMean", "eliteVarA")],
       burnIn = cfg$programme[c("nParents", "nCrosses", "nProgeny")])
}

## Number of plots in the elite nursery an arm grows from cycle 2 on.
nursery_size <- function(cfg, progeny) cfg$programme$nCrosses * progeny
