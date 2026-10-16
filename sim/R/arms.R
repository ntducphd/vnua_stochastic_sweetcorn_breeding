## =============================================================================
## arms.R -- the programmes compared (Section 2.2) as strategy objects
##
## An arm is a list of class "sim_arm" with
##   name, kind      "closed", "inflow" or "elite"
##   offer(ctx)      step 4 (inflow) or 5 (elite control): the pool of candidates the arm
##                   offers, with GEBV and true values and the introgression stages grown
##   admit(off, ctx) which members of that pool join the parents (indices, in order)
##   nurseryProgeny  progeny per cross of the elite nursery the arm grows at step 6
## Arms differ only through these four parts; run_cycle() is common to all of them.
##
## The cycle context `ctx` (built in cycle.R) carries t, cfg, scenario, sp, founder, the
## nursery X_t, its GEBV and true values, the GEBV ranking `rank`, the parents (first n_p
## of the ranking), b_t, b*_t and predict(pop).
## =============================================================================

new_arm <- function(name, kind, offer, admit, nurseryProgeny) {
  structure(list(name = name, kind = kind, offer = offer, admit = admit,
                 nurseryProgeny = as.integer(nurseryProgeny)), class = "sim_arm")
}

print.sim_arm <- function(x, ...) {
  cat(sprintf("<arm %s (%s), nursery progeny per cross %d>\n", x$name, x$kind, x$nurseryProgeny))
  invisible(x)
}

## -----------------------------------------------------------------------------
## Offers
## -----------------------------------------------------------------------------
empty_offer <- function() list(kind = "none", pop = NULL, gebv = numeric(0), gv = numeric(0),
                               stages = list())

## A donor-derived pool: GEBV from the current model, true values from AlphaSimR.
donor_offer <- function(ctx, pool, stages) {
  if (is.null(pool) || nInd(pool) == 0)
    return(list(kind = "donor", pop = NULL, gebv = numeric(0), gv = numeric(0), stages = stages))
  list(kind = "donor", pop = pool, gebv = ctx$predict(pool), gv = pool@gv[, 1], stages = stages)
}

## [R] draw donors and elite mates (each with replacement), [R] cross them (one F1 per
## cross, meiosis and phenotypes) and [R] self every F1 into `progenyPerF1` progeny.
## Donor founders are the female parents and elite parents the males.
donor_cross_and_self <- function(ctx, nCrosses, progenyPerF1) {
  donors <- ctx$founder$donors
  donorIdx <- sample.int(nInd(donors), nCrosses, replace = TRUE)
  mateIdx <- ctx$parents[sample.int(length(ctx$parents), nCrosses, replace = TRUE)]
  f1 <- makeCross2(donors, ctx$nursery, crossPlan = cbind(donorIdx, mateIdx),
                   nProgeny = 1, simParam = ctx$sp)
  self(f1, nProgeny = progenyPerF1, simParam = ctx$sp)
}

sweet_subset <- function(pop, ctx) which(is_sweet(pop, ctx$founder, ctx$sp))

## DIRECT: donor x elite F2, field-grown; the sweet homozygotes are offered.
offer_direct <- function(f2PerCross) {
  function(ctx) {
    f2 <- donor_cross_and_self(ctx, ctx$cfg$introgression$nDonorCrosses, f2PerCross)
    s <- sweet_subset(f2, ctx)
    donor_offer(ctx, f2[s], list(stage_record("F2", nInd(f2), length(s), byMarker = FALSE)))
  }
}

## MAS: the same cross, but the locus is scored on `seedlingsPerCross` seedlings per F1
## and only the sweet homozygotes are grown; all of them are offered.
offer_mas <- function(seedlingsPerCross) {
  function(ctx) {
    f2 <- donor_cross_and_self(ctx, ctx$cfg$introgression$nDonorCrosses, seedlingsPerCross)
    s <- sweet_subset(f2, ctx)
    donor_offer(ctx, f2[s], list(stage_record("MAS seedlings", nInd(f2), length(s), byMarker = TRUE)))
  }
}

## BRIDGE: nothing from the first cross is offered. Its `families` sweet F2 of highest
## GEBV (all of them when fewer are sweet) are crossed to [R] elite mates drawn with
## replacement, [R] each F1 is selfed into `progeny` lines, and the sweet homozygotes of
## that second step (all of them: both parents are sweet homozygotes) are offered.
offer_bridge <- function(f2PerCross, families, progeny) {
  function(ctx) {
    f2 <- donor_cross_and_self(ctx, ctx$cfg$introgression$nDonorCrosses, f2PerCross)
    s1 <- sweet_subset(f2, ctx)
    stage1 <- stage_record("F2", nInd(f2), length(s1), byMarker = FALSE)
    if (!length(s1)) return(donor_offer(ctx, NULL, list(stage1)))
    sweetF2 <- f2[s1]
    keep <- top_n(ctx$predict(sweetF2), families)
    bridgeParents <- sweetF2[keep]
    nFam <- nInd(bridgeParents)
    mateIdx <- ctx$parents[sample.int(length(ctx$parents), nFam, replace = TRUE)]
    f1 <- makeCross2(bridgeParents, ctx$nursery, crossPlan = cbind(seq_len(nFam), mateIdx),
                     nProgeny = 1, simParam = ctx$sp)
    step2 <- self(f1, nProgeny = progeny, simParam = ctx$sp)
    s2 <- sweet_subset(step2, ctx)
    stage2 <- stage_record("bridge progeny", nInd(step2), length(s2), byMarker = FALSE)
    donor_offer(ctx, step2[s2], list(stage1, stage2))
  }
}

## Elite controls offer the elite nursery itself; their supply counts (equations 8, 9)
## are therefore taken over the whole nursery (Section 2.1.3).
offer_elite <- function(ctx) {
  list(kind = "elite", pop = ctx$nursery, gebv = ctx$gebv, gv = ctx$gv, stages = list())
}

## -----------------------------------------------------------------------------
## Admission rules
## -----------------------------------------------------------------------------
admit_none <- function(offer, ctx) integer(0)

## Equation (6): at most kappa candidates with GEBV >= b_t, largest first.
admit_threshold_cap <- function(gebv, b, cap) {
  ok <- which(gebv >= b)
  if (!length(ok) || cap <= 0) return(integer(0))
  ok[top_n(gebv[ok], cap)]
}

admit_by_threshold <- function(offer, ctx) {
  admit_threshold_cap(offer$gebv, ctx$b, ctx$cfg$programme$cap)
}

## Elite controls: the schedule[t] elite candidates ranked just below the parent cut,
## admitted without the threshold.
admit_next_ranked <- function(schedule) {
  force(schedule)
  function(offer, ctx) {
    k <- as.integer(schedule[ctx$t])
    if (is.na(k) || k <= 0L) return(integer(0))
    ctx$rank[ctx$cfg$programme$nParents + seq_len(k)]
  }
}

## -----------------------------------------------------------------------------
## Factory: the arm called `name` in `scenario`
## -----------------------------------------------------------------------------
## schedule: admissions per cycle, required for ELITE2M_<arm> (the inflow arm's nAdmitted
##           in the same replicate and scenario).
build_arm <- function(name, cfg, scenario, schedule = NULL) {
  P <- cfg$programme; I <- cfg$introgression; R <- cfg$resourceMatching
  rm <- isTRUE(scenario$resourceMatch)
  fullProgeny <- P$nProgeny
  armProgeny <- if (rm) R$nurseryProgeny else P$nProgeny   # inflow arms and their ELITE2M

  if (name == "CLOSED")
    return(new_arm(name, "closed", function(ctx) empty_offer(), admit_none, fullProgeny))

  if (name == "ELITE2")   # fills the cap in every cycle; grows no introgression nursery
    return(new_arm(name, "elite", offer_elite, admit_next_ranked(rep(P$cap, P$nCycles)),
                   fullProgeny))

  if (name == "DIRECT")
    return(new_arm(name, "inflow",
                   offer_direct(if (rm) R$directF2PerCross else I$f2PerCross),
                   admit_by_threshold, armProgeny))

  if (name == "MAS")
    return(new_arm(name, "inflow",
                   offer_mas(if (rm) R$masSeedlingsPerCross else I$masScreenFactor * I$f2PerCross),
                   admit_by_threshold, armProgeny))

  if (name == "BRIDGE")
    return(new_arm(name, "inflow",
                   offer_bridge(if (rm) R$bridgeF2PerCross else I$f2PerCross,
                                I$bridgeFamilies,
                                if (rm) R$bridgeProgeny else I$bridgeProgeny),
                   admit_by_threshold, armProgeny))

  if (startsWith(name, MATCHED_PREFIX)) {
    target <- substring(name, nchar(MATCHED_PREFIX) + 1)
    if (!target %in% INFLOW_ARMS) stop("ELITE2M must name an inflow arm: ", name)
    if (is.null(schedule) || length(schedule) != P$nCycles || anyNA(schedule) || any(schedule < 0))
      stop(sprintf("%s needs a schedule of %d non-negative admission counts", name, P$nCycles))
    return(new_arm(name, "elite", offer_elite, admit_next_ranked(as.integer(schedule)), armProgeny))
  }
  stop("unknown arm: ", name)
}
