## =============================================================================
## outcomes.R -- the per-cycle record (DATA_DICTIONARY.md, "Per-cycle records")
##
## Population outcomes are measured at step 3 of Algorithm 1 on the unselected nursery
## X_t, before any admission. The admission columns of the same row describe what steps
## 4-5 of cycle t offered and admitted, and plotsGrown counts the plots of cycle t (the
## nursery X_t plus that cycle's introgression nurseries).
## =============================================================================

## Column order of every per-cycle CSV written by the engine.
RECORD_COLUMNS <- c(
  "cycle", "year", "meanG", "gain", "genicVar", "genicPct", "totalVarG", "accuracy",
  "donorPct", "pctSweet", "nAdmitted", "nCand", "nClearBar", "nClearTrue", "accDonor",
  "biasDonor", "admitTrueGap", "gapCycle", "nRedraw", "nF2screened", "nSweetUsable",
  "nProgeny", "plotsGrown", "strategy", "scenario", "visualSort", "resourceMatch",
  "rep", "seed", "RVersion", "AlphaSimR", "ompThreads", "locusA", "h2", "nAdmit",
  "nDonorCross", "realisedGap", "donorMean", "baseMean", "locusChr", "maxFreqDiff",
  "cullDonor")

## Columns that describe genetics and selection (everything a random draw can change).
## Two runs are genetically identical when these agree; plotsGrown is accounting only.
GENETIC_COLUMNS <- c("meanG", "gain", "genicVar", "genicPct", "totalVarG", "accuracy",
                     "donorPct", "pctSweet", "nAdmitted", "nCand", "nClearBar", "nClearTrue",
                     "accDonor", "biasDonor", "admitTrueGap", "gapCycle", "nF2screened",
                     "nSweetUsable", "nProgeny")

## Step 3: outcomes of the unselected nursery.
measure_nursery <- function(X, snpX, gebv, gv, founder, sp, t, cfg, progeny) {
  sc <- founder$scale
  V <- genic_variance(X, sp)
  meanG <- to_tha(mean(gv), sc)
  list(cycle = t,
       year = t * cfg$programme$yearsPerCycle,
       meanG = meanG,
       gain = meanG - founder$base$meanG,
       genicVar = var_to_tha(V, sc),
       genicPct = 100 * V / founder$base$genicVar,
       totalVarG = var_to_tha(pop_var(gv), sc),
       accuracy = safe_cor(gebv, gv),
       donorPct = donor_allele_pct(snpX, founder$diagnostic),
       pctSweet = 100 * mean(is_sweet(X, founder, sp)),
       gapCycle = meanG - founder$base$donorMean,
       nProgeny = progeny,
       nurseryPlots = nInd(X))
}

## One introgression nursery (or screen) of a cycle, for the plot accounting.
##   n         progeny screened (F2, second-step progeny or MAS seedlings)
##   nSweet    sweet homozygotes among them
##   byMarker  TRUE when the locus is scored on seedlings, so only sweet plants are grown
stage_record <- function(label, n, nSweet, byMarker) {
  list(label = label, n = as.integer(n), nSweet = as.integer(nSweet), byMarker = byMarker)
}

## Plots an introgression stage occupies: every progeny under field discard, only the
## sweet homozygotes under visual sorting or marker screening. Visual sorting changes this
## count and nothing else.
stage_plots <- function(stage, visualSort) {
  if (stage$byMarker || visualSort) stage$nSweet else stage$n
}

## Steps 4-5: the offered pool, the uncensored supply (equations 8 and 9), the predictor
## inside a donor-derived pool (equation 10) and the admitted lines' true merit.
summarise_offer <- function(offer, admitted, ctx, nurseryPlots) {
  gebv <- offer$gebv; gv <- offer$gv
  donorPool <- identical(offer$kind, "donor")
  screened <- sum(vapply(offer$stages, `[[`, integer(1), "n"))
  sweet <- sum(vapply(offer$stages, `[[`, integer(1), "nSweet"))
  intro <- sum(vapply(offer$stages, stage_plots, numeric(1), visualSort = ctx$scenario$visualSort))
  list(nAdmitted = length(admitted),
       nCand = length(gebv),
       nClearBar = sum(gebv >= ctx$b),
       nClearTrue = sum(gv >= ctx$bTrue),
       accDonor = if (donorPool) safe_cor(gebv, gv) else NA_real_,
       biasDonor = if (donorPool && length(gebv)) mean(gebv - gv) else NA_real_,
       admitTrueGap = if (length(admitted)) mean(gv[admitted]) - ctx$bTrue else NA_real_,
       nF2screened = as.integer(screened),
       nSweetUsable = as.integer(sweet),
       plotsGrown = as.integer(nurseryPlots + intro))
}

## Stack the cycle rows of one run and add the identification and founder columns.
assemble_records <- function(rows, founder, cfg, armName, scenario, seed, rep) {
  df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  meta <- list(
    strategy = armName, scenario = scenario$name,
    visualSort = scenario$visualSort, resourceMatch = scenario$resourceMatch,
    rep = rep, seed = seed,
    RVersion = software_versions()$R,
    AlphaSimR = software_versions()$AlphaSimR,
    ompThreads = Sys.getenv("OMP_NUM_THREADS", unset = NA_character_),
    locusA = cfg$sugar$additive, h2 = cfg$founders$h2,
    nAdmit = cfg$programme$cap, nDonorCross = cfg$introgression$nDonorCrosses,
    realisedGap = founder$base$realisedGap, donorMean = founder$base$donorMean,
    baseMean = founder$base$meanG, locusChr = founder$locus$chr,
    maxFreqDiff = founder$locus$maxFreqDiff, cullDonor = founder$cull$cullDonor,
    nRedraw = founder$nRedraw)
  for (k in names(meta)) df[[k]] <- meta[[k]]
  df[, RECORD_COLUMNS]
}
