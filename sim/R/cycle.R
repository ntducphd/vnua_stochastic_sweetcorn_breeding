## =============================================================================
## cycle.R -- one cycle of Algorithm 1
##
## Steps marked [R] draw from R's generator, in this order and no other:
##   1 [R] training candidates drawn from X_t; appended to the archive (latest records kept)
##   2     ridge regression on the archive, GEBV of X_t, the n_p parents, b_t (equation 5)
##   3     outcomes of the unselected X_t
##   4 [R] inflow arms only: donors and mates, crosses, selfs (meioses and phenotypes);
##         BRIDGE's second step; admission by equation (6)
##   5     elite controls only: the scheduled next-ranked elite candidates
##   6 [R] random mating of parents + admitted lines into X_{t+1} (meioses and phenotypes)
## Steps 2, 3 and 5 draw nothing, so every arm reaches step 4 of cycle 1 with the same
## generator state, training sample and parents.
## =============================================================================

## state: list(nursery = Pop, nurseryProgeny = int, archive = archive)
## env:   list(cfg, scenario, sp, founder)
run_cycle <- function(state, t, arm, env, trace = FALSE) {
  cfg <- env$cfg; sp <- env$sp; P <- cfg$programme
  X <- state$nursery
  tr <- if (trace) list(rngStart = rng_state()) else NULL

  ## 1 [R] training sample; its phenotypes were drawn when X_t was created
  trainIdx <- sample.int(nInd(X), P$nTrain)
  snpX <- pullSnpGeno(X, simParam = sp)
  state$archive <- archive_add(state$archive, snpX[trainIdx, , drop = FALSE],
                               X@pheno[trainIdx, 1])

  ## 2 model, GEBV, parents, threshold
  model <- fit_marker_model(state$archive)
  gebv <- predict_gebv(model, snpX)
  gv <- X@gv[, 1]
  rank <- top_n(gebv, length(gebv))
  parents <- rank[seq_len(P$nParents)]
  ctx <- list(t = t, cfg = cfg, scenario = env$scenario, sp = sp, founder = env$founder,
              nursery = X, gebv = gebv, gv = gv, rank = rank, parents = parents,
              b = mean(gebv[parents]), bTrue = mean(gv[parents]),
              predict = function(pop) predict_gebv(model, pullSnpGeno(pop, simParam = sp)))

  ## 3 outcomes before any admission
  rec <- measure_nursery(X, snpX, gebv, gv, env$founder, sp, t, cfg, state$nurseryProgeny)

  ## 4-5 offer and admission
  if (trace) tr$rngBeforeOffer <- rng_state()
  offer <- arm$offer(ctx)
  admitted <- arm$admit(offer, ctx)
  rec <- c(rec, summarise_offer(offer, admitted, ctx, rec$nurseryPlots))

  ## 6 [R] random mating of the parents and the admitted lines (parents first, by rank)
  if (trace) tr$rngBeforeMating <- rng_state()
  pool <- X[parents]
  if (length(admitted)) pool <- c(pool, offer$pop[admitted])
  state$nursery <- randCross(pool, nCrosses = P$nCrosses, nProgeny = arm$nurseryProgeny,
                             simParam = sp)
  state$nurseryProgeny <- arm$nurseryProgeny

  if (trace) {
    tr$trainIdx <- trainIdx
    tr$rank <- rank
    tr$parents <- parents
    tr$admittedIds <- if (length(admitted)) offer$pop@id[admitted] else character(0)
    tr$nextNurseryIds <- state$nursery@id
    tr$b <- ctx$b
  }
  list(state = state, record = rec, trace = tr)
}
