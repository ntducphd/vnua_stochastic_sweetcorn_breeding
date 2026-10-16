## =============================================================================
## founders.R -- the founder object F of a replicate (Section 2.1.5)
##
## F comprises the coalescent draw, the choice of the endosperm locus and the cull at it,
## the draw of QTL, markers and effects, and the burn-in. build_founders() makes it under
## its own seed; load_or_build_founders() reuses a saved file when there is one, because
## under AlphaSimR 2.1.0 runMacs does not reproduce from a seed and a saved founder file
## therefore *defines* the replicate. The programme seed is set only after F exists, so
## a programme run is the same whether F was just built or read from disk.
##
## Order of the draws in build_founders() (all from R's generator except the coalescent):
##   1 runMacs (R supplies one seed per chromosome; MaCS itself is not seed-controlled)
##   2 locus choice and cull (deterministic), repeated from 1 if a pool has too few survivors
##   3 QTL, effects and markers (addTraitA, addSnpChip), redrawn while either set holds the locus
##   4 phenotypes of the kept elite founders, then per burn-in cycle: selection on phenotype
##     (deterministic) and random mating (meioses and the new generation's phenotypes)
## =============================================================================

FOUNDER_CLASS <- "sweetcorn_founders"

## -----------------------------------------------------------------------------
## Building blocks
## -----------------------------------------------------------------------------

## runMacs with split = s places the two pools in the first and second half of the
## individuals. The endosperm locus is the segregating site with the largest absolute
## difference in allele frequency between the halves (first such site on ties). The elite
## pool is the half in which allele 1 is the more frequent, so that allele 1 is the sweet
## allele as equation (1) requires.
choose_endosperm_locus <- function(haplo, nPerPool) {
  G <- pullSegSiteGeno(haplo)
  a <- seq_len(nPerPool); b <- nPerPool + seq_len(nPerPool)
  pA <- colMeans(G[a, , drop = FALSE]) / 2
  pB <- colMeans(G[b, , drop = FALSE]) / 2
  diff <- abs(pA - pB)
  j <- which.max(diff)
  eliteIsA <- pA[j] >= pB[j]
  name <- colnames(G)[j]
  list(name = name,
       chr = as.integer(sub("_.*$", "", name)),
       site = as.integer(sub("^.*_", "", name)),
       index = j,
       maxFreqDiff = unname(diff[j]),
       freqSweetElite = unname(if (eliteIsA) pA[j] else pB[j]),
       freqSweetDonor = unname(if (eliteIsA) pB[j] else pA[j]),
       eliteRows = if (eliteIsA) a else b,
       donorRows = if (eliteIsA) b else a,
       dosage = as.integer(G[, j]))
}

## Elite founders must be homozygous sweet (dosage 2), donor founders homozygous
## non-sweet (dosage 0). The first survivors in each pool are kept.
cull_at_locus <- function(locus, nKeep) {
  eliteOK <- locus$eliteRows[locus$dosage[locus$eliteRows] == 2L]
  donorOK <- locus$donorRows[locus$dosage[locus$donorRows] == 0L]
  n <- length(locus$eliteRows)
  list(ok = length(eliteOK) >= nKeep && length(donorOK) >= nKeep,
       eliteKeep = utils::head(eliteOK, nKeep),
       donorKeep = utils::head(donorOK, nKeep),
       nEliteSurvivors = length(eliteOK), nDonorSurvivors = length(donorOK),
       cullElite = 1 - length(eliteOK) / n,
       cullDonor = 1 - length(donorOK) / n)
}

## Draw the QTL (with additive effects scaled to mean 0 and variance 1 in the combined
## founder population) and then the markers from the sites not chosen as QTL, so the two
## sets are disjoint. A fresh SimParam is used for every draw; the draw is repeated until
## neither set contains the endosperm locus. nRedraw counts the rejected draws.
draw_qtl_and_markers <- function(haplo, genome, locusName, maxRedraws) {
  nRedraw <- 0L
  repeat {
    sp <- SimParam$new(haplo)
    sp$nThreads <- 1L
    sp$addTraitA(nQtlPerChr = genome$nQtlPerChr, mean = 0, var = 1)
    sp$addSnpChip(nSnpPerChr = genome$nSnpPerChr)
    qtl <- loci_names(sp$traits[[1]])
    snp <- loci_names(sp$snpChips[[1]])
    if (!(locusName %in% qtl) && !(locusName %in% snp)) break
    nRedraw <- nRedraw + 1L
    if (nRedraw > maxRedraws) stop(sprintf("no valid QTL/marker draw in %d redraws", maxRedraws))
  }
  list(sp = sp, nRedraw = nRedraw, qtl = qtl, snp = snp)
}

## Burn-in by phenotypic selection in the elite pool only: in each cycle the n_p
## individuals of highest phenotype are mated at random into n_c crosses of n_o progeny,
## whose phenotypes AlphaSimR draws as they are created. The last generation is the
## cycle-1 elite nursery X_1.
burn_in <- function(elite, cycles, programme, sp) {
  pop <- setPheno(elite, simParam = sp)
  for (b in seq_len(cycles)) {
    sel <- top_n(pop@pheno[, 1], programme$nParents)
    pop <- randCross(pop[sel], nCrosses = programme$nCrosses,
                     nProgeny = programme$nProgeny, simParam = sp)
  }
  pop
}

## -----------------------------------------------------------------------------
## The founder object
## -----------------------------------------------------------------------------
build_founders <- function(cfg, seed, verbose = TRUE) {
  g <- cfg$genome; fo <- cfg$founders
  nSim <- fo$nSimulatedPerPool
  t0 <- Sys.time()
  seed_rng(seed)

  ## 1-2  coalescent draw, endosperm locus, cull
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    haplo <- runMacs(nInd = 2L * nSim, nChr = g$nChr, segSites = g$segSites,
                     species = g$species, split = g$splitGenerations, nThreads = 1L)
    locus <- choose_endosperm_locus(haplo, nSim)
    cull <- cull_at_locus(locus, fo$nKeptPerPool)
    if (cull$ok) break
    if (verbose) log_msg("founders %d: cull left %d elite / %d donor survivors (need %d); new coalescent draw",
                         seed, cull$nEliteSurvivors, cull$nDonorSurvivors, fo$nKeptPerPool)
    if (attempt >= fo$maxCoalescentAttempts)
      stop(sprintf("founders %d: too few survivors of the cull in %d coalescent draws", seed, attempt))
  }

  ## 3  QTL, effects, markers
  draw <- draw_qtl_and_markers(haplo, g, locus$name, fo$maxRedraws)
  sp <- draw$sp
  ## The trait is purely additive, so additive and genetic variance coincide; both use the
  ## divisor N, the convention under which AlphaSimR scaled the founders to variance 1.
  combined <- newPop(haplo, simParam = sp)            # no residual variance yet: no draw
  varGFounders <- pop_var(combined@gv[, 1])
  varE <- varGFounders * (1 - fo$h2) / fo$h2           # held fixed from here on
  sp$setVarE(varE = varE)                              # phenotypes drawn at every new population
  elite <- combined[cull$eliteKeep]
  donors <- combined[cull$donorKeep]

  ## t/ha scale from the kept elite founders
  scale <- list(eliteMean = fo$eliteMean,
                g0 = mean(elite@gv[, 1]),
                sdA0 = sqrt(pop_var(elite@gv[, 1]) / fo$eliteVarA),
                varE = varE, varGFounders = varGFounders)

  ## 4  burn-in
  nursery <- burn_in(elite, fo$burnInCycles, cfg$programme, sp)

  ## base quantities shared by every arm and scenario of the replicate
  snpNursery <- pullSnpGeno(nursery, simParam = sp)
  snpDonor <- pullSnpGeno(donors, simParam = sp)
  baseG <- mean(nursery@gv[, 1])
  varNursery <- pop_var(nursery@gv[, 1])
  base <- list(
    meanG = to_tha(baseG, scale),
    genicVar = genic_variance(nursery, sp),
    donorMean = to_tha(mean(donors@gv[, 1]), scale),
    sdA = sqrt(var_to_tha(varNursery, scale)),
    h2Cycle1 = varNursery / (varNursery + varE))
  base$realisedGap <- base$meanG - base$donorMean

  f <- list(
    seed = seed,
    configName = cfg$name,
    config = founder_config(cfg),
    software = software_versions(),
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    buildSeconds = seconds_since(t0),
    SP = sp,
    nursery = nursery,
    nurseryProgeny = cfg$programme$nProgeny,
    donors = donors,
    locus = locus[c("name", "chr", "site", "index", "maxFreqDiff", "freqSweetElite",
                    "freqSweetDonor")],
    qtl = draw$qtl, snp = draw$snp,
    nRedraw = draw$nRedraw,
    nCoalescentDraws = attempt,
    cull = cull[c("cullElite", "cullDonor", "nEliteSurvivors", "nDonorSurvivors")],
    keptRows = list(elite = cull$eliteKeep, donor = cull$donorKeep),
    eliteFounderGv = elite@gv[, 1],                    # internal scale, before burn-in
    scale = scale,
    base = base,
    diagnostic = diagnostic_markers(snpNursery, snpDonor))
  class(f) <- FOUNDER_CLASS
  if (verbose) log_msg("founders %d: built in %.0f s (locus %s, %d redraws, donor cull %.1f%%, gap %.2f t/ha)",
                       seed, f$buildSeconds, f$locus$name, f$nRedraw, 100 * f$cull$cullDonor,
                       f$base$realisedGap)
  f
}

## One-line description of a founder object.
print.sweetcorn_founders <- function(x, ...) {
  cat(sprintf("<founders seed %d, %s: locus %s (chr %d, max freq diff %.3f), %d redraws, cull elite %.1f%% donor %.1f%%,\n  base %.2f t/ha, donor %.2f t/ha, gap %.2f t/ha; AlphaSimR %s, R %s>\n",
              x$seed, x$configName, x$locus$name, x$locus$chr, x$locus$maxFreqDiff, x$nRedraw,
              100 * x$cull$cullElite, 100 * x$cull$cullDonor, x$base$meanG, x$base$donorMean,
              x$base$realisedGap, x$software$AlphaSimR, x$software$R))
  invisible(x)
}

## -----------------------------------------------------------------------------
## Saving, reloading and checksums
## -----------------------------------------------------------------------------
founder_path <- function(dir, seed) file.path(dir, sprintf("founder_%d.rds", seed))

## Write atomically (temporary file, then rename) with a SHA-256 sidecar, and refresh
## the directory's checksum list.
save_founders <- function(f, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".partial")
  saveRDS(f, tmp)
  if (!file.rename(tmp, path)) stop("could not move ", tmp, " to ", path)
  writeLines(sprintf("%s  %s", sha256_file(path), basename(path)), paste0(path, ".sha256"))
  write_founder_checksums(dirname(path))
  invisible(path)
}

sha256_file <- function(path) unname(tools::sha256sum(path))

## SHA256SUMS lists every founder file of a directory ("<hash>  <file>").
write_founder_checksums <- function(dir) {
  files <- sort(list.files(dir, pattern = "^founder_.*\\.rds$"))
  lines <- if (length(files)) sprintf("%s  %s", sha256_file(file.path(dir, files)), files) else character(0)
  tmp <- file.path(dir, "SHA256SUMS.partial")
  writeLines(lines, tmp)
  file.rename(tmp, file.path(dir, "SHA256SUMS"))
  invisible(lines)
}

## TRUE when every file in SHA256SUMS still has its recorded hash.
verify_founder_checksums <- function(dir) {
  sums <- file.path(dir, "SHA256SUMS")
  if (!file.exists(sums)) return(NA)
  x <- read.table(sums, colClasses = "character", col.names = c("hash", "file"))
  all(sha256_file(file.path(dir, x$file)) == x$hash)
}

## Reload a founder file and check that it belongs to this seed and configuration and that
## its content matches its checksum sidecar.
read_founders <- function(path, cfg = NULL, seed = NULL) {
  side <- paste0(path, ".sha256")
  if (file.exists(side)) {
    rec <- strsplit(readLines(side, n = 1), "  ", fixed = TRUE)[[1]][1]
    if (!identical(rec, sha256_file(path))) stop("checksum mismatch for ", path)
  }
  f <- readRDS(path)
  if (!inherits(f, FOUNDER_CLASS)) stop(path, " is not a founder object of this engine")
  if (!is.null(seed) && !identical(as.numeric(f$seed), as.numeric(seed)))
    stop(sprintf("%s holds seed %s, not %s", path, f$seed, seed))
  if (!is.null(cfg) && !identical(f$config, founder_config(cfg)))
    stop(sprintf("%s was built under a different founder configuration than '%s'", path, cfg$name))
  f
}

## The replicate's founder object: from disk when saved, otherwise built and saved.
load_or_build_founders <- function(cfg, seed, dir = NULL, verbose = TRUE) {
  if (!is.null(dir)) {
    path <- founder_path(dir, seed)
    if (file.exists(path)) {
      f <- read_founders(path, cfg, seed)
      attr(f, "source") <- "file"
      if (verbose) log_msg("founders %d: read %s", seed, path)
      return(f)
    }
  }
  f <- build_founders(cfg, seed, verbose = verbose)
  if (!is.null(dir)) save_founders(f, founder_path(dir, seed))
  attr(f, "source") <- "built"
  f
}
