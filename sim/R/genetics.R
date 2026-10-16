## =============================================================================
## genetics.R -- reading genotypes and expressing genetic quantities
##
## The engine keeps every genetic value on AlphaSimR's internal scale and converts to
## t/ha only when a record is written (equation in Section 2.1.1):
##     t/ha = eliteMean + (g - g0) / sdA0,
## with g0 and sdA0^2 the mean and additive variance of the kept elite founders.
## =============================================================================

## Copies of the sweet allele (allele 1) at the endosperm locus, equation (1).
locus_dosage <- function(pop, founder, sp) {
  if (is.null(pop) || nInd(pop) == 0) return(integer(0))
  as.integer(pullMarkerGeno(pop, founder$locus$name, simParam = sp)[, 1])
}

## Marketability m_i = 1{X_i,l0 = 2}, equation (3). Never selected on.
is_sweet <- function(pop, founder, sp) locus_dosage(pop, founder, sp) == 2L

## Internal genetic value -> t/ha, and internal variance -> (t/ha)^2.
to_tha <- function(g, scale) scale$eliteMean + (g - scale$g0) / scale$sdA0
var_to_tha <- function(v, scale) v / scale$sdA0^2

## Genic variance, equation (11): sum over QTL of 2 p (1 - p) alpha^2 (internal scale).
genic_variance <- function(pop, sp) {
  Q <- pullQtlGeno(pop, simParam = sp)
  p <- colMeans(Q) / 2
  alpha <- sp$traits[[1]]@addEff
  sum(2 * p * (1 - p) * alpha^2)
}

## Names ("chr_site") of the loci of an AlphaSimR LociMap (trait or SNP chip).
loci_names <- function(lociMap) {
  paste0(rep(seq_along(lociMap@lociPerChr), lociMap@lociPerChr), "_", lociMap@lociLoc)
}

## -----------------------------------------------------------------------------
## Donor-allele content at the markers (record column donorPct).
## A marker is diagnostic when the post-burn-in elite nursery X_1 is fixed for one allele
## and the frozen donor pool carries the other. The content of a population is the summed
## frequency of the diagnostic (donor-only) alleles divided by their summed frequency in
## the donor pool, in percent: 0 for any population descending from X_1 alone, 100 for the
## donor pool itself, about 50 for donor x elite F2 (in expectation, before selection).
## -----------------------------------------------------------------------------
diagnostic_markers <- function(snpElite, snpDonor) {
  n <- nrow(snpElite)
  s <- colSums(snpElite)
  fixed0 <- s == 0
  fixed2 <- s == 2 * n
  qDonor1 <- colMeans(snpDonor) / 2                     # donor frequency of allele 1
  donorFreq <- ifelse(fixed0, qDonor1, 1 - qDonor1)     # frequency of the allele X_1 lacks
  keep <- (fixed0 | fixed2) & donorFreq > 0
  list(column = which(keep), allele = ifelse(fixed0[keep], 1L, 0L),
       donorFreq = donorFreq[keep])
}

donor_allele_pct <- function(snp, diag) {
  if (!length(diag$column)) return(NA_real_)
  p1 <- colMeans(snp[, diag$column, drop = FALSE]) / 2
  pDiag <- ifelse(diag$allele == 1L, p1, 1 - p1)
  100 * sum(pDiag) / sum(diag$donorFreq)
}
