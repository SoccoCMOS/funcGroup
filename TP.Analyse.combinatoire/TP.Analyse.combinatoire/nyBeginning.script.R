#==============================================================================
#==============================================================================
#
#     myBEGINNING.R = set of instructions to run at first
#
#==============================================================================
#==============================================================================


source("myCombinat.R")
source("mySeparating.R")


#==============================================================================
# Multiplicative decomposition of assembly performances
#==============================================================================
res    <- multiplicative.decomposition(mOccur, Fobs, rm.mono = TRUE)
size   <- res$size
mOccur <- res$mOccur
Fobs   <- res$Fobs
Fmono  <- res$Fmono
Fscale <- res$Fscale

fobs   <- Fobs/Fscale
alpha  <- res$alpha
beta   <- res$beta
trans  <- res$trans
nbAss  <- dim(mOccur)[1]
nbElt  <- dim(mOccur)[2]


nbclassesMax <- nbElt



#==============================================================================
# Compute the hierarchical tree
#==============================================================================
if (length(unique(xpr)) > 1) {

  alphaElt <- sdalphaElt <- betaElt <-
    sdbetaElt <- matrix(0, nrow = nbElt, ncol = length(setXpr))
  rownames(alphaElt) <- rownames(betaElt) <- colnames(mOccur)


  for (tim in seq_along(unique(xpr))) {
    indXpr <- which(xpr == unique(xpr)[tim])

    for (elt in seq_len(nbElt)) {
      index           <- which(mOccur[indXpr, elt] == 1)

      alphaElt[elt, tim]   <- gmean(alpha[index])
      sdalphaElt[elt, tim] <- gsd(alpha[index])

      betaElt[elt, tim]    <- amean(beta[index])
      sdbetaElt[elt, tim]  <- asd(beta[index])
    }
  }

} else {

  alphaElt <- betaElt <- fobsElt <-
    sdalphaElt <- sdbetaElt <- sdfobsElt <- numeric(nbElt)
  names(alphaElt) <- names(betaElt) <- colnames(mOccur)

  for (elt in seq_len(nbElt)) {
    index           <- which(mOccur[ , elt] == 1)

    alphaElt[elt]   <- gmean(alpha[index])
    sdalphaElt[elt] <- gsd(alpha[index])

    betaElt[elt]    <- amean(beta[index])
    sdbetaElt[elt]  <- asd(beta[index])
  }
  fobsElt   <- alphaElt * betaElt
  sdfobsElt <- (sdalphaElt / alphaElt + sdbetaElt / betaElt) * fobsElt
}

dot       <- cbind(alphaElt, betaElt)
dot.dist  <- dist(x = dot, method = "euclidean", diag = FALSE, upper = FALSE)
Rtree     <- hclust(d = dot.dist, method = "ward.D")


tree.div   <- as.mytree(Rtree)
mAffectElt <- mAffect.elements(Rtree)
mAssMotifs <- mAffect.motifs(tree.div, mOccur)


#==============================================================================
#==============================================================================
#
#     END of FILE
#
#==============================================================================
#==============================================================================
