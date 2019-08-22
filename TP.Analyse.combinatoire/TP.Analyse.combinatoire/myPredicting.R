#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                               myPREDICTING.R
#
#                         Benoit JAILLARD, summer 2016
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


source("myStats.R")


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#       List of functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# predict.amean.bymot          <- function(assMotifs, mOccur, fct)
# predict.xpr.amean.bymot      <- function(assMotifs, mOccur, fct, xpr)
#
# amean.bymot.LOO              <- function(fctMot)
# predict.amean.bymot.LOO      <- function(assMotifs, mOccur, fct)
# predict.xpr.amean.bymot.LOO  <- function(assMotifs, mOccur, fct, xpr)
#
# amean.bymot.jack             <- function(fctMot, jack)
# predict.amean.bymot.jack     <- function(assMotifs, mOccur, fct, jack)
# predict.xpr.amean.bymot.jack <- function(assMotifs, mOccur, fct, jack, xpr)
#
# predict.gmean.bymot          <- function(assMotifs, mOccur, fct)
# predict.xpr.gmean.bymot      <- function(assMotifs, mOccur, fct, xpr)
#
# gmean.bymot.LOO              <- function(fctMot)
# predict.gmean.bymot.LOO      <- function(assMotifs, mOccur, fct)
# predict.xpr.gmean.bymot.LOO  <- function(assMotifs, mOccur, fct, xpr)
#
# gmean.bymot.jack             <- function(fctMot, jack)
# predict.gmean.bymot.jack     <- function(assMotifs, mOccur, fct, jack)
# predict.xpr.gmean.bymot.jack <- function(assMotifs, mOccur, fct, jack, xpr)
#
#
# predict.amean.byelt          <- function(assMotifs, mOccur, fct)
# predict.xpr.amean.byelt      <- function(assMotifs, mOccur, fct, xpr)
#
# amean.byelt.LOO              <- function(fctMot)
# predict.amean.byelt.LOO      <- function(assMotifs, mOccur, fct)
# predict.xpr.amean.byelt.LOO  <- function(assMotifs, mOccur, fct, xpr)
#
# amean.byelt.jack             <- function(fctMot, jack)
# predict.amean.byelt.jack     <- function(assMotifs, mOccur, fct, jack)
# predict.xpr.amean.byelt.jack <- function(assMotifs, mOccur, fct, jack, xpr)
#
# predict.gmean.byelt          <- function(assMotifs, mOccur, fct)
# predict.xpr.gmean.byelt      <- function(assMotifs, mOccur, fct, xpr)
#
# gmean.byelt.LOO              <- function(fctMot)
# predict.gmean.byelt.LOO      <- function(assMotifs, mOccur, fct)
# predict.xpr.gmean.byelt.LOO  <- function(assMotifs, mOccur, fct, xpr)
#
# gmean.byelt.jack             <- function(fctMot, jack)
# predict.gmean.byelt.jack     <- function(assMotifs, mOccur, fct, jack)
# predict.xpr.gmean.byelt.jack <- function(assMotifs, mOccur, fct, jack, xpr)
#
#
# predict.cal     <- function(assMotifs, mOccur, fct, xpr,
#                             opt.mean = "amean",
#                             opt.mod  = "bymot")
#
# predict.prd     <- function(assMotifs, mOccur, fct, xpr,
#                             opt.mean = "amean",
#                             opt.mod  = "bymot",
#                             opt.jack = FALSE,   jack = c(2, 5))
#
#
# compute.bias     <- function(tree.prd, fct)
# compute.mStats   <- function(mCal, mPrd, fct, nbK)
# compute.tree.Stats <- function(mCal, mPrd, mStats, fct, nbK)
#
#
# predict.function <- function(tree.cal, mOccur, fct, xpr,
#                              opt.mean = "amean",
#                              opt.mod  = "byelt",
#                              opt.jack = FALSE,  jack = c(2,5))
#
# predict.twin     <- function(tree.cal, mOccur, alpha, beta, xpr,
#                              fscale    = 1,
#                              titre     = "",
#                              opt.alpha = "gmean",
#                              opt.beta  = "amean",
#                              opt.mod   = "byelt",
#                              opt.jack  = FALSE,    jack = c(2,5))
#
# plot.mStats      <- function(mStats, nbElt, titre = "")
#
#
# plot.prediction.simple <- function(Fprd, Fobs, assMotifs, nbcl,
#                                    titre    = "",
#                                    opt.mean = "amean",
#                                    opt.aov  = FALSE,
#                                    pvalue   = 0.025)
#
# plot.prediction.LOO    <- function(Fcal, Fprd, Fobs, assMotifs, nbcl,
#                                    titre    = "",
#                                    opt.mean = "amean",
#                                    opt.aov  = FALSE,
#                                    pvalue   = 0.025)
#
# add.ass.names          <- function(Fcal, Fprd, Fobs,
#                                    assMotifs, AssNames)
#
# plot.prediction        <- function(res,
#                                    xpr      = rep(1, length(fct)),
#                                    titre    = "",
#                                    opt.aov  = FALSE,   pvalue = 0.025,
#
#                                    clu.stat = FALSE,
#                                    clu.cal  = FALSE,
#                                    clu.prd  = FALSE,
#                                    tre.stat = FALSE,
#                                    tre.prd  = FALSE,
#                                    tre.pub  = FALSE,
#                                    tre.best = FALSE,
#                                    tre.opt  = TRUE,
#                                    ass.loc  = FALSE,
#                                    opt.all  = FALSE)
#



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#             CROSS-VALIDATION PREDICTIONS BY MOTIF
#                    (Non-linear model only)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      bymot = by motif in a whole (WITHOUT taking into account
#                                                       species contribution)
#      by including all assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.bymot <- function(assMotifs, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments (xpr)
#      amean = by using arithmetic mean
#      bymot = by motif in a whole (WITHOUT taking into account
#                                                         species contribution)
#      by including all assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.bymot <- function(assMotifs, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.amean.bymot(assMotifs[indXpr], fct[indXpr] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      LOO   = by excluding the assembly to predict    (Leave One Out)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean.bymot.LOO <- function(fctMot) {

  nbass <- length(fctMot)

  if (nbass > 1) {
    fctPrd <- numeric(nbass)
    for (ind in seq_len(nbass)) fctPrd[ind] <- amean(fctMot[-ind])
  } else {
    fctPrd <- NA
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.bymot.LOO <- function(assMotifs, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean.bymot.LOO(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments (xpr)
#      amean = by using arithmetic mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      LOO   = by excluding the assembly to predict    (Leave One Out)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.bymot.LOO <- function(assMotifs, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.amean.bymot.LOO(assMotifs[indXpr], fct[indXpr] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      Jack   = by excluding jack[1] assemblies to predict    (Jackknife)
#               jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean.bymot.jack <- function(fctMot, jack) {

  nbass <- length(fctMot)

  if (nbass > jack[1] * jack[2]) {

    fctPrd <- numeric(nbass)
    index  <- sample.int(nbass)
    size   <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fctPrd[indjack] <- amean(fctMot[-indjack])
    }

    indjack         <- index[(ind * size + 1):nbass]
    fctPrd[indjack] <- amean(fctMot[-indjack])

  } else {

    fctPrd <- amean.bymot.LOO(fctMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.bymot.jack <- function(assMotifs, fct, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean.bymot.jack(fct[indMot], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments (xpr)
#      amean = by using arithmetric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      Jack   = by excluding jack[1] assemblies to predict    (Jackknife)
#               jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.bymot.jack <- function(assMotifs, fct, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.amean.bymot.jack(assMotifs[indXpr],
                                               fct[indXpr], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      by including all the assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.bymot <- function(assMotifs, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      by including all the assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.bymot <- function(assMotifs, fct, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.gmean.bymot(assMotifs[indXpr], fct[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      LOO   = by excluding the assembly to predict    (Leave One Out)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean.bymot.LOO <- function(fctMot) {

  nbass <- length(fctMot)

  if (nbass > 1) {
    fctPrd <- numeric(nbass)
    for (ind in seq_len(nbass)) fctPrd[ind] <- gmean(fctMot[-ind])
  } else {
    fctPrd <- NA
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.bymot.LOO <- function(assMotifs, fct) {

  fctPrd    <- numeric(length(assMotifs))

  setMot    <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean.bymot.LOO(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      LOO   = by excluding the assembly to predict    (Leave One Out)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.bymot.LOO <- function(assMotifs, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.gmean.bymot.LOO(assMotifs[indXpr], fct[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      Jack   = by excluding jack[1] assemblies to predict    (Jackknife)
#               jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean.bymot.jack <- function(fctMot, jack) {

  nbass <- length(fctMot)

  if (nbass > jack[1] * jack[2]) {

    fctPrd <- numeric(nbass)
    index  <- sample.int(nbass)
    size   <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fctPrd[indjack] <- gmean(fctMot[-indjack])
    }

    indjack         <- index[(ind * size + 1):nbass]
    fctPrd[indjack] <- gmean(fctMot[-indjack])

  } else {

    fctPrd <- gmean.bymot.LOO(fctMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.bymot.jack <- function(assMotifs, fct, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean.bymot.jack(fct[indMot], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      xpr   = by using several experiments
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      Jack  = by excluding jack[1] assemblies to predict    (Jackknife)
#              jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.bymot.jack <- function(assMotifs, fct, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.gmean.bymot.jack(assMotifs[indXpr],
                                               fct[indXpr], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#          CROSS-VALIDATION PREDICTION BY ELEMENTS INSIDE EACH MOTIF
#                 (Non-linear model + local linear model)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Old versions

predict.amean.byelt.old1 <- function(assMotifs, mOccur, fct) {

  fctPrd  <- numeric(length(assMotifs))

  setMot  <- unique(assMotifs)
  mfct <- matrix(NA, nrow = length(setMot), ncol = dim(mOccur)[2])
  #  colnames(mfct) <- colnames(mOccur)
  #  rownames(mfct) <- setMot

  for (mot in seq_along(setMot)) {

    motif  <- setMot[mot]
    indMot <- which(assMotifs == motif)
    setElt <- unique(which((mOccur[indMot, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(mOccur[indMot, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- amean(fct[indMot[indElt]])
    }
  }


  msk     <- logical(length(fct))
  sizeAss <- apply(mOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setMot)) {

    motif  <- setMot[mot]
    indMot <- which(assMotifs == motif)

    if (length(indMot) > 0) {
      setElt <- unique(which((mOccur[indMot, , drop = FALSE] == 1),
                             arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setElt)) {
        element <- setElt[elt]
        indElt  <- which(mOccur[indMot, element] == 1)

        if (length(indElt) > 0) {
          index <- indMot[indElt]
          fctPrd[index] <- fctPrd[index] + mfct[motif, element]
        }
      }

      fctPrd[indMot] <- fctPrd[indMot] / sizeAss[indMot]
      msk[indMot]     <- TRUE
    }
  }

  fctPrd[msk == FALSE] <- NA

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.byelt.old2 <- function(assMotifs, mOccur, fct) {

  fctPrd  <- numeric(length(assMotifs))

  setMot  <- unique(assMotifs)
  sizeAss <- apply(mOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setMot)) {
    indMot <- which(assMotifs == setMot[mot])

    setElt <- unique(which((mOccur[indMot, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])
    for (elt in seq_along(setElt)) {
      indElt <- which(mOccur[indMot, setElt[elt]] == 1)

      if (length(indElt) > 0) {
        index         <- indMot[indElt]
        fctPrd[index] <- fctPrd[index] + amean(fct[index])
      }
    }

    fctPrd[indMot] <- fctPrd[indMot] / sizeAss[indMot]
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean.byelt <- function(fctMot, mOccurMot) {

  fctPrd <- numeric(length(fctMot))

  setElt <- unique(which((mOccurMot[ , , drop = FALSE] == 1),
                          arr.ind = TRUE)[ , 2])

  for (elt in seq_along(setElt)) {
    indElt         <- which(mOccurMot[ , setElt[elt]] == 1)
    fctPrd[indElt] <- fctPrd[indElt] + amean(fctMot[indElt])
  }

  fctPrd <- fctPrd / apply(mOccurMot, MARGIN = 1, FUN = sum)

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.byelt <- function(assMotifs, mOccur, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])

    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean.byelt(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- fct[indMot]
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.byelt <- function(assMotifs, mOccur, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.amean.byelt(assMotifs[indXpr],
                                          mOccur[indXpr, ],
                                          fct[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      LOO   = by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# old version

predict.amean.byelt.LOO.old <- function(assMotifs, mOccur, fct) {

  fctPrd    <- numeric(length(assMotifs))
  fctPrd[ ] <- NA

  vfct <- numeric(dim(mOccur)[2])
  #  names(mfct) <- colnames(mOccur)

  setMot    <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])

    if (length(indMot) > 1)
      for (ind in seq_along(indMot)) {

        vfct[] <- NA
        indOth <- indMot[-ind]
        setElt <- which(mOccur[indMot[ind], ] != 0)

        for (elt in seq_along(setElt)) {
          element <- setElt[elt]
          indElt  <- which(mOccur[indOth, element] == 1)
          if (length(indElt) > 0)
               vfct[element] <- amean(fct[indOth[indElt]])
        }

        fctPrd[indMot[ind]] <- amean(vfct[setElt])
      }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean.byelt.LOO <- function(fctMot, mOccurMot) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccurMot[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])

#      if (length(indElt) > 0) {
#        vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])
#    } else {
#        vfct[setElt[elt]] <- amean(fctMot[indOth])
#    }

    }

    fctPrd[ind] <- amean(vfct[setElt])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.byelt.LOO <- function(assMotifs, mOccur, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {
    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean.byelt.LOO(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      LOO   = by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.byelt.LOO <- function(assMotifs, mOccur, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict.amean.byelt.LOO(assMotifs[indXpr],
                                              mOccur[indXpr, ],
                                              fct[indXpr] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      Jack  = by excluding jack[1] assemblies to predict    (Jackknife)
#              jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# old version

amean.byelt.jack.old <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack <- index[(ind - 1) * size + (1:size)]

      vfct[]  <- NA
      indOth  <- seq_len(nbass)[-indjack]
      setElt  <- unique(which((mOccurMot[indjack, , drop = FALSE] == 1),
                               arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setElt)) {
        indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
        if (length(indElt) > 0)
          vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])
      }

      fctPrd[indjack] <- amean(vfct[setElt])
    }

    indjack <- index[(ind * size + 1):nbass]

    vfct[]  <- NA
    indOth  <- seq_len(nbass)[-indjack]
    setElt  <- unique(which((mOccurMot[indjack, , drop = FALSE] == 1),
                            arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])
    }

    fctPrd[indjack] <- amean(vfct[setElt])

  } else {

    fctPrd[] <- amean.byelt.LOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean.byelt.jack <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

      tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

    tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)

  } else {

    fctPrd[ ] <- amean.byelt.LOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.byelt.jack <- function(assMotifs, mOccur, fct, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean.byelt.jack(fct[indMot], mOccur[indMot, ], jack)
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      Jack  = by excluding jack[1] assemblies to predict    (Jackknife)
#              jack[2] is the number of quartiles
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.amean.byelt.jack <- function(assMotifs, mOccur, fct, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict.amean.byelt.jack(assMotifs[index],
                                              mOccur[index, ],
                                              fct[index],  jack  )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Old version

predict.gmean.byelt.old <- function(assMotifs, mOccur, fct) {

  fctPrd    <- numeric(length(assMotifs))
  fctPrd[ ] <- 1

  setMot <- unique(assMotifs)
  mfct   <- matrix(NA, nrow = length(setMot), ncol = dim(mOccur)[2])
#  colnames(mfct) <- colnames(mOccur)
#  rownames(mfct) <- setMot

  for (mot in seq_along(setMot)) {
    motif  <- setMot[mot]
    indMot <- which(assMotifs == motif)

    setElt <- unique(which((mOccur[indMot, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element              <- setElt[elt]
      indElt               <- which(mOccur[indMot, element] == 1)
      mfct[motif, element] <- gmean(fct[indMot[indElt]])
    }
  }


  msk        <- logical(length(fct))
  sizeAss    <- apply(mOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setMot)) {
    motif  <- setMot[mot]
    indMot <- which(assMotifs == motif)

    setElt <- unique(which((mOccur[indMot, , drop = FALSE] == 1),
                            arr.ind = TRUE)[ , 2])

    if (length(indMot) > 0) {
      for (elt in seq_along(setElt)) {
        element <- setElt[elt]
        indElt  <- which(mOccur[indMot, element] == 1)

        if (length(indElt) > 0) {
          index         <- indMot[indElt]
          fctPrd[index] <- fctPrd[index] * mfct[motif, element]
        }
      }

      fctPrd[indMot] <- fctPrd[indMot] ^ (1/sizeAss[indMot])
      msk[indMot]    <- TRUE
    }
  }

  fctPrd[msk == FALSE] <- NA

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean.byelt <- function(fctMot, mOccurMot) {

  fctPrd    <- numeric(length(fctMot))
  fctPrd[ ] <- 1

  setElt    <- unique(which((mOccurMot[ , , drop = FALSE] == 1),
                             arr.ind = TRUE)[ , 2])
  for (elt in seq_along(setElt)) {
    indElt <- which(mOccurMot[ , setElt[elt]] == 1)

    if (length(indElt) > 0)
      fctPrd[indElt] <- fctPrd[indElt] * gmean(fctMot[indElt])
  }

  fctPrd <- fctPrd ^ (1/apply(mOccurMot, MARGIN = 1, FUN = sum))

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.byelt <- function(assMotifs, mOccur, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean.byelt(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- fct[indMot]
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.byelt <- function(assMotifs, mOccur, fct, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict.gmean.byelt(assMotifs[index],
                                         mOccur[index, ],
                                         fct[index] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      LOO = by excluding the assembly to predict    (Leave One Out)
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Old version

predict.gmean.byelt.LOO.old <- function(assMotifs, mOccur, fct) {

  fctPrd    <- numeric(length(assMotifs))
  fctPrd[ ] <- NA

  mfct <- numeric(dim(mOccur)[2])
  #  names(mfct) <- colnames(mOccur)

  setMot    <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])

    if (length(indMot) > 1)
      for (ind in seq_along(indMot)) {

        indOth <- indMot[-ind]
        setElt <- which(mOccur[indMot[ind], ] != 0)

        mfct[ ] <- NA
        for (elt in seq_along(setElt)) {

          element <- setElt[elt]
          indElt  <- which(mOccur[indOth, element] == 1)
          if (length(indElt) > 0)
            mfct[element] <- gmean(fct[indOth[indElt]])
        }

        fctPrd[indMot[ind]] <- gmean(mfct[setElt])
      }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean.byelt.LOO <- function(fctMot, mOccurMot) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccurMot[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- gmean(fctMot[indOth[indElt]])

#      if (length(indElt) > 0) {
#        vfct[setElt[elt]] <- gmean(fctMot[indOth[indElt]])
#    } else {
#        vfct[setElt[elt]] <- gmean(fctMot[indOth])
#    }

    }

    fctPrd[ind] <- gmean(vfct[setElt])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.byelt.LOO <- function(assMotifs, mOccur, fct) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean.byelt.LOO(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      LOO = by excluding the assembly to predict    (Leave One Out)
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.byelt.LOO <- function(assMotifs, mOccur, fct, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr    <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict.gmean.byelt.LOO(assMotifs[index],
                                             mOccur[index, ],
                                             fct[index] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      Jack  = by excluding jack[1] assemblies to predict    (Jackknife)
#              jack[2] is the number of quartiles
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# old version

gmean.byelt.jack.old <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack <- index[(ind - 1) * size + (1:size)]

      vfct[] <- NA
      indOth <- seq_len(nbass)[-indjack]
      setElt <- unique(which((mOccurMot[indjack, , drop = FALSE] == 1),
                              arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setElt)) {
        indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
        if (length(indElt) > 0)
          vfct[setElt[elt]] <- gmean(fctMot[indOth[indElt]])
      }

      fctPrd[indjack] <- gmean(vfct[setElt])
    }

    indjack <- index[(ind * size + 1):nbass]

    vfct[] <- NA
    indOth <- seq_len(nbass)[-indjack]
    setElt <- unique(which((mOccurMot[indjack, , drop = FALSE] == 1),
                            arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- gmean(fctMot[indOth[indElt]])
    }

    fctPrd[indjack] <- gmean(vfct[setElt])

  } else {

    fctPrd[] <- gmean.byelt.LOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean.byelt.jack <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

      tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

    tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)

  } else {

    fctPrd[] <- gmean.byelt.LOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.byelt.jack <- function(assMotifs, mOccur, fct, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean.byelt.jack(fct[indMot], mOccur[indMot, ], jack)
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed
#      x     = by using several experiments
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      Jack  = by excluding jack[1] assemblies to predict    (Jackknife)
#              jack[2] is the number of quartiles
#      for any Fonction fct (Fobs for instance)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.xpr.gmean.byelt.jack <- function(assMotifs, mOccur, fct, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict.gmean.byelt.jack(assMotifs[index],
                                              mOccur[index, ],
                                              fct[index],  jack  )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#             GENERAL FUNCTIONS OF PREDICTIon BY CROSS-VALIDATION
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed by including all the assemblies to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.cal <- function(assMotifs, mOccur, fct, xpr,
                        opt.mean = "amean",
                        opt.mod  = "bymot") {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.mod == "byelt") optmod <- "byelt"

  optxpr <- ""
  if (length(unique(xpr)) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot. =
             predict.amean.bymot(assMotifs, fct),
           amean.bymot.xpr =
             predict.xpr.amean.bymot(assMotifs, fct, xpr),
           gmean.bymot. =
             predict.gmean.bymot(assMotifs, fct),
           gmean.bymot.xpr =
             predict.xpr.gmean.bymot(assMotifs, fct, xpr),

           amean.byelt. =
             predict.amean.byelt(assMotifs, mOccur, fct),
           amean.byelt.xpr =
             predict.xpr.amean.byelt(assMotifs, mOccur, fct, xpr),
           gmean.byelt. =
             predict.gmean.byelt(assMotifs, mOccur, fct),
           gmean.byelt.xpr =
             predict.xpr.gmean.byelt(assMotifs, mOccur, fct, xpr)  )
    )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed by excluding (LOO) the assembly to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.prd <- function(assMotifs, mOccur, fct, xpr,
                        opt.mean = "amean",
                        opt.mod  = "bymot",
                        opt.jack = FALSE,   jack = c(2, 5)) {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.mod == "byelt") optmod <- "byelt"

  optjack <- ""
  if (opt.jack == TRUE) optjack <- "jack"

  optxpr <- ""
  if (length(unique(xpr)) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optjack, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot.. =
             predict.amean.bymot.LOO(assMotifs, fct),
           amean.bymot..xpr =
             predict.xpr.amean.bymot.LOO(assMotifs, fct, xpr),
           amean.bymot.jack. =
             predict.amean.bymot.jack(assMotifs, fct, jack),
           amean.bymot.jack.xpr =
             predict.xpr.amean.bymot.jack(assMotifs, fct, jack, xpr),

           gmean.bymot.. =
             predict.gmean.bymot.LOO(assMotifs, fct),
           gmean.bymot..xpr =
             predict.xpr.gmean.bymot.LOO(assMotifs, fct, xpr),
           gmean.bymot.jack. =
             predict.gmean.bymot.jack(assMotifs, fct, jack),
           gmean.bymot.jack.xpr =
             predict.xpr.gmean.bymot.jack(assMotifs, fct, jack, xpr),

           amean.byelt.. =
             predict.amean.byelt.LOO(assMotifs, mOccur, fct),
           amean.byelt..xpr =
             predict.xpr.amean.byelt.LOO(assMotifs, mOccur, fct, xpr),
           amean.byelt.jack. =
             predict.amean.byelt.jack(assMotifs, mOccur, fct, jack),
           amean.byelt.jack.xpr =
             predict.xpr.amean.byelt.jack(assMotifs, mOccur, fct, jack, xpr),

           gmean.byelt.. =
             predict.gmean.byelt.LOO(assMotifs, mOccur, fct),
           gmean.byelt..xpr =
             predict.xpr.gmean.byelt.LOO(assMotifs, mOccur, fct, xpr),
           gmean.byelt.jack. =
             predict.gmean.byelt.jack(assMotifs, mOccur, fct, jack),
           gmean.byelt.jack.xpr =
             predict.xpr.gmean.byelt.jack(assMotifs, mOccur, fct, jack, xpr) )
    )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#           FUNCTIONS OF COMPUTATION OF VARIOUS PARAMETERS
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute the matrix 3 x nbclusters of R2cal, R2prd, and missing values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute.mStats <- function(mCal, mPrd, fct, nbK) {

  nbclMax <- dim(mCal)[1]

  tmp    <- c("missing", "R2cal", "R2prd", "pcal", "pprd", "AIC", "AICc")
  mStats <- matrix(0, nrow = nbclMax, ncol = length(tmp))
  rownames(mStats) <- seq_len(nbclMax)
  colnames(mStats) <- tmp

  for (nbcl in seq_len(nbclMax)) {
    mStats[nbcl, "missing"] <-
      length(na.action(na.omit(mPrd[nbcl, ]))) / length(mCal[nbcl, ])
    mStats[nbcl, "R2cal"] <- R2mse(mCal[nbcl, ], fct)
    mStats[nbcl, "R2prd"] <- R2mse(mPrd[nbcl,],  fct)
    mStats[nbcl, "AIC"]   <- AIC(  mCal[nbcl, ], fct, nbK[nbcl])
    mStats[nbcl, "AICc"]  <- AICc( mCal[nbcl, ], fct, nbK[nbcl])

    if (nbK[nbcl] > 1) {
      mStats[nbcl, "pcal"] <- pmse( mCal[nbcl, ], fct, nbK[nbcl])
      mStats[nbcl, "pprd"] <- pmse( mPrd[nbcl,],  fct, nbK[nbcl])
    }
  }

  return(mStats)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Compute and plot the best Prediction
#                        as a concatenation of the best predictions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute.tree.Stats <- function(mCal, mPrd, mStats, fct, nbK) {

  nbclMax <- dim(mCal)[1]
  nbAss   <- dim(mCal)[2]

  tCal <- tPrd <- tNbcl <- matrix(NA, nrow = nbclMax, ncol = nbAss)
  rownames(tCal) <- rownames(tPrd) <- rownames(tNbcl) <- rownames(mCal)
  colnames(tCal) <- colnames(tPrd) <- colnames(tNbcl) <- colnames(mCal)

  for (nbcl in seq_len(nbclMax)) {
    for (ass in seq_len(nbAss)) {
      last  <- min(length(na.omit(mPrd[1:nbcl, ass])),
                   length(na.omit(mStats[1:nbcl, "R2prd"])),
                   nbcl)
      index <- first.argmax(mStats[1:last, "R2prd"])

      tNbcl[nbcl, ass] <- index
      tPrd[nbcl, ass]  <- mPrd[index, ass]
      tCal[nbcl, ass]  <- mCal[index, ass]
    }
  }
  tStats <- compute.mStats(tCal, tPrd, fct, nbK)


  if (dim(tStats)[1] > 1) for (nbcl in 2:nbclMax)
    if (tStats[nbcl, "R2prd"] <= tStats[nbcl - 1, "R2prd"]) {
      tStats[nbcl, "R2prd"] <- tStats[nbcl - 1, "R2prd"]
      tNbcl[nbcl, ]         <- tNbcl[nbcl - 1, ]
      tPrd[nbcl, ]          <- tPrd[nbcl - 1, ]
      tCal[nbcl, ]          <- tCal[nbcl - 1, ]
  }
  tStats <- compute.mStats(tCal, tPrd, fct, nbK)

  res        <- list(tCal, tPrd, tStats, tNbcl)
  names(res) <- c("tCal", "tPrd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Compute the statistiques of each motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute.motif.Stats <- function(tCal, tPrd, mMotifs) {

  nbclMax <- dim(tCal)[1]

  # build the general structure of motif.Stats

  uTab <- list()
  for (nbcl in 1:nbclMax) uTab[[nbcl]] <- table(mMotifs[nbcl, ])


  # compute the values of motif.Stats

  uMean <- uSd <- uRmse <- uR2 <- uSlope <- uTab

  for (nbcl in 1:nbclMax)
    for (motif in names(uTab[[nbcl]])) {

      index <- which(mMotifs[nbcl, ] == motif)

      uMean[[nbcl]][motif] <- amean(tPrd[nbcl, index])

      if (length(index) > 1) {

        uSd[[nbcl]][motif]    <- asd(tPrd[nbcl, index])
        uRmse[[nbcl]][motif]  <- rmse(tPrd[nbcl, index], tCal[nbcl, index])
        uR2[[nbcl]][motif]    <- R2mse(tPrd[nbcl, index], tCal[nbcl, index])
        uSlope[[nbcl]][motif] <-
          lm(tPrd[nbcl, index] ~ tCal[nbcl, index])$coef[2]

      } else {

        oldMotif <- mMotifs[nbcl - 1, index]

        uSd[[nbcl]][motif]    <- uSd[[nbcl - 1]][oldMotif]
        uRmse[[nbcl]][motif]  <- uRmse[[nbcl - 1]][oldMotif]
        uR2[[nbcl]][motif]    <- uR2[[nbcl - 1]][oldMotif]
        uSlope[[nbcl]][motif] <- uSlope[[nbcl - 1]][oldMotif]
      }
    }

  res        <- list(uTab, uMean, uSd, uRmse, uR2, uSlope)
  names(res) <- c("uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute assembly functioning and associated statistiques
#          using the Clustering model
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.function <- function(tree.cal, mOccur, fct,
                             xpr = rep(1, length(fct)),
                             # options for computing
                             opt.mean = "amean",
                             opt.mod  = "byelt",
                             opt.jack = FALSE,  jack = c(2,5)) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  names(fct) <- names(xpr) <- rownames(mOccur)

  nbElt      <- dim(mOccur)[2]
  nbAss      <- length(fct)

  mAssMotifs <- mAffect.motifs(tree.cal, mOccur)

  # Compute the cross-validation predictions
  mCal <- mPrd <- matrix(NA, nrow = nbElt, ncol = nbAss)
  rownames(mCal) <- rownames(mPrd) <- seq_len(nbElt)
  colnames(mCal) <- colnames(mPrd) <- rownames(mOccur)

  for (nbcl in seq_len(nbElt)) {
    assMotifs    <- mAssMotifs[nbcl, ]
    mCal[nbcl, ] <-
      predict.cal(assMotifs, mOccur, fct, xpr, opt.mean, opt.mod)
    mPrd[nbcl, ] <-
      predict.prd(assMotifs, mOccur, fct, xpr,
                  opt.mean, opt.mod, opt.jack, jack)
  }

  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mAssMotifs[nbcl, ]))

  mStats   <- compute.mStats(mCal, mPrd, fct, nbK)

  # Compute the Global Predictions for all the number of clusters
  tree.prd <- compute.tree.Stats(mCal, mPrd, mStats, fct, nbK)

  uStats   <- compute.motif.Stats(tree.prd$tCal, tree.prd$tPrd, mAssMotifs)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
              mCal, mPrd, mStats, mAssMotifs,
              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl,
              uStats$uTab, uStats$uMean, uStats$uSd,
              uStats$uRmse, uStats$uR2, uStats$uSlope)


  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
                  "mCal", "mPrd", "mStats", "mMotifs",
                  "tCal", "tPrd", "tStats", "tNbcl",
                  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute assembly functioning and associated statistiques
#          using the Clustering model
#          for Interaction effet (alpha), Composition effect (beta),
#            and assembly Functioning as Alpha * Beta * Fscale
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.twin <- function(tree.cal, mOccur, alpha, beta,
                         xpr       = rep(1, length(fct)),
                         fscale    = 1,
                         titre     = "",
                         opt.alpha = "gmean",
                         opt.beta  = "amean",
                         opt.mod   = "byelt",
                         opt.jack  = FALSE,    jack = c(2,5) ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  fct     <- alpha * beta * fscale

  mAssMotifs <- mAffect.motifs(tree.cal, mOccur)

  nbElt   <- dim(mOccur)[2]
  nbAss   <- length(fct)
  nbElt <- dim(mAssMotifs)[1]

  # Compute the cross-validation predictions
  mCal <- mPrd <- matrix(NA, nrow = nbElt, ncol = nbAss)
  rownames(mCal) <- rownames(mPrd) <- seq_len(nbElt)

  for (nbcl in seq_len(nbElt)) {
    assMotifs <- mAssMotifs[nbcl, ]

    mCal[nbcl, ] <-
      predict.cal(assMotifs, mOccur, alpha, xpr, opt.alpha, opt.mod) *
      predict.cal(assMotifs, mOccur, beta, xpr, opt.beta, opt.mod) *
      fscale

    mPrd[nbcl, ] <-
      predict.prd(assMotifs, mOccur, alpha, xpr,
                  opt.alpha, opt.mod, opt.jack, jack) *
      predict.prd(assMotifs, mOccur, beta, xpr,
                  opt.beta, opt.mod, opt.jack, jack) *
      fscale
  }

  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mAssMotifs[nbcl, ]))
  mStats <- compute.mStats(mCal, mPrd, fct, nbK)

  # Compute the Global Predictions for all the number of clusters
  tree.prd <- compute.tree.Stats(mCal, mPrd, mStats, fct, nbK)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  opt.mean <- paste(opt.beta)

  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
              mCal, mPrd, mStats, mAssMotifs,
              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl)

  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
                  "mCal", "mPrd", "mStats", "mMotifs",
                  "tCal", "tPrd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                      PLOTTING FUNCTIONS
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot the matrix 3 x nbclusters of R2cal, R2prd, and missing values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.mStats <- function(mStats, nbElt, titre = "") {

  nbclMax <- dim(mStats)[1]

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a first graph without any pvaluesx
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-

  plot(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
       xlim = c(1, nbElt),   ylim = c(0,1),
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "R2 (in black), E (in red)",
       xlab = "number of clusters",
       main = titre)

  # plot R2 calibration
  points(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
         type = "b", pch = 1, cex = 2, bg = "white", col = "black")

  # plot R2 prediction with pvalues
  points(x = seq_len(nbclMax), y = mStats[ ,"R2prd"],
         type = "b", pch = 1, cex = 2, bg = "white", col = "red3")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a second graph without any pvalues
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  plot(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
       xlim = c(1, nbElt), ylim = c(0,1),
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "Predicting ratio",
       xlab = "number of clusters",
       main = titre)

  # plot predictiong ratio
  points(x = seq_len(nbclMax), y = 1 - mStats[ ,"missing"],
         type = "b", pch = 0, cex = 2, bg = "white", col = "blue")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a third graph with AIC and AICc
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  ylim <- c(min(mStats[ , "AIC"], mStats[ , "AICc"], na.rm = TRUE),
            max(mStats[1, "AIC"], mStats[1, "AICc"], na.rm = TRUE))

  # plot AIC
  plot(x = seq_len(nbclMax), y = mStats[ ,"AIC"],
       xlim = c(1, nbElt), ylim = ylim,
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "AIC", xlab = "number of clusters",
       main = titre)

  fct  <- na.omit(mStats[ , "AIC"])
  abline(v = first.argmin(fct)[1], col = "green3")
  points(x = seq_along(fct), y = fct,
         type = "b", pch = 1, cex = 2, bg = "white", col = "blue3")

  # plot AICc
  plot(x = seq_len(nbclMax), y = mStats[ ,"AIC"],
       xlim = c(1, nbElt), ylim = ylim,
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "AICc", xlab = "number of clusters",
       main = titre)

  fct  <- na.omit(mStats[ , "AICc"])
  abline(v = first.argmin(fct)[1], col = "green3")
  points(x = seq_along(fct), y = fct,
         type = "b", pch = 1, cex = 2, bg = "white", col = "red3")

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot Predictions vs Observations
#
#  Inputs : Fprd : vector of modelled or predicted values
#           Fobs : vector of observed values
#           assMotifs : vector of motifs of which belong the assemblies
#
#  Options : titre : main titre of Figure
#            opt.reg : technical information on the quality of the regression
#            opt.aov : variance analysis of assembly function by motif
#            pvalue : threshold for the aov analysis
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction.simple <- function(Fprd, Fobs, assMotifs, nbcl,
                                   titre    = "",
                                   xylim    = range(Fobs),
                                   opt.mean = "amean",
                                   opt.aov  = FALSE,
                                   pvalue   = 0.025) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs and plot the figure
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbAss <- length(Fobs)
  tmp   <- mean.fct(Fobs, opt.mean)

  index <- na.action(na.omit(Fprd))
  if (length(index) == nbAss) stop("the vector Fprd cannot be null")

  if (length(index) > 0) {
    Fprd      <- Fprd[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
  }

  plot(x = Fobs, xlab = "Observations", xlim = xylim,
       y = Fprd, ylab = "Predictions",  ylim = xylim,
       main = titre,
       type = "n", tck = 0.02, las = 1)

  abline(v = tmp, h = tmp, lty = "dotted", col = "blue")

  lines(x = xylim, y = xylim, lty = "solid", col = "red")

  points(x = Fobs, y = Fprd,
         pch = figures[assMotifs], col = couleurs[assMotifs],
         bg = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Adding of various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75*width
  posX2 <- xylim[1] + 1.02*width
  posX3 <- xylim[1] - 0.02*width
  posX4 <- xylim[1] + 0.02*width
  posX5 <- xylim[1] + 0.25*width

  posY1 <- xylim[1] + 0.05*width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02*width
  posY4 <- xylim[1] + 0.98*width


  # predicting ratio
  tmp <- paste("predicting = ", length(na.omit(Fprd)), "/", nbAss, sep = "")
  text(x = posX1, y = posY1, labels = tmp, col = "red")

  # R2 value
  if ((nbAss - length(index)) > 1) {
    tmp <- paste("R2 = ", signif(R2mse(Fprd, Fobs), digits = 3), sep = "")
    text(x = posX1, y = posY2, labels = tmp, col = "red")
  }

  # Number of clusters
  text(paste("Nb clusters = ", nbcl, sep = ""),
       x = posX5, y = posY4, col = "red")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Optional informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # option variance analysis

  if (opt.aov == TRUE) {
    setMot <- sort(unique(assMotifs))

    text(x = posX2, y = posY3, labels = "pred", col = "black")

    test <- test.posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(rownames(test) == motif)
        text(x = posX2, y = test[index, "means"],
             labels = as.character(test[index, "groups"]),
             col =  couleurs[motif], font = 3)
      }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot Predictions vs Observations
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction.LOO <- function(Fcal, Fprd, Fobs, assMotifs, nbcl,
                                xylim    = range(Fobs),
                                titre    = "",
                                opt.mean = "amean",
                                opt.aov  = FALSE,
                                pvalue   = 0.025) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs and plot the figure
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbAss <- length(Fobs)
#  xylim <- c(0.7, 2.5) # pour alpha
#  xylim <- c(0.45, 1.73) # pour beta
  tmp   <- mean.fct(Fobs, opt.mean)

  index <- na.action(na.omit(Fprd))
  if (length(index) == nbAss) return(FALSE)

  if (length(index) > 0) {
    Fprd      <- Fprd[-index]
    Fcal      <- Fcal[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
  }

  plot(x = Fobs, xlab = "Observations", xlim = xylim,
       y = Fcal, ylab = "Predictions",  ylim = xylim,
       main = titre,
       las = 1, type = "n", tck = 0.02)

  for (elt in seq_along(Fobs))
    lines(x = c(Fobs[elt], Fobs[elt]), y = c(Fprd[elt], Fcal[elt]),
          col = couleurs[assMotifs[elt]], lty = "solid")

  abline(v = tmp, h = tmp, lty = "dotted", col = "blue")

  lines(x = xylim, y = xylim, lty = "solid", col = "red")

  points(x = Fobs, y = Fcal,
         pch = figures[assMotifs], col = couleurs[assMotifs],
         bg = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Adding of various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75*width
  posX2 <- xylim[1] + 1.02*width
  posX3 <- xylim[1] - 0.02*width
  posX4 <- xylim[1] + 0.02*width
  posX5 <- xylim[1] + 0.25*width

  posY1 <- xylim[1] + 0.05*width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02*width
  posY4 <- xylim[1] + 0.98*width


  # Predicting assemblies

  tmp <- paste("predicting = ", length(na.omit(Fprd)), "/", nbAss, sep = "")
  text(x = posX1, y = posY1, labels = tmp, col = "red")

  #  Number of clusters and R2 value

  if ((nbAss - length(index)) > 1) {
    tmp <- paste0("R2 = ",  signif(R2mse(Fcal, Fobs), digits = 3),
                 "  E = ", signif(R2mse(Fprd, Fobs), digits = 3))
    text(x = posX1, y = posY2, labels = tmp, col = "red")

    tmp <- paste0("Nb clusters = ", nbcl, "  E/R2 = ",
                  signif(R2mse(Fprd, Fobs)/R2mse(Fcal, Fobs), digits = 3))
    text(tmp, x = posX5, y = posY4, col = "red")
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Optional informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # option variance analysis

  if (opt.aov == TRUE) {
    setMot <- sort(unique(assMotifs))

    text(x = posX3, y = posY3, labels = "cal", col = "black")
    test   <- test.posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(rownames(test) == motif)
        text(x = posX3, y = test[index, "means"],
             labels = as.character(test[index, "groups"]),
             col =  couleurs[motif], font = 3)
      }

    text(x = posX2, y = posY3, labels = "prd", col = "black")
    test   <- test.posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(rownames(test) == motif)
        text(x = posX2, y = test[index, "means"],
             labels = as.character(test[index, "groups"]),
             col =  couleurs[motif], font = 3)
      }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Add the names of assemblies
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

add.ass.names <- function(Fcal, Fprd, Fobs,
                          assMotifs, AssNames) {

  index <- na.action(na.omit(Fprd))
  if (length(index) == length(Fobs)) return(FALSE)

  if (length(index) != 0) {
    Fprd      <- Fprd[-index]
    Fcal      <- Fcal[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
    AssNames  <- AssNames[-index]
  }

  thre   <- max(Fobs) - (max(Fobs) - min(Fobs)) / 5
  index1 <- which(Fobs <= thre)
  index2 <- which(Fobs >  thre)

  if (length(index1) != 0)
    text(x = Fobs[index1], y = Fcal[index1], labels = AssNames[index1],
         col = couleurs[assMotifs[index1]], pos = 4)          # to the right of

  if (length(index2) != 0)
    text(x = Fobs[index2], y = Fcal[index2], labels = AssNames[index2],
         col = couleurs[assMotifs[index2]], pos = 2)          # to the left of
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot all the necessary predicting figures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction <- function(res,
                            xpr.who  = 5,
                            nbOpt    = 0,

                            titre    = "",

                            clu.stat = FALSE,
                            clu.cal  = FALSE,
                            clu.prd  = FALSE,

                            tre.stat = FALSE,
                            tre.prd  = FALSE,
                            tre.pub  = FALSE,
                            tre.best = FALSE,
                            tre.opt  = TRUE,
                            tre.calvsprd = FALSE,

                            ass.loc  = FALSE,

                            opt.aov  = FALSE, pvalue = 0.025,

                            opt.all  = FALSE  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.all == TRUE) {
    clu.stat <- clu.cal <- clu.prd <-
    tre.stat <- tre.prd <- tre.pub <- tre.best <- tre.opt <-
    tre.calvsprd <- ass.loc <- TRUE
  }

  opt.xpr <- FALSE
  setXpr  <- unique(res$xpr)
  indxpr  <- seq_along(setXpr)
  if (length(indxpr) > 1) {
    opt.xpr <- TRUE
    setAss  <- res$names[which(res$xpr == setXpr[1])]
    pas     <- length(setAss)
  }

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Plot analysis cluster by cluster
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt <- dim(res$mStats)[1]
  nbAss <- length(res$fct)

  nbMax <- first.argmax(res$tStats[ ,"R2prd"])
  if (nbOpt == 0) { nbOpt <- first.argmin(res$tStats[ ,"AICc"])
  } else {
    if (nbOpt > nbMax) nbOpt <- nbMax
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Stats
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.stat == TRUE) {
    titre2 <- paste("R2 calibration and prediction:", titre,
                    paste(res$opt.mean, res$opt.mod, sep = "/"), sep = " ")
    plot.mStats(res$mStats, nbElt, titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Calibrations for all the number of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.cal == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Calibration with", nbcl, "clusters", sep = " ")
      plot.prediction.simple(res$mCal[nbcl, ], res$fct,
                             res$mMotifs[nbcl, ], nbcl,
                             xylim    = range(res$fct),
                             opt.mean = res$opt.mean,
                             opt.aov  = opt.aov,
                             pvalue   = pvalue,
                             titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Predictions for all the number of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.prd == TRUE)
    for (nbcl in seq_len(nbMax))
      if (length(na.omit(res$mPrd[nbcl, ])) > 0) {
        titre2 <- paste(titre, "Prediction with", nbcl, "clusters", sep = " ")
        plot.prediction.LOO(res$mCal[nbcl, ], res$mPrd[nbcl, ], res$fct,
                            res$mMotifs[nbcl, ], nbcl,
                            opt.mean = res$opt.mean,
                            opt.aov  = opt.aov,
                            pvalue   = pvalue,
                            titre    = titre2)
      }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Stats of Tree predictions
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.stat == TRUE) {
    titre2 <- paste("R2 tree calibration and tree prediction:", titre,
                    paste(res$opt.mean, res$opt.mod, sep = "/"), sep = " ")
    plot.mStats(res$tStats, nbElt, titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #    using different colors for predictions with different cluster numbers
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.prd == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          res$tNbcl[nbcl, ], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      index <- which(res$tNbcl[nbcl, ] == (nbcl - 1))
      if (length(index) != 0)
        add.ass.names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index],
                      res$tNbcl[nbcl, index], res$names[index])
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #          with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.pub == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          rep(1, nbAss), nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the best Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.best == TRUE) {
    nbcl   <- nbMax
    titre2 <- paste(titre, "Best Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum (according to AICc) Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions versus Tree Calibrations for all numbers of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.calvsprd == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction vs Calibration with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.simple(res$tPrd[nbcl, ], res$tCal[nbcl, ],
                             res$mMotifs[nbcl, ], nbcl,
                             xylim    = range(res$fct),
                             opt.mean = res$opt.mean,
                             opt.aov  = opt.aov,
                             pvalue   = pvalue,
                             titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by experiment with the names of Assemblages
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.xpr == TRUE) {
    nbcl   <- nbOpt

    for (ipr in seq_along(indxpr)) {
      titre2 <- paste0(titre,
                       " Tree Predictions for the Experiment ", indxpr[ipr])
      index  <- which(res$xpr == setXpr[ipr])
      plot.prediction.LOO(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                          res$fct[index], res$mMotifs[nbcl, index],
                          nbcl, xylim = range(res$fct),
                          opt.mean = res$opt.mean,
                          titre    = titre2)
      if (ass.loc == TRUE)
        add.ass.names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index], res$mMotifs[nbcl, index],
                      res$names[index])
    }
  } else {
    if (ass.loc == TRUE) {
      nbcl   <- nbOpt
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          res$mMotifs[nbcl, ], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      add.ass.names(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                    res$mMotifs[nbcl, ], res$names)
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by Assemblage over experiments
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.xpr == TRUE) {

    nbcl   <- nbOpt
    if (is.numeric(xpr.who) == TRUE) {
        setIndex <- sort(sample(1:pas, size = xpr.who[1]))
    } else {
        setIndex <- NULL
        for (elt in seq_along(xpr.who))
          setIndex <- c(setIndex, which(res$names == xpr.who[elt])[1])
    }

    for (ass in setIndex) {
      titre2 <-
        paste0(titre, " Tree Prediction of the Assemblage '",
               setAss[ass], "'")
      index <- ass + ((1:length(indxpr)) - 1) * pas
      plot.prediction.LOO(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                          res$fct[index], res$mMotifs[nbcl, index],
                          nbcl, xylim = range(res$fct),
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      if (ass.loc == TRUE)
        add.ass.names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index], res$mMotifs[nbcl, index],
                      res$xpr[index])
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

write.prediction <- function(filename, res) {

#  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
#                  "mCal", "mPrd", "mStats", "mMotifs",
#                  "tCal", "tPrd", "tStats", "tNbcl",
#                  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")
#
#  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
#              mCal, mPrd, mStats, mAssMotifs,
#              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl,
#              uStats$uTab, uStats$uMean, uStats$uSd,
#              uStats$uRmse, uStats$uR2, uStats$uSlope)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  nbElt <- dim(res$mStats)[1]
  nbAss <- length(res$fct)

  # nbElt, nbAss, "opt.mean", "opt.mod": 1 word
  tmp           <- matrix(c(nbElt, nbAss, res$opt.mean, res$opt.mod),
                          nrow = 1, byrow = FALSE)
  colnames(tmp) <- c("nbElt", "nbAss", "opt.mean", "opt.mod")
  write.table(x = tmp,
              file = paste(filename, "pred.str", "csv", sep = "."),
              append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

  # "names", "fct", "xpr" : nbAss words
  tmp           <- rbind(res$names, res$fct, res$xpr)
  write.table(x = tmp,
              file = paste(filename, "pred.vct", "csv", sep = "."),
              append = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")

  # "mCal", "mPrd", "mMotifs", "tCal", "tPrd", "tNbcl" : nbElt x nbAss
  tmp <- rbind(res$mCal, res$mPrd, res$mMotifs,
               res$tCal, res$tPrd, res$tNbcl)
  colnames(tmp) <- res$names
  write.table(x = tmp,
              file = paste(filename, "pred.mat", "csv", sep = "."),
              append = FALSE, col.names = FALSE, row.names = FALSE, sep = ",")

  # "mStats", "tStats" : nbElt x nbStat
  tmp <- rbind(res$mStats, res$tStats)
  write.table(x = tmp,
              file = paste(filename, "pred.stat", "csv", sep = "."),
              append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read.prediction <- function(filename) {

  #  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
  #                  "mCal", "mPrd", "mStats", "mMotifs",
  #                  "tCal", "tPrd", "tStats", "tNbcl",
  #                  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")
  #
  #  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
  #              mCal, mPrd, mStats, mAssMotifs,
  #              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl,
  #              uStats$uTab, uStats$uMean, uStats$uSd,
  #              uStats$uRmse, uStats$uR2, uStats$uSlope)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read options
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp  <- read.table(file = paste(filename, "pred.str", "csv", sep = "."),
                     header = TRUE, sep = ",")

  nbElt    <- as.integer(  tmp[1, 1])
  nbAss    <- as.integer(  tmp[1, 2])
  opt.mean <- as.character(tmp[1, 3])
  opt.mod  <- as.character(tmp[1, 4])


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read vectors
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  tmp  <- read.table(file = paste(filename, "pred.vct", "csv", sep = "."),
                     header = FALSE, sep = ",")

  name. <- as.character(as.vector(unlist(tmp[1, ])))
  fct   <- as.numeric(  as.vector(unlist(tmp[2, ])))
  xpr   <- as.character(as.vector(unlist(tmp[3, ])))
  names(fct) <- names(xpr) <- name.


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read matrices
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  row      <- 1

  tmp  <- read.table(file = paste(filename, "pred.mat", "csv", sep = "."),
                     header = FALSE, sep = ",")

  mCal     <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  mPrd     <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  mMotifs  <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  tCal     <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  tPrd     <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  tNbcl    <- as.matrix(tmp[row:(row + nbElt - 1), ])

  colnames(mCal) <- colnames(mPrd) <- colnames(mMotifs) <- name.
  colnames(tCal) <- colnames(tPrd) <- colnames(tNbcl)   <- name.

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Read Stats
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  row      <- 1

  tmp  <- read.table(file = paste(filename, "pred.stat", "csv", sep = "."),
                     header = TRUE, sep = ",")

  mStats   <- as.matrix(tmp[row:(row + nbElt - 1), ])
  row      <- row + nbElt

  tStats   <- as.matrix(tmp[row:(row + nbElt - 1), ])


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # A FINIR
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  res <- list(name., fct, xpr, opt.mean, opt.mod,
              mCal, mPrd, mStats, mMotifs,
              tCal, tPrd, tStats, tNbcl)
#  ,
#              uTab, uMean, uSd, uRmse, uR2, uSlope)

  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
                  "mCal", "mPrd", "mStats", "mMotifs",
                  "tCal", "tPrd", "tStats", "tNbcl")
#  ,
#                  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#           FUNCTIONS OF PREDICTION of A SUPPLEMENTARY DATASET
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      amean = by using arithmetic mean
#      bymot = by motif in a whole (WITHOUT taking into account
#                                                       species contribution)
#      by including all assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.amean.bymot.sup <- function(supMotifs, appMotifs, appFct) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- amean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      by including all the assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.gmean.bymot.sup <- function(supMotifs, appMotifs, appFct) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- gmean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# version 1

predict.amean.byelt.sup.old <- function(supMotifs, supOccur,
                                        appMotifs, appOccur, appFct) {

  setSupMot <- unique(supMotifs)
  setAppMot <- unique(appMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2])
  colnames(mfct) <- colnames(appOccur)
  rownames(mfct) <- setMot

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- amean(appFct[indApp[indElt]])
    }
  }


  supFct  <- numeric(length(supMotifs))
  msk     <- logical(length(supFct))
  sizeSup <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {
    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)
    indAppMot <- which(appMotifs == motif)

    if ( (length(indSupMot) > 0) & (length(indAppMot) > 0) ) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)
        indAppElt <- which(appOccur[indAppMot, element] == 1)

        if ( (length(indSupElt) > 0) & (length(indAppElt) > 0) ) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] + mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] / sizeSup[indSupMot]
      msk[indSupMot]    <- TRUE
    }
  }

  supFct[msk == FALSE] <- NA

  return(supFct)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# version 2  without msk

predict.amean.byelt.sup <- function(supMotifs, supOccur,
                                    appMotifs, appOccur, appFct) {

  setSupMot <- unique(supMotifs)
  setAppMot <- unique(appMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2])
  colnames(mfct) <- colnames(appOccur)
  rownames(mfct) <- setMot

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- amean(appFct[indApp[indElt]])
    }
  }


  supFct  <- numeric(length(supMotifs))
  sizeSup <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {
    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] + mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] / sizeSup[indSupMot]
    }
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#version 1

predict.gmean.byelt.sup.old <- function(supMotifs, supOccur,
                                        appMotifs, appOccur, appFct) {

  setSupMot <- unique(supMotifs)
  setAppMot <- unique(appMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2])
  colnames(mfct) <- colnames(appOccur)
  rownames(mfct) <- setMot

  for (mot in seq_along(setAppMot)) {
    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- gmean(appFct[indApp[indElt]])
    }
  }


  supFct   <- numeric(length(supMotifs))
  supFct[] <- 1
  msk      <- logical(length(supFct))
  sizeSup  <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {

    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)
    indAppMot <- which(appMotifs == motif)

    if ( (length(indSupMot) > 0) & (length(indAppMot) > 0) ) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)
        indAppElt <- which(appOccur[indAppMot, element] == 1)

        if ( (length(indSupElt) > 0) & (length(indAppElt) > 0) ) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] * mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] ^ (1/sizeSup[indSupMot])
      msk[indSupMot]    <- TRUE
    }
  }

  supFct[msk == FALSE] <- NA

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#version 2 without msk

predict.gmean.byelt.sup <- function(supMotifs, supOccur,
                                    appMotifs, appOccur, appFct) {

  setSupMot <- unique(supMotifs)
  setAppMot <- unique(appMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2])
  colnames(mfct) <- colnames(appOccur)
  rownames(mfct) <- setMot

  for (mot in seq_along(setAppMot)) {
    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- gmean(appFct[indApp[indElt]])
    }
  }


  supFct   <- numeric(length(supMotifs))
  supFct[] <- 1
  sizeSup  <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {
    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] * mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] ^ (1/sizeSup[indSupMot])
    }
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed by excluding (LOO) the assembly to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.sup <- function(supMotifs, supOccur,
                        appMotifs, appOccur, appFct,
                        opt.mean = "amean",
                        opt.mod  = "bymot"  ) {

  option <- paste(opt.mean, opt.mod, sep = ".")

  return(
    switch(option,
           amean.bymot =
             predict.amean.bymot.sup(supMotifs, appMotifs, appFct),
           gmean.bymot =
             predict.gmean.bymot.sup(supMotifs, appMotifs, appFct),

           amean.byelt =
             predict.amean.byelt.sup(supMotifs, supOccur,
                                     appMotifs, appOccur, appFct),
           gmean.byelt =
             predict.gmean.byelt.sup(supMotifs, supOccur,
                                     appMotifs, appOccur, appFct)
    )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute assembly functioning and associated statistiques
#          using the Clustering model
#          for Supplementary Assemblies by knowing their elemental Composition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict.function.sup <- function(tree, res.prd,
                                 supOccur, appOccur, appFct,

                                 # options for computing
                                 opt.mean = "amean",
                                 opt.mod  = "byelt"    ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #  Cas 1 où on suppose le nombre d'éléments identiques
  #      dans apprentissage et supplémentaires'
  #  Prévoir le cas 2 où le nombre d'éléments est différent
  #      dans apprentissage et supplémentaire'

  nbElt      <- dim(appOccur)[2]    # ????

  nbAppAss   <- dim(appOccur)[1]
  nbSupAss   <- dim(supOccur)[1]

  mAssMotifs <- mAffect.motifs(tree, rbind(appOccur, supOccur))
  mAppMotifs <- mAssMotifs[ , 1:nbAppAss]
  mSupMotifs <- mAssMotifs[ , (nbAppAss + 1):(nbAppAss + nbSupAss)]

  # Compute the raw predictions
  mSup <- mSd <- tNbcl <- matrix(NA, nrow = nbElt, ncol = nbSupAss)
  rownames(mSup) <- rownames(mSd) <- rownames(tNbcl) <- seq_len(nbElt)
  colnames(mSup) <- colnames(mSd) <- colnames(tNbcl) <- rownames(supOccur)

  for (nbcl in seq_len(nbElt)) {

    mSup[nbcl, ] <- predict.sup(mSupMotifs[nbcl, ], supOccur,
                                mAppMotifs[nbcl, ], appOccur, appFct,
                                opt.mean, opt.mod)

    mSd[nbcl, ] <- res.prd$uRmse[[nbcl]][mSupMotifs[nbcl, ]]

    index <- which(mSup[nbcl, ] > appFct)
    mSd[nbcl, index] <- -mSd[nbcl, index]
  }


  # Compute the associated statistiques
  tNbcl[1, ] <- 1
  for (nbcl in 2:nbElt) {
    tNbcl[nbcl, ] <- nbcl
    tNbcl[nbcl, is.na(mSup[nbcl, ])] <- tNbcl[(nbcl - 1), is.na(mSup[nbcl, ])]
  }

  tSup <- mSup
  for (nbcl in 2:nbElt)
    tSup[nbcl, is.na(mSup[nbcl, ])] <- tSup[(nbcl - 1), is.na(mSup[nbcl, ])]


  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mSupMotifs[nbcl, ]))

  mStats <- compute.mStats(mSup, mSup, appFct, nbK)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  #                  , "tError" A vérifier

  res        <- list(rownames(supOccur), appFct, opt.mean, opt.mod,
                     mSupMotifs, tSup, mSd, tNbcl, mStats)
  names(res) <- c("names", "fct", "opt.mean", "opt.mod",
                  "mMotifs", "tSup", "tSd", "tNbcl", "tStats")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot all the necessary predicting figures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction.sup <- function(res, nbOpt,

                                titre    = "",

                                sup.prd  = FALSE,
                                sup.pub  = FALSE,
                                sup.opt  = TRUE,
                                ass.loc  = FALSE,
                                opt.all  = FALSE  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.all == TRUE)
      sup.prd <- sup.pub <- sup.opt <- ass.loc <- TRUE

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Plot analysis cluster by cluster
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt  <- dim(res$tStats)[1]
  nbAss  <- length(res$fct)

  indFct <- which(!is.na(res$fct))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Prediction of Supplementary data with",
                    nbcl, "clusters", sep = " ")
    plot.prediction.LOO(res$tSup[nbcl, indFct],
                        res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                        res$fct[indFct],
                        rep(1, nbAss)[indFct], nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #    using different colors for predictions with different cluster numbers
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.prd == TRUE)
    for (nbcl in seq_len(nbOpt)) {

      titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.LOO(res$tSup[nbcl, indFct],
                          res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                          res$fct[indFct],
                          res$tNbcl[nbcl, indFct], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      index <- which(res$tNbcl[nbcl, indFct] == (nbcl - 1))
      if (length(index) != 0) {
        index <- indFct[index]
        add.ass.names(res$tSup[nbcl, index],
                      res$tSup[nbcl, index] + res$tSd[nbcl, index],
                      res$fct[index],
                      res$tNbcl[nbcl, index], res$names[index])
      }
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #          with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.pub == TRUE)
    for (nbcl in seq_len(nbOpt)) {

      titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                      "clusters", sep = " ")
      plot.prediction.LOO(res$tSup[nbcl, indFct],
                          res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                          res$fct[indFct],
                          rep(1, nbAss)[indFct], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by experiment with the names of Assemblages
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (ass.loc == TRUE) {

    nbcl   <- nbOpt
    titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                    "clusters", sep = " ")

    plot.prediction.LOO(res$tSup[nbcl, indFct],
                        res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                        res$fct[indFct],
                        res$mMotifs[nbcl, indFct], nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)

    add.ass.names(res$tSup[nbcl, indFct],
                  res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                  res$fct[indFct],
                  res$mMotifs[nbcl, indFct], res$names[indFct])
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   END of FILE myPREDICTING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
