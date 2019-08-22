#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                 myCOMBINAT.R
#
#  set of functions for combinatorial analysis of community diversity effects
#
#                         Benoit JAILLARD, summer 2016
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


source("myStats.R")
source("myCombinat.R")



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                               mySEPARATING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot a vector of points with or without std bars
#
#  Inputs : x, y : vectors of coordinates
#           figure, couleur, nom : informations to plot the vectors
#           opt.std : with (TRUE) or without( FALSE) error-bars
#           scale.log : plot in log-scale (TRUE) or not (FALSE).
#           If TRUE, the vectors x and y must be input as log2(x) and log2(y)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

points.sd <- function(x, y,
                      figure, couleur, nom,
                      opt.std = TRUE,
                      scale.log = FALSE) {

  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # plot at first the error-bars
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if ((opt.std == TRUE) & (length(x) > 1)) {
    dx  <- sd(x)
    dy  <- sd(y)

    if (scale.log == FALSE) {
      if (dx > EPSILON) arrows(x0 = mx - dx, y0 = my,
                          x1 = mx + dx, y1 = my,
                          col = couleur, length = 0.1, angle = 90, code = 3,
                          lwd = 1, lty = "solid")
      if (dy > EPSILON) arrows(x0 = mx, y0 = my - dy,
                          x1 = mx, y1 = my + dy,
                          col = couleur, length = 0.1, angle = 90, code = 3,
                          lwd = 1, lty = "solid")
    } else {

      if ((dx > EPSILON) & ((mx - dx) > EPSILON))
        arrows(x0 = log2(mx - dx), y0 = log2(my),
               x1 = log2(mx + dx), y1 = log2(my),
               col = couleur, length = 0.1, angle = 90, code = 3,
               lwd = 1, lty = "solid")
      if ((dy > EPSILON) & ((my - dy) > EPSILON))
        arrows(x0 = log2(mx), y0 = log2(my - dy),
               x1 = log2(mx), y1 = log2(my + dy),
               col = couleur, length = 0.1, angle = 90, code = 3,
               lwd = 1, lty = "solid")
    }
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # plot second the mean points with a white background
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (scale.log == FALSE) {
    points(x = mx, y = my,
           type = "p", pch = figure, col = couleur, bg = "white", cex = 4)
    text(x = mx, y = my,
         labels = nom, pos = 4, col = couleur, bg = "white")

  } else {

    points(x = log2(mx), y = log2(my),
           type = "p", pch = figure, col = couleur, bg = "white", cex = 4)
    text(x = log2(mx), y = log2(my),
         labels = nom, pos = 4, col = couleur, bg = "white")
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot alpha x beta
#
#  Inputs : x, y : alpha, beta vectors of coordinates
#           figure, couleur, nom : informations to plot the vectors
#           opt.std : with (TRUE) or without( FALSE) error-bars
#           scale.log : plot in log-scale (TRUE) or not (FALSE).
#           If TRUE, the vectors x and y must be input as x and y
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.alphaVSbeta <- function(x, y, xlab, ylab,
                             xlim = range(x), ylim = range(y),
                             titre = "",
                             scale.log = FALSE) {

  R2 <- cor2(x, y)

  if (scale.log == FALSE) {

    plot(x = x, xlim = xlim, xlab = xlab,
         y = y, ylim = ylim, ylab = ylab,
         main = titre,
         las = 1, type = "n", tck = 0.02)

    abline(h = 1, col = "black")
    abline(v = mean(x), h = mean(y), col = "blue")

    points(x = x, y = y, cex = 2, pch = 21, bg = "white")

    text(x = min(x) + 3/4*(max(x) - min(x)), y = min(y),
         labels = paste("R2 = ", signif(R2, digits = 3), sep = ""),
         pos = 3, col = "blue")

  } else {

    log2x <- log2(x)
    log2y <- log2(y)

    plot(x = log2x, xlim = log2(xlim),
         xlab = paste("log2(", xlab, ")", sep = ""),
         y = log2y, ylim = log2(ylim),
         ylab = paste("log2(", ylab, ")", sep = ""),
         main = paste("log2(", titre, ")", sep = ""),
         las = 1, type = "n", tck = 0.02)

    abline(h = 0, v = 0, col = "black")
    abline(v = log2(mean(x)), h = log2(mean(y)), col = "blue")

    points(x = log2x, y = log2y, cex = 2, pch = 21, bg = "white")

    text(x = min(log2x) + 3/4*(max(log2x) - min(log2x)), y = min(log2y),
         labels = paste("R2 = ", signif(R2, digits = 3), sep = ""),
         pos = 3, col = "blue")
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Remove species assemblies that were in several exemplar
#     and return the right mat (with only one assemblage)
#     and fct (as mean of fucntions of double assemblages)
#
#   Inputs  : mOccur and Fobs : matrix of occurrence and function
##  Outputs : mOccur and Fobs, without doublons
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm.dual.assemblies <- function(mat, fct,
                               xpr      = rep(1, length(fct)),
                               opt.mean = "amean") {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setXpr <- unique(xpr)
  setAss <- which(xpr == setXpr[1])
  pas    <- length(setAss)
  mut    <- mat[setAss, ]

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  AssMotifs  <- affect.motifs(c(1:dim(mut)[2]), mut)
  index      <- table(AssMotifs)
  bool       <- !logical(length(AssMotifs))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setDoublon <- which(index > 1)
  for (i in seq_along(setDoublon)) {
    doublon <- which(AssMotifs == setDoublon[i])

    for (ipr in seq_along(setXpr)) {
      off <- pas * (ipr - 1)
      fct[doublon[1] + off] <- mean.fct(fct[doublon + off], opt.mean)
      bool[doublon[2:length(doublon)] + off] <- FALSE
    }
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  res        <- list(mat[bool, ], fct[bool])
  names(res) <- c("mat", "fct")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Remove one or more species from data
#     and return the index.row and index.col to select
#
#     Inputs : mOccur and Fobs : matrix of occurrence and function
#             elements : a vector of names of elements to remove.
#                        Names of elements must be in colnames of mOccur.
#     Outputs : mOccur and Fobs, without the elements to remove
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm.elements <- function(mOccur, fct, elements) {

  index.row <- c(1:dim(mOccur)[1])
  index.col <- c(1:dim(mOccur)[2])

  for (elt in seq_along(elements)) {
    index.row <- setdiff(index.row, which(mOccur[ ,elements[elt]] == 1))
    index.col <- setdiff(index.col, which(elements[elt] == colnames(mOccur)))
  }

  mOccur <- mOccur[index.row, index.col]
  fct    <- fct[index.row]

  res        <- list(mOccur, fct)
  names(res) <- c("mOccur", "fct")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Decompose the observed performance into a.inter and a.comp ratios
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

multiplicative.decomposition <- function(mOccur, Fobs,
                                         rm.mono = TRUE) {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check if all monocultures are observed
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  size      <- apply(mOccur, MARGIN = 1, FUN = sum)
  mono      <- which(size == 1)
  if (length(mono) == 0) stop("There is no monoculture in mOccur")

  elements  <- which(apply(mOccur[mono,], MARGIN = 2, FUN = sum) != 0)

  if (length(elements) < dim(mOccur)[2]) {
    missing <- setdiff(colnames(mOccur), colnames(mOccur)[elements])
    res     <- rm.elements(mOccur, Fobs, missing)

    mOccur  <- res$mOccur
    Fobs    <- res$fct
    size    <- apply(mOccur, MARGIN = 1, FUN = sum)
    mono    <- which(size == 1)
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Look for and compute the monoculture performance
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  Fmono <-
    apply(mOccur[mono, ] * Fobs[mono], MARGIN = 2, FUN = sum) /
    apply(mOccur[mono, ], MARGIN = 2, FUN = sum)

  Fref    <- (mOccur %*% Fmono) / size
  Fscale  <- mean(Fmono)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Compute the interaction and composition effects
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  alpha     <- as.vector(Fobs/Fref)
  beta      <- as.vector(Fref/Fscale)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Remove the monocultures
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (rm.mono == TRUE) {
    index  <- which(size != 1)
    size   <- size[index]
    mOccur <- mOccur[index, ]
    Fobs   <- Fobs[index]
    alpha  <- alpha[index]
    beta   <- beta[index]
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Identify the transgressive performances
  # 1: transgressive under-yielding;  2: under-yielding
  # 3: over-yielding ;                4: transgressive over-yielding
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  trans <- numeric(dim(mOccur)[1])
  tmp   <- t( t(mOccur) * Fmono )
  for (i in 1:dim(mOccur)[1]) {
    ttp         <- tmp[i, (tmp[i,] != 0)]
    if (Fobs[i] <  mean(ttp)) trans[i] <- 2
    if (Fobs[i] <  min(ttp))  trans[i] <- 1
    if (Fobs[i] >= mean(ttp)) trans[i] <- 3
    if (Fobs[i] >  max(ttp))  trans[i] <- 4
  }


  res        <- list(size, trans, mOccur,
              Fobs, alpha, beta, Fmono, Fscale)
  names(res) <- c("size", "trans", "mOccur",
                  "Fobs", "alpha", "beta", "Fmono", "Fscale")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the histogram of fct-values
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ghisto <- function(x, namex, scale.log = FALSE) {

  if (scale.log == FALSE) {

    hist(x, xlab = namex, ylab = "Density", prob = TRUE, main = namex,
         las = 1, breaks = 20)
    abline(v = 1, col = "black")
    abline(v = gmean(x), col = "blue")

  } else {

    titre <- paste("log2(", namex, ")", sep = "")
    hist(log2(x), xlab = titre, ylab = "Density", prob = TRUE, main = titre,
         las = 1, breaks = 20)
    abline(v = 0, col = "black")
    abline(v = log2(gmean(x)), col = "blue")
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the dynamic trajectories of elements on alpha x beta biplot
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

scex <- function(n) {

  p1 <- c(1, 1)
  p2 <- c(n, 5)

  a   <- (p2[2] - p1[2]) / (p2[1] - p1[1])
  b   <-  p2[2] - a * p2[1]
  pas <- (p2[1] - p1[1]) / (n - 1)
  x   <- p1[1] + ((1:n) - 1) * pas
  res <- a * x + b

  return(res)
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
plot.dynamics <- function(alpha, beta, mOccur, xpr,
                          figs = figures[1], cols = couleurs[1] ) {

  nbElt  <- dim(mOccur)[2]
  setXpr <- unique(xpr)

  mAlphaElt <- mBetaElt <- msdAlphaElt <- msdBetaElt <-
    matrix(0, nrow = length(setXpr), ncol = nbElt)
  colnames(mAlphaElt) <- colnames(mBetaElt) <-
    colnames(msdAlphaElt) <- colnames(msdBetaElt) <- colnames(mOccur)

  for (tim in seq_along(setXpr))
    for (elt in 1:nbElt) {

      index <- which( (xpr == setXpr[tim]) & (mOccur[ , elt] == 1) )

      mAlphaElt[tim, elt]   <- gmean(alpha[index])
      msdAlphaElt[tim, elt] <- gsd(alpha[index])

      mBetaElt[tim, elt]    <- amean(beta[index])
      msdBetaElt[tim, elt]  <- asd(beta[index])
  }


  xlim <- c(0.95, 1.05) * range(mBetaElt)
  ylim <- c(0.95, 1.05) * range(mAlphaElt)

  if (length(figs) == 1) figs <- rep(figs, nbElt)
  if (length(cols) == 1) cols <- rep(cols, nbElt)

  # A figure for all Elements at different Times
  plot(x = mBetaElt[1, ],  xlab = "beta",  xlim = xlim,
       y = mAlphaElt[1, ], ylab = "alpha", ylim = ylim,
       main = "All elements",
       type = "n", tck = 0.02, las = 1)
  abline(h = 1, v = 1, lty = "solid", col = "blue")

  for (elt in 1:nbElt) {
    points(x = mBetaElt[ , elt],
           y = mAlphaElt[ , elt],
           type = "b", pch = figures[figs[elt]],
           col = couleurs[cols[elt]], bg = "white",
           cex = scex(length(setXpr)) )

    posT <- 3                                               # over the point
    if (length(setXpr) > 1)
       if (mAlphaElt[2, elt] > mAlphaElt[1, elt]) posT <- 1  # under the point

    text(x = mBetaElt[1, elt], y = mAlphaElt[1, elt],
         labels = colnames(mOccur)[elt],
         col = couleurs[cols[elt]], pos = posT)
  }


  # A figure by Element at different Times
  for (elt in 1:nbElt) {

    titre <- colnames(mOccur)[elt]
    plot(x = mBetaElt[1, ],  xlab = "beta",  xlim = xlim,
         y = mAlphaElt[1, ], ylab = "alpha", ylim = ylim,
         main = titre,
         type = "n", tck = 0.02, las = 1)
    abline(h = 1, v = 1, lty = "solid", col = "blue")

    points(x = mBetaElt[ , elt],
           y = mAlphaElt[ , elt],
           type = "b", pch = figures[figs[elt]],
           col = couleurs[cols[elt]], bg = "white",
           cex = scex(length(setXpr)) )
  }


  # A figure by Assembly at different Times
  nbAss <- length(which(xpr == setXpr[1]))
  index <- (seq_along(setXpr) - 1) * nbAss

  xlim <- c(0.9, 1.1) * range(beta)
  ylim <- c(0.9, 1.1) * range(alpha)

  figs <- check.symbol(figs, nbAss)
  cols <- check.symbol(cols, nbAss)

  size <- apply(mOccur, MARGIN = 1, FUN = sum)
  figs <- size
  cols <- size

  for (ass in 1:nbAss) {

    titre <- rownames(mOccur)[ass]
    plot(x = beta[ass + index],  xlab = "beta",  xlim = xlim,
         y = alpha[ass + index], ylab = "alpha", ylim = ylim,
         main = titre,
         type = "n", tck = 0.02, las = 1)
    abline(h = 1, v = 1, lty = "solid", col = "blue")

    points(x = beta[ass + index],
           y = alpha[ass + index],
           type = "b", pch = figures[figs[ass]],
           col = couleurs[cols[ass]], bg = "white",
           cex = scex(length(setXpr)) )
  }

}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                         END of FILE mySEPARATING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
