#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#           Set of Functions et Procedures
#                  for manipulating Matrix of Abundances
#
#                     (Benoit JAILLARD, February 2018)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  Compute the Shannon index
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Shannon.index <- function(x) {
  n   <- sum(x)
  s   <- length(x)
  res <- sum(x/n * log2(n/x))

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  Compute the Piélou index
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Pielou.index <- function(x) {
  res <- Shannon.index(x) / log2(length(x))

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Shorten the names of elements
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

shorten.names <- function(full.names,
                          lmax = 4, separator = "") {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Evaluate the shortest possible length of names
  #            according the number of doublons, that should be lower than 9
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  cut <- 1
  while (max(table(substr(full.names, 1, cut))) > 9) cut <- cut + 1
  cut <- cut - 1
  if (cut < lmax) cut <- lmax

  short.names <- substr(full.names, 1, cut)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Build short names as : radical + sep + 1-9, for distinguishing the doublons
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  doublon <- table(short.names)
  for (dbl in 1:length(doublon)) {
    index <- which(short.names == names(doublon[dbl]))
    if (length(index) > 1)
      short.names[index] <- paste(names(doublon[dbl]),
                                  c(1:length(index)), sep = separator)
  }

  return(short.names)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Convert a numeric abundance matrix into a binary occurrence matrix
#
#  Each species is treated separately.
#  The distribution of each species on the whole experiment is computed,
#    then partitioned in nbLev intervals centred on the opt option.
#    Opt = "amean", "gmean" or "median"
#    lev = the level of threshold in the interval scale (lev <= nbLev).
#
#  Inputs : mat   : the matrix to binarise
#           lev   : the level to cut the matrix
#           nbLev : the number of intervals that matrix.
#                   it is 1 for cutting the at the median
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

segment.vector <- function(v, nbSeg,
                           option = "median"       # "amean", "gmean"
                           ) {

  sv <- sort(v, index.return = TRUE)
  sx <- sv$x - min(sv$x)

  # option = "median"
  sw <- sx
  sw[sw != 0] <- 1
  sv.median  <- cumsum(sw)

    # option = "amean"
  sv.amean   <- cumsum(sx)

  # option = "gmean"
  sv.gmean   <- c(rep(0, sum(sx == 0)), cumsum(log(sx[sx != 0])))

  order <- switch(option,
                  median = sv.median,
                  amean  = sv.amean,
                  gmean  = sv.gmean    )

  thres                <- numeric(nbSeg + 1)
  thres[2:(nbSeg + 1)] <- seq_len(nbSeg) * round(max(order)/nbSeg)
  thres[nbSeg + 1]     <- max(order)


  index  <- numeric(length(v))
  for (i in 1:nbSeg)
    index[sv$ix[which((order > thres[i]) & (order <= thres[i + 1]))]] <- i

  return(index)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

binarise.matrix <- function(mAbund, mOccur,
                            offset = 0, pOccur = 0.30,
                            option = "median"    # "amean", "gmean"
                            ) {

  # clean the matrix for non-occurring or too scarse element

  nbOcc  <- apply(mOccur, MARGIN = 2, FUN = sum)
  mOccur <- mOccur[ , nbOcc > offset]

  nbAss  <- dim(mOccur)[1]
  nbElt  <- dim(mOccur)[2]
  nbOcc  <- apply(mOccur, MARGIN = 2, FUN = sum)


  # determine the number of classes
  #      for the probability to be observed are close to pOccur

  nbLev  <- round(nbOcc/nbAss / pOccur)
  nbLev[nbLev == 0] <- 1


  # split the observations in regard to the number of classes
  mat <- matrix(0, nrow = nbAss, ncol = max(nbLev))
  colnames(mat) <- c(1:max(nbLev))
  rownames(mat) <- rownames(mAbund)

  resOccur <- resAbund <- NULL
  for (elt in 1:nbElt) {

    mat[,] <- 0
    nbcol  <- 1:nbLev[elt]

    index  <- segment.vector(mAbund[ , elt], nbLev[elt],
                             option = "median" )

    for (i in nbcol) mat[(index == i), i] <- mAbund[(index == i), elt]

    colnames(mat)[nbcol] <- paste(colnames(mAbund)[elt], nbcol, sep = ".")
    resAbund <- cbind(resAbund, mat[ , nbcol, drop = FALSE])

    mat[mat != 0] <- 1
    colnames(mat)[nbcol] <- paste(colnames(mOccur)[elt], nbcol, sep = ".")
    resOccur <- cbind(resOccur, mat[ , nbcol, drop = FALSE])
  }

  res        <- list(resAbund, resOccur, nbLev)
  names(res) <- c("mAbund", "mOccur", "nbLev")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  Plot a matrix-image of co-occurrence
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

co.occur <- function(mOccur) {

  nbElt <- dim(mOccur)[2]
  res   <- matrix(0, nrow = nbElt, ncol = nbElt)
  rownames(res) <- colnames(res) <- colnames(mOccur)

  for (elt in seq_len(nbElt)) {
    index      <- which(mOccur[ ,elt] == TRUE)
    res[elt, ] <- apply(mOccur[index, ], MARGIN = 2, FUN = sum) / length(index)
  }

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.co.occurrence <- function(mOccur,
                               order = c(1:dim(mOccur)[2]),
                               wdw   = c(1:dim(mOccur)[2])   ) {

  co.mOccur <- co.occur(mOccur)[order[wdw], order[wdw]]

  image(z = t(co.mOccur[rev(wdw), ]),  col = grey((9:0)/9),
        axes = FALSE, las = 1)

  coord <- seq(0, 1, by = 1/(dim(co.mOccur)[1] - 1) )

  axis(side = 3, at = coord,
       labels = colnames(mOccur)[order[wdw]],
       las = 2, tick = FALSE)
  axis(side = 2, at = coord,
       labels = colnames(mOccur)[rev(order[wdw])],
       las = 2, tick = FALSE)

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Generic function to plot a value vector versus sorted elements
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.fofelements <- function(v,
                             ylim   = range(v),
                             ytitre = "title",

                             index = seq_along(v),
                             wdw   = seq_along(v),

                             cols  = rep(couleurs[1], length(v)),
                             figs  = rep(figures[1], length(v))  ) {

  plot(y = v[index][wdw], ylim = ylim, ylab = ytitre,
       x = wdw, xlab = "elements",
       col = cols[index][wdw], pch = figs[index][wdw],
       las = 1, cex = 2)

  axis(at = seq_along(v)[wdw], side = 3,
       labels = names(v)[index][wdw],
       tick = FALSE, las = 2, col = "black", pos = max(v))

  abline(h = mean(v), col = "red")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# plot the occurrence of split or not elements

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.occurrence <- function(mOccur,
                            mOccur.split = mOccur,
                            nbLev = rep(1, dim(mOccur)[2]),

                            order = c(1:dim(mOccur)[2]),
                            windw = c(1:dim(mOccur)[2]),

                            cols  = rep(couleurs[1], length(order)),
                            figs  = rep(figures[1],  length(order))  ) {

  if (windw[length(windw)] > dim(mOccur)[2]) windw <- windw[1]:dim(mOccur)[2]

  # plot the whole occurrence of elements

  vfreq <- apply(mOccur, MARGIN = 2, FUN = sum) / dim(mOccur)[1]
  plot.fofelements(v = vfreq,
                   ytitre = "frequency of occurrence of split elements",
                   ylim   = c(0, 1),
                   index  = order, wdw = windw,
                   cols   = cols,  figs = figs)


  # plot the partial occurrence of split elements

  if ( (dim(mOccur.split)[2] > dim(mOccur)[2]) &
       (length(nbLev) == dim(mOccur)[2]) &
       (max(nbLev) > 1) ) {

    cumfreq     <- numeric(length(vfreq))
    vfreq.split <- apply(mOccur.split, MARGIN = 2, FUN = sum) /
                   dim(mOccur.split)[1]

    lev <- 1
    for (lev in 1:(max(nbLev) - 1)) {

      index        <- cumsum(nbLev) - nbLev + lev
      who          <- which(nbLev > lev)
      cumfreq[who] <- cumfreq[who] + vfreq.split[index[who]]
      points(x = who[windw], y = cumfreq[who][windw],
             col = couleurs[1 + lev], pch = figures[1 + lev], cex = 2)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# plot the abundance distribution
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.abundance <- function(mAbund, mOccur,
                           order = c(1:dim(mOccur)[2]),
                           windw = c(1:dim(mOccur)[2]),
                           scale.log = FALSE,

                           cols  = rep(couleurs[1], length(order)),
                           figs  = rep(figures[1],  length(order))  ) {

  if (scale.log == TRUE) {
    titre <- "log10(abundance)"
    fct   <- function(x) { log10(median(x[(x != 0)])) }
#    fct   <- function(x) { log10(mean(x[(x != 0)])) }
  } else {
    titre <- "abundance"
    fct   <- function(x) { median(x[(x != 0)]) }
#    fct   <- function(x) { mean(x[(x != 0)]) }
  }

  # compute the mean abundances for each element
  tmp <- apply(mAbund, MARGIN = 2, FUN = fct)
  names(tmp) <- colnames(mOccur)

  plot.fofelements(v = tmp, ytitre = titre,
                   index = order, wdw = windw,
                   cols  = cols, figs = figs)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Domaine de mOccur.split
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.occurring <- function(nbOpt, xwdw = c(1:50),
                           tree.cal, res.prd,
                           mAbund, mOccur, nbLev,
                           mAbund.split, mOccur.split) {

  affectElt <- cut.tree(tree.cal, nbOpt)

  col.tree  <- couleurs[shift.affectElt(affectElt)]
  fig.tree  <- figures[ shift.affectElt(affectElt)]

  # sort the elements by decreasing mean abundances
  #    tmp       <- apply(mAbund, MARGIN = 2, FUN = sum)
  tmp       <- apply(mAbund, MARGIN = 2, FUN = median)
  ind.abund <- sort(tmp, index.return = TRUE, decreasing = TRUE)$ix

  ind.tree  <- sort.tree(tree.cal, index.return = TRUE)$ix

  tmp       <- sort.split.names(affectElt, ind.tree,
                                colnames(mOccur), nbLev)
  ind.elt   <- tmp$index.elt
  ind.split <- tmp$index.split


  # plot the occurrence frequency sorted by decreasing value
  #       cluster by cluster
  plot.occurrence(mOccur, order = ind.abund)  #

  plot.occurrence(mOccur, order = ind.abund,
                  mOccur.split, nbLev)         #

  # plot the occurrence frequency sorted according to tree clustering
  #       cluster by cluster
  plot.occurrence(mOccur.split,
                  order = ind.tree, cols = col.tree, figs = fig.tree)

  plot.occurrence(mOccur.split,
                  order = ind.tree, cols = col.tree, figs = fig.tree,
                  windw = xwdw     )

  #       element by element
  plot.occurrence(mOccur.split,
                  order = ind.split, cols = col.tree, figs = fig.tree)

  plot.occurrence(mOccur.split,
                  order = ind.split, cols = col.tree, figs = fig.tree,
                  windw = xwdw)


  # plot the element abundaces sorted by median values
  plot.abundance(mAbund, mOccur,
                 order = ind.elt, scale.log = FALSE)

  plot.abundance(mAbund, mOccur,
                 order = ind.elt, scale.log = TRUE)

  # plot the element abundances sorted according to tree clustering
  #       cluster by cluster
  plot.abundance(mAbund.split, mOccur.split,
                 order = ind.tree,
                 cols  = col.tree, figs = fig.tree,
                 scale.log = TRUE)

  plot.abundance(mAbund.split, mOccur.split,
                 order = ind.tree,
                 cols  = col.tree, figs = fig.tree,
                 scale.log = TRUE,
                 windw = xwdw     )

  #       element by element
  plot.abundance(mAbund.split, mOccur.split,
                 order = ind.split,
                 cols  = col.tree, figs = fig.tree,
                 scale.log = TRUE)

  plot.abundance(mAbund.split, mOccur.split,
                 order = ind.split,
                 cols  = col.tree, figs = fig.tree,
                 scale.log = TRUE,
                 windw = xwdw     )
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
unique.labelling <- function(v, data) {

  res <- character(dim(data)[1])
  for (i in 1:dim(data)[1]) {
    tmp <- NULL
    for (j in 1:length(v)) tmp <- paste(tmp, data[i, v[j]], sep = ".")
    res[i] <- substr(tmp, 2, nchar(tmp))
  }

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Re-concatène les abondances qui appartiennent au même cluster

which.str <- function(str, chr) {
  i   <- 1
  while ( (substr(str, i, i) != chr) & (i < nchar(str)) ) i <- i + 1
  if (i == nchar(str)) i <- 0
  return(i)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# seconde version

auPoint <- FALSE
if (auPoint == TRUE) {

  affectElt <- cut.tree(tree.cal, nbOpt)

  nbElt     <- dim(mOccur.split)[2]
  nbAss     <- dim(mOccur.split)[1]
  posPnt    <- integer(nbElt)
  for (elt in 1:nbElt) posPnt[elt] <- which.str(names(affectElt)[elt], ".")

  setRad  <- substr(names(affectElt), 1, posPnt - 1)
  setRad[posPnt == 0] <- names(affectElt)[posPnt == 0]

  setFix  <- substr(names(affectElt), posPnt + 1, nchar(names(affectElt)))
  setFix[posPnt == 0] <- 0

  mat  <- mOccur.split

  for (nbcl in 1:length(unique(affectElt))) {

    indClu  <- which(affectElt == nbcl)
    doublon <- table(setRad[indClu])

    for (elt in 1:length(doublon)) {

      indRad <- which(setRad == names(doublon)[elt])
      indFix <- which(indRad %in% indClu)
      index  <- indRad[indFix]

      nblev  <- nbLev[which(names(nbLev) == names(doublon)[elt])]

      # compact the redundant vectors of Abundance
      mat[ , index[1]] <- apply(mOccur.split[ , index, drop = FALSE],
                                MARGIN = 1, FUN = sum)

      # generate a compact name
      if (setFix[index][1] == 0) {
        colnames(mat)[index[1]] <- setRad[index][1]
      } else {
        tmp <- ""
        for (i in 1:length(indFix)) tmp <- paste0(tmp, setFix[index][i])
        colnames(mat)[index[1]] <- paste(setRad[index][1], tmp, sep = ".")
        colnames(mat)[index[-1]] <- c("X")
      }
    }
  }

  indexOK  <- which(colnames(mat) != "X")
  resOccur <- mat[ , indexOK]

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   END of MAIN
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

