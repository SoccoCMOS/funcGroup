#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                 myPLOT.R
#
#                set of functions for plotting the results
#                of elements or motifs sorting and clustering
#
#                         Benoit JAILLARD, summer 2016
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


source("myStats.R")


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   List of functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


# plot.histogram    <- function(y, clusters, titre, pas = 0,
#                               opt.mean = "amean")
#
# sort.motifs       <- function(fct, AssMotifs, AssNoms)
#                               titre    = "",
#                               pvalue   = 0.05,
#                               opt.hor  = TRUE  )
#
# sort.tree                <- function(X, index.return = FALSE)
# plot.tree                <- function(tree, col = "black", titre = "")
#
# delstr.begin             <- function(v, n)
# delstr.end               <- function(v, n)
# concat.byline            <- function(v, nbchar = 70)
# plot.bypage              <- function(tab, nbline = 25)
# plot.motifs.content      <- function(fct, AssMotifs, AssNoms)
# plot.clusters.content    <- function(affectElt, mOccur)
#
# plot.tree.all            <- function(tree, Fobs, mOccur,
#                                      titre  = "",
#                                     tre.cal = FALSE,
#                                     tre.prd = FALSE,
#                                     tre.pub = FALSE,
#                                     tre.leg = FALSE,
#                                     mot.cal = FALSE, pvalue = 0.025,
#                                     mot.pub = FALSE,
#                                     mot.leg = FALSE,
#                                     opt.all = FALSE)





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot the clusters of assemblies as histograms
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.histogram <- function(y, clusters, titre,
                           pas = 0,
                           opt.mean = "amean") {

  ymin <- floor(min(y))
  ymax <- 1 + floor(max(y))

  if (pas == 0) pas <- (ymax - ymin)/20
  breaks <- seq(ymin, ymax, by = pas)

  hist(y, prob = TRUE, las = 1, ylab = "Density",
       main = titre, breaks = breaks)

  for (iter in seq_along(table(clusters)))
    if (length(y[(clusters == iter)]) > 1) {

      tmp <- density(y[(clusters == iter)])
      lines(x = tmp$x, y = tmp$y * length(y[(clusters == iter)]) / length(y),
            col = couleurs[iter])

      abline(v = mean.fct(y[(clusters == iter)], opt.mean),
             col = couleurs[iter])

      rug(y[(clusters == iter)], col = couleurs[iter])
    }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot the clusters of assemblies as a boxplot
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort.motifs <- function(fct, AssMotifs, AssNoms,
                        opt.aov = TRUE,  pvalue = 0.05,
                        opt.dec = TRUE     ) {

  # opt.dec = TRUE by default because test.posthoc() is "decreasing"

  sres  <- test.posthoc(x = fct, clusters = AssMotifs, pvalue = pvalue)

  indAss <- unique(AssMotifs)
  indNom <- unique(AssNoms)

  snoms  <- character(length(sres$motifs))
  nombre <- integer(length(sres$motifs))

  for (mot in seq_along(sres$motifs)) {
    motif       <- sres$motifs[mot]
    snoms[mot]  <- indNom[which(indAss == motif)]
    nombre[mot] <- length(which(AssMotifs == motif))
  }

  sres <- cbind(snoms, nombre, sres)
  colnames(sres)[1:5] <- c("nom", "nombre", "motif", "mean", "group")

  if (opt.dec == FALSE)
    for (j in seq_len(dim(sres)[2])) sres[ ,j] <- rev(sres[ ,j])

  rownames(sres) <- c(1:dim(sres)[1])

  return(sres)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.box.motifs <- function(fct, AssMotifs, sres,

                            ordre   = NULL,
                            ylim    = range(fct),
                            titre   = "",
                            opt.aov = TRUE, pvalue = 0.05,
                            opt.hor = TRUE) {

  # re-organize the motifs to appear as continuously increasing or decreasing
  ass.sort <- integer(length(fct))
  if (length(ordre) == 0) {

    for (mot in seq_along(sres$motif)) {
      motif <- sres$motif[mot]
      ass.sort[which(AssMotifs == motif)] <- mot
    }

  } else {

    index <- integer(length(ordre))
    for (mot in seq_along(ordre)) {
      motif      <- ordre[mot]
      index[mot] <- which(sres$motif == motif)
      ass.sort[which(AssMotifs == motif)] <- mot
    }
    sres <- sres[index, ]
  }


  # horizontal plotting of boxes is the most complete possibly
  if (opt.hor == TRUE) {

    boxplot(fct ~ ass.sort,
            names = sres$nom,
            ylab  = "Assembly motifs",
            xlab  = "Observed Function", ylim = ylim,
            horizontal = opt.hor, las = 1)

    points(y = seq_along(sres$motif), x = sres$mean,
           pch = 0, col = couleurs[seq_along(sres$motif)], bg = "white")

    lines(x = rep(amean(fct), 2),
          y = range(ass.sort), col = "red", lty = "solid")

    text(y = seq_along(sres$motif), x = min(fct), pos = 2,
         col    = couleurs[seq_along(sres$motif)],
         labels = as.character(sres$nombre) )

    text(y = seq_along(sres$motif), x = max(fct), pos = 2,
         col    = couleurs[seq_along(sres$motif)],
         labels = as.character(sres$group),
         font   = 3)

  } else {

    # vertical plotting of boxes

    boxplot(fct ~ ass.sort,
            names = sres$nom,
            xlab  = "Assembly motifs",
            ylab  = "Observed Function", ylim = ylim,
            horizontal = opt.hor, las = 2)

    points(x = seq_along(sres$motif), y = sres$mean,
           pch = 0, col = couleurs[seq_along(sres$motif)], bg = "white")

    lines(y = rep(amean(fct), 2),
          x = range(ass.sort), col = "red", lty = "solid")

    axis(side = 3, at = seq_along(sres$motif),
         labels = as.character(sres$group),
         col    = couleurs[seq_along(sres$motif)],
         tick   = FALSE, las = 2, font = 3 )

  }

  title(titre)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot the clusters of assemblies that contain a species as a boxplot
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

reverse.matrix <- function(mat) {

  for (j in seq_len(dim(mat)[2])) mat[ ,j] <- rev(mat[ ,j])
  return(mat)

}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

sort.elements <- function(fct, mOccur, cols,
                          pvalue  = pvalue,
                          opt.dec = TRUE  ) {

  # opt.dec = TRUE by default because test.posthoc() is "decreasing"

  nbElt  <- dim(mOccur)[2]
  if (length(cols) < nbElt) cols <- check.symbol(cols, nbElt)

  numElt <- integer(nbElt)
  fctElt <- motElt <- NULL

  for (elt in seq_len(nbElt)) {

    index       <- which(mOccur[ , elt] == 1)
    numElt[elt] <- length(index)
    fctElt      <- c(fctElt, fct[index])
    motElt      <- c(motElt, rep(elt, length(index)))
  }

  sres   <- test.posthoc(x = fctElt, clusters = motElt, pvalue = pvalue)
  nomElt <- colnames(mOccur)[sres$motifs]
  numElt <- numElt[sres$motifs]
  coll   <- cols[sres$motifs]

  sres   <- cbind(nomElt, coll, numElt, sres)
  colnames(sres)[1:6] <- c("nom", "color", "nombre", "index", "mean", "group")

  if (opt.dec == FALSE) sres <- reverse.matrix(sres)

  rownames(sres) <- c(1:(dim(sres)[1]))

  ind.inturn <- index.inturn(sres$index)
  sres       <- cbind(ind.inturn, sres)

  return(sres)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Get together the abundances and occurrences / seconde version

sort.split.names <- function(affectElt, ind.tree,
                             short.names, nbLev) {

  nbElt.split <- length(affectElt)

  posPnt <- integer(nbElt.split)
  for (elt in 1:nbElt.split)
    posPnt[elt] <- which.str(names(affectElt)[elt], ".")

  setRad  <- substr(names(affectElt), 1, posPnt - 1)
  setRad[posPnt == 0] <- names(affectElt)[posPnt == 0]

  sort.short.names <- unique(setRad[ind.tree])

  index <- integer(length(short.names))
  for (elt in seq_along(index))
    index[elt] <- which(short.names == sort.short.names[elt])

  # index tel que : s.short.names <- short.names[index]

  index.split        <- integer(nbElt.split)
  names(index.split) <- names(affectElt)

  count <- 0
  for (elt in seq_along(index)) {

    element          <- index[elt]
    ind.split        <- sum(nbLev[1:element]) - nbLev[element] + 1
    who              <- ind.split + 1:nbLev[element] - 1

    index.split[who] <- count + 1:nbLev[element]
    count            <- count + nbLev[element]
  }

  index.split  <- sort(index.split, index.return = TRUE)$ix

  res        <- list(index, index.split)
  names(res) <- c("index.elt", "index.split")

  return(res)
}


#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.box.elements <- function(fct, mOccur, sres,

                              ordre   = NULL,
                              ylim    = range(fct),
                              elt.wdw = seq_len(dim(mOccur)[2]),

                              titre   = "",
                              opt.aov = TRUE, pvalue = 0.05,
                              opt.hor = TRUE  ) {

  nbElt  <- dim(mOccur)[2]

  # re-organize the motifs to appear as continuously increasing or decreasing
  fctElt <- motElt <- NULL
  for (elt in seq_len(nbElt)) {
    element <- sres$index[elt]
    indElt  <- which(mOccur[ , element] == 1)
    fctElt  <- c(fctElt, fct[indElt])
    motElt  <- c(motElt, rep(elt, length(indElt)))
  }
  elt.sort <- motElt


  # sort the sres-table according input order
  fres <- sres
  if (length(ordre) != 0) {
    index <- integer(length(ordre))
    for (elt in seq_along(ordre)) {
      element    <- ordre[elt]
      index[elt] <- which(sres$index == element)
    }
    fres <- sres[index, ]
  }


  # horizontal plotting of boxes is the most complete possibly
  indxy <- compact.index(elt.wdw)

  if (opt.hor == TRUE) {

    boxplot(fctElt ~ elt.sort,
            subset = (elt.sort %in% elt.wdw),
            names  = fres$nom[elt.wdw],
            ylab   = "Element inside communities",
            xlab   = "Observed Function", ylim = ylim,
            horizontal = opt.hor, las = 1)

    points(y   = indxy,
           x   = fres$mean[elt.wdw],
           col = as.character(fres$color[elt.wdw]),
           pch = 0, bg = "white"  )

    lines(x = rep(amean(fct), 2),
          y = range(indxy),
          col = "red", lty = "solid")

    text(labels = as.character(fres$nombre[elt.wdw]),
         y      = indxy,
         x      = min(fct),
         col    = as.character(fres$color[elt.wdw]),
         pos    = 2  )

    text(labels = as.character(fres$group[elt.wdw]),
         y      = indxy,
         x      = max(fct), pos = 2,
         col    = as.character(fres$color[elt.wdw]),
         font   = 3)

  } else {

    # vertical plotting of boxes

    boxplot(fctElt ~ elt.sort,
            subset = (elt.sort %in% elt.wdw),
            names  = fres$nom[elt.wdw],
            xlab   = "Element inside communities",
            ylab   = "Observed Function", ylim = ylim,
            horizontal = opt.hor, las = 2)

    points(x   = indxy,
           y   = fres$mean[elt.wdw],
           col = as.character(fres$color[elt.wdw]),
           pch = 0, bg = "white"  )

    lines(y = rep(amean(fct), 2),
          x = range(indxy), col = "red", lty = "solid")

    axis(labels = as.character(fres$group[elt.wdw]),
         at     = indxy,
         col    = as.character(fres$color[elt.wdw]),
         side   = 3,
         tick   = FALSE, las = 2, font = 3 )

  }

  title(titre)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot the contents of clusters of elements or motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

delstr.begin <- function(v, n) {
  return(substr(v, n + 1, nchar(v)))
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

delstr.end <- function(v, n) {
  return(substr(v, 1, nchar(v) - n))
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# thres = nb of characters by line

concat.byline <- function(v, nbchar = 70) {

  res   <- NULL

  i <- 1
  while (i <= length(v)) {

    a <- ""
    while ((nchar(a) < nbchar) & (i <= length(v)) ) {
      a <- paste(a, v[i], sep = ", ")
      i <- i + 1
    }
    a <- paste(a, ",", sep = "")

    res <- c(res, delstr.begin(a, 2))
  }

  return(res)
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  nbline = number of lines by page

plot.bypage <- function(tab, nbline = 25) {

  page <- 0
  ltab <- length(tab)
  while (nbline * page < ltab) page <- page + 1

  for (p in 1:page) {

    tmp <- tab[intersect((p - 1) * nbline + 1:nbline, seq_along(tab))]

    plot(x = seq_len(nbline), y = seq_len(nbline),
         xlab = "",          ylab = "",
         type = "n", axes = FALSE)

    axis(side = 4, at = seq_along(tmp), rev(tmp),
         pos = 1, tick = FALSE, las = 2, col = "black",
         family = "mono")
  }
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.clusters.content <- function(affectElt, mOccur,
                                  opt.sort = TRUE,
                                  nbline   = 25) {

  noms   <- names.assembly(affectElt, mOccur)

  setClu <- table(noms$clusters)
  nbrow  <- length(setClu)

  # compute the number of necessary columns
  tmp   <- list()

  for (i in seq_len(nbrow))
    tmp[[i]] <- noms$elements[which(noms$clusters == names(setClu)[i])]

  if (opt.sort == TRUE)
    for (i in seq_len(nbrow)) tmp[[i]] <- sort(tmp[[i]])

  # adjust the length of rownames
  lmax <- max(nchar(names(setClu)))
  for (i in seq_len(nbrow)) {
    a <- names(setClu)[i]
    while (nchar(a) < lmax) a <- paste(" ", a, sep = "")
    names(tmp[[i]])[1] <- a
  }

  blanc <- NULL
  for (k in 1:(lmax + 5)) blanc <- paste(blanc, " ", sep = "")

  # concat the matrix by row
  tab <- c("Clusters content : (sorted by decreasing effect on Fobs)", "")
  for (i in seq_len(nbrow)) {

    tpp    <- concat.byline(tmp[[i]])
    tpp[1] <- paste(names(tmp[[i]])[1], " = { ", tpp[1], sep = "")
    if (length(tpp) > 1)
      for (j in 2:length(tpp)) tpp[j] <- paste(blanc, tpp[j], sep = "")
    tpp[length(tpp)] <- paste(delstr.end(tpp[length(tpp)], 1), " }", sep = "")

    tab  <- c(tab, tpp)
  }

  # plot the legend
  plot.bypage(tab, nbline)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

plot.motifs.content <- function(fct, AssMotifs, AssNoms,
                                opt.sort = TRUE  ) {

  tmp     <- sort.motifs(fct, AssMotifs, AssNoms)
  tmp$nom <- as.character(tmp$nom)

  # compute the number of necessary columns
  ttp <- list()
  for (i in seq_along(tmp$nom))
    ttp[[i]] <- rownames(mOccur)[which(AssNoms == tmp$nom[i])]

  if (opt.sort == TRUE)
    for (i in seq_along(tmp$nom)) ttp[[i]] <- sort(ttp[[i]])

  # adjust the length of rownames
  lmax <- max(nchar(tmp$nom))
  for (i in seq_along(tmp$nom)) {
    a <- tmp$nom[i]
    while (nchar(a) < lmax) a <- paste(" ", a, sep = "")
    names(ttp[[i]])[1] <- a
  }

  blanc <- NULL
  for (k in 1:(lmax + 5)) blanc <- paste(blanc, " ", sep = "")

  # concat the matrix by row
  tab <- c("", "Motifs content : (sorted by decreasing mean of Fobs)", "")
  for (i in seq_along(tmp$nom)) {

    tpp    <- concat.byline(ttp[[i]])
    tpp[1] <- paste(names(ttp[[i]])[1], " = { ", tpp[1], sep = "")
    if (length(tpp) > 1)
      for (j in 2:length(tpp)) tpp[j] <- paste(blanc, tpp[j], sep = "")
    tpp[length(tpp)] <- paste(delstr.end(tpp[length(tpp)], 1), " }", sep = "")

    tab  <- c(tab, tpp)
  }
  tab  <- c(tab, "")

  # plot the legend
  plot.bypage(tab)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute Statistiques for each Assembly Motif
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute.mStats.byMotif <- function(res, nbOpt, AssNoms) {

  setMots <- unique(res$mMotifs[nbOpt, ])
  nbMots  <- length(setMots)

  motR2cal <- motR2prd <- motSlope <- matrix(0, nrow = nbOpt, ncol = nbMots)
  colnames(motR2cal) <- colnames(motR2prd) <-
                        colnames(motSlope) <- unique(AssNoms)

  nbcl <- 1
  motR2cal[nbcl, ] <- R2mse(res$tCal[nbcl, ], res$fct)
  motR2prd[nbcl, ] <- R2mse(res$tPrd[nbcl, ], res$fct)
  motSlope[nbcl, ] <- lm(res$tCal[nbcl, ] ~ res$fct)$coefficients[2]

  if (nbOpt > 1) for (nbcl in 2:nbOpt) {

    motR2cal[nbcl, ] <- motR2cal[nbcl - 1, ]
    motR2prd[nbcl, ] <- motR2prd[nbcl - 1, ]
    motSlope[nbcl, ] <- motSlope[nbcl - 1, ]

    indAss <- which(res$tNbcl[nbcl, ] == nbcl)
    indMot <- unique(res$mMotifs[nbcl, indAss])

    for (mot in 1:length(indMot)) {

      motif   <- indMot[mot]

      indAss2 <- which(res$mMotifs[nbcl, ] == motif)
      indMot2 <- unique(res$mMotifs[nbOpt, indAss2])

      motR2cal[nbcl, indMot2] <-
        R2mse(res$tCal[nbcl, indAss2], res$fct[indAss2])
      motR2prd[nbcl, indMot2] <-
        R2mse(res$tPrd[nbcl, indAss2], res$fct[indAss2])
      motSlope[nbcl, indMot2] <-
        lm(res$tCal[nbcl, indAss2] ~ res$fct[indAss2])$coefficients[2]
    }
  }

  resMot        <- list(motR2cal, motR2prd, motSlope)
  names(resMot) <- c("R2cal", "R2prd", "Slope")

  return(resMot)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot all the necessary tree figures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.clustering <- function(tree, res, mOccur,

                            short.names = colnames(mOccur),
                            nbLev       = rep(1, length(short.names)),

                            mSignifElt  = NULL,

                            nbOpt    = NULL,
                            nbmin    = 3,

                            col.tree = "black",
                            xwdw     = c(1:20),
                            ywdw     = c(1:15),
                            titre    = "",

                            tre.cal  = FALSE,
                            tre.prd  = FALSE,
                            tre.pub  = FALSE,
                            tre.leg  = FALSE,
                            tre.detail = FALSE,

                            mot.cal  = FALSE,
                            mot.prd  = FALSE,
                            mot.pub  = FALSE,
                            mot.leg  = FALSE,
                            mot.tre  = FALSE,

                            elt.cal  = FALSE,

                            opt.aov  = FALSE, pvalue = 0.05,
                            opt.all  = FALSE  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  # Check opt.all

  nbElt <- dim(tree$aff)[1]

  if (opt.all == TRUE) {
      tre.cal <- tre.prd <- tre.pub <- tre.leg <- tre.detail <-
      mot.cal <- mot.prd <- mot.pub <- mot.leg <- mot.tre <-
      elt.cal <- TRUE
  }

  # Check tree

  tree.prd     <- tree
  tree.prd$cor <- res$tStats[ ,"R2cal"]                            #res.prd

  nbMax <- first.argmax(res$tStats[ ,"R2prd"])
  if (length(nbOpt) == 0) {
    nbOpt <- first.argmin(res$tStats[ ,"AICc"])
  } else {
    if (nbOpt > nbMax) nbOpt <- nbMax
  }

  cutline   <- amean(tree.prd$cor[(nbOpt - 1):nbOpt])

  affectElt <- cut.tree(tree, nbOpt)
  AssMotifs <- affect.motifs(affectElt, mOccur)
  AssNoms   <- names.assembly(affectElt, mOccur)$motifs
  nbAss     <- length(AssMotifs)


  # Compute different ordre-indices

  tmp       <- sort.tree(tree, index.return = TRUE)
  nom.tree  <- tmp$noms
  ind.tree  <- tmp$ix

  tmp       <- sort.split.names(affectElt, ind.tree, short.names, nbLev)
  ind.elt   <- tmp$index.elt
  ind.split <- tmp$index.split

  lwdw <- length(which(affectElt > 1))
  if (length(xwdw) > length(affectElt)) xwdw <- seq_along(affectElt)
  if (length(xwdw) < lwdw + nbmin)      xwdw <- seq_len(lwdw + nbmin)

  twdw      <- logical(length(affectElt))
  twdw[ind.tree][xwdw] <- TRUE
  stree     <- simplify.tree(tree, twdw)


  # Compute color vectors

    if ((length(col.tree) == 1) & (col.tree[1] == "black")) {
    cols        <- couleurs[shift.affectElt(affectElt)]
    names(cols) <- names(affectElt)
  } else {
    cols        <- check.symbol(col.tree, nbElt)[1:nbElt]
    names(cols) <- names(affectElt)
  }


  # check opt.xpr

  opt.xpr <- FALSE
  setXpr  <- unique(res$xpr)
  if (length(setXpr) != 1) {
    opt.xpr <- TRUE
    setAss  <- res$names[which(res$xpr == setXpr[1])]
    pas     <- length(setAss)
  }

  # check mSignifElt

  if (length(mSignifElt) == 0) {
    tre.detail <- FALSE
    mSignifElt <- matrix(TRUE, nrow = nbElt, ncol = nbElt)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the whole calibration Tree as a diagram
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.cal == TRUE) {

    # plot full tree

    if (tre.detail == TRUE)
      for (nb in 2:(nbOpt - 1)) {
        sig <- mSignifElt[nb, ]
        plot.tree(tree.cal, sig,
                  col   = cols,
                  titre = paste0(titre, "   nb = ", nb)  )
        plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                              opt.sort = TRUE)
      }

    sig <- mSignifElt[nbOpt, ]
    plot.tree(tree.cal, sig,
              col = cols,
              titre = titre)
    lines(x = c(0, nbElt + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, nbElt + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")
    plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                          opt.sort = TRUE)

    # plot head of tree

    stree <- simplify.tree(tree.cal, twdw)
    if (tre.detail == TRUE)
      for (nb in 2:(nbOpt - 1)) {
        sig <- mSignifElt[nb, twdw]
        plot.tree(stree, sig,
                  col   = cols[twdw],
                  titre = paste0(titre, "   nb = ", nb)  )
        plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                              opt.sort = TRUE)
      }

    sig <- mSignifElt[nbOpt, twdw]
    plot.tree(stree, sig,
              col   = cols[twdw],
              titre = paste0(titre, " x-window")  )
    lines(x = c(0, length(twdw) + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, length(twdw) + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")
    plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                          opt.sort = TRUE)

  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the whole prediction Tree as a diagram
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.prd == TRUE) {

    # plot full tree

    if (tre.detail == TRUE)
      for (nb in 2:(nbOpt - 1)) {
        sig <- mSignifElt[nb, ]
        plot.tree(tree.prd, sig,
                  col   = cols,
                  titre = paste0(titre, "   nb = ", nb)  )
        plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                              opt.sort = TRUE)
      }

    sig <- mSignifElt[nbOpt, ]
    plot.tree(tree.prd, sig,
              col = cols,
              titre = titre)
    lines(x = c(0, nbElt + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, nbElt + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")
    plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                          opt.sort = TRUE)

    # plot head of tree

    stree <- simplify.tree(tree.prd, twdw)
    if (tre.detail == TRUE)
      for (nb in 2:(nbOpt - 1)) {
        sig <- mSignifElt[nb, twdw]
        plot.tree(stree, sig,
                  col   = cols[twdw],
                  titre = paste0(titre, "   nb = ", nb)  )
        plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                              opt.sort = TRUE)
      }

    sig <- mSignifElt[nbOpt, twdw]
    plot.tree(stree, sig,
              col   = cols[twdw],
              titre = paste0(titre, " x-window")  )
    lines(x = c(0, length(twdw) + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, length(twdw) + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")
    plot.clusters.content(affectElt[sig], mOccur[ , sig, drop = FALSE],
                          opt.sort = TRUE)

  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the prediction Tree as a monochrom diagram for publication
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.pub == TRUE) {

    plot.tree(tree.prd, mSignifElt[nbOpt, ],
              col = "black", titre = titre)
    lines(x = c(0, nbElt + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, nbElt + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")

    stree <- simplify.tree(tree.prd, twdw)
    plot.tree(stree, mSignifElt[nbOpt, twdw],
              col   = "black",
              titre = paste0(titre, " x-window")   )
    lines(x = c(0, length(twdw) + 1),
          y = rep(cutline, 2), col = "blue", lty = "dashed")
    lines(x = c(0, length(twdw) + 1),
          y = rep(res$tStats[nbOpt, "R2prd"], 2), col = "red")
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the motifs as boxplots
  # If opt.xpr is TRUE, plot a boxplot by experiment
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  opt.aov        <- TRUE
  opt.decreasing <- FALSE
  opt.horizontal <- TRUE


  if (mot.cal == TRUE) {

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    #  Function drawn is the observed function
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre2  <- paste(titre, "calibration", sep = ".")

    fct     <- res$fct
    ylim    <- range(fct)

    ref.all <- sort.motifs(fct, AssMotifs, AssNoms,
                           opt.aov = opt.aov, pvalue = pvalue,
                           opt.dec = opt.decreasing  )
    ordre   <- ref.all$motif

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    if (opt.xpr == TRUE) {

      titre3  <- paste(titre2, "all.Xpr", sep = ".")
      plot.box.motifs(fct, AssMotifs, ref.all,
                      ordre   = ordre,
                      ylim    = ylim,
                      titre   = titre3,
                      opt.aov = opt.aov, pvalue = pvalue,
                      opt.hor = opt.horizontal  )

      ref.xpr <- list()
      for (ipr in seq_along(setXpr)) {

        titre3         <- paste(titre2, setXpr[ipr], sep = ".")
        index          <- which(res$xpr == setXpr[ipr])
        ref.xpr[[ipr]] <- sort.motifs(fct[index], AssMotifs[index],
                                      AssNoms[index],
                                      opt.aov = opt.aov, pvalue = pvalue,
                                      opt.dec = opt.decreasing  )
        plot.box.motifs(fct[index], AssMotifs[index], ref.xpr[[ipr]],
                        ordre   = ordre,
                        ylim    = ylim,
                        titre   = titre3,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = opt.horizontal  )
      }

    } else {

      plot.box.motifs(fct, AssMotifs, ref.all,
                      ordre   = ordre,
                      ylim    = ylim,
                      titre   = titre2,
                      opt.aov = opt.aov, pvalue = pvalue,
                      opt.hor = opt.horizontal  )
    }


    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    #  Function drawn is the function predicted by the clustering model
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    if (mot.prd == TRUE) {

      titre2 <- paste(titre, "prediction", sep = ".")

      fct     <- res$tPrd[nbOpt, ]
      ref.all <- sort.motifs(fct, AssMotifs, AssNoms,
                             opt.aov = opt.aov, pvalue = pvalue,
                             opt.dec = opt.decreasing  )

      #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
      if (opt.xpr == TRUE) {

        titre3 <- paste(titre2, "all.Xpr", sep = ".")
        plot.box.motifs(fct, AssMotifs, sres = ref.all,
                        ordre   = ordre,
                        ylim    = ylim,
                        titre   = titre3,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = opt.horizontal  )

        for (ipr in seq_along(setXpr)) {

          titre3         <- paste(titre2, setXpr[ipr], sep = ".")
          index          <- which(res$xpr == setXpr[ipr])
          ref.xpr[[ipr]] <- sort.motifs(fct[index], AssMotifs[index],
                                        AssNoms[index],
                                        opt.aov = opt.aov, pvalue = pvalue,
                                        opt.dec = opt.decreasing  )
          plot.box.motifs(fct[index], AssMotifs[index], ref.xpr[[ipr]],
                          ordre   = ordre,
                          ylim    = ylim,
                          titre   = titre3,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = opt.horizontal  )
        }

      } else {

        plot.box.motifs(fct, AssMotifs, ref.all,
                        ordre   = ordre,
                        ylim    = ylim,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = opt.horizontal  )
      }
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the assembly content of each motif
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (mot.leg == TRUE) {

    if (opt.xpr == TRUE) {

      index <- which(res$xpr == setXpr[1])
      plot.motifs.content(res$fct[index], AssMotifs[index], AssNoms[index])

    } else {

      plot.motifs.content(res$fct, AssMotifs, AssNoms)

    }
  }




  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted by increasing effect on assembly function
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = TRUE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  opt.aov        <- TRUE
  opt.decreasing <- FALSE
  opt.horizontal <- TRUE


  if (elt.cal == TRUE) {

    fct     <- res$fct
    ylim    <- range(fct)

    ref.all <- sort.elements(fct, mOccur, cols,
                             pvalue = pvalue, opt.dec = FALSE  )
    ordre   <- ref.all$index

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")
      plot.box.elements(fct, mOccur, ref.all,
                        ordre   = ordre,
                        ylim    = ylim,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      ref.xpr <- list()
      for (ipr in seq_along(setXpr)) {

        titre2         <- paste(titre, setXpr[ipr], sep = ".")
        index          <- which(res$xpr == setXpr[ipr])
        ref.xpr[[ipr]] <- sort.elements(fct[index],
                                    mOccur[index, , drop = FALSE], cols,
                                    pvalue = pvalue, opt.dec = FALSE  )
        plot.box.elements(fct[index], mOccur[index, , drop = FALSE],
                          sres    = ref.xpr[[ipr]],
                          ordre   = ordre,
                          ylim    = ylim,
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }


      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(fct, mOccur, ref.all,
                        ordre   = ordre,
                        elt.wdw = rev(nbElt - ywdw + 1),
                        ylim    = ylim,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      for (ipr in seq_along(setXpr)) {
        titre2 <- paste(titre, "all.Xpr.zoom   Higher part",
                        setXpr[ipr], sep = ".")
        index <- which(res$xpr == setXpr[ipr])
        plot.box.elements(fct[index], mOccur[index, , drop = FALSE],
                          sres    = ref.xpr[[ipr]],
                          ordre   = ordre,
                          elt.wdw = rev(nbElt - ywdw + 1),
                          ylim    = ylim,
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(fct, mOccur, ref.all,
                        ordre   = ordre,
                        ylim    = ylim,
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom", sep = ".")
      plot.box.elements(fct, mOccur, ref.all,
                        ordre   = ordre,
                        ylim    = ylim,
                        elt.wdw = rev(nbElt - ywdw + 1),
                        titre   = paste(titre2, "higher part", sep = " "),
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted as the clustering tree
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = TRUE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (elt.cal == TRUE) {

    indOrd  <- ref.all$ind.inturn[rev(ind.tree)]

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      for (ipr in seq_along(setXpr)) {

        titre2     <- paste(titre, setXpr[ipr], sep = ".")
        index      <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }

      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = rev(nbElt - ywdw + 1),
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      for (ipr in seq_along(setXpr)) {

        titre2 <- paste(titre, "all.Xpr.zoom   Higher part",
                        setXpr[ipr], sep = ".")
        index  <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          elt.wdw = rev(nbElt - ywdw + 1),
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom  Higher part", sep = ".")

      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = rev(nbElt - ywdw + 1),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted by element as the clustering tree
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = TRUE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (elt.cal == TRUE) {

    indOrd  <- ref.all$ind.inturn[rev(ind.split)]

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      for (ipr in seq_along(setXpr)) {

        titre2     <- paste(titre, setXpr[ipr], sep = ".")
        index      <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }

      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = rev(nbElt - ywdw + 1),
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      for (ipr in seq_along(setXpr)) {

        titre2 <- paste(titre, "all.Xpr.zoom   Higher part",
                        setXpr[ipr], sep = ".")
        index  <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          elt.wdw = rev(nbElt - ywdw + 1),
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = TRUE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )

      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom  Higher part", sep = ".")

      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = rev(nbElt - ywdw + 1),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = TRUE  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted by increasing effect on assembly function
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = FALSE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  opt.aov        <- TRUE
  opt.decreasing <- TRUE
  opt.horizontal <- FALSE


  if (elt.cal == TRUE) {

    ref.all <- sort.elements(res$fct, mOccur, cols,
                             pvalue = pvalue, opt.dec = TRUE  )

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index,
                        sres    = ref.all,
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      ref.xpr <- list()
      for (ipr in seq_along(setXpr)) {

        titre2     <- paste(titre, setXpr[ipr], sep = ".")
        index      <- which(res$xpr == setXpr[ipr])

        ref.xpr[[ipr]] <- sort.elements(res$fct[index],
                                        mOccur[index, , drop = FALSE], cols,
                                        pvalue = pvalue, opt.dec = TRUE  )

        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index,
                          sres    = ref.xpr[[ipr]],
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }

      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index,
                        sres    = ref.all,
                        elt.wdw = ywdw,
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      for (ipr in seq_along(setXpr)) {
        titre2  <- paste(titre, "all.Xpr.zoom   Higher part",
                         setXpr[ipr], sep = ".")
        index   <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index,
                          sres    = ref.xpr[[ipr]],
                          elt.wdw = ywdw,
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index,
                        sres    = ref.all,
                        ylim    = range(res$fct),
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )


      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom  Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index,
                        sres    = ref.all,
                        elt.wdw = ywdw,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted as clustering tree
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = FALSE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (elt.cal == TRUE) {

    indOrd  <- ref.all$ind.inturn[ind.tree]

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")

      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      for (ipr in seq_along(setXpr)) {

        titre2     <- paste(titre, setXpr[ipr], sep = ".")
        index      <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }


      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = ywdw,
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      for (ipr in seq_along(setXpr)) {
        titre2  <- paste(titre, "all.Xpr.zoom   Higher part",
                         setXpr[ipr], sep = ".")
        index   <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          elt.wdw = ywdw,
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )


      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom  Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = ywdw,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the element performances
  #            sorted by element as the clustering tree
  #                                 as a boxplot for comparison with motifs
  #                  opt.horizontal = FALSE
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (elt.cal == TRUE) {

    indOrd  <- ref.all$ind.inturn[ind.split]

    if (opt.xpr == TRUE) {

      # plot all the figure
      titre2  <- paste(titre, "all.Xpr", sep = ".")

      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      for (ipr in seq_along(setXpr)) {

        titre2     <- paste(titre, setXpr[ipr], sep = ".")
        index      <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }


      # plot only the higher part of the figure
      titre2 <- paste(titre, "all.Xpr.zoom   Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = ywdw,
                        ylim    = range(res$fct),
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )

      for (ipr in seq_along(setXpr)) {
        titre2  <- paste(titre, "all.Xpr.zoom   Higher part",
                         setXpr[ipr], sep = ".")
        index   <- which(res$xpr == setXpr[ipr])
        plot.box.elements(res$fct[index], mOccur[index, , drop = FALSE],
                          ordre   = ref.all$index[indOrd],
                          sres    = ref.xpr[[ipr]][indOrd, ],
                          elt.wdw = ywdw,
                          ylim    = range(res$fct),
                          titre   = titre2,
                          opt.aov = opt.aov, pvalue = pvalue,
                          opt.hor = FALSE  )
      }

    } else {

      # plot all the figure
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        ylim    = range(res$fct),
                        titre   = titre,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )


      # plot only the higher part of the figure
      titre2 <- paste(titre, "zoom  Higher part", sep = ".")
      plot.box.elements(res$fct, mOccur,
                        ordre   = ref.all$index[indOrd],
                        sres    = ref.all[indOrd, ],
                        elt.wdw = ywdw,
                        titre   = titre2,
                        opt.aov = opt.aov, pvalue = pvalue,
                        opt.hor = FALSE  )
    }
  }



  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Compute Statistiques for each Assembly Motif
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#  mStats <- compute.mStats.byMotif(res, nbOpt, AssNoms)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the hierarchical tree of motifs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  mot.tre <- FALSE

  if (mot.tre == TRUE) {

    index <- which(res$xpr == setXpr[1])
    tree.mot <- cluster.motifs(affectElt, mOccur[index, , drop = FALSE],
                               Fobs[index],
                               opt.meth = "divisive")
    plot.tree(tree.mot, col = "black", titre)
  }

}





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                             END of file myPLOT.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
