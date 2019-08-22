#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#        SCRIPT for Complete Standard Treatment of Dynamic DataSet
#
#                             Benoît JAILLARD
#                              Janvier 2018
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
affectElt <- rep(1, nbElt)


dat <- 1
for (dat in seq_along(setData)) {

  ALPHA <- BETA <- fOBS <- FOBS <- FOBSF <-
           FMONO <- FSCALE <- MOCCUR <- XPR <- NULL

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nom    <- setData[dat]

  nomfic <- paste(rad, nom, sep = "/")
  if (!file.exists(nomfic)) dir.create(nomfic)

  filename <- paste(rad, paste(nom, "csv", sep = "."), sep = "/")
  data     <- read.table(filename, header = TRUE, sep = ",")

  setXpr  <- colnames(data)[(1 + nbElt + indXpr[1]):
                              (1 + nbElt + indXpr[length(indXpr)])]


#  filename <- paste(rad, paste(paste("eco", as.character(dat), sep = ""),
#                               "csv", sep = "."), sep = "/")
#  trait    <- as.vector(unlist(read.table(filename, header = FALSE, sep = ",")))
#  tmp      <- sort(trait, decreasing = TRUE, index.return = TRUE)
#  trait    <- as.vector(tmp$x)

  ipr <- indXpr[1]
  for (ipr in seq_along(setXpr)) {

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Matrices
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # Read of raw data
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    mOccur           <- as.matrix(data[ , 1 + (1:nbElt)])
    rownames(mOccur) <- as.character(data[ , 1])

    Fobs   <- as.vector(as.numeric(unlist(data[ , (1 + nbElt) + ipr])))
    size   <- apply(mOccur, MARGIN = 1, FUN = sum)

    mOccur <- mOccur[(size != 0),]
    Fobs   <- Fobs[(size != 0)]
    size   <- size[(size != 0)]

    tmp    <- rm.dual.assemblies(mOccur, Fobs)
    mOccur <- tmp$mat
    Fobs   <- tmp$fct

    nbAss    <- length(Fobs)
    xpr      <- rep(setXpr[ipr], nbAss)

    figures  <- check.symbol(figures,  max(nbAss, nbElt))
    couleurs <- check.symbol(couleurs, max(nbAss, nbElt))

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

    source("nyBeginning.script.R")

    nbAss  <- length(Fobs)
    xpr    <- rep(setXpr[ipr], nbAss)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("Separating", setXpr[ipr], sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)

    source("nySeparating.script.R")

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("Predicting", setXpr[ipr], sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)

    source("nyPredicting.script.R")

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("Fobs", setXpr[ipr], sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)
    open.txt(enregistre, filename)

    tree.cal <- cluster.elements(mOccur, Fobs, opt.mean = "gmean")
    write.tree(filename, tree.cal)
    tree.cal <- read.tree(filename)

    res.prd  <- predict.function(tree.cal, mOccur, Fobs, opt.mean = "gmean")

    nbOpt    <- first.argmin(res.prd$tStats[, "AICc"])

    plot.prediction(res.prd,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    tre.prd = TRUE, col = couleurs[gf])

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("alpha", setXpr[ipr], sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)
    open.txt(enregistre, filename)

    tree.cal <- cluster.elements(mOccur, alpha, opt.mean = "gmean")
    write.tree(filename, tree.cal)
    tree.cal <- read.tree(filename)

    res.prd <- predict.function(tree.cal, mOccur, alpha, opt.mean = "gmean")

    plot.prediction(res.prd,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    tre.prd = TRUE, col = couleurs[gf])

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("beta", setXpr[ipr], sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)
    open.txt(enregistre, filename)

    tree.cal <- cluster.elements(mOccur, beta, opt.mean = "amean")
    write.tree(filename, tree.cal)
    tree.cal <- read.tree(filename)

    res.prd  <- predict.function(tree.cal, mOccur, beta, opt.mean = "amean")

    plot.prediction(res.prd,
                    nbOpt   = 3,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    nbOpt   = 3,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    nbOpt   = 3,
                    titre   = titre,
                    opt.aov = TRUE, pvalue = 0.001,
                    tre.prd = TRUE, col = couleurs[gf])

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # Preparing plotting All experiments over time

    ALPHA     <- c(ALPHA, alpha)
    BETA      <- c(BETA,  beta)
    FOBSF     <- c(FOBSF, Fobs)
    FOBS      <- c(FOBS,  Fobs / amean(Fobs))
    fOBS      <- c(fOBS,  fobs)
    MOCCUR    <- rbind(MOCCUR, mOccur)
    FMONO     <- c(FMONO, Fmono)
    FSCALE    <- c(FSCALE, rep(Fscale, nbAss))
    XPR       <- c(XPR,   xpr)
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #    All experiments over time
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  alpha  <- ALPHA
  beta   <- BETA
  Fobs   <- FOBSF
  fobs   <- fOBS
  mOccur <- MOCCUR
  Fmono  <- FMONO
  Fscale <- unique(FSCALE)
  xpr    <- XPR

  #  tmp              <- paste(rownames(mOccur), paste0("t", XPR), sep = ".")
#  rownames(MOCCUR) <- c(rownames(mOccur),
#                        tmp[(length(rownames(mOccur)) + 1):length(XPR)])
  nbAss            <- length(XPR)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.Separating")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)

  source("nySeparating.script.R")

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.Predicting")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)

  source("nyPredicting.script.R")

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.FOBSF")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(MOCCUR, FOBS, xpr = XPR, opt.mean = "gmean")
  write.tree(filename, tree.cal)
  tree.cal <- read.tree(filename)

  res.prd <- predict.function(tree.cal, MOCCUR, FOBS,
                              xpr = XPR, opt.mean = "gmean")

  plot.prediction(res.prd,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  tre.prd = TRUE, col = couleurs[gf])

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.fOBS")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(MOCCUR, fOBS, xpr = XPR, opt.mean = "gmean")
  write.tree(filename, tree.cal)
  tree.cal <- read.tree(filename)

  res.prd <- predict.function(tree.cal, MOCCUR, fOBS,
                              xpr = XPR, opt.mean = "gmean")

  plot.prediction(res.prd,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  tre.prd = TRUE, col = couleurs[gf])

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.ALPHA")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(MOCCUR, ALPHA, xpr = XPR, opt.mean = "gmean")
  write.tree(filename, tree.cal)
  tree.cal <- read.tree(filename)

  res.prd <- predict.function(tree.cal, MOCCUR, ALPHA,
                              xpr = XPR, opt.mean = "gmean")

  plot.prediction(res.prd,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  tre.prd = TRUE, col = couleurs[gf])

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.BETA")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(MOCCUR, BETA, xpr = XPR, opt.mean = "amean")
  write.tree(filename, tree.cal)
  tree.cal <- read.tree(filename)

  res.prd <- predict.function(tree.cal, MOCCUR, BETA,
                              xpr = XPR, opt.mean = "amean")

  plot.prediction(res.prd,
                  nbOpt   = 4,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  nbOpt   = 4,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, MOCCUR,
                  nbOpt   = 4,
                  titre   = titre,
                  opt.aov = TRUE, pvalue = 0.001,
                  tre.prd = TRUE, col = couleurs[gf])

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #    Dynamics
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  titre <- c("all.Dynamics")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)

  plot.dynamics(ALPHA, BETA, MOCCUR, xpr = XPR,
                figs = gf, cols = gf)

#  plot(x = amean.Fmono, xlab = "Fmono",
#       y = amean.Fobs, ylab = "Fobs",
#       pch = figures[1], col = "black", bg = "white",
#       type = "p", cex = 2, tck = 0.02, las = 1)

#  points.sd(x = Fmono, y = Fobs,
#           figure = figures[1], couleur = "black",
#           nom = "", opt.std = TRUE, scale.log = FALSE)


  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  END of FILE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

