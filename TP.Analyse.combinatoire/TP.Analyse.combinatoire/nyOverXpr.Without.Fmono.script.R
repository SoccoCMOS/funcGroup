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


dat <- 1
for (dat in seq_along(setData)) {

  FOBS <- FOBSF <- MOCCUR <- TIM <- NULL

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nom    <- setData[dat]

  nomfic <- paste(rad, nom, sep = "/")
  if (!file.exists(nomfic)) dir.create(nomfic)

  filename <- paste(rad, paste(nom, "csv", sep = "."), sep = "/")
  data     <- read.table(filename, header = TRUE, sep = ",")

#  filename <- paste(rad, paste(paste("eco", as.character(dat), sep = ""),
#                               "csv", sep = "."), sep = "/")
#  trait    <- as.vector(unlist(read.table(filename, header = FALSE, sep = ",")))
#  tmp      <- sort(trait, decreasing = TRUE, index.return = TRUE)
#  trait    <- as.vector(tmp$x)

  tim <- 1
  for (tim in seq_len(nbTime)) {

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Matrices
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # Read of raw data
    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    mOccur           <- as.matrix(data[ , 1 + (1:nbElt)])
    rownames(mOccur) <- as.character(data[ , 1])

    Fobs   <- as.vector(as.numeric(unlist(data[ , (1 + nbElt) + tim])))
    size   <- apply(mOccur, MARGIN = 1, FUN = sum)

    figures  <- check.symbol(figures, length(Fobs))
    couleurs <- check.symbol(couleurs, length(Fobs))

    mOccur <- mOccur[(size != 0),]
    Fobs   <- Fobs[(size != 0)]
    size   <- size[(size != 0)]

    mOccur <- mOccur[(size != 1),]
    Fobs   <- Fobs[(size != 1)]
    size   <- size[(size != 1)]

    nbAss <- dim(mOccur)[1]

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    titre <- paste("Fobs", paste0("t", tim), sep = ".")

    filename <- paste(rad, nom, titre, sep = "/")
    open.pdf(enregistre, filename)
    open.txt(enregistre, filename)

    tree.cal <- cluster.elements(mOccur, Fobs)

    res.prd  <- predict.function(tree.cal, mOccur, Fobs)

    plot.prediction(res.prd,
                    titre   = titre,
                    opt.aov = FALSE, pvalue = 0.05,
                    opt.all = TRUE)

    plot.clustering(tree.cal, res.prd, mOccur,
                    titre   = titre,
                    opt.aov = FALSE, pvalue = 0.05,
                    opt.all = TRUE)

    close.pdf(enregistre)

    #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
    # Preparing plotting All experiments over time

    FOBSF     <- c(FOBSF, Fobs)
    FOBS      <- c(FOBS,  Fobs / amean(Fobs))
    MOCCUR    <- rbind(MOCCUR, mOccur)
    TIM       <- c(TIM,   rep(tim, nbAss))

  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #    All experiments over time
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  FobsF  <- FOBSF
  Fobs   <- FOBS
  mOccur <- MOCCUR

  nbAss <- dim(mOccur)[1]


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.FobsF.raw")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(mOccur, FobsF, xpr = TIM)

  res.prd <- predict.function(tree.cal, mOccur, FobsF, xpr = TIM)

  plot.prediction(res.prd,
                  titre   = titre,
                  opt.aov = FALSE, pvalue = 0.05,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, mOccur,
#                  Fmono   = Fmono,
                  titre   = titre,
                  opt.aov = FALSE, pvalue = 0.05,
                  opt.all = TRUE)

  close.pdf(enregistre)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  titre <- c("all.FobsF")

  filename <- paste(rad, nom, titre, sep = "/")
  open.pdf(enregistre, filename)
  open.txt(enregistre, filename)

  tree.cal <- cluster.elements(mOccur, Fobs, xpr = TIM)

  res.prd <- predict.function(tree.cal, mOccur, Fobs, xpr = TIM)

  plot.prediction(res.prd,
                  titre   = titre,
                  opt.aov = FALSE, pvalue = 0.05,
                  opt.all = TRUE)

  plot.clustering(tree.cal, res.prd, mOccur,
                  titre   = titre,
                  opt.aov = FALSE, pvalue = 0.05,
                  opt.all = TRUE)

  close.pdf(enregistre)

 #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  END of FILE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
