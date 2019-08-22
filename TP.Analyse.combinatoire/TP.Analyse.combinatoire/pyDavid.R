#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Analysis year-by-year of datset from Tilman
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  Loading of usefull files and data
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list = ls())


source("myLibraries.R")

source("myTools.R")
source("myStats.R")
source("myTree.R")
source("myCombinat.R")
source("myPlot.R")

source("mySeparating.R")
source("myPredicting.R")
source("myClustering.R")





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#    MAIN
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rad       <- c("David")

setData   <- c("f120.unsorted.moy")

nbElt     <- 16
nbXpr     <- 11

gf        <- c(1,3,2,4,1,3,3,2,1,2,1,4,2,3,4,4)

enregistre <- TRUE
#enregistre <- FALSE


indXpr <- c(1:8)
source("nyOverXpr.withFmono.script.1.R")


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  END of FILE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxm



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  Conventional species clustering
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


affectGf <- gf
names(affectGf) <- colnames(MOCCUR)
nbOpt    <- length(unique(affectGf))

opt.aov <- TRUE
pvalue  <- 0.001
opt.all <- TRUE

predict.function.gf <- function(affectElt, mOccur, Fobs,
                                xpr = rep(1, length(Fobs)),
                                # options for computing
                                opt.mean = "amean",
                                opt.mod  = "byelt",
                                opt.jack = FALSE,  jack = c(2,5)) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt     <- dim(mOccur)[2]
  nbAss     <- length(Fobs)

  assMotifs <- affect.motifs(affectElt, mOccur)
  assNoms   <- names.assembly(affectElt, mOccur)$motifs
  nbAss     <- length(assMotifs)

  # Compute the cross-validation predictions
  Fcal <- predict.cal(assMotifs, mOccur, Fobs, xpr, opt.mean, opt.mod)
  Fprd <- predict.prd(assMotifs, mOccur, Fobs, xpr,
                      opt.mean, opt.mod, opt.jack, jack)

  # Compute the associated statistiques
  nbK <- length(unique(assMotifs))

  tmp    <- c("missing", "R2cal", "R2prd", "AIC", "AICc", "pcal", "pprd")
  vStats <- numeric(length(tmp))
  names(vStats) <- tmp

  vStats["missing"] <- length(na.action(na.omit(Fprd))) / nbAss
  vStats["R2cal"]   <- R2mse(Fcal, Fobs)
  vStats["R2prd"]   <- R2mse(Fprd, Fobs)
  vStats["AIC"]     <- AIC(  Fcal, Fobs, nbK)
  vStats["AICc"]    <- AICc( Fcal, Fobs, nbK)

  if (nbK > 1) {
    vStats["pcal"]  <- pmse( Fcal, Fobs, nbK)
    vStats["pprd"]  <- pmse( Fprd, Fobs, nbK)
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  res <- list(rownames(mOccur), Fobs, xpr, opt.mean, opt.mod,
              Fcal, Fprd, vStats, assMotifs, assNoms)

  names(res) <- c("names", "Fobs", "xpr", "opt.mean", "opt.mod",
                  "Fcal", "Fprd", "vStats", "vMotifs", "vNames")

  return(res)
}

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#   ALPHA
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

titre <- "all.ALPHA.gf"
Fobs <- ALPHA
opt.mean <- "gmean"

filename <- paste(rad, nom, titre, sep = "/")
open.pdf(enregistre, filename)

res.prd <- predict.function.gf(affectGf, MOCCUR, Fobs,
                               xpr = XPR,
                               # options for computing
                               opt.mean = opt.mean,
                               opt.mod  = "byelt")

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction.LOO(res.prd$Fcal, res.prd$Fprd, res.prd$Fobs, res.prd$vMotifs,
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)


plot.prediction.LOO(res.prd$Fcal, res.prd$Fprd, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.LOO(res.prd$Fprd, res.prd$Fcal, res.prd$Fobs, res.prd$vMotifs,
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.LOO(res.prd$Fprd, res.prd$Fcal, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.simple(res.prd$Fprd, res.prd$Fobs, res.prd$vMotifs,
                       nbOpt,
                       titre = titre, opt.mean = opt.mean,
                       opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.simple(res.prd$Fprd, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                       nbOpt,
                       titre = titre, opt.mean = opt.mean,
                       opt.aov  = TRUE, pvalue = 0.001)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
opt.aov        <- TRUE
opt.decreasing <- FALSE
opt.horizontal <- TRUE

assMotifs <- res.prd$vMotifs
assNoms   <- res.prd$vNames

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  Function drawn is the observed function
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre2  <- paste(titre, "calibration", sep = ".")

fct     <- res.prd$Fobs
ylim    <- range(fct)

ref.all <- sort.motifs(fct, assMotifs, assNoms,
                       opt.aov = opt.aov, pvalue = pvalue,
                       opt.dec = opt.decreasing  )
ordre   <- ref.all$motif

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre3  <- paste(titre2, "all.Xpr", sep = ".")
plot.box.motifs(fct, assMotifs, ref.all,
                ordre   = ordre,
                ylim    = ylim,
                titre   = titre3,
                opt.aov = opt.aov, pvalue = pvalue,
                opt.hor = opt.horizontal  )

ref.xpr <- list()
for (ipr in seq_along(setXpr)) {

  titre3         <- paste(titre2, setXpr[ipr], sep = ".")
  index          <- which(xpr == setXpr[ipr])
  ref.xpr[[ipr]] <- sort.motifs(fct[index], assMotifs[index],
                                assNoms[index],
                                opt.aov = opt.aov, pvalue = pvalue,
                                opt.dec = opt.decreasing  )

  plot.box.motifs(fct[index], assMotifs[index], ref.xpr[[ipr]],
                  ordre   = ordre,
                  ylim    = ylim,
                  titre   = titre3,
                  opt.aov = opt.aov, pvalue = pvalue,
                  opt.hor = opt.horizontal  )
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  Function drawn is the function predicted by the clustering model
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

titre2 <- paste(titre, "prediction", sep = ".")

fct     <- res.prd$Fprd
fct[na.action(na.omit(fct))] <- amean(fct)
ref.all <- sort.motifs(fct, assMotifs, assNoms,
                       opt.aov = opt.aov, pvalue = pvalue,
                       opt.dec = opt.decreasing  )


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre3 <- paste(titre2, "all.Xpr", sep = ".")
plot.box.motifs(fct, assMotifs, ref.all,
                ordre   = ordre,
                ylim    = ylim,
                titre   = titre3,
                opt.aov = opt.aov, pvalue = pvalue,
                opt.hor = opt.horizontal  )

for (ipr in seq_along(setXpr)) {

  titre3         <- paste(titre2, setXpr[ipr], sep = ".")
  index          <- which(xpr == setXpr[ipr])
  ref.xpr[[ipr]] <- sort.motifs(fct[index], assMotifs[index],
                                assNoms[index],
                                opt.aov = opt.aov, pvalue = pvalue,
                                opt.dec = opt.decreasing  )

  plot.box.motifs(fct[index], assMotifs[index], ref.xpr[[ipr]],
                  ordre   = ordre,
                  ylim    = ylim,
                  titre   = titre3,
                  opt.aov = opt.aov, pvalue = pvalue,
                  opt.hor = opt.horizontal  )
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Plot the assembly content of each motif
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

index <- which(xpr == setXpr[1])
plot.motifs.content(res.prd$Fobs[index], assMotifs[index], assNoms[index])

close.pdf(enregistre)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#   BETA
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

titre <- "all.BETA.gf"
Fobs  <- BETA
opt.mean <- "amean"

filename <- paste(rad, nom, titre, sep = "/")
open.pdf(enregistre, filename)

res.prd <- predict.function.gf(affectGf, MOCCUR, Fobs,
                               xpr = XPR,
                               # options for computing
                               opt.mean = opt.mean,
                               opt.mod  = "byelt")

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.prediction.LOO(res.prd$Fcal, res.prd$Fprd, res.prd$Fobs, res.prd$vMotifs,
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)


plot.prediction.LOO(res.prd$Fcal, res.prd$Fprd, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.LOO(res.prd$Fprd, res.prd$Fcal, res.prd$Fobs, res.prd$vMotifs,
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.LOO(res.prd$Fprd, res.prd$Fcal, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                    nbOpt,
                    titre = titre, opt.mean = opt.mean,
                    opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.simple(res.prd$Fprd, res.prd$Fobs, res.prd$vMotifs,
                       nbOpt,
                       titre = titre, opt.mean = opt.mean,
                       opt.aov  = TRUE, pvalue = 0.001)

plot.prediction.simple(res.prd$Fprd, res.prd$Fobs, rep(1, length(res.prd$vMotifs)),
                       nbOpt,
                       titre = titre, opt.mean = opt.mean,
                       opt.aov  = TRUE, pvalue = 0.001)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
opt.aov        <- TRUE
opt.decreasing <- FALSE
opt.horizontal <- TRUE

assMotifs <- res.prd$vMotifs
assNoms   <- res.prd$vNames

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  Function drawn is the observed function
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre2  <- paste(titre, "calibration", sep = ".")

fct     <- res.prd$Fobs
ylim    <- range(fct)

ref.all <- sort.motifs(fct, assMotifs, assNoms,
                       opt.aov = opt.aov, pvalue = pvalue,
                       opt.dec = opt.decreasing  )
ordre   <- ref.all$motif

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre3  <- paste(titre2, "all.Xpr", sep = ".")
plot.box.motifs(fct, assMotifs, ref.all,
                ordre   = ordre,
                ylim    = ylim,
                titre   = titre3,
                opt.aov = opt.aov, pvalue = pvalue,
                opt.hor = opt.horizontal  )

ref.xpr <- list()
for (ipr in seq_along(setXpr)) {

  titre3         <- paste(titre2, setXpr[ipr], sep = ".")
  index          <- which(xpr == setXpr[ipr])
  ref.xpr[[ipr]] <- sort.motifs(fct[index], assMotifs[index],
                                assNoms[index],
                                opt.aov = opt.aov, pvalue = pvalue,
                                opt.dec = opt.decreasing  )

  plot.box.motifs(fct[index], assMotifs[index], ref.xpr[[ipr]],
                  ordre   = ordre,
                  ylim    = ylim,
                  titre   = titre3,
                  opt.aov = opt.aov, pvalue = pvalue,
                  opt.hor = opt.horizontal  )
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#  Function drawn is the function predicted by the clustering model
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

titre2 <- paste(titre, "prediction", sep = ".")

fct     <- res.prd$Fprd
fct[na.action(na.omit(fct))] <- amean(fct)
ref.all <- sort.motifs(fct, assMotifs, assNoms,
                       opt.aov = opt.aov, pvalue = pvalue,
                       opt.dec = opt.decreasing  )


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre3 <- paste(titre2, "all.Xpr", sep = ".")
plot.box.motifs(fct, assMotifs, ref.all,
                ordre   = ordre,
                ylim    = ylim,
                titre   = titre3,
                opt.aov = opt.aov, pvalue = pvalue,
                opt.hor = opt.horizontal  )

for (ipr in seq_along(setXpr)) {

  titre3         <- paste(titre2, setXpr[ipr], sep = ".")
  index          <- which(xpr == setXpr[ipr])
  ref.xpr[[ipr]] <- sort.motifs(fct[index], assMotifs[index],
                                assNoms[index],
                                opt.aov = opt.aov, pvalue = pvalue,
                                opt.dec = opt.decreasing  )

  plot.box.motifs(fct[index], assMotifs[index], ref.xpr[[ipr]],
                  ordre   = ordre,
                  ylim    = ylim,
                  titre   = titre3,
                  opt.aov = opt.aov, pvalue = pvalue,
                  opt.hor = opt.horizontal  )
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Plot the assembly content of each motif
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

index <- which(xpr == setXpr[1])
plot.motifs.content(res.prd$Fobs[index], assMotifs[index], assNoms[index])

close.pdf(enregistre)



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

affectElt <- cut.tree(tree.cal, 4)

affectElt

mean(FMONO)
sd(FMONO)

mean(FMONO[affectElt == 1])
sd(FMONO[affectElt == 1])
range(FMONO[affectElt == 1])

mean(FMONO[affectElt == 2])
sd(FMONO[affectElt == 2])
range(FMONO[affectElt == 2])

mean(FMONO[affectElt == 3])
sd(FMONO[affectElt == 3])
range(FMONO[affectElt == 3])

mean(FMONO[affectElt == 4])
sd(FMONO[affectElt == 4])
range(FMONO[affectElt == 4])

affectAOV <- rep(c(1:16), 8)
summary(lm(FMONO ~ affectAOV))
test.posthoc(FMONO, affectAOV, pvalue = 0.001)

affectAOV <- rep(affectElt, 8)
summary(aov(FMONO ~ affectAOV)) # p = 0.00277
test.posthoc(FMONO, affectAOV, pvalue = 0.001) # 4 (198.3) et 2 (164.2) en a,
#                                                3 (71.9)  et 1 (66.2)  en b



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre <- "all.gf.number"
filename <- paste(rad, nom, titre, sep = "/")
open.pdf(enregistre, filename)

size      <- apply(MOCCUR, MARGIN = 1, FUN = sum)

titre     <- c("all.BETA")
nbOpt     <- 4
filename  <- paste(rad, nom, titre, sep = "/")
tree.cal  <- read.tree(filename)
affectElt <- cut.tree(tree.cal, nbOpt)
AssMotifs <- affect.motifs(affectElt, MOCCUR)
AssNoms   <- names.assembly(affectElt, MOCCUR)$motifs

boxplot(BETA ~ nchar(AssNoms), main = "BETA", ylim = range(BETA))
abline(h = 1, col = "red")
moyen <- numeric(nbOpt)
for (i in seq_len(nbOpt)) moyen[i] <- mean(BETA[nchar(AssNoms) == i])
points(y = moyen, x = seq_len(nbOpt), pch = 0)

summary(aov(BETA ~ AssNoms))        # p < 10-16
summary(aov(BETA ~ nchar(AssNoms))) # p = 0.77
table(nchar(AssNoms))
number.motifs(affectElt, MOCCUR)
summary(aov(BETA ~ XPR))        # p = 0.947


affectElt[] <- gf
AssMotifs   <- affect.motifs(affectElt, MOCCUR)
AssNoms     <- names.assembly(affectElt, MOCCUR)$motifs

boxplot(BETA ~ size, main = "BETA", ylim = range(BETA))
abline(h = 1, col = "red")
moyen <- numeric(nbOpt)
for (i in seq_len(nbOpt)) moyen[i] <- mean(BETA[size == names(table(size))[i]])
points(y = moyen, x = seq_len(nbOpt), pch = 0)

summary(aov(BETA ~ AssNoms))        # p < 10-16
summary(aov(BETA ~ size))           # p = 0.042
table(size)

nbElt  <- dim(MOCCUR)[2]
fctElt <- motElt <- NULL
for (elt in seq_len(nbElt)) {
  element <- elt
  indElt  <- which(mOccur[ , element] == 1)
  fctElt  <- c(fctElt, BETA[indElt])
  motElt  <- c(motElt, rep(elt, length(indElt)))
}
summary(aov(fctElt ~ motElt)) # p = 0.0585


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
titre    <- c("all.ALPHA")
nbOpt     <- 3
filename <- paste(rad, nom, titre, sep = "/")
tree.cal <- read.tree(filename)
affectElt <- cut.tree(tree.cal, nbOpt)
AssMotifs <- affect.motifs(affectElt, MOCCUR)
AssNoms   <- names.assembly(affectElt, MOCCUR)$motifs

boxplot(ALPHA ~ nchar(AssNoms), main = "ALPHA", ylim = range(ALPHA))
abline(h = 1, col = "red")
moyen <- numeric(nbOpt)
for (i in seq_len(nbOpt)) moyen[i] <- mean(ALPHA[nchar(AssNoms) == i])
points(y = moyen, x = seq_len(nbOpt), pch = 0)

summary(aov(ALPHA ~ nchar(AssNoms))) # p < 10-16  a, b , c
test.posthoc(ALPHA, nchar(AssNoms), pvalue = 0.001)
summary(aov(ALPHA ~ AssNoms))        # p < 10-16
table(nchar(AssNoms))
summary(aov(ALPHA ~ XPR))        # p = 0.021

boxplot(ALPHA ~ size, main = "ALPHA", ylim = range(ALPHA))
abline(h = 1, col = "red")
moyen <- numeric(length(table(size)))
for (i in seq_len(length(table(size)))) moyen[i] <- mean(ALPHA[size == names(table(size))[i]])
points(y = moyen, x = seq_len(length(table(size))), pch = 0)

summary(aov(ALPHA ~ size)) # p < 10-16
test.posthoc(ALPHA, nchar(AssNoms), pvalue = 0.001)

#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
close.pdf(enregistre)

