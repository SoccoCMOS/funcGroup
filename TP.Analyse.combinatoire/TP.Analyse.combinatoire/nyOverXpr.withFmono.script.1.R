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
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

nom <- setData[dat]

nomfic <- paste(rad, nom, sep = "/")
if (!file.exists(nomfic)) dir.create(nomfic)

filename <- paste(rad, paste(nom, "csv", sep = "."), sep = "/")
data     <- read.table(filename, header = TRUE, sep = ",")

setXpr  <- colnames(data)[(1 + nbElt + indXpr[1]):
                            (1 + nbElt + indXpr[length(indXpr)])]


ipr <- indXpr[5]
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


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  END of FILE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

