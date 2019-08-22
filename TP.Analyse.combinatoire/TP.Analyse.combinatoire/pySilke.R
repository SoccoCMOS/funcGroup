#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#     Combinatorial analysis
#     of an experiment of dynamic microorganism growth
#     on different combinations of 3 substrates
#
#     Data from Silke LANGENHEDER et al., 2010
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list = ls())

setwd("C:/Users/Jaillard/Mes Documents/1.BEF.Relation/0.R.code/Y1")
setwd("C:/Users/Jaillard/Mes Documents/Y1")

library(clusterCrit)

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

rad       <- c("Silke")

setData   <- c("glucose",
               "xylose",
               "galactose",
               "glucose.xylose",
               "glucose.galactose",
               "xylose.galactose",
               "glucose.xylose.galactose")

nbElt     <- 6
nbXpr     <- 6
indXpr    <- c(1:nbXpr)

gf        <- seq(1, nbElt)

enregistre <- TRUE
#enregistre <- FALSE


source("nyOverXpr.withFmono.script.R")


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  END of FILE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
