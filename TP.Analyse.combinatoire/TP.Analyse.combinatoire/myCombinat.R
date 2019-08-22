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



library(plyr)


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  List of functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# compact.index         <- function(v)

# affect.motifs         <- function(affectElt, mOccur)
# mAffect.elements      <- function(Rtree)
# mAffect.motifs        <- function(tree, mOccur)

# names.assembly        <- function(affectElt, mOccur)
# number.motifs         <- function(nbcl)

# delstr.begin          <- function(v, n)
# delstr.end            <- function(v, n)
# concat.byline         <- function(v, nbchar = 70)
# plot.bypage           <- function(tab, nbline = 25)
# plot.clusters.content <- function(affectElt, mOccur)
# plot.motifs.content   <- function(fct, AssMotifs, noms.motifs

# offset.combinations   <- function(nbcl)
# table.combinations    <- function(nbcl)
# table.size.motifs     <- function(affectElt, mOccur)

# label.elements        <- function(partition)
# set.partition         <- function(partition)
# generate.moccur       <- function(partition)

# hypergeo              <- function(partition, sampling)
# rsize                 <- function(partition, somme, level, v, res)
# fsize                 <- function(partition, somme)
# gsize                 <- function(partition, somme)
# table.size            <- function(partition)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return a vector continuously indexed from 1 to max(index)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compact.index <- function(v) {

  res <- v
  set <- sort(unique(v))
  for (i in seq_along(set)) res[v == set[i]] <- i

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Cluster the assemblages by assembly motifs
#
#  Inputs  : affectElt and mat Occurrence
#  Outputs : the vector of motifs to which belong the assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

affect.motifs <- function(affectElt, mOccur) {

  nbAss <- dim(mOccur)[1]

  EltClasses <- t( t(mOccur) * compact.index(affectElt) )
  rownames(EltClasses) <- seq_len(nbAss)

  myfct1  <- function(x) { x[x != 0] }
  tmp1    <- alply(EltClasses[ , , drop = FALSE], .margins = 1, .fun = myfct1)

  myfct2  <- function(x) { sort(unique(x)) }
  tmp2    <- llply(tmp1, .fun = myfct2)

  size    <- unlist(llply(tmp2, .fun = length))
  setSize <- sort(unique(size))

  motif     <- 1
  AssMotifs <- integer(nbAss)
  for (siz in seq_along(setSize)) {

    index <- which(size == setSize[siz])

    while (length(index) > 0) {
      ind <- 1
      ref <- tmp2[[index[ind]]]
      for (i in 2:length(index))
        if (identical(ref, tmp2[[index[i]]])) ind <- c(ind, i)

      AssMotifs[index[ind]] <- motif

      motif <- motif + 1
      index <- index[-ind]
    }
  }

  return(AssMotifs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Cluster the assemblages by assembly motifs
#
#  Inputs  : a tree and mat Occurrence
#  Outputs : the matrix of motifs to which belong the assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mAffect.motifs <- function(tree, mOccur) {

  nbclTot <- length(tree$cor)
  nbclMax <- which(tree$cor == max(tree$cor))[1]

  mAssMotifs <- matrix(as.integer(0), nrow = nbclTot, ncol = dim(mOccur)[1])
  for (nbcl in seq_len(nbclMax)) {
    affectElt          <- cut.tree(tree, nbcl)
    mAssMotifs[nbcl, ] <- affect.motifs(affectElt, mOccur)
  }

  if (nbclMax < nbclTot)
    for (nbcl in (nbclMax + 1):nbclTot)
      mAssMotifs[nbcl, ] <- mAssMotifs[nbclMax, ]

  colnames(mAssMotifs) <- rownames(mOccur)
  rownames(mAssMotifs) <- c(1:nbclTot)

  return(mAssMotifs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# FONCTION alternative, OK pour affectMot non fourni (on ne s'en sert pas...)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mAffect.elements <- function(Rtree) {

  nbclMax    <- length(Rtree$order)
  mAffectElt <- matrix(as.integer(0), nrow = nbclMax, ncol = nbclMax)
  for (nbcl in seq_len(nbclMax)) mAffectElt[nbcl, ] <- cutree(Rtree, nbcl)

  colnames(mAffectElt) <- Rtree$labels
  rownames(mAffectElt) <- c(1:nbclMax)

  return(mAffectElt)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the names (in minus letters) of elements, clusters of elements,
#              and all possible assembly motifs associated to element clusters
#
#  Input:  mOccur    : matrix of occurrence
#          affectElt : affectation of species into classes
#  Output: names.elements         : names of elements
#          names.element.clusters : names of clusters of elements
#          names.assembly.motifs  : names of elemental assembly motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

names.assembly <- function(affectElt, mOccur) {

  nbAss <- dim(mOccur)[1]
  nbElt <- dim(mOccur)[2]

  # label the element clusters
  affectElt <- shift.affectElt(affectElt)

  nomAffect <- myLetters[affectElt]

  # label the assembly motifs
  EltClasses <- t( t(mOccur) * affectElt )
  rownames(EltClasses) <- seq_len(nbAss)

  myfct1  <- function(x) { x[x != 0] }
  tmp1    <- alply(EltClasses[ , , drop = FALSE], .margins = 1, .fun = myfct1)

  myfct2  <- function(x) { sort(unique(x)) }
  tmp2    <- llply(tmp1, .fun = myfct2)

  nomMotif <- character(nbAss)
  for (ass in seq_len(nbAss))
    for (elt in seq_along(tmp2[[ass]]))
      nomMotif[ass] <- paste(nomMotif[ass],
                             myLetters[tmp2[[ass]][elt]], sep = "")

  noms        <- list(colnames(mOccur), nomAffect, rownames(mOccur), nomMotif)
  names(noms) <- c("elements", "clusters", "assemblies", "motifs")

  return(noms)
}


#system.time(for (i in 1:100000) affect_motif.name(x)    )

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the number of possible motifs
#  Input:  affectElt and mOccur
#  Output: a vector of number of motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

number.motifs <- function(affectElt, mOccur) {

  partition <- table(affectElt)

  noms      <- names.assembly(affectElt, mOccur)
  nomMotif  <- unique(noms$motifs)
  nbMots    <- length(nomMotif)

  nbMotif   <- integer(nbMots)
  nbMotif[] <- 1
  names(nbMotif) <- nomMotif

  for (mot in 1:nbMots) {

    index <- NULL
    for (i in 1:nchar(nomMotif[mot]))
      index <- c(index, which(myLetters == substr(nomMotif[mot], i, i)))

    for (i in seq_along(index))
      nbMotif[mot] <- nbMotif[mot] *
                        sum(choose(partition[index[i]], 1:partition[index[i]]))
  }

  return(nbMotif)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

offset.combinations <- function(nbcl)
{
  offset <- integer(nbcl)
  off    <- choose(nbcl, 1:nbcl)
  for (i in 2:nbcl) offset[i] <- sum(off[1:(i - 1)])

  return(offset)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

table.combinations <- function(nbcl)
{
  tabComb <- list()
  for (i in 1:nbcl)
    tabComb <- c(tabComb, combinat::combn(nbcl, i, simplify = FALSE))

  return(tabComb)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

table.size.motifs <- function(affectElt, mOccur) {

  nbElt         <- length(affectElt)

  AssMotifs     <- affect.motifs(affectElt, mOccur)
  nbMots        <- length(unique(AssMotifs))

  mat           <- matrix(0, nrow = nbElt, ncol = nbMots)
  colnames(mat) <- c(1:nbMots)
  rownames(mat) <- c(1:nbElt)

  size          <- apply(mOccur, MARGIN = 1, FUN = sum)
  for (i in 1:nbElt) {
    tmp <- table(AssMotifs[size == i])
    mat[i, names(tmp) ] <- tmp
  }

  noms          <- names.assembly(affectElt, mOccur)
  colnames(mat) <- unique(noms$motifs)

  return(mat)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# name each element by class
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

label.elements <- function(partition) {

  v1 <- v2 <- NULL
  for (i in seq_along(partition)) {
    v1 <- c(v1, rep(letters[i], partition[i]))
    v2 <- c(v2, seq(1:partition[i]))
  }

  res <- paste0(v1, v2)

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# name and index each element belonging to a partition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

set.partition <- function(partition) {

  affectElt <- NULL
  for (i in seq_along(partition))
    affectElt <- c(affectElt, rep(i, partition[i]))

  names(affectElt) <- label.elements(partition)

  return(affectElt)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# generate an occurring matrix (mOccur)
#        belonging all the possible combinations assembled with affectElt
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

generate.moccur <- function(partition) {

  nbElt  <- sum(partition)

  mOccur <- NULL
  for (j in 1:nbElt) {

    nbComb <- choose(nbElt, j)
    b      <- t(combn(nbElt, j))
    a      <- matrix(0, nrow = nbComb, ncol = nbElt)
    for (i in seq_len(nbComb)) a[i, b[i, 1:j]] <- 1
    mOccur <- rbind(mOccur, a)
  }

  colnames(mOccur) <- label.elements(partition)

  return(mOccur)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# sample mOccur from a full-mOccur matrix
#   for  nbAss assemblages,
#   with nbMot different motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# first version
sample.moccur <- function(nbAss, affectElt, full.mOccur, opt.size = NULL) {

  index <- seq_len(dim(full.mOccur)[1])


  if (length(opt.size) > 0) {
    nbMot <- length(opt.size)
    size  <- apply(full.mOccur, MARGIN = 1, FUN = sum)
    index <- which(size %in% opt.size)
  }

  tmp.mOccur <- full.mOccur[index, ]


  tmp.nbAss <- dim(tmp.mOccur)[1]
  if (nbAss > tmp.nbAss) nbAss <- tmp.nbAss

  full.nbMot <- 2^length(unique(affectElt)) - 1
  if (nbMot > full.nbMot) nbMot <- full.nbMot


  AssMotifs <- 1
  while (length(unique(AssMotifs)) < nbMot) {
    index     <- sample(seq_len(tmp.nbAss), size = nbAss, replace = FALSE)
    mOccur    <- tmp.mOccur[index, ]
    AssMotifs <- affect.motifs(affectElt, mOccur)
  }

  return(mOccur)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sample.moccur.bysize <- function(nbAss, affectElt, full.mOccur,
                                 opt.size = NULL) {

  # select a subset of full.mOccur

  nbMot     <- 2
  index     <- index2 <- seq_len(dim(full.mOccur)[1])
  full.size <- aaply(full.mOccur, .margins = 1, .fun = sum)

  if (length(opt.size) > 0) index <- which(full.size %in% opt.size)
  if (length(index) == 0) stop("no assemblage corresponds to your request")

  tmp.mOccur <- full.mOccur[index, ]
  tmp.size   <- full.size[index]

  # checking the parameters

  tmp.nbAss <- dim(tmp.mOccur)[1]
  if (nbAss > tmp.nbAss) nbAss <- tmp.nbAss

  # loocking for an adequat sampling

  motif <- 1
  while (length(unique(motif)) < nbMot) {
    index  <- sample(seq_len(tmp.nbAss), size = nbAss, replace = FALSE)
    mOccur <- tmp.mOccur[index, ]
    motif  <- affect.motifs(affectElt, mOccur)
  }

  return(mOccur)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sample.moccur.bysize <- function(nbAss, AssMotifs, full.mOccur,
                                 opt.size = NULL) {

  # select a subset of full.mOccur

  nbMot     <- 2
  index     <- seq_len(dim(full.mOccur)[1])
  full.size <- aaply(full.mOccur, .margins = 1, .fun = sum)

  if (length(opt.size) > 0) index <- which(full.size %in% opt.size)
  if (length(index) == 0) stop("no assemblage corresponds to your request")

  tmp.mOccur <- full.mOccur[index, ]
  tmp.size   <- full.size[index]
  tmp.motif  <- AssMotifs[index]

  # checking the parameters

  tmp.nbAss <- dim(tmp.mOccur)[1]
  if (nbAss > tmp.nbAss) nbAss <- tmp.nbAss

  # loocking for an adequat sampling

  motif <- 1
  while (length(unique(motif)) < nbMot) {
    index  <- sample(seq_len(tmp.nbAss), size = nbAss, replace = FALSE)
    mOccur <- tmp.mOccur[index, ]
    motif  <- tmp.motif[index]
  }

  return(index)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sample.moccur.eachmotif <- function(nbAss, affectElt, full.mOccur) {

  # select a subset of full.mOccur

  AssMotifs <- affect.motifs(affectElt, full.mOccur)
  setMots   <- unique(AssMotifs)
  nbMots    <- length(setMots)

  index <- NULL
  for (mot in seq_len(nbMots)) {
    motif  <- setMots[mot]
    ind    <- which(AssMotifs == motif)
    if (nbAss > length(ind)) {
      snbAss <- length(ind)
    } else {
      snbAss <- nbAss
    }
    indtmp <- sample(x = ind, size = snbAss, replace = FALSE)
    index  <- c(index, indtmp)
  }
  mOccur <- full.mOccur[index, ]

  return(mOccur)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sample.moccur.eachmotif <- function(nbAss, AssMotifs, full.mOccur) {

  # select a subset of full.mOccur

  setMots <- unique(AssMotifs)
  nbMots  <- length(setMots)

  index <- NULL
  for (mot in seq_len(nbMots)) {
    motif  <- setMots[mot]
    ind    <- which(AssMotifs == motif)
    if (nbAss > length(ind)) { snbAss <- length(ind)
    } else {                   snbAss <- nbAss
    }
    indtmp <- sample(x = ind, size = snbAss, replace = FALSE)
    index  <- c(index, indtmp)
  }

  return(index)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# generate a random vector of performance (Fobs) of same size
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

generate.fobs <- function(AssMotifs, amplitude, cv,
                          opt.law = "norm", dev = borne, coef = correct) {

  nbAss   <- length(AssMotifs)
  setMots <- unique(AssMotifs)
  nbMots  <- length(setMots)

  if (opt.law == "unif")
    noise <- sqrt(3) * runif(nbAss, min = -1, max = +1)
  if (opt.law == "norm")
    noise <- coef    * rnorm.trunc(nbAss, min = -dev, max = +dev)

  perform <- runif(nbMots, min = 1, max = amplitude)

  Fobs <- numeric(nbAss)
  for (mot in 1:nbMots) {
    motif       <- setMots[mot]
    index       <- which(AssMotifs == motif)
    Fobs[index] <- perform[motif] * (1 + cv * noise[index])
  }

  return(Fobs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Hyper-geometrical law
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

hypergeo <- function(partition, sampling) {

  res <- 1
  for (i in seq_along(partition))
    res <- prod(res, choose(partition[i], sampling[i]))

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Recursive determining the different no-empty combinations
#          of size = somme from a set = partition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rsize <- function(partition, somme, level, v, res) {

  for (index in 1:partition[level]) {

    if ( (level < length(partition)) & (sum(v) < somme) ) {
      v[level] <- index
      res      <- rsize(partition, somme, level + 1, v, res)

    } else {
      v[level] <- somme - sum(v)
      if (v[level] > 0) res <- rbind(res, v)

      break
    }
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Function calling the recursive
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

fsize <- function(partition, somme) {

  res <- 0
  if (somme >= length(partition)) {
    res <- NULL
    res <- rsize(partition, somme,
                 level = 1, v = rep(0, length(partition)), res)
    rownames(res) <- seq_len(dim(res)[1])
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  recursive determining the different sampling to obtain a given sum
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gsize <- function(partition, somme) {

  tmp1 <- fsize(partition, somme)

  tmp2 <- 0
  if (is.matrix(tmp1))
    for (i in seq_len(dim(tmp1)[1]))
      tmp2 <- sum(tmp2, hypergeo(partition, tmp1[i, ]) )

  return(tmp2)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  recursive determining the different sampling to obtain a given sum
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

table.size <- function(partition) {

  nbElt         <- sum(partition)
  nbMot         <- 2^length(partition) - 1

  mat           <- matrix(as.integer(0), nrow = nbElt, ncol = nbMot)
  colnames(mat) <- seq_len(nbMot)
  rownames(mat) <- seq_len(nbElt)

  iter <- 0
  v    <- integer(nbElt)
  nom  <- character(nbMot)
  for (size in 1:length(partition)) {
    tmp1 <- as.matrix(combinat::combn(myLetters[seq_along(partition)], size))
    tmp2 <- as.matrix(combinat::combn(partition, size))
    for (mot in seq_len(dim(tmp1)[2])) {

      iter <- iter + 1
      tmp  <- NULL
      for (i in 1:size) tmp <- paste(tmp, tmp1[i, mot], sep = "")
      nom[iter] <- tmp

      for (i in 1:nbElt) v[i] <- gsize(tmp2[ , mot], i)
      mat[ , iter] <- v
    }
  }

  colnames(mat) <- nom

  return(mat)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                          END of FILE myCOMBINAT.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



