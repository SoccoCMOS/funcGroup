#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                 myTREE.R
#
#     set of functions for manipulating hierarchical trees
#
#                         Benoit JAILLARD, summer 2016
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


source("myStats.R")

library(plyr)
library(partitions)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   List of functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# simplify.tree            <- function(tree, vBool)
# sort.tree                <- function(tree, index.return = FALSE)
# plot.tree                <- function(tree,
#                                      signifElt = tree$aff[dim(tree$aff)[2],],
#                                      col = "black", titre = "")
# rm.tree.level            <- function(tree, level)
# rm.stats.level           <- function(stats, level)
# check.tree               <- function(tree, stats)
# permute                  <- function(mat, x, y)
# reverse                  <- function(v)
# next.partitions          <- function(partition)
# is.identical.partition   <- function(partition, mPartitions)
# is.connected.partition   <- function(parent, child)
# is.tree.partition        <- function(mPartitions)
# cut.tree                 <- function(res, nbcl)
# as.mytree                <- function(Rtree)
# as.Rtree                 <- function(mytree)
# write.tree               <- function(filename, tree)
# read.tree                <- function(filename)


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Simplify a tree by keeping only the "TRUE" elements
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

simplify.tree <- function(tree, vBool) {

  tree$aff <- tree$aff[ , vBool]

  fct   <- function(x) { length(unique(x)) }
  laff  <- apply(tree$aff, MARGIN = 1, FUN = fct)
  ind1  <- which(table(laff) > 1)

  index <- NULL
  for (elt in seq_along(ind1)) index <- c(index, which(laff == ind1[elt])[-1])

  if (length(index) > 0) {

    tree$cor <- tree$cor[-index]
    tree$aff <- tree$aff[-index, ]

    for (elt in 2:dim(tree$aff)[2]) {

      l1 <- tree$aff[elt - 1, ]
      l2 <- tree$aff[elt, ]

      for (clu in 1:elt) {

        if (length(which(l2 == clu)) == 0) {

          index <- which(l1 == clu)
          for (i in seq_along(index)) {
            elt2  <- elt
            while (   (elt2 <= dim(tree$aff)[2])
                   && (tree$aff[elt2, index[i]] == l2[index[i]]) ) {
              tree$aff[elt2, index[i]] <- clu
              elt2 <- elt2 + 1
            }
          }

        }
        tree$aff[elt, ] <- compact.index(tree$aff[elt, ])
      }
    }

    rownames(tree$aff) <- seq_len(sum(vBool != 0))
  }

  return(tree)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Sort the resulting file of tree
#
#  Input:  X : matrix of leaves : affectation of species to classes
#          Y : vector of distance
#          Z : matrix of distance from the best element
#  Output: files of leaves sorted for plotting
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort.tree <- function(tree, index.return = FALSE) {

  X      <- shift.affectElt(tree$aff)

  nbline <- dim(X)[1]
  nbitem <- dim(X)[2]


  # Separate first the sorted from unsorted elements
  mask <- ordre <- seq_len(nbitem)
  for (lin in 2:nbline) {
    index       <- sort(x = X[lin, mask], index.return = TRUE)
    index       <- index$ix + nbitem - length(mask)
    X[ , mask]  <- X[ , index]
    ordre[mask] <- ordre[index]
    mask        <- mask[X[lin, mask] == max(X)]
  }

  # Sort the leaves of tree
  new <- old <- seq_len(nbitem)
  for (lin in seq_len(nbline)) {
    for (item in seq_len(nbitem))
      old[item] <- which(unique(X[lin,]) == X[lin, item])
    new             <- sort(old, index.return = TRUE)
    X[ , seq_len(nbitem)]  <- X[ , new$ix]
    ordre[seq_len(nbitem)] <- ordre[new$ix]
  }

  # Re-numeration of leaves
  old <- integer(nbline)
  for (lin in seq_len(nbline)) {
    v  <- setdiff(unique(X[lin, ]), unique(old))
    if (length(v) == 1) old[lin] <- v
  }

  newX <- matrix(0, nrow = nbline, ncol = nbitem)
  for (lin in seq_len(nbline))
    for (j in seq_len(lin)) newX[lin, which(old[j] == X[lin, ])] <- j

  colnames(newX) <- colnames(X)[ordre]

  if (index.return == TRUE) {
    res        <- list(newX, ordre, colnames(newX))
    names(res) <- c("x", "ix", "noms")
  } else {
    res <- newX
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#          Plot a Divisive or Agglomerative Hierarchical Clustering
#    obtained with the options ["divisive" and "tree"] or "agglomerative"
#
#  Inputs: X: matrix of classes of elements
#          Y: vector of distance
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot.tree <- function(tree,
                      signifElt = NULL,
                      col       = "black", titre = "") {

  if (length(signifElt) == 0) signifElt <- rep(TRUE, length(tree$cor))

  vpch            <- numeric(dim(tree$aff)[2])
  vpch[ ]         <- 1
  vpch[signifElt] <- 19

  if (length(tree$cor) > 1) {
    tmpx <- sort.tree(tree, index.return = TRUE)
    X    <- tmpx$x
    vpch <- vpch[tmpx$ix]
    if (length(col) > 1) col <- col[tmpx$ix]
  } else {
    X    <- shift.affectElt(tree$aff)
    vpch <- vpch[1]
    if (length(col) > 1) col <- col[1]
  }

  Y   <- tree$cor
  Y[Y < 0] <- 0


  nbline <- dim(X)[1]
  nbitem <- dim(X)[2] + 1
  YY     <- c(0,Y)

  xx     <- matrix(0, nrow = nbline, ncol = nbline)
  xy     <- NULL

  # plot the framework
  plot(x = seq_len(nbitem), y = seq(1/nbitem, 1, 1/nbitem),
       xlab = "element",    ylab = "R2-value",
       xlim = c(1, nbitem), ylim = c(0, 1),
       type = "n", tck = 0.02, las = 1)

  # plot the vertical lines
  for (lin in 1:nbline)
    for (leave in 1:lin) {
      xx[lin, leave] <- 1 + mean(which(X[lin, ] == leave))
      lines(x = rep(xx[lin, leave], 2), y = c(YY[lin], YY[lin + 1]),
            lty = "solid")
    }

  # plot the horizontal jonction lines
  if (nbline > 1) for (lin in 2:nbline) {
    XX <- which(xx[lin, ] != xx[lin - 1, ])
    if (length(XX) == 1) XX <- rep(XX, 2)
    lines(x = xx[lin, XX], y = rep(YY[lin], 2), lty = "solid" )
  }

  # plot the final horizontal jonction lines
  setLeave <- sort(unique(X[nbline, ]))
  for (leave in seq_along(setLeave)) {
    tmp <- which(X[nbline, ] == setLeave[leave])
    lines(x = c(min(tmp), max(tmp)) + 1,
          y = rep(YY[nbline + 1], 2),
          lty = "solid" )
  }

  # plot the symbols of the most likely partition
  xpos <- 1 + c(1:(nbitem - 1))
  ypos <- YY[length(YY)] + 0.02

  if (length(col) > 1) {
    points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch,
           cex = 2, col = col, bg = "white")
  } else {
    points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch, bg = "white")
  }

  #  plot the names of elements
  axis(side = 3, at = 2:nbitem, labels = colnames(X),
       tick = FALSE, las = 2, pos = ypos)

  title(titre)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Delete a level in a clustering tree
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm.tree.level <- function(tree, level) {

  if ((level < 2) | (level > dim(tree$aff)[1]))
    stop(" 'level' should be comprized between 2 and dim(tree)[1] ")

  levMax <- dim(tree$aff)[1]

  tree$aff[tree$aff >= (level + 1)] <- tree$aff[tree$aff >= (level + 1)] - 1

  tree$aff[level:(levMax - 1), ] <- tree$aff[(level + 1):levMax, ]

  index <- which(table(tree$aff[levMax, ]) > 1)
  index <- which(tree$aff[levMax, ] == index)[2]
  tree$aff[levMax, index] <- levMax

  tree$cor[level:(levMax - 2)] <- tree$cor[(level + 1):(levMax - 1)]

  return(tree)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm.stats.level <- function(stats, level) { return(stats[-level, ]) }




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Check a clustering tree by deleting all unuseful levels
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

check.tree <- function(tree, stats) {

  tree$cor <- stats[ ,"R2cal"]

  while
    (first.optimum(tree$cor, "max") != which(tree$cor == max(tree$cor))[1]) {

    level <- first.optimum(tree$cor, "max") + 1
    tree  <- rm.tree.level(tree, level)
    stats <- rm.stats.level(stats, level)
    }

  res <- list(tree, stats)
  names(res) <- c("tree", "stats")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# permute two values into a matrix
# return the right matrix
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

permute <- function(mat, x, y) {

  mask.x <- which(mat == x)
  mask.y <- which(mat == y)

  mat[mask.x] <- y
  mat[mask.y] <- x

  return(mat)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# permute two values into a matrix
# return the right matrix
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

reverse <- function(v) {

  index <- unique(v)
  mask1 <- which(v == index[1])
  mask2 <- which(v == index[2])

  v[mask1] <- index[2]
  v[mask2] <- index[1]

  return(v)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# return all the possible next connected partitions
#                              issued from a given partition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

next.partitions <- function(partition) {

  sizeparts <- table(partition)
  lastclass <- length(names(sizeparts))

  sizeparts <- sizeparts[sizeparts > 1]
  nomclass  <- as.integer(names(sizeparts))

  newpartition <- NULL
  for (i in 1:length(sizeparts)){

    p     <- t(partitions::parts(sizeparts[i]))
    index <- which(apply((p != 0), MARGIN = 1, FUN = sum) == 2)
    p     <- p[index, c(1:2), drop = FALSE]

    mask  <- which(partition == nomclass[i])
    for (j in 1:dim(p)[1]) {

      res <- t(partitions::setparts(p[j, ]))
      res[res == 2] <- lastclass + 1
      res[res == 1] <- nomclass[i]

      newpart         <- matrix(partition, nrow = dim(res)[1],
                                ncol = length(partition), byrow = TRUE)
      newpart[, mask] <- res
      newpartition    <- rbind(newpartition, newpart)
    }
  }

  res <- newpartition
  colnames(res) <- names(partition)

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# look for two identical partitions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

is.identical.part <- function(part1, part2) {

  crit  <- sort(table(part1), decreasing = TRUE)
  nbclu <- length(crit)
  perm  <- matrix(unlist(combinat::permn(as.integer(names(crit)))),
                  ncol = nbclu, byrow = TRUE)

  mask <- list(logical)
  for (i in 1:nbclu) mask[[i]] <- which(part1 == i)

  res  <- FALSE
  for (per in 1:dim(perm)[1]) {

    part <- part1
    for (i in 1:nbclu) part[mask[[i]]] <- perm[per, i]

    mtmp <- rbind(part2, part)
    tmp  <- apply(mtmp[ , , drop = FALSE], MARGIN = 1,
                  FUN = identical, y = mtmp[dim(mtmp)[1], ])
    res  <- (length(setdiff(which(tmp == TRUE), dim(mtmp)[1])) > 0)
    if (res == TRUE) break
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# look for identical partitions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

is.identical.partition.other <- function(partition, mPartitions) {

  # coerce "mPartitions" to be a matrix
  if (!is.matrix(mPartitions))
    mPartitions <- matrix(mPartitions, nrow = 1, ncol = length(partition))

  # keep only the same type (cardinal) of partitions
  res   <- FALSE
  crit1 <- sort(table(partition), decreasing = TRUE)

  tmp   <- t(apply(mPartitions[ , , drop = FALSE], MARGIN = 1, FUN = table))
  crit2 <- t(apply(tmp[ , , drop = FALSE], MARGIN = 1,
                   FUN = sort, decreasing = TRUE))
  lcrit2 <- unique(apply(crit2[,, drop = FALSE], MARGIN = 1, FUN = length))

  if (length(crit1) != lcrit2) {
    return(res)
  } else {
    tmp2 <- rbind(crit2, crit1)
  }

  index <- apply(tmp2[ , , drop = FALSE], MARGIN = 1,
                 FUN = identical, y = tmp2[dim(tmp2)[1], ])
  mPart <- mPartitions[setdiff(which(index == TRUE), dim(tmp2)[1]),
                       , drop = FALSE ]

  # look for possible permutations of set notation
  nbclu <- length(crit1)
  perm  <- matrix(unlist(combinat::permn(as.integer(names(crit1)))),
                  ncol = nbclu, byrow = TRUE)

  mask <- list(logical)
  for (i in 1:nbclu) mask[[i]] <- which(partition == i)

  for (per in 1:dim(perm)[1]) {

    part <- partition
    for (i in 1:nbclu) part[mask[[i]]] <- perm[per, i]

    mtmp <- rbind(mPart, part)
    tmp  <- apply(mtmp[ , , drop = FALSE], MARGIN = 1,
                  FUN = identical, y = mtmp[dim(mtmp)[1], ])
    res  <- (length(setdiff(which(tmp == TRUE), dim(mtmp)[1])) > 0)
    if (res == TRUE) break
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# look for if a given partition is possibly chid partition
#                                        connected to another parent partition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

is.connected.partition <- function(parent, child) {

  res <- FALSE

  crit1 <- table(parent)
  crit2 <- table(child)
  if (length(crit2) != length(crit1) + 1) return(res)

  children <- next.partitions(parent)
  res      <- is.identical.partition(child, children)

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# look for if a given set of partitions is possibly a tree partition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

is.tree.partition <- function(mPartitions) {

  tmp    <- apply(mPartitions[ , , drop = FALSE], MARGIN = 1, FUN = table)
  crit2  <- lapply(tmp, FUN = sort, decreasing = TRUE)
  lcrit2 <- unlist(unique(lapply(crit2, FUN = length)))
  index  <- sort(lcrit2, index.return = TRUE)
  mPartitions <- mPartitions[index$ix, ]

  res <- FALSE
  if (!identical(lcrit2, c(lcrit2[1]:(lcrit2[1] + length(lcrit2) - 1))))
    return(res)

  nbClu <- dim(mPartitions)[1]
  res   <- logical(nbClu - 1)
  for (clu in 1:(nbClu - 1))
    res[clu] <- is.connected.partition(mPartitions[clu, ],
                                       mPartitions[clu + 1, ])

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Cut a tree
#    and output a vector of affectation as a raw vector
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cut.tree <- function(res, nbcl) { compact.index(res$aff[nbcl, ]) }



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Function to transform R-tree into matrix-tree
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

as.mytree <- function(Rtree) {

  height <- rev(c(0, Rtree$height))
  height <- 1 - height / max(height)
  nbLev  <- length(height)

  mytree <- matrix(0, nrow = nbLev, ncol = nbLev)
  rownames(mytree) <- colnames(mytree) <- seq_len(nbLev)

  for (nbcl in seq_len(nbLev)) mytree[nbcl, ] <- cutree(Rtree, nbcl)

  mask2 <- list()
  nbcl  <- 1
  for (i in seq_len(nbcl)) mask2[[i]] <- which(mytree[nbcl, ] == i)

  bool  <- matrix(FALSE, nrow = nbLev, ncol = nbLev)

  for (nbcl in 2:nbLev) {
    mask1 <- mask2
    for (i in seq_len(nbcl)) mask2[[i]] <- which(mytree[nbcl, ] == i)

    for (i in seq_len(nbcl - 1))
      for (j in seq_len(nbcl))
        bool[i, j] <- identical(mask1[[i]], mask2[[j]])

      index <- which(apply(bool, MARGIN = 2, FUN = sum) == 0)

      mytree[nbcl, ] <- mytree[nbcl - 1, ]
      mytree[nbcl, mask2[[ index[2] ]]] <- nbcl
  }

#  mytree <- mytree[ , Rtree$order]

  res <- list(mytree,  height)
  names(res) <- c("aff", "cor")

  return(res)
}

#plot(dot.arbre, hang = -1)
#mytree <- as.mytree(dot.arbre)
#plot.tree(mytree$aff, mytree$cor, "essai")



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Function to transform matrix-tree into R-tree  INACHEVEE !!!
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

as.Rtree <- function(mytree) {

  aff <- mytree$aff
  cor <- mytree$cor

  nbLev <- dim(aff)[1] - 1
  nbElt <- dim(aff)[2]

  tmp <- res$merge <- matrix(0, nrow = nbLev, ncol = 2)
  aff <- sort.tree(aff, index.return = TRUE)

  mask <- !logical(nbElt)
  for (lev in nbLev:1) {
    index <- which(aff$x[lev + 1, ] != aff$x[lev, ])
    leaf1 <- aff$x[lev + 1, index]
    leaf2 <- aff$x[lev,     index]
    line  <- range(leaf1, leaf2)

    if (mask[line[1]] == TRUE) {
      mask[line[1]] <- FALSE
      line[1]       <- - line[1]
    }

    if (mask[line[2]] == TRUE) {
      mask[line[2]] <- FALSE
      line[2]       <- - line[2]
    }

    if (!( (line[1] < 0) & (line[2] < 0))) {
      line <- range(line)
    } else {
      line <- rev(range(line))
    }

    tmp[lev, ] <- line
  }

  tmp <- tmp[rev(seq_len(nbLev)), ]

  v <- as.vector(t(tmp))
  i <- 0
  mask <- which(v > 0)
  w <- v[v > 0]
  index <- numeric(length(v))
  while (length(w) > 0) {
    i <- i + 1
    index[i] <- w[1]
    w <- w[w != index[i]]
  }
  index <- index[index > 0]

  w <- numeric(length(mask))
  for (i in 1:length(index)) w[v[mask] == index[i]] <- i
  v[mask] <- w

  tmp2 <- matrix(v, nrow = nbLev, ncol = 2, byrow = TRUE)

  res$merge <- tmp2

  res$height <- cor[1:nbLev]

  res$order  <- aff$ix

  res$labels <- colnames(mytree$aff)

  res$method <- c("ward.D")

  res$call   <- c("hclust(d = dot.dist, method = 'ward.D')")

  res$dist.method <- c("euclidean")

  res <- list(res$merge, res$height, res$order, res$labels,
              res$method, res$call, res$dist.method)
  names(res) <- c("merge", "height", "order", "labels",
                  "method", "call", "dist.method")

  attr(res, "class") <- c("hclust")

  return(res)
}

#res <- as.Rtree(res.Fobs$aff)
#plot(res)


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

write.tree <- function(filename, tree) {

  if ((names(tree)[1] == "aff") & (names(tree)[2] == "cor")) {

    tmp           <- cbind(tree$cor, tree$aff)
    colnames(tmp) <- c("cor", colnames(tree$aff))
    write.table(x = tmp, file = paste(filename, "tree", "csv", sep = "."),
                append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

  } else {

    write.table(x = tmp, file = paste(filename, "tree", "csv", sep = "."),
                append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")

  }

}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

read.tree <- function(filename) {

  tmp  <- read.table(file = paste(filename, "tree", "csv", sep = "."),
                     header = TRUE, sep = ",")
  .aff <- as.matrix(tmp[ , 2:dim(tmp)[2]])
  .cor <- as.vector(tmp[ , 1])

  res  <- list(.aff, .cor)
  names(res) <- c("aff", "cor")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                        END of file myTREE.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
