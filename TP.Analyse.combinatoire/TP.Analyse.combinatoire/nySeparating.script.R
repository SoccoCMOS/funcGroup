#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#     mySEPARATING.R = set of instructions to analyse a diversity dataset
#
#   SCRIPT for Separating the Interaction effect from the Composition effect
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Parameters for graph plotting
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
opt.size <- TRUE
opt.gf   <- FALSE
opt.log  <- FALSE
opt.aov  <- FALSE


higth <- max(alpha) - min(alpha)
width <- max(beta)  - min(beta)
posX1 <- min(beta)  - 0.02*width
posX2 <- min(beta)  + 1.02*width
posY1 <- min(alpha) - 0.02*higth
posY2 <- min(alpha) + 1.02*higth



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot Fmono
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if (length(unique(xpr)) > 1) {

  speFmono <- matrix(0, nrow = length(unique(xpr)), ncol = nbElt)
  index    <- c(0, (1:(length(unique(xpr)) - 1)) * nbElt)
  for (elt in 1:nbElt) speFmono[ , elt] <- Fmono[elt + index]

  mFmono <- apply(speFmono, MARGIN = 2, FUN = amean)
  sFmono <- apply(speFmono, MARGIN = 2, FUN = asd)
  index  <- sort(mFmono, decreasing = TRUE, index.return = TRUE)

  plot(x = seq_len(nbElt), xlab = "elements",
       y = mFmono[index$ix], ylab = "performance in monoculture",
       ylim = c(0, max(mFmono + sFmono, mFmono - sFmono)),
       type = "n", las = 1, cex = 2, tck = 0.02)

  for (elt in seq_along(index$ix))
    points.sd(x = rep(elt, length(unique(xpr))),
              y = speFmono[ , index$ix[elt]],
              figure = figures[1], couleur = "black",
              nom = "",
              opt.std = TRUE,
              scale.log = FALSE)

  axis(seq_len(nbElt), colnames(mOccur)[index$ix],
       side = 3, tick = FALSE, las = 2, col = "black")


  plot(x = seq_along(unique(xpr)), xlab = "time",
       y = Fscale, ylab = "performance in monoculture",
       ylim = c(75, 160),
       type = "b", las = 1, cex = 2, tck = 0.02)

} else {
  index <- sort(Fmono, decreasing = TRUE, index.return = TRUE)

  plot(x = seq_len(nbElt), xlab = "elements",
       y = Fmono[index$ix], ylab = "performance in monoculture",
       type = "p", las = 1, cex = 2, tck = 0.02)

  axis(seq_len(nbElt), colnames(mOccur)[index$ix],
       side = 3, tick = FALSE, las = 2, col = "black")
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot beta versus Fmono
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ymin <- min(betaElt - sdbetaElt)
ymax <- max(betaElt + sdbetaElt)
fct  <- Fmono

plot(x = fct,   xlab = "performance in monoculture",
     y = betaElt, ylab = "beta", ylim = c(ymin, ymax),
     type = "n", las = 1, cex = 2, tck = 0.02)

arrows(x0 = fct, y0 = betaElt - sdbetaElt,
       x1 = fct, y1 = betaElt + sdbetaElt,
       col = "black", length = 0.1, angle = 90, code = 3,
       lwd = 1, lty = "solid")

points(x = fct,  y = betaElt,
       type = "p", las = 1, cex = 2, bg = "white")

R2 <- cor2(betaElt, fct)

text(x = min(fct) + 3/4*(max(fct) - min(fct)), y = min(betaElt),
       labels = paste("R2 = ", signif(R2, digits = 3), sep = ""),
       pos = 3, col = "blue")


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot fobs
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "fobs"
ghisto(fobs, titre, scale.log = FALSE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot alpha
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "alpha"

ghisto(alpha, titre, scale.log = FALSE)
plot.alphaVSbeta(x = fobs, xlab = "fobs",
                 y = alpha, ylab = titre,
                 titre = titre, scale.log = FALSE)

plot.alphaVSbeta(y = fobs, ylab = "fobs",
                 x = alpha, xlab = titre,
                 titre = titre, scale.log = FALSE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot beta
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "beta"

ghisto(beta, titre, scale.log = FALSE)
plot.alphaVSbeta(x = fobs, xlab = "fobs",
                 y = beta, ylab = titre,
                 titre = titre, scale.log = FALSE)

plot.alphaVSbeta(y = fobs, ylab = "fobs",
                 x = beta, xlab = titre,
                 titre = titre, scale.log = FALSE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot log2(alpha) vs log2(beta)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot.alphaVSbeta(x = beta,  xlab = "beta",
                 y = alpha, ylab = "alpha",
                 titre = "alpha vs beta", scale.log = FALSE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot alpha vs beta WITH TRANSGRESSIVE ASSEMBLAGES in RED
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if (opt.aov == TRUE) {
  talpha <- test.posthoc(alpha, trans, pvalue)
  tbeta  <- test.posthoc(beta,  trans, pvalue)
}
setNoms <- c("trans-under", "under", "over", "trans-over")

plot(x = beta,  xlab = "beta",
     y = alpha, ylab = "alpha",
     main = "alpha vs beta: transgressive overyielding",
     type = "n", tck = 0.02, las = 1)

abline(v = mean(beta), h = mean(alpha), col = "blue")
text(x = posX1, y = posY1,labels = "beta", col = "black", font = 3)
text(x = posX2, y = posY2,labels = "alpha", col = "black", font = 3)

setTrans <- sort(unique(trans), decreasing = TRUE)
for (i in 1:length(setTrans)) {
  j <- setTrans[i]
  points(x = beta[(trans == j)], y = alpha[(trans == j)],
         type = "p", pch = figures[i], col = couleurs[i],
         bg = "white", cex = 2)

  points.sd(x = beta[(trans == j)],  y = alpha[(trans == j)],
            figure = figures[i], couleur = couleurs[i], nom = setNoms[j],
            opt.std = TRUE, scale.log = FALSE)

  if (opt.aov == TRUE) if (is.list(talpha)) {
    index <- which(talpha[ , "motifs"] == j)
    text(x = posX2,
         y = as.numeric(talpha[index, "means"]),
         labels = as.character(talpha[index, "groups"]),
         col =  couleurs[i], font = 3)
  }

  if (opt.aov == TRUE) if (is.list(tbeta)) {
    index <- which(tbeta[ , "motifs"] == j)
    text(x = as.numeric(tbeta[index, "means"]), y = posY1,
         labels = as.character(tbeta[index, "groups"]),
         col =  couleurs[i], font = 3)
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the mean performance of assemblages containing a given element
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
for (elt in seq_len(nbElt)) {

  index <- which(mOccur[ ,elt] == 1)

  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       main = paste("Element = ", colnames(mOccur)[elt], sep = ""),
       type = "n", tck = 0.02, las = 1)
  points(x = beta[-index], y = alpha[-index],
         type = "p", pch = figures[1], col = "black", bg = "white", cex = 2)
  points(x = beta[index], y = alpha[index],
         type = "p", pch = figures[2], col = "blue", bg = "white", cex = 2)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[2], couleur = "blue", nom = colnames(mOccur)[elt],
            opt.std = TRUE, scale.log = FALSE)
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITH RAW DATA
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = beta,  xlab = "beta",
     y = alpha, ylab = "alpha",
     main = "All elements",
     pch = figures[1], col = "black", bg = "white",
     type = "p", cex = 2, tck = 0.02, las = 1)

abline(h = 0, v = 0, col = "black")
abline(h = gmean(alpha), v = mean(beta), col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[ ,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = FALSE, scale.log = FALSE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITH RAW DATA WITH SD
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = beta,  xlab = "beta",
     y = alpha, ylab = "alpha",
     main = "All elements",
     pch = figures[1], col = "black", bg = "white",
     type = "p", cex = 2, tck = 0.02, las = 1)

abline(h = 0, v = 0, col = "black")
abline(h = gmean(alpha), v = mean(beta), col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = TRUE, scale.log = FALSE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITHOUT RAW DATA
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = beta,  xlab = "beta",
     y = alpha, ylab = "alpha",
     #     xlim = c(-0.2,0.4), ylim = c(0.4,0.7),
     main = "All elements",
     type = "n", bg = "white", tck = 0.02, las = 1)
abline(h = 0, v = 0, col = "black")
abline(h = gmean(alpha), v = gmean(beta), col = "blue")

R2.mean <- summary(lm(formula = alpha ~ beta))$"r.squared"
text(x = 0.8*max(beta), y = min(alpha),
     labels = paste("R2.mean = ", signif(R2.mean, digits = 3),sep = ""),
     pos = 3, col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = FALSE, scale.log = FALSE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all elements with standard deviation WITHOUT RAW DATA WITH SD
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = beta,  xlab = "beta",
     y = alpha, ylab = "alpha",
     #     xlim = c(-0.2,0.4), ylim = c(0.4,0.7),
     main = "All elements",
     type = "n", bg = "white", tck = 0.02, las = 1)
abline(h = 0, v = 0, col = "black")
abline(h = gmean(alpha), v = mean(beta), col = "blue")

R2.mean <- summary(lm(formula = alpha ~ beta))$"r.squared"
text(x = 0.5 * max(beta), y = min(alpha),
     labels = paste("R2.mean = ", signif(R2.mean, digits = 3), sep = ""),
     pos = 3, col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = TRUE, scale.log = FALSE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the hierarchical tree
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(Rtree, xlab = "hierarchical tree",
     main = "Agregation Euclidien/Ward", hang = -1, las = 1)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances CLUSTERED BY ASSEMBLY MOTIF
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
for (nbclasses in seq_len(nbclassesMax)) {
  # Plot one species
  affectElt <- mAffectElt[nbclasses, ]
  AssMotifs <- mAssMotifs[nbclasses, ]

  setMots   <- sort(unique(AssMotifs))
  nbMots    <- length(setMots)

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }
  noms <- names.assembly(affectElt, mOccur)
  plot.clusters.content(affectElt, mOccur)


  # only raw data without either mean neither standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    text(x = mean(beta[index]), y = gmean(alpha[index]),
         col = couleurs[motif],
         labels = noms$motifs[motif], pos = 3, font = 3)
  }


  # without standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = FALSE, scale.log = FALSE)
  }


  # with standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")
  text(x = posX1, y = posY1,labels = "beta", col = "black", font = 3)
  text(x = posX2, y = posY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = TRUE, scale.log = FALSE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = posX2, y = as.numeric(talpha[index2, "means"]),
           labels = as.character(talpha[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = as.numeric(tbeta[index2, "means"]), y = posY1,
           labels = as.character(tbeta[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }


  }


  # Standard deviation only without raw data
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")
  text(x = posX1, y = posY1,labels = "beta", col = "black", font = 3)
  text(x = posX2, y = posY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = TRUE, scale.log = FALSE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = posX2, y = as.numeric(talpha[index2, "means"]),
           labels = as.character(talpha[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = as.numeric(tbeta[index2, "means"]), y = posY1,
           labels = as.character(tbeta[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

  }


  # with Alpha and Beta by motif
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[mot], col = couleurs[mot], bg = "white",
           type = "p", cex = 2)

    for (elt in seq_len(nbElt)) {
      index2 <- which(mOccur[index,elt] == 1)
      titre  <- paste(noms$motifs[motif],
                      colnames(mOccur)[elt], sep = ".")
      points.sd(x = beta[index[index2]], y = alpha[index[index2]],
                figure = figures[elt], couleur = couleurs[mot],
                nom = titre,
                opt.std = FALSE, scale.log = FALSE)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances clustered by the NUMBER of SPECIES
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if (opt.size == TRUE) {

  mainTitre <- "Assembly size"

  size      <- apply(mOccur, MARGIN = 1, FUN = sum)
  AssMotifs <- numeric(length(size))
  for (i in seq_along(unique(size))) AssMotifs[(size == unique(size)[i])] <- i

  setMots <- sort(unique(AssMotifs))
  nbMots  <- length(setMots)
  noms    <- as.character(sort(unique(size)))

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # without standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif],
           type = "p", bg = "white", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms[motif],
              opt.std = FALSE, scale.log = FALSE)
  }


  # with standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = amean(beta), col = "blue")
  text(x = posX1, y = posY1, labels = "beta", col = "black", font = 3)
  text(x = posX2, y = posY2, labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms[motif],
              opt.std = TRUE, scale.log = FALSE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = posX2, y = as.numeric(talpha[index2, "means"]),
           labels = as.character(talpha[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = as.numeric(tbeta[index2, "means"]), y = posY1,
           labels = as.character(tbeta[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances clustered by A PRIORI FUNCTIONAL GROUPS
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if (opt.gf == TRUE) {

  mainTitre <- "Functional Groups"

  affectElt <- gf
  names(affectElt) <- colnames(mOccur)

  AssMotifs <- affect.motifs(affectElt, mOccur)
  setMots   <- sort(unique(AssMotifs))
  nbMots    <- length(setMots)

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }

  noms <- names.assembly(affectElt, mOccur)
  plot.clusters.content(affectElt, mOccur)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # without standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       main = mainTitre,
       type = "n", bg = "white", cex = 2, tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = FALSE, scale.log = FALSE)
  }


  # with standard deviation
  plot(x = beta,  xlab = "beta",
       y = alpha, ylab = "alpha",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = gmean(alpha), v = mean(beta), col = "blue")
  text(x = posX1, y = posY1,labels = "beta", col = "black", font = 3)
  text(x = posX2, y = posY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = beta[index], y = alpha[index],
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = TRUE, scale.log = FALSE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = posX2, y = as.numeric(talpha[index2, "means"]),
           labels = as.character(talpha[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = as.numeric(tbeta[index2, "means"]), y = posY1,
           labels = as.character(tbeta[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# SCALE LOG == TRUE  (opt.log == TRUE)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if (opt.log == TRUE) {


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Parameters for graph plotting
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
loghigth <- max(log2(alpha)) - min(log2(alpha))
logwidth <- max(log2(beta))  - min(log2(beta))
logposX1 <- min(log2(beta))  - 0.02*logwidth
logposX2 <- min(log2(beta))  + 1.02*logwidth
logposY1 <- min(log2(alpha)) - 0.02*loghigth
logposY2 <- min(log2(alpha)) + 1.02*loghigth



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot fobs
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "fobs"
ghisto(fobs, titre, scale.log = TRUE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot alpha
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "alpha"

ghisto(alpha, titre, scale.log = TRUE)
plot.alphaVSbeta(x = fobs,  xlab = "fobs",
                 y = alpha, ylab = titre,
                 titre = titre, scale.log = TRUE)

plot.alphaVSbeta(y = fobs,  ylab = "fobs",
                 x = alpha, xlab = titre,
                 titre = titre, scale.log = TRUE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot beta
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
titre <- "beta"

ghisto(beta, titre, scale.log = TRUE)
plot.alphaVSbeta(x = fobs, xlab = "fobs",
                 y = beta, ylab = titre,
                 titre = titre, scale.log = TRUE)

plot.alphaVSbeta(y = fobs, ylab = "fobs",
                 x = beta, xlab = titre,
                 titre = titre, scale.log = TRUE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot log2(alpha) vs log2(beta)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot.alphaVSbeta(x = beta,  xlab = "beta",
                 y = alpha, ylab = "alpha",
                 titre = "alpha vs beta", scale.log = TRUE)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot log2(alpha) vs log2(beta) WITH TRANSGRESSIVE ASSEMBLAGES in RED
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if (opt.aov == TRUE) {
  talpha <- test.posthoc(alpha, trans, pvalue)
  tbeta  <- test.posthoc(beta,  trans, pvalue)
}
setNoms <- c("trans-under", "under", "over", "trans-over")

plot(x = log2(beta),  xlab = "log2(beta)",
     y = log2(alpha), ylab = "log2(alpha)",
     main = "log2(alpha) vs log2(beta): transgressive overyielding",
     type = "n", tck = 0.02, las = 1)

abline(v = log2(mean(beta)), h = log2(mean(alpha)), col = "blue")
text(x = logposX1, y = logposY1,labels = "beta", col = "black", font = 3)
text(x = logposX2, y = logposY2,labels = "alpha", col = "black", font = 3)

setTrans <- sort(unique(trans), decreasing = TRUE)
for (i in 1:length(setTrans)) {
  j <- setTrans[i]
  points(x = log2(beta[(trans == j)]), y = log2(alpha[(trans == j)]),
         type = "p", pch = figures[i], col = couleurs[i],
         bg = "white", cex = 2)

  points.sd(x = beta[(trans == j)],  y = alpha[(trans == j)],
            figure = figures[i], couleur = couleurs[i], nom = setNoms[j],
            opt.std = TRUE, scale.log = TRUE)

  if (opt.aov == TRUE) if (is.list(talpha)) {
    index <- which(talpha[ , "motifs"] == j)
    text(x = logposX2, y = log2(as.numeric(talpha[index, "means"])),
         labels = as.character(talpha[index, "M"]),
         col =  couleurs[i], font = 3)
  }

  if (opt.aov == TRUE) if (is.list(tbeta)) {
    index <- which(tbeta[ , "motifs"] == j)
    text(x = log2(as.numeric(tbeta[index, "means"])), y = logposY1,
         labels = as.character(tbeta[index, "M"]),
         col =  couleurs[i], font = 3)
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the mean performance of assemblages containing a given element in LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
for (elt in seq_len(nbElt)) {
  index <- which(mOccur[,elt] == 1)

  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       main = paste("Element = ", colnames(mOccur)[elt], sep = ""),
       type = "n", tck = 0.02, las = 1)
  points(x = log2(beta[-index]), y = log2(alpha[-index]),
         type = "p", pch = figures[1], col = "black", bg = "white", cex = 2)
  points(x = log2(beta[index]), y = log2(alpha[index]),
         type = "p", pch = figures[2], col = "blue", bg = "white", cex = 2)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[2], couleur = "blue", nom = colnames(mOccur)[elt],
            opt.std = TRUE, scale.log = TRUE)
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(mean(beta)), col = "blue")
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITH RAW DATA  IN LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = log2(beta),  xlab = "log2(beta)",
     y = log2(alpha), ylab = "log2(alpha)",
     main = "All elements",
     pch = figures[1], col = "black", bg = "white",
     type = "p", cex = 2, tck = 0.02, las = 1)

abline(h = 0, v = 0, col = "black")
abline(h = log2(gmean(alpha)), v = log2(mean(beta)), col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = FALSE, scale.log = TRUE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITH RAW DATA WITH SD  IN LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = log2(beta),  xlab = "log2(beta)",
     y = log2(alpha), ylab = "log2(alpha)",
     main = "All elements",
     pch = figures[1], col = "black", bg = "white",
     type = "p", cex = 2, tck = 0.02, las = 1)

abline(h = 0, v = 0, col = "black")
abline(h = log2(gmean(alpha)), v = log2(mean(beta)), col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = log2(beta[index]), y = log2(alpha[index]),
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = TRUE, scale.log = FALSE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all species without standard deviation WITHOUT RAW DATA  IN LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = log2(beta),  xlab = "log2(beta)",
     y = log2(alpha), ylab = "log2(alpha)",
     #     xlim = c(-0.2,0.4), ylim = c(0.4,0.7),
     main = "All elements",
     type = "n", bg = "white", tck = 0.02, las = 1)
abline(h = 0, v = 0, col = "black")
abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

R2.mean <- summary(lm(formula = alpha ~ beta))$"r.squared"
text(x = 0.8*max(log2(beta)), y = min(log2(alpha)),
     labels = paste("R2.mean = ", signif(R2.mean, digits = 3),sep = ""),
     pos = 3, col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = FALSE, scale.log = TRUE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot all elements with standard deviation WITHOUT RAW DATA WITH SD IN LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(x = log2(beta),  xlab = "log2(beta)",
     y = log2(alpha), ylab = "log2(alpha)",
     #     xlim = c(-0.2,0.4), ylim = c(0.4,0.7),
     main = "All elements",
     type = "n", bg = "white", tck = 0.02, las = 1)
abline(h = 0, v = 0, col = "black")
abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

R2.mean <- summary(lm(formula = alpha ~ beta))$"r.squared"
text(x = 0.5 * max(log2(beta)), y = min(log2(alpha)),
     labels = paste("R2.mean = ", signif(R2.mean, digits = 3), sep = ""),
     pos = 3, col = "blue")

for (elt in 1:nbElt) {
  index <- which(mOccur[,elt] == 1)
  points.sd(x = beta[index], y = alpha[index],
            figure = figures[elt], couleur = couleurs[elt],
            nom = colnames(mOccur)[elt], opt.std = TRUE, scale.log = TRUE)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the hierarchical tree
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
plot(dot.arbre, xlab = "hierarchical tree",
     main = "Agregation Euclidien/Ward", hang = -1, las = 1)




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances CLUSTERED BY ASSEMBLY MOTIF  IN LOG
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
for (nbclasses in seq_len(nbclassesMax)) {
  # Plot one species
  affectElt <- mAffectElt[nbclasses, ]
  AssMotifs <- mAssMotifs[nbclasses, ]

  setMots   <- sort(unique(AssMotifs))
  nbMots    <- length(setMots)

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }

  noms <- names.assembly(affectElt, mOccur)
  plot.clusters.content(affectElt, mOccur)


  # only raw data without either mean neither standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    text(x = log2(gmean(beta[index])), y = log2(gmean(alpha[index])),
         col = couleurs[motif],
         labels = noms$motifs[motif], pos = 3, font = 3)
  }


  # without standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = FALSE, scale.log = TRUE)
  }


  # with standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")
  text(x = logposX1, y = logposY1,labels = "beta", col = "black", font = 3)
  text(x = logposX2, y = logposY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = TRUE, scale.log = TRUE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = logposX2, y = log2(as.numeric(talpha[index2, "means"])),
           labels = as.character(talpha[index2, "M"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = log2(as.numeric(tbeta[index2, "means"])), y = logposY1,
           labels = as.character(tbeta[index2, "M"]),
           col =  couleurs[motif], font = 3)
    }

  }


  # with Alpha and Beta by motif
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1,
       main = paste("number of classes = ", nbclasses, sep = ""))
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[mot], col = couleurs[mot], bg = "white",
           type = "p", cex = 2)

    for (elt in seq_len(nbElt)) {
      index2 <- which(mOccur[index,elt] == 1)
      titre  <- paste(noms$motifs[motif],
                      colnames(mOccur)[elt], sep = ".")
      points.sd(x = beta[index[index2]], y = alpha[index[index2]],
                figure = figures[elt], couleur = couleurs[mot],
                nom = titre,
                opt.std = FALSE, scale.log = TRUE)
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances clustered by the NUMBER of SPECIES
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
if (opt.size == TRUE) {

  mainTitre <- "Assembly size"

  size      <- apply(mOccur, MARGIN = 1, FUN = sum)
  AssMotifs <- numeric(length(size))
  for (i in seq_along(unique(size))) AssMotifs[(size == unique(size)[i])] <- i

  setMots <- sort(unique(AssMotifs))
  nbMots  <- length(setMots)
  noms    <- as.character(sort(unique(size)))

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # without standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

  for (mot in 1:nbMots) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif],
           type = "p", bg = "white", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms[motif],
              opt.std = FALSE, scale.log = TRUE)
  }


  # with standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")
  text(x = logposX1, y = logposY1,labels = "beta", col = "black", font = 3)
  text(x = logposX2, y = logposY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms[motif],
              opt.std = TRUE, scale.log = TRUE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = logposX2, y = log2(as.numeric(talpha[index2, "means"])),
           labels = as.character(talpha[index2, "M"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = log2(as.numeric(tbeta[index2, "means"])), y = logposY1,
           labels = as.character(tbeta[index2, "M"]),
           col =  couleurs[motif], font = 3)
    }
  }
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Plot the assembly performances clustered by A PRIORI FUNCTIONAL GROUPS
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

if (opt.gf == TRUE) {

  mainTitre <- "Functional Groups"

  affectElt <- gf
  AssMotifs <- affect.motifs(affectElt, mOccur)
  setMots   <- sort(unique(AssMotifs))
  nbMots    <- length(setMots)

  if (opt.aov == TRUE) {
    talpha <- test.posthoc(alpha, AssMotifs, pvalue)
    tbeta  <- test.posthoc(beta,  AssMotifs, pvalue)
  }

  noms <- names.assembly(affectElt, mOccur)
  plot.clusters.content(affectElt, mOccur)


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # without standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       main = mainTitre,
       type = "n", bg = "white", cex = 2, tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = FALSE, scale.log = TRUE)
  }


  # with standard deviation
  plot(x = log2(beta),  xlab = "log2(beta)",
       y = log2(alpha), ylab = "log2(alpha)",
       main = mainTitre,
       type = "n", cex = 2, bg = "white", tck = 0.02, las = 1)
  abline(h = 0, v = 0, col = "black")
  abline(h = log2(gmean(alpha)), v = log2(gmean(beta)), col = "blue")
  text(x = logposX1, y = logposY1,labels = "beta", col = "black", font = 3)
  text(x = logposX2, y = logposY2,labels = "alpha", col = "black", font = 3)

  for (mot in seq_len(nbMots)) {
    motif <- setMots[mot]
    index <- which(AssMotifs == motif)
    points(x = log2(beta[index]), y = log2(alpha[index]),
           pch = figures[motif], col = couleurs[motif], bg = "white",
           type = "p", cex = 2)
    points.sd(x = beta[index], y = alpha[index],
              figure = figures[motif], couleur = couleurs[motif],
              nom = noms$motifs[motif],
              opt.std = TRUE, scale.log = TRUE)

    if (opt.aov == TRUE) if (is.list(talpha)) {
      index2 <- which(talpha[ , "motifs"] == motif)
      text(x = logposX2, y = log2(as.numeric(talpha[index2, "means"])),
           labels = as.character(talpha[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }

    if (opt.aov == TRUE) if (is.list(tbeta)) {
      index2 <- which(tbeta[ , "motifs"] == motif)
      text(x = log2(as.numeric(tbeta[index2, "means"])), y = logposY1,
           labels = as.character(tbeta[index2, "groups"]),
           col =  couleurs[motif], font = 3)
    }
  }
}

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# END of SCALE LOG == TRUE  (opt.log == TRUE)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# END of SCRIPT
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


