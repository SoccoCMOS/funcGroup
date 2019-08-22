
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                    MYSTATS.R
#
#                 Set of mathematical and statistical functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Clear Memory, load files and global variables
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


library(multcompView)


EPSILON <- 1e-15





#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# List of functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx








#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the arithmetic mean and std
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the arithmetic mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

amean <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the geometric standard deviation
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

asd <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  mu <- amean(x)
  return( ifelse(length(x), sqrt( sum((x - mu) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the weighted mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wamean <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the weighted standard deviation
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wasd <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  res <- NA
  if (length(x) > 0) {
    tmp <- sum(x * w) / sum(w)
    res <- sum( w * ((x - tmp) ^ 2) ) / sum(w)
  }

  return( sqrt(res) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the geometric mean and std
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the geometric mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

#  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the geometric standard deviation
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gsd <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  return( ifelse(length(x),
                 exp( sqrt( sum( log(x/gmean(x)) ^ 2 ) / length(x) ) ),
                 NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the weighted geometric mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wgmean <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the geometric standard deviation
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the harmonic mean and std
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the harmonic mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

hmean <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  return( ifelse(length(x), length(x) / sum(1/x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the weight harmonic mean value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

whmean <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  return( ifelse(length(x), sum(w) / sum(w/x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the mean Prediction computed
#                    by excluding (LOO) the assembly to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mean.fct <- function(fct, opt.mean) {

  return(switch(opt.mean,
                amean = amean(fct),
                gmean = gmean(fct)
  ))
}

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Output the R2 protected for var(X)=0 or var(Y)=0
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the biased variance
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

var2 <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  res <- NA
  if (length(x) > 0) {
    mu  <- mean(x)
    res <- sum((x - mu) ^ 2) / length(x)
  }

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the moment of 3rd order
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

var3 <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  res <- NA
  if (length(x) > 0) {
    mu  <- mean(x)
    res <- sum((x - mu) ^ 3) / length(x)
  }

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the moment of 'th order
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

var4 <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  res <- NA
  if (length(x) > 0) {
    mu  <- mean(x)
    res <- sum((x - mu) ^ 4) / length(x)
  }

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Pearson' R2 of a linear regression
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cor2 <- function(x, y, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(y))
    x     <- x[index]
    y     <- y[index]
  }

  res <- NA
  if (length(x) > 0)
  {
    z1 <- var(x)
    z2 <- var(y)
    z3 <- cov(x,y)

    if ((z1 > 0) && (z2 > 0)) res <- (z3 ^ 2) / (z1 * z2)
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Functions
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

r2z  <- function(r) { return( log( (1 + r)/(1 - r) ) / 2 ) }
z2r  <- function(z) { return( (exp(2*z) - 1) / (exp(2*z) + 1) ) }



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#
#  Output the Prediction error
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Residual Sum of Square (RSS for "Residual") of Predicted values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss <- function(fcal, fct) {

  res     <- sum( (fcal - fct) ^ 2, na.rm = TRUE )
  if (res < EPSILON) res <- 0

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Mean Square Error of Predicted values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mse <- function(pred, obs, na.rm = TRUE)
{
  if (na.rm == TRUE) {
    index <- !(is.na(pred) | is.na(obs))
    pred  <- pred[index]
    obs   <- obs[index]
  }

  return( ifelse(length(obs), sum((pred - obs) ^ 2) / length(obs), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the ROOT Mean Square Error of Predicted values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rmse <- function(pred, obs, na.rm = TRUE) {

  return( sqrt(mse(pred, obs, na.rm)) )

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Coefficient of Variation of ROOT Mean Square Error
#       of Predicted values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cvrmse <- function(pred, obs, na.rm = TRUE) {

  return( rmse(pred, obs, na.rm) / amean(obs) )

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Coefficient of Determination R2 of Predicted values
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2mse <- function(pred, obs, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(pred) | is.na(obs))
    pred  <- pred[index]
    obs   <- obs[index]
  }

  res <- NA
  if (length(obs) > 1) {
    tss   <- sum((obs - amean(obs)) ^ 2)
    rss   <- sum((obs - pred) ^ 2)
    res   <- (tss - rss) / tss
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Probability associated to a R2mse
#     ( with "r" for "Residual", "e" for "Explained" and "t" for "Total")
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pmse <- function(prd, obs, nbK, na.rm = TRUE) {

  if (nbK == 1) stop("nbK should be higher than one")

  if (na.rm == TRUE) {
    index <- !(is.na(prd) | is.na(obs))
    prd   <- prd[index]
    obs   <- obs[index]
  }

  res <- 1
  if (length(obs) > 1) {
    tss   <- sum((obs - amean(obs)) ^ 2)
    rss   <- sum((obs - prd) ^ 2)
    ess   <- tss - rss

    nbfdt <- length(obs)
    nbfde <- nbK - 1
    nbfdr <- nbfdt - nbfde

    if ( (nbfdr > 0) && (nbfde > 0) ) {
       Fratio <- (ess / nbfde) / (rss / nbfdr)
       res    <- 1 - pf(Fratio, nbfde, nbfdr)
    }
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the associated CVmse when knowing a Coefficient of Determination R2
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2toCVmse <- function(R2, obs) {

  return( sqrt( (1 - R2) * var2(obs) ) / amean(obs) )

}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the associated CVmse when knowing a Coefficient of Determination R2
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2topmse <- function(R2, obs, nbK, na.rm = TRUE) {

  if (nbK == 1) stop("nbK should be higher than one")

  if (na.rm == TRUE) obs <- obs[!is.na(obs)]

  res <- 1
  if (length(obs) > 1) {
    tss <- sum((obs - amean(obs)) ^ 2)
    ess <- R2 * tss
    rss <- tss - ess

    nbfdt <- length(obs)
    nbfde <- nbK - 1
    nbfdr <- nbfdt - nbfde

    if ( (nbfdr > 0) && (nbfde > 0) ) {
      Fratio <- (ess / nbfde) / (rss / nbfdr)
      res    <- 1 - pf(Fratio, nbfde, nbfdr)
    }
  }

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the AIC computed from the RSS value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

AIC <- function(fprd, fobs, K) {

  res <- NA
  n   <- length(fobs)
  RSS <- rss(fprd, fobs)
  if (RSS > EPSILON) {
    res <- n * log(rss(fprd, fobs) / n) + 2 * K
    if (res == -Inf) res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the AICc computed from the RSS value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

AICc <- function(fprd, fobs, K) {

  res <- NA
  n   <- length(fobs)
  RSS <- rss(fprd, fobs)
  if ( (RSS > EPSILON) & (n + 1 > K) ) {
    res <- n * log(RSS / n) + 2 * K * (n + 2) / (n + 1 - K)
    if (res == -Inf) res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the Square BIAS
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

bias2 <- function(pred, obs, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(pred) | is.na(obs))
    pred  <- pred[index]
    obs   <- obs[index]
  }

  return( ifelse(length(obs) - 1,
                 mse(pred, obs, na.rm) - var2(pred, na.rm), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the index of the (first) minimum (or maximum) on a vector
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmin <- function(v) {

  arg <- which(v == min(v, na.rm = TRUE))

  return(arg)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first.argmin <- function(v) {

  v   <- v[!is.na(v)]

  arg <- 1
  while ((v[arg + 1] < v[arg]) & (arg < length(v))) arg <- arg + 1

  return(arg)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmax <- function(v) {

  arg <- which(v == max(v, na.rm = TRUE))

  return(arg)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first.argmax <- function(v) {

  v   <- v[!is.na(v)]

  arg <- 1
  while ((v[arg + 1] > v[arg]) & (arg < length(v))) arg <- arg + 1

  return(arg)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the number of Stirling of second kind
#           = the number of possible partitions into k custers among n elements
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

stirling <- function(n) {

  res        <- integer(n)
  names(res) <- seq_len(n)

  for (k in seq_len(n)) {
    tmp <- 0
    for (j in 0:k) tmp <- tmp + (-1) ^ (k - j) * choose(k, j) * j ^ n
    res[k] <- tmp / factorial(k)
  }

  return(res)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Variance analysis with HSD posthoc test
#
#  Inputs : x        = vector of data,
#           clusters = vector of factors
#           pvalue   = value of 1er risk of first species
#  Output : results
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test.posthoc <- function(x, clusters, pvalue = 0.01) {

  # check the inputs
  if (!is.factor(clusters)) clusters <- as.factor(clusters)

  if (sum(is.na(x)) > 0) {
    x        <- x[!is.na(x)]
    clusters <- clusters[!is.na(x)]
  }

  levels  <- sort(unique(clusters))
  nblevel <- length(levels)

  res   <- NA
  index <- which(table(clusters) > 1)
  if (length(index) > 1) {

    # compute the mean values of groups
    means  <- numeric(nblevel)
    for (i in seq_len(nblevel)) means[i] <- amean(x[clusters == levels[i]])
    tmp    <- sort(means, decreasing = TRUE, index.return = TRUE)

    # analyse the variances
    model  <- aov(x ~ clusters)
    test   <- TukeyHSD(x = model, "clusters", conf.level = (1 - pvalue))
    vlet   <- multcompLetters4(model, test)

    motifs <- levels[tmp$ix]
    means  <- tmp$x
    groups <- (vlet$clusters)$Letters
    mlet   <- (vlet$clusters)$LetterMatrix

    res    <- data.frame(motifs, means, groups, mlet)
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Drawing of the linear regression "Predicted.values vs Observed.values"
#
#  Inputs : reg = regression model,
#           pvalue = value of 1er risk of first species
#  Output : plots
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

confidence.intervals <- function(xlim, reg, pvalue) {

  x   <- (reg$model)[,2]

  n   <- length(x)
  xm  <- amean(x)
  xss <- sum(x ^ 2) - sum(x) ^ 2 / n
  tss <- qt(1 - pvalue/2, (n - 2))

  xx  <- seq(xlim[1], xlim[2], length = 20)

  inter.confidence <-
    tss * sqrt(summary(reg)$sigma ^ 2 * (1/n + (xx - xm) ^ 2 / xss))

  inter.prediction <-
    tss * sqrt(summary(reg)$sigma ^ 2 * (1 + 1/n + (xx - xm) ^ 2/xss))

  res        <- list(xx, inter.confidence, inter.prediction)
  names(res) <- c("x", "confidence", "prediction")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Test if two dependent R2
#  Inputs: v1, v2 : the two dependent data-vectors or R2 to compare
#          v3     : the vector or R2 between the two dependent data-vectors
#          n      : length of vectors
#     v1, v2 and V3 can be r.squared between data-vectors or data-vectors
#  Output: the p.value
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test.dependent.R2 <- function(v1, v2, v3, n = length(v1)) {

  if ((length(v1) > 1) || (length(v2) > 1) || (length(v3) > 1))
  {
    if  ((length(v1) != n) || (length(v2) != n) || (length(v3) != n))
      stop("The vectors v1, v2 and v3 should have a length of n")

    r13 <- abs(cor(v1, v3, method = "pearson"))
    r23 <- abs(cor(v2, v3, method = "pearson"))
    r12 <- abs(cor(v1, v2, method = "pearson"))

  } else {
    r13 <- sqrt(v1)
    r23 <- sqrt(v2)
    r12 <- sqrt(v3)
  }

  R      <- diag(1,3)
  R[2,1] <- R[1,2] <- r13
  R[3,1] <- R[1,3] <- r23
  R[2,3] <- R[3,2] <- r12

  if ((r13 == r23) || (r12 == 1)) {
    p.value <- 1
  } else {
    tobs <- abs(r13 - r23) *
      sqrt((1 + r12) /
           (2*abs(det(R)) / (n - 3) + (r13 + r23) ^ 2 * (1 - r12) ^ 3 /
                            (4*(n - 1))))
    p.value <- 2 * (1 - pt(tobs, (n - 3)))
  }

  return(p.value)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute the differences between successive R2 are significant
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test.dependent.R2mse <- function(pred2, pred1, obs) {

  R1 <- R2mse(pred2, obs)
  R2 <- R2mse(pred1, obs)
  R3 <- R2mse(pred2, pred1)

  test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
  test2 <- (is.na(R1) || is.na(R2) || is.na(R3))

  res  <- NA
  if ( !((test1 || test2) == TRUE) )
    res <- test.dependent.R2(R1, R2, R3, length(obs))

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute the differences between successive R2 are significant
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pvalue.dependent.R2mse <- function(mpred, obs) {

  nbline <- dim(mpred)[1]
  res    <- numeric(nbline)

  for (lin in 2:nbline) {

    index <- na.action(na.omit(mpred[lin,]))
    if (length(index) > 0) {

      pred  <- mpred[lin,     -index]
      ppred <- mpred[lin - 1, -index]
      pobs   <-   obs[-index]

    } else {

      pred  <- mpred[lin,    ]
      ppred <- mpred[lin - 1,]
      pobs   <- obs
    }

    R1 <- R2mse(pred,  pobs)
    R2 <- R2mse(ppred, pobs)
    R3 <- R2mse(ppred, pred)

    test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
    test2 <- (is.na(R1) || is.na(R2) || is.na(R3))

    res[lin] <- NA
    if (!((test1 || test2) == TRUE))
      res[lin] <- test.dependent.R2(R1, R2, R3, length(pobs))
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return a R2-star for a clustering
#
#  "adjusted" returns the index proposed by Theil (1978) for multi-regression
#  "modified" returns the index proposed by GoldBerger (1991)
#  "Calenski" returns the index proposed by Calinski-Harabasz (1974)
#                               for clustering
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

R2Star.cluster <- function(R2, nbcl, nbelt, opt.optima) {

  switch(opt.optima,
         "adjusted": { 1 - (1 - R2) * (nbelt - 1) / (nbelt - nbcl) },
         "modified": { (1 - (nbcl - 1) / (nbelt - 1)) * R2         },
         "calinski": { R2 / (1 - R2) * (nbcl - 1) / (nbelt - nbcl) }
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

parsimony <- function(nbAss, nbcl) {

  return( ifelse(!sum(nbcl == 0), sum(nbAss/nbcl) / sum(nbAss), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

normalize.variable <- function(fct, option = "Johnson") {

  if (option == "BoxCox") {

    library(forecast)

    lambda <- BoxCox.lambda(fct)
    res    <- BoxCox(fct, lambda)
  }

  if (option == "Johnson") {

    library(Johnson)

    tmp <- RE.Johnson(fct)
    res <- tmp$transformed
  }

  return(res)

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  as "rnorm" function, but with values comprised between min and max
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rnorm.trunc <- function(n, mean = 0, sd = 1, min, max) {

  res <- NULL
  while (length(res) < n) {
    x     <- rnorm(n - length(res), mean, sd)
    index <- which( (x >= min) & (x <= max) )
    x     <- x[index]
    res   <- c(res, x)
  }

  return(res)
}






#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  END of FILE
#
##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

