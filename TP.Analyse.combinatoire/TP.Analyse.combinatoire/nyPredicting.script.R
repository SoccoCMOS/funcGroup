#==============================================================================
#==============================================================================
#
#                               myPREDICTING.R
#
#==============================================================================
#==============================================================================


#==============================================================================
# Clear Memory, load files and global variables
#==============================================================================


#==============================================================================
#==============================================================================
#
#   SCRIPT for Predicting performance and Plotting the Results
#
#==============================================================================
#==============================================================================

opt.R2   <- TRUE
opt.cal  <- FALSE
opt.prd  <- FALSE
opt.glb  <- FALSE
opt.pub  <- FALSE
opt.aov  <- FALSE

opt.size <- FALSE
opt.gf   <- FALSE
opt.box  <- FALSE
opt.all  <- FALSE


opt.size <- TRUE
opt.gf   <- TRUE
#opt.box <- TRUE
#opt.all <- TRUE


opta <- "gmean"
optb <- "amean"

#==============================================================================
# Plot the predicted ALPHA, BETA and FOBS by ASSEMBLY MOTIFS
#==============================================================================

res.alpha <- predict.function(tree.div, mOccur, alpha,
                              opt.mean = opta,
                              opt.mod  = "byelt")

plot.prediction(res.alpha, titre = "alpha",
                opt.aov = FALSE, pvalue = 0.05,
                opt.all = TRUE)


res.beta <- predict.function(tree.div, mOccur, beta,
                             opt.mean = optb,
                             opt.mod  = "byelt")

plot.prediction(res.beta, titre = "beta",
                opt.aov = FALSE, pvalue = 0.05,
                opt.all = TRUE)


res.abfobs <- predict.twin(tree.div, mOccur, alpha, beta,
                           fscale = 1,
                           opt.alpha = opta,
                           opt.beta  = optb,
                           opt.mod  = "byelt")

plot.prediction(res.abfobs, titre = "ab = fobs",
                opt.aov = FALSE, pvalue = 0.05,
                opt.all = TRUE)




#==============================================================================
#==============================================================================
# END of SCRIPT
#==============================================================================
#==============================================================================



