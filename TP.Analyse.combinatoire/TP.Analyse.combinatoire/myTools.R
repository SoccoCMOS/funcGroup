
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#           Set of Functions et Procedures facilitating Input-Output
#
#                     (Benoit JAILLARD, February 2011)
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


# what           <- function(questions = "", answers = "")
# check.symbol   <- function(symbol, long)

# open.pdf       <- function(opt.prn = FALSE, filename = "", nbGraph = 1)
# close.pdf      <- function(opt.prn = FALSE)

# open.txt       <- function(opt.prn = FALSE, filename= "")
# write.txt      <- function(opt.prn = FALSE, filename = "", txt = "")



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  loading of usefull libraries and local sources
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

library(Cairo)



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  defining usefull constantes
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# to avoid program failure with error message "figure margins too large"
par(mar = c(1, 1, 1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

figures  <- c(21, 22, 24, 23, 25)
# 21 = "circle", 22 = "square", 23 = "diamond",
# 24 = triangle point-up", 25 = "triangle point-down"

couleurs <- c("red3", "blue2", "orange2", "turquoise3", "magenta", "green4",
              "pink", "violet", "salmon4", "skyblue2", "sienna3", "olivedrab3")

myLetters <- c(letters, LETTERS)

tmp  <-  hsv(h = seq(2/3, 1, (1-2/3)/21), s = 1.0, v = 1.0, alpha = 1)
# for hue (couleur), saturation (intensité), value (luminosité)
# h = 2/3 "blue" h = 1 = "red"



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#   Functions and Procedures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Function : to interact with a locutor
#  Input: a list of questions
#  Output: a list of answers
#          (by convention, 0 for No or False, and 1 for Yes or T
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

what <- function(questions = "", answers = "") {

  nbQuestions <- length(questions)
  questions   <- as.data.frame(questions)
  reponses    <- as.data.frame(answers)

  questions   <- cbind(questions, reponses)
  colnames(questions) <- c("My questions", "Your answers")

  comments    <- as.data.frame(matrix(
                     c(" ", "To answer to a question:", "To close the window:",
                       " ", "TRUE/FALSE=1/0",          "click on X"),
                       nrow = 3, ncol = 2))
  colnames(comments) <- c("My questions", "Your answers")
  questions   <- rbind(questions, comments)

  reponses    <- edit(questions)

  res         <- list(reponses[1:nbQuestions, 2])
  names(res)  <- c("answers")

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Return the index of "chr" in the string "str"
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


which.str <- function(str, chr) {

  i   <- 1
  while ( (substr(str, i, i) != chr) & (i < nchar(str)) ) i <- i + 1
  if (i == nchar(str)) i <- 0

  return(i)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Function : input to generate vectors of colours and symbols
#                             as many long as the number of studied communities
#  Input:  an index for sorting elements according a criterion
#  Output: the index that allows to reorganize in turn the elements-
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

index.inturn <- function(index) {

  res <- index
  for (i in seq_along(index)) res[index[i]] <- i

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Function : to generate vectors of colours and symbols
#                             as many long as the number of studied communities
#  Input: base symbol vector, and a vector of elements to plot
#  Output: a vector of symbols of length > length(vector)-
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

check.symbol <- function(symbol, long) {

  while (length(symbol) < long) symbol <- c(symbol, symbol)

  return(symbol)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Function : to generate vectors of colours and symbols
#                             as many long as the number of studied communities
#  Input: base symbol vector, and a vector of elements to plot
#  Output: a vector of symbols of length > length(vector)-
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

check.myLetters <- function(long) {

  iter   <- 0
  symbol <- myLetters
  while (length(symbol) < long) {
    iter   <- iter + 1
    symbol <- c(symbol, paste0(myLetters, iter))
  }

  return(symbol)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Procedure: to open a graphic window (if option is FALSE)
#                        or a PDF-file (if option is TRUE)
#
#  Input : option to record or not a PDF-file in the current directory
#          nbGraph = number of graphs (between 1 and 4) to record
#          filename = title of graph AND name of PDF-file
#  Output: no
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

open.pdf <- function(opt.prn = FALSE, filename = "", nbGraph = 1) {

  if (!is.logical(opt.prn))  stop("the 'opt.prn' should be logical")

  if (opt.prn == TRUE) {

    if (length(filename) == 0) stop("'filename' should be a non-empty string")
    if (nbGraph < 1)           stop("the option 'nbGraph' should be no null")

    nb      <- floor(sqrt(nbGraph))
    lb      <- ifelse(nbGraph %% (nb ^ 2), nb + 1, nb)
    layout(matrix(seq_len(lb ^ 2), ncol = lb, nrow = lb, byrow = TRUE))

    Cairo::CairoPDF(file = paste(filename, "pdf", sep = "."),
                    onefile = TRUE, encoding = "default",
                    height = 8, width = 11,
                    family = "Helvetica", pointsize = 12,
                    bg = "transparent")
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Procedure: to Close a graphic window,
#             which was previously open by "open.windows"
#
#  Input:  option=TRUE to save on disk, FALSE otherwise
#  Output: no
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

close.pdf <- function(opt.prn = FALSE) {

  if (opt.prn == TRUE) dev.off()
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Procedure: to Open a String CSV file
#
#  Input:  option=TRUE to save on disk, FALSE otherwise
#  Output: no
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

open.txt <- function(opt.prn = FALSE, filename = "") {

  if (opt.prn == TRUE) {

    if (length(filename) == 0) stop("'filename' should be a non-empty string")

    write(x      = c(as.character(filename), "\n"),
          file   = paste(filename, "csv", sep = "."),
          append = FALSE)
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  Procedure: to Open a String CSV file
#
#  Input:  option=TRUE to save on disk, FALSE otherwise
#  Output: no
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

write.txt <- function(opt.prn = FALSE, filename = "", txt = "") {

  if (opt.prn == TRUE) {

    if (length(filename) == 0) stop("'filename' should be a non-empty string")

    write( x      = c(txt, "\n"),
           file   = paste(filename, "csv", sep = "."),
           append = TRUE)
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#  boxplot with mean
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

box.plot <- function(y, x, main,
                     ylim = range(y),
                     xlim = c(names(table(x)[1]),
                              names(table(x)[length(table(x))])) ) {

  xxlim <- c(as.integer(xlim[1]) - 1.5, as.integer(xlim[2]) - 0.5)
  boxplot(y ~ x, main = main, ylim = ylim, xlim = xxlim,
          tck = 0.02, las = 1)

  tab   <- table(x)
  moyen <- numeric(length(tab))
  for (i in seq_along(tab)) moyen[i] <- mean(y[(x == names(tab)[i])])

  points(y = moyen, x = seq_along(tab), pch = 0)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#
#   END of FILE myTOOLS.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx=

