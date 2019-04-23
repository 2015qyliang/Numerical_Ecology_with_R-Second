screestick <- function(ev, las = 1, positive = TRUE, 
                       evmean = FALSE, relative = FALSE) 
{

# Description

# A function to draw a screeplot of a vector of eigenvalues
# with overlay of the values predicted by the broken stick model
#
# This function mimics the default output of function 
# screeplot.cca {vegan} but it works on a vector of eigenvalues, 
# which makes it independent of the function and package that has 
# computed the ordination.

# Arguments

# ev       vector of eigenvalues
# las      orientation of the x axis labels
# ww       width of bars in barplot
# ss       space between bars in barplot
# evmean   logical: should the mean of the eigenvalues be overlaid?
# relative logical: should the eigenvalues be divided by their
#          sum to represent relative proportions of variation?

# Author
# Daniel Borcard, Universite de Montreal, September 2017
# Licence: GPL-2

require(vegan)

if(is.vector(ev) == FALSE) 
  stop("Please provide a *vector* of eigenvalues", call. = FALSE)

evs <- sort(ev, decreasing = TRUE)

if(sum(abs(ev - evs)) > 0)
{
  cat("The values were not in decreasing order. They have been sorted.\n")
}

if(positive == TRUE)
{
  ntemp <- length(evs)
  evs <- evs[evs > 0]
  n <- length(evs)

  if(ntemp > n) 
  {
    cat("Warning: presence of", ntemp - n, "negative eigenvalues. \n")
    cat("Only the", n, "positive eigenvalues are considered.\n")
  }
}
else
{
  n <- length(evs)
}

if(relative == TRUE)
{
evs <- evs/sum(evs)
}

names(evs) <- paste("EV", 1 : n)

# Broken stick model
broken <- bstick(n) * sum(evs)

if(relative == TRUE)
  {
  ylabel = "Relative eigenvalue"
  }
else
  {
  ylabel = "Eigenvalue"
  }

# Scree plot
barplot(evs,
        ylim = c(0, max(c(evs, broken))),
        main = "Eigenvalues and broken stick model",
        ylab = ylabel,
        las = las)

# Overlay the lines and points representing the broken stick model.
# The bars of the bar plot have width ww = 1 and are preceded and
# separated by an interval ss = 0.2. There are n bars. Argument 
# 'absc' below is computed to ensure that the line breaks and the 
# corresponding points of the broken stick model are aligned with 
# the centers of the bars of the scree plot.
# ww and ss might be added as arguments in a later version.
ww <- 1
ss <- 0.2
ws <- ww + ss
n05 <- ww/2
absc <- seq((n05 + ss), ((n * ws) - n05), by = ws)
lines(
  absc, 
  broken,
  type = "o",
  pch = 1,
  col = "red"
)

if(evmean == TRUE) 
  {
    abline (h = mean(evs), col = "darkgray")
legend(
  "topright",
  legend = c("Broken Stick", "Mean of eigenvalues"),
  pch = c(1, NA_integer_),
  lty = 1, 
  col = c("red", "gray"),
  bty = "n")
  }
else
  {
legend(
  "topright", 
  "Broken Stick", 
  pch = 1, 
  lty = 1, 
  col = "red", 
  bty = "n")
  }
}
