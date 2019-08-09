scalog <- function(res, np = 999, alpha = c(0.05, 0.01, 0.001), cex=2)
{

# A function to compute a scalogram (Legendre and Legendre 2012, 
# p. 864) representing the eigenvalues of an RDA with a series of 
# spatial eigenfunctions (e.g. dbMEM) as explanatory variables.
# The eigenfunctions must be placed in decreasing order and they
# must be orthogonal to one another.
# In case of RDA the R^2 is variance, in case of CCA it is inertia.

# Arguments

# res    An RDA or CCA result object produced by vegan::rda().
#        The RDA or CCA must have been computed with the formula
#        interface.
# np     number of permutations in the RDA test
# alpha  probability thresholds for the color coding

# License: GPL-2 
# Author: Daniel Borcard, 2017

test <- anova(res, by = "terms", permutations = how(nperm = np))
inert <- test[,2]
variance <- inert[-length(inert)] / sum(inert)
signi <- test$"Pr(>F)"[-length(test$"Pr(>F)")]
n <- length(variance)

if(class(res)[1] == "rda") 
  ylabel <- expression(italic(R)^{2})
else
  ylabel <- "Inertia"
   
plot(
  1 : n, 
  variance, 
  type = "n", 
  main = "Scalogram", 
  xlab = "Eigenfunction number", 
  ylab = ylabel
)
lines(1 : n, variance)
alpha <- c(1, sort(alpha, decreasing = TRUE))
colour <- c("white", "greenyellow", "orange", "red")
for(i in 1 : 4)
{
  points(
    (1 : n)[signi <= alpha[i]], 
    variance[signi <= alpha[i]], 
    pch = 22, 
    cex = cex,
    bg = colour[i]
  )
}
legend(
  "topright",
  fill = c("red", "orange", "greenyellow"),
  legend = c(paste("p <= ", alpha[4]), paste("p <= ", alpha[3]), 
             paste("p <= ", alpha[2]))
)
invisible(test)
}

