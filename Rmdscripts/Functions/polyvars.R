polyvars <- function(X, degr = 2, raw = FALSE) 
{

# A function computing polynomials of vectors within a matrix.
# Contrary to function poly() on which it is based, this function 
# only computes polynomials separately for each vector of the provided matrix,
# e.g. x, x^2, x^3, and not combinations such as xy, x^2y and so on.
#
# Author: Daniel Borcard, December 2014, March 2017
# License: GPL2
#
# Usage
# -----
# polymatrix(X = rawdatamatrix, degr = 3, raw = FALSE)
#
# Arguments
# ---------
#    X: a matrix or data frame containing quantitative variables
#
#    degr: the degree to which the variables must be raised. Default: 2
#
#    raw: logical; if TRUE raw polynomials are computed directly from 
#         the raw variables. If FALSE (default), orthogonal polynomials 
#         are computed.
#
# Value
# -----
# A data frame containing the polynomials. In the output matrix, each
# variable appears in turn, followed by its polynomial terms, e.g.
# v1, v2_square, v2, v2_square, and so on.
#
# Details
# -------
# When raw = FALSE, the function computes orthogonal polynomial terms 
# of each variable separately. This means that in the resulting matrix
# the polynomial terms of each variable are orthogonal, but that they
# are not orthogonal to the terms of the other variables.

class.verif <- apply(X, 2, class)
if (any(class.verif == "factor") | any(class.verif == "character") == TRUE)
stop("No factor or character variables allowed.", call. = FALSE)

## Store or assign variable names
if(!is.null(colnames(X)))
{
   var.names <- colnames(X)
}
else
{
   var.names <- paste("v", 1 : ncol(X), sep = "")
}

## Compute polynomial terms
X.poly <- matrix(0, nrow(X), ncol(X) * degr)
for(i in 0: (ncol(X) - 1)) 
{
   toto <- poly(X[, (i + 1)], degr, raw=raw)
   X.poly[,(i * degr + 1) : ((i + 1) * degr)] <- toto
 }
                         
if((ncol(X) * degr) > (nrow(X) - 1) ) 
{
cat("\n------------------------------------------------------------------")
cat("\nWARNING: the number of polynomial terms is equal to or larger than")
cat("\nthe number of observations.")
cat("\n------------------------------------------------------------------\n")
}

## Create new column names
indices <- rep(1 : degr, ncol(X))
tmp <- vector(length = ncol(X.poly))
for(j in 1 : ncol(X))
{
   tmp[(j * degr - degr + 1) : (j * degr)] <- rep(var.names[j], degr)
}
var.poly <- paste(tmp, indices, sep = ".")
colnames(X.poly) <- var.poly

X.poly.df <- as.data.frame(X.poly)

X.poly.df

}


## Examples

## Construction of a fictitious matrix of 5 observations and 4 variables:
# env <- matrix(1:20, 5)

## Computation of orthogonal polynomials of degree 3:
# env.ortho.deg3 <- polymatrix(env, degr = 3)

## Computation of a matrix of raw polynomials of degree 4:
# env.raw.deg4 <- polymatrix(env, degr = 4, raw = TRUE)

