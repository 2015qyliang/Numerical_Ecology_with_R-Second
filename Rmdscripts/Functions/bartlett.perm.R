bartlett.perm <- function(y, fact, centr="MEDIAN", nperm=999, alpha=0.05)

# Computation of parametric, permutational and bootstrap versions of the
# Bartlett test of homogeneity of variances.
#
# The data are centred to their within-group medians (default) or means 
# before the tests.
#
# Prior to the computation of the test of homogeneity of variances,
# a Shapiro-Wilk test of normality of residuals is computed. If the residuals
# are not normally distributed, a warning is issued because this nonnormality
# influences the type I error of the parametric test, which will very likely 
# have an incorrect type I error.
#
# USAGE
# bartlett.perm(y, fact, centr, nperm, alpha)
#
# ARGUMENTS
# y       a numeric vector of data values
# fact    a vector or factor object giving the group for the corresponding
#         elements of y
# centr   should the data, within groups, be centred on their medians ("MEDIAN")
#         or on their means ("MEAN")?
# nperm   number of permutations
# alpha   level of rejection of the H0 hypothesis of normality of residuals
#         in the Shapiro-Wilk test
#
# RESULT
# Bartlett          Bartlett's K-squared test statistic
# Param.prob        Parametric probability (P-value) of the test
# Permut.prob       Permutational probability (P-value) of the test
# Bootstrap.prob    Bootstrap probability (P-value) of the test
#
# DETAILS
#
# Centring the groups on their median or mean is very important for permutation
# and bootstrap tests to be correct when the groups do not share the same 
# position. Permuting groups with unequal mean or median artificially increases 
# the within-group variances of the permuted data.
#
#                         Daniel Borcard
#                         Universite de Montreal
#                         1 February 2016

# ------------------------------------------------------
# EXAMPLE
#
# Species abundance type data:
# y1 <- log1p(round(rlnorm(5,0.2,2)))
# y2 <- log1p(round(rlnorm(5,1,2)))
# y3 <- log1p(round(rlnorm(5,2,5)))
# yy <- c(y1,y2,y3)
#
# Factor
# fac <- gl(3,5, labels=c("groupe1","groupe2","groupe3"))
#
# Bartlett test with centring on the group medians
# bartlett.perm(yy, fac, centr="MEDIAN", nperm=999, alpha=0.05)
# ------------------------------------------------------


{

fact <- as.factor(fact)

normal <- shapiro.test(resid(lm(y~fact)))
if(normal[2]<=alpha){
cat("\n-------------------------------------------------------------------")
cat("\nWARNING") 
cat("\nThe residuals of the ANOVA model are not normally distributed.")
cat("\nThis is likely to change the rate of type I error of the test")
cat("\nof homogeneity of variances.")
cat("\n-------------------------------------------------------------------")
cat("\n")
}

# Trap for groups with 0 dispersion
y.sd <- tapply(y, fact, sd)
if(any(y.sd==0)) {
cat("\n-------------------------------------------------------------------")
cat("\nPROBLEM ENCOUNTERED") 
cat("\nOne or more groups have zero variance. Please chek and correct.")
cat("\nThe computations are meaningless if a group is made of observations")
cat("\nthat all have the same value.")
cat("\n-------------------------------------------------------------------")
cat("\n")
stop
}


CENTRE <- c("MEDIAN", "MEAN")
    centr <- match.arg(centr, CENTRE)


# Within-group centring of data

if(centr == "MEDIAN"){
  meds <- tapply(y, fact, median, na.rm=TRUE)
  y.c <- y - meds[fact]
                   }
else{
  means <- tapply(y, fact, mean, na.rm=TRUE)
  y.c <- y - means[fact]
    }


bart <- bartlett.test(y.c,fact)

# Permutation tests

cat("\nPermutations running...")
cat("\n")

compt.perm <- 1

for(i in 1:nperm) {

yprime <- sample(y.c)

bart.perm <- bartlett.test(yprime,fact)
if(bart.perm[[1]] >= bart[[1]]){
   compt.perm=compt.perm+1}

}

# Bootstrap tests
# Difference with permutation test: resampling is done with replacement

cat("\nBootstrap running...")
cat("\n")
cat("\n")

compt.boot <- 1

for(i in 1:nperm) {

yboot <- sample(y.c, replace=TRUE)

bart.boot <- bartlett.test(yboot,fact)
if(bart.boot[[1]] >= bart[[1]]){
   compt.boot=compt.boot+1}

}


Result <- matrix(0,1,4)
colnames(Result) <- c("Statistic", "Param.prob", "Permut.prob", "Bootstrap.prob")
rownames(Result) <- "Bartlett" 

Result[1,1] <- round(bart[[1]],4)
Result[1,2] <- round(bart[[3]],4)
Result[1,3] <- compt.perm/(nperm+1)
Result[1,4] <- compt.boot/(nperm+1)

Result

}
