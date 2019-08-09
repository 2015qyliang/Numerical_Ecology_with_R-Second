'quickMEM' <- 
	function(Y, XY, stand=FALSE, thresh=NULL, method="fwd", myspat=NULL, alpha=0.05, rangexy=FALSE, detrend=TRUE, perm.max=999)
{

#              *** Quick exploratory dbMEM analysis ***
#                  Version 1.0.0 (21 july 2016)
#
#   This function replaces quickPCNM, which computed "classical" PCNM eigen-
#   vectors after Borcard and Legendre (2002). This new function, quickMEM, 
#   computes distance-based Moran's eigenvector maps (dbMEM) as described by 
#   Dray et al. (2006) and Legendre and Legendre (2012). The eigenvectors are 
#   the same, but dbMEM eigenvalues are directly proportional to Moran'I 
#   coefficient of spatial correlation of the corresponding eigenvectors. 
#   This makes it easy to identify the eigenvectors modeling positive spatial
#   correlation: they are the ones whose eigenvalue is larger than Moran's I
#   expectaion, which is E(I) = -1/(n-1) where n is the number of
#   observations (in practice, E(I) is close to 0 when n is large).
#
#   This function performs a dbMEM analysis for quick assessment of spatial 
#   structures, either on the basis of geographical coordinates (which can be 
#   one- or two-dimensional) or wich a set of pre-computed complex spatial 
#   variables (e.g. dbMEM with negative spatial corrrelation, MEM, AEM, 
#   polynomials...). 
#   If ONLY geographical coordinates are provided, the function automatically
#   computes dbMEM variables. Only dbMEM with positive eigenvalues (and,
#   therefore, Moran's I larger than its expectation) are computed.
#   Geographical coordinates must be provided even if complex spatial variables
#   are given as well.
#   If MEM or other complex variables are provided in the myspat argument, 
#   the function runs the spatial analysis directly on these variables.
#   Be careful NOT to provide enough spatial variables to saturate the model,
#   i.e., the number of spatial variables must be smaller than n-1. Otherwise
#   the global test will be nonsignificant and the analysis will be stopped.
#
#   The RDA is computed on a covariance matrix. If one wishes an RDA on
#   a correlation matrix instead, the dependent variables have to be
#   standardized to zero mean and unit variance with argument stand=TRUE.
#
#   Y:   a numeric vector, matrix or data frame containing response data 
#        ready for the analysis (i.e. pretransformed if necessary)
#   XY: a matrix containing the geographical (x or x-y) coordinates 
#        of the sites.
#   stand: should the response data be standardized to mean 0 and unit variance?
#       If TRUE then the RDAs are computed on a correlation matrix. Note that 
#       this prior standardization has the exact same effect as requesting an 
#       rda() on unstandardized response variables with argument scale=TRUE. 
#   thresh : user-provided truncation distance (for cases where the otherwise
#        automatically computed largest value of minimum spanning tree
#        is not appropriate).
#   method: specifies the method of selection of dbMEM variables: 
#     "none": no selection. All the spatial variables are used. Sometimes
#        useful when the user provides precomputed spatial variables.
#     "fwd": regression-based forward selection. Default, double stopping
#        criterion (alpha level and adjusted R2), strongly recommended choice.
#   myspat: an optional user=provided matrix containing precomputed spatial 
#        variables like MEM, AEM or polynomials.
#   alpha: level of significance of the tests (default = 0.05).
#   rangexy: if set to TRUE, rescales the spatial coordinates to a minimum of
#       0 on both axes and a maximum of 10 on the axis wih the largest range,
#       without altering the ratio between the x and y coordinates. If only
#       one coordinate is provided, rescales it into a range of [0;10].
#       Useful when the x-y coordinates are provided in ridiculously large or
#       small values.
#   detrend: if TRUE (default), prior ot the dbMEM analysis, the function tests 
#       if a linear trend is present and significant at the level given by 
#       "alpha", and, if yes, runs a linear detrending of the response data. 
#       It is useless to waste sine-shaped dbMEM variables to model linear
#       trends.
#   perm.max: possibility of limiting the maximum number of permutations in the
#       RDA tests.
#       
#   The function provides spatial plots of the several first canonical 
#        axes (thanks to Guillaume Blanchet for improving the plots!)
#        and a list of results containing the dbMEM variables, RDA results, 
#        a global test of the RDA and, if relevant, tests of the canonical axes.
#   Several summary diagnostics are displayed on-screen.
#
#   Call (examples): 
#
#   data(mite)
#   data(mite.xy)
#   mite.h <- decostand(mite,"hellinger")
#       A. Usual application with defaut settings (works well in most
#             situations):
#          exampleA <- quickMEM(mite.h, mite.xy)
#       B. Same as A but with precomputed third-order spatial polynomials,
#             a modified alpha level and rescaling of coordinates:
#          xy.poly3 <- poly(as.matrix(mite.xy), degree=3)
#          exampleB <- quickMEM(mite.h, mite.xy, myspat=xy.poly3, alpha=0.01,
#                      rangexy=TRUE)
#       C. Fit a plane (first degree polynomial) on the data:
#          exampleC <- quickMEM(mite.h, mite.xy, myspat=mite.xy, method="none",
#                      detrend=FALSE)
#
#   If you want to run the function several times to compare results, don't
#   forget to ask for a new graphical windows. Otherwise the active window
#   created during the first run will be overwritten.
#   When the run is completed type 'summary(name_of_object)' to get the RDA 
#   results.
#
#
#                                        Daniel Borcard
#                                        Universite de Montreal
#                                        July 2016
#                                        Licence: GPL-2

# require(adegraphics)
require(adespatial)
require(vegan)

a <- system.time({

Y <- as.matrix(Y)
   if(stand) Y <- scale(Y)
XY <- as.matrix(XY)
n <- nrow(Y)
epsilon <- sqrt(.Machine$double.eps)
if(!is.null(thresh)) thresh <- thresh+epsilon

# ----------------------------------------------------------------------------

if(rangexy==TRUE) {

if(ncol(XY)==1){
XY <- (XY-min(XY))/((max(XY)-min(XY))*0.1)
               }
else{
mini <- apply(XY,2,"min")
maxi <- apply(XY,2,"max")
xy.trans <- sweep(XY,2,mini)
range.max <- 0.1*(max((maxi[1]-mini[1]),(maxi[2]-mini[2])))
XY <- as.matrix(xy.trans/range.max)
    }
                  }

# ----------------------------------------------------------------------------

if (is.null(myspat)) {

### Building the dbMEM variables
#   Using function dbmem {adespatial}

# Problem if user-provided threshold is too small

if(!is.null(thresh)){
thresh.min <- give.thresh(dist(XY))
if(thresh < thresh.min+epsilon) {
cat("\n ------------------------------------------------------------------")
cat("\n User-provided truncation threshold, 'thresh', too small.") 
cat("\n Discontinuous groups would be created.")
cat("\n The smallest value holding all observations together,", thresh.min,",")
cat("\n will be used instead.")
cat("\n ------------------------------------------------------------------\n")
thresh <- thresh.min
}
}

dbMEM.comput <- dbmem(XY, thresh=thresh, MEM.autocor="positive", silent=FALSE)
dbMEMbase <- as.data.frame(dbMEM.comput)

# Compensation for internal loss of object by R if necessary
assign("dbMEMbase", dbMEMbase, envir=.GlobalEnv)
                     }
                     
else {
dbMEMbase <<- as.data.frame(myspat)
dmin <- "implicit in file of spatial variables"
meanMEM <- sum(apply(dbMEMbase,2,mean))
if(abs(meanMEM) > epsilon ) {
cat("\n ------------------------------------------------------------------")
cat("\n WARNING: the user-provided spatial variables are not centred")
cat("\n to zero mean. Are you sure that they are correct?")
cat("\n ------------------------------------------------------------------")
                            }
sumcorMEM <- sum(cor(dbMEMbase))
if(abs(sumcorMEM) > ncol(dbMEMbase)+epsilon ) {
cat("\n ------------------------------------------------------------------")
cat("\n WARNING: the user-provided spatial variables are not orthogonal")
cat("\n to one another. Are you sure that they are correct?")
cat("\n ------------------------------------------------------------------")
cat("\n")
                                              }
     }
nb.ev <- ncol(dbMEMbase)
ev <- attributes(dbMEMbase)$values
if(nb.ev >= n-1)  stop(" Too many explanatory variables. They should be less than n-1", call.=FALSE)

# ----------------------------------------------------------------------------

### RDA "response matrix x dbMEM"

## Preliminary step: linear detrending of response data if trend is significant
if(detrend==TRUE){
   trace.correl <- sum(diag(cor(Y,XY)))
   if(abs(trace.correl) < epsilon){
     Y.det <- Y
     # Compensation for internal loss of object by R if necessary
     assign("Y.det", Y.det, envir=.GlobalEnv)  
     temp2.test <- matrix(rep(1,5),1)
                                  }
   else{  
   temp2 <- rda(Y,XY)
   temp2.test <- anova(temp2,permutations=how(nperm=perm.max))
      if(temp2.test[1, ncol(temp2.test)] <= alpha) {
         Y.det <- resid(lm(Y~XY))
         # Compensation for internal loss of object by R if necessary
         assign("Y.det", Y.det, envir=.GlobalEnv)   
       }                 
      else {
         Y.det <- Y
         # Compensation for internal loss of object by R if necessary
         assign("Y.det", Y.det, envir=.GlobalEnv)   
   temp2.test <- matrix(rep(1,5),1)
           }
       }
                 }
else {
cat("\n User asked for no detrending\n\n")
   Y.det <- Y
     # Compensation for internal loss of object by R if necessary
     assign("Y.det", Y.det, envir=.GlobalEnv)   
     }

## RDA with all dbMEM variables(complete model)

mod1 <- rda(Y.det~.,data=dbMEMbase)
global.test <- anova(mod1,permutations=how(nperm=perm.max))
mod1.sum <- summary(mod1,scaling=1)
R2glob <- RsquareAdj(mod1)[[1]]
R2glob.a <- RsquareAdj(mod1)[[2]]

if(global.test[1, ncol(global.test)] >= alpha) {
   cat("\n ------------------------------------------------------------------")
   cat("\n *** Procedure stopped ***")
   cat("\n p-value of global test: ",global.test[1,ncol(global.test)])
   cat("\n No significant spatial structure detected by global dbMEM analysis.")
   cat("\n Selection of dbMEM variables would lead to spurious model.")
   cat("\n ------------------------------------------------------------------","\n")
   stop(call.=FALSE)
                                               }

# ----------------------------------------------------------------------------

else {

METHODS <- c("none", "fwd")
    method <- match.arg(method, METHODS)

if(method == "none"){ 
mod <- mod1
mod.test <- global.test
mod.sum <- mod1.sum
R2glob <- R2glob
R2glob.a <- R2glob.a
vars.sign <- c(1:ncol(dbMEMbase))
nb.sig.ev <- length(vars.sign)
                    }

else{

  if(method == "fwd"){   # Open  fwd

  ## Regression-based forward selection of dbMEM variables and RDA on 
  ## significant dbMEM variables
   
  # If there is only one response variable that is normally distributed, save
  # time by replacing permutation tests by parametric tests. Otherwise and for
  # a multivariate response matrix, forward selection with the adjusted R2 as 
  # aditional stopping criterion.
    if(ncol(Y)==1){
    norm.test <- shapiro.test(Y.det)
       if(norm.test[2] > 0.05) {
cat("\n ------------------------------------------------------------------")
cat("\n Only one, normally distributed response variable found. Parametric")
cat("\n forward selection on standardized response variable is run.")
cat("\n ------------------------------------------------------------------","\n")
       Y.det <-scale(Y.det)
       fwd.sel <- forward.sel.par(Y.det,dbMEMbase,alpha=alpha,adjR2thresh=R2glob.a)
                               } 
	   else {
cat("\n ------------------------------------------------------------------")
cat("\n The only response variable is not normally distributed.")
cat("\n Permutational forward selection is run.")
cat("\n ------------------------------------------------------------------","\n")
       fwd.sel <- forward.sel(Y.det,dbMEMbase,alpha=alpha,adjR2thresh=R2glob.a)
            }
                  } 
	else{
  fwd.sel <- forward.sel(Y.det,dbMEMbase,alpha=alpha,adjR2thresh=R2glob.a)
        }

  nb.sig.ev <- nrow(fwd.sel)
  vars.sign <- sort(fwd.sel[,2])
                     }   # 1.2.1 close  fwd
    }

dbMEMred <- as.data.frame(dbMEMbase[,c(vars.sign)])
# compensation for internal loss of object by R if necessary
assign("dbMEMred",dbMEMred,envir=.GlobalEnv)  
mod <- rda(Y.det~.,data=dbMEMred)
mod.sum <- summary(mod,scaling=1)
mod.test <- anova(mod)

R2 <- RsquareAdj(mod)[[1]]
R2adj <- RsquareAdj(mod)[[2]]

     }

# At this point, this function looses the track of one or the other object 
# defined higher. The workaround is to export several internal objects from the 
# function during the run and assign it to the global R environment (e.g. 
# Line 281). these objects are thus present in the main R workXY and should 
# be deleted prior to another quickMEM run.

if(ncol(Y)==1 || nb.sig.ev == 1) {
   nb.ax <- 1                    }
 else {
   mod.axes.test <- anova(mod,by="axis",cutoff=0.10)
   ## Count the significant axes 
   nb.ax <- length(which(mod.axes.test[,ncol(mod.axes.test)]<=alpha))

      }         
# Compensation for internal loss of object by R if necessary
# assign("nb.ax",nb.ax,envir=.GlobalEnv)  

# ----------------------------------------------------------------------------

## Plot of significant axes
fitted.scores <- scores(mod,display="lc",choices=1:nb.ax)
# Compensation for internal loss of object by R if necessary
assign("fitted.scores",fitted.scores,envir=.GlobalEnv)  
 par(mfrow=c(round(sqrt(nb.ax)),ceiling(sqrt(nb.ax))))

# Works but doomed (when s.value disappears from ade4):
# if(ncol(XY)==2){
#   for(i in 1:nb.ax){
#     ade4::s.value(XY,fitted.scores[,i],sub=paste("Axis ",i))
#                    }
#                } 
               
# Instead of s.value, call a slightly modified version (at the bottom 
# of this function)
               
 if(ncol(XY)==2){
  for(i in 1:nb.ax){
    sr.value(XY,fitted.scores[,i],sub=paste("Axis ",i))
                   }
               }

# Works but arranges the graphs automatically in an unaesthetic fashion
# in the graphical window:
# Compensation for internal loss of object by R if necessary
# assign("XY",XY,envir=.GlobalEnv)  
#   if(ncol(XY)==2){
#   s.value(XY, fitted.scores[,1:nb.ax])
#                  }

  else {
    for(i in 1:nb.ax){
    plot(XY,fitted.scores[,i],type="l",ylab="Fitted site scores")
                     }
      }

# ----------------------------------------------------------------------------

## Screen output
if(detrend==TRUE){
   if(temp2.test[1, ncol(temp2.test)] <= alpha) {
      cat("\n -------------------------------------------------------")
      cat("\n A significant linear trend has been found in the response data.")
      cat("\n The response data have been detrended prior to dbMEM analysis.")
                                                }
    else{
      cat("\n -------------------------------------------------------")
      cat("\n No significant linear trend has been found in the response data.")
      cat("\n The data have NOT been detrended prior to dbMEM analysis.")     
                                                }
                 }


cat("\n-------------------------------------------------------\n")
cat(" ",nb.ev,"dbMEM eigenvectors have been produced","\n")
cat("  R2 of global model = ",round(R2glob,4),"\n")
cat("  Adjusted R2 of global model = ",round(R2glob.a,4),"\n")
if(method != "none") {
if(nb.sig.ev==1){
cat(" ",nb.sig.ev,"dbMEM eigenvector has been selected","\n")
   if(nb.ax==1){
   cat("  The final model has",nb.ax,"significant canonical axis","\n")
               }
   else{
   cat("  The final model has",nb.ax,"significant canonical axes","\n")
       }
                }
else{
cat(" ", nb.sig.ev," dbMEM eigenvectors have been selected","\n")}
cat("  R2 of minimum (final) model = ",round(R2,4),"                    ","\n")
cat("  Adjusted R2 of minimum (final) model = ",round(R2adj,4),"        ","\n")
   if(nb.ax==1){
   cat("  The final model has",nb.ax,"significant canonical axis","\n")
               }
   else{
   cat("  The final model has",nb.ax,"significant canonical axes","\n")
       }
                     }
cat("---------------------------------------------------------")
cat("\n")

})
a[3] <- sprintf("%2f",a[3])
cat(" Time to compute quickMEM =",a[3]," sec",'\n')

## Extraction of results to be returned

if(method == "none") {
   if(ncol(Y)>1 && nb.sig.ev > 1) {
      table <- list(dbMEMbase, ev[1:nb.ev], mod, mod.test, mod.axes.test)
      names(table) <- c("dbMEM","eigenvalues","RDA","RDA_test","RDA_axes_tests")
                                  }
   else {
      table <- list(dbMEMbase,ev[1:nb.ev],mod,mod.test)
      names(table) <- c("dbMEM","eigenvalues","RDA","RDA_test")   
        }
                     } 
else {
   if(ncol(Y)>1 && nb.sig.ev > 1) {
      table <- list(dbMEMbase, ev[1:nb.ev], fwd.sel, dbMEMred, mod, mod.test, mod.axes.test)
      names(table) <- c("dbMEM","eigenvalues","fwd.sel","dbMEM_red_model","RDA","RDA_test","RDA_axes_tests")
                                  } 
   else {
      table <- list(dbMEMbase,ev[1:nb.ev], fwd.sel, dbMEMred,mod, mod.test)
      names(table) <- c("dbMEM","eigenvalues","fwd.sel","dbMEM_red_model","RDA","RDA_test")
        }
     }
return(table)
}



sr.value <- function (dfxy, z, xax = 1, yax = 2, method = c("bubble",
	"greylevel"), zmax = NULL, csize = 1, cpoint = 0, pch = 20,
	clegend = 0.75, neig = NULL, cneig = 1, xlim = NULL, ylim = NULL,
	grid = TRUE, addaxes = TRUE, cgrid = 0.75, include.origin = TRUE,
	origin = c(0, 0), sub = "", csub = 1, possub = "topleft",
	pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE)
#
# Slightly modified version of ade4's s.value() graphical function.
# Draws round instead of square bubbles in some plots when argument 
# "bubble" is called.
#
# License: GPL-2
# Author of the original function s.value: Daniel Chessel
# Modification: Francois Gillet, 25 August 2012
#
{
	dfxy <- data.frame(dfxy)
	if (length(z) != nrow(dfxy))
		stop(paste("Non equal row numbers", nrow(dfxy), length(z)))
	opar <- par(mar = par("mar"))
	on.exit(par(opar))
	par(mar = c(0.1, 0.1, 0.1, 0.1))
	coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax,
		xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
		cgrid = cgrid, include.origin = include.origin, origin = origin,
		sub = sub, csub = csub, possub = possub, pixmap = pixmap,
		contour = contour, area = area, add.plot = add.plot)
	if (!is.null(neig))
	{
		if (is.null(class(neig))) neig <- NULL
		if (class(neig) != "neig") neig <- NULL
		deg <- attr(neig, "degrees")
		if (length(deg) != length(coo$x)) neig <- NULL
	}
	if (!is.null(neig))
	{
		fun <- function(x, coo)
		{
			segments(coo$x[x[1]], coo$y[x[1]], coo$x[x[2]], coo$y[x[2]],
				lwd = par("lwd") * cneig)
		}
		apply(unclass(neig), 1, fun, coo = coo)
	}
	method <- method[1]
	if (method == "greylevel")
	{
		br0 <- pretty(z, 6)
		nborn <- length(br0)
		coeff <- diff(par("usr")[1:2])/15
		numclass <- cut.default(z, br0, include = TRUE, lab = FALSE)
		valgris <- seq(1, 0, le = (nborn - 1))
		h <- csize * coeff
		for (i in 1:(nrow(dfxy)))
			{
				symbols(coo$x[i], coo$y[i], circles = h/2, 
					bg = gray(valgris[numclass[i]]),
					add = TRUE, inch = FALSE)
			}
		scatterutil.legend.circle.grey(br0, valgris, h/2, clegend)
		if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
	}
	else if (method == "bubble")
	{
		coeff <- diff(par("usr")[1:2])/15
		sq <- sqrt(abs(z))
		if (is.null(zmax)) zmax <- max(abs(z))
		w1 <- sqrt(zmax)
		sq <- csize * coeff * sq/w1
		for (i in 1:(nrow(dfxy)))
		{
			if (sign(z[i]) >= 0)
			{
				symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "black", 
					fg = "white", add = TRUE, inch = FALSE)
			}
			else
			{
				symbols(coo$x[i], coo$y[i], circles = sq[i]/2, bg = "white", 
					fg = "black", add = TRUE, inch = FALSE)
			}
		}
		br0 <- pretty(z, 4)
		l0 <- length(br0)
		br0 <- (br0[1:(l0 - 1)] + br0[2:l0])/2
		sq0 <- sqrt(abs(br0))
		sq0 <- csize * coeff * sq0/w1
		sig0 <- sign(br0)
		if (clegend > 0) scatterutil.legend.bw.circle(br0, sq0, sig0, clegend)
		if (cpoint > 0) points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
	}
	else if (method == "circlesize") print("not yet implemented")
	if (!add.plot) box()
	invisible(match.call())
}



scatterutil.legend.bw.circle <- function (br0, sq0, sig0, clegend)
{
	br0 <- round(br0, dig = 6)
	cha <- as.character(br0[1])
	for (i in (2:(length(br0)))) cha <- paste(cha, br0[i], sep = " ")
	cex0 <- par("cex") * clegend
	yh <- max(c(strheight(cha, cex = cex0), sq0))
	h <- strheight(cha, cex = cex0)
	y0 <- par("usr")[3] + yh/2 + h/2
	ltot <- strwidth(cha, cex = cex0) + sum(sq0) + h
	rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, 
		par("usr")[1] + ltot + h/4, y0 + yh/2 + h/4, col = "white")
	x0 <- par("usr")[1] + h/2
	for (i in (1:(length(sq0))))
	{
		cha <- br0[i]
		cha <- paste(" ", cha, sep = "")
		xh <- strwidth(cha, cex = cex0)
		text(x0 + xh/2, y0, cha, cex = cex0)
		z0 <- sq0[i]
		x0 <- x0 + xh + z0/2
		if (sig0[i] >= 0)
			symbols(x0, y0, circles = z0/2, bg = "black", fg = "white",
				add = TRUE, inch = FALSE)
		else symbols(x0, y0, circles = z0/2, bg = "white", fg = "black",
			add = TRUE, inch = FALSE)
		x0 <- x0 + z0/2
	}
	invisible()
}



scatterutil.legend.circle.grey <- function (br0, valgris, h, clegend)
{
	if (clegend <= 0) return(invisible())
	br0 <- round(br0, dig = 6)
	nborn <- length(br0)
	cex0 <- par("cex") * clegend
	x0 <- par("usr")[1] + h
	x1 <- x0
	for (i in (2:(nborn)))
	{
		x1 <- x1 + h
		cha <- br0[i]
		cha <- paste(cha, "]", sep = "")
		xh <- strwidth(cha, cex = cex0)
		if (i == (nborn)) break
		x1 <- x1 + xh + h
	}
	yh <- max(strheight(paste(br0), cex = cex0), h)
	y0 <- par("usr")[3] + yh/2 + h/2
	rect(par("usr")[1] + h/4, y0 - yh/2 - h/4, x1 - h/4, y0 + yh/2 + h/4, 
		col = "white")
	x0 <- par("usr")[1] + h
	for (i in (2:(nborn)))
	{
		symbols(x0, y0, circles = h/2, bg = gray(valgris[i - 1]), add = TRUE, 
			inch = FALSE)
		x0 <- x0 + h
		cha <- br0[i]
		if (cha < 1e-05) cha <- round(cha, dig = 3)
		cha <- paste(cha, "]", sep = "")
		xh <- strwidth(cha, cex = cex0)
		if (i == (nborn)) break
		text(x0 + xh/2, y0, cha, cex = cex0)
		x0 <- x0 + xh + h
	}
	invisible()
}
