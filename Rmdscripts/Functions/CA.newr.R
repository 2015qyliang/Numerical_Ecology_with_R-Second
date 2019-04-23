`CA.newr` <- 
   function(Y, use.svd=TRUE, cumfit.obj=TRUE, cumfit.var=TRUE)
#
# Compute correspondence analysis (CA).
# Data table Y must contain frequencies or equivalent.
#
# use.svd: decomposition is done by svd (default). It can also be done by eigen.
#          The signs of axes may differ between the two methods.
# color  : color of the species symbols and labels in the biplots.
#
# For exercise, you may choose 'Table_9.11.txt' (small example from Chapter 9
# of the manual) or 'Spiders_28x12_spe.txt' (larger data set, real data).
#
#           Pierre Legendre, Universite de Montreal, January 2008
{
# Begin internal functions
sq.length <- function(vec) sum(vec^2)
#
'cumul.fit.va' <- function(Fhat,n,p,k,var.names)
	# Compute the table of "Cumulative fit per variable" 
{
	sp.var <- diag(var(Y))
	res <- matrix(NA,p,k)
	for(i in 1:p) {
	   res[i,] <- cumsum(Fhat[i,]^2)/sq.length(Fhat[i,])
	   }
	rownames(res) <- var.names
	colnames(res) <- paste("Cum.axis",1:k,sep=".")
	res
}
#
'cumul.fit.ob' <- function(F,n,p,k,obj.names)
	# Compute the table of "Cumulative fit of the objects" 
{
	res <- matrix(NA,n,k)
	for(i in 1:n) {
	   res[i,] <- cumsum(F[i,]^2)/sq.length(F[i,])
	   }
	rownames(res) <- obj.names
	colnames(res) <- paste("Cum.axis",1:k,sep=".")
	res
}
# End internal functions
#
Y = as.matrix(Y)
if(min(Y) < 0) stop("Negative values not allowed in CA")
#
# Calculate basic parameters of Y
n = nrow(Y)
p = ncol(Y)
#
# Save the row and column names
site.names = rownames(Y)
sp.names = colnames(Y)
#
# Construct the Qbar matrix (contributions to chi-square)
# Numerical ecology (1998), equations 9.31 and 9.32
fi. = matrix(apply(Y,1,sum),n,1)
f.j = matrix(apply(Y,2,sum),1,p)
f.  = sum(fi.)
pi. = as.vector(fi./f.)
p.j = as.vector(f.j/f.)
E = (fi. %*% f.j)/f.
Qbar = (Y - E) * E^(-0.5) / sqrt(f.)
inertia = sum(Qbar^2)
#
if(use.svd) {
   # Analyse Qbar by 'svd'
   svd.res = svd(Qbar)
   k = length(which(svd.res$d > 1e-8))
   values = svd.res$d[1:k]^2
   U = svd.res$v[,1:k]
   Uhat = svd.res$u[,1:k]
   } else {
   # Alternative analysis or Qbar by 'eigen'
   Qbar = as.matrix(Qbar)
   QprQ.eig = eigen( t(Qbar) %*% Qbar )
   k = length(which(QprQ.eig$values > 1e-16))
   values = QprQ.eig$values[1:k]
   U = QprQ.eig$vectors[,1:k]
   Uhat = Qbar %*% U %*% diag(values^(-0.5))
   }
#
rel.values = values/inertia
cum.rel <- cumsum(rel.values)
#
# Construct matrices V, Vhat, F, and Fhat for biplots, scalings 1 and 2
V = diag(p.j^(-0.5)) %*% U
Vhat = diag(pi.^(-0.5)) %*% Uhat
F = Vhat %*% diag(values^(0.5))
Fhat = V %*% diag(values^(0.5))
#
# Matrices for biplot, scaling = 3 (Symmetric scaling in Canoco)
spec3 = V %*% diag(values^(0.25))                # Species scores
site3 = Vhat %*% diag(values^(0.25))             # Site scores
#
   if(cumfit.var) {
      cfit.spe <- cumul.fit.va(Fhat,n,p,k,sp.names)
      } else {
      cfit.spe <- NULL
      }
#
   if(cumfit.obj) {
      cfit.obj <- cumul.fit.ob(F,n,p,k,site.names)
      } else {
      cfit.obj <- NULL
      }
#
rownames(U) <- rownames(V) <- rownames(spec3) <- rownames(Fhat) <- sp.names
rownames(Uhat) <- rownames(F) <- rownames(Vhat) <- rownames(site3) <- site.names
ax.names <- paste("Axis",1:k,sep="")
colnames(U) <- colnames(Uhat) <- colnames(V) <- colnames(spec3) <- colnames(Vhat) <- colnames(site3) <- colnames(F) <- colnames(Fhat) <- ax.names
#
general <- list(inertia=inertia, values=values, rel.values=rel.values, cum.rel=cum.rel)
scaling1 <- list(species=V, sites=F)
scaling2 <- list(species=Fhat, sites=Vhat)
scaling3 <- list(species=spec3, sites=site3)
scaling4 <- list(species=Fhat, sites=F)
fit <- list(cumulfit.spe=cfit.spe, cumulfit.obj=cfit.obj)
other <- list(U=U, Uhat=Uhat, F=F, Fhat=Fhat, site.names=site.names, sp.names=sp.names, Qbar=Qbar, call=match.call() )
#
out <- list(general=general, scaling1=scaling1, scaling2=scaling2, scaling3=scaling3, scaling4=scaling4, fit=fit, other=other)
class(out) <- "CA.newr"
out
}

`print.CA` <-
    function(x, kk=5, ...)
{
if (!inherits(x, "CA.newr")) stop("Object of class 'CA.newr' expected")
    cat("\nCorrespondence Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$other$call),'\n')
    cat("\nTotal inertia in matrix Qbar: ",x$general$inertia,'\n')
    cat("\nEigenvalues",'\n')
    cat(x$general$values,'\n')
    cat("\nRelative eigenvalues",'\n')
    cat(x$general$rel.values,'\n')
    cat("\nCumulative relative eigenvalues",'\n')
    cat(x$general$cum.rel,'\n')
    kk <- min(length(x$general$values), kk)
    if(!is.null(x$fit$cumulfit.spe)) {
       cat("\nCumulative fit per species (",kk,"axes)",'\n')
       print.default(x$fit$cumulfit.spe[,1:kk], digits=5)
       }
    if(!is.null(x$fit$cumulfit.obj)) {
       cat("\nCumulative fit of the objects (",kk,"axes)",'\n')
       print.default(x$fit$cumulfit.obj[,1:kk], digits=5)
       }
    cat('\n')
    invisible(x) 
}

`biplot.CA` <-
function(x, xax=1, yax=2, scaling=1, aspect=1, cex=1, color.sites="black", color.sp="red",...)
# xax and yax determine the axes that will be plotted.
# Use aspect=NA to remove the effect of parameter 'asp' in the biplot.
{
if (!inherits(x, "CA.newr")) stop("Object of class 'CA.newr' expected")
if(length(x$general$values) < 2) stop("There is a single eigenvalue. No plot can be produced.")
#
sp.names = x$other$sp.names
si.names = x$other$site.names
#

if(scaling == 1) {

# The sites are at the centroids (barycentres) of the species
# This projection preserves the chi-square distance among the sites
type = "scaling type 1"
sp = x$scaling1$species
si = x$scaling1$sites

} else if(scaling == 2) {

# The species are at the centroids (barycentres) of the sites
# This projection preserves the chi-square distance among the species
type = "scaling type 2"
sp = x$scaling2$species
si = x$scaling2$sites

} else if(scaling == 3) {

# Biplot, scaling = 3 (Symmetric scaling in Canoco)
type = "scaling type 3"
sp = x$scaling3$species
si = x$scaling3$sites

} else if(scaling == 4) {

# For contingency tables --
# Preserves the chi-square distance among the rows and among the columns
type = "scaling type 4"
sp = x$scaling2$species
si = x$scaling1$sites

} else { 

stop("Program stopped: error in scaling type")
}

# Find the limits of the axes
sp.range = apply(sp[,c(xax,yax)],2,range)
si.range = apply(si[,c(xax,yax)],2,range)

# Biplot: plot 'si' for sites, 'sp' for species
ran.si = si.range[2,] - si.range[1,]
ran.sp = sp.range[2,] - sp.range[1,]
ran.x = max(ran.si[1], ran.sp[1])
xmin = min(sp.range[1,1], si.range[1,1]) - ran.x/8
xmax = max(sp.range[2,1], si.range[2,1]) + ran.x/3
ymin = min(sp.range[1,2], si.range[1,2])
ymax = max(sp.range[2,2], si.range[2,2])
#
plot(si[,c(xax,yax)], asp=aspect, pch=20, cex=cex, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=paste("CA axis",xax), ylab=paste("CA axis",yax), col=color.sites)
text(si[,c(xax,yax)], labels=si.names, cex=cex, pos=4, offset=0.5, col=color.sites)
points(sp[,c(xax,yax)], pch=22, cex=cex, col=color.sp)
text(sp[,c(xax,yax)], labels=sp.names, cex=cex, pos=4, offset=0.5, col=color.sp)
title(main = c("CA biplot",type), family="serif")
#
invisible()
}
