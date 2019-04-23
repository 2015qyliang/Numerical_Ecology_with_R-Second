plot.lda = function(lda.out, groups, colour.vec=NULL, plot.sites=1, plot.centroids=0, xax=1, yax=2, plot.env=TRUE, plot.ell=TRUE, title="LDA predicted classes", mul.coef=2, pos.names=NULL, col.env="black", xlim=NULL, ylim=NULL) 
### 
# lda.out : Output file of function lda() in {MASS}.
# groups  : Vector listing the group number for each object (factor or numeric).
# colour.vec : personal vector with colour names, to be used for the groups in
#    the plot instead of the standard colours. Example vector with 7 colours:
#    col.vec = c("gray60","bisque","brown4","red","blue","darkgreen","orange4")
#    Function colors() makes 657 colours available to R users.
# plot.sites: 0 = Do not plot the sites.
#             1 = plot symbols for the sites
#             2 = print the site names
# plot.centroids: 0 = Do not plot the group centroids.
#                 1 = Plot the group centroids (triangles) with assigned 
#                     group colours.
# xax : The axis to be used for abscissa of the plot.
# yax : The axis to be used for ordinate of the plot.
# plot.env=TRUE: plot the explanatory variables to the plot. FALSE: Do not 
#                plot them.
# plot.ell=TRUE: plot the 95% coverage ellipses of the groups. FALSE: Do not
#                plot them.
# title : Allows user to customize the title that will print above the plot.
# mul.coef : Multiplication factor for the length of the variable arrows. Some
#            trial and error is needed here.
# pos.names : Offset the names of the binary variables: 1=bottom, 2=left, 
#             3=top, 4=right.
# col.env="black" : Colour for the environmental variable arrows and names.
# xlim, ylim : Vectors describing the minimum and maximum values of the plotted region.
#
# Authors:: Daniel Borcard and Pierre Legendre
{
if(class(lda.out) != "lda") stop("File lda.out was not produced by function lda of MASS")
if(min(summary(as.factor(groups))) < 2) stop("There is at least one group with less than 2 observations")
library(ellipse)
coef = lda.out$scaling   # Standardized discriminant function coefficients
k <- ncol(coef)  # Number of canonical axes
lev <- length(levels(as.factor(groups)))  # Number of groups
# print(c(k,lev))
Fp <- predict(lda.out)$x
centre <- matrix(NA,lev,k)
for(i in 1:lev) { # Compute the group centroids in LDA space
#	centre[i,] <- apply(Fp[groups==i,], 2, mean) }  
	centre[i,] <- apply(Fp[groups==levels(as.factor(gr))[i],], 2, mean) }  
# print(centre)
class.num <- as.numeric(predict(lda.out)$class) # Assignment of sites to classes
# print(class.num)
if(xax > k) stop("Their are not enough canonical axes; change the xax value")
if(yax > k) stop("Their are not enough canonical axes; change the yax value")
xlab=paste("LDA axis",xax," ")
ylab=paste("LDA axis",yax," ")
plot(Fp[,xax], Fp[,yax], type="n", main=title, xlab=xlab, ylab=ylab, xlim, ylim, asp=1)
if(plot.sites==1) {   # Plot symbols for the sites
	if(is.null(colour.vec)) {
		points(Fp[,xax], Fp[,yax], pch=21, bg=class.num+1)
		} else {
		colour.sel <- colour.vec[class.num]
		points(Fp[,xax], Fp[,yax], pch=21, bg=colour.sel)
		}
	} else if(plot.sites==2) {
	if(is.null(colour.vec)) {
		text(Fp[,xax], Fp[,yax], row.names(Fp), col=class.num+1)
		} else {
		colour.sel <- colour.vec[class.num]
		text(Fp[,xax], Fp[,yax], row.names(Fp), col=colour.sel)
		}
	}	
if(plot.centroids) {
	if(is.null(colour.vec)) {
		points(centre[,xax], centre[,yax], pch=24, bg=(1:lev)+1, cex=1.5)
		} else {
		colour.sel <- colour.vec[1:lev]
		# print(colour.sel)
		points(centre[,xax], centre[,yax], pch=24, bg=colour.sel, cex=1.5)
		}
	}
abline(v=0, lty="dotted")
abline(h=0, lty="dotted")
# Draw 95% ellipses around the groups
if(plot.ell) {
	for(i in 1:length(levels(as.factor(groups)))) { 
#		cov <- cov(Fp[groups==i,c(xax,yax)])
		cov <- cov(Fp[groups==levels(as.factor(gr))[i],
		c(xax,yax)])
#		centre <- apply(Fp[groups==i,c(xax,yax)], 2, mean)
		centre <- apply(Fp[groups==levels(as.factor(gr))[i],
		c(xax,yax)], 2, mean)
		lines(ellipse(cov, centre=centre, level=0.95))
		}
	}
if(plot.env) { 
	arrows(x0=0, y0=0, x1=coef[,xax]*mul.coef, y1=coef[,yax]*mul.coef, 
		col=col.env, code=2, lty=1, length=0.1, lwd=1) 	
	if(!is.null(rownames(coef))) {
		text(x=coef[,xax]*mul.coef, y=coef[,yax]*mul.coef, 
			rownames(coef), col=col.env, cex=1, pos=pos.names) }
	}
}
