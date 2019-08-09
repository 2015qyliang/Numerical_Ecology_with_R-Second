# # coldiss()
# # Color plots of a dissimilarity matrix, without and with ordering
# #
# # License: GPL-2
# # Author: Francois Gillet, 23 August 2012 - rev. 07 June 2016
# #
# 
# coldiss <- function(D,
#                       nc = 4,
#                       byrank = TRUE,
#                       diag = FALSE) {
#   require(gclus)
#   
#   D <- as.dist(as.matrix(D))
#   
#   if (max(D) > 1)
#     D <- D / max(D)
#   
#   if (byrank) {
#     spe.color <- dmat.color(1 - D, cm.colors(nc))
#   }
#   else {
#     spe.color <- dmat.color(1 - D, byrank = FALSE, cm.colors(nc))
#   }
#   
#   spe.o <- order.single(1 - D)
#   speo.color <- spe.color[spe.o, spe.o]
#   
#   op <- par(mfrow = c(1, 2), pty = "s")
#   
#   if (diag) {
#     plotcolors(
#       spe.color,
#       rlabels = attributes(D)$Labels,
#       main = "Dissimilarity Matrix",
#       dlabels = attributes(D)$Labels
#     )
#     plotcolors(
#       speo.color,
#       rlabels = attributes(D)$Labels[spe.o],
#       main = "Ordered Dissimilarity Matrix",
#       dlabels = attributes(D)$Labels[spe.o]
#     )
#   }
#   else {
#     plotcolors(spe.color, rlabels = attributes(D)$Labels,
#                main = "Dissimilarity Matrix")
#     plotcolors(speo.color,
#                rlabels = attributes(D)$Labels[spe.o],
#                main = "Ordered Dissimilarity Matrix")
#   }
#   
#   par(op)
# }
# 
# # Usage:
# # coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# # If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# # nc 							number of colours (classes)
# # byrank = TRUE		equal-sized classes
# # byrank = FALSE	equal-length intervals
# # diag = TRUE			print object labels also on the diagonal
# 
# # Example:
# # coldiss(spe.dj, nc = 9, byrank = FALSE, diag = TRUE)


#' Color plots of a dissimilarity matrix, without and with ordering
#'
#' This function generates color plots of dissimilarity matrices, with our
#' without ordering. It was supplied as supplementary data to "Numerical
#' ecology with R" Boccard & Legendre, 2012. And authored by Francois Gillet.
#'
#' @param D dissimilarity matrix. If D is not a dissimilarity matrix (max(D)>1),
#' then D is divided by max(D)
#' @param nc number of colors (classes, defaults to 4)
#' @param byrank specifiy either equal-sized classes (TRUE) or equal-length
#' intervals (FALSE). Defaults to TRUE
#' @param diag print object labels on the diagonal (defaults to FALSE)
#' @importFrom gclus dmat.color plotcolors order.single
#' @importFrom graphics par
#' @importFrom grDevices cm.colors
#' @keywords heatmap
#' @examples
#' ## Short example
#'
#' # coldiss(spe.dj, nc=9, byrank=F, diag=T)
#' # Original author: Francois Gillet, 23 aug 2012
#'
#' @export


coldiss <- function(D, nc = 4, byrank = TRUE, diag = FALSE)
{
  #require(gclus)
  
  if (max(D)>1) D <- D/max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1-D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1-D, byrank=FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1-D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow=c(1,2), pty="s")
  
  if (diag) {
    plotcolors(spe.color, rlabels=attributes(D)$Labels,
               main="Dissimilarity Matrix",
               dlabels=attributes(D)$Labels)
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o],
               main="Ordered Dissimilarity Matrix",
               dlabels=attributes(D)$Labels[spe.o])
  }
  else {
    plotcolors(spe.color, rlabels=attributes(D)$Labels,
               main="Dissimilarity Matrix")
    plotcolors(speo.color, rlabels=attributes(D)$Labels[spe.o],
               main="Ordered Dissimilarity Matrix")
  }
  par(op)
}