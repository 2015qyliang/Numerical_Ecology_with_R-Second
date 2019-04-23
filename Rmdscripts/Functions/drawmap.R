# Map of the Doubs river (see Chapter 2)

drawmap <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", tcol = "white") {

    # Draw the Doubs river
    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "x coordinate (km)",
      ylab = "y coordinate (km)"
    )
    lines(xy, col = "light blue")
    text(65, 20, "Upstream", cex = 1.2)
    text(15, 32, "Downstream", cex = 1.2)
    
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
#        pch = i + 20,
        pch = i + 15,
        cex = 3,
        col = i + 1,
        bg = i + 1
      )
    }
    text(xy,
         row.names(xy),
         cex = 0.8,
         col = tcol,
         font = 2)
    legend(
      "bottomright",
      paste("Cluster", 1:k),
#      pch = (1:k) + 20,
      pch = (1:k) + 15,
      col = 2:(k + 1),
      pt.bg = 2:(k + 1),
      pt.cex = 2,
      bty = "n"
    )
    
  }