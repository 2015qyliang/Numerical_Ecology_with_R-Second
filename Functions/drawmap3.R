# Map of the Doubs river (see Chapter 4)

drawmap3 <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", 
           colors = palette()[-1], pch = 21, tcol = "black") {

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
        pch = pch,
        cex = 3,
        col = "white",
        bg = colors[i]
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
      pch = 22,
      col = "white",
      pt.bg = colors,
      pt.cex = 2,
      bty = "n"
    )
    
  }