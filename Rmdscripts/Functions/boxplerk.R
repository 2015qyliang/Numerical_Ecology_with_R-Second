### Function to perform (unpaired) Kruskal-Wallis or paired Wilcoxon test with post-hoc
### multiple comparisons and boxplots with letters (using agricolae::kruskal)
### Francois Gillet, 15.10.2017 (adapted to agricolae 1.2-8)

# Arguments:
# X = numeric response variable
# Y = factor (groups, treatments...)
# p.adj = correction of p-values for multiple comparisons after Kruskal-Wallis test
# (default=none, bonferroni, holm...)
# paired = if TRUE and if 2 groups only, then Wilcoxon paired test is performed
# (objects are supposed to be ordered twice the same way in the data frame)

# Value:
# A summary table ($comparison) is printed and boxplots are drawn with result
# of the test ($p.value) and, if it is significant, of post-hoc tests as letters
# (in decreasing order)

# # Example:
# library(agricolae)
# source("boxplerk.R")
# library(stats)
# data(InsectSprays)
# boxplerk(
#   InsectSprays$count,
#   InsectSprays$spray,
#   ylab = "count",
#   xlab = "spray",
#   bcol = "bisque",
#   p.adj = "holm",
#   las = 1
# )


boxplerk <-
  function(X,
           Y,
           main = NULL,
           xlab = NULL,
           ylab = NULL,
           bcol = "bisque",
           p.adj = "none",
           cexy = 1,
           varwidth = TRUE,
           las = 1,
           paired = FALSE)
  {
    aa <- levels(as.factor(Y))
    an <- as.character(c(1:length(aa)))
    tt1 <- matrix(nrow = length(aa), ncol = 7)
    
    for (i in 1:length(aa))
    {
      temp <- X[Y == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- median(temp, na.rm = TRUE)
      tt1[i, 7] <- length(temp)
    }
    
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
    
    boxplot(
      X ~ Y,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth
    )
    require(agricolae)
    Yn <- factor(Y, labels = an)
    comp <- kruskal(X, Yn, p.adj = p.adj)
    sig <- "ns"
    
    if (paired == TRUE & length(aa) == 2)
    {
      coms <- wilcox.test(X ~ Yn, paired = TRUE)
      pp <- coms$p.value
    }
    else
    {
      pp <- comp$statistics$p.chisq
    }
    
    if (pp <= 0.1)
      sig <- "."
    if (pp <= 0.05)
      sig <- "*"
    if (pp <= 0.01)
      sig <- "**"
    if (pp <= 0.001)
      sig <- "***"
    
    gror <- comp$groups[order(rownames(comp$groups)), ]
    tt1$rank <- gror$X
    tt1$group <- gror$groups
    mtext(
      sig,
      side = 3,
      line = 0.5,
      adj = 0,
      cex = 2,
      font = 1
    )
    if (pp <= 0.1)
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 0.5,
        cex = 1,
        font = 4
      )
    
    list(comparison = tt1, p.value = pp)
    
  }
