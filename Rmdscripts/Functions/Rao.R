###################################################################################################################################
# 	The Rao function computes alpha, gamma and beta-components for taxonomic, functional and phylogenetic diversity with the Rao index
# 	The script integrates two functions: "Qdecomp", by Villeger & Mouillot (J Ecol, 2008) modified by Wilfried Thuiller, and "disc", by S. Pavoine, in the package ade4.
# 	For a regional assemblage of C local communities gamma = mean(alpha) + beta, where:
#  	gamma is the diversity of the regional pool
#  	alpha are the diversities of the local communities
#  	beta is the turn over between local communities
#  	diversity is estimated with the Rao quadratic entropy index (Rao 1982)
#
# INPUTS:
#	- "abundances": matrix of abundances (c x s) of the s species for the c local communities (or samples)
#	- "dfunct": matrix (s x s) or dist object with pairwise functional trait distances between the s species
#	- "dphyl": as dfunct but for phylogenetic distances
#	- "weight": defining if the correction by Villeger & Mouillot (J Ecol, 2008) is applied or not
#	- "Jost": defining if the correction Jost correction is applied (this paper and Jost 2007)
#	- "structure": a data frame containing the name of the group to which samples belong see
#      NA are not allowed in 'functdist'
#      NA are automatically replaced by 0 in 'abundances'
#
# OUTPUTS:
#	- The results are organized for Taxonomic diversity ($TD), Functional diversity ($FD) and phylogenetical diversity ($PD). Beta and gamma diversities are calculated for the whole data set and for each pair of samples ("Pairwise_samples")
#	- "$Richness_per_plot"(number of species per sample)
#	- "$Relative_abundance" (species relative abundances per plot)
#	- "$Pi" (species regional relative abundance)
#	- "$Wc" (weigthing factor),
#	- "$Mean_Alpha" (mean aplpha diversity; for taxonomic diversity the Simpson index is calculated)
#	- "$Alpha" (alpha diversity for each sample; for taxonomic diversity the Simpson index is calculated)
#	- "$Gamma" (gamma diversity; for taxonomic diversity the Simpson index is calculated)
#	- "$Beta_add" (Gamma-Mean_Alpha)
#	- "$Beta_prop" (Beta_add*100/Gamma)
#	- "$Pairwise_samples$Alpha" (mean alpha for each pair of samples)
#	- "$Pairwise_samples$Gamma" (gamma for each pair of samples)
#	- "$Pairwise_samples$Beta_add" (beta for each pair of samples as Gamma-Mean_Alpha)
#	- "$Pairwise_samples$Beta_prop" (beta for each pair of samples as Beta_add*100/Gamma)

#####################################################################################################################################



Rao <-
  function(sample,
           dfunc,
           dphyl,
           weight = F,
           Jost = F,
           structure = NULL)   {
    library(ade4)
    
    ####function Qdecomp by Villeger & Mouillot (J Ecol, 2008) modified by Wilfried Thuiller #####
    
    Qdecomp = function(functdist, abundances, w = TRUE) {
      # number and names of local communities
      c <- dim(abundances)[1]
      namescomm <- row.names(abundances)
      abundances <- as.matrix(abundances)
      
      # if necessary, transformation of functdist into matrix object
      if (is.matrix(functdist) == F)
        functdist <- as.matrix(functdist)
      
      # checking 'abundances' and 'functdist' dimensions
      if (dim(functdist)[1] != dim(functdist)[2])
        stop("error : 'functdist' has different number of rows and columns")
      if (dim(abundances)[2] != dim(functdist)[1])
        stop("error : different number of species in 'functdist' and 'abundances' ")
      
      # checking NA absence in 'functdist'
      if (length(which(is.na(functdist) == T)) != 0)
        stop("error : NA in 'functdist'")
      
      # replacement of NA by 0 in abundances
      if (is.na(sum(abundances)) == T)  {
        for (i in 1:dim(abundances)[1])
          for (j in 1:dim(abundances)[2])
          {
            if (is.na(abundances[i, j]) == T)
              abundances[i, j] <- 0
          } # end of i j
      } # end of if
      
      #  species richness and total abundances in local communities
      abloc <- apply(abundances, 1, sum)
      nbsploc <- apply(abundances, 1, function(x) {
        length(which(x > 0))
      })
      
      # relative abundances inside each local community
      locabrel <- abundances / abloc
      
      # alpha diversity
      Qalpha = apply(locabrel, 1, function(x)
        t(x) %*%  functdist %*% x)
      
      #Wc
      Wc = abloc / sum(abloc)
      
      # abundance-weighted mean alpha
      mQalpha <- as.numeric(Qalpha %*% abloc / sum(abloc))
      
      #Villeger's correction
      if (w == T) {
        # abundance-weighted mean alpha
        mQalpha <- as.numeric(Qalpha %*% abloc / sum(abloc))
        totabrel <- apply(abundances, 2, sum) / sum(abundances)
        Qalpha = Qalpha * Wc
      }
      
      # Rao's original definition: mean of Pi
      else {
        mQalpha <- mean(Qalpha)
        totabrel <- apply(locabrel, 2, mean)
      }
      
      
      # gamma diversity
      Qgamma <- (totabrel %*% functdist %*% totabrel) [1]
      
      # beta diversity
      Qbeta <- as.numeric(Qgamma - mQalpha)
      
      # standardized beta diversity
      Qbetastd <- as.numeric(Qbeta / Qgamma)
      
      # list of results
      resQ <-
        list(
          Richness_per_plot = nbsploc,
          Relative_abundance = locabrel,
          Pi = totabrel,
          Wc = Wc,
          Species_abundance_per_plot = abloc,
          Alpha = Qalpha,
          Mean_alpha = mQalpha,
          Gamma = Qgamma,
          Beta = Qbeta,
          Standardize_Beta = Qbetastd
        )
      
      return(resQ)
      
    }
    
    
    ###########function disc originally from S. Pavoine####
    
    disc = function (samples,
                     dis = NULL,
                     structures = NULL,
                     Jost = F)
    {
      if (!inherits(samples, "data.frame"))
        stop("Non convenient samples")
      if (any(samples < 0))
        stop("Negative value in samples")
      if (any(apply(samples, 2, sum) < 1e-16))
        stop("Empty samples")
      if (!is.null(dis)) {
        if (!inherits(dis, "dist"))
          stop("Object of class 'dist' expected for distance")
        # if (!is.euclid(dis))
        #stop("Euclidean property is expected for distance")
        dis <- as.matrix(dis)
        if (nrow(samples) != nrow(dis))
          stop("Non convenient samples")
      }
      if (!is.null(structures)) {
        if (!inherits(structures, "data.frame"))
          stop("Non convenient structures")
        m <- match(apply(structures, 2, function(x)
          length(x)),
          ncol(samples), 0)
        if (length(m[m == 1]) != ncol(structures))
          stop("Non convenient structures")
        m <-
          match(tapply(1:ncol(structures), as.factor(1:ncol(structures)),
                       function(x)
                         is.factor(structures[, x])),
                TRUE,
                0)
        if (length(m[m == 1]) != ncol(structures))
          stop("Non convenient structures")
      }
      Structutil <- function(dp2, Np, unit, Jost) {
        if (!is.null(unit)) {
          modunit <- model.matrix( ~ -1 + unit)
          sumcol <- apply(Np, 2, sum)
          Ng <- modunit * sumcol
          lesnoms <- levels(unit)
        }
        else {
          Ng <- as.matrix(Np)
          lesnoms <- colnames(Np)
        }
        sumcol <- apply(Ng, 2, sum)
        Lg <- t(t(Ng) / sumcol)
        colnames(Lg) <- lesnoms
        Pg <- as.matrix(apply(Ng, 2, sum) / nbhaplotypes)
        rownames(Pg) <- lesnoms
        deltag <- as.matrix(apply(Lg, 2, function(x)
          t(x) %*%
            dp2 %*% x))
        ug <- matrix(1, ncol(Lg), 1)
        if (Jost) {
          #dp2 <- as.matrix(as.dist(dfunct01))
          deltag <-
            as.matrix(apply(Lg, 2, function(x)
              t(x) %*% dp2 %*% x))
          X = t(Lg) %*% dp2 %*% Lg
          alpha = 1 / 2 * (deltag %*% t(ug) + ug %*% t(deltag))
          Gam = (X + alpha) / 2
          alpha = 1 / (1 - alpha) #Jost correction
          Gam = 1 / (1 - Gam)  #Jost correction
          Beta_add = Gam - alpha
          Beta_mult = 100 * (Gam - alpha) / Gam
        }
        else {
          deltag <- as.matrix(apply(Lg, 2, function(x)
            t(x) %*% dp2 %*% x))
          X = t(Lg) %*% dp2 %*% Lg
          alpha = 1 / 2 * (deltag %*% t(ug) + ug %*% t(deltag))
          Gam = (X + alpha) / 2
          Beta_add = Gam - alpha
          Beta_mult = 100 * (Gam - alpha) / Gam
        }
        colnames(Beta_add) <- lesnoms
        rownames(Beta_add) <- lesnoms
        return(
          list(
            Beta_add = as.dist(Beta_add),
            Beta_mult = as.dist(Beta_mult),
            Gamma = as.dist(Gam),
            Alpha = as.dist(alpha),
            Ng = Ng,
            Pg = Pg
          )
        )
      }
      Diss <- function(dis,
                       nbhaplotypes,
                       samples,
                       structures,
                       Jost) {
        structutil <- list(0)
        structutil[[1]] <-
          Structutil(dp2 = dis, Np = samples, NULL, Jost)
        diss <-
          list(
            structutil[[1]]$Alpha,
            structutil[[1]]$Gamma,
            structutil[[1]]$Beta_add,
            structutil[[1]]$Beta_mult
          )
        if (!is.null(structures)) {
          for (i in 1:length(structures)) {
            structutil[[i + 1]] <-
              Structutil(as.matrix(structutil[[1]]$Beta_add),
                         structutil[[1]]$Ng,
                         structures[, i],
                         Jost)
          }
          diss <-
            c(diss, tapply(1:length(structures), factor(1:length(structures)),
                           function(x)
                             as.dist(structutil[[x + 1]]$Beta_add)))
        }
        return(diss)
      }
      nbhaplotypes <- sum(samples)
      diss <- Diss(dis, nbhaplotypes, samples, structures, Jost)
      if (!is.null(structures)) {
        names(diss) <-
          c("Alpha", "Gamma", "Beta_add", "Beta_prop", "Beta_region")
        return(diss)
      }
      names(diss) <- c("Alpha", "Gamma", "Beta_add", "Beta_prop")
      return(diss)
    }
    
    
    
    
    
    
    TD <- FD <- PD <- NULL
    
    #Taxonomic diversity
    dS <-
      matrix(1, nrow(sample), nrow(sample)) - diag(rep(1, nrow(sample)))
    temp_qdec <-
      Qdecomp(dS, t(sample), w = weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
    TD$Richness_per_plot = temp_qdec$Richness_per_plot
    TD$Relative_abundance = temp_qdec$Relative_abundance
    TD$Pi = temp_qdec$Pi
    TD$Wc = temp_qdec$Wc
    if (Jost) {
      TD$Mean_Alpha = 1 / (1 - temp_qdec$Mean_alpha)
      TD$Alpha = 1 / (1 - temp_qdec$Alpha)
      TD$Gamma = 1 / (1 - temp_qdec$Gamma)
      TD$Beta_add = (TD$Gamma - TD$Mean_Alpha)
      TD$Beta_prop = 100 * TD$Beta_add / TD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      TD$Pairwise_samples <-
        disc(as.data.frame(sample),
             as.dist(dS),
             structure = structure,
             Jost = Jost)
    }
    else {
      TD$Mean_Alpha = temp_qdec$Mean_alpha
      TD$Alpha = temp_qdec$Alpha
      TD$Gamma = temp_qdec$Gamma
      TD$Beta_add = (TD$Gamma - TD$Mean_Alpha)
      TD$Beta_prop = 100 * TD$Beta_add / TD$Gamma
      #Call the disc function for alpha, gamma and beta estimations for each pair of samples
      TD$Pairwise_samples <-
        disc(as.data.frame(sample),
             as.dist(dS),
             structure = structure,
             Jost = Jost)
    }
    
    #Functional diversity estimation
    if (!is.null(dfunc)) {
      FD <- list()
      if (Jost) {
        if (max(dfunc) > 1)
          dfunc <-
            dfunc / max(dfunc)   #Make sure the distance are between 0 and 1 for the Jost correction
        temp_qdec <-
          Qdecomp(dfunc, t(sample), w = weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
        #  FD$Alpha = 1/(1-temp_qdec$Alpha)
        #  FD$Mean_Alpha = mean(FD$Alpha)
        FD$Mean_Alpha = 1 / (1 - temp_qdec$Mean_alpha)
        FD$Alpha = 1 / (1 - temp_qdec$Alpha)
        FD$Gamma = 1 / (1 - temp_qdec$Gamma)
        FD$Beta_add = (FD$Gamma - FD$Mean_Alpha)
        FD$Beta_prop = 100 * FD$Beta_add / FD$Gamma
        #Call the disc function for alpha, gamma and beta estimations for each pair of samples
        FD$Pairwise_samples <-
          disc(
            as.data.frame(sample),
            as.dist(dfunc),
            structure = structure,
            Jost = Jost
          )
      }
      else {
        temp_qdec <-
          Qdecomp(dfunc, t(sample), w = weight) #Call the Qdecomp function for alpha, gamma and beta estimations.
        FD$Mean_Alpha = temp_qdec$Mean_alpha
        FD$Alpha = temp_qdec$Alpha
        FD$Gamma = temp_qdec$Gamma
        FD$Beta_add = (FD$Gamma - FD$Mean_Alpha)
        FD$Beta_prop = 100 * FD$Beta_add / FD$Gamma
        #FD$Beta =  temp_qdec$Beta#
        #Call the disc function for alpha, gamma and beta estimations for each pair of samples
        FD$Pairwise_samples <-
          disc(
            as.data.frame(sample),
            as.dist(dfunc),
            structure = structure,
            Jost = Jost
          )
      }
    }
    #Phylogenetic diversity estimation
    if (!is.null(dphyl)) {
      PD <- list()
      if (Jost) {
        if (max(dphyl) > 1)
          dphyl <-
            dphyl / max(dphyl)   #Make sure the distance are between 0 and 1 for the Jost correction
        temp_qdec <-
          Qdecomp(dphyl, t(sample), w = weight)   #Call the Qdecomp function for alpha, gamma and beta estimations.
        PD$Mean_Alpha = 1 / (1 - temp_qdec$Mean_alpha)
        PD$Alpha = 1 / (1 - temp_qdec$Alpha)
        PD$Gamma = 1 / (1 - temp_qdec$Gamma)
        PD$Beta_add = (PD$Gamma - PD$Mean_Alpha)
        PD$Beta_prop = 100 * PD$Beta_add / PD$Gamma
        #Call the disc function for alpha, gamma and beta estimations for each pair of samples
        PD$Pairwise_samples <-
          disc(
            as.data.frame(sample),
            as.dist(dphyl),
            structure = structure,
            Jost = Jost
          )
      }
      else {
        temp_qdec <-
          Qdecomp(dphyl, t(sample), w = weight)  #Call the Qdecomp function for alpha, gamma and beta estimations.
        PD$Mean_Alpha = temp_qdec$Mean_alpha
        PD$Alpha = temp_qdec$Alpha
        PD$Gamma = temp_qdec$Gamma
        PD$Beta_add = (PD$Gamma - PD$Mean_Alpha)
        PD$Beta_prop = 100 * PD$Beta_add / PD$Gamma
        #PD$Beta =  temp_qdec$Beta
        #Call the disc function for alpha, gamma and beta estimations for each pair of samples
        PD$Pairwise_samples <-
          disc(
            as.data.frame(sample),
            as.dist(dphyl),
            structure = structure,
            Jost = Jost
          )
      }
      
      
      
      
    }
    out <- list(TD, FD, PD)
    names(out) <- c("TD", "FD", "PD")
    return(out)
    
  }
