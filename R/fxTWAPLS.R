#' Funtion to get the frequency of the climate value, which will be used to 
#'     provide fx correction for WA-PLS and TWA-PLS
#'
#' @param x the modern climate values 
#' @param bin binwidth to get the frequency of the modern climate values
#'
#' @return the frequency of the modern climate values
#' @export
#'
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                          "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                          package = "fxTWAPLS", 
#'                                          mustWork = TRUE))
#'     
#'     # Extract taxa
#'     taxaColMin <- which(colnames(modern_pollen) == "Abies")
#'     taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
#'     taxa <- modern_pollen[, taxaColMin:taxaColMax]
#'     
#'     # Get the frequency of each climate variable fx
#'     fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#'     fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#'     fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' }
fx <- function(x, bin) {
  pbin <- round((max(x) - min(x)) / bin, digits = 0)
  bin <- (max(x) - min(x)) / pbin
  hist <- hist(x, breaks = seq(min(x), max(x), by = bin))
  xbin <- seq(min(x) + bin / 2, max(x) - bin / 2, by = bin)
  counts <- hist[["counts"]]
  fx <- rep(NA, length(x))
  for (i in 1:length(x)) {
    fx[i] <- counts[which.min(abs(x[i] - xbin))]
  }
  if (any(fx == 0)) {
    print("Some x have a count of 0!")
  }
  plot(fx ~ x)
  return(fx)
}

#' WALPS training function, which can choose to perform fx correction
#' 
#' @importFrom stats lm
#' 
#' @param modern_taxa the modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate the modern climate value at each sampling site
#' @param nPLS the number of components to be extracted
#' @param usefx boolean flag on whether or not use fx correction.
#' @param fx the frequency of the climate value for fx correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#'
#' @return a list of the training results, which will be used by the predict 
#'     function. fit is the fitted value of modern training result.
#' @export
#'
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                          "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                          package = "fxTWAPLS", 
#'                                          mustWork = TRUE))
#'     
#'     # Extract taxa
#'     taxaColMin <- which(colnames(modern_pollen) == "Abies")
#'     taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
#'     taxa <- modern_pollen[, taxaColMin:taxaColMax]
#'     
#'     # Get the frequency of each climate variable fx
#'     fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#'     fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#'     fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#'     # MTCO
#'     fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#'     fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
#'                                     modern_pollen$Tmin, 
#'                                     nPLS = 5, 
#'                                     usefx = TRUE, 
#'                                     fx = fx_Tmin)
#' }
#' 
#' @seealso \code{\link{fx}} and \code{\link{TWAPLS.w}}
WAPLS.w <- function(modern_taxa, 
                    modern_climate, 
                    nPLS = 5, 
                    usefx = FALSE, 
                    fx = NA) {
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  # Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  # u of each component, standardized the same way as r
  u_sd <- matrix(NA, nc, nPLS) 
  optimum <- matrix(NA, nc, nPLS) # u updated
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  
  pls <- 1
  # Step 1. Take the centred environmental variable(xi) as initial site scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate new species scores (uk* by weighted averaging of the site scores)
  u[, pls] <- t(y) %*% x / sumi_yik # uk = sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  
  # Step 3. Calculate new site scores (ri) by weighted averaging of the species scores
  r[, pls] <- y %*% u[, pls] / sumk_yik # xi=sumk_yik*uk/sumk_yik; 1*nsite
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable (xJ on the components obtained so 
  # far using weights and take the fitted values as current estimates
  if(usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = sumk_yik / Ytottot)
  } else{
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site scores (rJ.
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new species scores (uk* by weighted averaging of the site scores)
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    
    # Step 3. Calculate new site scores (ri) by weighted averaging of the species scores
    r[, pls] <- y %*% u[, pls] / sumk_yik # xi=sumk_yik*uk/sumk_yik; 1*nsite
    
    # Step 4. For second and higher components, make the new site scores (ri) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized score as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable on the components obtained so 
    # far using weights and take the fitted values as current estimates 
    if(usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = sumk_yik / Ytottot)
    } else{
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, modern_climate, colnames(modern_taxa), optimum,comp, u, z, s, orth, alpha, mean(modern_climate), nPLS)
  names(list) <- c(c("fit", "x", "taxon_name", "optimum", "comp", "u", "z", "s", "orth", "alpha", "meanx", "nPLS"))
  return(list)
}

#' TWALPS training function, which can choose to perform fx correction
#' 
#' @importFrom stats lm
#' 
#' @param modern_taxa the modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate the modern climate value at each sampling site
#' @param nPLS the number of components to be extracted
#' @param usefx boolean flag on whether or not use fx correction.
#' @param fx the frequency of the climate value for fx correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#'
#' @return a list of the training results, which will be used by the predict 
#'     function. fit is the fitted value of modern training result.
#' @export
#'
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                          "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                          package = "fxTWAPLS", 
#'                                          mustWork = TRUE))
#'     
#'     # Extract taxa
#'     taxaColMin <- which(colnames(modern_pollen) == "Abies")
#'     taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
#'     taxa <- modern_pollen[, taxaColMin:taxaColMax]
#'     
#'     # Get the frequency of each climate variable fx
#'     fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#'     fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#'     fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#' 
#'     # MTCO
#'     fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#'     fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
#'                                       modern_pollen$Tmin, 
#'                                       nPLS = 5, 
#'                                       usefx = TRUE, 
#'                                       fx = fx_Tmin)
#' }
#' 
#' @seealso \code{\link{fx}} and \code{\link{WAPLS.w}}
TWAPLS.w <- function(modern_taxa,
                     modern_climate,
                     nPLS = 5,
                     usefx = FALSE,
                     fx = NA){
  # Step 0. Centre the environmental variable by subtracting the weighted mean
  x <- modern_climate
  y <- modern_taxa
  y <- as.matrix(y)
  nc <- ncol(modern_taxa)
  nr <- nrow(modern_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  #Define some matrix to store the values
  u <- matrix(NA, nc, nPLS) # u of each component
  u_sd <- matrix(NA, nc, nPLS) # u of each component, standardized the same way as r
  optimum <- matrix(NA, nc, nPLS) # u updated
  t <- matrix(NA, nc, nPLS) # tolerance
  r <- matrix(NA, nr, nPLS) # site score
  z <- matrix(NA, 1, nPLS) # standardize
  s <- matrix(NA, 1, nPLS) # standardize
  orth <- list() # store orthogonalization parameters
  alpha <- list() # store regression coefficients
  comp <- matrix(NA, nr, nPLS) # each component
  fit <- matrix(NA, nr, nPLS) # current estimate
  
  pls <- 1
  # Step 1. Take the centred environmental variable (xi) as initial site scores (ri). 
  r[, pls] <- x - mean(x)
  
  # Step 2. Calculate uk and tk
  u[, pls] <- t(y) %*% x / sumi_yik # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  n2 <- matrix(NA, nc, 1)
  for (k in 1:nc) {
    t[k, pls] <- sqrt(sum(y[, k] * (x - u[k, pls]) ^ 2) / sumi_yik[k])
    n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
    t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
  }
  
  # Step 3. Calculate new site scores (ri)
  r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2)) #xi; 1*nsite
  
  # Step 4. For the first axis go to Step 5.
  
  # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[, pls] <- mean(r[, pls], na.rm = TRUE)
  s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  
  # Step 6. Take the standardized score as the new component
  comp[, pls] <- r[, pls]
  
  # Step 7. Regress the environmental variable on the components obtained so far
  # using weights and take the fitted values as current estimates 
  if (usefx == FALSE) {
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = sumk_yik / Ytottot)
  } else{
    lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
  }
  
  fit[, pls] <- lm[["fitted.values"]]
  alpha[[pls]] <- lm[["coefficients"]]
  u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
  optimum[, pls] <- alpha[[pls]][1] + u_sd[, pls] * alpha[[pls]][2]
  
  for (pls in 2:nPLS) {
    # Go to Step 2 with the residuals of the regression as the new site scores (ri).
    r[, pls] <- lm[["residuals"]]
    
    # Step 2. Calculate new uk and tk
    u[, pls] <- t(y) %*% r[, pls] / sumi_yik # uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    n2 <- matrix(NA, nc, 1)
    for (k in 1:nc) {
      t[k, pls] <- sqrt(sum(y[, k] * (r[, pls] - u[k, pls]) ^ 2) / sumi_yik[k])
      n2[k] <- 1 / sum((y[, k] / sum(y[, k])) ^ 2)
      t[k, pls] <- t[k, pls] / sqrt(1 - 1 / n2[k])
    }
    
    # Step 3. Calculate new site scores (r;) by weighted averaging of the species scores, i.e. new
    r[,pls] <- (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
    
    # Step 4. For second and higher components, make the new site scores (r;) 
    # uncorrelated with the previous components by orthogonalization 
    # (Ter Braak, 1987 : Table 5 .2b)
    v <- rep(NA, pls - 1)
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      v[pls - j] <- sum(sumk_yik * fi * xi) / Ytottot
      xinew <- xi - v[pls - j] * fi
    }
    orth[[pls]] <- v
    # plot(xinew~r[,pls]);abline(0,1)
    r[, pls] <- xinew
    
    # Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[, pls] <- mean(r[, pls], na.rm = TRUE)
    s[, pls] <- sqrt(sum((r[, pls] - z[, pls]) ^ 2, na.rm = TRUE) / Ytottot)
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    
    # Step 6. Take the standardized scores as the new component
    comp[, pls] <- r[, pls]
    
    # Step 7. Regress the environmental variable (xJ on the components obtained 
    # so far using weights and take the fitted values as current estimates 
    if (usefx == FALSE) {
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = sumk_yik / Ytottot)
    } else{
      lm <- MASS::rlm(modern_climate ~ comp[, 1:pls], weights = 1 / fx ^ 2)
    }
    
    fit[, pls] <- lm[["fitted.values"]]
    alpha[[pls]] <- lm[["coefficients"]]
    
    u_sd[, pls] <- (u[, pls] - z[, pls]) / s[, pls]
    optimum[, pls] <-
      alpha[[pls]][1] + u_sd[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, modern_climate, colnames(modern_taxa), optimum, comp, u, t, z, s, orth, alpha, mean(modern_climate), nPLS)
  names(list) <- c(c("fit", "x", "taxon_name", "optimum", "comp", "u", "t", "z", "s", "orth", "alpha", "meanx", "nPLS"))
  return(list)
}

#' WAPLS predict function
#'
#' @param WAPLSoutput the output of the \code{\link{WAPLS.w}} training function, 
#'     either with or without fx correction
#' @param fossil_taxa fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents a 
#'     taxon.
#'
#' @return a list of the reconstruction results. fit is the fitted value.
#' @export
#'
# @examples
WAPLS.predict.w <- function(WAPLSoutput, fossil_taxa) {
  y <- fossil_taxa
  y <- as.matrix(y)
  nc <- ncol(fossil_taxa)
  nr <- nrow(fossil_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  nPLS <- WAPLSoutput[["nPLS"]]
  meanx <- WAPLSoutput[["meanx"]]
  u <- WAPLSoutput[["u"]]
  z <- WAPLSoutput[["z"]]
  s <- WAPLSoutput[["s"]]
  orth <- WAPLSoutput[["orth"]]
  alpha <- WAPLSoutput[["alpha"]]
  
  if (nc != nrow(u)) {
    print("Number of taxa doesn't match!")
  }
  if (all(colnames(fossil_taxa) == WAPLSoutput[["taxon_name"]]) == FALSE) {
    print("Taxa don't match!")
  }
  
  # Define some matrix to store the values
  fit <- matrix(NA, nr, nPLS)
  r <- matrix(NA, nr, nPLS)
  comp <- matrix(NA, nr, nPLS)
  
  pls <- 1
  # xi=sumk_yik*uk/sumk_yik; 1*nsite
  r[, pls] <- y %*% u[, pls] / sumk_yik # xi=sumk_yik*uk/sumk_yik; 1*nsite
  # standardize the same way
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  # multiply the same regression coefficients
  comp[, pls] <- r[, pls]
  fit[, 1] <- alpha[[pls]][1] + comp[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # xi = sumk_yik*uk/sumk_yik; 1*nsite
    r[, pls] <- y %*% u[, pls] / sumk_yik # xi = sumk_yik*uk/sumk_yik; 1*nsite
    # orthoganlization the same way
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      xinew <- xi - orth[[pls]][pls - j] * fi
    }
    r[, pls] <- xinew
    
    # standardize the same way
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    # multiply the same regression coefficients
    comp[, pls] <- r[, pls]
    fit[, pls] <-
      alpha[[pls]][1] + comp[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, nPLS)
  names(list) <- c(c("fit", "nPLS"))
  return(list)
}

#' TWAPLS predict function
#'
#' @param TWAPLSoutput the output of the \code{\link{TWAPLS.w}} training 
#'     function, either with or without fx correction
#' @param fossil_taxa fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents 
#'     a taxon.
#'
#' @return a list of the reconstruction results. fit is the fitted value.
#' @export
#'
# @examples
TWAPLS.predict.w <- function(TWAPLSoutput, fossil_taxa) {
  y <- fossil_taxa
  y <- as.matrix(y)
  nc <- ncol(fossil_taxa)
  nr <- nrow(fossil_taxa)
  Ytottot <- sum(y)
  sumk_yik <- rowSums(y)
  sumi_yik <- colSums(y)
  
  nPLS <- TWAPLSoutput[["nPLS"]]
  meanx <- TWAPLSoutput[["meanx"]]
  u <- TWAPLSoutput[["u"]]
  t <- TWAPLSoutput[["t"]]
  z <- TWAPLSoutput[["z"]]
  s <- TWAPLSoutput[["s"]]
  orth <- TWAPLSoutput[["orth"]]
  alpha <- TWAPLSoutput[["alpha"]]
  
  if (nc != nrow(u)) {
    print("Number of taxa doesn't match!")
  }
  if (all(colnames(fossil_taxa) == TWAPLSoutput[["taxon_name"]]) == FALSE) {
    print("Taxa don't match!")
  }
  
  # Define some matrix to store the values
  fit <- matrix(NA, nr, nPLS)
  r <- matrix(NA, nr, nPLS)
  comp <- matrix(NA, nr, nPLS)
  
  pls <- 1
  # xi=sumk_yik*uk/sumk_yik; 1*nsite
  r[, pls] <- (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2)) #xi; 1*nsite
  # standardize the same way
  r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
  # multiply the same regression coefficients
  comp[, pls] <- r[, pls]
  fit[, 1] <- alpha[[pls]][1] + comp[, pls] * alpha[[pls]][2]
  
  for(pls in 2:nPLS) {
    # xi=sumk_yik*uk/sumk_yik; 1*nsite
    r[, pls] = (y %*% (u[, pls] / t[, pls] ^ 2)) / (y %*% (1 / t[, pls] ^ 2)) # xi; 1*nsite
    # orthoganlization the same way
    for (j in 1:(pls - 1)) {
      fi <- r[, pls - j]
      xi <- r[, pls]
      xinew <- xi - orth[[pls]][pls - j] * fi
    }
    r[, pls] <- xinew
    
    # standardize the same way
    r[, pls] <- (r[, pls] - z[, pls]) / s[, pls]
    # multiply the same regression coefficients
    comp[, pls] <- r[, pls]
    fit[, pls] <-
      alpha[[pls]][1] + comp[, 1:pls] %*% as.matrix(alpha[[pls]][2:(pls + 1)])
  }
  
  list <- list(fit, nPLS)
  names(list) <- c(c("fit", "nPLS"))
  return(list)
}

#' Calculate Sample Specific Errors
#'
#' @param modern_taxa the modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate the modern climate value at each sampling site
#' @param fossil_taxa fossil taxa abundance data to reconstruct past climates, 
#'     each row represents a site to be reconstructed, each column represents a 
#'     taxon.
#' @param trainfun training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}
#' @param predictfun predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}
#' @param nboot the number of bootstrap cycles you want to use
#' @param nPLS the number of components to be extracted
#' @param nsig the significant number of components to use to reconstruct past 
#'     climates, this can be obtained from the cross-validation results.
#' @param usefx boolean flag on whether or not use fx correction.
#' @param fx the frequency of the climate value for fx correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param seed seed for reproducibility
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purpouses only.
#' @param test_it number of iterations to use in the test mode
#'
#' @return the bootstrapped standard error for each site
#' @export
#'
# @examples
sse.sample <- function(modern_taxa,
                       modern_climate,
                       fossil_taxa,
                       trainfun,
                       predictfun,
                       nboot,
                       nPLS,
                       nsig,
                       usefx,
                       fx,
                       cpus = 4,
                       seed = NULL,
                       test_mode = FALSE,
                       test_it = 5) {
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  # `%dorng%` <- doRNG::`%dorng%`
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Make list of row numbers by sampling with 
  # replacement
  k_samples <- replicate(nboot, sample(1:nrow(modern_taxa),
                                         size = nrow(modern_taxa),
                                         replace = TRUE))
  
  # Create list of indices to loop through
  idx <- 1:nboot
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  xboot <- foreach::foreach(i = idx,
                            .combine = cbind) %dopar% {
                              # Extract list of row numbers by sampling with 
                              # replacement
                              k <- k_samples[, i]
                              # k <- sample(1:nrow(modern_taxa),
                              #             size = nrow(modern_taxa),
                              #             replace = TRUE)
                              
                              # Reorganise modern_taxa obs in k order
                              modern_taxa <- modern_taxa[k, ] 
                              modern_climate <- modern_climate[k]
                              col_not0 <- which(colSums(modern_taxa) > 0)
                              # Strip out zero-sum cols
                              modern_taxa <- modern_taxa[, col_not0]
                              # Apply train function, with modern_climate also 
                              # in k order
                              mod <- trainfun(modern_taxa, modern_climate)
                              # Make reconstruction
                              predictfun(mod, 
                                         fossil_taxa[, col_not0])$fit[, nsig]
                            }
  parallel::stopCluster(cl) # Stop cluster
  
  avg.xboot <- rowMeans(xboot, na.rm = TRUE)
  v1 <- boot.mean.square <- rowMeans((xboot - avg.xboot) ^ 2 , na.rm = TRUE)
  return(sqrt(v1))
}

#' Leave one out cross validation as 
#' rioja (\url{https://cran.r-project.org/package=rioja})
#' 
#' @importFrom foreach `%dopar%`
#' 
#' @param modern_taxa the modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate the modern climate value at each sampling site
#' @param nPLS the number of components to be extracted
#' @param trainfun training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}
#' @param predictfun predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}
#' @param usefx boolean flag on whether or not use fx correction.
#' @param fx the frequency of the climate value for fx correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purpouses only.
#' @param test_it number of iterations to use in the test mode
#'
#' @return leave-one-out cross validation results
#' @export
#'
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                          "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                          package = "fxTWAPLS", 
#'                                          mustWork = TRUE))
#'     
#'     # Extract taxa
#'     taxaColMin <- which(colnames(modern_pollen) == "Abies")
#'     taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
#'     taxa <- modern_pollen[, taxaColMin:taxaColMax]
#'     
#'     # Get the frequency of each climate variable fx
#'     fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
#'     fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
#'     fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)
#'     
#'     # MTCO
#'     ## fx
#'     fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
#'     
#'     ## LOOCV
#'     ### without fx
#'     cv_Tmin <- fxTWAPLS::cv.w(taxa,
#'                               modern_pollen$Tmin,
#'                               nPLS = 5,
#'                               fxTWAPLS::WAPLS.w,
#'                               fxTWAPLS::WAPLS.predict.w,
#'                               cpus = 1)
#'     ### with fx
#'     cv_f_Tmin <- fxTWAPLS::cv.w(taxa,
#'                                 modern_pollen$Tmin,
#'                                 nPLS = 5,
#'                                 fxTWAPLS::WAPLS.w,
#'                                 fxTWAPLS::WAPLS.predict.w,
#'                                 usefx = TRUE,
#'                                 fx = fx_Tmin,
#'                                 cpus = CPUS)  
#' }
cv.w <- function(modern_taxa,
                 modern_climate,
                 nPLS = 5,
                 trainfun,
                 predictfun,
                 usefx = FALSE,
                 fx = NA,
                 cpus = 4,
                 test_mode = FALSE,
                 test_it = 5) {
  x <- modern_climate
  y <- modern_taxa
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)

  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- 1:length(x)
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  all.cv.out <- foreach::foreach(i = idx,
                                 .combine = rbind, #comb_pb(max(idx)),
                                 .verbose = FALSE) %dopar% {
                                   fit <- trainfun(y[-i, ], x[-i], nPLS, usefx, fx[-i])
                                   xnew <- predictfun(fit, y[i, ])[["fit"]]
                                   data.frame(x[i], xnew)
                                 }
  parallel::stopCluster(cl) # Stop cluster
  colnames(all.cv.out) <- c("test.x", paste0("comp", 1:nPLS))
  return(all.cv.out)
}

#' Leave out cross validation with pseudo (geographically and climatically close) 
#' sites removed from the training set.
#' Get the ones which are both geographically and climatically close, and could 
#' therefore results in pseudo-replication.
#'
#' @param dist distance matrix which contains the distance from other sites.
#' @param x the modern climate values
#' @param cpus number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purpouses only.
#' @param test_it number of iterations to use in the test mode
#' 
#' @return the geographically and climatically close sites to each test site.
#' @export
#'
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                           "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                           package = "fxTWAPLS", 
#'                                           mustWork = TRUE))
#'     point <- modern_pollen[, c("Long", "Lat")]
#'     dist <- get_distance(point, cpus = 1)
#'     pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, 
#'                                         modern_pollen$Tmin, 
#'                                         cpus = 1)
#' }
#' @seealso \code{\link{get_distance}}
get_pseudo <- function(dist, x, cpus = 4, test_mode = FALSE, test_it = 5) {
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- 1:length(x)
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  pseudo <- foreach::foreach(i = idx) %dopar% {
    which(dist[i, ] < 50000 & abs(x - x[i]) < 0.02 * (max(x) - min(x)))
  }
  parallel::stopCluster(cl) # Stop cluster
  return(pseudo)
}

#' Pseudo removed leave out cross validation
#'
#' @param modern_taxa the modern taxa abundance data, each row represents a 
#'     sampling site, each column represents a taxon.
#' @param modern_climate the modern climate value at each sampling site
#' @param nPLS the number of components to be extracted
#' @param trainfun training function you want to use, either 
#'     \code{\link{WAPLS.w}} or \code{\link{TWAPLS.w}}
#' @param predictfun predict function you want to use: if \code{trainfun} is 
#'     \code{\link{WAPLS.w}}, then this should be \code{\link{WAPLS.predict.w}}; 
#'     if \code{trainfun} is \code{\link{TWAPLS.w}}, then this should be 
#'     \code{\link{TWAPLS.predict.w}}
#' @param pseudo the geographically and climatically close sites to each test 
#'     site, obtained from \code{\link{get_pseudo}} function
#' @param usefx boolean flag on whether or not use fx correction.
#' @param fx the frequency of the climate value for fx correction: if 
#'     \code{usefx} is FALSE, this should be \code{NA}; otherwise, this should 
#'     be obtained from the \code{\link{fx}} function.
#' @param cpus number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purpouses only.
#' @param test_it number of iterations to use in the test mode
#'
#' @return leave-one-out cross validation results
#' @export
#'
# @examples
cv.pr.w <- function(modern_taxa,
                    modern_climate,
                    nPLS = 5,
                    trainfun,
                    predictfun,
                    pseudo,
                    usefx = FALSE,
                    fx = NA,
                    cpus = 4,
                    test_mode = TRUE,
                    test_it = 5) {
  x <- modern_climate
  y <- modern_taxa
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- 1:length(x)
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  all.cv.out <- foreach::foreach(i = idx,
                                 .combine = rbind) %dopar% {
                                   leave <- unlist(pseudo[i])
                                   fit <- trainfun(y[-leave, ], x[-leave], nPLS, usefx, fx[-leave])
                                   xnew <- predictfun(fit, y[i, ])[["fit"]]
                                   data.frame(x[i], xnew)
                                 }
  parallel::stopCluster(cl) # Stop cluster
  
  # assign column names to all.cv.out
  colnames(all.cv.out) <- c("test.x", paste0("comp", 1:nPLS))
  
  return(all.cv.out)
}


#' Random t-test 
#' 
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats rbinom
#' 
#' @param cvoutput cross-validation output either from \code{\link{cv.w}} or 
#'     \code{\link{cv.pr.w}}
#' @param n.perm the number of permutation times to get the p value, which 
#'     assesses whether using the current number of components is significantly 
#'     different from using one less.
#'
#' @return a matrix of the statistics of the cross-validation results
#' @export
#'
# @examples
rand.t.test.w <- function(cvoutput, n.perm = 999) {
  ncomp <- ncol(cvoutput) - 1
  output <- matrix(NA, ncomp, 11)
  colnames(output) <- c("R2",
                        "Avg.Bias",
                        "Max.Bias",
                        "Min.Bias",
                        "RMSEP",
                        "delta.RMSEP",
                        "p",
                        "Compre.b0",
                        "Compre.b1",
                        "Compre.b0.se",
                        "Compre.b1.se")
  
  for (i in 1:ncomp) {
    cv.x <- cvoutput[, 1]
    cv.i <- cvoutput[, 1 + i]
    output[i, "RMSEP"] <- sqrt(mean((cv.i - cv.x) ^ 2))
    output[i, "R2"] <- cor(cv.i, cv.x) ^ 2
    output[i, "Avg.Bias"] <- mean(cv.i - cv.x)
    output[i, "Max.Bias"] <- max(abs(cv.i - cv.x))
    output[i, "Min.Bias"] <- min(abs(cv.i - cv.x))
    output[i, c("Compre.b0", "Compre.b1")] <-
      summary(lm(cv.i ~ cv.x))[["coefficients"]][, "Estimate"]
    output[i, c("Compre.b0.se", "Compre.b1.se")] <-
      summary(lm(cv.i ~ cv.x))[["coefficients"]][, "Std. Error"]
  }
  # get delta.RMSEP
  for(i in 1:ncomp) {
    if (i == 1) {
      rmsep.null <- sqrt(mean((cv.x - mean(cv.x)) ^ 2))
      output[i, "delta.RMSEP"] <-
        (output[i, "RMSEP"] - rmsep.null) * 100 / rmsep.null
    } else{
      output[i, "delta.RMSEP"] <-
        (output[i, "RMSEP"] - output[i - 1, "RMSEP"]) * 100 / output[i - 1, "RMSEP"]
    }
  }
  # get p-value, which describes whether using the number of components now has 
  # a significant difference than using one less
  e0 <- cv.x - mean(cv.x)
  e <- cbind(e0, cvoutput[, 2:ncol(cvoutput)] - cv.x)
  
  t.res <- vector("numeric", ncomp)
  t <- vector("numeric", n.perm + 1)
  t.res[] <- NA
  n <- nrow(e)
  for (i in 1:ncomp) {
    d <- e[, i] ^ 2 - e[, i + 1] ^ 2
    t[1] <- mean(d, na.rm = TRUE)
    for (j in 1:n.perm) {
      sig <- 2 * rbinom(n, 1, 0.5) - 1
      t[j + 1] <- mean(d * sig, na.rm = TRUE)
    }
    t.res[i] <- sum(t >= t[1]) / (n.perm + 1)
  }
  output[, "p"] <- t.res
  
  print(output)
  return(output)
}

#' Get the distance between points
#' 
#' @importFrom foreach `%dopar%` 
#' 
#' @param point each row represents a sampling site, the first column is 
#'     longitude and the second column is latitude, both in decimal format
#' @param cpus number of CPUs for simultaneous iterations to execute, check
#'     \code{parallel::detectCores()} for available CPUs on your machine.
#' @param test_mode boolean flag to execute the function with a limited number
#'     of iterations, \code{test_it}, for testing purpouses only.
#' @param test_it number of iterations to use in the test mode
#'    
#' @return distance matrix, the value at the ith row, means the distance between 
#'     the ith sampling site and the whole sampling sites
#' @export
#' 
#' @examples
#' \dontrun{
#'     # Load modern pollen data
#'     modern_pollen <- read.csv(system.file("extdata", 
#'                                           "Modern_Pollen_gdd_alpha_Tmin.csv", 
#'                                           package = "fxTWAPLS", 
#'                                           mustWork = TRUE))
#'     point <- modern_pollen[, c("Long", "Lat")]
#'     dist <- get_distance(point, cpus = 1)
#' }
#' 
#' @seealso \code{\link{get_pseudo}}
get_distance <- function(point, cpus = 4, test_mode = FALSE, test_it = 5) {
  colnames(point) <- c("Long", "Lat")
  tictoc::tic("Distance between points")
  
  # Check the number of CPUs does not exceed the availability
  avail_cpus <- parallel::detectCores() - 1
  cpus <- ifelse(cpus > avail_cpus, avail_cpus, cpus)
  
  # Start parallel backend
  cl <- parallel::makeCluster(cpus, setup_strategy = "sequential")
  doParallel::registerDoParallel(cl)
  
  # Load binary operator for backend
  `%dopar%` <- foreach::`%dopar%`
  
  # Create list of indices to loop through
  idx <- 1:nrow(point)
  # Reduce the list of indices, if test_mode = TRUE
  if (test_mode) {
    idx <- 1:test_it
  }
  dist <- foreach::foreach(i = idx,
                           .combine = rbind) %dopar% {
                             tmp <- rep(0, nrow(point))
                             lon1 <- point[i, "Long"]
                             lat1 <- point[i, "Lat"]
                             for (j in 1:nrow(point)) {
                               lon2 <- point[j, "Long"]
                               lat2 <- point[j, "Lat"]
                               tmp[j] <- geosphere::distm(c(lon1, lat1),
                                                          c(lon2, lat2),
                                                          fun = geosphere::distHaversine)
                             }
                             tmp
                           }
  
  parallel::stopCluster(cl) # Stop cluster
  tictoc::toc()
  return(dist)
}

#' Plot the training results
#' 
#' @param train_output Training output, can be the output of WAPLS, WAPLS with 
#'     fx correction, TWAPLS, or TWAPLS with fx correction
#' @param col choose which column of the fitted value to plot, in other words, 
#'     how many number of components you want to use
#'     
#' @export
#' 
# @examples
plot_train <- function(train_output, col) {
  x <- train_output[["x"]]
  fitted <- train_output[["fit"]][, col]
  plotdata <- cbind.data.frame(x, fitted)
  
  max <- max(fitted, x)
  min <- min(fitted, x)

  # plot the fitted curve, the black line is the 1:1 line, the red line is the 
  # linear regression line to fitted and x, which shows the overall compression
  ggplot2::ggplot(plotdata, aes(x, fitted)) + 
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
    ggplot2::geom_abline(slope = 1, intercept = 0) + 
    ggplot2::xlim(min, max) + ggplot2::ylim(min, max) +
    ggplot2::geom_smooth(method = 'lm',
                         formula = y ~ x,
                         color = 'red')
}

#' Plot the residuals of the training results
#' 
#' @param train_output Training output, can be the output of WAPLS, WAPLS with 
#'     fx correction, TWAPLS, or TWAPLS with fx correction
#' @param col choose which column of the fitted value to plot, in other words, 
#'     how many number of components you want to use
#' 
#' @export
#' 
# @examples
plot_residuals <- function(train_output, col) {
  x <- train_output[["x"]]
  residuals <- train_output[["fit"]][, col] - train_output[["x"]]
  plotdata <- cbind.data.frame(x, residuals)
  
  maxr <- max(abs(residuals))

  # plot the residuals, the black line is 0 line, the red line is the locally 
  # estimated scatterplot smoothing, which shows the degree of local compression
  ggplot2::ggplot(plotdata, aes(x, residuals)) + 
    ggplot2::geom_point(size = 0.4) + ggplot2::theme_bw() +
    ggplot2::geom_abline(slope = 0, intercept = 0) + 
    ggplot2::xlim(min(x), max(x)) + ggplot2::ylim(-maxr, maxr) +
    ggplot2::geom_smooth(method = 'loess',
                         color = 'red',
                         se = FALSE)
}