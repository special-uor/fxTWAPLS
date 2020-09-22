# if(!require(matrixStats)){install.packages("matrixStats");library(matrixStats)}

#' Funtion to get the density of x
#'
#' @param x 
#' @param bin 
#'
#' @return
#' @export
#'
# @examples
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

# Define WAPLS and TWAPLS training funtions fit represents the fitted value
#' WALPS training function fit
#' 
#' @importFrom stats lm
#' 
#' @param modern_taxa 
#' @param modern_climate 
#' @param nPLS 
#' @param usefx 
#' @param fx 
#'
#' @return
#' @export
#'
# @examples
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
  # if(!require(MASS)){install.packages("MASS");library(MASS)}
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

#' TWALPS training function fit
#' @importFrom stats lm
#' @param modern_taxa 
#' @param modern_climate 
#' @param nPLS 
#' @param usefx 
#' @param fx 
#'
#' @return
#' @export
#'
# @examples
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
  # if(!require(MASS)){install.packages("MASS");library(MASS)}
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

# Define WAPLS and TWAPLS predict funtions ----------------------------------------
# fit represents the fitted value
#' WAPLS predict function
#'
#' @param WAPLSoutput 
#' @param fossil_taxa 
#'
#' @return
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
#' @param TWAPLSoutput 
#' @param fossil_taxa 
#'
#' @return
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
#' @param modern_taxa 
#' @param modern_climate 
#' @param fossil_taxa 
#' @param trainfun 
#' @param predictfun 
#' @param nboot 
#' @param nPLS 
#' @param nsig 
#' @param usefx 
#' @param fx 
#'
#' @return
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
                       fx) {
  # Make NA filled list of names for each taxon 
  predr <- rep(NA, nrow(fossil_taxa))
  
  # Make many sets, run WAPLS 
  xboot <-
    replicate(nboot, { # Do this n times...
      k <- sample(1:nrow(modern_taxa),
                  size = nrow(modern_taxa),
                  replace = TRUE) # Make list of row numbers by sampling with replacement
      modern_taxa <- modern_taxa[k, ] # Reorganise modern_taxa obs in k order
      modern_climate <- modern_climate[k]
      col_not0 <- which(colSums(modern_taxa) > 0)
      modern_taxa <- modern_taxa[, col_not0] # Strip out zero-sum cols
      
      mod <- trainfun(modern_taxa, modern_climate) # Apply WAPLS, with modern_climate also in k order
      pred <- as.data.frame(predictfun(mod, fossil_taxa[, col_not0]))[, nsig] # Make reconstruction
      predr <- pred
    })
  
  avg.xboot <- rowMeans(xboot, na.rm = TRUE)
  v1 <- boot.mean.square <- rowMeans((xboot - avg.xboot) ^ 2 , na.rm = TRUE)
  return(sqrt(v1))
}

#' Leave one out cross validation as rioja
#'
#' @param modern_taxa 
#' @param modern_climate 
#' @param nPLS 
#' @param trainfun 
#' @param predictfun 
#' @param usefx 
#' @param fx 
#'
#' @return
#' @export
#'
# @examples
cv.w <- function(modern_taxa,
                 modern_climate,
                 nPLS = 5,
                 trainfun,
                 predictfun,
                 usefx = FALSE,
                 fx = NA) {
  x <- modern_climate
  y <- modern_taxa

  all.cv.out <- data.frame(matrix(nrow = nrow(modern_taxa), ncol = nPLS + 1))
  for (i in 1:length(x)) {
    if (i %% 100 == 0) {
      print(i)
    } #show progress of the calculation
    fit <- trainfun(y[-i, ], x[-i], nPLS, usefx, fx[-i])
    xnew <- predictfun(fit, y[i, ])[["fit"]]
    all.cv.out[i, ] <- data.frame(x[i], xnew)
    
  }
  for(k in 1:ncol(all.cv.out)) {
    if (k == 1) {
      colnames(all.cv.out)[k] <- "test.x"
    } else{
      colnames(all.cv.out)[k] <- paste("comp", k - 1, sep = "")
    }
  } #assign column names to all.cv.out
  
  return(all.cv.out)
}

#' Leave out cross validation with pseudo (geographically and climatically close) 
#' sites removed from the training set.
#' Get the ones which are both geographically and climatically close, and could 
#' therefore results in pseudo-replication.
#'
#' @param dist 
#' @param x 
#'
#' @return
#' @export
#'
# @examples
get_pseduo <- function(dist, x) {
  pseduo <- as.list(matrix(nrow = nrow(dist)))
  for (i in 1:nrow(dist)) {
    if (i %% 100 == 0) {
      print(i)
    }
    pseduo[[i]] <-
      which(dist[i, ] < 50000 & abs(x - x[i]) < 0.02 * (max(x) - min(x)))
  }
  return(pseduo)
}

#Pseudo removed leave out cross validation
cv.pr.w<-function(modern_taxa,modern_climate,nPLS=5,trainfun,predictfun,pseduo,usefx=FALSE,fx=NA){
  x<-modern_climate
  y<-modern_taxa
  all.cv.out<-data.frame(matrix(nrow=nrow(modern_taxa),ncol=nPLS+1))
  for(i in 1:length(x)) {
    if(i%%100==0){print(i)} #show progress of the calculation
    leave<-unlist(pseduo[i])
    fit<-trainfun(y[-leave,],x[-leave],nPLS,usefx,fx[-leave])
    xnew <- predictfun(fit, y[i,])[["fit"]]
    all.cv.out[i,]<-data.frame(x[i],xnew)
    
  }
  
  for(k in 1:ncol(all.cv.out)){
    if(k==1){
      colnames(all.cv.out)[k]<-"test.x"
    }else{
      colnames(all.cv.out)[k]<-paste("comp",k-1,sep="")
    }}#assign column names to all.cv.out
  
  return(all.cv.out)
  
}


# Random t-test ---------------------------------------------
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats rbinom
rand.t.test.w<-function(cvoutput,n.perm=999){
  ncomp<-ncol(cvoutput)-1
  output<-matrix(NA,ncomp,11)
  colnames(output)<-c("R2","Avg.Bias","Max.Bias","Min.Bias","RMSEP","delta.RMSEP","p","Compre.b0","Compre.b1","Compre.b0.se","Compre.b1.se")
  
  for(i in 1:ncomp){
    cv.x<-cvoutput[,1]
    cv.i<-cvoutput[,1+i]
    output[i,"RMSEP"]<-sqrt(mean((cv.i - cv.x)^2))
    output[i,"R2"]<-cor(cv.i,cv.x)^2
    output[i,"Avg.Bias"]<-mean(cv.i - cv.x)
    output[i,"Max.Bias"]<-max(abs(cv.i - cv.x))
    output[i,"Min.Bias"]<-min(abs(cv.i - cv.x))
    output[i,c("Compre.b0","Compre.b1")]<-summary(lm(cv.i~cv.x))[["coefficients"]][,"Estimate"]
    output[i,c("Compre.b0.se","Compre.b1.se")]<-summary(lm(cv.i~cv.x))[["coefficients"]][,"Std. Error"]
  }
  #get delta.RMSEP
  for(i in 1:ncomp){
    if(i==1){
      rmsep.null<-sqrt(mean((cv.x - mean(cv.x))^2))
      output[i,"delta.RMSEP"]<-(output[i,"RMSEP"]-rmsep.null)*100/rmsep.null
    }else{
      output[i,"delta.RMSEP"]<-(output[i,"RMSEP"]-output[i-1,"RMSEP"])*100/output[i-1,"RMSEP"]
    }
  }
  #get p-value, which describes whether using the number of components now has a significant difference than using one less
  e0 <- cv.x - mean(cv.x)
  e <- cbind(e0, cvoutput[,2:ncol(cvoutput)]-cv.x)
  
  t.res <- vector("numeric",ncomp)
  t <- vector("numeric",n.perm+1)
  t.res[] <- NA
  n <- nrow(e)
  for (i in 1:ncomp) {
    d <- e[, i]^2 - e[, i+1]^2
    t[1] <- mean(d, na.rm=TRUE)
    for (j in 1:n.perm) {
      sig <- 2 * rbinom(n, 1, 0.5) - 1
      t[j+1] <- mean(d * sig, na.rm=TRUE)
    }
    t.res[i] <- sum(t >= t[1]) / (n.perm+1)
  }
  output[,"p"]<-t.res
  
  print(output)
  return(output)
  
}