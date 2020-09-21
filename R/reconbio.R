if(!require(matrixStats)){install.packages("matrixStats");library(matrixStats)}

# Funtions to get the density of x --------------------------------------
fx<-function(x,bin){
  pbin<-round((max(x)-min(x))/bin,digits = 0)
  bin<-(max(x)-min(x))/pbin
  hist<-hist(x,breaks = seq(min(x),max(x),by=bin))
  xbin<-seq(min(x)+bin/2,max(x)-bin/2,by=bin)
  counts<-hist[["counts"]]
  fx<-rep(NA,length(x))
  for(i in 1:length(x)){
    fx[i]<-counts[which.min(abs(x[i]-xbin))]
  }
  if(any(fx==0)){print("Some x have a count of 0!")}
  plot(fx~x)
  return(fx)
}
# Define WAPLS and TWAPLS training funtions ----------------------------------------
# fit represents the fitted value
WAPLS.w<-function(modern_taxa,modern_climate,nPLS=5,usefx=FALSE,fx=NA){
  
  #Step 0. Centre the environmental variable by subtracting the weighted mean
  x<-modern_climate
  y<-modern_taxa
  y<-as.matrix(y)
  nc<-ncol(modern_taxa);nr<-nrow(modern_taxa)
  Ytottot<-sum(y)
  sumk_yik<-rowSums(y)
  sumi_yik<-colSums(y)
  
  #Define some matrix to store the values
  u<-matrix(NA,nc,nPLS); #u of each component
  u_sd<-matrix(NA,nc,nPLS); #u of each component, standardized the same way as r
  optimum<-matrix(NA,nc,nPLS);  #u updated
  r<-matrix(NA,nr,nPLS);  #site score
  z<-matrix(NA,1,nPLS);  #standardize
  s<-matrix(NA,1,nPLS);  #standardize
  orth<-list();#store orthogonalization parameters
  alpha<-list();#store regression coefficients
  comp<-matrix(NA,nr,nPLS) #each component
  fit<-matrix(NA,nr,nPLS) #current estimate
  
  
  pls<-1
  #Step 1. Take the centred environmental variable(xi) as initial site scores (ri). 
  r[,pls]<-x-mean(x)
  
  #Step 2. Calculate new species scores (uk* by weighted averaging of the site scores)
  u[,pls] = t(y)%*%x / sumi_yik; #uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  
  #Step 3. Calculate new site scores (ri) by weighted averaging of the species scores
  r[,pls] = y%*%u[,pls]/ sumk_yik; #xi=sumk_yik*uk/sumk_yik; 1*nsite
  
  #Step 4. For the first axis go to Step 5.
  
  #Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[,pls]<-mean(r[,pls],na.rm=TRUE)
  s[,pls]<-sqrt(sum((r[,pls]-z[,pls])^2,na.rm=TRUE)/Ytottot)
  r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
  
  #Step 6. Take the standardized score as the new component
  comp[,pls]<-r[,pls]
  
  #Step 7. Regress the environmental variable (xJ on the components obtained so far using weights and take the fitted values as current estimates 
  if(!require(MASS)){install.packages("MASS");library(MASS)}
  if(usefx==FALSE){
    lm<-rlm(modern_climate~comp[,1:pls],weights = sumk_yik/Ytottot )
  }else{
    lm<-rlm(modern_climate~comp[,1:pls],weights = 1/fx^2  )
  }
  
  fit[,pls]<-lm[["fitted.values"]]
  alpha[[pls]]<-lm[["coefficients"]]
  u_sd[,pls]<-(u[,pls]-z[,pls])/s[,pls]
  optimum[,pls]<-alpha[[pls]][1]+u_sd[,pls]*alpha[[pls]][2]
  
  
  for(pls in 2:nPLS){
    #Go to Step 2 with the residuals of the regression as the new site scores (rJ.
    r[,pls]<-lm[["residuals"]]
    
    #Step 2. Calculate new species scores (uk* by weighted averaging of the site scores)
    u[,pls] = t(y)%*%r[,pls] / sumi_yik; #uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    
    #Step 3. Calculate new site scores (ri) by weighted averaging of the species scores
    r[,pls] = y%*%u[,pls]/ sumk_yik; #xi=sumk_yik*uk/sumk_yik; 1*nsite
    
    #Step 4. For second and higher components, make the new site scores (ri) uncorrelated with the previous components by orthogonalization (Ter Braak, 1987 : Table 5 .2b)
    v<-rep(NA,pls-1)
    for(j in 1:(pls-1)){
      fi<-r[,pls-j]
      xi<-r[,pls]
      v[pls-j]<-sum(sumk_yik*fi*xi)/Ytottot
      xinew<-xi-v[pls-j]*fi
    }
    orth[[pls]]<-v
    r[,pls]<-xinew
    
    
    
    #Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[,pls]<-mean(r[,pls],na.rm=TRUE)
    s[,pls]<-sqrt(sum((r[,pls]-z[,pls])^2,na.rm=TRUE)/Ytottot)
    r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
    
    #Step 6. Take the standardized score as the new component
    comp[,pls]<-r[,pls]
    
    #Step 7. Regress the environmental variable on the components obtained so far using weights and take the fitted values as current estimates 
    if(usefx==FALSE){
      lm<-rlm(modern_climate~comp[,1:pls],weights = sumk_yik/Ytottot )
    }else{
      lm<-rlm(modern_climate~comp[,1:pls],weights = 1/fx^2 )
    }
    
    fit[,pls]<-lm[["fitted.values"]]
    alpha[[pls]]<-lm[["coefficients"]]
    
    u_sd[,pls]<-(u[,pls]-z[,pls])/s[,pls]
    optimum[,pls]<-alpha[[pls]][1]+u_sd[,1:pls]%*%as.matrix(alpha[[pls]][2:(pls+1)])
    
  }
  
  list<-list(fit,modern_climate,colnames(modern_taxa),optimum,comp,u,z,s,orth,alpha,mean(modern_climate),nPLS)
  names(list)<-c(c("fit","x","taxon_name","optimum","comp","u","z","s","orth","alpha","meanx","nPLS"))
  return(list)
  
}
TWAPLS.w<-function(modern_taxa,modern_climate,nPLS=5,usefx=FALSE,fx=NA){
  
  #Step 0. Centre the environmental variable by subtracting the weighted mean
  x<-modern_climate
  y<-modern_taxa
  y<-as.matrix(y)
  nc<-ncol(modern_taxa);nr<-nrow(modern_taxa)
  Ytottot<-sum(y)
  sumk_yik<-rowSums(y)
  sumi_yik<-colSums(y)
  
  #Define some matrix to store the values
  u<-matrix(NA,nc,nPLS); #u of each component
  u_sd<-matrix(NA,nc,nPLS); #u of each component, standardized the same way as r
  optimum<-matrix(NA,nc,nPLS);  #u updated
  t<-matrix(NA,nc,nPLS); #tolerance
  r<-matrix(NA,nr,nPLS);  #site score
  z<-matrix(NA,1,nPLS);  #standardize
  s<-matrix(NA,1,nPLS);  #standardize
  orth<-list();#store orthogonalization parameters
  alpha<-list();#store regression coefficients
  comp<-matrix(NA,nr,nPLS) #each component
  fit<-matrix(NA,nr,nPLS) #current estimate
  
  pls<-1
  #Step 1. Take the centred environmental variable (xi) as initial site scores (ri). 
  r[,pls]<-x-mean(x)
  
  #Step 2. Calculate uk and tk
  u[,pls] = t(y)%*%x / sumi_yik; #uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
  n2<-matrix(NA,nc,1)
  for(k in 1:nc){
    t[k,pls] = sqrt(sum(y[,k]*(x-u[k,pls])^2)/sumi_yik[k])
    n2[k]<-1/sum((y[,k]/sum(y[,k]))^2)
    t[k,pls]<-t[k,pls]/sqrt(1-1/n2[k])
  }
  
  #Step 3. Calculate new site scores (ri)
  r[,pls] = (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
  
  #Step 4. For the first axis go to Step 5.
  
  #Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
  z[,pls]<-mean(r[,pls],na.rm=TRUE)
  s[,pls]<-sqrt(sum((r[,pls]-z[,pls])^2,na.rm=TRUE)/Ytottot)
  r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
  
  #Step 6. Take the standardized score as the new component
  comp[,pls]<-r[,pls]
  
  #Step 7. Regress the environmental variable on the components obtained so far using weights and take the fitted values as current estimates 
  if(!require(MASS)){install.packages("MASS");library(MASS)}
  if(usefx==FALSE){
    lm<-rlm(modern_climate~comp[,1:pls],weights = sumk_yik/Ytottot )
  }else{
    lm<-rlm(modern_climate~comp[,1:pls],weights = 1/fx^2 )
  }
  
  fit[,pls]<-lm[["fitted.values"]]
  alpha[[pls]]<-lm[["coefficients"]]
  u_sd[,pls]<-(u[,pls]-z[,pls])/s[,pls]
  optimum[,pls]<-alpha[[pls]][1]+u_sd[,pls]*alpha[[pls]][2]
  
  
  for(pls in 2:nPLS){
    #Go to Step 2 with the residuals of the regression as the new site scores (ri).
    r[,pls]<-lm[["residuals"]]
    
    #Step 2. Calculate new uk and tk
    u[,pls] = t(y)%*%r[,pls] / sumi_yik; #uk=sumi_yik*xi/sumi_yik; 1*nmodern_taxa
    n2<-matrix(NA,nc,1)
    for(k in 1:nc){
      t[k,pls] = sqrt(sum(y[,k]*(r[,pls]-u[k,pls])^2)/sumi_yik[k])
      n2[k]<-1/sum((y[,k]/sum(y[,k]))^2)
      t[k,pls]<-t[k,pls]/sqrt(1-1/n2[k])
    }
    
    #Step 3. Calculate new site scores (r;) by weighted averaging of the species scores, i.e. new
    r[,pls] = (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
    
    #Step 4. For second and higher components, make the new site scores (r;) uncorrelated with the previous components by orthogonalization (Ter Braak, 1987 : Table 5 .2b)
    v<-rep(NA,pls-1)
    for(j in 1:(pls-1)){
      fi<-r[,pls-j]
      xi<-r[,pls]
      v[pls-j]<-sum(sumk_yik*fi*xi)/Ytottot
      xinew<-xi-v[pls-j]*fi
    }
    orth[[pls]]<-v
    #plot(xinew~r[,pls]);abline(0,1)
    r[,pls]<-xinew
    
    
    
    #Step 5. Standardize the new site scores (ri) ter braak 1987 5.2.c
    z[,pls]<-mean(r[,pls],na.rm=TRUE)
    s[,pls]<-sqrt(sum((r[,pls]-z[,pls])^2,na.rm=TRUE)/Ytottot)
    r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
    
    #Step 6. Take the standardized scores as the new component
    comp[,pls]<-r[,pls]
    
    #Step 7. Regress the environmental variable (xJ on the components obtained so far using weights and take the fitted values as current estimates 
    if(usefx==FALSE){
      lm<-rlm(modern_climate~comp[,1:pls],weights = sumk_yik/Ytottot )
    }else{
      lm<-rlm(modern_climate~comp[,1:pls],weights = 1/fx^2 )
    }
    
    fit[,pls]<-lm[["fitted.values"]]
    alpha[[pls]]<-lm[["coefficients"]]
    
    u_sd[,pls]<-(u[,pls]-z[,pls])/s[,pls]
    optimum[,pls]<-alpha[[pls]][1]+u_sd[,1:pls]%*%as.matrix(alpha[[pls]][2:(pls+1)])
    
  }
  
  list<-list(fit,modern_climate,colnames(modern_taxa),optimum,comp,u,t,z,s,orth,alpha,mean(modern_climate),nPLS)
  names(list)<-c(c("fit","x","taxon_name","optimum","comp","u","t","z","s","orth","alpha","meanx","nPLS"))
  return(list)
  
}

# Define WAPLS and TWAPLS predict funtions ----------------------------------------
# fit represents the fitted value
WAPLS.predict.w<-function(WAPLSoutput,fossil_taxa){
  y<-fossil_taxa
  y<-as.matrix(y)
  nc<-ncol(fossil_taxa);nr<-nrow(fossil_taxa)
  Ytottot<-sum(y)
  sumk_yik<-rowSums(y)
  sumi_yik<-colSums(y)
  
  nPLS<-WAPLSoutput[["nPLS"]]
  meanx<-WAPLSoutput[["meanx"]]  
  u<-WAPLSoutput[["u"]]
  z<-WAPLSoutput[["z"]]
  s<-WAPLSoutput[["s"]]
  orth<-WAPLSoutput[["orth"]]
  alpha<-WAPLSoutput[["alpha"]]
  
  
  if(nc!=nrow(u)){
    print("Number of taxa doesn't match!")
  }
  if(all(colnames(fossil_taxa)==WAPLSoutput[["taxon_name"]])==FALSE){
    print("Taxa don't match!")
  }
  
  #Define some matrix to store the values
  fit<-matrix(NA,nr,nPLS)
  r<-matrix(NA,nr,nPLS)
  comp<-matrix(NA,nr,nPLS)
  
  pls=1
  #xi=sumk_yik*uk/sumk_yik; 1*nsite
  r[,pls] = y%*%u[,pls]/ sumk_yik; #xi=sumk_yik*uk/sumk_yik; 1*nsite
  #standardize the same way
  r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
  #multiply the same regression coefficients
  comp[,pls]<-r[,pls]
  fit[,1]<-alpha[[pls]][1]+comp[,pls]*alpha[[pls]][2]
  
  for(pls in 2:nPLS){
    #xi=sumk_yik*uk/sumk_yik; 1*nsite
    r[,pls] = y%*%u[,pls]/ sumk_yik; #xi=sumk_yik*uk/sumk_yik; 1*nsite
    #orthoganlization the same way
    for(j in 1:(pls-1)){
      fi<-r[,pls-j]
      xi<-r[,pls]
      xinew<-xi-orth[[pls]][pls-j]*fi
    }
    r[,pls]<-xinew
    
    #standardize the same way
    r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
    #multiply the same regression coefficients
    comp[,pls]<-r[,pls]
    fit[,pls]<-alpha[[pls]][1]+comp[,1:pls]%*%as.matrix(alpha[[pls]][2:(pls+1)])
  }
  
  list<-list(fit,nPLS)
  names(list)<-c(c("fit","nPLS"))
  return(list)
  
}
TWAPLS.predict.w<-function(TWAPLSoutput,fossil_taxa){
  y<-fossil_taxa
  y<-as.matrix(y)
  nc<-ncol(fossil_taxa);nr<-nrow(fossil_taxa)
  Ytottot<-sum(y)
  sumk_yik<-rowSums(y)
  sumi_yik<-colSums(y)
  
  nPLS<-TWAPLSoutput[["nPLS"]]
  meanx<-TWAPLSoutput[["meanx"]]  
  u<-TWAPLSoutput[["u"]]
  t<-TWAPLSoutput[["t"]]
  z<-TWAPLSoutput[["z"]]
  s<-TWAPLSoutput[["s"]]
  orth<-TWAPLSoutput[["orth"]]
  alpha<-TWAPLSoutput[["alpha"]]
  
  
  if(nc!=nrow(u)){
    print("Number of taxa doesn't match!")
  }
  if(all(colnames(fossil_taxa)==TWAPLSoutput[["taxon_name"]])==FALSE){
    print("Taxa don't match!")
  }
  
  #Define some matrix to store the values
  fit<-matrix(NA,nr,nPLS)
  r<-matrix(NA,nr,nPLS)
  comp<-matrix(NA,nr,nPLS)
  
  pls=1
  #xi=sumk_yik*uk/sumk_yik; 1*nsite
  r[,pls] = (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
  #standardize the same way
  r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
  #multiply the same regression coefficients
  comp[,pls]<-r[,pls]
  fit[,1]<-alpha[[pls]][1]+comp[,pls]*alpha[[pls]][2]
  
  for(pls in 2:nPLS){
    #xi=sumk_yik*uk/sumk_yik; 1*nsite
    r[,pls] = (y%*%(u[,pls]/t[,pls]^2))/(y%*%(1/t[,pls]^2)); #xi; 1*nsite
    #orthoganlization the same way
    for(j in 1:(pls-1)){
      fi<-r[,pls-j]
      xi<-r[,pls]
      xinew<-xi-orth[[pls]][pls-j]*fi
    }
    r[,pls]<-xinew
    
    #standardize the same way
    r[,pls]<-(r[,pls]-z[,pls])/s[,pls]
    #multiply the same regression coefficients
    comp[,pls]<-r[,pls]
    fit[,pls]<-alpha[[pls]][1]+comp[,1:pls]%*%as.matrix(alpha[[pls]][2:(pls+1)])
  }
  
  list<-list(fit,nPLS)
  names(list)<-c(c("fit","nPLS"))
  return(list)
  
}

#Sample specific errors
sse.sample<-function(modern_taxa,modern_climate,fossil_taxa,trainfun,predictfun,nboot,nPLS,nsig,usefx,fx){
  # Make NA filled list of names for each taxon 
  predr<-rep(NA, nrow(fossil_taxa))
  
  # Make many sets, run WAPLS 
  xboot<-replicate(nboot, {                                         # Do this n times...
    k<-sample(1:nrow(modern_taxa), size=nrow(modern_taxa), replace=TRUE) # Make list of row numbers by sampling with replacement
    modern_taxa<-modern_taxa[k,]                                         # Reorganise modern_taxa obs in k order
    modern_climate<-modern_climate[k]
    col_not0<-which(colSums(modern_taxa)>0)
    modern_taxa<-modern_taxa[,col_not0]                            # Strip out zero-sum cols
    
    mod<-trainfun(modern_taxa, modern_climate)                                # Apply WAPLS, with modern_climate also in k order
    pred<-as.data.frame(predictfun(mod,fossil_taxa[,col_not0]))[,nsig]      # Make reconstruction
    predr<-pred
  })
  
  avg.xboot<-rowMeans(xboot, na.rm=TRUE)               
  v1<-boot.mean.square<- rowMeans((xboot-avg.xboot)^2 , na.rm=TRUE )  
  return(sqrt(v1))
  
}


# Leave one out cross validation as rioja---------------------------------------------
cv.w<-function(modern_taxa,modern_climate,nPLS=5,trainfun,predictfun,usefx=FALSE,fx=NA){
  x<-modern_climate
  y<-modern_taxa

  all.cv.out<-data.frame(matrix(nrow=nrow(modern_taxa),ncol=nPLS+1))
  for(i in 1:length(x)) {
    if(i%%100==0){print(i)} #show progress of the calculation
    fit<-trainfun(y[-i,],x[-i],nPLS,usefx,fx[-i])
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

# Leave out cross validation with pseudo (geographically and climatically close) sites removed from the training set ---------------------------------------------
#Get the ones which are both geographically and climatically close, and could therefore results in pseudo-replication
get_pseduo<-function(dist,x){
  pseduo<-as.list(matrix(nrow=nrow(dist)))
  for(i in 1:nrow(dist)){
    if(i%%100==0){print(i)}
    pseduo[[i]]<-which(dist[i,]<50000&abs(x-x[i])<0.02*(max(x)-min(x)))
    
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
#Get and plot the training results -------------------------
train.sig<-function(modern_climate,fit,nsig, fit_t,nsig_t, fit_f,nsig_f, fit_tf,nsig_tf){
  train<-cbind.data.frame(modern_climate,fit[["fit"]][,nsig],fit_t[["fit"]][,nsig_t],
                         fit_f[["fit"]][,nsig_f],fit_tf[["fit"]][,nsig_tf])
  colnames(train)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  return(train)
}
plot.train.sig<-function(train_sig_Tmin,train_sig_gdd,train_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #MTCO
  max_Tmin<-max(train_sig_Tmin);min_Tmin<-min(train_sig_Tmin)
  p1_Tmin<-ggplot(train_sig_Tmin,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(a)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("WA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
    
  p2_Tmin<-ggplot(train_sig_Tmin,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(b)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("TWA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
    
  p3_Tmin<-ggplot(train_sig_Tmin,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(c)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("WA-PLS with fx  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
    
  p4_Tmin<-ggplot(train_sig_Tmin,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_Tmin, x =min_Tmin,label="(d)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression(paste("TWA-PLS with fx  MTCO"," ", (degree~C))), x = expression(paste("MTCO"," ", (degree~C))))+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_Tmin,max_Tmin),breaks=c(-50,-25,0,25),labels = function(x) sprintf("%g", x))
    
  

  #GDD0
  max_gdd<-max(train_sig_gdd);min_gdd<-min(train_sig_gdd)
  p1_gdd<-ggplot(train_sig_gdd,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(e)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('WA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
    
  p2_gdd<-ggplot(train_sig_gdd,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(f)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('TWA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
    
  p3_gdd<-ggplot(train_sig_gdd,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(g)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('WA-PLS with fx  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
    
  p4_gdd<-ggplot(train_sig_gdd,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_gdd, x =min_gdd,label="(h)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= bquote('TWA-PLS with fx  '~ GDD[0]), x = bquote(GDD[0]))+
    scale_y_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_gdd,max_gdd),breaks=c(0,3000,6000,9000),labels = function(x) sprintf("%g", x))
    
  
  #alpha
  max_alpha<-max(train_sig_alpha);min_alpha<-min(train_sig_alpha)
  p1_alpha<-ggplot(train_sig_alpha,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(i)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("WA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
    
  p2_alpha<-ggplot(train_sig_alpha,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(j)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("TWA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
    
  p3_alpha<-ggplot(train_sig_alpha,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(k)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("WA-PLS with fx   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
    
  p4_alpha<-ggplot(train_sig_alpha,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_alpha, x =min_alpha,label="(l)")+geom_abline(slope=1, intercept=0)+
    geom_smooth(method='lm', formula= y~x,color='red')+
    labs(y= expression("TWA-PLS with fx   "*alpha), x = expression(alpha))+
    scale_y_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_alpha,max_alpha),breaks=c(0,0.5,1,1.5),labels = function(x) sprintf("%0.2f", x))
    
  p<-ggarrange(p1_Tmin,p1_gdd,p1_alpha,
               p2_Tmin,p2_gdd,p2_alpha,
               p3_Tmin,p3_gdd,p3_alpha,
               p4_Tmin,p4_gdd,p4_alpha,  ncol = 3)
  ggsave(file="Training results.jpeg",p,width=9.5,height=10)
  
}
plot.train.residual.sig<-function(train_sig_Tmin,train_sig_gdd,train_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  
  #Get the residuals
  train_residual_sig_Tmin<-cbind.data.frame(train_sig_Tmin[,"x"],train_sig_Tmin[,-1]-train_sig_Tmin[,"x"])
  colnames(train_residual_sig_Tmin)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  train_residual_sig_gdd<-cbind.data.frame(train_sig_gdd[,"x"],train_sig_gdd[,-1]-train_sig_gdd[,"x"])
  colnames(train_residual_sig_gdd)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  train_residual_sig_alpha<-cbind.data.frame(train_sig_alpha[,"x"],train_sig_alpha[,-1]-train_sig_alpha[,"x"])
  colnames(train_residual_sig_alpha)<-c("x","WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  
  #MTCO
  max_x_Tmin<-max(train_residual_sig_Tmin[,1]);min_x_Tmin<-min(train_residual_sig_Tmin[,1])
  max_y_Tmin<-max(train_residual_sig_Tmin[,-1]);min_y_Tmin<-min(train_residual_sig_Tmin[,-1])
  
  p1_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(a)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("WA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p2_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(b)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("TWA-PLS  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p3_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(c)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("WA-PLS with fx  MTCO"," ", (degree~C))), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  p4_Tmin<-ggplot(train_residual_sig_Tmin,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_Tmin, x = min_x_Tmin,label="(d)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression(paste("TWA-PLS with fx  MTCO"," ", (degree~C))), x = expression(paste("MTCO"," ", (degree~C))))+
    scale_y_continuous(limits=c(min_y_Tmin,max_y_Tmin),breaks=c(-30,-20,-10,0,10,20,30),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_Tmin,max_x_Tmin),breaks=c(-40,-30,-20,-10,0,10,20),labels = function(x) sprintf("%g", x))
  
  
  
  #GDD0
  max_x_gdd<-max(train_residual_sig_gdd[,1]);min_x_gdd<-min(train_residual_sig_gdd[,1])
  max_y_gdd<-max(train_residual_sig_gdd[,-1]);min_y_gdd<-min(train_residual_sig_gdd[,-1])
  
  p1_gdd<-ggplot(train_residual_sig_gdd,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(e)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('WA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p2_gdd<-ggplot(train_residual_sig_gdd,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(f)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('TWA-PLS  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p3_gdd<-ggplot(train_residual_sig_gdd,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(g)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('WA-PLS with fx  '~ GDD[0]), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  p4_gdd<-ggplot(train_residual_sig_gdd,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_gdd, x =min_x_gdd,label="(h)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= bquote('TWA-PLS with fx  '~ GDD[0]), x = bquote(GDD[0]))+
    scale_y_continuous(limits=c(min_y_gdd,max_y_gdd),breaks=c(-4000,-3000,-2000,-1000,0,1000,2000,3000,4000),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(limits=c(min_x_gdd,max_x_gdd),breaks=c(0,2000,4000,6000,8000),labels = function(x) sprintf("%g", x))
  
  
  #alpha
  max_x_alpha<-max(train_residual_sig_alpha[,1]);min_x_alpha<-min(train_residual_sig_alpha[,1])
  max_y_alpha<-max(train_residual_sig_alpha[,-1]);min_y_alpha<-min(train_residual_sig_alpha[,-1])
  
  p1_alpha<-ggplot(train_residual_sig_alpha,aes(x,WAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text",y= max_y_alpha, x =min_x_alpha,label="(i)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("WA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p2_alpha<-ggplot(train_residual_sig_alpha,aes(x,TWAPLS))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(j)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("TWA-PLS   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p3_alpha<-ggplot(train_residual_sig_alpha,aes(x,WAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(k)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("WA-PLS with fx   "*alpha), x = NULL)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p4_alpha<-ggplot(train_residual_sig_alpha,aes(x,TWAPLS.fx))+geom_point(size=0.4)+theme_bw()+
    annotate("text", y= max_y_alpha, x =min_x_alpha,label="(l)")+geom_abline(slope=0, intercept=0)+
    geom_smooth(method='loess',color='red', se = FALSE)+
    labs(y= expression("TWA-PLS with fx   "*alpha), x = expression(alpha))+
    scale_y_continuous(limits=c(min_y_alpha,max_y_alpha),breaks=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(limits=c(min_x_alpha,max_x_alpha),breaks=c(0,0.2,0.4,0.6,0.8,1.0,1.2),labels = function(x) sprintf("%0.2f", x))
  
  p<-ggarrange(p1_Tmin,p1_gdd,p1_alpha,
               p2_Tmin,p2_gdd,p2_alpha,
               p3_Tmin,p3_gdd,p3_alpha,
               p4_Tmin,p4_gdd,p4_alpha,  ncol = 3)
  ggsave(file="Training residuals.jpeg",p,width=9.5,height=10)
  
}

#Get and plot the core results -------------------------
core.sig<-function(fossil,nsig, fossil_t,nsig_t, fossil_f,nsig_f, fossil_tf,nsig_tf){
  core<-cbind.data.frame(fossil[["fit"]][,nsig],fossil_t[["fit"]][,nsig_t],
                          fossil_f[["fit"]][,nsig_f],fossil_tf[["fit"]][,nsig_tf])
  colnames(core)<-c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")
  return(core)
}
plot.core.sig<-function(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
  
  xbreak<-2000*c(seq(0,6))
  
  #MTCO
  plotsite_Tmin<-core_sig_Tmin[which(core_sig_Tmin$Site==sitename),]
  max_Tmin<-max(plotsite_Tmin[,-c(1,2)],na.rm=TRUE);max_Tmin<-max(max_Tmin,core_modern_Tmin,zero_Tmin)
  min_Tmin<-min(plotsite_Tmin[,-c(1,2)],na.rm=TRUE);min_Tmin<-min(min_Tmin,core_modern_Tmin,zero_Tmin)
 
  p_Tmin<-ggplot(plotsite_Tmin)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    labs(y= expression(paste("Reconstructed MTCO"," ", (degree~C))), x = NULL)+
    annotate("text", y= max_Tmin, x = (-1000),label="(a)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_Tmin,yend = core_modern_Tmin)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_Tmin,linetype="dashed")
  
  #GDD0
  plotsite_gdd<-core_sig_gdd[which(core_sig_gdd$Site==sitename),]
  max_gdd<-max(plotsite_gdd[,-c(1,2)],na.rm=TRUE);max_gdd<-max(max_gdd,core_modern_gdd,zero_gdd)
  min_gdd<-min(plotsite_gdd[,-c(1,2)],na.rm=TRUE);min_gdd<-min(min_gdd,core_modern_gdd,zero_gdd)
  
  p_gdd<-ggplot(plotsite_gdd)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    labs(y= bquote('Reconstructed '~ GDD[0]), x = NULL)+
    annotate("text", y= max_gdd, x = (-1000),label="(b)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_gdd,yend = core_modern_gdd)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_gdd,linetype="dashed")
  
  #alpha
  plotsite_alpha<-core_sig_alpha[which(core_sig_alpha$Site==sitename),]
  max_alpha<-max(plotsite_alpha[,-c(1,2)],na.rm=TRUE);max_alpha<-max(max_alpha,core_modern_alpha,zero_alpha)
  min_alpha<-min(plotsite_alpha[,-c(1,2)],na.rm=TRUE);min_alpha<-min(min_alpha,core_modern_alpha,zero_alpha)
  
  p_alpha<-ggplot(plotsite_alpha)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    labs(y= expression("Reconstructed "*alpha), x = "Age (yr BP)")+
    annotate("text", y= max_alpha, x = (-1000),label="(c)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_alpha,yend = core_modern_alpha)+
    theme_bw()+theme(axis.title.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_alpha,linetype="dashed")
  
  #Range of reconstruction
  range_WAPLS<-c(max(plotsite_Tmin[,"WAPLS"],na.rm=TRUE)-min(plotsite_Tmin[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_gdd[,"WAPLS"],na.rm=TRUE)-min(plotsite_gdd[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_alpha[,"WAPLS"],na.rm=TRUE)-min(plotsite_alpha[,"WAPLS"],na.rm=TRUE))
  #range_TWAPLS<-c(max(plotsite_Tmin[,"TWAPLS"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS"],na.rm=TRUE),
               #  max(plotsite_gdd[,"TWAPLS"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS"],na.rm=TRUE),
                # max(plotsite_alpha[,"TWAPLS"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS"],na.rm=TRUE))
 # range_WAPLS.fx<-c(max(plotsite_Tmin[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"WAPLS.fx"],na.rm=TRUE),
                    # max(plotsite_gdd[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"WAPLS.fx"],na.rm=TRUE),
                    # max(plotsite_alpha[,"WAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"WAPLS.fx"],na.rm=TRUE))
  range_TWAPLS.fx<-c(max(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE),
                 max(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE),
                 max(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE))
  #range<-rbind.data.frame(range_WAPLS,range_TWAPLS,range_WAPLS.fx,range_TWAPLS.fx)
  #colnames(range)<-c("range_Tmin","range_gdd","range_alpha");  
  #range$method<-c("WA-PLS","TWA-PLS","WA-PLS with fx","TWA-PLS with fx");
  #range$method<-factor(range$method,levels=c("WA-PLS","TWA-PLS","WA-PLS with fx","TWA-PLS with fx"))
  range<-rbind.data.frame(range_WAPLS,range_TWAPLS.fx)
  colnames(range)<-c("range_Tmin","range_gdd","range_alpha");  
  range$method<-c("WA-PLS","TWA-PLS with fx");
  range$method<-factor(range$method,levels=c("WA-PLS","TWA-PLS with fx"))
  
  p_Tmin_range<-ggplot(range)+geom_col(aes(method,range_Tmin),width=0.2)+ylim(0,1.2*max(range$range_Tmin))+
    geom_text(aes(method,range_Tmin),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_Tmin,digits = 2)))+
    ylab(expression(paste("Range of MTCO"," ", (degree~C))))+
    annotate("text", y= 1.2*max(range$range_Tmin), x = 0.6,label="(d)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_gdd_range<-ggplot(range)+geom_col(aes(method,range_gdd),width=0.2)+ylim(0,1.2*max(range$range_gdd))+
    geom_text(aes(method,range_gdd),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_gdd,digits = 0)))+
    ylab(bquote('Range of'~ GDD[0]))+annotate("text", y= 1.2*max(range$range_gdd), x = 0.6,label="(e)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_alpha_range<-ggplot(range)+geom_col(aes(method,range_alpha),width=0.2)+xlab("Method")+ylim(0,1.2*max(range$range_alpha))+
    geom_text(aes(method,range_alpha),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_alpha,digits = 2)))+
    labs(y=expression("Range of "*alpha),x=NULL)+annotate("text", y= 1.2*max(range$range_alpha), x = 0.6,label="(f)")+theme_bw()
  
  pleft<-ggarrange(p_Tmin,p_gdd,p_alpha,ncol=1)
  pright<-ggarrange(p_Tmin_range,p_gdd_range,p_alpha_range,ncol=1)
  
  p<-arrangeGrob(pleft ,pright, ncol = 7,nrow=1,layout_matrix = cbind(1,1,1,1,1,2,2))
  ggsave(file=paste(sitename,"reconstruction results using the last significant number of components.jpeg"),p,width=10,height=9)
  
}
plot.core.sig.sse<-function(sitename, core_sig_Tmin,core_modern_Tmin, core_sig_gdd,core_modern_gdd, core_sig_alpha,core_modern_alpha,zero_Tmin,zero_gdd,zero_alpha,sse_core_sig_Tmin,sse_core_sig_gdd,sse_core_sig_alpha){
  if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
  if(!require(egg)){ install.packages("egg");library(egg)}
  if(!require(gridExtra)){install.packages("gridExtra");library(gridExtra)}
  
  xbreak<-2000*c(seq(0,6))
  
  #MTCO
  plotsite_Tmin<-cbind.data.frame(core_sig_Tmin[which(core_sig_Tmin$Site==sitename),],sse_core_sig_Tmin[which(core_sig_Tmin$Site==sitename),])
  max_Tmin<-max(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_Tmin<-max(max_Tmin,core_modern_Tmin,zero_Tmin)
  min_Tmin<-min(plotsite_Tmin[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_Tmin<-min(min_Tmin,core_modern_Tmin,zero_Tmin)
  
  p_Tmin<-ggplot(plotsite_Tmin)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= expression(paste("Reconstructed MTCO"," ", (degree~C))), x = NULL)+
    annotate("text", y= max_Tmin, x = (-1000),label="(a)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_Tmin,yend = core_modern_Tmin)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_Tmin,max_Tmin),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_Tmin,linetype="dashed")
  
  #GDD0
  plotsite_gdd<-cbind.data.frame(core_sig_gdd[which(core_sig_gdd$Site==sitename),],sse_core_sig_gdd[which(core_sig_gdd$Site==sitename),])
  max_gdd<-max(plotsite_gdd[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_gdd<-max(max_gdd,core_modern_gdd,zero_gdd)
  min_gdd<-min(plotsite_gdd[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_gdd<-min(min_gdd,core_modern_gdd,zero_gdd)
  
  p_gdd<-ggplot(plotsite_gdd)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= bquote('Reconstructed '~ GDD[0]), x = NULL)+
    annotate("text", y= max_gdd, x = (-1000),label="(b)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_gdd,yend = core_modern_gdd)+
    theme_bw()+theme(axis.title.x=element_blank(),axis.text.x=element_blank())+
    scale_y_continuous(limits=c(min_gdd,max_gdd),labels = function(x) sprintf("%g", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_gdd,linetype="dashed")
  
  #alpha
  plotsite_alpha<-cbind.data.frame(core_sig_alpha[which(core_sig_alpha$Site==sitename),],sse_core_sig_alpha[which(core_sig_alpha$Site==sitename),])
  max_alpha<-max(plotsite_alpha[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);max_alpha<-max(max_alpha,core_modern_alpha,zero_alpha)
  min_alpha<-min(plotsite_alpha[,c("WAPLS","TWAPLS","WAPLS.fx","TWAPLS.fx")],na.rm=TRUE);min_alpha<-min(min_alpha,core_modern_alpha,zero_alpha)
  
  p_alpha<-ggplot(plotsite_alpha)+
    geom_point(size=0.4,aes(Age.cal.BP,WAPLS))+geom_line(aes(Age.cal.BP,WAPLS))+
    geom_ribbon(aes(ymin=WAPLS-1.96*sse_WAPLS, ymax=WAPLS+1.96*sse_WAPLS,x=Age.cal.BP),alpha=0.22,fill="black")+
    geom_point(col="red",size=0.4,aes(Age.cal.BP,TWAPLS.fx))+geom_line(col="red",aes(Age.cal.BP,TWAPLS.fx))+
    geom_ribbon(aes(ymin=TWAPLS.fx-sse_TWAPLS.fx, ymax=TWAPLS.fx+sse_TWAPLS.fx,x=Age.cal.BP),alpha=0.36,fill="red")+
    labs(y= expression("Reconstructed "*alpha), x = "Age (yr BP)")+
    annotate("text", y= max_alpha, x = (-1000),label="(c)")+ 
    geom_segment(x = (-100),xend=(-500), y = core_modern_alpha,yend = core_modern_alpha)+
    theme_bw()+theme(axis.title.x=element_blank())+
    scale_y_continuous(limits=c(min_alpha,max_alpha),labels = function(x) sprintf("%0.2f", x))+
    scale_x_continuous(breaks = xbreak,limits=c(-1000,12000))+
    geom_hline(yintercept = zero_alpha,linetype="dashed")
  
  #Range of reconstruction
  range_WAPLS<-c(max(plotsite_Tmin[,"WAPLS"],na.rm=TRUE)-min(plotsite_Tmin[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_gdd[,"WAPLS"],na.rm=TRUE)-min(plotsite_gdd[,"WAPLS"],na.rm=TRUE),
                 max(plotsite_alpha[,"WAPLS"],na.rm=TRUE)-min(plotsite_alpha[,"WAPLS"],na.rm=TRUE))
  range_TWAPLS.fx<-c(max(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_Tmin[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_gdd[,"TWAPLS.fx"],na.rm=TRUE),
                     max(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE)-min(plotsite_alpha[,"TWAPLS.fx"],na.rm=TRUE))
  range<-rbind.data.frame(range_WAPLS,range_TWAPLS.fx)
  colnames(range)<-c("range_Tmin","range_gdd","range_alpha");  
  range$method<-c("WA-PLS","TWA-PLS with fx");
  range$method<-factor(range$method,levels=c("WA-PLS","TWA-PLS with fx"))
  
  p_Tmin_range<-ggplot(range)+geom_col(aes(method,range_Tmin),width=0.2)+ylim(0,1.2*max(range$range_Tmin))+
    geom_text(aes(method,range_Tmin),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_Tmin,digits = 2)))+
    ylab(expression(paste("Range of MTCO"," ", (degree~C))))+
    annotate("text", y= 1.2*max(range$range_Tmin), x = 0.6,label="(d)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_gdd_range<-ggplot(range)+geom_col(aes(method,range_gdd),width=0.2)+ylim(0,1.2*max(range$range_gdd))+
    geom_text(aes(method,range_gdd),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_gdd,digits = 0)))+
    ylab(bquote('Range of'~ GDD[0]))+annotate("text", y= 1.2*max(range$range_gdd), x = 0.6,label="(e)")+theme_bw()+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank())
  
  p_alpha_range<-ggplot(range)+geom_col(aes(method,range_alpha),width=0.2)+xlab("Method")+ylim(0,1.2*max(range$range_alpha))+
    geom_text(aes(method,range_alpha),position = "identity",vjust=(-1),label=sprintf("%0.2f", round(range$range_alpha,digits = 2)))+
    labs(y=expression("Range of "*alpha),x=NULL)+annotate("text", y= 1.2*max(range$range_alpha), x = 0.6,label="(f)")+theme_bw()

  p1<-ggarrange(p_Tmin,p_gdd,p_alpha,ncol=1)
  p2<-ggarrange(p_Tmin_range,p_gdd_range,p_alpha_range,ncol=1)
  
  p<-arrangeGrob(p1 ,p2, ncol = 7,nrow=1,layout_matrix = cbind(1,1,1,1,1,2,2))
  ggsave(file=paste(sitename,"reconstruction results using the last significant number of components.jpeg"),p,width=10,height=9)
  
}





