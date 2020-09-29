################################################################################
#############################  Set up environment  #############################
################################################################################
rm(list = ls())
install.packages("remotes")
remotes::install_github("special-uor/fxTWAPLS")
if(!require(here)){install.packages("here");library(here)}
if(!require(tictoc)){install.packages("tictoc");library(tictoc)}

################################################################################
##################################  Training  ##################################
################################################################################
dir.create(here::here("proceedings"), FALSE)
setwd(here::here("proceedings"))
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))

taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

dir.create(here::here("training"), FALSE)
setwd(here::here("training"))

# Get the frequency of each climate variable fx
fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)

# In step 7 of training, climate variable is regressed to the components 
# obtained,
# MTCO
fit_Tmin <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_t_Tmin <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$Tmin, nPLS = 5)
fit_f_Tmin <- fxTWAPLS::WAPLS.w(taxa, 
                                modern_pollen$Tmin, 
                                nPLS = 5, 
                                usefx = TRUE, 
                                fx = fx_Tmin)
fit_tf_Tmin <- fxTWAPLS::TWAPLS.w(taxa, 
                                  modern_pollen$Tmin, 
                                  nPLS = 5, 
                                  usefx = TRUE, 
                                  fx = fx_Tmin)

# GDD0
fit_gdd <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$gdd, nPLS = 5)
fit_t_gdd <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$gdd, nPLS = 5)
fit_f_gdd <- fxTWAPLS::WAPLS.w(taxa, 
                               modern_pollen$gdd, 
                               nPLS = 4, 
                               usefx = TRUE, 
                               fx = fx_gdd)
fit_tf_gdd <- fxTWAPLS::TWAPLS.w(taxa, 
                                 modern_pollen$gdd, 
                                 nPLS = 5, 
                                 usefx = TRUE, 
                                 fx = fx_gdd)

# alpha
fit_alpha <- fxTWAPLS::WAPLS.w(taxa, modern_pollen$alpha, nPLS = 5)
fit_t_alpha <- fxTWAPLS::TWAPLS.w(taxa, modern_pollen$alpha, nPLS = 5)
fit_f_alpha <- fxTWAPLS::WAPLS.w(taxa, 
                                 modern_pollen$alpha, 
                                 nPLS = 4, 
                                 usefx = TRUE, 
                                 fx = fx_alpha)
fit_tf_alpha <- fxTWAPLS::TWAPLS.w(taxa, 
                                   modern_pollen$alpha, 
                                   nPLS = 5, 
                                   usefx = TRUE, 
                                   fx = fx_alpha)

################################################################################
######################  LOOCV same as rioja (Table S5.2)  ######################
################################################################################
dir.create(here::here("training/loocv_as_rioja"), FALSE)
setwd(here::here("training/loocv_as_rioja"))

# Set the number of CPUS
## Find the number of available CPUs with the following command
print(parallel::detectCores())
CPUS <- 15

# MTCO
tictoc::tic("MTCO")
tictoc::tic("Tmin")
cv_Tmin2 <- fxTWAPLS::cv.w(taxa,
                           modern_pollen$Tmin,
                           nPLS = 5,
                           fxTWAPLS::WAPLS.w,
                           fxTWAPLS::WAPLS.predict.w,
                           cpus = CPUS)
write.csv(cv_Tmin, "cv_Tmin.csv")
tictoc::toc()
tictoc::tic("t_Tmin")
cv_t_Tmin <- fxTWAPLS::cv.w(taxa, 
                            modern_pollen$Tmin, 
                            nPLS = 5, 
                            fxTWAPLS::TWAPLS.w, 
                            fxTWAPLS::TWAPLS.predict.w, 
                            cpus = CPUS)
write.csv(cv_t_Tmin, "cv_t_Tmin.csv")
tictoc::toc()
tictoc::tic("f_Tmin")
cv_f_Tmin <- fxTWAPLS::cv.w(taxa,
                            modern_pollen$Tmin,
                            nPLS = 5,
                            fxTWAPLS::WAPLS.w,
                            fxTWAPLS::WAPLS.predict.w,
                            usefx = TRUE,
                            fx = fx_Tmin,
                            cpus = CPUS)
write.csv(cv_f_Tmin, "cv_f_Tmin.csv")
tictoc::toc()
tictoc::tic("tf_Tmin")
cv_tf_Tmin <- fxTWAPLS::cv.w(taxa,
                             modern_pollen$Tmin,
                             nPLS = 5,
                             fxTWAPLS::TWAPLS.w,
                             fxTWAPLS::TWAPLS.predict.w,
                             usefx = TRUE,
                             fx = fx_Tmin,
                             cpus = CPUS)
write.csv(cv_tf_Tmin, "cv_tf_Tmin.csv")
tictoc::toc()
rand_Tmin <- fxTWAPLS::rand.t.test.w(cv_Tmin, n.perm = 999)
write.csv(rand_Tmin, "rand_Tmin.csv")
rand_t_Tmin <- fxTWAPLS::rand.t.test.w(cv_t_Tmin, n.perm = 999)
write.csv(rand_t_Tmin, "rand_t_Tmin.csv")
rand_f_Tmin <- fxTWAPLS::rand.t.test.w(cv_f_Tmin, n.perm = 999)
write.csv(rand_f_Tmin, "rand_f_Tmin.csv")
rand_tf_Tmin <- fxTWAPLS::rand.t.test.w(cv_tf_Tmin, n.perm = 999)
write.csv(rand_tf_Tmin, "rand_tf_Tmin.csv")
tictoc::toc()

## Load data
cv_Tmin <- read.csv("cv_Tmin.csv", row.names = 1)
cv_t_Tmin <- read.csv("cv_t_Tmin.csv", row.names = 1)
cv_f_Tmin <- read.csv("cv_f_Tmin.csv", row.names = 1)
cv_tf_Tmin <- read.csv("cv_tf_Tmin.csv", row.names = 1)
rand_Tmin <- read.csv("rand_Tmin.csv", row.names = 1)
rand_t_Tmin <- read.csv("rand_t_Tmin.csv", row.names = 1)
rand_f_Tmin <- read.csv("rand_f_Tmin.csv", row.names = 1)
rand_tf_Tmin <- read.csv("rand_tf_Tmin.csv", row.names = 1)

# GDD0
tictoc::tic("GDD0")
tictoc::tic("gdd")
cv_gdd <- fxTWAPLS::cv.w(taxa, 
                         modern_pollen$gdd, 
                         nPLS = 5, 
                         fxTWAPLS::WAPLS.w, 
                         fxTWAPLS::WAPLS.predict.w, 
                         cpus = CPUS)
write.csv(cv_gdd, "cv_gdd.csv")
tictoc::toc()
tictoc::tic("t_gdd")
cv_t_gdd <- fxTWAPLS::cv.w(taxa, 
                           modern_pollen$gdd, 
                           nPLS = 5, 
                           fxTWAPLS::TWAPLS.w, 
                           fxTWAPLS::TWAPLS.predict.w, 
                           cpus = CPUS)
write.csv(cv_t_gdd, "cv_t_gdd.csv")
tictoc::toc()
tictoc::tic("f_gdd")
cv_f_gdd <- fxTWAPLS::cv.w(taxa,
                           modern_pollen$gdd,
                           nPLS = 4,
                           fxTWAPLS::WAPLS.w,
                           fxTWAPLS::WAPLS.predict.w,
                           usefx = TRUE,
                           fx = fx_gdd,
                           cpus = CPUS)
write.csv(cv_f_gdd, "cv_f_gdd.csv")
tictoc::toc()
tictoc::tic("tf_gdd")
cv_tf_gdd <- fxTWAPLS::cv.w(taxa,
                            modern_pollen$gdd,
                            nPLS = 5,
                            fxTWAPLS::TWAPLS.w,
                            fxTWAPLS::TWAPLS.predict.w,
                            usefx = TRUE,
                            fx = fx_gdd,
                            cpus = CPUS)
write.csv(cv_tf_gdd, "cv_tf_gdd.csv")
tictoc::toc()

rand_gdd <- fxTWAPLS::rand.t.test.w(cv_gdd, n.perm = 999)
write.csv(rand_gdd, "rand_gdd.csv")
rand_t_gdd <- fxTWAPLS::rand.t.test.w(cv_t_gdd, n.perm = 999)
write.csv(rand_t_gdd, "rand_t_gdd.csv")
rand_f_gdd <- fxTWAPLS::rand.t.test.w(cv_f_gdd, n.perm = 999)
write.csv(rand_f_gdd, "rand_f_gdd.csv")
rand_tf_gdd <- fxTWAPLS::rand.t.test.w(cv_tf_gdd, n.perm = 999)
write.csv(rand_tf_gdd, "rand_tf_gdd.csv")
tictoc::toc()

## Load data
cv_gdd <- read.csv("cv_gdd.csv", row.names = 1)
cv_t_gdd <- read.csv("cv_t_gdd.csv", row.names = 1)
cv_f_gdd <- read.csv("cv_f_gdd.csv", row.names = 1)
cv_tf_gdd <- read.csv("cv_tf_gdd.csv", row.names = 1)
rand_gdd <- read.csv("rand_gdd.csv", row.names = 1)
rand_t_gdd <- read.csv("rand_t_gdd.csv", row.names = 1)
rand_f_gdd <- read.csv("rand_f_gdd.csv", row.names = 1)
rand_tf_gdd <- read.csv("rand_tf_gdd.csv", row.names = 1)

# alpha
tictoc::tic("alpha") # global
tictoc::tic("alpha") # local
cv_alpha <- fxTWAPLS::cv.w(taxa,
                           modern_pollen$alpha,
                           nPLS = 5,
                           fxTWAPLS::WAPLS.w,
                           fxTWAPLS::WAPLS.predict.w,
                           cpus = CPUS)
write.csv(cv_alpha, "cv_alpha.csv")
tictoc::toc()
tictoc::tic("t_alpha")
cv_t_alpha <- fxTWAPLS::cv.w(taxa,
                             modern_pollen$alpha,
                             nPLS = 5,
                             fxTWAPLS::TWAPLS.w,
                             fxTWAPLS::TWAPLS.predict.w,
                             cpus = CPUS)
write.csv(cv_t_alpha, "cv_t_alpha.csv")
tictoc::toc()
tictoc::tic("f_alpha")
cv_f_alpha <- fxTWAPLS::cv.w(taxa,
                             modern_pollen$alpha,
                             nPLS = 4,
                             fxTWAPLS::WAPLS.w,
                             fxTWAPLS::WAPLS.predict.w,
                             usefx = TRUE,
                             fx = fx_alpha,
                             cpus = CPUS)
write.csv(cv_f_alpha, "cv_f_alpha.csv")
tictoc::toc()
tictoc::tic("tf_alpha")
cv_tf_alpha <- fxTWAPLS::cv.w(taxa,
                              modern_pollen$alpha,
                              nPLS = 5,
                              fxTWAPLS::TWAPLS.w,
                              fxTWAPLS::TWAPLS.predict.w,
                              usefx = TRUE,
                              fx = fx_alpha,
                              cpus = CPUS)
write.csv(cv_tf_alpha, "cv_tf_alpha.csv")
tictoc::toc()

rand_alpha <- fxTWAPLS::rand.t.test.w(cv_alpha, n.perm = 999)
write.csv(rand_alpha, "rand_alpha.csv")
rand_t_alpha <- fxTWAPLS::rand.t.test.w(cv_t_alpha, n.perm = 999)
write.csv(rand_t_alpha, "rand_t_alpha.csv")
rand_f_alpha <- fxTWAPLS::rand.t.test.w(cv_f_alpha, n.perm = 999)
write.csv(rand_f_alpha, "rand_f_alpha.csv")
rand_tf_alpha <- fxTWAPLS::rand.t.test.w(cv_tf_alpha, n.perm = 999)
write.csv(rand_tf_alpha, "rand_tf_alpha.csv")
tictoc::toc()

## Load data
cv_alpha <- read.csv("cv_alpha.csv", row.names = 1)
cv_t_alpha <- read.csv("cv_t_alpha.csv", row.names = 1)
cv_f_alpha <- read.csv("cv_f_alpha.csv", row.names = 1)
cv_tf_alpha <- read.csv("cv_tf_alpha.csv", row.names = 1)
rand_alpha <- read.csv("rand_alpha.csv", row.names = 1)
rand_t_alpha <- read.csv("rand_t_alpha.csv", row.names = 1)
rand_f_alpha <- read.csv("rand_f_alpha.csv", row.names = 1)
rand_tf_alpha <- read.csv("rand_tf_alpha.csv", row.names = 1)

# Put the results together
rand_rioja <- rbind.data.frame(rand_Tmin,
                               rand_t_Tmin,
                               rand_f_Tmin,
                               rand_tf_Tmin,
                               rand_gdd,
                               rand_t_gdd,
                               rand_f_gdd,
                               rand_tf_gdd,
                               rand_alpha,
                               rand_t_alpha,
                               rand_f_alpha,
                               rand_tf_alpha)
    
write.csv(rand_rioja, here::here("rand_rioja.csv"))

################################################################################
## Cross validation with pseudo sites removed from the training set (Table 2) ##
################################################################################
dir.create(here::here("training/cv_pseudo_removed"), FALSE)
setwd(here::here("training/cv_pseudo_removed"))

# Calculate the distance between points
point <- modern_pollen[, c("Long", "Lat")]
dist <- fxTWAPLS::get_distance(point, cpus = CPUS)
write.csv(dist, "distance.csv")

# Get the geographically and climatically close sites
dist <- read.csv("distance.csv", row.names = 1)
tictoc::tic("Pseudo")
tictoc::tic("Tmin")
pseudo_Tmin <- fxTWAPLS::get_pseudo(dist, modern_pollen$Tmin, cpus = CPUS)
tictoc::toc()
tictoc::tic("gdd")
pseudo_gdd <- fxTWAPLS::get_pseudo(dist, modern_pollen$gdd, cpus = CPUS)
tictoc::toc()
tictoc::tic("alpha")
pseudo_alpha <- fxTWAPLS::get_pseudo(dist, modern_pollen$alpha, cpus = CPUS)
tictoc::toc()

if(!require(rlist)){install.packages("rlist");library(rlist)} 
rlist::list.save(pseudo_Tmin, 'pseudo_Tmin.rdata')
rlist::list.save(pseudo_gdd, 'pseudo_gdd.rdata')
rlist::list.save(pseudo_alpha, 'pseudo_alpha.rdata')

# Pseudo removed leave out cross validation
pseudo_Tmin <- rlist::list.load('pseudo_Tmin.rdata')
pseudo_gdd <- rlist::list.load('pseudo_gdd.rdata')
pseudo_alpha <- rlist::list.load('pseudo_alpha.rdata')

# MTCO
tictoc::tic("Pseudo MTCO")
tictoc::tic("Tmin")
cv_pr_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                modern_pollen$Tmin,
                                nPLS = 5,
                                fxTWAPLS::WAPLS.w,
                                fxTWAPLS::WAPLS.predict.w,
                                pseudo_Tmin,
                                cpus = CPUS)
write.csv(cv_pr_Tmin, "cv_pr_Tmin.csv")
tictoc::toc()
tictoc::tic("t_Tmin")
cv_pr_t_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                  modern_pollen$Tmin,
                                  nPLS = 5,
                                  fxTWAPLS::TWAPLS.w,
                                  fxTWAPLS::TWAPLS.predict.w,
                                  pseudo_Tmin,
                                  cpus = CPUS)
write.csv(cv_pr_t_Tmin, "cv_pr_t_Tmin.csv")
tictoc::toc()
tictoc::tic("f_Tmin")
cv_pr_f_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                  modern_pollen$Tmin,
                                  nPLS = 5,
                                  fxTWAPLS::WAPLS.w,
                                  fxTWAPLS::WAPLS.predict.w,
                                  usefx = TRUE,
                                  fx = fx_Tmin,
                                  pseudo_Tmin,
                                  cpus = CPUS)
write.csv(cv_pr_f_Tmin, "cv_pr_f_Tmin.csv")
tictoc::toc()
tictoc::tic("tf_Tmin")
cv_pr_tf_Tmin <- fxTWAPLS::cv.pr.w(taxa,
                                   modern_pollen$Tmin,
                                   nPLS = 5,
                                   fxTWAPLS::TWAPLS.w,
                                   fxTWAPLS::TWAPLS.predict.w,
                                   usefx = TRUE,
                                   fx = fx_Tmin,
                                   pseudo_Tmin,
                                   cpus = CPUS)
write.csv(cv_pr_tf_Tmin, "cv_pr_tf_Tmin.csv")
tictoc::toc()

rand_pr_Tmin <- fxTWAPLS::rand.t.test.w(cv_pr_Tmin, n.perm = 999)
write.csv(rand_pr_Tmin, "rand_pr_Tmin.csv")
rand_pr_t_Tmin <- fxTWAPLS::rand.t.test.w(cv_pr_t_Tmin, n.perm = 999)
write.csv(rand_pr_t_Tmin, "rand_pr_t_Tmin.csv")
rand_pr_f_Tmin <- fxTWAPLS::rand.t.test.w(cv_pr_f_Tmin, n.perm = 999)
write.csv(rand_pr_f_Tmin, "rand_pr_f_Tmin.csv")
rand_pr_tf_Tmin <- fxTWAPLS::rand.t.test.w(cv_pr_tf_Tmin, n.perm = 999)
write.csv(rand_pr_tf_Tmin, "rand_pr_tf_Tmin.csv")
tictoc::toc()

# GDD0
tictoc::tic("Pseudo GDD0")
tictoc::tic("gdd")
cv_pr_gdd <- fxTWAPLS::cv.pr.w(taxa,
                               modern_pollen$gdd,
                               nPLS = 5,
                               fxTWAPLS::WAPLS.w,
                               fxTWAPLS::WAPLS.predict.w,
                               pseudo_gdd,
                               cpus = CPUS)
write.csv(cv_pr_gdd, "cv_pr_gdd.csv")
tictoc::toc()
tictoc::tic("t_gdd")
cv_pr_t_gdd <- fxTWAPLS::cv.pr.w(taxa,
                                 modern_pollen$gdd,
                                 nPLS = 5,
                                 fxTWAPLS::TWAPLS.w,
                                 fxTWAPLS::TWAPLS.predict.w,
                                 pseudo_gdd,
                                 cpus = CPUS)
write.csv(cv_pr_t_gdd, "cv_pr_t_gdd.csv")
tictoc::toc()
tictoc::tic("f_gdd")
cv_pr_f_gdd <- fxTWAPLS::cv.pr.w(taxa,
                                 modern_pollen$gdd,
                                 nPLS = 4,
                                 fxTWAPLS::WAPLS.w,
                                 fxTWAPLS::WAPLS.predict.w,
                                 usefx = TRUE,
                                 fx = fx_gdd,
                                 pseudo_gdd,
                                 cpus = CPUS)
write.csv(cv_pr_f_gdd, "cv_pr_f_gdd.csv")
tictoc::toc()
tictoc::tic("tf_gdd")
cv_pr_tf_gdd <- fxTWAPLS::cv.pr.w(taxa,
                                  modern_pollen$gdd,
                                  nPLS = 5,
                                  fxTWAPLS::TWAPLS.w,
                                  fxTWAPLS::TWAPLS.predict.w,
                                  usefx = TRUE,
                                  fx = fx_gdd,
                                  pseudo_gdd,
                                  cpus = CPUS)
write.csv(cv_pr_tf_gdd, "cv_pr_tf_gdd.csv")
tictoc::toc()

rand_pr_gdd <- fxTWAPLS::rand.t.test.w(cv_pr_gdd, n.perm = 999)
write.csv(rand_pr_gdd, "rand_pr_gdd.csv")
rand_pr_t_gdd <- fxTWAPLS::rand.t.test.w(cv_pr_t_gdd, n.perm = 999)
write.csv(rand_pr_t_gdd, "rand_pr_t_gdd.csv")
rand_pr_f_gdd <- fxTWAPLS::rand.t.test.w(cv_pr_f_gdd, n.perm = 999)
write.csv(rand_pr_f_gdd, "rand_pr_f_gdd.csv")
rand_pr_tf_gdd <- fxTWAPLS::rand.t.test.w(cv_pr_tf_gdd, n.perm = 999)
write.csv(rand_pr_tf_gdd, "rand_pr_tf_gdd.csv")
tictoc::toc()

# alpha
tictoc::tic("Pseudo alpha")
tictoc::tic("alpha")
cv_pr_alpha <- fxTWAPLS::cv.pr.w(taxa,
                                 modern_pollen$alpha,
                                 nPLS = 5,
                                 fxTWAPLS::WAPLS.w,
                                 fxTWAPLS::WAPLS.predict.w,
                                 pseudo_alpha,
                                 cpus = CPUS)
write.csv(cv_pr_alpha, "cv_pr_alpha.csv")
tictoc::toc()
tictoc::tic("t_alpha")
cv_pr_t_alpha <- fxTWAPLS::cv.pr.w(taxa,
                                   modern_pollen$alpha,
                                   nPLS = 5,
                                   fxTWAPLS::TWAPLS.w,
                                   fxTWAPLS::TWAPLS.predict.w,
                                   pseudo_alpha,
                                   cpus = CPUS)
write.csv(cv_pr_t_alpha, "cv_pr_t_alpha.csv")
tictoc::toc()
tictoc::tic("f_alpha")
cv_pr_f_alpha <- fxTWAPLS::cv.pr.w(taxa,
                                   modern_pollen$alpha,
                                   nPLS = 4,
                                   fxTWAPLS::WAPLS.w,
                                   fxTWAPLS::WAPLS.predict.w,
                                   usefx = TRUE,
                                   fx = fx_alpha,
                                   pseudo_alpha,
                                   cpus = CPUS)
write.csv(cv_pr_f_alpha, "cv_pr_f_alpha.csv")
tictoc::toc()
tictoc::tic("tf_alpha")
cv_pr_tf_alpha <- fxTWAPLS::cv.pr.w(taxa,
                                    modern_pollen$alpha,
                                    nPLS = 5,
                                    fxTWAPLS::TWAPLS.w,
                                    fxTWAPLS::TWAPLS.predict.w,
                                    usefx = TRUE,
                                    fx = fx_alpha,
                                    pseudo_alpha,
                                    cpus = CPUS)
write.csv(cv_pr_tf_alpha, "cv_pr_tf_alpha.csv")
tictoc::toc()

rand_pr_alpha <- fxTWAPLS::rand.t.test.w(cv_pr_alpha, n.perm = 999)
write.csv(rand_pr_alpha, "rand_pr_alpha.csv")
rand_pr_t_alpha <- fxTWAPLS::rand.t.test.w(cv_pr_t_alpha, n.perm = 999)
write.csv(rand_pr_t_alpha, "rand_pr_t_alpha.csv")
rand_pr_f_alpha <- fxTWAPLS::rand.t.test.w(cv_pr_f_alpha, n.perm = 999)
write.csv(rand_pr_f_alpha, "rand_pr_f_alpha.csv")
rand_pr_tf_alpha <- fxTWAPLS::rand.t.test.w(cv_pr_tf_alpha, n.perm = 999)
write.csv(rand_pr_tf_alpha, "rand_pr_tf_alpha.csv")
tictoc::toc()

# Put the results together
rand_pseudo_removed <- rbind.data.frame(rand_pr_Tmin,
                                        rand_pr_t_Tmin,
                                        rand_pr_f_Tmin,
                                        rand_pr_tf_Tmin,
                                        rand_pr_gdd,
                                        rand_pr_t_gdd,
                                        rand_pr_f_gdd,
                                        rand_pr_tf_gdd,
                                        rand_pr_alpha,
                                        rand_pr_t_alpha,
                                        rand_pr_f_alpha,
                                        rand_pr_tf_alpha)

write.csv(rand_pseudo_removed, "rand_pseudo_removed.csv")

################################################################################
################################ Reconstruction ################################
################################################################################
Holocene <- read.csv(system.file("extdata", "Holocene.csv", 
                                 package = "fxTWAPLS", 
                                 mustWork = TRUE), 
                     row.names = 1)
taxaColMin <- which(colnames(Holocene) == "Abies")
taxaColMax <- which(colnames(Holocene) == "Zygophyllaceae")
core <- Holocene[, taxaColMin:taxaColMax]

dir.create(here::here("proceedings/core"), FALSE)
setwd(here::here("proceedings/core"))

# MTCO
fossil_Tmin <- fxTWAPLS::WAPLS.predict.w(fit_Tmin, core)
fossil_t_Tmin <- fxTWAPLS::TWAPLS.predict.w(fit_t_Tmin, core)
fossil_f_Tmin <- fxTWAPLS::WAPLS.predict.w(fit_f_Tmin, core)
fossil_tf_Tmin <- fxTWAPLS::TWAPLS.predict.w(fit_tf_Tmin, core)

# GDD0
fossil_gdd <- fxTWAPLS::WAPLS.predict.w(fit_gdd, core)
fossil_t_gdd <- fxTWAPLS::TWAPLS.predict.w(fit_t_gdd, core)
fossil_f_gdd <- fxTWAPLS::WAPLS.predict.w(fit_f_gdd, core)
fossil_tf_gdd <- fxTWAPLS::TWAPLS.predict.w(fit_tf_gdd, core)

# alpha
fossil_alpha <- fxTWAPLS::WAPLS.predict.w(fit_alpha, core)
fossil_t_alpha <- fxTWAPLS::TWAPLS.predict.w(fit_t_alpha, core)
fossil_f_alpha <- fxTWAPLS::WAPLS.predict.w(fit_f_alpha, core)
fossil_tf_alpha <- fxTWAPLS::TWAPLS.predict.w(fit_tf_alpha, core)

#Get the sample specific errors, use the default nboot=100 as in rioja
# MTCO
sse_Tmin_WAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                       modern_climate = modern_pollen$Tmin,
                                       fossil_taxa = core,
                                       trainfun = fxTWAPLS::WAPLS.w,
                                       predictfun = fxTWAPLS::WAPLS.predict.w,
                                       nboot = 100,
                                       nPLS = 5,
                                       nsig = 3,
                                       usefx = FALSE,
                                       fx = NA)

sse_Tmin_TWAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                        modern_climate = modern_pollen$Tmin,
                                        fossil_taxa = core,
                                        trainfun = fxTWAPLS::TWAPLS.w,
                                        predictfun = fxTWAPLS::TWAPLS.predict.w,
                                        nboot = 100,
                                        nPLS = 5,
                                        nsig = 3,
                                        usefx = FALSE,
                                        fx = NA)

sse_Tmin_WAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                          modern_climate = modern_pollen$Tmin,
                                          fossil_taxa = core,
                                          trainfun = fxTWAPLS::WAPLS.w,
                                          predictfun = fxTWAPLS::WAPLS.predict.w,
                                          nboot = 100,
                                          nPLS = 5,
                                          nsig = 3,
                                          usefx = TRUE,
                                          fx = fx_Tmin)

sse_Tmin_TWAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                           modern_climate = modern_pollen$Tmin,
                                           fossil_taxa = core,
                                           trainfun = fxTWAPLS::TWAPLS.w,
                                           predictfun = fxTWAPLS::TWAPLS.predict.w,
                                           nboot = 100,
                                           nPLS = 5,
                                           nsig = 4,
                                           usefx = TRUE,
                                           fx = fx_Tmin)

sse_core_sig_Tmin <- cbind.data.frame(sse_Tmin_WAPLS,
                                      sse_Tmin_TWAPLS,
                                      sse_Tmin_WAPLS.fx,
                                      sse_Tmin_TWAPLS.fx)
colnames(sse_core_sig_Tmin) <- c("sse_WAPLS", 
                                 "sse_TWAPLS", 
                                 "sse_WAPLS.fx", 
                                 "sse_TWAPLS.fx")

# GDD0
sse_gdd_WAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa, 
                                      modern_climate = modern_pollen$gdd, 
                                      fossil_taxa = core, 
                                      trainfun = fxTWAPLS::WAPLS.w, 
                                      predictfun = fxTWAPLS::WAPLS.predict.w, 
                                      nboot = 100, 
                                      nPLS = 5, 
                                      nsig = 2, 
                                      usefx = FALSE, 
                                      fx = NA)

sse_gdd_TWAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                       modern_climate = modern_pollen$gdd,
                                       fossil_taxa = core,
                                       trainfun = fxTWAPLS::TWAPLS.w,
                                       predictfun = fxTWAPLS::TWAPLS.predict.w,
                                       nboot = 100,
                                       nPLS = 5,
                                       nsig = 2,
                                       usefx = FALSE,
                                       fx = NA)
sse_gdd_WAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                         modern_climate = modern_pollen$gdd,
                                         fossil_taxa = core,
                                         trainfun = fxTWAPLS::WAPLS.w,
                                         predictfun = fxTWAPLS::WAPLS.predict.w,
                                         nboot = 100,
                                         nPLS = 4,
                                         nsig = 2,
                                         usefx = TRUE,
                                         fx = fx_gdd)

sse_gdd_TWAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                          modern_climate = modern_pollen$gdd,
                                          fossil_taxa = core,
                                          trainfun = fxTWAPLS::TWAPLS.w,
                                          predictfun = fxTWAPLS::TWAPLS.predict.w,
                                          nboot = 100,
                                          nPLS = 5,
                                          nsig = 3,
                                          usefx = TRUE,
                                          fx = fx_gdd)

sse_core_sig_gdd <- cbind.data.frame(sse_gdd_WAPLS,
                                     sse_gdd_TWAPLS,
                                     sse_gdd_WAPLS.fx,
                                     sse_gdd_TWAPLS.fx)
colnames(sse_core_sig_gdd) <- c("sse_WAPLS", 
                                "sse_TWAPLS", 
                                "sse_WAPLS.fx", 
                                "sse_TWAPLS.fx")

# alpha
sse_alpha_WAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                        modern_climate = modern_pollen$alpha,
                                        fossil_taxa = core,
                                        trainfun = fxTWAPLS::WAPLS.w,
                                        predictfun = fxTWAPLS::WAPLS.predict.w,
                                        nboot = 100,
                                        nPLS = 5,
                                        nsig = 3,
                                        usefx = FALSE,
                                        fx = NA)

sse_alpha_TWAPLS <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                         modern_climate = modern_pollen$alpha,
                                         fossil_taxa = core,
                                         trainfun = fxTWAPLS::TWAPLS.w,
                                         predictfun = fxTWAPLS::TWAPLS.predict.w,
                                         nboot = 100,
                                         nPLS = 5,
                                         nsig = 4,
                                         usefx = FALSE,
                                         fx = NA)

sse_alpha_WAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                           modern_climate = modern_pollen$alpha,
                                           fossil_taxa = core,
                                           trainfun = fxTWAPLS::WAPLS.w,
                                           predictfun = fxTWAPLS::WAPLS.predict.w,
                                           nboot = 100,
                                           nPLS = 4,
                                           nsig = 2,
                                           usefx = TRUE,
                                           fx = fx_alpha)

sse_alpha_TWAPLS.fx <- fxTWAPLS::sse.sample(modern_taxa = taxa,
                                            modern_climate = modern_pollen$alpha,
                                            fossil_taxa = core,
                                            trainfun = fxTWAPLS::TWAPLS.w,
                                            predictfun = fxTWAPLS::TWAPLS.predict.w,
                                            nboot = 100,
                                            nPLS = 5,
                                            nsig = 3,
                                            usefx = TRUE,
                                            fx = fx_alpha)

sse_core_sig_alpha <- cbind.data.frame(sse_alpha_WAPLS,
                                       sse_alpha_TWAPLS,
                                       sse_alpha_WAPLS.fx,
                                       sse_alpha_TWAPLS.fx)
colnames(sse_core_sig_alpha) <- c("sse_WAPLS", 
                                  "sse_TWAPLS", 
                                  "sse_WAPLS.fx", 
                                  "sse_TWAPLS.fx")
