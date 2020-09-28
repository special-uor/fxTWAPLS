##################################### cv.w #####################################
times_df <- fxTWAPLS::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = fxTWAPLS::cv.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = fxTWAPLS::WAPLS.w, 
                                    predictfun = fxTWAPLS::WAPLS.predict.w)

################################### cv.pr.w ####################################
times_df <- fxTWAPLS::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = fxTWAPLS::cv.pr.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = fxTWAPLS::WAPLS.w, 
                                    predictfun = fxTWAPLS::WAPLS.predict.w,
                                    pseduo = pseduo_Tmin)

################################## get_pseduo ##################################
times_df <- fxTWAPLS::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = fxTWAPLS::get_pseduo, 
                                    plot = TRUE, 
                                    dist = dist, 
                                    x = modern_pollen$Tmin)
