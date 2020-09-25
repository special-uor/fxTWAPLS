##################################### cv.w #####################################
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::cv.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = reconbio::WAPLS.w, 
                                    predictfun = reconbio::WAPLS.predict.w)

################################### cv.pr.w ####################################
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::cv.pr.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = reconbio::WAPLS.w, 
                                    predictfun = reconbio::WAPLS.predict.w,
                                    pseduo = pseduo_Tmin)

################################## get_pseduo ##################################
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::get_pseduo, 
                                    plot = TRUE, 
                                    dist = dist, 
                                    x = modern_pollen$Tmin)
