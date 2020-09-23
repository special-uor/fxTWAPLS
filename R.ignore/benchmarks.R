##################################### cv.w #####################################
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::cv.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = reconbio::WAPLS.w, 
                                    predictfun = reconbio::WAPLS.predict.w)

################################## get_pseduo ##################################
tictoc::tic("Tmin")
pseduo_Tmin <- reconbio::get_pseduo(dist, modern_pollen$Tmin, cpus = 15)
tictoc::toc()

out <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12), 
                               FUN = reconbio::get_pseduo, 
                               plot = TRUE, 
                               dist = dist, 
                               x = modern_pollen$Tmin)
