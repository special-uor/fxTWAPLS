# Load modern pollen data
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))

# Extract taxa
taxaColMin <- which(colnames(modern_pollen) == "Abies")
taxaColMax <- which(colnames(modern_pollen) == "Zygophyllaceae")
taxa <- modern_pollen[, taxaColMin:taxaColMax]

# Get the frequency of each climate variable fx
fx_Tmin <- fxTWAPLS::fx(modern_pollen$Tmin, bin = 0.02)
fx_gdd <- fxTWAPLS::fx(modern_pollen$gdd, bin = 20)
fx_alpha <- fxTWAPLS::fx(modern_pollen$alpha, bin = 0.002)

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

test_that("fx works", {
  expect_equal(length(fx_Tmin), 6458)
  expect_equal(length(fx_gdd), 6458)
  expect_equal(length(fx_alpha), 6458)
})

test_that("WPALS training function works", {
  expected_names <- c("fit", "x", "taxon_name", "optimum", "comp", "u", 
                      "z", "s", "orth", "alpha", "meanx", "nPLS")
  expect_equal(names(fit_Tmin), expected_names)
  expect_equal(names(fit_f_Tmin), expected_names)
})

test_that("TWPALS training function works", {
  expected_names <- c("fit", "x", "taxon_name", "optimum", "comp", "u", 
                      "t", "z", "s", "orth", "alpha", "meanx", "nPLS")
  expect_equal(names(fit_t_Tmin), expected_names)
  expect_equal(names(fit_tf_Tmin), expected_names)
})

test_that("LOOCV as in rioja works", {
  # MTCO
  # cv_Tmin <- fxTWAPLS::cv.w(taxa[1:500, ],
  #                           modern_pollen$Tmin[1:500],
  #                           nPLS = 5,
  #                           fxTWAPLS::WAPLS.w,
  #                           fxTWAPLS::WAPLS.predict.w,
  #                           cpus = 4)
})
