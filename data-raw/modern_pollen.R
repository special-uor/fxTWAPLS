## code to prepare `modern_pollen` dataset
modern_pollen <- read.csv(system.file("extdata", 
                                      "Modern_Pollen_gdd_alpha_Tmin.csv", 
                                      package = "fxTWAPLS", 
                                      mustWork = TRUE))
# Removed accented characters and special symbols 
## From variable names
cln_names <- sapply(names(modern_pollen), fxTWAPLS:::cln_str, keep = "\\._-")
colnames(modern_pollen) <- cln_names
# From entities
cln_entities <- sapply(modern_pollen$Entity.name, 
                       fxTWAPLS:::cln_str, keep = "\\.\\(\\)_-")
modern_pollen$Entity.name <- cln_entities
usethis::use_data(modern_pollen, overwrite = TRUE)
