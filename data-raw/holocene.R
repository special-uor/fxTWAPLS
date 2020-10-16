## code to prepare `holocene` dataset
holocene <- read.csv(system.file("extdata", 
                                 "Holocene.csv",
                                 package = "fxTWAPLS",
                                 mustWork = TRUE),
                     row.names = 1)
# Removed accented characters and special symbols from variable names
cln_names <- sapply(names(holocene), fxTWAPLS:::cln_str, keep = "\\._-")
colnames(holocene) <- cln_names
usethis::use_data(holocene, overwrite = TRUE)
