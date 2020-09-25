test_that("hexagonal logo works", {
  hex_logo(output = "hex_logo.png")
  expect_true(file.exists("hex_logo.png"))
  expect_false(dir.exists("hex_logo.png"))
  expect_gt(file.size("hex_logo.png"), 0)
  file.remove("hex_logo.png")
  expect_false(file.exists("hex_logo.png"))
})

test_that("parallel benchmark works", {
  # Define toy function that sleeps for (60/cpus) seconds
  a <- function(cpus) {Sys.sleep(4/cpus)}
  times_df <- par_benchmark(c(1, 2, 4), a, quiet = TRUE)
  expect_equal(floor(times_df$times), c(4, 2, 1))
  expect_output(par_benchmark(c(4), a, quiet = FALSE))
  expect_output(par_benchmark(c(4), a, quiet = FALSE, plot = TRUE))
  expect_true(file.exists("Rplots.pdf"))
  expect_false(dir.exists("Rplots.pdf"))
  # expect_gt(file.size("Rplots.pdf"), 0)
  file.remove("Rplots.pdf")
  expect_false(file.exists("Rplots.pdf"))
})
