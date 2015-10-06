library(BioCompToolsR)
context("sampler statistics")

test_that("Fst OF testing", {

  expect_equal(FstOF( t(matrix(c(0.5,0.5))), c(0.5,0.5) ), 0.0)

})
