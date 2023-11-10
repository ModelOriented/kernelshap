test_that("rowpaste() does what it should", {
  M <- cbind(c(0, 0), c(1, 0), c(1, 1))
  expect_equal(rowpaste(M), c("011", "001"))
})

test_that("shapley_weights() does what it should", {
  expect_equal(shapley_weights(5, 2), factorial(2) * factorial(5 - 2 - 1) / factorial(5))
})
