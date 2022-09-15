test_that("Sum of kernel weights is 1", {
  expect_equal(sum(kernel_weights(5)), 1)
  expect_equal(sum(kernel_weights(10)), 1)
})

test_that("Sum of kernel weights is 1, even for subset of domain", {
  expect_equal(sum(kernel_weights(10, S = 2:5)), 1)
})

p <- 10
m <- 100

test_that("Random z have right output dim and the sums are between 1 and p-1", {
  Z <- sample_Z(p, m = m)
  
  expect_equal(dim(Z), c(m, p))
  expect_true(all(rowSums(Z) %in% 1:(p-1)))
})

test_that("Random z have right output dim and the sums are in subset S", {
  S <- 2:3
  Z <- sample_Z(p, m = m, S = S)
  
  expect_equal(dim(Z), c(m, p))
  expect_true(all(rowSums(Z) %in% S))
})

test_that("Sampling input looks ok (deg = 0)", {
  input <- input_sampling(p, m = m, deg = 0, paired = TRUE)

  expect_equal(dim(input$Z), c(m, p))
  expect_equal(sum(input$w), 1)
  expect_equal(dim(input$A), c(p, p))
  expect_equal(diag(input$A), rep(0.5, p))
})

test_that("Sampling input looks ok (deg = 0, unpaired)", {
  input <- input_sampling(p, m = m, deg = 0, paired = FALSE)
  
  expect_equal(dim(input$Z), c(m, p))
  expect_equal(sum(input$w), 1)
  expect_equal(dim(input$A), c(p, p))
 # expect_equal(diag(input$A), rep(0.5, p)) # This is not TRUE
})

test_that("Sampling input looks ok (deg = 1)", {
  input <- input_sampling(p, m = m, deg = 1, paired = TRUE)
  
  expect_equal(dim(input$Z), c(m, p))
  expect_true(sum(input$w) < 1)
  expect_equal(dim(input$A), c(p, p))
  expect_true(all(diag(input$A) < 0.5))
})

test_that("Sampling input looks ok (deg = 2)", {
  input <- input_sampling(p, m = m, deg = 2, paired = TRUE)
  
  expect_equal(dim(input$Z), c(m, p))
  expect_true(sum(input$w) < 1)
  expect_equal(dim(input$A), c(p, p))
  expect_true(all(diag(input$A) < 0.5))
})

