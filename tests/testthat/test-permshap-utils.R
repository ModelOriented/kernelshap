test_that("rowpaste() does what it should", {
  M <- cbind(c(0, 0), c(1, 0), c(1, 1))
  expect_equal(rowpaste(M), c("011", "001"))
})

test_that("shapley_weights() does what it should", {
  expect_equal(shapley_weights(5, 2), factorial(2) * factorial(5 - 2 - 1) / factorial(5))
})

test_that("balanced_chains() returns output of correct structure", {
  p <- 5L
  result <- balanced_chains(p)
  expect_equal(length(result), p)
  expect_equal(sapply(result, `[[`, 1L), 1:p)
  expect_equal(sapply(result, sort), .row(c(p, p)))
})

test_that("exact_Z_balanced() returns correct matrix", {
  expected <- rbind(
    c(FALSE, TRUE, TRUE),
    c(TRUE, FALSE, TRUE),
    c(TRUE, TRUE, FALSE),
    c(TRUE, FALSE, FALSE),
    c(FALSE, TRUE, FALSE),
    c(FALSE, FALSE, TRUE)
  )
  colnames(expected) <- 1:3
  expect_equal(exact_Z_balanced(3, 1:3), expected)
})

test_that("sample_Z_from_chain() provides correct matrices", {
  J <- c(1, 3, 2, 4, 5)
  nms <- letters[seq_along(J)]

  result <- sample_Z_from_chain(J, nms)

  # Switch off in the order of J
  expected_forward <- rbind(
    c(FALSE, TRUE, FALSE, TRUE, TRUE),
    c(FALSE, FALSE, FALSE, TRUE, TRUE)
  )
  expected_backward <- !expected_forward
  expected <- rbind(expected_forward, expected_backward)
  colnames(expected) <- nms
  expect_equal(result, expected)

  # Large p
  p <- 10L
  chain <- 1L:p
  result <- sample_Z_from_chain(chain, letters[chain])
  expect_equal(dim(result), c(2L * (p - 3L), p))
  expect_equal(result[1L:(p - 3L), ], !result[(p - 2L):(2L * (p - 3L)), ])
})

test_that("init_vzj() provides correct matrices", {
  p <- 4L
  v0 <- cbind(1, 2)
  v1 <- cbind(a = 2, b = 5)
  expected_1 <- c(2, 0, 0, 0, 1, 0, 0, 0, 2)
  expected_2 <- c(5, 0, 0, 0, 2, 0, 0, 0, 5)
  expect_equal(init_vzj(p, v0 = v0, v1 = v1), cbind(a = expected_1, b = expected_2))

  # One-dim
  expect_equal(
    init_vzj(p, v0 = v0[, 1L, drop = FALSE], v1 = v1[, 1L, drop = FALSE]),
    cbind(a = expected_1)
  )
})
