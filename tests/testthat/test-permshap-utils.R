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

test_that("sample_Z_from_chain() provides correct matrices", {
  J <- c(1, 3, 2)
  nms <- c("a", "b", "c")

  result <- sample_Z_from_chain(J, nms)

  # Switch off in the order of J
  expected_forward <- rbind(
    c(FALSE, TRUE, TRUE),
    c(FALSE, TRUE, FALSE)
  )

  # Swith on in the order of J
  expected_backward <- rbind(
    c(TRUE, FALSE, FALSE),
    c(TRUE, FALSE, TRUE)
  )
  expected <- rbind(expected_forward, expected_backward)
  colnames(expected) <- nms
  expect_equal(result, expected)

  # Large p
  p <- 10L
  chain <- 1L:p
  result <- sample_Z_from_chain(chain, letters[chain])
  expect_equal(dim(result), c(2L * (p - 1L), p))
  expect_equal(result[1L:(p - 1L), ], !result[p:(2L * (p - 1L)), ])
})

test_that("pad_vz() works", {
  # One-dimensional output
  v0 <- matrix(0, dimnames = list("a", "x"))
  v1 <- matrix(1, dimnames = list(NULL, "x"))
  vz <- cbind(c(2, 3, 3, 2))
  expect_equal(pad_vz(vz, v0, v1), expected = cbind(x = c(1, 2, 3, 0, 3, 2, 1)))

  # Two-dimensional output
  v0 <- rbind(c(0, 0))
  dimnames(v0) <- list("a", c("x", "y"))
  v1 <- rbind(c(1, 1))
  colnames(v1) <- c("x", "y")
  vz <- cbind(x = c(2, 3, 3, 2), y = c(3, 4, 4, 3))
  expected <- cbind(
    x = c(1, 2, 3, 0, 3, 2, 1),
    y = c(1, 3, 4, 0, 4, 3, 1)
  )
  expect_equal(pad_vz(vz, v0, v1), expected = expected)
})
