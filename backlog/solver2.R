# https://github.com/iancovert/shapley-regression/blob/master/shapreg/shapley.py
solver2 <- function(A, b, constraint) {
  A_inv_1 <- solve(A, matrix(1, nrow = nrow(A)))
  A_inv_b <- solve(A, b)
  num <- colSums(A_inv_b) - constraint
  A_inv_b - A_inv_1 %*% num / sum(A_inv_1)
}
