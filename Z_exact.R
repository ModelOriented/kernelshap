# This function creates z vectors and replicates them to match Kernel SHAP weights
create_exact_Z <- function(p) {
  stopifnot(p >= 2L, p <= 5)

  # Create all z vectors
  Zi <- do.call(expand.grid, replicate(p, 0:1, simplify = FALSE))
  Zi <- Zi[2:(nrow(Zi) - 1), ]
  rs <- unname(rowSums(Zi))
  
  # Replicate in a way to end up with correct distribution of sums
  S <- 1:(p - 1)
  probs <- (p - 1) / (choose(p, S) * S * (p - S))
  probs <- probs / sum(probs)
  fct <- probs / as.vector(table(rs))
  fct <- fct / min(fct)  # Does only work for p <= 6
  reps <- round(fct[rs])
  
  Zi <- Zi[rep(1:nrow(Zi), times = reps), , drop = FALSE]
  Zi <- as.matrix(Zi)
  dimnames(Zi) <- NULL
  
  check <- as.vector(prop.table(table(rowSums(Zi))))
  stopifnot(all.equal(check, probs))
  Zi
}

Z_exact <- c(list(matrix(NA)), lapply(2:5, create_exact_Z))
