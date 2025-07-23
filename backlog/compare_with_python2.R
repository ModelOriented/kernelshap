library(kernelshap)

n <- 100

X <- data.frame(
  x1 = seq(1:n) / 100,
  x2 = log(1:n),
  x3 = sqrt(1:n),
  x4 = sin(1:n),
  x5 = (seq(1:n) / 100)^2,
  x6 = cos(1:n)
)
head(X)

pf <- function(model, newdata) {
  x <- newdata
  x[, 1] * x[, 2] * x[, 3] * x[, 4] + x[, 5] + x[, 6]
}
ks <- kernelshap(pf, head(X), bg_X = X, pred_fun = pf)
ks # -1.196216 -1.241848 -0.9567848 3.879420 -0.33825  0.5456252
es <- permshap(pf, head(X), bg_X = X, pred_fun = pf)
es # -1.196216 -1.241848 -0.9567848 3.879420 -0.33825  0.5456252

set.seed(10)
kss <- kernelshap(
  pf,
  head(X, 1),
  bg_X = X,
  pred_fun = pf,
  hybrid_degree = 0,
  exact = F,
  m = 9000,
  max_iter = 100,
  tol = 0.0005
)
kss # -1.203484 -1.241578 -0.9653531 3.889111 -0.3360786 0.5493293

set.seed(2)
ksh <- kernelshap(
  pf,
  head(X, 1),
  bg_X = X,
  pred_fun = pf,
  hybrid_degree = 1,
  exact = FALSE,
  max_iter = 10000,
  tol = 0.0005
)
ksh # -1.195273 -1.240459 -0.9546994 3.884469 -0.3382316 0.5361404

set.seed(1)
ksh2 <- kernelshap(
  pf,
  head(X, 1),
  bg_X = X,
  pred_fun = pf,
  hybrid_degree = 2,
  exact = FALSE,
  m = 10000,
  max_iter = 10000,
  tol = 0.0001
)
ksh2 # -1.195725 -1.241576 -0.9567066 3.87909 -0.3385894 0.5454539

set.seed(1)
ps0 <- permshap(
  pf,
  head(X, 1),
  bg_X = X,
  pred_fun = pf,
  exact = FALSE,
  low_memory = FALSE,
  max_iter = 20000,
  tol = 0.0001
)
ps0 # -1.20929 -1.225766 -0.9602608 3.879889 -0.33825 0.5456252
