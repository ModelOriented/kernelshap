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
ps <- permshap(pf, head(X), bg_X = X, pred_fun = pf)
ps # -1.196216 -1.241848 -0.9567848 3.879420 -0.33825  0.5456252

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
kss # -1.198078 -1.246508 -0.9580638 3.877532 -0.3241824 0.541247

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
ksh # -1.191981 -1.240656 -0.9516264 3.86776 -0.3342143 0.5426642

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
ksh2 # 1.195976 -1.241107 -0.9565121 3.878891 -0.3384621 0.5451118

set.seed(1)
pss <- permshap(
  pf,
  head(X, 1),
  bg_X = X,
  pred_fun = pf,
  exact = FALSE,
  max_iter = 40000,
  tol = 0.0001
)
pss # -1.222608 -1.252001 -0.9312635 3.890444 -0.33825 0.5456252 non-convergence
