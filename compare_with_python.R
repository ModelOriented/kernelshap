library(ggplot2)
library(kernelshap)

# Turn ordinal factors into unordered
ord <- c("clarity", "color", "cut")
diamonds[, ord] <- lapply(diamonds[ord], factor, ordered = FALSE)

# Fit model
fit <- lm(log(price) ~ log(carat) * (clarity + color + cut), data = diamonds)

# Subset of 120 diamonds used as background data
bg_X <- diamonds[seq(1, nrow(diamonds), 450), ]

# Subset of 1018 diamonds to explain
X_small <- diamonds[seq(1, nrow(diamonds), 53), c("carat", ord)]

# Exact KernelSHAP (5s)
system.time(
  ks <- kernelshap(fit, X_small, bg_X = bg_X)  
)
ks

# SHAP values of first 2 observations:
#          carat     clarity     color        cut
# [1,] -2.050074 -0.28048747 0.1281222 0.01587382
# [2,] -2.085838  0.04050415 0.1283010 0.03731644

# Pure sampling version takes a bit longer (13 seconds)
system.time(
  ks2 <- kernelshap(fit, X_small, bg_X = bg_X, exact = FALSE, hybrid_degree = 0)  
)
ks2

# SHAP values of first 2 observations:
#          carat     clarity     color        cut
# [1,] -2.050074 -0.28048747 0.1281222 0.01587382
# [2,] -2.085838  0.04050415 0.1283010 0.03731644

# Using parallel backend
library("doFuture")

registerDoFuture()
plan(multisession, workers = 2)  # Windows
# plan(multicore, workers = 2)   # Linux, macOS, Solaris

# 3 seconds in the second call
system.time(
  ks3 <- kernelshap(fit, X_small, bg_X = bg_X, parallel = TRUE)  
)
ks3

library(shapviz)

sv <- shapviz(ks)
sv_dependence(sv, "carat")


# More features (but non-sensical model)
# Fit model
fit <- lm(
  log(price) ~ log(carat) * (clarity + color + cut) + x + y + z + table + depth, 
  data = diamonds
)

# Subset of 1018 diamonds to explain
X_small <- diamonds[seq(1, nrow(diamonds), 53), setdiff(names(diamonds), "price")]

# Exact KernelSHAP on X_small, using X_small as background data 
# (58/67(?) seconds for exact, 25/18 for hybrid deg 2, 16/9 for hybrid deg 1, 
# 26/17 for pure sampling; second number with 2 parallel sessions on Windows)
system.time(
  ks <- kernelshap(fit, X_small, bg_X = bg_X)  
)
ks

# SHAP values of first 2 observations:
#          carat        cut     color     clarity         depth         table          x           y            z
# [1,] -1.842799 0.01424231 0.1266108 -0.27033874 -0.0007084443  0.0017787647 -0.1720782 0.001330275 -0.006445693
# [2,] -1.876709 0.03856957 0.1266546  0.03932912 -0.0004202636 -0.0004871776 -0.1739880 0.001397792 -0.006560624

#========================
# The same in Python
#========================

import numpy as np
import pandas as pd
from plotnine.data import diamonds
from statsmodels.formula.api import ols
from shap import KernelExplainer

# Turn categoricals into integers because, inconveniently, kernel SHAP
# requires numpy array as input
ord = ["clarity", "color", "cut"]
# x = ["carat"] + ord + ["table", "depth", "x", "y", "z"]
x = ["carat"] + ord
diamonds[ord] = diamonds[ord].apply(lambda x: x.cat.codes)
X = diamonds[x].to_numpy()

# Fit model with interactions and dummy variables
fit = ols(
  "np.log(price) ~ np.log(carat) * (C(clarity) + C(cut) + C(color))", # + x + y + z + table + depth", 
  data=diamonds
).fit()

# Background data (120 rows)
bg_X = X[0:len(X):450]

# Define subset of 1018 diamonds to explain
X_small = X[0:len(X):53]

# Calculate KernelSHAP values
ks = KernelExplainer(
  model=lambda X: fit.predict(pd.DataFrame(X, columns=x)), 
  data = bg_X
)
sv = ks.shap_values(X_small)  # 11 minutes
sv[0:2]

# Output for four features (exact calculations, 1 minute)
# array([[-2.05007406, -0.28048747,  0.12812216,  0.01587382],
#        [-2.0858379 ,  0.04050415,  0.12830103,  0.03731644]])


# Output for nine features (exact calculations, 13 minutes)
# array([[-1.84279897e+00, -2.70338744e-01,  1.26610769e-01,
#         1.42423108e-02,  1.77876470e-03, -7.08444295e-04,
#         -1.72078182e-01,  1.33027467e-03, -6.44569296e-03],
#        [-1.87670887e+00,  3.93291219e-02,  1.26654599e-01,
#         3.85695742e-02, -4.87177593e-04, -4.20263565e-04,
#         -1.73988040e-01,  1.39779179e-03, -6.56062359e-03]])