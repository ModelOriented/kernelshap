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

# Exact KernelSHAP (4s)
system.time(
  ks <- kernelshap(fit, X_small, bg_X = bg_X)  
)
ks

# SHAP values of first 2 observations:
#          carat     clarity     color        cut
# [1,] -2.050074 -0.28048747 0.1281222 0.01587382
# [2,] -2.085838  0.04050415 0.1283010 0.03731644

# Sampling version takes a bit longer (6 seconds)
system.time(
  ks2 <- kernelshap(fit, X_small, bg_X = bg_X, sampling_strategy = "paired")  
)
ks2

# SHAP values of first 2 observations:
#          carat     clarity     color        cut
# [1,] -2.050074 -0.28048747 0.1281222 0.01587382
# [2,] -2.085838  0.04050415 0.1283010 0.03731644

# Using parallel backend (1 second from the second call on)
library("doFuture")

registerDoFuture()
plan(multisession, workers = 2)  # Windows
# plan(multicore, workers = 3)   # Linux, macOS, Solaris

# 3 seconds
system.time(
  ks3 <- kernelshap(fit, X_small, bg_X = bg_X, parallel = TRUE)  
)
ks3

library(shapviz)
sv <- shapviz(ks)
sv_dependence(sv, "carat", "auto")


# More features (but non-sensical model)
# Fit model
fit <- lm(
  log(price) ~ log(carat) * (clarity + color + cut) + x + y + z + table + depth, 
  data = diamonds
)

# Subset of 1018 diamonds to explain
X_small <- diamonds[seq(1, nrow(diamonds), 53), setdiff(names(diamonds), "price")]

# Exact KernelSHAP on X_small, using X_small as background data 
# (56 seconds for exact, 30s for hybrid, 24s for hybrid parallel)
system.time(
  ks <- kernelshap(fit, X_small, bg_X = bg_X, parallel=T)  
)
ks

# SHAP values of first 2 observations:
#   carat        cut     color     clarity         depth         table          x           y            z
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
x = ["carat"] + ord
diamonds[ord] = diamonds[ord].apply(lambda x: x.cat.codes)
X = diamonds[x].to_numpy()

# Fit model with interactions and dummy variables
fit = ols(
  "np.log(price) ~ np.log(carat) * (C(clarity) + C(cut) + C(color))", 
  data=diamonds
).fit()

# Define subset of 270 diamonds
X_small = X[0:len(X):200]

# Calculate KernelSHAP values for X_small, using X_small as 
# background data. Predict function has to map numpy array 
# to pd.Dataframe
ks = KernelExplainer(
  model=lambda X: fit.predict(pd.DataFrame(X, columns=x)), 
  data = X_small
)
sv = ks.shap_values(X_small)  # 30 seconds
sv[0:2]

# array([[-2.09368837, -0.28757279,  0.11651242,  0.01496767],
#        [ 0.01148493, -0.1191795 ,  0.11157978, -0.02016471]])
