library(ggplot2)
library(kernelshap)

# Turn ordinal factors into unordered
ord <- c("clarity", "color", "cut")
diamonds[, ord] <- lapply(diamonds[ord], factor, ordered = FALSE)

# Fit model
fit <- lm(log(price) ~ log(carat) * (clarity + color + cut), data = diamonds)

# Subset of 270 diamonds
X_small <- diamonds[seq(1, nrow(diamonds), 200), c("carat", ord)]

# Exact KernelSHAP on X_small, using X_small as background data (3 seconds)
system.time(
  ks <- kernelshap(fit, X_small, bg_X = X_small)  
)
ks

# SHAP values of first 2 observations:
#             carat    clarity     color         cut
# [1,] -2.09368837 -0.2875728 0.1165124  0.01496767
# [2,]  0.01148493 -0.1191795 0.1115798 -0.02016471

# Sampling version takes a bit longer (10 seconds)
system.time(
  ks2 <- kernelshap(fit, X_small, bg_X = X_small, exact = FALSE)  
)
ks2

# SHAP values of first 2 observations:
#             carat    clarity     color         cut
# [1,] -2.09368837 -0.2875728 0.1165124  0.01496767
# [2,]  0.01148493 -0.1191795 0.1115798 -0.02016471

library(shapviz)
sv <- shapviz(ks)
sv_dependence(sv, "carat", "auto")

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
