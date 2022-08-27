library(ggplot2)
library(kernelshap)

# Turn ordinal factors into unordered
ord <- c("clarity", "color", "cut")
diamonds[, ord] <- lapply(diamonds[ord], factor, ordered = FALSE)

# Fit model
fit <- lm(log(price) ~ log(carat) * clarity + cut + color, data = diamonds)

# Apply Kernel SHAP to first two diamonds using every 500th diamond as background data
x <- c("carat", ord)
bg_data <- diamonds[seq(1, nrow(diamonds), 500), ]
ks <- kernelshap(diamonds[1:2, x], function(z) predict(fit, z), bg_X = bg_data)
ks

# SHAP values of first 2 observations:
#          carat    clarity    color         cut
# [1,] -2.180243 -0.2542866 0.109037 0.024558363
# [2,] -2.364556 -0.1029186 0.109037 0.002575027

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
  "np.log(price) ~ np.log(carat) * C(clarity) + C(cut) + C(color)", 
  data=diamonds
).fit()

# Calculate SHAP values for first two obs using every 500th diamond
# as background data. Predict function has to map numpy array to pd.Dataframe
ks = KernelExplainer(
  model=lambda X: fit.predict(pd.DataFrame(X, columns=x)), 
  data = X[0:len(X):500]
)
sv = ks.shap_values(X[0:2])
sv

# array([[-2.18024332, -0.2542866 ,  0.10903704,  0.02455836],
#        [-2.3645562 , -0.10291862,  0.10903704,  0.00257503]])

