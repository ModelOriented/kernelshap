library(ggplot2)
library(kernelshap)

fit <- lm(log(price) ~ log(carat) + I(as.integer(cut)-1) + I(as.integer(color)-1) + I(as.integer(clarity)-1), data = diamonds)
x <- c("carat", "cut", "color", "clarity")
ks <- kernelshap(head(diamonds[x]), function(z) predict(fit, z), bg_X = head(diamonds[x], 100))
ks

# SHAP values of first 2 observations:
#           carat        cut     color    clarity
# [1,] -0.4992360 0.05224771 0.1512496 -0.2424388
# [2,] -0.6700215 0.02038935 0.1512496 -0.1193734

#========================
# The same in Python
#========================

# import numpy as np
# import pandas as pd
# from sklearn.preprocessing import OrdinalEncoder, FunctionTransformer
# from sklearn.compose import ColumnTransformer
# from sklearn.linear_model import LinearRegression
# from plotnine.data import diamonds
# 
# x_ord = ["cut", "color", "clarity"]
# ord_levels = [diamonds[x].cat.categories.to_list() for x in x_ord]
# ord_levels
# 
# preprocessor = ColumnTransformer(
#   transformers=[
#     ("log", FunctionTransformer(np.log), ["carat"]),
#     ("encoder", OrdinalEncoder(categories=ord_levels), x_ord)
#   ],
#   remainder="drop"
# )
# 
# X = preprocessor.fit_transform(diamonds)
# y = np.log(diamonds.price)
# 
# ols = LinearRegression()
# ols.fit(X, y)
# 
# ols.coef_, ols.intercept_
# # (array([ 1.87734554, 0.03185836, -0.0779637 , 0.12306539]), 8.262516487229426)
# 
# from shap import KernelExplainer
# 
# ks = KernelExplainer(model=lambda X: ols.predict(X), data = X[0:100])
# sv = ks.shap_values(X[0:2])
# sv
# 
# # array([[-0.49923603,  0.05224771,  0.15124958, -0.24243881],
# #       [-0.6700215 ,  0.02038935,  0.15124958, -0.11937342]])
