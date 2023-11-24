library(kernelshap)
library(ranger)

differences <- numeric(4)

set.seed(1)

for (depth in 1:4) {
  fit <- ranger(
    Sepal.Length ~ ., 
    mtry = 3,
    data = iris, 
    max.depth = depth
  )
  ps <- permshap(fit, iris[2:5], bg_X = iris)
  ks <- kernelshap(fit, iris[2:5], bg_X = iris)
  differences[depth] <- mean(abs(ks$S - ps$S))
}

differences  # for tree depth 1, 2, 3, 4
# 5.053249e-17 9.046443e-17 2.387905e-04 4.403375e-04

# SHAP values of first two rows with tree depth 4
ps
#      Sepal.Width Petal.Length Petal.Width      Species
# [1,]  0.11377616   -0.7130647  -0.1956012 -0.004437022
# [2,] -0.06852539   -0.7596562  -0.2259017 -0.006575266

ks
#      Sepal.Width Petal.Length Petal.Width      Species
# [1,]  0.11463191   -0.7125194  -0.1951810 -0.006258208
# [2,] -0.06828866   -0.7597391  -0.2259833 -0.006647530

  
# larger data, more features
library(xgboost)
library(shapviz)

colnames(miami) <- tolower(colnames(miami))
miami$log_ocean <- log(miami$ocean_dist)
x <- c("log_ocean", "tot_lvg_area", "lnd_sqfoot", "structure_quality", "age", "month_sold")

# Train/valid split
set.seed(1)
ix <- sample(nrow(miami), 0.8 * nrow(miami))

y_train <- log(miami$sale_prc[ix])
y_valid <- log(miami$sale_prc[-ix])
X_train <- data.matrix(miami[ix, x])
X_valid <- data.matrix(miami[-ix, x])

dtrain <- xgb.DMatrix(X_train, label = y_train)
dvalid <- xgb.DMatrix(X_valid, label = y_valid)

# Fit via early stopping (depth 1 to 3)
differences <- numeric(3)

for (i in 1:3) {
  fit <- xgb.train(
    params = list(learning_rate = 0.15, objective = "reg:squarederror", max_depth = i),
    data = dtrain,
    watchlist = list(valid = dvalid),
    early_stopping_rounds = 20,
    nrounds = 1000,
    callbacks = list(cb.print.evaluation(period = 100))
  )
  ps <- permshap(fit, X = head(X_valid, 500), bg_X = head(X_valid, 500))
  ks <- kernelshap(fit, X = head(X_valid, 500), bg_X = head(X_valid, 500))
  differences[i] <- mean(abs(ks$S - ps$S))
}
differences
# 2.904010e-09 5.158383e-09 6.586577e-04

ps
# SHAP values of first observations:
# log_ocean tot_lvg_area lnd_sqfoot structure_quality        age  month_sold
# 0.2224359   0.04941044  0.1266136         0.1360166 0.01036866 0.005557032
# 0.3674484   0.01045079  0.1192187         0.1180312 0.01426247 0.005465283

ks
# SHAP values of first observations:
# log_ocean tot_lvg_area lnd_sqfoot structure_quality        age  month_sold
# 0.2245202  0.049520308  0.1266020         0.1349770 0.01142703 0.003355770
# 0.3697167  0.009575195  0.1198201         0.1168738 0.01544061 0.003450425
