library(ranger)
library(survival)
library(kernelshap)

set.seed(1)

fit <- ranger(Surv(time, status) ~ ., data = veteran, num.trees = 20)
fit2 <- ranger(time ~ . - status, data = veteran, num.trees = 20)
fit3 <- ranger(time ~ . - status, data = veteran, quantreg = TRUE, num.trees = 20)
fit4 <- ranger(status ~ . - time, data = veteran, probability = TRUE, num.trees = 20)

xvars <- setdiff(colnames(veteran), c("time", "status"))

kernelshap(fit, head(veteran), feature_names = xvars, bg_X = veteran)
kernelshap(fit, head(veteran), feature_names = xvars, bg_X = veteran, survival = "prob")
kernelshap(fit2, head(veteran), feature_names = xvars, bg_X = veteran)
kernelshap(fit3, head(veteran), feature_names = xvars, bg_X = veteran, type = "quantiles")
kernelshap(fit4, head(veteran), feature_names = xvars, bg_X = veteran)

permshap(fit, head(veteran), feature_names = xvars, bg_X = veteran)
permshap(fit, head(veteran), feature_names = xvars, bg_X = veteran, survival = "prob")
permshap(fit2, head(veteran), feature_names = xvars, bg_X = veteran)
permshap(fit3, head(veteran), feature_names = xvars, bg_X = veteran, type = "quantiles")
permshap(fit4, head(veteran), feature_names = xvars, bg_X = veteran)

# Sampling
permshap(fit, head(veteran), feature_names = xvars, bg_X = veteran, exact = FALSE)
permshap(fit, head(veteran), feature_names = xvars, bg_X = veteran, exact = FALSE, low_memory = TRUE)
permshap(fit, head(veteran), feature_names = xvars, bg_X = veteran, survival = "prob", exact = FALSE)
permshap(fit2, head(veteran), feature_names = xvars, bg_X = veteran, exact = FALSE)
permshap(fit3, head(veteran), feature_names = xvars, bg_X = veteran, type = "quantiles", exact = FALSE)
permshap(fit4, head(veteran), feature_names = xvars, bg_X = veteran, exact = FALSE)
permshap(fit4, head(veteran), feature_names = xvars, bg_X = veteran, exact = FALSE, low_memory = TRUE)
