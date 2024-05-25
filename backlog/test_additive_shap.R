library(mgcv)
library(gam)
library(survival)
library(splines)
library(testthat)

formulas_ok <- list(
  Sepal.Length ~ Sepal.Width + Petal.Width + Species,
  Sepal.Length ~ log(Sepal.Width) + poly(Petal.Width, 2) + ns(Petal.Length, 2),
  Sepal.Length ~ log(Sepal.Width) + poly(Sepal.Width, 2)
)

formulas_bad <- list(
  Sepal.Length ~ Species * Petal.Length,
  Sepal.Length ~ Species + Petal.Length + Species:Petal.Length,
  Sepal.Length ~ log(Petal.Length / Petal.Width)
)

models <- list(mgcv::gam, mgcv::bam, gam::gam)

for (formula in formulas_ok) {
  for (model in models) {
    fit <- model(formula, data = iris)
    s <- additive_shap(fit, head(iris), verbose = FALSE)
  }
}

for (formula in formulas_bad) {
  for (model in models) {
    fit <- model(formula, data = iris)
    expect_error(s <- additive_shap(fit, head(iris), verbose = FALSE))
  }
}

# Survival
iris$s <- rep(1, nrow(iris))
formulas_ok <- list(
  Surv(Sepal.Length, s) ~ Sepal.Width + Petal.Width + Species,
  Surv(Sepal.Length, s) ~ log(Sepal.Width) + poly(Petal.Width, 2) + ns(Petal.Length, 2),
  Surv(Sepal.Length, s) ~ log(Sepal.Width) + poly(Sepal.Width, 2)
)

formulas_bad <- list(
  Surv(Sepal.Length, s) ~ Species * Petal.Length,
  Surv(Sepal.Length, s) ~ Species + Petal.Length + Species:Petal.Length,
  Surv(Sepal.Length, s) ~ log(Petal.Length / Petal.Width)
)

models <- list(survival::coxph, survival::survreg)

for (formula in formulas_ok) {
  for (model in models) {
    fit <- model(formula, data = iris)
    s <- additive_shap(fit, head(iris), verbose = FALSE)
  }
}

for (formula in formulas_bad) {
  for (model in models) {
    fit <- model(formula, data = iris)
    expect_error(s <- additive_shap(fit, head(iris), verbose = FALSE))
  }
}
