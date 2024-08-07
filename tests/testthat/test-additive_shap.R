test_that("Additive formulas give same as agnostic SHAP with full training data as bg data", {
  formulas <- list(
    Sepal.Length ~ .,
    Sepal.Length ~ log(Sepal.Width) + poly(Sepal.Width, 2) + Petal.Length,
    form <- Sepal.Length ~ log(Sepal.Width) + Species + poly(Petal.Length, 2)
  )
  xvars <- list(
    setdiff(colnames(iris), "Sepal.Length"),
    c("Sepal.Width", "Petal.Length"),
    xvars <- c("Sepal.Width", "Petal.Length", "Species")
  )
  
  for (j in seq_along(formulas)) {
    fit <- list(
      lm = lm(formulas[[j]], data = iris),
      glm = glm(formulas[[j]], data = iris, family = quasipoisson)
    )
    
    shap1 <- lapply(fit, additive_shap, head(iris), verbose = FALSE)
    shap2 <- lapply(
      fit, permshap, head(iris), bg_X = iris, verbose = FALSE, feature_names = xvars[[j]]
    )
    shap3 <- lapply(
      fit, kernelshap, head(iris), bg_X = iris, verbose = FALSE, feature_names = xvars[[j]]
    )
    
    for (i in seq_along(fit)) {
      expect_equal(shap1[[i]]$S, shap2[[i]]$S)
      expect_equal(shap1[[i]]$S, shap3[[i]]$S)
    }
  }
})

test_that("formulas with more than one covariate per term fail", {
  formulas_bad <- list(
    Sepal.Length ~ Species * Petal.Length,
    Sepal.Length ~ Species + Petal.Length + Species:Petal.Length,
    Sepal.Length ~ log(Petal.Length / Petal.Width)
  )
  
  for (formula in formulas_bad) {
    fit <- list(
      lm = lm(formula, data = iris),
      glm = glm(formula, data = iris, family = quasipoisson)
    )
    for (f in fit)
      expect_error(additive_shap(f, head(iris), verbose = FALSE))
  }  
})

