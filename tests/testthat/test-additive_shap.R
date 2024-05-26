test_that("simple additive formula gives same as permshap()", {
  form <- Sepal.Length ~ .
  fit_lm <- lm(form, data = iris)
  fit_glm <- glm(form, data = iris, family = quasipoisson)
  
  s_add_lm <- additive_shap(fit_lm, head(iris), verbose = FALSE)
  s_add_glm <- additive_shap(fit_glm, head(iris), verbose = FALSE)
  
  X <- head(iris[-1L])
  s_perm_lm <- permshap(fit_lm, X = X, bg_X = iris, verbose = FALSE)
  s_perm_glm <- permshap(
    fit_glm, X = X, bg_X = iris, verbose = FALSE
  )
  expect_equal(s_add_lm$S, s_perm_lm$S)
  expect_equal(s_add_glm$S, s_perm_glm$S)
  expect_equal(s_add_lm$predictions, unname(predict(fit_lm, newdata = X)))
  expect_equal(s_add_glm$predictions, unname(predict(fit_glm, newdata = X)))
})

test_that("formula where feature appears in two terms gives same as permshap()", {
  form <- Sepal.Length ~ log(Sepal.Width) + poly(Sepal.Width, 2) + Petal.Length
  fit_lm <- lm(form, data = iris)
  fit_glm <- glm(form, data = iris, family = quasipoisson)
  
  s_add_lm <- additive_shap(fit_lm, head(iris), verbose = FALSE)
  s_add_glm <- additive_shap(fit_glm, head(iris), verbose = FALSE)
  
  X <- head(iris[2:3])
  s_perm_lm <- permshap(fit_lm, X = X, bg_X = iris, verbose = FALSE)
  s_perm_glm <- permshap(
    fit_glm, X = X, bg_X = iris, verbose = FALSE
  )
  expect_equal(s_add_lm$S, s_perm_lm$S)
  expect_equal(s_add_glm$S, s_perm_glm$S)
  expect_equal(s_add_lm$predictions, unname(predict(fit_lm, newdata = X)))
  expect_equal(s_add_glm$predictions, unname(predict(fit_glm, newdata = X)))
})

test_that("formula with complicated terms gives same as permshap()", {
  form <- Sepal.Length ~ 
    log(Sepal.Width) + Species + poly(Petal.Length, 2)
  
  fit_lm <- lm(form, data = iris)
  fit_glm <- glm(form, data = iris, family = quasipoisson)
  
  s_add_lm <- additive_shap(fit_lm, head(iris), verbose = FALSE)
  s_add_glm <- additive_shap(fit_glm, head(iris), verbose = FALSE)
  
  X <- head(iris[c(2, 3, 5)])
  s_perm_lm <- permshap(fit_lm, X = X, bg_X = iris, verbose = FALSE)
  s_perm_glm <- permshap(
    fit_glm, X = X, bg_X = iris, verbose = FALSE
  )
  expect_equal(s_add_lm$S, s_perm_lm$S)
  expect_equal(s_add_glm$S, s_perm_glm$S)
  expect_equal(s_add_lm$predictions, unname(predict(fit_lm, newdata = X)))
  expect_equal(s_add_glm$predictions, unname(predict(fit_glm, newdata = X)))
})

test_that("formulas with more than one covariate per term fail", {
  formulas_bad <- list(
    Sepal.Length ~ Species * Petal.Length,
    Sepal.Length ~ Species + Petal.Length + Species:Petal.Length,
    Sepal.Length ~ log(Petal.Length / Petal.Width)
  )
  
  for (formula in formulas_bad) {
    fit <- lm(formula, data = iris)
    expect_error(s <- additive_shap(fit, head(iris), verbose = FALSE))
    
    fit <- glm(formula, data = iris, family = quasipoisson)
    expect_error(s <- additive_shap(fit, head(iris), verbose = FALSE))
  }  
})
