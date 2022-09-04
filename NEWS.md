# kernelshap 0.1.900 DEVEL

## Breaking change

The interface of `kernelshap()` has been revised. Instead of specifying a predition function, it suffices now to pass the fitted model object. The default `pred_fun` is now `stats::predict`, which works in most cases. Some other cases are catched via model class ("ranger" and mlr3 "Learner"s). The `pred_fun` can be overwritten by a function of the form `function(object, X, ...)`.

Example: Logistic regression with predictions on logit scale

```
kernelshap(fit, X, bg_X)
```

Example: Logistic regression with predictions on probability scale

```
kernelshap(fit, X, bg_X, type = "response")
```

Example: Log-linear regression to be evaluated on original scale.
Here, the default predict function needs to be overwritten:

```
kernelshap(fit, X, bg_X, pred_fun = function(m, X) exp(predict(m, X)))
```

## Major improvements

- The package now supports multidimensional predictions. Hurray!
- Major speed-up.
- For $2 \le p \le 5$ features, the algorithm now returns exact Kernel SHAP values. (For $p = 1$, *exact Shapley values* are returned.)
- Besides `matrix`, `data.frames`, and `tibbles`, the package now also accepts `data.tables`.

## User visible changes

- `kernelshap()` is less picky regarding the output structure of `pred_fun()`.
- `kernelshap()` is less picky about the column structure of the background data `bg_X`. It should simply contain the columns of `X` (but can have more or in different order). The old behaviour was to launch an error if `colnames(X) != colnames(bg_X)`.
- The default `m = "auto"` has been changed from `trunc(20 * sqrt(p))` to `max(trunc(20 * sqrt(p)), 5 * p`. This will have an effect for cases where the number of features $p > 16$. The change will imply more robust results for large p.
- There were too many "ks_*()" functions to extract elements of a "kernelshap" object. They are now all deprecated and replaced by `ks_extract(, what = "S")`.
- Added package "MASS" to dependencies, see above.

## Bug fixes

- Depending on $m$ and $p$, the matrix inversion required in the constrained least-squares solution could fail. It is now replaced by `MASS::ginv()`, the Moore-Penrose pseudoinverse using `svd()`.

# kernelshap 0.1.0

This is the initial release.
