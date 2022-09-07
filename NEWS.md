# kernelshap 0.2.0.900 DEVEL

## Major improvements

- Thanks to David Watson, exact calculations are now possible for $p>5$ features, although we do not recommend it for $p>8$. In the exact case with $p > 1$, `kernelshap()` does predicions on data sets with $(2^p-1)N$ rows, where $N$ is the number of rows of the background data. This explodes with large $p$.

## User visible changes

- The default for `exact` argument is now `"auto"` instead of `TRUE`. This will use exact
calculations up to $p=8$ features. For $p>8$, sampling is used.

## Bug fixes

- Predictions of `mgcv::gam()` would cause an error in `check_pred()` (they are 1D-arrays).
- Fixed small mistakes in Readme (mlr3 and mgcv).

# kernelshap 0.2.0

## Breaking change

The interface of `kernelshap()` has been revised. Instead of specifying a prediction function, it suffices now to pass the fitted model object. The default `pred_fun` is now `stats::predict`, which works in most cases. Some other cases are catched via model class ("ranger" and mlr3 "Learner"). The `pred_fun` can be overwritten by a function of the form `function(object, X, ...)`. Additional arguments to the prediction function are passed via `...` of `kernelshap()`.

Some examples:

- Logistic regression (logit scale): `kernelshap(fit, X, bg_X)`
- Logistic regression (probabilities): `kernelshap(fit, X, bg_X, type = "response")`
- Linear regression with logarithmic response, but evaluated on original scale: Here, the default predict function needs to be overwritten: `kernelshap(fit, X, bg_X, pred_fun = function(m, X) exp(predict(m, X)))`

## Major improvements

- `kernelshap()` has received a more intuitive interface, see breaking change above.
- The package now supports multidimensional predictions. Hurray!
- Thanks to David Watson, parallel computing is now supported. The user needs to set up the parallel backend before calling `kernelshap()`, e.g., using the "doFuture" package, and then set `parallel = TRUE`. Especially on Windows, sometimes not all global variables or packages are loaded in the parallel instances. These can be specified by `parallel_args`, a list of arguments passed to `foreach()`.
- Even without parallel computing, `kernelshap()` has become much faster.
- For $2 \le p \le 5$ features, the algorithm now returns exact Kernel SHAP values with respect to the given background data. (For $p = 1$, exact *Shapley values* are returned.)
- Direct handling of "tidymodels" models.

## User visible changes

- Besides `matrix`, `data.frame`s, and `tibble`s, the package now also accepts `data.table`s (if the prediction function can deal with them).
- `kernelshap()` is less picky regarding the output structure of `pred_fun()`.
- `kernelshap()` is less picky about the column structure of the background data `bg_X`. It should simply contain the columns of `X` (but can have more or in different order). The old behaviour was to launch an error if `colnames(X) != colnames(bg_X)`.
- The default `m = "auto"` has been changed from `trunc(20 * sqrt(p))` to `max(trunc(20 * sqrt(p)), 5 * p`. This will have an effect for cases where the number of features $p > 16$. The change will imply more robust results for large p.
- There were too many "ks_*()" functions to extract elements of a "kernelshap" object. They are now all deprecated and replaced by `ks_extract(, what = "S")`.
- Added "MASS", "doRNG", and "foreach" to dependencies.

## Bug fixes

- Depending on $m$ and $p$, the matrix inversion required in the constrained least-squares solution could fail. It is now replaced by `MASS::ginv()`, the Moore-Penrose pseudoinverse using `svd()`.

## New contributor

- David Watson

# kernelshap 0.1.0

This is the initial release.
