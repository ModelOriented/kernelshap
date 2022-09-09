# kernelshap 0.2.0.900 DEVEL

## Major improvements

- Thanks to David Watson, exact calculations are now possible for $p>5$ features. Since predicions are done on data sets with $(2^p-2)N$ rows, where $N$ is the number of rows of the background data, we do not recommend the exact algorithm for $p > 10$.
- The arguments `exact` and `paired_sampling` have been entangled and replaced by a new argument `sampling_strategy`. Its default is "auto", which will use `sampling_strategy = 'exact'` for up to eight features, and `sampling_strategy = 'hybrid'` for more features.
- Inspired by the Python implementation of Scott Lundberg, we have introduced a new sampling strategy "hybrid". Instead of randomly drawing $m$ binary "on-off" vectors $z$ from the Kernel SHAP distribution, we procede in two steps: In the first (exact) step, all vectors $z$ with $sum(z) \in \{1, p-1\}$ are listed (degree 1 hybrid). If $m$ is sufficiently large, also all $z$ with $sum(z) \in \{2, p-2\}$ are added (degree 2). In the second step, the remaining $z$ are sampled from the (renormalized) kernel SHAP distribution according to the paired strategy. $z$ vectors are weighted to match the Kernel SHAP distribution. Convergence behaviour of this strategy is expected to be much better compared to a pure sampling technique. The reason is that the Kernel SHAP distribution has very high weights for values $1$ or $p-1$ (and also for values $2$ and $p-2$). A pure sampling strategy would mainly pick $z$ vectors with $sum(z) \in \{1, p-1\}$, producing many duplicated $z$ and leaving only a few other $z$. The hybrid strategy corrects this inefficient behaviour by enumerating the $z$ with high weights and sampling the others.

## Not backward compatible changes

- Arguments `exact` and `paired_sampling` are depreciated in favour of the new argument `sampling_strategy`, see above. They will be removed in version 0.4.0.
- The new default, `sampling_strategy = "auto"`, will use exact calculations for $p \le 8$ and the new hybrid strategy (see above) otherwise.
- The default for `m` is now \code{NULL} instead of "auto". Also the value of `m` in this case has changed to `max(30, pp * (pp + 5), 4 * p)`, where `pp = min(p, 14)`. This is in line with the hybrid strategy.
- Instead of doubling `m` with paired sampling, we now double the default value of `m` if `m = "auto"`. As a consequence, paired and unpaired sampling will use the same number of on-off vectors per iteration, which improves comparability of the two approaches.

## Bug fixes

- Predictions of `mgcv::gam()` would cause an error in `check_pred()` (they are 1D-arrays).
- Fixed mistakes in the examples of the Readme (mlr3 and mgcv).

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
