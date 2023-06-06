# kernelshap 0.3.8

## API improvements

- Multi-output case: column names of predictions are now used as list names of the resulting `S` and `SE` lists.

## Maintenance

- Added explanation of sampling Kernel SHAP to help file

# kernelshap 0.3.7

## Maintenance

- Fixed problem in Latex math for MacOS.

# kernelshap 0.3.6

## Maintenance

- Improved help files and README

# kernelshap 0.3.5

## Maintenance

- New contributor: Przemyslaw Biecek - welcome on board!
- My new cozy home: https://github.com/ModelOriented/kernelshap
- Webpage created with "pkgdown"
- Introduced Github workflows
- More unit tests

## Small visible changes

- Removed the `ks_extract()` function. It was designed to extract objects like the matrix `S` of SHAP values from the resulting "kernelshap" object `x`. We feel that the standard extraction options (`x$S`, `x[["S"]]`, or `getElement(x, "S")`) are sufficient.
- Adding $(n \times K)$ matrix of predictions to the output, where $n$ is the number of rows in the explainer data `X`, and $K$ is the dimension of a single prediction (usually 1).
- Setting `verbose = FALSE` now does not suppress the warning on too large background data anymore. Use `suppressWarnings()` instead.

# kernelshap 0.3.4

## Documentation

- New logo
- Better package description
- Better README

# kernelshap 0.3.3

## Less dependencies

- Removed dependency "dorng". This might have an impact on the seeding if in parallel mode.
- Removed dependency "MASS"

# kernelshap 0.3.2

## Documentation

- Rewritten README and examples to better show the role of the background data.

## Bug fixes

- When `bg_X` contained more columns than `X`, unflexible prediction functions could fail when being applied to `bg_X`. 

# kernelshap 0.3.1

## Changes

- New argument `feature_names` allows to specify the features to calculate SHAP values for. The default equals to `colnames(X)`. This should be changed only in situations when `X` (the dataset to be explained) contains non-feature columns.
- The background dataset can now consist of a single row only. This is useful in situations with natural "off" value such as for image data or for models that can naturally deal with missing values.


# kernelshap 0.3.0

## Major improvements

### Exact calculations

Thanks to David Watson, exact calculations are now also possible for $p>5$ features. By default, the algorithm uses exact calculations for $p \le 8$ and a hybrid strategy otherwise, see the next section. At the same time, the exact algorithm became much more efficient.

A word of caution: Exact calculations mean to create $2^p-2$ on-off vectors $z$ (cheap step) and evaluating the model on a whopping $(2^p-2)N$ rows, where $N$ is the number of rows of the background data (expensive step). As this explodes with large $p$, we do not recommend the exact strategy for $p > 10$.

### Hybrid strategy

The iterative Kernel SHAP sampling algorithm of Covert and Lee (2021) [1] works by randomly sample $m$ on-off vectors $z$ so that their sum follows the SHAP Kernel weight distribution (renormalized to the range from $1$ to $p-1$). Based on these vectors, many predictions are formed. Then, Kernel SHAP values are derived as the solution of a constrained linear regression, see [1] for details. This is done multiple times until convergence.

A drawback of this strategy is that many (at least 75%) of the $z$ vectors will have $\sum z \in \{1, p-1\}$, producing many duplicates. Similarly, at least 92% of the mass will be used for the $p(p+1)$ possible vectors with $\sum z \in \{1, 2, p-1, p-2\}$ etc. This inefficiency can be fixed by a hybrid strategy, combining exact calculations with sampling. 
The hybrid algorithm has two steps:

1. Step 1 (exact part): There are $2p$ different on-off vectors $z$ with $\sum z \in \{1, p-1\}$, covering a large proportion of the Kernel SHAP distribution. The degree 1 hybrid will list those vectors and use them according to their weights in the upcoming calculations. Depending on $p$, we can also go a step further to a degree 2 hybrid by adding all $p(p-1)$ vectors with $\sum z \in \{2, p-2\}$ to the process etc. The necessary predictions are obtained along with other calculations similar to those in [1].
2. Step 2 (sampling part): The remaining weight is filled by sampling vectors $z$ according to Kernel SHAP weights renormalized to the values not yet covered by Step 1. Together with the results from Step 1 - correctly weighted - this now forms a complete iteration as in Covert and Lee (2021). The difference is that most mass is covered by exact calculations. Afterwards, the algorithm iterates until convergence. The output of Step 1 is reused in every iteration, leading to an extremely efficient strategy.

The default behaviour of `kernelshap()` is as follows:

- $p \le 8$: Exact Kernel SHAP (with respect to the background data)
- $9 \le p \le 16$: Degree 2 hybrid
- $p > 16$: Degree 1 hybrid
- $p = 1$: Exact Shapley values

It is also possible to use a pure sampling strategy, see Section "User visible changes" below. While this is usually not advisable compared to a hybrid approach, the options of `kernelshap()` allow to study different properties of Kernel SHAP and doing empirical research on the topic.

Kernel SHAP in the Python implementation "shap" uses a quite similar hybrid strategy, but without iterating. The new logic in the R package thus combines the efficiency of the Python implementation with the convergence monitoring of [1].

[1] Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.

## User visible changes

- The default value of `m` is reduced from $8p$ to $2p$ except when `hybrid_degree = 0` (pure sampling).
- The default value of `exact` is now `TRUE` for $p \le 8$ instead of $p \le 5$.
- A new argument `hybrid_degree` is introduced to control the exact part of the hybrid algorithm. The default is 2 for $4 \le p \le 16$ and degree 1 otherwise. Set to 0 to force a pure sampling strategy (not recommended but useful to demonstrate superiority of hybrid approaches).
- The default value of `tol` was reduced from 0.01 to 0.005.
- The default of `max_iter` was reduced from 250 to 100.
- The order of some of the arguments behind the first four has been changed.
- Paired sampling no longer duplicates `m`.
- Thanks to Mathias Ambuehl, the random sampling of z vectors is now fully vectorized.
- The output of `print()` is now more slim. 
- A new `summary()` function shows more infos.

## Other changes

- The resulting object now contains `m_exact` (the number of on-off vectors used for the exact part), `prop_exact` (proportion of mass treated in exact fashion), `exact` flag, and `txt` (the info message when starting the algorithm).

## Bug fixes

- Predictions of `mgcv::gam()` would cause an error in `check_pred()` (they are 1D-arrays).
- Fixed small mistakes in the examples of the README (mlr3 and mgcv).

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
