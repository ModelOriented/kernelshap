# kernelshap 0.2.0.900 DEVEL

## Major improvements

### Exact calculations

Thanks to David Watson, exact calculations are now possible for $p>5$ features. By default, the algorithm now uses exact calculations for $p \le 8$ and a hybrid strategy otherwise, see next section. Note that the exact case will create $2^p-2$ on-off vectors $z$. Consequently, predictions will be done on $(2^p-2)N$ rows, where $N$ is the number of rows of the background data. As this explodes with large $p$, we do not recommend it for $p > 10$.

### Hybrid strategy

The iterative Kernel SHAP sampling algorithm of Covert and Lee (2021) [1] works by randomly sample $m$ on-off vectors $z$ so that their sum follows the SHAP Kernel weight distribution (renormalized to the range from $1$ to $p-1$). Based on these vectors, many predictions are formed and Kernel SHAP values are derived, see [1] for details. This is done multiple times until convergence.

A drawback of this strategy is that many (at least 75%) of the $z$ vectors will have $\sum z \in \{1, p-1\}$, producing many duplicates. Similarly, at least 92% of the mass will be used for the $p(p+1)$ possible vectors with $\sum z \in \{1, 2, p-1, p-2\}$ etc. This inefficiency can be removed by a hybrid strategy, combining exact calculations with sampling. 

The hybrid algorithm works in two steps:

1. Step 1 (exact part): All z vectors with $\sum z \in \{1, p-1\}$ are created (degree 1 hybrid). Then, if `m_exact_max` is sufficiently large, also all vectors with $\sum z \in \{2, p-2\}$ are added (degree 2) etc. until at most `m_exact_max` on-off vectors are generated with total Kernel SHAP weight $w_{exact}$. The necessary predictions are obtained along with other calculations similar to those in [1].
2. Step 2 (sampling part): The remaining weight $1 - w_{exact}$ is filled by sampling vectors $z$ according to Kernel SHAP weights renormalized to the values not yet covered by Step 1. Together with the results from Step 1, correctly weighted, this now forms a complete iteration like in Covert and Lee (2021) with the difference that most mass is covered by exact calculations. Then, the algorithm iterates until convergence. The output of Step 1 is reused in every iteration, leading to an extremely efficient strategy that usually converges in the minimal number of iterations (2).

The default behaviour of `kernelshap()` is as follows:

- $p \le 8$: Exact Kernel SHAP (with respect to the background data)
- $9 \le p \le 16$: Degree 2 hybrid
- $p > 16$: Degree 1 hybrid
- $p = 1$: Exact Shapley values

It is also possible to use a pure sampling strategy, see Section "User visible changes" below.

Kernel SHAP in the Python implementation "shap" uses a quite similar hybrid strategy, but without iterating. The new logic in the R package thus combines the efficiency of the Python implementation with the convergence monitoring of [1].

[1] Ian Covert and Su-In Lee. Improving KernelSHAP: Practical Shapley Value Estimation Using Linear Regression. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:3457-3465, 2021.

## User visible changes

- The default value of `m` is now $\min(256, 8p)$`.
- The default value of `exact` is now `TRUE` if $p <= 8$ and no `hybrid_degree` is set.
- A new argument `hybrid_degree` is introduced to control the exact part of the hybrid algorithm. The default, `NULL`, ensures hybrid degree 2 up to $p\le 16$ and degree 1 for $p > 16$. Set to 0 to force a pure sampling strategy.
- The default value of `tol` was reduced from 0.01 to 0.005.
- The default of `max_iter` was reduced from 250 to 25.
- The order of some of the arguments behind the first four has been changed.
- Paired sampling no longer duplicates `m`.
- Thanks to Mathias Ambuehl, the random sampling of z vectors is now fully vectorized.

## Bug fixes

- Predictions of `mgcv::gam()` would cause an error in `check_pred()` (they are 1D-arrays).
- Fixed small mistakes in the examples of the README (mlr3 and mgcv).

## New contributor

- David Watson

# kernelshap 0.1.0

This is the initial release.
