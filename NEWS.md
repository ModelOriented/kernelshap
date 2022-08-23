# kernelshap 0.1.900 DEVEL

## Major improvements

- The package now supports multi-output predictions. Hurray!
- Major speed-up of data operations and minimization of calls to the prediction function.

## Major changes

- The default `m = "auto"` has been changed from `trunc(20 * sqrt(p))` to `max(trunc(20 * sqrt(p)), 5 * p`. This will have an effect for cases where the number of features $p > 16$. The change will imply more robust results for large p.

## Bug fixes

- Depending on $m$ and $p$, a matrix inversion required in the constrained least-squares solution could fail. It is now replaced by `MASS::ginv()`, the Moore-Penrose pseudoinverse.

## Minor changes

- There were too many "ks_*()" functions to extract elements of a "kernelshap" object. They are now all deprecated and replaced by `ks_extract(, what = "S")`.
- Added package "Mass" to dependencies, see above.

# kernelshap 0.1.0

This is the initial release.
