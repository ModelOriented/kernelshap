# kernelshap 0.2.0 DEVEL

## Major improvements

- The package now supports multi-output predictions. Hurray!
- Major speed-up of data operations.

## Minor changes

- Added "data.table" package to "Suggests".
- Major simplification of linear algebra.
- There were too many "ks_*()" functions to extract elements of a "kernelshap" object. They are now all deprecated and replaced by `ks_extract(, what = "S")`.

# kernelshap 0.1.0

This is the initial release.
