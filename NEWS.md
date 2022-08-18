# kernelshap 0.1.900 DEVEL

## Major improvements

- The package now supports multi-output predictions. Hurray!
- Major speed-up of data operations and minimization of calls to the prediction function.

## Minor changes

- There were too many "ks_*()" functions to extract elements of a "kernelshap" object. They are now all deprecated and replaced by `ks_extract(, what = "S")`.

# kernelshap 0.1.0

This is the initial release.
