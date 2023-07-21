Hello CRAN

My recent change introduced Latex errors in help files on MacOS and "oldrel". The reason are "text" instead of "textrm" tags in eqn{} blocks. I will need to fix this in two additional packages, but I will test it with this one.

https://cran.r-project.org/web/checks/check_results_kernelshap.html

## Checks

### Revdep

survex:

- OK: 1
- BROKEN: 0

### `check(manual = TRUE, cran = TRUE)`

❯ checking for future file timestamps ... NOTE
  unable to verify current time

❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  
### `check_win_devel()`

OK

### `check_rhub()` Notes

- kipping checking HTML validation: no command 'tidy' found
- Skipping checking math rendering: package 'V8' unavailable


