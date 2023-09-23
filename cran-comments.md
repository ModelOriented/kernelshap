# Resubmission II

I have now removed the non-registered DOI from the description file.

This was Uwe's comment:

   Found the following (possibly) invalid DOIs:
     DOI: 10.5555/3295222.3295230
       From: DESCRIPTION
       Status: 404
       Message: Not Found

and looking for the DOI shows it is not registered.
Even when going to
https://dl.acm.org/doi/10.5555/3295222.3295230
and clicking on the "Publisher Site" link leads us into nirvana.

So I guess the is ill registered and you have to revert the change. My 
apologoies, I had not checked whether the DOI given is registreed or not.

Best,


# Resubmission I

I stumbled over a DOI: Uwe gently pointed this out:

"
   The Description field contains
     <https://dl.acm.org/doi/10.5555/3295222.3295230>, and Covert and Lee
   Please use permanent DOI markup for linking to publications as in 
<doi:prefix/suffix>.
"

This resubmission fixes this.

# Original message

Hello CRAN

This is a small maintenance release only.

Thanks a lot

Michael

## Checks

### Revdep

survex 1.1.3                                                                             
- OK: 1
- BROKEN: 0

### Local check with innocent NOTE

❯ checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  
### `check_win_devel()` NOTE

- R Under development (unstable) (2023-09-11 r85126 ucrt)

### `check_rhub()` NOTES

-> Note sure where the 403 problem comes from. Is it relevant?

Found the following (possibly) invalid URLs:
  URL: https://dl.acm.org/doi/10.5555/3295222.3295230
    From: DESCRIPTION
    Status: 403
    Message: Forbidden
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable
