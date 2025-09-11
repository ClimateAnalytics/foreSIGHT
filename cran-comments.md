## Test environments
* R version 4.4.2 (2023-10-31) on Windows 11 [devtools::check(args = "--as-cran")
]
* R-devel on win-builder [devtools::check_win_devel()]

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTES [devtools::check(args = "--as-cran") only]:
* NOTE about size of LazyData database (30.6MB). Compression has been enabled using `LazyDataCompression: xz`. The large size is due to included datasets used in examples and vignettes.
* NOTE about large number of packages. All packages in imports are considered to be essential for this package. 

## Downstream dependencies
There are no known reverse dependencies.

## Submission summary
* This is a new major release (v2.0.0).
* Major new features include:
  - Enabled use of flexible time steps (hourly to annual) instead of just daily
  - Updated attribute names to reflect different time steps and aggregation periods
  - Changed format of reference climate data to support multiple time steps
  - Implemented sub-daily and monthly stochastic weather generators (SWGs)
  - Added post-processing routines for modifying SWG output
  - Created functionality for tied attributes, which match changes to perturbed attributes
  - Added diagnostic plots for exploring how perturbed attributes affect other attributes
  - Supported multivariable attributes
  - Supported parallel processing in scenario generation 
* All examples, vignettes, and tests pass across platforms.
