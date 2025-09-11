## Test environments
* R version 4.4.2 (2023-10-31) on Windows 11 [devtools::check(args = "--as-cran")]
* R-devel on win-builder [devtools::check_win_devel()]

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTES [devtools::check(args = "--as-cran") only]:

NOTE 1: Imports includes 25 non-default packages.
The package imports a larger number of packages due to its integration of multiple statistical, hydrological, and visualization tools to support a comprehensive climate scenario evaluation framework. Further reduction is being considered for future versions.

NOTE 2: Installed size is 29.2Mb (data: 26.4Mb, doc: 1.2Mb)
The majority of the package size is due to example datasets (in data/) used to illustrate and validate the package functionality. These datasets are important for reproducibility, examples, and vignettes. The size has been minimized where feasible without compromising usability.

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
