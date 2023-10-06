
Results from rhub::check_for_cran(platforms=c('windows-x86_64-devel','ubuntu-gcc-release','fedora-clang-devel'))

2 notes (both seem unavoidable)

─  Building package
─  Uploading package
─  Preparing build, see status at
   https://builder.r-hub.io/status/foreSIGHT_1.2.0.tar.gz-394bdd12530741f9a02f0e19632a9f97
   https://builder.r-hub.io/status/foreSIGHT_1.2.0.tar.gz-f92be62eb5374004a1034ce314b4bc84
   https://builder.r-hub.io/status/foreSIGHT_1.2.0.tar.gz-31f23149fe10499ca93cd6c5eee005c2
─  Build started
─  Creating new user
─  Downloading and unpacking package file
─  Querying package dependencies
─  Installing package dependencies
─  Running R CMD check
   setting _R_CHECK_FORCE_SUGGESTS_ to false
   setting R_COMPILE_AND_INSTALL_PACKAGES to never
   setting R_REMOTES_STANDALONE to true
   setting R_REMOTES_NO_ERRORS_FROM_WARNINGS to true
   setting _R_CHECK_FORCE_SUGGESTS_ to true
   setting _R_CHECK_CRAN_INCOMING_USE_ASPELL_ to true
   'getOption("repos")' replaces Bioconductor standard repositories, see
   'help("repositories", package = "BiocManager")' for details.
   Replacement repositories:
       CRAN: https://cloud.r-project.org
─  using log directory 'C:/Users/USERSALQidPpfB/foreSIGHT.Rcheck'
─  using R Under development (unstable) (2023-07-21 r84722 ucrt)
─  using platform: x86_64-w64-mingw32
─  R was compiled by
       gcc.exe (GCC) 12.2.0
       GNU Fortran (GCC) 12.2.0
─  running under: Windows Server 2022 x64 (build 20348)
─  using session charset: UTF-8
─  using option '--as-cran'
✔  checking for file 'foreSIGHT/DESCRIPTION'
─  this is package 'foreSIGHT' version '1.2.0'
─  checking CRAN incoming feasibility ... [11s] NOTE
   Maintainer: 'David McInerney <david.mcinerney@adelaide.edu.au>'
   
   New maintainer:
     David McInerney <david.mcinerney@adelaide.edu.au>
   Old maintainer(s):
     Bree Bennett <bree.bennett@adelaide.edu.au>
✔  checking package namespace information
✔  checking package dependencies
✔  checking if this is a source package
✔  checking if there is a namespace
✔  checking for executable files
✔  checking for hidden files and directories
✔  checking for portable file names
─  checking whether package 'foreSIGHT' can be installed ... [44s] OK
─  used C compiler: 'gcc.exe (GCC) 12.2.0'
─  used C++ compiler: 'G__~1.EXE (GCC) 12.2.0'
✔  checking installed package size
✔  checking package directory
✔  checking for future file timestamps
✔  checking 'build' directory
✔  checking DESCRIPTION meta-information
✔  checking top-level files
✔  checking for left-over files
✔  checking index information
✔  checking package subdirectories
✔  checking R files for non-ASCII characters
✔  checking R files for syntax errors
✔  checking whether the package can be loaded
✔  checking whether the package can be loaded with stated dependencies
✔  checking whether the package can be unloaded cleanly
✔  checking whether the namespace can be loaded with stated dependencies
✔  checking whether the namespace can be unloaded cleanly
✔  checking loading without being on the library search path (12.6s)
✔  checking startup messages can be suppressed
✔  checking use of S3 registration
✔  checking dependencies in R code
✔  checking S3 generic/method consistency
✔  checking replacement functions
✔  checking foreign function calls
─  checking R code for possible problems ... [30s] OK (41.5s)
✔  checking Rd files
✔  checking Rd metadata
✔  checking Rd line widths
✔  checking Rd cross-references
✔  checking for missing documentation entries
✔  checking for code/documentation mismatches
✔  checking Rd \usage sections (5.9s)
✔  checking Rd contents
✔  checking for unstated dependencies in examples
✔  checking contents of 'data' directory
✔  checking data for non-ASCII characters
✔  checking LazyData
✔  checking data for ASCII and uncompressed saves
✔  checking line endings in C/C++/Fortran sources/headers
✔  checking pragmas in C/C++ headers and code
✔  checking compilation flags used
✔  checking compiled code (18.5s)
✔  checking installed files from 'inst/doc'
✔  checking files in 'vignettes'
─  checking examples ... [17s] NOTE (19.6s)
   Examples with CPU (user + system) or elapsed time > 5s
                             user system elapsed
   plotPerformanceSpaceMulti 5.39   0.09    5.48
✔  checking for unstated dependencies in 'tests'
─  checking tests
✔  Running 'testthat.R' (13.8s)
✔  checking for unstated dependencies in vignettes
✔  checking package vignettes in 'inst/doc'
─  checking re-building of vignette outputs ... [112s] OK (1m 51.3s)
─  checking PDF version of manual ... [18s] OK (9.7s)
N  checking for non-standard things in the check directory (12.7s)
✔  checking HTML version of manual
   Found the following files/directories:
N  checking for detritus in the temp directory
   Found the following files/directories:
     ''NULL''
     'lastMiKTeXException'
   
─  Done with R CMD check
─  Cleaning up files and user