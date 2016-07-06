# Description #
Excursions is an R package that contains functions that compute probabilistic excursion sets, contour credibility regions, and simultaneous confidence bands for latent Gaussian random processes and fields. 

The theory of the methods used in the package are described in [Bolin and Lindgren (2014)](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12055/abstract). A manual for the package (version 1.1, so slightly outdated at the moment) can be found [here](http://www.math.chalmers.se/~bodavid/software/excursions/excursions_manual.pdf) and the code used in the manual can be downloaded [here](http://www.math.chalmers.se/~bodavid/software/excursions/code.zip). 

# Versions #
The current stable version of the package on CRAN is 2.2.1. The development version of the package contains new features and fixes that are not on CRAN. To install this version, see the instructions below. 

# Installation instructions #
The stable version of the package can be installed directly from CRAN, or by using the command
```
#!r

devtools::install_bitbucket("excursions","davidbolin",ref="stable")
```
in R. The development version can be installed using the command 
```
#!r

devtools::install_bitbucket("excursions","davidbolin",ref="default")
```
The development version of the package is called excursionsdevel, so it is loaded by running

```
#!r
library(excursionsdevel)
```

If you want to install the package using the install_bitbucket-method on Windows, you first need to install Rtools and add the paths to Rtools and gcc to the Windows PATH environment variable. This can be done for the current R session only using the commands
```
#!r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH=paste(c(gcc,rtools, Sys.getenv("PATH")),collapse=";"))
```
where the variables rtools and gcc need to be changed if Rtools is not installed directly on C:.

# Repository branch workflows #
The package version format is `major.minor.bugfix`. All regular development should be performed on the `default` branch, or on feature branches derived from `default`. The version on the `default` branch is the next version to be prepared for release. The version on the `release` branch is the bugfix version currently being prepared for stable release. Bugfixes may only be applied to the `release` branch, and must be merged into `default`.

   * Prepare a new release:
```
hg update release
hg merge default
## (Resolve any version number conflict in favour of the default version)
hg commit -m "Start new release"
hg update default
## Update the version number as major.(minor+1).0
hg commit -m "Start next development version"
```
  * Prepare a stable version:
```
hg update release
## Update the date in the DESCRIPTION
hg commit -m "Update release date"
## Perform CRAN checks, if unsuccessful then stop, and do a bugfix instead.
hg update stable
hg merge release
## (Resolve package name conflicts in favour of stable name (excursions), and README.md removed)
hg commit -m "New stable version"
hg tag vX.X.X
## X.X.X = major.minor.bugfix from DESCRIPTION
```
  * Do a bugfix:
```
hg update release
## Update the version number to major.minor.(bugfix+1)
## Fix the bug
hg commit
hg update default
hg merge release
## (Resolve version conflict in favour of the default version)
hg commit -m "Apply bugfix from release branch: ..."
## Optionally, prepare a stable version
```
  * Submit to CRAN
```
## If not already done, prepare a stable version
## Submit to CRAN
## If not accepted, do bugfixes and resubmit
## If accepted, do
hg update default
## Update CRAN version in README.md
hg commit -m "Update CRAN version in README.md"
```