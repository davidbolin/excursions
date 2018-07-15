# Description #

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/excursions)](https://cran.r-project.org/package=excursions)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/excursions)](https://cranlogs.r-pkg.org/badges/grand-total/excursions)

`excursions` is an R package that contains functions that compute probabilistic excursion sets, contour credibility regions, and simultaneous confidence bands for latent Gaussian random processes and fields.

The theory of the methods used in the package are described in the papers [Bolin and Lindgren (2015)](http://onlinelibrary.wiley.com/doi/10.1111/rssb.12055/abstract), [Bolin and Lindgren (2016)](http://www.tandfonline.com/doi/full/10.1080/10618600.2016.1228537) , and [Bolin et al (2015)](http://www3.stat.sinica.edu.tw/statistica/j25n1/J25N120/J25N120.html).


# Manual #
A manual for the package (version 2.2.2) can be found [here](http://www.math.chalmers.se/~bodavid/software/excursions/excursions_manual_v2.pdf) and the code used in the manual can be downloaded [here](http://www.math.chalmers.se/~bodavid/software/excursions/code.zip).

# Versions #
The development version of the package contains new features and fixes that are not on CRAN, see `NEWS.md`. To install this version, see the instructions below.

# Installation instructions #
The latest CRAN release of the package can be installed directly from CRAN with `install.packages("excursions")`.
The latest stable version (which is sometimes slightly more recent than the CRAN version), can be installed by using the command
```r
devtools::install_bitbucket("excursions", "davidbolin", ref = "release")
```
in R. The development version can be installed using the command
```r
devtools::install_bitbucket("excursions", "davidbolin", ref = "default")
```

If you want to install the package using the `devtools::install_bitbucket`-method on Windows, you first need to install `Rtools` and add the paths to `Rtools` and `gcc` to the Windows `PATH` environment variable. This can be done for the current R session only using the commands
```r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH = paste(c(gcc, rtools, Sys.getenv("PATH")), collapse = ";"))
```
where the variables `rtools` and `gcc` need to be changed if `Rtool`s is not installed directly on `C:`.

# Repository branch workflows #
The package version format for released versions is `major.minor.bugfix`. All regular development should be performed on the `default` branch, or on feature branches derived from `default`. On the `default` branch, the vestion number is `major.minor.bugfix.9000`, where the first three component reflect the latest released version with changes present in the `default` branch. The version on the `release` branch is the bugfix version currently being prepared for stable release. Bugfixes may only be applied to the `release` branch, and must be merged into `default`.

  * Prepare a new release:
```
hg update release
hg merge default
## Update the version number as major.(minor+1).0
## Update the version in NEWS.md
hg commit -m "Start new release"
hg update default
hg merge release
## Resolve/update the version number conflict in favour of the release version, with extra .9000
## Update the version in NEWS.md
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
## (Resolve version conflict in favour of the release version, with extra .9000)
hg commit -m "Apply bugfix from release branch: ..."
## Optionally, prepare a stable version
```
  * Submit to CRAN
```
## If not already done, prepare a stable version
## Perform CRAN checks, if unsuccessful then stop, and do bugfixes
## Submit to CRAN
## If not accepted then stop, and do bugfixes
## If accepted, do
hg update default
## Update CRAN version in README.md
hg commit -m "Update CRAN version in README.md"
```
