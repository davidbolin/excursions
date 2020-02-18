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
remotes::install_bitbucket("excursions/davidbolin", ref = "master")
```
in R. The development version can be installed using the command
```r
remotes::install_bitbucket("excursions/davidbolin", ref = "devel")
```

If you want to install the package using the `remotes::install_bitbucket`-method on Windows, you first need to install `Rtools` and add the paths to `Rtools` and `gcc` to the Windows `PATH` environment variable. This can be done for the current R session only using the commands
```r
rtools = "C:\\Rtools\\bin"
gcc = "C:\\Rtools\\gcc-4.6.3\\bin"
Sys.setenv(PATH = paste(c(gcc, rtools, Sys.getenv("PATH")), collapse = ";"))
```
where the variables `rtools` and `gcc` need to be changed if `Rtool`s is not installed directly on `C:`.

# Repository branch workflows #
The package version format for released versions is `major.minor.bugfix`. All regular development should be performed on the `devel` branch or in a feature branch, managed with `git flow feature`. On the `devel` branch, the vestion number is `major.minor.bugfix.9000`, where the first three components reflect the latest released version with changes present in the `default` branch. Bugfixes should be applied via the `git flow bugfix` and `git flow hotfix` methods, as indicated below. For `git flow` configuration, use `master` as the stable master branch, `devel` as the develop branch, and `v` as the version tag prefix. See [the `git flow` tutorial](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more information.

For non `master` and `devel` branches that collaborators need access to (e.g. release branches, feature branches, etc, use the `git flow publish` mechanism).

  * Prepare a new stable release with CRAN submission:
```
git flow release start major.(minor+1).0
## Update the DESCRIPTION version number as major.(minor+1).0
## Update the version in NEWS.md
## Commit the changes
## At this point, see the CRAN submission section below.
git flow release finish 'VERSION'
## Resolve/update the DESCRIPTION and NEWS.md version number conflict
## in favour of the released version, with extra .9000, e.g. with
## the help of  git mergetool
## Add a new version section in NEWS.md
## Commit the merge
```
  * Do a hotfix (branch from stable master; use bugfix for release branch bugfixes):
```
git flow hotfix start hotfix_branch_name
## Do the bugfix, update the verison number major.minor.(bugfix+1), and commit
## Optionally, do CRAN submission
git flow hotfix finish hotfix_branch_name
## Resolve merge conflicts (hopefully mostly due to version numbers)
```
  * CRAN submission
```
## Perform CRAN checks (usually on the release branch version)
## If unsuccessful then do bugfixes with increasing bugfix version, until ok
## Submit to CRAN
## If not accepted then do more bugfixes and repeat
```
