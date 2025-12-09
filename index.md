# The `excursions` package

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/excursions)](https://cran.r-project.org/package=excursions)
[![CRAN_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/excursions)](https://cranlogs.r-pkg.org/badges/grand-total/excursions)

The ability to find regions where a stochastic process exceeds a certain
level, or is significantly different from some reference level, is
important in several areas of application. A related problem is
uncertainty quantification of contour curves and more generally of
contour maps, which are often used to display estimates of continuous
surfaces. The number of contours used in a contour map should typically
reflect the uncertainty in the estimate, since one should be allowed to
draw many contours if the uncertainty of the estimated surface is low
and fewer contours if the uncertainty is high. The ability to quantify
the uncertainty in the contour map is important if one should be able to
choose the number of contours in a rigorous way.

`excursions` contains functions for solving these problems for latent
Gaussian models (LGMs), which is a large model class that is widely used
in applications. The computational methods are especially well-suited
for models where the latent field has Markov properties. Solving the
problems involves computing high-dimensional Gaussian integrals, which
can be done more efficiently if Markov properties can be utilized. With
the ability to efficiently compute Gaussian integrals, one can also
compute simultaneous credible bands for latent Gaussian processes, and
more generally for mixtures of Gaussian processes.

The package supports at least three ways of specifying the model that
should be analyzed. The standard method for purely Gaussian models is to
specify the model by providing the parameters of the Gaussian process.
For more general models, the input can either be given as Monte Carlo
simulations of the process or as the result from an analysis using the
[R-INLA](https://davidbolin.github.io/excursions//articles/inla.html "`INLA` vignette")
or
[inlabru](https://davidbolin.github.io/excursions//articles/inlabru.html "`inlabru` vignette")
packages.

See the [Getting
started](https://davidbolin.github.io/rSPDE//articles/excursions.html "An introduction to the excursions package")
vignette for an introduction to the main functions of the package. A
complete list of functions is contained in the [`excursions`
documentation](https://davidbolin.github.io/rSPDE/reference/index.html "`excursions` documentation.")
and a brief introduction to the theory is provided in the
[Methodology](https://davidbolin.github.io/rSPDE//articles/theory.html "Methodology")
vignette.

The `INLA` interface is introduced in the
[R-INLA](https://davidbolin.github.io/excursions//articles/inla.html "`INLA` vignette")
vignette and the `inlabru` interface is introduced in the
[inlabru](https://davidbolin.github.io/excursions//articles/inlabru.html "`inlabru` vignette")
vignette.

# Installation instructions

The latest CRAN release of the package can be installed directly from
CRAN with `install.packages("excursions")`. The latest stable version
(which is sometimes slightly more recent than the CRAN version), can be
installed by using the command

``` r
remotes::install_github("davidbolin/excursions", ref = "stable")
```

in R. The development version can be installed using the command

``` r
remotes::install_github("davidbolin/excursions", ref = "devel")
```
