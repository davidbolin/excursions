# Excursions: Excursion Sets and Contour Credibility Regions for Random Fields

`excursions` contains functions that compute probabilistic excursion
sets, contour credibility regions, contour avoiding regions, contour map
quality measures, and simultaneous confidence bands for latent Gaussian
random processes and fields. A detailed manual can be found in the paper
Bolin, D and Lindgren, F (2018) *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, 86(5), 1–20.

## Details

The main functions in the package fall into three different categories
described below.

**Excursion sets, contour credibility regions, and contour avoiding
regions**

The main functions for computing excursion sets, contour credibility
regions, and contour avoiding regions are

- [`excursions()`](https://davidbolin.github.io/excursions/reference/excursions.md)
  :

  The main function for Gaussian models.

- [`excursions.inla()`](https://davidbolin.github.io/excursions/reference/excursions.inla.md)
  :

  Interface for latent Gaussian models estimated using INLA.

- [`excursions.mc()`](https://davidbolin.github.io/excursions/reference/excursions.mc.md)
  :

  Function for analyzing models that have been estimated using Monte
  Carlo methods.

The output from the functions above provides a discrete domain estimate
of the regions. Based on this estimate, the function
[`continuous()`](https://davidbolin.github.io/excursions/reference/continuous.md)
computes a continuous domain estimate.

The main reference for these functions is Bolin, D. and Lindgren, F.
(2015) *Excursion and contour uncertainty regions for latent Gaussian
models*, JRSS-series B, vol 77, no 1, pp 85-106.

**Contour map quality measures**

The package provides several functions for computing contour maps and
their quality measures. These quality measures can be used to decide on
an appropriate number of contours to use for the contour map.

The main functions for computing contour maps and the corresponding
quality measures are

- [`contourmap()`](https://davidbolin.github.io/excursions/reference/contourmap.md)
  :

  The main function for Gaussian models.

- [`contourmap.inla()`](https://davidbolin.github.io/excursions/reference/contourmap.inla.md)
  :

  Interface for latent Gaussian models estimated using INLA.

- [`contourmap.mc()`](https://davidbolin.github.io/excursions/reference/contourmap.mc.md)
  :

  Function for analyzing models that have been estimated using Monte
  Carlo methods.

Other noteworthy functions relating to contourmaps are
[`tricontour()`](https://davidbolin.github.io/excursions/reference/tricontour.md)
and
[`tricontourmap()`](https://davidbolin.github.io/excursions/reference/tricontour.md),
which compute contour curves for functinos defined on triangulations, as
well as
[`contourmap.colors()`](https://davidbolin.github.io/excursions/reference/contourmap.colors.md)
which can be used to compute appropriate colors for displaying contour
maps.

The main reference for these functions is Bolin, D. and Lindgren, F.
(2017) *Quantifying the uncertainty of contour maps*, Journal of
Computational and Graphical Statistics, 26:3, 513-524.

**Simultaneous confidence bands**

The main functions for computing simultaneous confidence bands are

- [`simconf()`](https://davidbolin.github.io/excursions/reference/simconf.md)
  :

  Function for analyzing Gaussian models.

- [`simconf.inla()`](https://davidbolin.github.io/excursions/reference/simconf.inla.md)
  :

  Function for analyzing latent Gaussian models estimated using INLA.

- [`simconf.mc()`](https://davidbolin.github.io/excursions/reference/simconf.mc.md)
  :

  Function for analyzing models estimated using Monte Carlo methods.

- [`simconf.mixture()`](https://davidbolin.github.io/excursions/reference/simconf.mixture.md)
  :

  Function for analyzing Gaussian mixture models.

The main reference for these functions is Bolin et al. (2015)
*Statistical prediction of global sea level from global temperature*,
Statistica Sinica, Vol 25, pp 351-367.

## See also

Useful links:

- <https://github.com/davidbolin/excursions>

- <https://davidbolin.github.io/excursions/>

- Report bugs at <https://github.com/davidbolin/excursions/issues>

## Author

**Maintainer**: David Bolin <davidbolin@gmail.com>
([ORCID](https://orcid.org/0000-0003-2361-5465))

Authors:

- Finn Lindgren <finn.lindgren@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-5833-2011))

Other contributors:

- Suen Man Ho <M.H.Suen@sms.ed.ac.uk> \[contributor\]
