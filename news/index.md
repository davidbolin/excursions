# Changelog

## excursions 2.5.11

- Update link to external inla vignette in documentation

## excursions 2.5.10

- Corrections to various links in documentation

## excursions 2.5.9

- Add -DNPRINT to CAMD build, and remove fflush(stdout) use, to avoid
  fprint and stdout usage
- Fix issue with identity matrices, ensuring they are not treated as
  all-zero
- Remove obsolete and unneded includes of Rdefines.h and R_ext/PrtUtil.h

## excursions 2.5.8

CRAN release: 2023-11-30

- Minor updates to C code to avoid warning with gcc13 and clang17
  compilers

## excursions 2.5.7

CRAN release: 2023-10-23

- Update to make sure that max.threads properly limits the number of
  threads
- Minor update to the documentation

## excursions 2.5.6

- Remove rgdal suggest
- Update some examples to use fmesher instead of INLA

## excursions 2.5.5

CRAN release: 2023-01-25

- Add R compiler configuration extraction to configure.ac
- Add support for results computed with the new default `compact` mode
  in INLA

## excursions 2.5.4

CRAN release: 2022-10-14

- Remove local gsl code, and rely on SystemRequirements instead

## excursions 2.5.3

CRAN release: 2022-08-25

- Add support for INLA results computed with the `experimental` option.
- Avoid deprecated Matrix (\>=1.4-2) class coercion methods

## excursions 2.5.2

CRAN release: 2021-09-17

- Add `return.marginals.predictor = TRUE` option to tests and examples,
  for new INLA compatibility
- Move inclusion of `omp.h` to before R-related header includes, for
  clang version 13 compatibility

## excursions 2.5.1

CRAN release: 2021-01-21

- Update links to INLA in documentation

## excursions 2.5.0

- Remove dependency on INLA for tests on CRAN
- Fix bug for diagonal matrices of size 1x1
- Add fields URL and BugReports to DESCRIPTION

## excursions 2.4.5

CRAN release: 2020-02-26

- Fix bug to handle empty sets in continuous interpretations
- Fix compiler warnings; unused variable and mismatching Qinv
  declaration
- Update COPYRIGHT information

## excursions 2.4.4

CRAN release: 2018-08-31

- Set fixed buffer size for RngStream.c stream name
- Update Makevars to avoid OpenMP warning

## excursions 2.4.3

- Improved INLA backwards compatibility

## excursions 2.4.2

- Add support for general manifolds in continuous interpretation methods
- Bug fixes to continuous interpretation methods
- Update CITATION information with new JSS manuscript

## excursions 2.4.1

CRAN release: 2018-07-19

- Minor fixed to C code to avoid warnings.

## excursions 2.4.0

- Added a `NEWS.md` file to track changes to the package.

- Updated repository for INLA

- Add n.iter and seed options for
  [`contourmap()`](https://davidbolin.github.io/excursions/reference/contourmap.md)
  calculations.

- Add support for the QC method for
  [`contourmap.inla()`](https://davidbolin.github.io/excursions/reference/contourmap.inla.md).

- Updated documentation.

## excursions 2.3.6

- Fix bug affecting contourmap calculations, caused by the previous fix
  to
  [`gaussint()`](https://davidbolin.github.io/excursions/reference/gaussint.md).

## excursions 2.3.5

- Fixed reordering bug in
  [`excursions.variances()`](https://davidbolin.github.io/excursions/reference/excursions.variances.md)
  caused by the previous fix to
  [`gaussint()`](https://davidbolin.github.io/excursions/reference/gaussint.md).

## excursions 2.3.4

- Fixed reordering bug in
  [`gaussint()`](https://davidbolin.github.io/excursions/reference/gaussint.md)
  that affected standalone use only.

## excursions 2.3.3

CRAN release: 2018-04-03

- Removed use of deprecated `rBind` and `cBind` and add dependency on
  `R >= 3.2.0`

## excursions 2.3.2

CRAN release: 2018-01-28

- Added citation information to DESCRIPTION

## excursions 2.3.1

- Add whitespace after `-f` for make, for wider compatibility, see
  [http://pubs/opegroup.org/onlinepubs/9699919799/utilities/make.html](http://pubs/opegroup.org/onlinepubs/9699919799/utilities/make.md)

- Updated CITATION information

## excursions 2.3.0

CRAN release: 2018-01-27
