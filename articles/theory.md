# Definitions and computational methodology

Hierarchical models are of great importance in many areas of statistics.
In the simplest form, a hierarchical model has a likelihood distribution
$`\pi({\boldsymbol{\mathrm{Y}}}|{\boldsymbol{\mathrm{X}}}, {\boldsymbol{\mathrm{\theta}}})`$
for observed data $`{\boldsymbol{\mathrm{Y}}}`$, which is specified
conditionally on a latent process of interest,
$`{\boldsymbol{\mathrm{X}}}`$, which has a distribution
$`\pi({\boldsymbol{\mathrm{X}}}|{\boldsymbol{\mathrm{\theta}}})`$. For
Bayesian hierarchical models, one also specifies prior distributions for
the model parameters $`{\boldsymbol{\mathrm{\theta}}}`$. The most
important special case of these models are the LGMs, which are obtained
by assuming that
$`{\boldsymbol{\mathrm{X}}}|{\boldsymbol{\mathrm{\theta}}}`$ has a
Gaussian distribution (Rue et al. 2009). Numerous applications can be
studied using models of this form, and these are therefore the main
focus of the methods in `excursions`.

A statistical analysis using an LGM often concludes with reporting the
posterior mean
$`E({\boldsymbol{\mathrm{X}}}|{\boldsymbol{\mathrm{Y}}})`$ as a point
estimate of the latent field, possibly together with posterior variances
as a measure of uncertainty. In many applications, however, reporting
posterior means and variances are not enough. As stated in the
introduction, one may be interested in computing regions where the
latent field exceeds some given threshold, contour curves with their
associated uncertainty, or simultaneous confidence bands. In some
applications, only a contour map of the posterior mean is reported,
where the number of levels in the contour map should represent the
uncertainty in the estimate. These are quantities that can be computed
with `excursions` and we now define these in more detail before
outlining how they can be computed. For details we refer to (Bolin and
Lindgren 2015) and (Bolin and Lindgren 2017).

## Definitions

The main quantities that can be computed using `excursions` are (1)
excursion sets, (2) contour credible regions and level avoiding sets,
(3) excursion functions, (4) contour maps and their quality measures,
(5) simultaneous confidence bands. This section defines these in more
detail.

Throughout the section, $`X({\boldsymbol{\mathrm{s}}})`$ will denote a
stochastic process defined on some domain of interest, $`\Omega`$, which
we assume is open with a well-defined area $`|\Omega|<\infty`$. Since it
is not necessary for the definitions, we will not explicitly state the
dependency on the data, but simply allow $`X`$ to have some possible
non-stationary distribution. In practice, however, the distribution of
$`X`$ will typically be a posterior distribution conditionally on data,
$`X({\boldsymbol{\mathrm{s}}})|{\boldsymbol{\mathrm{Y}}}`$. For
frequentist models, the distribution of $`X`$ could also be
conditionally on for example a maximum likelihood estimate of the
parameters,
$`X({\boldsymbol{\mathrm{s}}})|{\boldsymbol{\mathrm{Y}}},\widehat{{\boldsymbol{\mathrm{\theta}}}}`$.

### Excursion sets

An excursion set is a set where the process
$`X({\boldsymbol{\mathrm{s}}})`$ exceeds or goes below some given level
of interest, $`u`$. A where $`X({\boldsymbol{\mathrm{s}}})>u`$ is
referred to as a positive excursion set, whereas a set where
$`X({\boldsymbol{\mathrm{s}}})<u`$ is referred to as a negative
excursion set. If
$`X({\boldsymbol{\mathrm{s}}})=f({\boldsymbol{\mathrm{s}}})`$ is a known
function, these sets can be computed directly as
$`A_u^+(f) = \{ {\boldsymbol{\mathrm{s}}}\in\Omega; f({\boldsymbol{\mathrm{s}}})>u \}`$
and
$`A_u^-(f) = \{ {\boldsymbol{\mathrm{s}}}\in\Omega; f({\boldsymbol{\mathrm{s}}})<u \}`$
respectively. If $`X({\boldsymbol{\mathrm{s}}})`$ is a latent random
process, one can only provide a region where it with some (high)
probability exceeds the level. More specifically, the positive level
$`u`$ excursion set with probability $`\alpha`$,
$`E_{{u,\alpha}}^{{+}}(X)`$, is defined as the largest set so that with
probability $`1-\alpha`$ the level $`u`$ is exceeded at all locations in
the set,
``` math
E_{{u,\alpha}}^{{+}} = \mathop{\mathrm{arg\,max}}_{D}\{|D|:P[D\subset A_u^+(X)]\geq 1 - \alpha\}.
```
Similarly, the negative $`u`$ excursion set with probability $`\alpha`$,
$`E_{{u,\alpha}}^{{-}}(X)`$, is defined as the largest set so that with
probability $`1-\alpha`$ the process is below the level $`u`$ at all
locations in the set. This set is obtained by replacing $`A_u^+(X)`$
with $`A_u^-(X)`$ in the equation above.

### Contour credible regions and level avoiding sets

For a function $`f`$, a contour curve (or set in general) of a level
$`u`$ is defined as the set of all level $`u`$ crossings. Formally, the
level curve is defined as
$`A_u^c(f) = \left(A_u^+(f)^o\cup A_u^-(f)^o\right)^c`$, where $`B^o`$
denotes the interior of the set $`B`$ and $`B^c`$ denotes the
complement. Note that $`A_u^c(f)`$ not only includes the set of
locations where $`f({\boldsymbol{\mathrm{s}}})=u`$, but also all
discontinuous level crossings.

For a latent random process $`X`$, one can only provide a credible
region for the contour curve. A level $`u`$ contour credibility region,
$`E_{{u,\alpha}}^{{c}}(X)`$, is defined as the smallest set such that
with probability $`1-\alpha`$ level $`u`$ crossings of $`X`$ are in the
set. This set can be seen as the complement of the level $`u`$ contour
avoiding set $`E_{{u,\alpha}}^{{}}(X)`$, which is defined as the largest
union $`M_{u,\alpha}^+ \cup M_{u,\alpha}^-`$, where jointly
$`X({\boldsymbol{\mathrm{s}}})>u`$ in $`M_{u,\alpha}^+`$ and
$`X({\boldsymbol{\mathrm{s}}})<u`$ in $`M_{u,\alpha}^-`$. Formally,
``` math
\begin{equation*}
(M_{u,\alpha}^+(X), M_{u,\alpha}^-(X)) = \mathop{\mathrm{arg\,max}}_{(D^+,D^-)}\{|D^-\cup D^+| :
P(D^- \subseteq A_u^-(X),\, D^+ \subseteq A_u^+(X))
\geq 1-\alpha\},
\end{equation*}
```
where the sets $`(D^+,D^-)`$ are open. The sets $`M_{u,\alpha}^+`$ and
$`M_{u,\alpha}^-`$ are denoted as the pair of level $`u`$ avoiding sets.
The notion of level avoiding sets can naturally be extended to multiple
levels $`u_1< u_2 < \cdots < u_k`$, which is needed when studying
contour maps. In this case, the multilevel contour avoiding set is
denoted $`C_{{\boldsymbol{\mathrm{u}}},\alpha}(X)`$ (For a formal
definition, see (Bolin and Lindgren 2017)).

### Excursion functions

introduced excursion functions as a tool for visualizing excursion sets
simultaneously for all values of $`\alpha`$. For a level $`u`$, the
positive and negative excursion functions are defined as
$`F_u^+({\boldsymbol{\mathrm{s}}}) = 1 - \inf\{\alpha ; {\boldsymbol{\mathrm{s}}}\in E_{{u,\alpha}}^{{+}} \}`$
and
$`F_u^-({\boldsymbol{\mathrm{s}}}) = 1 - \inf\{\alpha ; {\boldsymbol{\mathrm{s}}}\in E_{{u,\alpha}}^{{-}} \}`$,
respectively. Similarly, the contour avoidance function, and the contour
function are defined as
$`F_u({\boldsymbol{\mathrm{s}}}) = 1 -\inf\{\alpha ; {\boldsymbol{\mathrm{s}}}\in E_{{u,\alpha}}^{{}} \}`$
and
$`F_u^c({\boldsymbol{\mathrm{s}}}) = \sup\{\alpha; {\boldsymbol{\mathrm{s}}}\in E_{{u,\alpha}}^{{c}}\}`$,
respectively. Finally, for levels $`u_1< u_2 < \cdots < u_k`$, one can
define a contour map function as
$`F({\boldsymbol{\mathrm{s}}}) = \sup\{1-\alpha; {\boldsymbol{\mathrm{s}}}\in C_{u,\alpha}\}`$.

These functions take values between zero and one and each set
$`E_{{u,\alpha}}^{{\bullet}}`$ can be retrieved as the $`1-\alpha`$
excursion set of the function
$`F_u^{\bullet}({\boldsymbol{\mathrm{s}}})`$.

### Contour maps and their quality measures

For a function $`f({\boldsymbol{\mathrm{s}}})`$, a contour map $`C_f`$
with contour levels $`u_1 < u_2 < \ldots < u_K`$ is defined as the
collection of contour curves $`A_{u_1}^c(f), \ldots, A_{u_K}^c(f)`$ and
associated level sets
$`G_k = \{{\boldsymbol{\mathrm{s}}} : u_{k}< f({\boldsymbol{\mathrm{s}}}) < u_{k+1}\}`$,
for $`0\leq k\leq K`$, where one defines $`u_0 = -\infty`$ and
$`u_{K+1}=\infty`$. In practice, a process
$`X({\boldsymbol{\mathrm{s}}})`$ is often visualized using a contour map
of the posterior mean
$`E(X({\boldsymbol{\mathrm{s}}})|{\boldsymbol{\mathrm{Y}}})`$. The
contour maps is visualized either by just drawing the contour curves
labeled by their values, or by also visualizing each level set in a
specific color. The color for a set $`G_k`$ is typically chosen as the
color that corresponds to the level $`u_k^e = (u_k+u_{k+1})/2`$ in a
given color map.

In order to choose an appropriate number of contours, one must be able
to quantify the uncertainty of contour maps. The uncertainty can be
represented using a contour map quality measure $`P`$, which is a
function that takes values in $`[0, 1]`$. Here, $`P`$ should be chosen
in such a way that $`P \approx 1`$ indicates that the contour map, in
some sense, is appropriate as a description of the distribution of the
random field, whereas $`P \approx 0`$ should indicate that the contour
map is inappropriate.

An example of a contour map quality measure is the normalized integral
of the contour map function
``` math
\begin{equation}
P_0(X, C_f) = \frac1{|\Omega|}\int_{\Omega}F({\boldsymbol{\mathrm{s}}})d{\boldsymbol{\mathrm{s}}}.
\end{equation}
```
The most useful quality measure is denoted $`P_2`$ and is defined as the
simultaneous probability for the level crossings of
$`(u_1^e, \ldots , u_K^e)`$ all falling within their respective level
sets $`(G_1, \ldots, G_K)`$ .

An intuitively interpretable approach for choosing the number of
contours in a contour map is to find the largest K such that $`P_2`$ is
above some threshold. For a joint credibility of $`90\%`$, say, choose
the largest number of contours such that $`P_2 \geq 0.9`$. How this can
be done using `excursions` is illustrated in Section~.

### Simultaneous confidence bands

Especially for time series applications, the uncertainty in the latent
process is often visualized using pointwise confidence bands. A
pointwise confidence interval for $`X`$ at some location
$`{\boldsymbol{\mathrm{s}}}`$ is given by
$`[q_{\alpha/2}({\boldsymbol{\mathrm{s}}}),q_{1-\alpha/2}({\boldsymbol{\mathrm{s}}})]`$,
where $`q_{\alpha}({\boldsymbol{\mathrm{s}}})`$ denotes the
$`\alpha`$-quantile in the marginal distribution of
$`X({\boldsymbol{\mathrm{s}}})`$.

A problem with using pointwise confidence bands is that there is not
joint interpretation, and one is therefore often of interested in
computing simultaneous confidence bands. For a process
$`X({\boldsymbol{\mathrm{s}}}), {\boldsymbol{\mathrm{s}}}\in \Omega`$,
we define a simultaneous confidence band as the region
$`\{({\boldsymbol{\mathrm{s}}},y): {\boldsymbol{\mathrm{s}}}\in \Omega, q_{\rho}({\boldsymbol{\mathrm{s}}}) \leq y \leq q_{1-\rho}({\boldsymbol{\mathrm{s}}})\}`$.
Here $`\rho`$ is chosen such that
$`P(q_{\rho}({\boldsymbol{\mathrm{s}}})<X({\boldsymbol{\mathrm{s}}}) <q_{1-\rho}({\boldsymbol{\mathrm{s}}}), {\boldsymbol{\mathrm{s}}}\in \Omega)=1-\alpha`$.
Thus $`\alpha`$ controls the probability that the process is inside the
confidence band at all locations in $`\Omega`$.

## Computational methods

If the latent process $`X({\boldsymbol{\mathrm{s}}})`$ is defined on a
continuous domain $`\Omega`$, one has to use a discretization,
$`{\boldsymbol{\mathrm{x}}}`$, of the process in the statistical
analysis. The vector $`{\boldsymbol{\mathrm{x}}}`$ may be the process
evaluated at some locations of interest or more generally the weights in
a basis expansion
$`X({\boldsymbol{\mathrm{s}}}) = \sum_i \varphi_i({\boldsymbol{\mathrm{s}}}) x_i`$.
Computations in the package are in a first stage performed using the
distribution of $`{\boldsymbol{\mathrm{x}}}`$. If $`\Omega`$ is
continuous, computations in a second stage interpolates the results for
$`{\boldsymbol{\mathrm{x}}}`$ to $`\Omega`$. In this section, we briefly
outline the computational methods used in these two stages. As most
quantities of interest can be obtained using some excursion function, we
focus on how these are computed. As a representative example, we outline
how $`F_u^+({\boldsymbol{\mathrm{s}}})`$ is computed in the following
sections.

As before, let $`{\boldsymbol{\mathrm{Y}}}`$ and
$`{\boldsymbol{\mathrm{\theta}}}`$ be a vectors respectively containing
observations and model parameters. Computing an excursion function
$`{\boldsymbol{\mathrm{F}}}_u^+ = \{F_u^+(x_1), \ldots, F_u^+(x_n)\}`$
requires computing integrals of the posterior distribution for
$`{\boldsymbol{\mathrm{x}}}`$. To save computation time, it is assumed
that $`E_{{u,\alpha_1}}^{{+}} \subset E_{{u,\alpha_2}}^{{+}}`$ if
$`\alpha_1 > \alpha_2`$. This means that $`{\boldsymbol{\mathrm{F}}}_u`$
can be obtained by first reordering the nodes and then computing a
sequential integral. The reordering is in this case obtained by sorting
the marginal probabilities $`P(x_i > u)`$ (for other options, see (Bolin
and Lindgren 2015). After reordering, the $`i`$:th element of
$`{\boldsymbol{\mathrm{F}}}_u^+`$ is obtained as the integral
``` math
\int_u^{\infty}\pi({\boldsymbol{\mathrm{x}}}_{1:i}|{\boldsymbol{\mathrm{Y}}}) d{\boldsymbol{\mathrm{x}}}_{1:i}.
```
Using sequential importance sampling as described below, the whole
sequence of integrals can be obtained with the same cost as computing
only one integral with $`i=n`$, making the computation of
$`{\boldsymbol{\mathrm{F}}}_u^+`$ feasible also for large problems.

### Gaussian integrals

The basis for the computational methods in the package is the ability to
compute the required integral when the posterior distribution is
Gaussian. In this case, one should compute an integral
``` math
\begin{equation}\label{eq:markovint}
\frac{|{\boldsymbol{\mathrm{Q}}}|^{1/2}}{(2\pi)^{d/2}}\int_{{\boldsymbol{\mathrm{u}}}-{\boldsymbol{\mathrm{\mu}}}\leq {\boldsymbol{\mathrm{x}}}}\exp\left(-\frac{1}{2}{\boldsymbol{\mathrm{x}}}^{\top}{\boldsymbol{\mathrm{Q}}} {\boldsymbol{\mathrm{x}}}\right) \,\mathrm{d}{\boldsymbol{\mathrm{x}}}.
\end{equation}
```
Here $`{\boldsymbol{\mathrm{\mu}}}`$ and $`{\boldsymbol{\mathrm{Q}}}`$
are the posterior mean and posterior precision matrix respectively. To
take advantage of the possible sparsity of $`{\boldsymbol{\mathrm{Q}}}`$
if a Markovian model is used, the integral is rewritten as
``` math
\begin{equation}\label{eq:seqint}
\int_{a_d}^{\infty}\pi(x_d)
\int_{a_{d-1}}^{\infty}\pi(x_{d-1}|x_d)
\cdots
\int_{a_2}^{\infty}\pi(x_2|{\boldsymbol{\mathrm{x}}}_{3:d})
\int_{a_1}^{\infty} \pi(x_1|{\boldsymbol{\mathrm{x}}}_{2:d})
\,\mathrm{d}{\boldsymbol{\mathrm{x}}}
\end{equation}
```
where, if the problem has a Markov structure,
$`x_i|{\boldsymbol{\mathrm{x}}}_{i+1:d}`$ only depends on a few of the
elements in $`x_{i+1:d}`$ given by the Markov structure. The integral is
calculated using a sequential importance sampler by starting with the
integral of the last component $`\pi(x_d)`$ and then moving backward in
the indices, see for further details.

### Handling non-Gaussian data

Using the sequential importance sampler above,
$`{\boldsymbol{\mathrm{F}}}_u^+`$ can be computed for Gaussian models
with known parameters. For more complicated models, the latent Gaussian
structure has to be handled, and this can be done in different ways
depending on the accuracy that is needed. `excursions` currently
supports the following five methods described in detail in (Bolin and
Lindgren 2015): Empirical Bayes (`EB`), Quantile Correction (`QC`),
Numerical Integration (`NI`), Numerical integration with Quantile
Corrections (`NIQC`), and improved Numerical integration with Quantile
Corrections (`iNIQC`).

The `EB` method is the simplest method and is based on using a Gaussian
approximation of the posterior,
$`\pi({\boldsymbol{\mathrm{x}}}|{\boldsymbol{\mathrm{Y}}}) \approx \pi_G({\boldsymbol{\mathrm{x}}}|{\boldsymbol{\mathrm{Y}}},\widehat{{\boldsymbol{\mathrm{\theta}}}})`$.
The `QC` method uses the same Gaussian approximation but modifies the
limits in the integral to improve the accuracy. The three other methods
are intended for Bayesian models, where the posterior is obtained by
integrating over the parameters,
``` math
\begin{equation*}
\pi({\boldsymbol{\mathrm{x}}}\mid{\boldsymbol{\mathrm{Y}}}) = \int \pi({\boldsymbol{\mathrm{x}}}\mid{\boldsymbol{\mathrm{Y}}},{\boldsymbol{\mathrm{\theta}}})\pi({\boldsymbol{\mathrm{\theta}}}\mid{\boldsymbol{\mathrm{Y}}})d{\boldsymbol{\mathrm{\theta}}}.
\end{equation*}
```
The `NI` method approximates the integration with respect to the
parameters as in the INLA method, using a sum of representative
parameter configurations, and the `NIQC` and `iNIQC` methods combines
this with the `QC` method to improve the accuracy further.

In general, `EB` and `QC` are suitable for frequentist models and for
Bayesian models where the posterior distribution conditionally on the
parameters is approximately Gaussian. The methods are equivalent if the
posterior is Gaussian and in other cases `QC` is more accurate than
`EB`. For Bayesian models, the `NI` method is, in general, more accurate
than the `QC` method, and for non-Gaussian likelihoods, the `NIQC` and
`iNIQC` methods can be used for improved results. In general the
accuracy and computational cost of the methods are as follows

Accuracy: `EB` $`<`$`QC` $`<`$`NI` $`<`$`NIQC` $`<`$`iNIQC`

Computational cost: `EB` $`\approx`$`QC` $`<`$`NI` $`\approx`$`NIQC`
$`<`$`iNIQC`.

If the main purpose of the analysis is to construct excursion or contour
sets for low values of $`\alpha`$, we recommend using the `QC` method
for problems with Gaussian likelihoods and the `NIQC` method for
problems with non-Gaussian likelihoods. The increase in accuracy of the
`iNIQC` method is often small compared to the added computational cost.

### Continuous domain interpolations

For a continuous spatial domain, the excursion function
$`F_u^+({\boldsymbol{\mathrm{s}}})`$ can be approximated using
$`{\boldsymbol{\mathrm{F}}}_u^+`$ computed at discrete point locations.
The main idea is to interpolate $`{\boldsymbol{\mathrm{F}}}_u^+`$
assuming monotonicity of the random field between the discrete
computation locations. Specifically, assume that the values of
$`{\boldsymbol{\mathrm{F}}}_u^+`$ correspond to the values at the
vertices of some triangulated mesh. If the excursion set
$`E_{{u,\alpha}}^{{+}}(X)`$ should be computed for some specific value
of $`\alpha`$, one has to find the $`1-\alpha`$ contour for
$`F_u^+({\boldsymbol{\mathrm{s}}})`$. For interpolation to a location
$`{\boldsymbol{\mathrm{s}}}`$ within a specific triangle $`\mathcal{T}`$
with corners in
$`{\boldsymbol{\mathrm{s}}}_1, {\boldsymbol{\mathrm{s}}}_2,`$ and
$`{\boldsymbol{\mathrm{s}}}_3`$, `excursions` by default uses log-linear
interpolation,
$`F_u({\boldsymbol{\mathrm{s}}}) = \exp\{\sum_{k=1}^3 w_k\log[F_u({\boldsymbol{\mathrm{s}}}_k)]\}`$.
Here $`\{(w_1,w_2,w_3);\, w_1,w_2,w_3\geq 0,\, \sum_{k=1}^3 w_i=1\}`$
are the barycentric coordinates of $`{\boldsymbol{\mathrm{s}}}`$ within
the triangle.

Further technical details of the continuous domain construction are
given in (Bolin and Lindgren 2017). Studies of the resulting continuous
domain excursion sets in (Bolin and Lindgren 2017) indicate that the
log-linear interpolation method results in sets with coverage
probability on target or slightly above target for large target
probabilities.

## References

Bolin, David, and Finn Lindgren. 2015. “Excursion and Contour
Uncertainty Regions for Latent Gaussian Models.” *Journal of the Royal
Statistical Society B* 77: 85–106. <https://doi.org/10.1111/rssb.12055>.

Bolin, David, and Finn Lindgren. 2017. “Quantifying the Uncertainty of
Contour Maps.” *Journal of Computational and Graphical Statistics* 26
(3): 513–24. <https://doi.org/10.1080/10618600.2016.1228537>.

Rue, H, S Martino, and N Chopin. 2009. “Approximate Bayesian Inference
for Latent Gaussian Models Using Integrated Nested Laplace
Approximations (with Discussion).” *Journal of the Royal Statistical
Society B* 71 (2): 319–92.
<https://doi.org/10.1111/j.1467-9868.2008.00700.x>.
