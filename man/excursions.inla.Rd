\name{excursions.inla}
\alias{excursions.inla}
\title{
Excursion sets and contour credible regions for latent Gaussian
models
}
\description{
  Excursion sets and contour credible regions for latent Gaussian
models calculated using the INLA method.
}
\usage{
excursions.inla(result.inla,
                stack,
                name=NULL,
                tag=NULL,
                ind=NULL,
                method,
                alpha=1,
                F.limit,
                u,
                u.link = FALSE,
                type,
                n.iter=10000,
                verbose=0,
                max.threads=0,
                seed=NULL)
}
\arguments{
  \item{result.inla }{Result object from INLA call}
  \item{stack }{The stack object used in the INLA call.}
  \item{tag }{The tag of the component in the stack for which to do the calculation. This argument should only be used if a stack object is provided, use the name argument otherwise.}
  \item{name }{The name of the component for which to do the calculation. This argument should only be used if a stack object is not provided, use the tag argument otherwise.}
  \item{ind }{If only a part of a component should be used in the calculations, this argument specifies the indices for that part.}
  \item{method }{Method for handeling the latent Gaussian structure:
    \itemize{
      \item{'EB' }{Empirical Bayes}
      \item{'QC' }{Quantile correction}
      \item{'NI' }{Numerical integration}
      \item{'NIQC' }{Numerical integration with quantile correction}
      \item{'iNIQC' }{Improved integration with quantle correction}
    }
  }
  \item{alpha }{Error probability for the excursion set of interest. The default
   value is 1.}
  \item{F.limit }{Error probability for when to stop the calculation of the
    excursion function. The default value is \code{alpha}, and the value cannot
    be smaller than \code{alpha}. A smaller value of \code{F.limit} results in a
    smaller compontation time.}
  \item{u }{Excursion or contour level.}
  \item{u.link }{If u.link is TRUE, \code{u} is assumed to be in the scale of the
                data and is then transformed to the scale of the linear
                predictor (default FALSE)}
 \item{type}{Type of region:
    \itemize{
    \item{'>' }{positive excursions}
    \item{'<' }{negative excursions}
    \item{'!=' }{contour avoiding function}
    \item{'=' }{contour credibility function}}}
   	\item{n.iter}{Number or iterations in the MC sampler that is
	  used for approximating probabilities. The default value is
	  10000.}
   \item{verbose }{Set to TRUE for verbose mode (optional)}
   \item{max.threads }{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
   \item{seed }{Random seed (optional).}
}
\value{

\code{excursions.inla} returns an object of class "excurobj". This is a list that contains the following arguments:
  \item{E }{Excursion set, contour credible region, or contour avoiding
    set}
  \item{F }{The excursion function corresponding to the set \code{E} calculated
    for values up to \code{F.limit}}
  \item{G }{ Contour map set. \eqn{G=1} for all nodes where the \eqn{mu > u}.}
  \item{M }{ Contour avoiding set. \eqn{M=-1} for all non-significant nodes. \eqn{M=0} for nodes where the process is significantly below \code{u} and \eqn{M=1} for all nodes where the field is significantly above \code{u}. Which values that should be present depends on what type of set that is calculated.}
  \item{rho }{Marginal excursion probabilities}
  \item{mean }{Posterior mean}
  \item{vars }{Marginal variances}
  \item{meta }{A list containing various information about the calculation.}
}
\note{
  This function requires the \code{INLA} package, which is not a CRAN
  package.  See \url{http://www.r-inla.org/download} for easy
  installation instructions.
}
\author{
  David Bolin \email{davidbolin@gmail.com} and
  Finn Lindgren \email{finn.lindgren@gmail.com}
}

\references{
Bolin, D. and Lindgren, F. (2015) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, JRSS-series B, vol 77, no 1, pp 85-106.
}

\examples{
## In this example, we calculate the excursion function
## for a partially observed AR process.
\donttest{
if (require(INLA)) {
## Sample the process:
rho = 0.9
tau = 15
tau.e = 1
n = 100
x = 1:n
mu = 10*((x<n/2)*(x-n/2) + (x>=n/2)*(n/2-x)+n/4)/n
Q = tau*sparseMatrix(i=c(1:n, 2:n), j=c(1:n, 1:(n-1)),
                   x=c(1,rep(1+rho^2, n-2),1, rep(-rho, n-1)),
                   dims=c(n, n), symmetric=TRUE)
X = mu+solve(chol(Q), rnorm(n))

## measure the sampled process at n.obs random locations
## under Gaussian measurement noise.
n.obs = 50
obs.loc = sample(1:n,n.obs)
A = sparseMatrix(i=1:n.obs, j=obs.loc, x=rep(1, n.obs), dims=c(n.obs, n))
Y = as.vector(A \%*\% X + rnorm(n.obs)/sqrt(tau.e))

## Estimate the parameters using INLA
ef = list(c(list(ar=x),list(cov=mu)))
s.obs = inla.stack(data=list(y=Y), A=list(A), effects=ef, tag="obs")
s.pre = inla.stack(data=list(y=NA), A=list(1), effects=ef,tag="pred")
stack = inla.stack(s.obs,s.pre)
formula = y ~ -1 + cov + f(ar,model="ar1")
result = inla(formula=formula, family="normal", data = inla.stack.data(stack),
              control.predictor=list(A=inla.stack.A(stack),compute=TRUE),
              control.compute = list(config = TRUE))

## calculate the level 0 positive excursion function
res.qc =
excursions.inla(result, stack = stack, tag = 'pred', alpha=0.99, u=0,
                method='QC', type='>', max.threads=2)

## plot the excursion function and marginal probabilities
plot(res.qc$rho,type="l",
     main="marginal probabilities (black) and excursion function (red)")
lines(res.qc$F,col=2)
}
}
}