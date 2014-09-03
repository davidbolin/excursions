\name{excursions.inla}
\alias{excursions.inla}
\title{
  Calculate excursion sets based on INLA output
}
\description{
  Calculate excursion function or contour uncertainty function for
  models calculated using the INLA method.
}
\usage{
excursions.inla(result.inla, ind=NULL, method, alpha=1, u, type,
                n.iter=10000, verbose=0, max.threads=0)
}
\arguments{
  \item{result.inla}{Result object from INLA call}
  \item{alpha}{Error probability for when to stop the calculation of the
    excursion function. The default value is 1. If only excursion sets
    with small error probabilities are of interest, \code{alpha} can be
    decreased to decrease computation time. }
  \item{u}{Excursion or contour level.}
 \item{type}{Type of region:
    \itemize{
    \item{'>' }{positive excursions}
    \item{'<' }{negative excursions}
    \item{'!=' }{contour avoiding function}
    \item{'=' }{contour credibility function}}}
   	\item{n.iter}{Number or iterations in the MC sampler that is
	  used for approximating probabilities. The default value is
	  10000.}
 \item{ind}{Indices of nodes in the joint distribution for which the
   excursion function should be calculated}
   \item{verbose}{Set to TRUE for verbose mode (optional)}
   \item{max.threads}{Decides the number of threads the program can use. Set to 0 for using the maximum number of threads allowed by the system (default).}
  \item{method}{Method for handeling the latent Gaussian structure:
  \itemize{
  	\item{'EB' }{Empirical Bayes}
  	\item{'QC' }{Quantile correction}
  	\item{'NI' }{Numerical integration}
  	\item{'NIQC' }{Numerical integration with quantile correction}
  	\item{'iNIQC' }{Improved integration with quantle correction}
  }
  }
}
\value{
  A list:
  \item{F }{Excursion function calculated for values up to alpha}
  \item{rho }{Marginal excursion probabilities}
  \item{mean }{Posterior mean}
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
Bolin, D. and Lindgren, F. (2013) \emph{Excursion and contour uncertainty regions for latent Gaussian models}, Journal of the Royal Statistical Society: Series B (in press).
}

\examples{
## In this example, we calculate the excursion function 
## for a partially observed AR process.
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
iprd <- inla.stack.index(stack, "pred")$data
formula = y ~ -1 + cov + f(ar,model="ar1")
result = inla(formula=formula, family="normal", 
		data = inla.stack.data(stack),
		control.predictor=list(A=inla.stack.A(stack), 
			compute=TRUE,cdf=c(0)),
		control.compute = list(config = TRUE))

## calculate the level 0 positive excursion function
res.qc =
excursions.inla(result, ind=iprd, alpha=0.99, u=0, method='QC', type='>',
                max.threads=2)

## plot the excursion function and marginal probabilities
plot(res.qc$rho,type="l",
     main="marginal probabilities (black) and excursion function (red)")
lines(res.qc$F,col=2)
}
}
