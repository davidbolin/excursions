
integration.testdata1 <- function() {
  n = 11
  Q = Matrix(toeplitz(c(1, -0.1, rep(0, n-2))))
  mu = seq(-5, 5, length=n)
  L = chol(Q)
  seed = c(1750459768, 1598840523, 249795150, 
           2124039812, 1116456203, 1714779982)
  a = rep(-3,n) 
  b = rep(3,n)
  list(n=n,Q=Q,mu=mu,seed=seed,L=L,a=a,b=b)
}       


    