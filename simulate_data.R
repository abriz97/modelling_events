library(survival)
library(data.table)

# baseline hazard: Weibull
# h_0(t)=lambda rho t^(rho-1)

# N = sample size    
# lambda = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C

simulWeib <- function(N, lambda, rho, beta, rateC)
{
  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  # Weibull latent event times
  u <- runif(n=N)
  # Inverse cdf
  Tlat <- (- log(u) / (lambda * exp(x * beta)))^(1 / rho)

  # censoring times
  C <- rexp(n=N, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # data set
  data.frame(id=1:N,
             time=time,
             status=status,
             x=x)
}
