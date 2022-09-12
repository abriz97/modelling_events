library(survival)
library(data.table)
library(ggfortify)

# baseline hazard: Weibull
# h_0(t)=lambda rho t^(rho-1)

# N = sample size    
# lambda = scale parameter in h0()
# rho = shape parameter in h0()
# beta = fixed effect parameter
# rateC = rate parameter of the exponential distribution of C

# Simulate data
# _____________

simulWeib <- function(N, lambda, rho, beta, rateC)
{
  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  # Weibull latent event times
  u <- runif(n=N)
  # Inverse cdf
  T <- (- log(u) / (lambda * exp(x * beta)))^(1 / rho)

  # censoring times
  C <- rexp(n=N, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(T, C)
  delta <-  as.numeric(T <= C)

  # data set
  data.table(id=1:N,
             time=time,
             delta=delta,
             x=x)
}

simulWeib_diff_cens <- function(N, lambda, rho, beta, rateC)
{
  # covariate --> N Bernoulli trials
  x <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))

  # Weibull latent event times
  u <- runif(n=N)
  # Inverse cdf
  T <- (- log(u) / (lambda * exp(x * beta)))^(1 / rho)

  # censoring times:
  C_0 <- rexp(n=N, rate=rateC)
  C_1 <- rexp(n=N, rate= 1/2 * rateC)

  if(rho == 2)
  {
          cat(mean(T), ' - ')
          cat(1/(4*lambda) * sqrt(2*pi), '\n')
          thres <- 1/2*sqrt(lambda/pi)
  }else{
          thres <- mean(T)
  }

  # select which participants have long Survival times
  idx <- T > thres

  C <- C_0
  C[idx] <- C_1[idx]

  # follow-up times and event indicators
  time <- pmin(T, C)
  delta <-  as.numeric(T <= C)

  # data set
  data.table(id=1:N,
             time=time,
             delta=delta,
             x=x)
}

weibull_cdf <- function(t, rho, lambda)
        1 - exp( - lambda * t^rho)


# Test
# ____

N=100
BETA=-.6
LAMBDA=.01
RHO=2
RATEC=.001


for(k in 1:1e3)
{
        dat <- simulWeib(N=N, lambda=LAMBDA, rho=RHO, beta=BETA, rateC=RATEC)

        # Representation of survival data: + means failure time is larger than X
        
        if(0)
        {
                with(dat,
                        Surv(time,delta)
                ) -> tmp

                tmp1 <- data.table(
                        range = 1:700,
                        survival0 = 1- weibull_cdf(1:700, rho=1, lambda=.01 * exp(0) ),
                        survival1 = 1- weibull_cdf(1:700, rho=1, lambda=.01 * exp(-0.6) )
                )
                tmp2 <- tmp1[, .(range, survival1)]

                idx0 <- dat$x == 0
                idx1 <- dat$x == 1

                plot(tmp)
                lines(tmp1, col='blue')
                lines(tmp2, col='red')

                plot(tmp[idx0])
                lines(tmp1, col='blue')

                plot(tmp[idx1])
                lines(tmp2, col='red')
        }

        fit <- coxph(Surv(time, delta) ~ x, data=dat)

        betaHat[k] <- fit$coef
}

cat(
    "Mean Estimated beta:", mean(betaHat), '\n',
    "True Beta:", BETA, '\n'
)

hist(betaHat)
abline(v=BETA, col='red')
abline(v=0, col='blue')


# Dependent censoring

for(k in 1:1e3)
{
        dat <- simulWeib_diff_cens(N=N, lambda=LAMBDA, rho=RHO, beta=BETA, rateC=RATEC)

        # Representation of survival data: + means failure time is larger than X
        
        if(0)
        {
                with(dat,
                        Surv(time,delta)
                ) -> tmp

                tmp1 <- data.table(
                        range = 1:700,
                        survival0 = 1- weibull_cdf(1:700, rho=1, lambda=.01 * exp(0) ),
                        survival1 = 1- weibull_cdf(1:700, rho=1, lambda=.01 * exp(-0.6) )
                )
                tmp2 <- tmp1[, .(range, survival1)]

                idx0 <- dat$x == 0
                idx1 <- dat$x == 1

                plot(tmp)
                lines(tmp1, col='blue')
                lines(tmp2, col='red')

                plot(tmp[idx0])
                lines(tmp1, col='blue')

                plot(tmp[idx1])
                lines(tmp2, col='red')
        }

        fit <- coxph(Surv(time, delta) ~ x, data=dat)

        betaHat[k] <- fit$coef
}

cat(
    "Mean Estimated beta:", mean(betaHat), '\n',
    "True Beta:", BETA, '\n'
)

hist(betaHat)
abline(v=BETA, col='red')
abline(v=0, col='blue')



