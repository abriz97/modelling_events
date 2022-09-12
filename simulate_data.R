library(survival)

#set parameters
set.seed(1234)

n = 40000 #sample size


#functional relationship

lambda=0.000020 #constant baseline hazard 2 per 100000 per 1 unit time

b_haz <-function(t) #baseline hazard
  {
    lambda #constant hazard wrt time 
  }

x = cbind(hba1c=rnorm(n,2,.5)-2,age=rnorm(n,40,5)-40,duration=rnorm(n,10,2)-10)

B = c(1.1,1.2,1.3) # hazard ratios (model coefficients)

hist(x %*% B) #distribution of scores

haz <-function(t) #hazard function
  b_haz(t) * exp(x %*% B)

c_hf <-function(t) #cumulative hazards function
  exp(x %*% B) * lambda * t 

S <- function(t) #survival function
{
  exp(-c_hf(t))
}

S(.005)
S(1)
S(5)

#simulate censoring

time = rnorm(n,10,2)

S_prob = S(time)

#simulate events

event = ifelse(runif(1)>S_prob,1,0)

#model fit

km = survfit(Surv(time,event)~1,data=data.frame(x))

plot(km) #kaplan-meier plot

#Cox PH model

fit = coxph(Surv(time,event)~ hba1c+age+duration, data=data.frame(x))

summary(fit)            

cox.zph(fit)
