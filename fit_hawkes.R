library(hawkes)
library(data.table)
library(ggplot2)
library(stats) # for constrOptim
# https://cran.r-project.org/web/packages/hawkes/hawkes.pdf

if(dir.exists('/home/andrea/'))
{
        git.dir <- '~/git/modelling_events'
}else{
        git.dir <- ''
}
# TODO: make your own git.dir

outdir <- file.path(git.dir, 'results')
outdir_1D <- file.path(outdir, 'hawkes_1D')


# source(file.path(git.dir, 'visualisations.R'))
source(file.path(git.dir, 'simulate_hawkes.R'))
source(file.path(git.dir, 'hawkes_helpers.R'))


# Hawkes library only allows for exponential excitation functions
# meaning h(t) = alpha exp(-beta t) -> int h(t) = alpha/beta 
# so:
# E[ children caused ] = alpha/beta
# E[ time to children ] = 1/beta

LAMBDA0 = 2
ALPHA = .6
BETA = .81
HORIZON = 100

recreate_allal <- function(lambda0=LAMBDA0, alpha=ALPHA, beta=BETA, horizon=HORIZON, return.plot=TRUE)
{

        h_expdecay <- function(x, alpha, beta)
                alpha * exp(- beta * x )

        compensator <- function(x, pars, history)
        {
                out <- rep(NA, length(x))

                for(idx in seq_along(x))
                {
                        y <- x[idx]
                        
                        events_prior_y <- history[history <= y]
                        out[idx] <- pars$lambda0 + sum(h_expdecay( y - events_prior_y, pars$alpha, pars$beta))
                }
                out
        }

        integrated_compensator <- function(x, pars,  history)
        {
                # check entries are sorted
                stopifnot( all(sort(x, decreasing=FALSE) == x ))

                integrated_h <- function(x)
                        with(pars, alpha/beta * (1 - exp(-beta * x)))

                # intensities associated to events that happened
                out <- rep(NA, length(x))
                for(idx in seq_along(x))
                {
                        y <- x[idx]
                        base = pars$lambda0 * y
                        events_prior_y <- history[history <= y]
                        out[idx] <- base + sum(integrated_h(y - events_prior_y))
                }
                out
        }

        print_parameters <- function(x)
        {
                for (n in names(x))
                        cat(n, ': ', x[[n]], '\n')
                cat('\n\n')
        }


        # simulate data: one dimensional Hawkes process
        h_1 <- simulateHawkes(lambda0=lambda0,
                            alpha=alpha,
                            beta=beta,
                            horizon=horizon)

        # estimate parameters via ML
        hat_pars <- maximum_likelihood(h_1[[1]])
        print_parameters(hat_pars)
        
        # find the compensator
        pois_history <- integrated_compensator(h_1[[1]], pars=hat_pars, history=h_1[[1]])
        pois_difftimes <- diff(pois_history)
                
        # Make plots and tests
        if(return.plot)
        {
                p1 <- plot.trajectory.and.conditionalintensity(history=h_1[[1]])
                force(p1)
        }

        envelope <- plot.envelope.test(poisson_t = pois_history)
        p2 <- envelope$plot
        if(return.plot)
                force(p2)

        kstest <- make.cdfplot.kstest(poisson_deltat=pois_difftimes)
        if(return.plot)
        {
                p3 <- kstest$plot
                force(p3)
        }

        if(return.plot)
        {
                p4 <- plot.martingale.residuals(pois_history, history=h_1[[1]])
                force(p4)
        }

        out <- list(KStest = kstest$pvalue, envelopeout=envelope$test)
        if(return.plot)
        {
                p <- ggpubr::ggarrange(p1, p2, p3, p4)
                out1 <- list( plot = p)

                out <- append(out, out1)
        }
        return(out)
}


out <- recreate_allal()
filename=file.path(outdir_1D, 'allal_4diagnostics_example.png')
ggsave(out$plot, filename=filename, width=12, height=10, units='cm')



if(0)
{
#Multivariate Hawkes process
        lambda0<-c(0.2,0.2)
        alpha<-matrix(c(0.5,0,0,0.5),byrow=TRUE,nrow=2)
        beta<-c(0.7,0.7)
        horizon<-3600#one hour
        h<-simulateHawkes(lambda0,alpha,beta,horizon)
}

