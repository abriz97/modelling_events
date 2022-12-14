library(purrr)
library(hawkes)
library(data.table)
library(ggplot2)
library(stats) # for constrOptim
library(ggpubr)
# https://cran.r-project.org/web/packages/hawkes/hawkes.pdf


if(dir.exists('/home/andrea/'))
{
        git.dir <- '~/git/modelling_events'
}else{
        git.dir <- '~/Documents/CDT Modules/modelling_events'
}
# TODO: make your own git.dir

outdir <- file.path(git.dir, 'results')
outdir_1D <- file.path(outdir, 'hawkes_1D')
outdir_multi <- file.path(outdir, 'hawkes_multid')


# source(file.path(git.dir, 'visualisations.R'))
source(file.path(git.dir, 'simulate_hawkes.R'))
source(file.path(git.dir, 'functions', 'hawkes_helpers.R'))

# Hawkes library only allows for exponential excitation functions
# meaning h(t) = alpha exp(-beta t) -> int h(t) = alpha/beta 
# so:
# E[ children caused ] = alpha/beta
# E[ time to children ] = 1/beta

LAMBDA0 = 2
ALPHA = .6
BETA = .81
HORIZON = 100


# out <- recreate_allal()
# filename=file.path(outdir_1D, 'allal_4diagnostics_example.png')
# ggsave(out$plot, filename=filename, width=12, height=10, units='cm')

LAMBDA0 <- c(0.2, 0.3)
ALPHA <- matrix(c(.5, .3, 0, .1), byrow=T, 2)
BETA <- matrix(c(.6, .6, .6, .6), byrow=T, 2)
HORIZON <- 3600


recreate_multidimensional_allal <- function( lambda0=LAMBDA0, alpha=ALPHA, beta=BETA, horizon=HORIZON)
{
        h_expdecay <- function(x, pars, i,j)
                with(pars, alpha[i,j] * ( exp(-beta[i,j] * x)))

        compensator <- function(x, history, pars)
        {

                # get dimension of process
                D <- 1
                if(is.list(history))
                        D <- length(history)

                # Initialise output
                .init.with.baseline <- function(i)
                        rep(pars$lambda0[[i]], length(x))
                out <- lapply(1:D, .init.with.baseline)

                        
                # for each x0, find the events that happened, 
                # and then find the h-contributions to the rate by dimension
                for (idx in seq_along(x))
                {
                        x0 <- x[idx]
                        .select <- function(vec)
                                vec[vec <= x0]

                        events_prior_x0 <- lapply(history, .select)
                        diff_x0_epx0 <- lapply(events_prior_x0, function(a){x0 - a})

                        for(dim in 1:D)
                        {
                                .add.h.contribution.from.dim.i <- function(dim_i)
                                        sum(h_expdecay( diff_x0_epx0[[dim_i]], pars=hat_pars, j=dim, i=dim_i ))

                                tmp <- lapply(1:D, .add.h.contribution.from.dim.i )
                                out[[dim]][idx] <- purrr::reduce(tmp, `+`)
                        }
                }
                        
                return(out)
        }

        integrated_compensator <- function(pars,  history)
        {
                D <- length(history)

                integrated_h <- function(x, i,j)
                        with(pars, alpha[i,j]/beta[i,j] * (1 - exp(-beta[i,j] * x)))

                get.integrated.compensator <- function(dim)
                {

                        # intensities associated to events that happened
                        d <- history[[dim]]
                        out <- rep(NA, length(d))

                        for(idx in seq_along(d))
                        {
                                y <- d[idx]

                                # baseline intensity
                                base = pars$lambda0[[dim]] * y

                                .select <- function(vec)
                                        vec[vec <= y]
                                events_prior_y <- lapply(history, .select)

                                .f <- function(dim2)
                                {
                                        epy <- events_prior_y[[dim2]] 
                                        integrated_h( y - epy, i=dim2, j=dim)
        
                                }
                                exc_contributions <- unlist(lapply( 1:D, .f))
                        
                                out[idx] <- base + sum(exc_contributions)
                        }
                        out
                }

                out <- lapply(1:length(history), get.integrated.compensator)
        }

        print_parameters <- function(x)
        {
                for (n in names(x))
                        cat(n, ': ', round(x[[n]],2), '\n')
                cat('\n\n')
        }


        # simulate data: one dimensional Hawkes process
        true_pars <- list(
                lambda0 = lambda0,
                alpha = alpha,
                beta = beta,
                horizon = horizon)

        D <- max(
                length(lambda0),
                length(beta),
                sqrt(length(alpha))
        )
        
        
        h_1 <- with(true_pars,
                    simulateHawkes(lambda0=lambda0,
                                   alpha=alpha,
                                   beta=beta,
                                   horizon=horizon)
                    )
        plot_multi_hawkes(h_1, max=100)

        # estimate parameters via ML
        hat_pars <- maximum_likelihood_hawkes(h_1)

        print_parameters(true_pars)
        print_parameters(hat_pars)
        
        # find the compensator
        pois_history <- integrated_compensator( pars=hat_pars, history=h_1)
        pois_difftimes <- lapply(pois_history, diff)

        # Make plots and tests 
        if(return.plot)
        {
                p1 <- plot.trajectory.and.conditionalintensity(t=0, history=h_1, by_step=2)
                if(D == 1)
                        p1 <- list(plot=p1)
                force(p1)
        }

        envelope <- lapply(pois_history, plot.envelope.test)
        if(return.plot)
        {
                p2 <- lapply(envelope, `[[`, 'plot')
                force(p2)
        }

        kstest <- lapply(pois_difftimes, make.cdfplot.kstest)
        if(return.plot)
        {
                p3 <- lapply(kstest, `[[`, 'plot')
                force(p3)
        }

        if(return.plot)
        {
                .f <- function(i)
                        plot.martingale.residuals(pois_history[[i]], history=h_1[[i]])
                p4 <- lapply(1:D, .f) 
                force(p4)
        }

        out <- list(
                    KStest=sapply(kstest, `[[`, 'pvalue'),
                    envelopeout=sapply(envelope, `[[`, 'test')
        )
        
        if(return.plot)
        {
                .f <- function(i)
                {
                        
                        g <- ggarrange(p1[[i]], p2[[i]], p3[[i]], p4[[i]])
                        g <- annotate_figure(g, top=paste0('Dimension ', i ))

                        tmp <- paste0('allal_4diagnostics_', D, 'dims_dim', i, '.png')
                        filename=file.path(outdir_multi, tmp)
                        ggsave(g, filename=filename, width=12.5, height=10, units='cm')
                        g
                }
        
                ttheme(
                       base_size = 4,
                       base_style = "default",
                       base_colour = "black",
                       padding = unit(c(4, 4), "mm"),
                       colnames.style = colnames_style(size = 4),
                       rownames.style = rownames_style(size = 4),
                       tbody.style = tbody_style(size = 4)
                )


                summary_tabs <- .table.compare.params(PARS1=true_pars[1:3], PARS2=hat_pars[1:3])
                out1 <- list(plot=lapply(1:D, .f))
                out <- append(out, out1)
        }

        if(D == 2 & return.plot)
        {
                g <- ggarrange(
                        out$plot[[1]],
                        out$plot[[2]],
                        summary_tabs[[1]],
                        summary_tabs[[2]],
                        ncol=2, nrow=2,
                        heights=c(.9, .1))

                tmp <- paste0('allal_4diagnostics_', D, 'all_dims_estimates.png')
                filename=file.path(outdir_multi, tmp)
                ggsave(g, filename=filename, width=30, height=20, units='cm')
                g
        }

        out
}
