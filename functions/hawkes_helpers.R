plot_multi_hawkes <- function(lst, max=Inf)
{
        # lst <- copy(h)
        tmp <- as.data.table(lst)
        tmp <- melt(tmp)
        tmp <- tmp[value <= max]
        
        ggplot(tmp, aes(x=value, y=variable, color=variable)) +
                geom_point(pch='+') + 
                theme_bw() +
                labs(x='Event times', y='Dimension', color='') + 
                theme(legend.position='none')
}

maximum_likelihood_hawkes <- function( t )
{
        D <- 1
        if(is.list(t))
                D <- length(t)

        if(D==1)
                t <- t[[1]]

        # define likelihood

        idx_lambda0 <- seq(1, D, by=1)
        idx_alpha <- seq(D+1, D + D^2, by=1)
        idx_beta <- seq(D+D^2+1, D + 2*D^2, by=1)
        
        .f <- function(x)
        {
                .m <- function(x, idx)
                {
                        if (D == 1)
                        {
                                return(x[idx])
                        }else{
                                return(matrix(x[idx], byrow=TRUE, nrow=D))
                        }
                }

                likelihoodHawkes(lambda0 = x[idx_lambda0],
                                 alpha = .m(x, idx_alpha),
                                 beta = .m(x, idx_beta),
                                 history=t)
        }
        
        x0 <- c(
                rep(0.5, length(idx_lambda0)),
                as.vector(.6 * diag(D)),
                as.vector(0.8 * matrix(1, byrow=TRUE, nrow=D, ncol=D))
                )


        # perform constrained optimization if D == 1,
        # else can't be bothered (and is it possible?) to get ui and ci
        ui <- 0
        if(D == 1)
                ui <- c(0, -1, 1)
        
        
        out <- constrOptim(theta=x0,  f=.f, grad=NULL, ui=ui, ci=0-.0001)

        with(out,
             list(
                  lambda0=par[idx_lambda0],
                  alpha=matrix(par[idx_alpha], byrow=T, D),
                  beta=matrix(par[idx_beta], byrow=T, D ),
                  loglik=value
             )
        )
}

theme_text_sizes <- theme( 
                          text = element_text(size=6),
                          axis.text = element_text(size=5),
                          axis.title = element_text(size=6),
                          plot.title = element_text(size=8))

recreate_allal <- function(lambda0=LAMBDA0, alpha=ALPHA, beta=BETA, horizon=HORIZON, return.plot=TRUE)
{
        # Deprecated: the multidimensional function also works for 1D

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


make.cdfplot.kstest <- function(poisson_deltat)
{
        tmp <- ks.test(poisson_deltat, pexp) 
        lab <- data.table( lab=paste0('\t KS-test: \t\n\tp-value: ', round(tmp$p.value, 2), '\t\n' ))
        lab[, `:=` (xpos=Inf, ypos=-Inf, hjustvar=1, vjustvar=0)]


        p <- ggplot() + 
                stat_function( fun=pexp, color='black') +
                stat_ecdf( aes(x=poisson_deltat), geom='step', color='blue' ) +
                geom_hline(aes(yintercept=c(0,1)), linetype='dotted') +
                geom_text(data=lab, aes(label=lab, x=xpos, y=ypos, hjust=hjustvar, vjust=vjustvar), size=3) + 
                
                theme_bw() +
                theme_text_sizes +
                scale_y_continuous(expand=c(0, 0)) +
                scale_x_continuous(expand=c(0, 0)) +
                labs(x='Lambda(t)', y='CDF', title='Cumulative Distribution functions')

        list(plot=p, pvalue=tmp$p.value)
}

plot.trajectory.and.conditionalintensity <- function(t, history, by_step=.1)
{
        if(0)
        {
                history=h_1
                by_step=2
        }

        # t seems to be useless? maybe its not in 1D
        tmp <- data.table(x=seq(from=0, to=horizon, by=by_step))
        cols <- paste0('y', 1:D)
        tmp[, (cols) := compensator(x=x, pars=hat_pars, history=history)]
        
        .f <- function(col)
        {
                dim <- as.integer(gsub('y([0-9]+)$', '\\1', col))
                hstr <- history[[dim]]

                sel <- c('x', col)
                tmp1 <- tmp[, ..sel]
                setnames(tmp1, col, 'y')
                tmp1[, rescale:=length(hstr)/max(y)]
                tail(tmp1)

                p <- ggplot() + 
                        geom_line(data=tmp1, aes(x=x, y=y*rescale), color='grey80') + 
                        geom_step(aes(x=hstr,y=1:length(hstr)), color='red') +
                        theme_bw()  +
                        theme_text_sizes +
                        scale_y_continuous(expand=c(0, 0)) +
                        scale_x_continuous(expand=c(0, 0)) +
                        labs(x='t', y='', title='Trajectory and conditional intensity')  
                force(p)
                p
        }
        out <- lapply(cols, .f)

        if(length(out) == 1)
                out <- out[[1]]

        return(out)
}

plot.envelope.test <- function(poisson_t)
{

        # envelope test

        # the steps of the PP are Gamma(N, 1) random variables
        # so we can get quantiles

        .f <- function(n)
                data.table(
                        Nt=n,
                        cc_025=qgamma(c(.025), shape= n, rate=1),
                        cc_975=qgamma(c(.975), shape= n, rate=1)
                )

        drange <- lapply(1:length(poisson_t), .f)
        drange <- rbindlist(drange)
        drange[, t_hat := poisson_t ]

        drange[, DUMMY := (t_hat <= cc_025 | t_hat >= cc_975) ]
        test <- drange[, any(DUMMY)]

        p <- ggplot(drange, aes(y=Nt)) + 
                geom_ribbon(aes(xmin=cc_025, xmax=cc_975), fill='green', alpha=.4) +
                geom_step(aes(x=t_hat)) +
                theme_bw() +
                theme_text_sizes +
                scale_y_continuous(expand=c(0, 0)) +
                scale_x_continuous(expand=c(0, 0)) +
                labs(x='Lambda(t)', y='N_t', title='Envelope Test') 
        list(plot=p, test=test)
}

plot.martingale.residuals <- function(poisson_t, history)
{
        res <- data.table(
                          t = history,
                          N_t = seq_along(poisson_t),
                          Lambda_t = poisson_t)
        res[, R_t := N_t - Lambda_t ]

        p <- ggplot(res, aes(x=t, y=R_t)) + 
                geom_line() +  
                theme_bw() +
                theme_text_sizes +
                scale_y_continuous(expand=c(0, 0)) +
                scale_x_continuous(expand=c(0, 0)) +
                labs(x='time', y='R(t) = N(t) - Lambda(t)', title='Martingale Residuals') 
        p
}

.table.compare.params <- function(PARS1, PARS2)
{
        .get.contributions.to.dim.i <- function(i, lst)
        {
                out <- lapply(lst, function(inp){
                               if(is.matrix(inp)){
                                       round(inp[, i],2)
                               }else{
                                       round(inp[i], 2)
                               }
                          })

                n_times <- lapply(out, length)
                out <- Reduce(c, out)
                out

                # get names
                nms <- names( unlist(lapply(n_times, function(n) 1:n ) ))
                nms <- gsub('([a-z])([0-9])', '\\1_\\2', nms)
                nms <- paste0(nms, i)

                out1 <- copy(out)
                names(out1) <- nms
                t(out1)
        }

        true <- lapply(1:D, .get.contributions.to.dim.i, lst=PARS1)
        est <- lapply(1:D, .get.contributions.to.dim.i, lst=PARS2)

        .f <- function(dim)
        {
                out <- rbind(
                        true[[dim]],
                        est[[dim]])

                out <- cbind( c('True:', 'Est:'), out)
                out
        }
        
        out <- lapply(1:D, .f)
        out <- lapply(out, ggtexttable)
}
