maximum_likelihood_hawkes <- function( t )
{
        # define likelihood
        .f <- function(x)
                likelihoodHawkes(x[1], x[2], x[3], t)

        # define 
        out <- constrOptim(c(0.5, 0.5, 2),  .f, grad=NULL, ui=c(0, -1, 1), ci=0)

        with(out,
             list(
                  lambda0=par[1],
                  alpha=par[2],
                  beta=par[3],
                  loglik=value
             )
        )
}

theme_text_sizes <- theme( 
                          text = element_text(size=6),
                          axis.text = element_text(size=5),
                          axis.title = element_text(size=6),
                          plot.title = element_text(size=8))

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

        list(plot=p, pvalue=p)
}

plot.trajectory.and.conditionalintensity <- function(history)
{

        tmp <- data.table(x=seq(from=0, horizon, by=.1))
        tmp[, y:=compensator(x=x ,hat_pars, history=history)]
        tmp[, rescale:=length(history)/max(y)]

        p <- ggplot() + 
                geom_step(aes(x=history, y=1:length(history)), color='red') +
                geom_line(data=tmp, aes(x=x, y=y*rescale)) + 
                theme_bw()  +
                theme_text_sizes +
                scale_y_continuous(expand=c(0, 0)) +
                scale_x_continuous(expand=c(0, 0)) +
                labs(x='t', y='', title='Trajectory and conditional intensity')  
        p
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

plot.martingale.residuals <- function(poisson_t=pois, history)
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
