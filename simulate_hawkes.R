homogeneous_poisson_proc <- function(intensity,time_span){
  wait_times <- vector("list")
  wait_times[[1]] <- 0
  while(Reduce("+",wait_times) < time_span){
    samp <- rexp(n = 1, rate = intensity)
    wait_times <- append(wait_times,samp)
  }
  event_times <- cumsum(unlist(wait_times))
  return(event_times[2:(length(event_times)-1)])
}

trunc_unif_poisson_proc <- function(average_reprod_rate, reprod_period, time_span){
  mean_points <- average_reprod_rate*min(reprod_period, time_span)
  total_points <- rpois(1, mean_points)
  event_times <- runif(total_points, min = 0, max = min(reprod_period,time_span))
  return(sort(event_times))
}

exponential_poisson_proc <- function(average_intensity, exponential_rate, time_span){
  # intensity(t) = average_intensity*exponential_rate^exp(-exponential_rate*t)
  mean_points <- pexp(time_span, exponential_rate)*average_intensity
  total_points <- rpois(1, mean_points)
  # simulating exponentials conditional on being within timespan
  # use inverse CDF sampling and simulate uniforms on restricted interval (eps,1]
  # exponential = -log(uniform)/rate must be less than time_span
  # alternatively can generate on whole line then truncate afterwards
  unif_cutoff <- exp(-time_span*exponential_rate)
  conditional_unifs <- runif(total_points, min = unif_cutoff, max = 1)
  event_times <- -log(conditional_unifs)/exponential_rate
  return(sort(event_times))
}

# Simulating population model with immigration and birth:
# See Example 2 on page 5 of 'Bayesian Inference for Hawkes Processes' (Rasumussen)

# setting example parameters
end_time <- 100 # time span on which to consider the HP
background_intensity <- 1 # intensity at which immigrants arrive
mark_rate <- 1 # rate of the exponential dist used to sample reproduction period
average_reprod_rate <- 0.5 # rate of reproduction per unit time

simulate_population_hawkes <- function(
  end_time,background_intensity,mark_rate,average_reprod_rate,max_gen=20) {
# initial cluster from homogenous poisson process
  immigrant_times <- homogeneous_poisson_proc(background_intensity, time_span=end_time)
  immigrant_marks <- rexp(length(immigrant_times), rate = mark_rate)
  immigrant_gen <- rep_len(0, length(immigrant_times))
  times_marks <- cbind(immigrant_times,immigrant_marks, immigrant_gen)
  colnames(times_marks) <- c("time","mark","gen")
  times_marks <- data.table(times_marks)

  # descendant clusters generated sequentially
  gen <- 0
  start_index <- 1
  stop_index <- nrow(times_marks)
  # descendants of immigrants
  while ( (stop_index - start_index > 0) & (gen<max_gen) ){ # loop until no new descendants
    gen <- gen + 1 # in practice generation is not observable
    for ( i in start_index:stop_index ){ # Add in desecendants of current generation
      new_times <- trunc_unif_poisson_proc(average_reprod_rate =  average_reprod_rate,
      reprod_period = times_marks[i,2],
      time_span = end_time - times_marks[i,1])
      new_marks <- rexp(length(new_times),rate = mark_rate)
      new_gens <- rep_len(gen, length(new_times))
      new_times_marks <- matrix(c(new_times,new_marks,new_gens),ncol=3)
      times_marks <- rbind(times_marks, new_times_marks, use.names=FALSE)
    }
    # on next iteration update indices to consider next generation
    start_index <- stop_index + 1
    stop_index <- nrow(times_marks)
  }
  return(times_marks) # may wish to sort by time
}
