# simulators.R
# These are wrapper functions to simulate from the stochastic model,
# specifically designed to work with the synlik package for parameter
# estimation via MCMC and synthetic likelihood.

mod_simulator <- function(params, init, nsims, nstep, start, today,
                          extraArgs = NULL, ...) {
  
  tmp <- evaluate.model(params, init , nsims, nstep, start, today)
  return(tmp)
}


mod_wrap <- function(param, nsim, extraArgs = extraArgs, ...) {
  param <- as.list(exp(param))
  names(param) <- c( "beta0", "sigma", "b", "a0", "beta.factor")
  
  # Set some extraArgs variables as objects or objects part of param
  # for passing to evaluate.model. This is done because synlik only
  # wants parameters to be estimated in the param argument. (Well, there
  # are ways around that, but this is easier.)
  nObs <- extraArgs$nObs
  init <- extraArgs$init
  start <- extraArgs$start
  today <- extraArgs$today
  param$dt <- extraArgs$dt
  param$w <- extraArgs$w
  param$z <- extraArgs$z
  param$c <- extraArgs$c
  param$presymptomatic <- extraArgs$presymptomatic
  
  stopifnot( !is.null(nObs) )
  
  # Return simul, with nrow = nsims and ncol = nsteps
  simul <- t(
    sapply(
      mod_simulator(params = param, init , nsims = nsim, nstep, 
                    start, today, extraArgs = extraArgs), 
      function(x) x$C
      )
    )
  
  if(is.null(extraArgs$timesToObs)) extraArgs$timesToObs <- FALSE
  
  # Save only the time steps with indexes in extraArgs$timesToSave
  report_times <- seq(from = 1, to = ncol(simul), by = 1/extraArgs$dt)
  if( extraArgs$timesToObs ){ 
    first_obs <- which(seq(start, today + 28, by = 1) == head(extraArgs$obsDataDate, 1))
    last_obs <- which(seq(start, today + 28, by = 1) == tail(extraArgs$obsDataDate, 1))
    return( simul[ , report_times[first_obs:last_obs], drop = FALSE] )
  } else {
    return( simul )
  }
}