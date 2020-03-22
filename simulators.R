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
  nObs <- extraArgs$nObs
  init <- extraArgs$init
  start <- extraArgs$start
  today <- extraArgs$today
  stopifnot( !is.null(nObs) )
  
  # Return simul, with nrow = nsims and ncol = nsteps
  simul <- t(sapply(mod_simulator(params = param, init , nsims = nsim, nstep, start, today,
                         extraArgs = extraArgs), function(x) x$C))
  
  # Save only the time steps with indexes in extraArgs$timesToSave
  if( !is.null(extraArgs$timesToSave) ){ 
    return( simul[ , extraArgs$timesToSave, drop = FALSE] )
  } else {
    return( simul )
  }
}