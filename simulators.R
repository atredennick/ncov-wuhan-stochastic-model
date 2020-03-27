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
  names(param) <- extraArgs$paramNames
  
  # Set some extraArgs variables as objects or objects part of param
  # for passing to evaluate.model. This is done because synlik only
  # wants parameters to be estimated in the param argument. (Well, there
  # are ways around that, but this is easier.)
  if(!is.null(extraArgs$nObs)) nObs <- extraArgs$nObs
  init <- extraArgs$init
  init$I1 <- ceiling(param$I0)
  param$I0 <- NULL
  start <- extraArgs$start
  today <- extraArgs$today
  param$dt <- extraArgs$dt
  param$w <- extraArgs$w
  param$z <- extraArgs$z
  param$c <- extraArgs$c
  param$presymptomatic <- extraArgs$presymptomatic
  if(!("sigma" %in% extraArgs$paramNames)) param$sigma <- extraArgs$sigma
  param$b <- extraArgs$b
  param$a0 <- extraArgs$a0
  nstep <- extraArgs$nstep
  
  # stopifnot( !is.null(nObs) )
  
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
    first_obs <- which(seq(start, today + 28, by = 1) == (head(extraArgs$obsDataDate, 1) - 1))
    last_obs <- which(seq(start, today + 28, by = 1) == tail(extraArgs$obsDataDate, 1))
    s1 <- simul[ , report_times[first_obs:last_obs], drop = FALSE]
    s1 <- apply(s1, 1, diff)
    # idx <- which(colSums(s1) != 0)
    return( t(s1) )
  } else {
    if(nrow(simul) == 1) {
      return(diff(simul[ , report_times]))
    } else {
      return( t(apply(simul[ , report_times], 1, diff)) )
    }
  }
}