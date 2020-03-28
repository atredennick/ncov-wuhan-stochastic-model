
## ----create-sims---------------------------------------------------------
# First source the necessary functions
source("stochastic-model.R")  # the stochastic disease transmission model
source("simulators.R")  # simulator wrappers
source("style-definitions.R")  # colors

# Time spans for simulation
start <- as.Date("2019-12-01")
today <- Sys.Date()

# Unknown parameters
theta <- list(beta0 = 0.657, beta.factor = 1, sigma = 1/6.4, I0 = 1)

# Known parameters
extraArgs <- list(nstep = NULL, start = start, 
                  today = today, dt = 0.05, w = 40, z = 45, c = 1, 
                  presymptomatic = 0, timesToObs = FALSE, 
                  b = 0.143, a0 = 0.0446, paramNames = names(theta), 
                  nObs = NULL, obsDataDate = NULL)

# Initial conditions
init <- list(S=59002000, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
             I1 = NA, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
             H=0, Ru=0, C=0)
extraArgs$init <- init

# Simulate from the model
res <- mod_wrap(param = log(unlist(theta)), nsim = 2, extraArgs = extraArgs)

matplot(t(res), type = "l", col = col.cases.ci, lty = 1, 
        main = "Simulated trajectories of COVID-19",
        xlab = "Day", ylab = "Number of new cases")


## ----sim-data------------------------------------------------------------
obs_data <- mod_wrap(param = log(unlist(theta)), nsim = 1, extraArgs = extraArgs)
extraArgs$obsData <- obs_data
extraArgs$nObs <- length(obs_data)
plot(obs_data, type='b',  xlab='Day', col = col.cases,
     ylab='Number of new cases', main='Simulated COVID-19 cases in Hubei')


## ----plot-ma-------------------------------------------------------------
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
stat_times <- ma(diff(log(obs_data+10)), n = 5)
d0 <- min(which(stat_times > 0.1))
d1 <- d0 + min(which(stat_times[d0:length(stat_times)] < -0.01))
par(mfrow = c(1,2))
plot(stat_times, xlab = "Day", ylab = "Moving average of differences", 
     type = "b", cex = 0.5)
abline(h = 0, col = "grey", lty = 2)
abline(v = c(d0, d1), col = "red")
plot(obs_data, type='b',  xlab='Day', col = col.cases,
     ylab='Number of new cases', cex = 0.5)
abline(v = c(d0, d1), col = "red")


## ----summary-stats-------------------------------------------------------
covid_stats <- function(x, extraArgs, ...) {
  ## obsData is a vector of observed path
  ## x is a M by n.t matrix of paths, each row of 'x' is a replicate
  
  ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}
  stat_times <- ma(diff(log(extraArgs$obsData+10)), n = 5)
  d0 <- min(which(stat_times > 0.1))
  d1 <- d0 + min(which(stat_times[d0:length(stat_times)] < -0.01))
  
  obsData <- extraArgs$obsData
  
  stopifnot(is.vector(obsData), length(obsData) != 0)
  if (!is.matrix(x)) x <- matrix(x, 1, length(x))
  
  x <- log(x + 10)  # log transform
  
  # Max incidence
  X0 <- (apply(x, 1, function(x) max(x[1:(d0-1)])))
  X0 <- cbind(X0, (apply(x, 1, function(x) max(x[d0:d1]))))
  X0 <- cbind(X0, (apply(x, 1, function(x) max(x[(d1+1):length(x)]))))
  
  # Cumulative incidence
  X0 <- cbind(X0, (apply(x, 1, function(x) sum(x[1:(d0-1)]))))
  X0 <- cbind(X0, (apply(x, 1, function(x) sum(x[d0:d1]))))
  X0 <- cbind(X0, (apply(x, 1, function(x) sum(x[(d1+1):length(x)]))))
  
  # Max along whole trajectory
  X0 <- cbind(X0, apply(x, 1, max))
  
  # Day of max
  X0 <- cbind(X0, apply(x, 1, function(x) which.max(x)))
  
  # Exponential increase and decline
  X0 <- cbind(X0, apply(x, 1, max) / (apply(x, 1, function(x) which.max(x)) - d0))
  
  # Linear regression coefs
  X0 <- cbind(X0, apply(x, 1, function(x) as.numeric(coef(lm(x[1:(d0-1)] ~ seq_along(x[1:(d0-1)])))[2])))
  X0 <- cbind(X0, apply(x, 1, function(x) as.numeric(coef(lm(x[d0:d1] ~ seq_along(x[d0:d1])))[2])))
  X0 <- cbind(X0, apply(x, 1, function(x) as.numeric(coef(lm(x[(d1+1):length(x)] ~ seq_along(x[(d1+1):length(x)])))[2])))
  
  # Final epidemic size
  X0 <- cbind(X0, apply(x, 1, sum))
  return(X0)
}


## ----ex-stats------------------------------------------------------------
summaries <- covid_stats(x = res, extraArgs = extraArgs)
print(summaries)


## ----obs-stats-----------------------------------------------------------
obs_summaries <- covid_stats(x = obs_data, extraArgs = extraArgs)
print(obs_summaries)


## ----dist----------------------------------------------------------------
d <- function(xsim, xobs) {
  mean(abs(xobs - xsim))
}

d(summaries[2,], obs_summaries)

dsims <- 500
dist <- numeric(dsims)
for(i in 1:dsims) {
  res <- mod_wrap(param = log(unlist(theta)), nsim = 1, extraArgs = extraArgs)
  summaries <- covid_stats(x = res, extraArgs = extraArgs)
  dist[i] <- d(summaries, obs_summaries)
}
hist(dist, main = "Distribution of statistic distances at known parameters",
     xlab = "d")
abline(v = median(dist), col = "red")
epsilon <- as.numeric(round(quantile(dist, 0.8)))


## ----abc-----------------------------------------------------------------
get_prior_d <- function(theta, theta_prime) {
  d_theta <- dnorm(theta[1], log(0.657), 0.5) +
             dnorm(theta[2], log(1), 0.5) +
             dnorm(theta[3], log(1/6,4), 0.5) +
             dunif(theta[4], -5, 5)
  
  d_theta_prime <- dnorm(theta_prime[1], log(0.657), 0.5) +
                   dnorm(theta_prime[2], log(1), 0.5) +
                   dnorm(theta_prime[3], log(1/6,4), 0.5) +
                   dunif(theta_prime[4], -5, 5)
  return(min(1, exp(d_theta_prime - d_theta)))
}


run_abc <- function(nburn, nsamp, epsilon, f_stats, f_distance, f_model, f_prior,
                    init_theta, obs_data, ...) {
  results <- matrix(nrow = nsamp, ncol = length(init_theta))
  i <- 0
  theta <- init_theta
  while(i < (nburn+nsamp)) {
    theta_prime <- rnorm(length(theta), theta, sd = c(0.1, 0, 0, 0))
    # p <- f_prior(as.numeric(theta), theta_prime)
    # if(rbinom(1, 1, p) == 1) {
      m <- f_model(param = theta_prime, nsim = 1, extraArgs = extraArgs)
      mS <- f_stats(x = m, extraArgs = extraArgs)
      mO <- f_stats(x = obs_data, extraArgs = extraArgs)
      dtest <- f_distance(mS, mO)
      if(dtest < epsilon) {
        theta <- theta_prime
      }
    # }
#     propsal = c
#     accept with p exp(prior(propsal) – prior(current))
#     if accept do ABC step
    i <- i + 1
    if(i > nburn) results[(i - nburn), ] <- theta
  }
  return(results)
}


## ----test----------------------------------------------------------------
init_theta <- log(c(beta0 = 0.5, beta.factor = 1, sigma = 1/6.4, I0 = 1))
out_mcmc <- run_abc(nburn = 2000, nsamp = 5000, epsilon = epsilon, f_stats = covid_stats, 
                    f_distance = d, f_model = mod_wrap, f_prior = get_prior_d,
                    init_theta = init_theta, obs_data = obs_data, extraArgs)

saveRDS(object = out_mcmc, file = "mcmc-chains-abc.RDS")

# Plot the chains
# png("mcmc-chains.png", width = 8.5, height = 6, units = "in", res = 300)
# par(mfrow = c(2, 2))
# plot(exp(out_mcmc[,1]), type = "l", xlab = "Iteration",
#      ylab = expression(log(beta)), las = 1)
# abline(h = as.numeric(theta[1]), col = "red")
# plot(out_mcmc[,2], type = "l", xlab = "Iteration",
#      ylab = expression(log(xi)), las = 1)
# abline(h = log(as.numeric(theta[2])), col = "red")
# plot(out_mcmc[,3], type = "l", xlab = "Iteration",
#      ylab = expression(log(sigma)), las = 1)
# abline(h = log(as.numeric(theta[3])), col = "red")
# plot(out_mcmc[,4], type = "l", xlab = "Iteration",
#      ylab = expression(log(I[0])), las = 1)
# abline(h = as.numeric(theta[4]), col = "red")
# dev.off()

png("mcmc-chains.png", width = 8.5, height = 4, units = "in", res = 300)
par(mfrow = c(1,2))
plot(exp(out_mcmc[,1]), type = "l", xlab = "MCMC Iteration",
          ylab = expression(beta), las = 1, bty = "n", main = "Traceplot of ABC-MCMC")
abline(h = as.numeric(theta[1]), col = "red", lwd = 2)
hist(exp(out_mcmc[,1]), breaks = 15, col = "grey", xlab = expression(beta),
     las = 1, main = "Posterior Distribution")
abline(v = as.numeric(theta[1]), col = "red", lwd = 2)
dev.off()

