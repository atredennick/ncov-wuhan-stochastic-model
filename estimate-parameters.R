# estimate-parameters.R
# Script to estimate model parameters using MCMC and synthetic
# likelihood, using the synlik package.
#
# Author: Andrew Tredennick


# Clear the decks ---------------------------------------------------------

rm(list = ls(all.names = TRUE))


# Load packages -----------------------------------------------------------

library(synlik)


# Source scripts ----------------------------------------------------------

source("stochastic-model.R")
source("simulators.R")



# Test the simulator ------------------------------------------------------

params <- list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, 
               c=1, presymptomatic=1, dt=0.05)

init <-  list(S=10600000, E1=0, E2=0, E3=0, E4=0, E5=6, E6=0,
              I1 = 1, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
              H=0, Ru=0, C=0)

nsims <- 2
nstep <- NULL
start <- as.Date("2020-03-01")
today <- Sys.Date()

res <- mod_simulator(params, init, nsims, nstep, start, today, extraArgs = NULL)
plot(res[[1]]$C)

extraArgs <- list("init" = init, "nstep" = nstep, "start" = start, "today" = today, "nObs" = 12)
res <- mod_wrap(param = log(unlist(params)), nsim = nsims, extraArgs = extraArgs)
matplot(t(res), type = "l")



# Create "synlik" object --------------------------------------------------

covid_sl <- new("synlik",
                simulator = mod_wrap,
                param = log(unlist(params)),
                extraArgs = extraArgs)

res <- simulate(covid_sl, nsim = 2)
matplot(t(res), type = "l")



