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
source('style-definitions.R')
source("my-smcmc.R")


# Pull in data from GitHub ------------------------------------------------

start <- as.Date("2020-03-01")
state.data <- read.csv('https://raw.githubusercontent.com/CEIDatUGA/COVID-19-DATA/master/UScases_by_state_wikipedia.csv?token=ABQMQXYOBV7TFRDQXWX2YQC6PEFIU', stringsAsFactors = FALSE)
state.data <- state.data[1:(dim(state.data)[1]-1),]
state.data$Date <- as.Date(state.data$Date, format='%b %d')
state.data[is.na(state.data)] <- 0
georgia <- data.frame(date=state.data$Date, cases=state.data$GA)
georgia <- subset(georgia, date >= start)


# Create "synlik" object --------------------------------------------------

params <- list(beta0=0.6584, sigma=1/6.4, b=0.143, a0=1/1.5)

init <-  list(S=10600000, E1=0, E2=0, E3=0, E4=0, E5=6, E6=0,
              I1 = 40, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
              H=0, Ru=0, C=0)

nsims <- 2
nstep <- NULL
start <- as.Date("2020-03-01")
today <- Sys.Date()

extraArgs <- list("init" = init, "nstep" = nstep, "start" = start, 
                  "today" = today, "dt" = 0.05, "w" = 12, "z" = 12, "c" = 1, 
                  "presymptomatic" = 1, "timesToObs" = TRUE)
extraArgs$nObs <- nrow(georgia)

covid_sl <- new("synlik",
                simulator = mod_wrap,
                param = log(unlist(params)),
                extraArgs = extraArgs)

covid_sl@data <- cumsum(georgia$cases)
covid_sl@extraArgs$obsData <- cumsum(georgia$cases)

res <- simulate(covid_sl, nsim = 5)
plot(georgia$date, cumsum(georgia$cases), 
     type='h', lwd=5, lend='butt', xlab='', 
     col = col.cases,
     ylab='New case notifications', 
     main='COVID-19 cases in Georgia')
for(i in 1:nrow(res)) {
  resDF <- data.frame(date = georgia$date,
                      cases = res[i, ])
  lines(resDF$date, resDF$cases, col = col.cases.ci)
}


# Define prior distributions ----------------------------------------------

covid_priors_sl <- function(input, ...){ sum( input ) +
    dnorm (input[1], log(0.67), 1, log = TRUE) +
    dnorm (input[2], log(0.15), 1, log = TRUE) +
    dnorm (input[3], log(0.143), 1, log = TRUE) +
    dnorm (input[4], log(0.67), 1, log = TRUE)}


# Define stats for synthetic likelihood -----------------------------------

covid_stats <- function(x, extraArgs, ...) {
  obsData <- as.vector(extraArgs$obsData)
  stopifnot(length(obsData) != 0)
  if (!is.matrix(x)) 
    x <- matrix(x, 1, length(x))
  x <- t(x)
  X0 <- cbind(log(colMeans(x)+0.0001), log(colSums(x)+0.0001))
  X0
}

covid_sl@summaries <- covid_stats
# checkNorm(covid_sl, nsim = 50)  # works


# Run MCMC ----------------------------------------------------------------

covMat <- diag(rep(0.2, length(params)))

# slik(covid_sl, param  = unname(log(unlist(params))), nsim   = 10)  # works
# slice(object = covid_sl,
#       ranges = list("beta0" = seq(-0.5, -0.3, by = 0.01)),
#       param = log(unlist(params)),
#       nsim = 10)  # works

init_params <- c(beta0 = 0.3, sigma = 0.2, b = 0.15, a0 = 0.4)
covid_mcmc <- smcmc(covid_sl, 
                  initPar = log(init_params),
                  nsim = 100,
                  niter = 50, 
                  burn = 0,
                  priorFun = covid_priors_sl, 
                  propCov = covMat)

plot(covid_mcmc)

