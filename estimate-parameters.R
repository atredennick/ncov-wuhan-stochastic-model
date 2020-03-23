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



# Read in data ------------------------------------------------------------

today <- Sys.Date()
start <- as.Date('12/01/2019',format='%m/%d/%Y')
today.day <- today - start + 1
data <- read.csv('china-province-data - wikipedia-cases.csv')
data[,1] <- as.Date(as.character(data[,1]), format='%m-%d-%Y')
names(data)[1] <- 'date'
data.province <- data #head(data,-3)
data.province[is.na(data.province)] <- 0
data.province$total <- rowSums(data.province[,c(2:32)])
data.province$day <- data.province$date - start + 1  # Day 1 is given by "start"
data.province$cum.cases <- cumsum(data.province$total)

model_data <- data.province[, c("date", "Hubei")]
model_data$cumcases <- cumsum(model_data$Hubei)


# Create "synlik" object --------------------------------------------------

params <- list(beta0=0.657, sigma=1/6.4, b=0.143, a0=0.0446, beta.factor = 1)

init <-  list(S=59002000, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
              I1 = 1, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
              H=0, Ru=0, C=0)

nsims <- 2
nstep <- NULL
start <- as.Date("2019-12-01")
today <- Sys.Date()

extraArgs <- list("init" = init, "nstep" = nstep, "start" = start, 
                  "today" = today, "dt" = 0.05, "w" = 40, "z" = 45, "c" = 1, 
                  "presymptomatic" = 0, "timesToObs" = TRUE, 
                  "nObs" = nrow(model_data), "obsData" = model_data$cumcases,
                  "obsDataDate" = model_data$date)

covid_sl <- new("synlik",
                simulator = mod_wrap,
                param = log(unlist(params)),
                extraArgs = extraArgs)

covid_sl@data <- model_data$cumcases

res <- simulate(covid_sl, nsim = 5)
# matplot(t(res), typ= "l", log = "y", ylim = c(1, 800000))


plot(model_data$date, model_data$cumcases, 
     type='h', lwd=5, lend='butt', xlab='', 
     col = col.cases, ylim = c(1, 800000),
     ylab='Cumulative case notifications', 
     main='COVID-19 cases in Hubei', log = "y")
for(i in 1:nrow(res)) {
  resDF <- data.frame(date = model_data$date,
                      cases = res[i, ])
  lines(resDF$date, resDF$cases, col = col.cases.ci)
}


# Define prior distributions ----------------------------------------------

covid_priors_sl <- function(input, ...) { sum( input ) +
    dnorm (input[1], log(0.67), 1, log = TRUE) +
    dnorm (input[2], log(0.15), 1, log = TRUE) +
    dnorm (input[3], log(0.143), 1, log = TRUE) +
    dnorm (input[4], log(0.67), 1, log = TRUE) }


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

