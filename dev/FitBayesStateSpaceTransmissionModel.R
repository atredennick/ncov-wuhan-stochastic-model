# FitBayesStateSpaceTransmissionModel.R


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(rjags)


# Load data ---------------------------------------------------------------

today <- Sys.Date()
start <- as.Date('12/01/2019',format='%m/%d/%Y')
today.day <- today - start + 1
#data <- read.csv('china-province-data - cases.csv')
data <- read.csv('china-province-data - wikipedia-cases.csv')
data[,1] <- as.Date(as.character(data[,1]), format='%m-%d-%Y')
names(data)[1] <- 'date'

data.province <- data #head(data,-3)
#data.province$date <- as.Date(data.province$Date, format='%m-%d-%Y')
data.province[is.na(data.province)] <- 0
#data.province$total <- rowSums(data.province[,3:35])
data.province$total <- rowSums(data.province[,c(2:32)])
data.province$day <- data.province$date - start + 1  # Day 1 is given by "start"
data.province$cum.cases <- cumsum(data.province$total)


# Plot the data -----------------------------------------------------------

plot(data.province$dat, cumsum(data.province$Hubei), type='h', lwd=5, col='grey',
     lend='butt', xlab='Date', ylab='Cumulative reported cases')
lines(data.province$date, data.province$Hubei)
points(data.province$date, data.province$Hubei, pch = 19, col = "white", cex = 1.25)
points(data.province$date, data.province$Hubei, pch = 21, cex = 1.25)


# Set up data for state-space model ---------------------------------------

yCumInfDet <- cumsum(data.province$Hubei)  # observed number of cumulative, reported infections


# Define functions for time-varying parameters ----------------------------

# Recovery rate
gamma <- function(z = 45, b=0.143, a0=0.0446, t){
  # linear function
  # default parameters z = 50, b=1/7, a0=0.3 from Paige Miller analysis; revised 2/1/20 based on updated analysis
  #    z: time at start of intervention (notionally Jan 19, 50 days after notional start on Dec 1 and the date of expanded testing in Wuhan)
  #    b: intercept (positive)
  #    a0: slope parameter (units of days, positive) -- original analysis had set to 0.03, but this doesn't seem fast enough
  #    t: time in the model
  
  gamma <- ifelse(t<=z, gamma <- b, gamma <- b + a0*(t-z))
  return(gamma)
}

# Case notification rate
eta <- function(t) ifelse(t<=55,1/(-0.47*t+27.2),1)

fever.clinic <- as.Date('2020-01-09')
w <- fever.clinic-start+1
q <- function(t, w=40, q0=0.11, q1=0.98) ifelse(t<=w,q0,q1)

#Function for making cells of projection matrix that hold parameters = NA
NA_matrix <- function(nrow, ncol, nyr){
  A <- array(0, dim=c(nrow, ncol, nyr))
  for(i in 1:nrow) {
    imax <- min(ncol, i+1)
    for(j in i:imax) {
      A[i, j, ] <- NA
    }
  }
  
  A[7, 12, ] <- NA
  A[11, 16, ] <- NA
  A[11, 12, ] <- 0
  A[16, 17, ] <- 0
  A[15, 16, ] <- 0
  A[15, 17, ] <- NA
  A[16, 18, ] <- NA
  A[17, 18, ] <- 0
  
  return(A)
}


# Set up data list for JAGS -----------------------------------------------

nClass <- 18
nSimDat <- length(yCumInfDet)
dt <- 1
nSimProc <- (nSimDat) / dt
steps <- 1:nSimProc
iDatInProc <- steps[seq(1, nSimProc, 1/dt)]

z <- 45
b <- 0.143
a0 <- 0.0446
gamma <- 4 * gamma(z=z, b=b, a0=a0, t=1:nSimProc)  # multiplier 4 for pseudo stages
sigma <- 1/6.4
sigma <- 6 * sigma  # multiplier 6 for pseudo stages
eta <- eta(t = 1:nSimProc)
qt <- q(t = 1:nSimProc, w = 40)

datList <- list(yCases = yCumInfDet,
                y.I = diag(1, nClass),
                nClass = nClass,
                dt = dt,
                nSimDat = nSimDat,
                nSimProc = nSimProc,
                iDatInProc = iDatInProc,
                A = NA_matrix(nrow = nClass, ncol = nClass, nyr = nSimProc),
                beta0 = 0.657, 
                sigma = sigma,
                gamma = gamma,
                b = sigma,
                eta = eta,
                q = qt,
                c = 1,
                S1 = 59002000 - 40,
                totNPopDat = 59002000)

n.iter=1000
n.update=1000
n.adapt=1000
model= "./dev/COVID-StateSpace-JAGS.R"
jm=jags.model(model, data=datList,n.chains=1, n.adapt = n.adapt)
update(jm,progress.bar="text", n.iter=n.update)
outMCMC=coda.samples(jm,variable.names=c("iTest"),n.iter=n.iter)
quants <- summary(outMCMC)$quantiles



I <- quants[,3]

xb <- 1:length(I)
plot(xb, I,xlab='',ylab='Cases',col=1, log = "y", type = "l")
day <- data.province$date - start
day <- xb
lines(day, cumsum(data.province$Hubei), type='h', col='lightskyblue', lwd=3, lend='butt' )
lines(xb, quants[,3])
lines(xb, quants[,2], lty = 2)
lines(xb, quants[,4], lty = 2)




# plot(I, type = "l")
# lines(quants[,1], lty =2)
# lines(quants[,5], lty = 2)
