# Bayesian State-Space Model of COVID-19 transmission

var A[nClass, nClass, nSimProc]

model {
  # Priors and initial states
  totN[1] <- totNPopDat
  infectDet[1] <- 1
  infectUnDet[1] <- 0
  tau.p ~ dunif(0, 10)
  
  N[1, 1] ~ dunif(S1, S1+2)
  for(ii in 2:8){
    N[1, ii] ~ dunif(0, 1)
  }
  N[1, 9] ~ dunif(1, 2)
  for(ii in 10:nClass){
    N[1, ii] ~ dunif(0, 1)
  }

  
  for(t in 2:nSimProc) {
    # Rate of S -> I
    phi[t] <- beta0*infectDet[t-1] / totN[t-1]+beta0*c*infectUnDet[t-1] / totN[t-1]
    
    # Fill in transition matrix
    # susceptible to exposed
    A[1, 1, t] <- exp(-phi[t]*dt)  # S -> E
    A[1, 2, t] <- 1 - exp(-phi[t]*dt)  # S-> E
    
    # Transitions between exposed classes and out
    A[2, 2, t] <- exp(-sigma*dt)
    A[2, 3, t] <- 1 - exp(-sigma*dt)
    A[3, 3, t] <- exp(-sigma*dt)
    A[3, 4, t] <- 1 - exp(-sigma*dt)
    A[4, 4, t] <- exp(-sigma*dt)
    A[4, 5, t] <- 1 - exp(-sigma*dt)
    A[5, 5, t] <- exp(-sigma*dt)
    A[5, 6, t] <- 1 - exp(-sigma*dt)
    A[6, 6, t] <- exp(-sigma*dt)
    A[6, 7, t] <- 1 - exp(-sigma*dt)
    A[7, 7, t] <- exp(-sigma*dt)
    A[7, 8, t] <- (1 - exp(-sigma*dt)) * q[t]
    A[7, 12, t] <- (1 - exp(-sigma*dt)) * (1-q[t])
    
    # Transitions between and out of I (detected)
    A[8, 8, t] <- exp(-gamma[t]*dt)
    A[8, 9, t] <- 1 - exp(-gamma[t]*dt)
    A[9, 9, t] <- exp(-gamma[t]*dt)
    A[9, 10, t] <- 1 - exp(-gamma[t]*dt)
    A[10, 10, t] <- exp(-gamma[t]*dt)
    A[10, 11, t] <- 1 - exp(-gamma[t]*dt)
    A[11, 11, t] <- exp(-gamma[t]*dt)
    A[11, 16, t] <- 1 - exp(-gamma[t]*dt)
    
    # Transitions between and out of I (undetected)
    A[12, 12, t] <- exp(-b*dt)
    A[12, 13, t] <- 1 - exp(-b*dt)
    A[13, 13, t] <- exp(-b*dt)
    A[13, 14, t] <- 1 - exp(-b*dt)
    A[14, 14, t] <- exp(-b*dt)
    A[14, 15, t] <- 1 - exp(-b*dt)
    A[15, 15, t] <- exp(-b*dt)
    A[15, 17, t] <- 1 - exp(-b*dt)
    
    # Recoveries
    A[16, 16, t] <- exp(-eta[t]*dt)
    A[16, 18, t] <- 1 - exp(-eta[t]*dt)
    
    A[17, 17, t] <- exp(-b*dt)
    A[18, 18, t] <- exp(-eta[t]*dt)
    
    # Estimates of mean of posterior distribution
    mu[t, 1:nClass] <- A[,,t] %*% N[t-1, 1:nClass]
    
    # Get logs of deterministic predictions
    for(j in 1:nClass){
      log_mu[t,j] <- log(max(1, mu[t, j]))
    }
    
    log_N[t, 1:nClass] ~ dmnorm(log_mu[t,1:nClass], tau.p*y.I) 
    
    # Exponentiate log of population size
    for(j in 1:nClass){
      N[t,j] <- min(exp(log_N[t,j]),10000)
    }
    totN[t] <- sum(N[t,1:nClass])
    
    # Likelihood based on case incidence
    infectUnDet[t] <- sum(N[t, 13:16])  # undetected infections
    infectDet[t] <- sum(N[t, 9:12])  # detected infections
  }  # end time step
  
  for(lt in 2:nSimDat) {
    yCases[lt] ~ dpois(infectDet[iDatInProc[lt]])  # likelihood of data given the model
  }
}  # end model
/*End of Model*/