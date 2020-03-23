# stochastic-model.R
# Functions to simulate from the stochastic disease tranmission 
# model for COVID-19.
#
# Author: John Drake

onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  
  S <- x[2]                            #local variable for susceptibles
  
  E1 <- x[3]                           #exposed classes
  E2 <- x[4]
  E3 <- x[5]
  E4 <- x[6]
  E5 <- x[7]
  E6 <- x[8]
  
  I1 <- x[9]                           #detected infectious classes
  I2 <- x[10]
  I3 <- x[11]
  I4 <- x[12]
  
  Iu1 <- x[13]                           #undetected infectious classes
  Iu2 <- x[14]
  Iu3 <- x[15]
  Iu4 <- x[16]
  
  I.detected <-I1+I2+I3+I4
  I.undetected <- Iu1+Iu2+Iu3+Iu4
  I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infectious
  
  H <- x[17]                           #local variable for Isolated
  Ru <- x[18]                          #local variable for undetected recovereds
  
  C <- x[19]    # local variable for notifications
  
  N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
  t <- x[20]                          #get current time
  
  with(                                #use with to simplify code
    as.list(params), 
    {
      gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
      sigmai <- 6*sigma  # multiplier 6 for pseudo stages
      etat <- eta(t)     # case notification rate
      
      ## ATT UPDATE HERE TO INCLUDE ALL ARGUMENTS
      betat <- beta(t,w,beta0,beta.factor)   # time dependent transmissibility, presymptomatic=1 causes this transmissibility to apply to late stage latent cases as well
      
      rates <- as.numeric(c(betat*I.detected/N+betat*c*I.undetected/N+presymptomatic*betat*c*E6/N,                           # movements out of S
                            sigmai, sigmai, sigmai, sigmai, sigmai, sigmai,   # movements out of E
                            gammai, gammai, gammai, gammai,                   # movements out of I (detected)
                            b, b, b, b,                                       # movements out of I (undetected)
                            etat))
      
      states0 <- x[2:(length(x)-1)]
      
      # transition probabilities
      
      p <- matrix(0, nrow=length(rates),ncol=length(states0))                                     # matrix to hold transitions probs
      
      p[1,]  <- c(exp(-rates[1]*dt),1-exp(-rates[1]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # S-> E
      
      p[2,]  <- c(0, exp(-rates[2]*dt), 1-exp(-rates[2]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[3,]  <- c(0, 0, exp(-rates[3]*dt), 1-exp(-rates[3]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[4,]  <- c(0, 0, 0, exp(-rates[4]*dt), 1-exp(-rates[4]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[5,]  <- c(0, 0, 0, 0, exp(-rates[5]*dt), 1-exp(-rates[5]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[6,]  <- c(0, 0, 0, 0, 0, exp(-rates[6]*dt), 1-exp(-rates[6]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[7,]  <- c(0, 0, 0, 0, 0, 0, exp(-rates[7]*dt), (1-exp(-rates[7]*dt))*q(t,w), 0, 0, 0, (1-exp(-rates[7]*dt))*(1-q(t,w)), 0, 0, 0, 0, 0, 0)    # Transitions out of E
      
      p[8,]  <- c(0, 0, 0, 0, 0, 0, 0, exp(-rates[8]*dt), 1-exp(-rates[8]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I (detected)
      p[9,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[9]*dt), 1-exp(-rates[9]*dt), 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I
      p[10,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[10]*dt), 1-exp(-rates[10]*dt), 0, 0, 0, 0, 0, 0, 0)  # Transitions out of I
      p[11,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[11]*dt), 0, 0, 0, 0, 1-exp(-rates[11]*dt), 0, 0)  # Transitions out of I -> H
      
      p[12,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[12]*dt), 1-exp(-rates[12]*dt), 0, 0, 0, 0, 0)    # Transitions out of I (undetected)
      p[13,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[13]*dt), 1-exp(-rates[13]*dt), 0, 0, 0, 0)    # Transitions out of I
      p[14,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[14]*dt), 1-exp(-rates[14]*dt), 0, 0, 0)  # Transitions out of I
      p[15,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[15]*dt), 0, 1-exp(-rates[15]*dt), 0)  # Transitions out of I -> R_u
      
      
      p[16,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[16]*dt), 0, 1-exp(-rates[16]*dt))  # Transitions R_d -> to C (notification)
      
      # update states
      
      states1 <- matrix(0, nrow=length(rates),ncol=length(states0))                                # matrix to hold future states
      
      for(i in 1:length(rates)){
        states1[i,] <- t(rmultinom(1, states0[i], p[i,])) 
      }
      
      states1 <- colSums(states1)
      states1[17] <- states1[17]+Ru  #add formerly Recovered undetected cases
      states1[18] <- states1[18]+C  #add formerly notified cases
      
      return(x <- c(dt, states1, tail(x,1)+dt))
    }
  )
}

model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(round(nstep + 1), length(x))) #set up array to store results
  colnames(output) <- c("time","S",
                        "E1", "E2", "E3", "E4", "E5", "E6",
                        "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4",
                        "H", "Ru", "C", "cum.time") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- as.numeric(onestep(x,params))
  }
  output                                    #return output
}

## UPDATE WITH JOHN
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

# eta <- function(t, w=12) ifelse(t<=w,1/3,1/3)
eta <- function(t) ifelse(t<=55,1/(-0.47*t+27.2),1)  ## UPDATE WITH JOHN
q <- function(t, w=50, q0=0.111, q1=0.98) ifelse(t<=w,q0,q1)  ## UPDATE WITH JOHN

beta <- function(t, w=12, beta0=0.6584, beta.factor=1) ifelse(t<=w,beta0,beta0/beta.factor)

evaluate.model <- function(params, init , nsims, nstep, start, today){
  
  if(is.null(nstep)) nstep <- (as.numeric(today-start)+1+28)/params$dt #run simulation from start to current time plus four weeks
  
  xstart <- c(time=0, unlist(init), cum.time = 0) #initial conditions
  
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  
  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    
  }
  
  return(data)
}

get.curve <- function(x){
  out <- x[[i]]
  dt <- x$cum.time[2]-x$cum.time[1]
  report.times <- seq(1,(length(x$cum.time)), by=1/dt)
  case.notifications <- diff(x$C[report.times])
}