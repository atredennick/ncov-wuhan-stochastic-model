

## ----notifications-------------------------------------------------------
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

plot(data.province$date, data.province$Hubei, type='h', lwd=5, col='lightskyblue',
     lend='butt', xlab='Date', ylab='Case notifications', main='Case notifications (Hubei)')


## ----onestep-------------------------------------------------------------
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
  I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infections
  
  H <- x[17]                           #local variable for hospitalized
  Ru <- x[18]                          #local variable for undetected recovereds
  
  C <- x[19]    # local variable for notifications
 
  N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
  t <- x[20]                          #get current time
  
  
  with(                                #use with to simplify code
       as.list(params), 
       {
         gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
         sigmai <- 6*sigma  # multiplier 6 for pseudo stages
         etat <- eta(t)
        # betat <- beta(beta0, t)

         rates <- as.numeric(c(beta0*I.detected/N+beta0*c*I.undetected/N,                           # movements out of S
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


## ----model---------------------------------------------------------------
model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,length(x)))         #set up array to store results
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


## ----evaluate------------------------------------------------------------
evaluate.model <- function(params=list(beta0=0.657, sigma=1/6.4, z=45, b=0.143, a0=0.0446, w=40, c=1,  dt=0.05),
                                       init = list(S=59002000, E1=0, E2=0, E3=0, E4=0, E5=6, E6=0,
                                                   I1 = 1, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                                   H=0, Ru=0, C=0),
                           nsims=2, nstep=NULL, start=as.Date("2019-12-01"),today=Sys.Date()){
  
  #popn size of Wuhan 11081000
  #popn size of Hubei 59002000
  if(is.null(nstep)) nstep <- (as.numeric(today-start)+1+14)/params$dt #run simulation from start to current time plus two weeks

  xstart <- c(time=0, unlist(init), cum.time = 0) #initial conditions

  data <- vector(mode='list',length=nsims) #initialize list to store the output

  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    
  }

  return(data)
}


## ----plot-function-------------------------------------------------------
plot.model <- function(data, log='y'){
  
  # process data
  nsims <- length(data)
 
  for(i in 1:nsims) data[[i]]$I <- data[[i]]$I1 + data[[i]]$I2 + data[[i]]$I3 +
      data[[i]]$I4
  for(i in 1:nsims) data[[i]]$Iu <- data[[i]]$Iu1 + data[[i]]$Iu2 + data[[i]]$Iu3 +
      data[[i]]$Iu4
  for(i in 1:nsims) data[[i]]$E <- data[[i]]$E1 + data[[i]]$E2 + data[[i]]$E3 +
      data[[i]]$E4 + data[[i]]$E5 + data[[i]]$E6
  
  max.time<-data[[1]]$cum.time[max(which(data[[1]]$I>0))] #maximum time in first simulation
  max.y<-max(data[[1]]$C)       #find max total confirmed cases for plotting range
  
  # calculate means
  m1 <- m2 <- m3 <- m4 <- m5 <- matrix(nrow=length(data[[1]]$I), ncol=nsims)
  for(i in 1:nsims){
    m1[,i] <- data[[i]]$E
    m2[,i] <- data[[i]]$I+data[[i]]$Iu
   # m3[,i] <- data[[i]]$Iu
    m4[,i] <- data[[i]]$H
    m5[,i] <- data[[i]]$C
  }
  E.mean <- rowMeans(m1)
  I.mean <- rowMeans(m2)
 # Iu.mean <- rowMeans(m3)
  H.mean <- rowMeans(m4)
  C.mean <- rowMeans(m5)
  
  # colors
  E.col <- rgb(0,1,0,.25)
  I.col <- rgb(1,0,0,.25)
  Iu.col <- rgb(0.5, 0.5, 0, 0.25)
  H.col <- rgb(0,0,1,.25)
  C.col <- rgb(0,0,0,.25)
  E.mean.col <- rgb(0,1,0,1)
  I.mean.col <- rgb(1,0,0,1)
  Iu.mean.col <- rgb(0.5,0.5,0,1)
  H.mean.col <- rgb(0,0,1,1)
  C.mean.col <- rgb(0,0,0,1)
  
  #set up plot
  plot(I~cum.time,data=data[[1]],xlab='',ylab='Cases',col=1,
       xlim=c(0,max.time),ylim=c(1,max.y), type='n', lty=1, log=log, axes=FALSE) # set up plot
  
  # add data to plot
  day <- data.province$date - start
  lines(day, cumsum(data.province$Hubei), type='h', col='lightskyblue', lwd=3, lend='butt' )
  
  # plot spaghetti
  lines(E~cum.time,data=data[[1]], col=E.col, lty=1)
  lines(I+Iu~cum.time,data=data[[1]], col=I.col, lty=1)
 # lines(Iu~cum.time,data=data[[1]], col=Iu.col, lty=1)
  lines(H~cum.time,data=data[[1]], col=H.col, lty=1)
  lines(C~cum.time,data=data[[1]], col=C.col, lty=1, lwd=3)

 
  axis(1, at=seq(0,max.time,5), labels=format(start+seq(0,max.time,5), format= '%b %d'))
  axis(2)
  box()

  if(nsims > 1){
  for (k in 2:min(100,nsims)) {              #add multiple epidemics to plot
    lines(E~cum.time, data=data[[k]], col=E.col, type='l', lty=1)
    lines(I+Iu~cum.time, data=data[[k]], col=I.col, type='l', lty=1)
 #   lines(Iu~cum.time, data=data[[k]], col=Iu.col, type='l', lty=1)
    lines(H~cum.time, data=data[[k]], col=H.col, type='l', lty=1)
    lines(C~cum.time, data=data[[k]], col=C.col, type='l', lty=1, lwd=3)
  }
  
  # plot means
  lines(E.mean~cum.time, data=data[[k]], col=E.mean.col, lty=1)
  lines(I.mean~cum.time, data=data[[k]], col=I.mean.col, lty=1)  
#  lines(Iu.mean~cum.time, data=data[[k]], col=Iu.mean.col, lty=1)
  lines(H.mean~cum.time, data=data[[k]], col=H.mean.col, lty=1)
  lines(C.mean~cum.time, data=data[[k]], col=C.mean.col, lty=1)
  }

   legend('topleft', lty=c(1,1,1,1,1,1), lwd=c(1,1,1,1,3,3), bty='n', cex=0.75,
         col=c(E.col, I.col, H.col, C.col, 'lightskyblue'),
         legend=c('Latent cases', 'Infectious cases', 'Hospitalized', 
                  'Cumulative reported cases (Model)', 'Cumulative reported cases (Data)'))
}


## ----plot-function2------------------------------------------------------
plot.model2 <- function(data, log='y'){
  
  # process data
  nsims <- length(data)
 
  for(i in 1:nsims) data[[i]]$I <- data[[i]]$I1 + data[[i]]$I2 + data[[i]]$I3 +
      data[[i]]$I4
  for(i in 1:nsims) data[[i]]$Iu <- data[[i]]$Iu1 + data[[i]]$Iu2 + data[[i]]$Iu3 +
      data[[i]]$Iu4
  for(i in 1:nsims) data[[i]]$E <- data[[i]]$E1 + data[[i]]$E2 + data[[i]]$E3 +
      data[[i]]$E4 + data[[i]]$E5 + data[[i]]$E6
  
  max.time<-data[[1]]$cum.time[max(which(data[[1]]$I>0))] #maximum time in first simulation
  #max.y<-max(data[[1]]$Q)       #find max total confirmed cases for plotting range
  max.y <- max(unlist(lapply(data, FUN=function(x) max(x$C))))
  
  # calculate means
  m1 <- m2 <- m3 <- m4 <- m5 <- matrix(nrow=length(data[[1]]$I), ncol=nsims)
  for(i in 1:nsims){
    m1[,i] <- data[[i]]$E
    m2[,i] <- data[[i]]$I
    m3[,i] <- data[[i]]$Iu
    m4[,i] <- data[[i]]$H
    m5[,i] <- data[[i]]$C
  }
  
  observed <- m5
  unobserved <- m1 + m2 + m3 + m4
  
  obs.mean <- rowMeans(observed)
  unobs.mean <- rowMeans(unobserved)
  cum.time <- data[[1]]$cum.time
  
  # colors
  E.col <- rgb(0,1,0,.25)
  I.col <- rgb(1,0,0,.25)
  Iu.col <- rgb(0.5, 0.5, 0, 0.25)
  H.col <- rgb(0,0,1,.25)
  C.col <- rgb(0,0,0,.25)
  E.mean.col <- rgb(0,1,0,1)
  I.mean.col <- rgb(1,0,0,1)
  Iu.mean.col <- rgb(0.5,0.5,0,1)
  H.mean.col <- rgb(0,0,1,1)
  C.mean.col <- rgb(0,0,0,1)
  
  #set up plot
  plot(obs.mean~cum.time, xlab='',ylab='Cases',col=1,
       xlim=c(0,max.time),ylim=c(1,max.y), type='n', lty=1, log=log, axes=FALSE) # set up plot
  
  # add data to plot
  day <- data.province$date - start
  lines(day, cumsum(data.province$Hubei), type='h', col='lightskyblue', lwd=3, lend='butt' )
  
  # plot spaghetti
  lines(unobserved[,1]~cum.time, col=E.col, lty=1)
  lines(observed[,1]~cum.time, col=C.col, lty=1, lwd=1)
 
  axis(1, at=seq(0,max.time,5), labels=format(start+seq(0,max.time,5), format= '%b %d'))
  axis(2)
  box()

  if(nsims > 1){
  for (k in 2:min(100,nsims)) {              #add multiple epidemics to plot
    lines(unobserved[,k]~cum.time,  col=E.col, type='l', lty=1)
    lines(observed[,k]~cum.time, col=C.col, type='l', lty=1, lwd=1)
  }
  
  # plot means
  lines(unobs.mean~cum.time, col=E.mean.col, lty=1)
  lines(obs.mean~cum.time, col=C.mean.col, lty=1)
  }

   legend('topleft', lty=c(1,1,1), lwd=c(1,1,3), bty='n', cex=0.75,
         col=c(E.col, C.col, 'lightskyblue'),
         legend=c('Unobserved cases', 'Cumulative reported cases (Model)', 'Cumulative reported cases (Data)'))
}


## ----plot-summary--------------------------------------------------------
plot.summary <- function(data){
   # process data
  nsims <- length(data)
 
  for(i in 1:nsims) data[[i]]$I <- data[[i]]$I1 + data[[i]]$I2 + data[[i]]$I3 +
      data[[i]]$I4
  for(i in 1:nsims) data[[i]]$Iu <- data[[i]]$Iu1 + data[[i]]$Iu2 + data[[i]]$Iu3 +
      data[[i]]$Iu4
  for(i in 1:nsims) data[[i]]$E <- data[[i]]$E1 + data[[i]]$E2 + data[[i]]$E3 +
      data[[i]]$E4 + data[[i]]$E5 + data[[i]]$E6
  
  max.time<-data[[1]]$cum.time[max(which(data[[1]]$I>0))] #maximum time in first simulation
  max.y<-max(data[[1]]$I)       #find max total confirmed cases for plotting range
  
  # calculate means
  m1 <- m2 <- m3 <- m4 <- m5 <- matrix(nrow=length(data[[1]]$I), ncol=nsims)
  for(i in 1:nsims){
    m1[,i] <- data[[i]]$E
    m2[,i] <- data[[i]]$I
    m3[,i] <- data[[i]]$Iu
    m4[,i] <- data[[i]]$H
    m5[,i] <- data[[i]]$C
  }
  E.mean <- rowMeans(m1)
  I.mean <- rowMeans(m2)
  Iu.mean <- rowMeans(m3)
  H.mean <- rowMeans(m4)
  C.mean <- rowMeans(m5)
  
  # colors
  E.col <- rgb(0,1,0,.25)
  I.col <- rgb(1,0,0,.25)
  Iu.col <- rgb(0.5, 0.5, 0, 0.25)
  H.col <- rgb(0,0,1,.25)
  C.col <- rgb(0,0,0,.25)
  E.mean.col <- rgb(0,1,0,1)
  I.mean.col <- rgb(1,0,0,1)
  Iu.mean.col <- rgb(0.5,0.5,0,1)
  H.mean.col <- rgb(0,0,1,1)
  C.mean.col <- rgb(0,0,0,1)
  
  #set up plot
  j <- max(which(data[[1]]$time < today.day))
  
  h1 <- hist(log10(m1[j,]+m2[j,]+m3[j,]+m4[j,]), plot=FALSE)  # total unnotifed cases in the community is E+I+Iu+H
  h2 <- hist(log10(m5[j,]), plot=FALSE)         # total case notifications
  
  xlim <- c(floor(min(c(h1$mids, h2$mids))),ceiling(max(c(h1$mids, h2$mids))))
  ylim <- c(0, ceiling(max(c(h1$density, h2$density))))
  
  plot(xlim, ylim, type='n', axes=FALSE, xlab='Outbreak size', ylab='Probability density', main=paste('Total observed and unobserved cases as of', today))
  lines(h1$mids, h1$density, type='l', lwd=3, lend='butt', col=rgb(0.5,0.5,0, 0.5))
  lines(h2$mids, h2$density, type='l', lwd=3, lend='butt', col=C.col)
  legend('topleft', lty=1, col = c(rgb(0.5,0.5,0, 0.5), C.col), legend=c('Unobserved cases', 'Case notifications'), bty='n')
  axis(1, seq(min(xlim), max(xlim), by=1), labels=formatC(10^seq(min(xlim), max(xlim), by=1), format='d', big.mark=","))
  box()
}



## ----gamma---------------------------------------------------------------
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


## ----recovery-rate-plot--------------------------------------------------
t <- seq(0, 80)
g <- gamma(t=t)
par(mfrow=c(1,2))
plot(t,g, type='l', xlab='Time (days)', ylab='Recovery rate', ylim=c(0,2.5), main=today)
abline(h=1/7, lty=2, col='red')
abline(h=0.33, lty=2, col='blue')
abline(v=today.day, lty=2, col='darkgrey')
legend('topleft', col=c('black','red','blue','darkgrey'), lty=c(1,2,2,2), bty='n',
       legend=c('Recovery rate','Baseline recovery','Transmissibility','Today'))
plot(t, 1/g, type='l', xlab='Time (days)', ylab='Infectious period (days)', main=today)
abline(h=7, lty=2, col='red')
abline(h=1/0.33, lty=2, col='blue')
abline(v=today.day, lty=2, col='darkgrey')


## ----eta-----------------------------------------------------------------
eta <- function(t) ifelse(t<=55,1/(-0.47*t+27.2),1)
plot(t, eta(t), type='l', xlab='', ylab='Case notification rate')
abline(v=today.day, lty=2, col='darkgrey')
text(today.day, 0.8, labels = paste('Today:',today), pos=4, col='darkgrey', cex=0.8)


## ----q-function----------------------------------------------------------
fever.clinic <- as.Date('2020-01-09')
w <- fever.clinic-start+1
q <- function(t, w=40, q0=0.11, q1=0.98) ifelse(t<=w,q0,q1)


## ----sim1----------------------------------------------------------------
start <- as.Date('12/01/2019',format='%m/%d/%Y')
set.seed(1292020)                #set seed
out1 <- evaluate.model(init = list(S=59002000,
                                 E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start)
plot.model(out1, log='y')
plot.model2(out1, log='y')
#plot.summary(out1)


## ----sim2----------------------------------------------------------------
start <- as.Date('11/15/2019',format='%m/%d/%Y')
set.seed(1292020)                #set seed
out2 <- evaluate.model(params=list(beta0=0.657, sigma=1/6.4, z=56, b=0.143, a0=0.0446, w=56, c=1, dt=0.05),
                       init = list(S=59002000,
                                 E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start, today=today)
plot.model2(out2)
#plot.summary(out3)


## ----sim3----------------------------------------------------------------
start <- as.Date('12/15/2019',format='%m/%d/%Y')
set.seed(1252020)                #set seed
out3 <- evaluate.model(params=list(beta0=0.657, sigma=1/6.4, z=30, b=0.143, a0=0.0446, w=26, c=1, dt=0.05),
                       init = list(S=59002000,
                                 E1=5, E2=5, E3=5, E4=5, E5=5, E6=5,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start, today=today)
plot.model2(out3)


## ----sim4----------------------------------------------------------------
start <- as.Date('12/01/2019',format='%m/%d/%Y')
set.seed(1252020)                #set seed
out4 <- evaluate.model(params=list(beta0=0.657, sigma=1/6.4, z=3000, b=0.143, a0=0.0446, w=40, c=1, dt=0.05),
                       init = list(S=59002000,
                                 E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start, today=today)
plot.model2(out4)


## ----sim5----------------------------------------------------------------
start <- as.Date('12/01/2019',format='%m/%d/%Y')
set.seed(1292020)                #set seed
out5 <- evaluate.model(params=list(beta0=0.657, sigma=1/6.4, z=45, b=0.143, a0=0.0446, w=40, c=0.52,  dt=0.05), init = list(S=59002000,
                                 E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start)
#plot.model(out5, log='y')
plot.model2(out5, log='y')


## ----onestep2------------------------------------------------------------
#this version of the model allows latent infections to transmit
onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  
  S <- x[2]                            #local variable for susceptibles
  
  E1 <- x[3]                           #exposed classes
  E2 <- x[4]
  E3 <- x[5]
  E4 <- x[6]
  E5 <- x[7]
  E6 <- x[8]
  
  E <- E1+E2+E3+E4+E5+E6
  
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
  I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infections
  
  H <- x[17]                           #local variable for hospitalized
  Ru <- x[18]                          #local variable for undetected recovereds
  
  C <- x[19]    # local variable for notifications
 
  N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
  t <- x[20]                          #get current time
  
  
  with(                                #use with to simplify code
       as.list(params), 
       {
         gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
         sigmai <- 6*sigma  # multiplier 6 for pseudo stages
         etat <- eta(t)
        # betat <- beta(beta0, t)

         rates <- as.numeric(c(beta0*I.detected/N+beta0*c*I.undetected/N+beta0*c*E/N,                           # movements out of S
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


## ----sim6----------------------------------------------------------------
start <- as.Date('12/01/2019',format='%m/%d/%Y')
set.seed(1292020)                #set seed
out6 <- evaluate.model(init = list(S=59002000,
                                 E1=0, E2=0, E3=0, E4=0, E5=0, E6=0,
                                 I1=1, I2=0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0, 
                                 R=0, Ru=0, Q=0),
                     nsims=25, nstep=NULL, start=start)
#plot.model(out5, log='y')
plot.model2(out6, log='y')

