# R-code for the trend tests considered in the manuscript
# New Tests for Trend in Time Censored Recurrent Event Data
# by Bo Henry Lindqvist and Jan Terje Kvaløy

# LR test for time truncated data
LRtid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  LR <- (sqrt(12)/(CV*tau*sqrt(n)))*(sum(tvec)-n*tau/2)
  return(LR)
}


# Extended LR for time truncation, not normalised
ELRtid0 <- function(tvec,tau,sigma="s",cv=1,a=0.5){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  #  konst <- 1/(CV*tau*sqrt(n)*sqrt((1/12)-a^2*(1-a)^2))
  konst <- 1/(CV*tau*sqrt(n))
  sumledd <- sum(abs(tvec-a*tau))
  sisteledd <- (0.5-a*(1-a))*n*tau  
  ELR <- konst*(sumledd-sisteledd)
  return(ELR)
}


# Extended LR for time truncation, normalised
ELRtid1 <- function(tvec,tau,sigma="s",cv=1,a=0.5){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  konst <- 1/(CV*tau*sqrt(n)*sqrt((1/12)-a^2*(1-a)^2))
  sumledd <- sum(abs(tvec-a*tau))
  sisteledd <- (0.5-a*(1-a))*n*tau  
  ELR <- konst*(sumledd-sisteledd)
  return(ELR)
}



# Adaptive Extended LR for time truncation based on each of: 
# sup_a |ELR(a)| 
# sup_a ELR(a) 
# int_a ELR(a)
# int_a |ELR(a)|
# Not normalised
AELR0tid <- function(tvec,tau,sigma="s",cv=1, agrid=seq(0,1,by=0.01)){
  n <- length(agrid)
  ELRs <- numeric(n)
  for(i in 1:n)
    ELRs[i] <- ELRtid0(tvec,tau,sigma="s",cv=1, a=agrid[i])
  AELRsup1 <- max(ELRs) # sup test
  AELRsup2 <- max(abs(ELRs)) # sup test2
  asup <- agrid[which.max(abs(ELRs))]
  delta <- 1/n 
  AELRint1 <- delta*sum(ELRs) # int1 test
  AELRint2 <- delta*sum(abs(ELRs)) # int2 test
  return(list(AELR0sup1=AELRsup1, AELR0sup2=AELRsup2, 
              AELR0int1=AELRint1, AELR0int2=AELRint2, asup=asup))
}


# Adaptive Extended LR for time truncation based on each of: 
# sup_a |ELR(a)| 
# sup_a ELR(a) 
# int_a ELR(a)
# int_a |ELR(a)|
# Normalised
AELR1tid <- function(tvec,tau,sigma="s",cv=1, agrid=seq(0,1,by=0.01)){
  n <- length(agrid)
  ELRs <- numeric(n)
  for(i in 1:n)
    ELRs[i] <- ELRtid1(tvec,tau,sigma="s",cv=1, a=agrid[i])
  AELRsup1 <- max(ELRs) # sup test
  AELRsup2 <- max(abs(ELRs)) # sup test2
  asup <- agrid[which.max(abs(ELRs))]
  delta <- 1/n 
  AELRint1 <- delta*sum(ELRs) # int1 test
  AELRint2 <- delta*sum(abs(ELRs)) # int2 test
  return(list(AELR1sup1=AELRsup1, AELR1sup2=AELRsup2, 
              AELR1int1=AELRint1, AELR1int2=AELRint2, asup=asup))
}



# ILR1
INTLR1tid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  konst <- sqrt(45)/(CV*sqrt(n))
  sumledd <- sum(tvec)/tau-sum(tvec^2)/(2*tau^2)
  sisteledd <- n/3  
  INT1 <- konst*(sumledd-sisteledd)
  return(INT1)
}


# ILR2
INTLR2tid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  konst <- sqrt(45)/(CV*sqrt(n))
  sumledd <- sum(tvec^2)/(2*tau^2)
  sisteledd <- n/6  
  INT2 <- konst*(sumledd-sisteledd)
  return(INT2)
}

# Function to calculate int_0^a V_tg
calcintVtg <- function(tvec,tau,sigma="s",cv=1,a)
{
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }  
  konst <- 1/(CV*tau*sqrt(n))
  sumledd <- sum(a*tau-pmin(tvec,a*tau))
  sisteledd <- 0.5*a^2*tau*n
  intVtg <- konst*(sumledd-sisteledd)
  return(intVtg)
}

# Integrate (Vtg(a))^2  and  sup|Vtg(a)|
INTCvMKS1tid <- function(tvec,tau,sigma="s",cv=1, agrid=seq(0,1,by=0.01)){
  n <- length(agrid)
  intVtgs <- numeric(n)
  for(i in 1:n)
    intVtgs[i] <- calcintVtg(tvec,tau,sigma="s",cv=1, a=agrid[i])
  delta <- 1/n 
  ICvM <- delta*sum(intVtgs^2)
  IKS <- max(abs(intVtgs))
  return(list(ICvM=ICvM, IKS=IKS))
}



# KS for time truncation
KStid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  t1vec <- n*c(0,tvec)/tau
  t2vec <- n*c(tvec,tau)/tau
  nvec <- 0:n
  KS <- (1/(CV*sqrt(n)))*max(c(max(abs(nvec-t1vec))),c(max(abs(nvec-t2vec))))
  return(KS)
}



# CV for time truncation
CVtid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  konst <- 1/(CV^2*n)
  indseq <- 0:(n-1)
  sumledd <- sum(indseq^2*xvec/tau)-n*sum(indseq*(tvec^2-c(0,tvec[1:(n-1)])^2)/tau^2)
  sisteledd <- n^2/3+n^2*(tvec[n]^2/tau^2-tvec[n]/tau)
  CV <- konst*(sumledd+sisteledd)
  return(CV)
}


# AD for time truncation
ADtid <- function(tvec,tau,sigma="s",cv=1){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="c"){
    xvec2 <- c(xvec,tau-tvec[n])
    mu <- tau/n
    sdest <- sqrt(sum(xvec2^2)/n-mu^2)
    CV <- sdest/mu
  }
  if(sigma=="l"){
    xdiff <- xvec[2:n]-xvec[1:(n-1)]
    sdest <- sqrt(sum(xdiff^2)/(2*(n-1)))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
  if(sigma=="fixCV"){
    CV <- cv
  }
  if(n==1)
    return((log(tau^2/((tau-tvec[1])*tvec[1]))-1)/CV^2)
  konst <- 1/(CV^2*n)
  iseq <- 1:(n-1)
  iseq2 <- iseq^2
  Niseq <- (n-1):1
  Niseq2 <- Niseq^2
  tip <- tvec[2:n]
  ti <- tvec[1:(n-1)]
  lntfrac <- log(tip/ti)
  lntautfrac <- log((tau-ti)/(tau-tip))
  sumledd <- sum(Niseq2*lntautfrac+iseq2*lntfrac)
  sisteledd <- n^2*(log(tau/(tau-tvec[1]))+log(tau/tvec[n])-1)
  AD <- konst*(sumledd+sisteledd)
  return(AD)
}