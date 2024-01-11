# Code to simulate the asymptotic null distribution of the 
# ICvM, IKS and SELR1 tests

library(e1071)

# ELR(a), not normalised
ELR0a <- function(x,a){
  n <- length(x)
  delta <- 1/n
  an <- ceiling(a*n)
  int1 <- ifelse(an>0,delta*sum(x[1:an]),0)
  int2 <- ifelse(an<n,delta*sum(x[(an+1):n]),0)
  ELRa <- int1-int2
  return(list(ELRa=ELRa,int1=int1,int2=int2))
}

# ELR(a), normalised
ELR1a <- function(x,a){
  n <- length(x)
  delta <- 1/n
  an <- ceiling(a*n)
  int1 <- ifelse(an>0,delta*sum(x[1:an]),0)
  int2 <- ifelse(an<n,delta*sum(x[(an+1):n]),0)
  ELRa <- int1-int2
  ELRa <- ELRa*(1/sqrt(1/12-(a^2)*(1-a)^2))
  return(list(ELRa=ELRa,int1=int1,int2=int2))
}

# Function to calculate: 
# sup_a ELR(a) 
# sup_a |ELR(a)| 
# int_a ELR(a)
# int_a |ELR(a)|
# Not normalised
AELR0s <- function(x,aseq){
  n <- length(aseq)
  delta <- 1/n
  ELRs <- numeric(n)
  ELR1s <- numeric(n)
  ELR2s <- numeric(n)
  for(i in 1:n){
    ELRresa <- ELR0a(x,aseq[i])
    ELRs[i] <- ELRresa$ELRa
    ELR1s[i] <- ELRresa$int1
    ELR2s[i] <- ELRresa$int2
  }
  AELRsup <- max(ELRs)
  AELRsup2 <- max(abs(ELRs))
  amax <- aseq[which.max(abs(ELRs))]
  AELRint1 <- delta*sum(ELRs)
  AELRint2 <- delta*sum(abs(ELRs))
  INTLR1 <- delta*sum(ELR1s)
  INTLR2 <- delta*sum(ELR2s)
  ICvM <- delta*sum(ELR1s^2)
  IKS <- max(abs(ELR1s))
  return(list(AELR0sup1=AELRsup,AELR0sup2=AELRsup2,amax=amax,AELR0int1=AELRint1,AELR0int2=AELRint2,
              INT0LR1=INTLR1, INT0LR2=INTLR2, INT0CvM=ICvM, INT0KS=IKS))
}


# Function to calculate: 
# sup_a ELR(a) 
# sup_a |ELR(a)| 
# int_a ELR(a)
# int_a |ELR(a)|
# Normalised
AELR1s <- function(x,aseq){
  n <- length(aseq)
  delta <- 1/n
  ELRs <- numeric(n)
  ELR1s <- numeric(n)
  ELR2s <- numeric(n)
  for(i in 1:n){
    ELRresa <- ELR1a(x,aseq[i])
    ELRs[i] <- ELRresa$ELRa
    ELR1s[i] <- ELRresa$int1
    ELR2s[i] <- ELRresa$int2
  }
  AELRsup <- max(ELRs)
  AELRsup2 <- max(abs(ELRs))
  amax <- aseq[which.max(abs(ELRs))]
  AELRint1 <- delta*sum(ELRs)
  AELRint2 <- delta*sum(abs(ELRs))
  INTLR1 <- delta*sum(ELR1s)
  INTLR2 <- delta*sum(ELR2s)
  ICvM <- delta*sum(ELR1s^2)
  IKS <- max(abs(ELR1s))
  return(list(AELR1sup1=AELRsup,AELR1sup2=AELRsup2,amax=amax,AELR1int1=AELRint1,AELR1int2=AELRint2,
              INT1LR1=INTLR1, INT1LR2=INTLR2, INT1CvM=ICvM, INT1KS=IKS))
}


# Function to simulate distributions
simAELRdists <- function(nBBs=100,deltaBBs=0.01,aseq=seq(0,1,by=0.01)){
  # nBBs = number of Brownian Bridges to be simulated
  # deltaBBs = step size in the simulated BBs 
  # aseq = the sequence of a's to search through
  amaxs <- numeric(nBBs)
  AELR0sup1s <- numeric(nBBs)
  AELR0sup2s <- numeric(nBBs)
  AELR0int1s <- numeric(nBBs)
  AELR0int2s <- numeric(nBBs)
  INT0LR1s <- numeric(nBBs)
  INT0LR2s <- numeric(nBBs)
  INT0CvMs <- numeric(nBBs)
  INT0KSs <- numeric(nBBs)
  AELR1sup1s <- numeric(nBBs)
  AELR1sup2s <- numeric(nBBs)
  AELR1int1s <- numeric(nBBs)
  AELR1int2s <- numeric(nBBs)
  INT1LR1s <- numeric(nBBs)
  INT1LR2s <- numeric(nBBs)
  INT1CvMs <- numeric(nBBs)
  INT1KSs <- numeric(nBBs)
  for(i in 1:nBBs){
    x <- as.vector(rbridge(frequency = 1/deltaBBs))
    simAELR0s <- AELR0s(x,aseq)
    amaxs[i] <- simAELR0s$amax 
    AELR0sup1s[i] <- simAELR0s$AELR0sup1 
    AELR0sup2s[i] <- simAELR0s$AELR0sup2 
    AELR0int1s[i] <- simAELR0s$AELR0int1 
    AELR0int2s[i] <- simAELR0s$AELR0int2 
    INT0LR1s[i] <- simAELR0s$INT0LR1 
    INT0LR2s[i] <- simAELR0s$INT0LR2 
    INT0CvMs[i] <- simAELR0s$INT0CvM 
    INT0KSs[i] <- simAELR0s$INT0KS 
    simAELR1s <- AELR1s(x,aseq)
    AELR1sup1s[i] <- simAELR1s$AELR1sup1 
    AELR1sup2s[i] <- simAELR1s$AELR1sup2 
    AELR1int1s[i] <- simAELR1s$AELR1int1 
    AELR1int2s[i] <- simAELR1s$AELR1int2 
    INT1LR1s[i] <- simAELR1s$INT1LR1 
    INT1LR2s[i] <- simAELR1s$INT1LR2 
    INT1CvMs[i] <- simAELR1s$INT1CvM 
    INT1KSs[i] <- simAELR1s$INT1KS 
  }
  return(list(asups=amaxs,
              AELR0sup1s=AELR0sup1s,AELR0sup2s=AELR0sup2s,AELR0int1s=AELR0int1s,AELR0int2s=AELR0int2s,
              INT0LR1s=INT0LR1s, INT0LR2s=INT0LR2s, INT0CvMs=INT0CvMs, INT0KSs=INT0KSs,
              AELR1sup1s=AELR1sup1s,AELR1sup2s=AELR1sup2s,AELR1int1s=AELR1int1s,AELR1int2s=AELR1int2s,
              INT1LR1s=INT1LR1s, INT1LR2s=INT1LR2s, INT1CvMs=INT1CvMs, INT1KSs=INT1KSs
  ))
}

AELRsim <- simAELRdists(nBBs=10000,deltaBBs=0.001,aseq=seq(0,1,by=0.001))
AELR0sup1dist <- AELRsim$AELR0sup1s
AELR0sup2dist <- AELRsim$AELR0sup2s
AELR0int1dist <- AELRsim$AELR0int1s
AELR0int2dist <- AELRsim$AELR0int2s
INT0LR1dist <- AELRsim$INT0LR1s
INT0LR2dist <- AELRsim$INT0LR2s
INT0CvMdist <- AELRsim$INT0CvMs
INT0KSdist <- AELRsim$INT0KSs
AELR1sup1dist <- AELRsim$AELR1sup1s
AELR1sup2dist <- AELRsim$AELR1sup2s
AELR1int1dist <- AELRsim$AELR1int1s
AELR1int2dist <- AELRsim$AELR1int2s
INT1LR1dist <- AELRsim$INT1LR1s
INT1LR2dist <- AELRsim$INT1LR2s
INT1CvMdist <- AELRsim$INT1CvMs
INT1KSdist <- AELRsim$INT1KSs


# Plotting some of the asymptotic distributions 
par(mfrow=c(3,1))
hist(AELR1sup1dist, probability = TRUE,nclass=50, col="burlywood")
abline(v=quantile(AELR1sup1dist, probs=c(0.95)), col="blue", lty=2)
abline(v=quantile(AELR1sup1dist, probs=c(0.99)), col="red", lty=2)
quantile(AELR1sup1dist, probs=c(0.5,0.9,0.95,0.99))


hist(INT1CvMdist, probability = TRUE,nclass=50, col="burlywood")
abline(v=quantile(INT1CvMdist, probs=c(0.95)), col="blue", lty=2)
abline(v=quantile(INT1CvMdist, probs=c(0.99)), col="red", lty=2)
quantile(INT1CvMdist, probs=c(0.5,0.9,0.95,0.99))


hist(INT1KSdist, probability = TRUE,nclass=50, col="burlywood")
abline(v=quantile(INT1KSdist, probs=c(0.95)), col="blue", lty=2)
abline(v=quantile(INT1KSdist, probs=c(0.99)), col="red", lty=2)
quantile(INT1KSdist, probs=c(0.5,0.9,0.95,0.99))
