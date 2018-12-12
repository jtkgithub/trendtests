### Functions for calculating trend tests for single process time censored data

# Most of the tests are described in: 
#  Anonymous (2018), A Class of Tests for Trend in Time Censored Recurrent Event Data, 
# See also: 
#  Lawless, Cigsar and Cook (2012), Testing for monotone trend in recurrent event processes,
#  Technometrics 54: 147-158.


## All trend test functions take the following arguments as input: 
# tvec - a vector of times to event
# tau - the censoring time
# sigma - how to calculate the standard deviation, options are:
#          - "s" for the usual sample standard deviation
#          - "c" for the estimator in equation (10) in Anonymous (2018)
#          - "l" for the estimator in equation (11) in Anonymous (2018)
#          - "fixCV" for forcing a specific value of CV
# cv - a specific value of cv used only if sigma="fixCV"

# Library used for calculating p-values for CvM and AD tests
# install.packages("goftest")
library(goftest)


# Function to calculate estimated CV
findCV <- function(tvec,tau,sigma,cv){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  if(sigma=="s"){
    sdest <- sqrt(var(xvec))
    mu <- mean(xvec)
    CV <- sdest/mu
  }
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
  if(sigma=="fixCV"){
    CV <- cv
  }
  return(CV)
}  
  
## Lewis Robinson test
  
# Calculate the test observator
LRtestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n=length(tvec)
  LR <- (sqrt(12)/(CV*tau*sqrt(n)))*(sum(tvec)-n*tau/2)
  return(LR)
}

# Run the test
LRtest <- function(tvec,tau,sigma="s",cv=1){
  LR <- LRtestobs(tvec,tau,sigma,cv)
  pLR <- 2*pnorm(-abs(LR))
  cat("Lewis-Robinson test for trend in time censored processes.\n")
  cat(paste("Test statistic: LR =",round(LR,digits=3),"\n"))
  cat(paste("p-value =",round(pLR,digits=5)),"\n")
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}


## Kolmogorov-Smirnov test

# Function to calculate p-values for the Kolmogorov-Smirnov test
pKS <- function(x,kmax=1000) 
  1-(sqrt(2*pi)/x)*sum(exp(-(2*(1:kmax)-1)^2*pi^2/(8*x^2)))

# Calculate the test observator
KStestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n=length(tvec)
  xvec=diff(c(0,tvec))
  t1vec <- n*c(0,tvec)/tau
  t2vec <- n*c(tvec,tau)/tau
  nvec <- 0:n
  KS <- (1/(CV*sqrt(n)))*max(c(max(abs(nvec-t1vec))),c(max(abs(nvec-t2vec))))
  return(KS)
}  

# Run the test
KStest <- function(tvec,tau,sigma="s",cv=1){
  KS <- KStestobs(tvec,tau,sigma,cv)
  pKS <- pKS(KS)
  cat("Kolmogorov-Smirnov test for trend in time censored processes.\n")
  cat(paste("Test statistic: KS =",round(KS,digits=3),"\n"))
  cat(paste("p-value =",round(pKS,digits=5)),"\n")
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}





## Cramer von-Mises test

# Calculate the test observator
CvMtestobs <- function(tvec,tau,sigma,cv){
  CV <- findCV(tvec,tau,sigma,cv)
  n=length(tvec)
  xvec=diff(c(0,tvec))
  konst <- 1/(CV^2*n)
  indseq <- 0:(n-1)
  sumledd <- sum(indseq^2*xvec/tau)-n*sum(indseq*(tvec^2-c(0,tvec[1:(n-1)])^2)/tau^2)
  sisteledd <- n^2/3+n^2*(tvec[n]^2/tau^2-tvec[n]/tau)
  CV <- konst*(sumledd+sisteledd)
  return(CV)
}  

# Run the test
CvMtest <- function(tvec,tau,sigma="s",cv=1){
  CvM <- CvMtestobs(tvec,tau,sigma,cv)
  pCvM <- 1-pCvM(q=CvM)
  cat("Cramer von-Mises test for trend in time censored processes.\n")
  cat(paste("Test statistic: CvM =",round(CvM,digits=3),"\n"))
  cat(paste("p-value =",round(pCvM,digits=5)),"\n")
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}



## Anderson-Darling test

# Calculate the test observator
ADtestobs <- function(tvec,tau,sigma="s",cv=1){
  CV <- findCV(tvec,tau,sigma,cv)
  n=length(tvec)
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

# Run the test
ADtest <- function(tvec,tau,sigma="s",cv=1){
  AD <- ADtestobs(tvec,tau,sigma,cv)
  pAD <- 1-pAD(q=AD)
  cat("Anderson-Darling test for trend in time censored processes.\n")
  cat(paste("Test statistic: AD =",round(AD,digits=3),"\n"))
  cat(paste("p-value =",round(pAD,digits=5)),"\n")
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}



## Extended Lewis-Robins test 

# Calculate the test observator
ELRtestobs <- function(tvec,tau,sigma,cv,a){
  CV <- findCV(tvec,tau,sigma,cv)
  n=length(tvec)
  konst <- 1/(CV*tau*sqrt(n)*sqrt((1/12)-a^2*(1-a)^2))
  sumledd <- sum(abs(tvec-a*tau))
  sisteledd <- (0.5-a*(1-a))*n*tau  
  ELR <- konst*(sumledd-sisteledd)
  return(ELR)
}

# Run the test
ELRtest <- function(tvec,tau,sigma="s",cv=1, a=0.5){
  ELR <- ELRtestobs(tvec,tau,sigma,cv,a)
  pELR <- 2*pnorm(-abs(ELR))
  cat("Extended Lewis-Robinson test for trend in time censored processes.\n")
  cat(paste("Test statistic: ELR =",round(ELR,digits=3),"\n"))
  cat(paste("p-value =",round(pELR,digits=5)),"\n")
  cat(paste("Using a =",round(a,digits=3),"\n"))
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}





## Linear rank test from Lawless, Cigsar and Cook (2012) 

# Calculate the test observator
LinRanktestobs <- function(tvec){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  Rvec=order(xvec)
  evec=vector(length=n)
  U=0
  VarU=0
  for(j in 1:n){
    rseq=n+1-1:Rvec[j]
    evec[j]=sum(1/rseq)
  }
  ebar=mean(evec)
  jvec=1:n-(n+1)/2
  U=sum(evec*jvec)
  VarU=sum(jvec^2)*sum((evec-ebar)^2/(n-1))
  LinRank <- U/sqrt(VarU)
  list(LinRank=LinRank,V=U,VarV=VarU)
}

# Run the test 
LinRanktest <- function(tvec){
  LinRank <- LinRanktestobs(tvec)$LinRank
  pLinRank <- 2*pnorm(-abs(LinRank))
  cat("Linear rank test.\n")
  cat(paste("Test statistic: LinRank =",round(LinRank,digits=3),"\n"))
  cat(paste("p-value =",round(pLinRank,digits=5)),"\n")
}




## Calculate p-values for all  tests by permutation

# B is the numer of permutations

Permtests <- function(tvec,tau,B=10000,sigma="s", cv=1, a=0.5){
  n=length(tvec)
  xvec=diff(c(0,tvec))
  LR0 <- LRtestobs(tvec,tau,sigma,cv)
  KS0 <- KStestobs(tvec,tau,sigma,cv)
  CvM0 <- CvMtestobs(tvec,tau,sigma,cv)
  AD0 <- ADtestobs(tvec,tau,sigma,cv)
  ELR0 <- ELRtestobs(tvec,tau,sigma,cv,a)
  LinRank0 <- LinRanktestobs(tvec)$LinRank
  LRpvec <- vector(length=B)
  KSpvec <- vector(length=B)  
  CvMpvec <- vector(length=B)  
  ADpvec <- vector(length=B)  
  ELRpvec <- vector(length=B)
  LinRankpvec <- vector(length=B)  
  for(i in 1:B){
    xvecp <- sample(xvec)
    tvecp <- cumsum(xvecp)
    LRpvec[i] <- LRtestobs(tvecp,tau,sigma,cv)
    KSpvec[i] <- KStestobs(tvecp,tau,sigma,cv)
    CvMpvec[i] <- CvMtestobs(tvecp,tau,sigma,cv)
    ADpvec[i] <- ADtestobs(tvecp,tau,sigma,cv)
    ELRpvec[i] <- ELRtestobs(tvecp,tau,sigma,cv,a)
    LinRankpvec[i] <- LinRanktestobs(tvecp)$LinRank
  }
  pLR <- sum(abs(LRpvec)>abs(LR0))/B
  pKS <- sum(KSpvec>KS0)/B
  pCvM <- sum(CvMpvec>CvM0)/B
  pAD <- sum(ADpvec>AD0)/B
  pLinRank <- sum(abs(LinRankpvec)>abs(LinRank0))/B
  pELR <- sum(abs(ELRpvec)>abs(ELR0))/B
  cat("\np-values of the trend tests for time censored processes\n") 
  cat(paste("computed by permutation using B =", B,"repetitions\n\n"))
  cat("Lewis-Robinson test.\n")
  cat(paste("Test statistic: LR =",round(LR0,digits=3),"\n"))
  cat(paste("p-value =",round(pLR,digits=5)),"\n\n")
  cat("Kolmogorov-Smirnov test.\n")
  cat(paste("Test statistic: KS =",round(KS0,digits=3),"\n"))
  cat(paste("p-value =",round(pKS,digits=5)),"\n\n")
  cat("Cramer von-Mises test.\n")
  cat(paste("Test statistic: CvM =",round(CvM0,digits=3),"\n"))
  cat(paste("p-value =",round(pCvM,digits=5)),"\n\n")
  cat("Anderson-Darling test.\n")
  cat(paste("Test statistic: AD =",round(AD0,digits=3),"\n"))
  cat(paste("p-value =",round(pAD,digits=5)),"\n\n")
  cat("Linear rank test.\n")
  cat(paste("Test statistic: LinRank =",round(LinRank0,digits=3),"\n"))
  cat(paste("p-value =",round(pLinRank,digits=5)),"\n\n")
  cat("Extended Lewis-Robinson test.\n")
  cat(paste("Test statistic: ELR =",round(ELR0,digits=3),"\n"))
  cat(paste("p-value =",round(pELR,digits=5)),"\n")
  cat(paste("Using a =",a,"\n\n"))
  CV <- findCV(tvec,tau,sigma,cv)
  cat(paste("CV =", round(CV,digits=3) ,"\n"))
}

