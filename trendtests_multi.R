### Functions for calculating trend tests for multiple process time censored data

# Most of the tests are described in: 
#  Kvaløy and Lindqvist (2018), A Class of Tests for Trend in Time Censored Recurrent Event Data, 
#  arXiv:1802.08339 [stat.ME], https://arxiv.org/abs/1802.08339
# See also: 
#  Lawless, Cigsar and Cook (2012), Testing for monotone trend in recurrent event processes,
#  Technometrics 54: 147-158.


## All trend test functions take the following arguments as input: 
# tlist - a data list with the following elements
#           - tvec, a list of the time to event vectors for each process
#           - tauvec, a vector of the censoring times for each process
#           - m, the number of processes
# weights - specification of weight with the following options:
#            - "equal" for equal weights
#            - "tau" for weights proportional to the length of each observation interval 
#            - "sqrtN" for weights proportional to sqrt(number of observations in process)
#            - "sqrtNCV" for weights proportional to sqrt(number of observations in process)/CV
#            - "CVntau" for weights proportional to cv*tau*sqrtN 
# sigma - how to calculate the standard deviation, options are:
#          - "s" for the usual sample standard deviation
#          - "c" for the estimator in equation (10) in Kvaløy & Lindqvist (2018)
#          - "l" for the estimator in equation (11) in Kvaløy & Lindqvist (2018)
#          - "fixCV" for forcing a specific value of CV
# cv - a specific value of cv, used only if sigma="fixCV"


# Need parts of the single process code
source("trendtests.R")


# Function to calculate the weights
findwvec <- function(tlist,weights,sigma,cv){
  wvec <- vector(length=tlist$m)
  if(weights=="equal")
    wvec <- rep(1/sqrt(tlist$m),tlist$m)
  if(weights=="tau"){
    wvec <- tlist$tauvec/sqrt(sum(tlist$tauvec^2))
  }
  if(weights=="sqrtN"){
    for(i in 1:tlist$m){
      n <- length(tlist$tvec[[i]])
      wvec[i] <- sqrt(n)
    }
    wvec <- wvec/sqrt(sum(wvec^2))
  }
  if(weights=="sqrtNCV"){
    for(i in 1:tlist$m){
      n=length(tlist$tvec[[i]])
      xvec=diff(c(0,tlist$tvec[[i]]))
      mu <- mean(xvec)
      if(sigma=="l"){
        xdiff <- xvec[2:n]-xvec[1:(n-1)]
        vest <- sum(xdiff^2)/(2*(n-1))
      }
      else
        vest <- var(xvec)
      CV <- sqrt(vest)/mu
      n <- length(tlist$tvec[[i]])
      wvec[i] <- sqrt(n)/CV
    }
    wvec <- wvec/sqrt(sum(wvec^2))
  }
  if(weights=="CVntau"){
    for(i in 1:tlist$m){
      n <- length(tlist$tvec[[i]])
      CV <- findCV(tvec=tlist$tvec[[i]],tau=tlist$tauvec[i],sigma,cv)
      wvec[i] <- CV*tlist$tauvec[i]*sqrt(n)
    }
    wvec <- wvec/sqrt(sum(wvec^2))
  }
  return(wvec)
}



## Lewis Robinson test

# Calculate the test observator
LRtestobs_multi <- function(tlist,weights,sigma,cv){
  wvec <- findwvec(tlist,weights, sigma, cv)
  LR <- 0
  for(i in 1:tlist$m)
    LR <- LR+wvec[i]*LRtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv) 
  return(LR)
}

# Run the test
LRtest_multi <- function(tlist,weights="sqrtNCV",sigma="s",cv=1){
  LR <- LRtestobs_multi(tlist,weights,sigma,cv)
  pLR <- 2*pnorm(-abs(LR))
  cat("Lewis-Robinson test for trend in multiple time censored processes.\n")
  cat(paste("Test statistic: LR =",round(LR,digits=3),"\n"))
  cat(paste("p-value =",round(pLR,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Extended Lewis Robinson test

# Calculate the test observator
ELRtestobs_multi <- function(tlist,weights,sigma,cv,avec){
  wvec <- findwvec(tlist,weights, sigma, cv)
  ELR <- 0
  for(i in 1:tlist$m)
    ELR <- ELR+wvec[i]*ELRtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv,avec[i]) 
  return(ELR)
}

# Run the test
ELRtest_multi <- function(tlist,weights="ntau",sigma="s",cv=1, avec){
  ELR <- ELRtestobs_multi(tlist,weights,sigma,cv,avec)
  pELR <- 2*pnorm(-abs(ELR))
  cat("Extended Lewis-Robinson test for trend in multiple time censored processes.\n")
  cat(paste("Test statistic: ELR =",round(ELR,digits=3),"\n"))
  cat(paste("p-value =",round(pELR,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
  cat(paste("Using avec: "))
  cat(round(avec,digits=3))
  cat("\n")
}


## Cramer von-Mises test, version 1 using simulated p-value

# Calculate the test observator and p-value
CvMtestobs_multi1 <- function(tlist,weights,sigma,cv,Npsim){
  wvec <- findwvec(tlist,weights,sigma,cv) 
  CvM0 <- 0
  for(i in 1:tlist$m)
    CvM0 <- CvM0+wvec[i]*CvMtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv)
  CvMvec <- wvec%*%matrix(qCvM(p=runif(tlist$m*Npsim)),nrow=tlist$m)
  pvalue <- mean(CvMvec>CvM0)
  list(CvM=CvM0,pCvM=pvalue)
}

# Run the test
CvMtest_multi1 <- function(tlist,weights="tau",sigma="s",cv=1,Npsim=1000){
  CvM <- CvMtestobs_multi1(tlist,weights,sigma,cv,Npsim)
  cat("Cramer von-Mises test for trend in multiple time censored processes, version 1.\n")
  cat(paste("Test statistic: CvM =",round(CvM$CvM,digits=3),"\n"))
  cat(paste("p-value =",round(CvM$pCvM,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Cramer von-Mises test, version 2 using normal approximation

# Calculate the test observator 
CvMtestobs_multi2 <- function(tlist,weights,sigma,cv){
  wvec <- findwvec(tlist,weights,sigma,cv) 
  CvM <- 0
  for(i in 1:tlist$m)
    CvM <- CvM+wvec[i]*CvMtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv)
  zCvM <- (CvM-sum(wvec)*0.166665)/sqrt(sum(wvec^2)*0.02224)
  return(zCvM)
}

# Run the test
CvMtest_multi2 <- function(tlist,weights="tau",sigma="s",cv=1){
  zCvM <- CvMtestobs_multi2(tlist,weights,sigma,cv)
  pzCvM <- 2*pnorm(-abs(zCvM))
  cat("Cramer von-Mises test for trend in multiple time censored processes, version 2.\n")
  cat(paste("Test statistic: zCvM =",round(zCvM,digits=3),"\n"))
  cat(paste("p-value =",round(pzCvM,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Anderson-Darling test, version 1 using simulated p-value

# Calculate the test observator and p-value

ADtestobs_multi1 <- function(tlist,weights,sigma,cv,Npsim){
  wvec <- findwvec(tlist,weights,sigma,cv) 
  AD0 <- 0
  for(i in 1:tlist$m)
    AD0 <- AD0+wvec[i]*ADtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv)
  ADvec <- wvec%*%matrix(qAD(p=runif(tlist$m*Npsim),n=Inf, fast=TRUE),nrow=tlist$m)
  pvalue <- mean(ADvec>AD0)
  list(AD=AD0,pAD=pvalue)
}

# Run the test
ADtest_multi1 <- function(tlist,weights="tau",sigma="s",cv=1,Npsim=1000){
  AD <- ADtestobs_multi1(tlist,weights,sigma,cv,Npsim)
  cat("Anderson-Darling test for trend in multiple time censored processes, version 1.\n")
  cat(paste("Test statistic: AD =",round(AD$AD,digits=3),"\n"))
  cat(paste("p-value =",round(AD$pAD,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Anderson-Darling test, version 2 using normal approximation

# Calculate the test observator 
ADtestobs_multi2 <- function(tlist,weights,sigma,cv){
  wvec <- findwvec(tlist,weights,sigma,cv) 
  AD <- 0
  ADvar <- 0
  for(i in 1:tlist$m){
    AD <- AD+wvec[i]*ADtestobs(tlist$tvec[[i]],tlist$tauvec[i],sigma,cv)
    ADvar <- ADvar+wvec[i]^2*(2*(pi^2-9)/3+(10-pi^2)/length(tlist$tvec[[i]]))
  }
  zAD <- (AD-sum(wvec))/sqrt(ADvar)
  return(zAD)
}

# Run the test
ADtest_multi2 <- function(tlist,weights="tau",sigma="s",cv=1){
  zAD <- ADtestobs_multi2(tlist,weights,sigma,cv)
  pzAD <- 2*pnorm(-abs(zAD))
  cat("Anderson-Darling test for trend in multiple time censored processes, version 2.\n")
  cat(paste("Test statistic: zAD =",round(zAD,digits=3),"\n"))
  cat(paste("p-value =",round(pzAD,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Linear rank test from Lawless, Cigsar and Cook (2012) 

# Calculate the test observator
LinRanktestobs_multi <- function(tlist){
  Vvec <- vector(length=tlist$m)
  VarVvec <- vector(length=tlist$m)
  for(j in 1:tlist$m){
    LinRanktest <- LinRanktestobs(tlist$tvec[[j]])
    Vvec[j] <- LinRanktest$V 
    VarVvec[j] <- LinRanktest$VarV 
  }
  mLinRank <- sum(Vvec)/sqrt(sum(VarVvec))
  return(mLinRank)
}

# Run the test 
LinRanktest_multi <- function(tlist){
  LinRank <- LinRanktestobs_multi(tlist)
  pLinRank <- 2*pnorm(-abs(LinRank))
  cat("Linear rank test for trend in multiple processes.\n")
  cat(paste("Test statistic: LinRank =",round(LinRank,digits=3),"\n"))
  cat(paste("p-value =",round(pLinRank,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



## Generalized Laplace test from Lawless, Cigsar and Cook (2012) 

# Calculate the test observator
GLtestobs_multi <- function(tlist){
  Ulist <- vector(length=tlist$m)
  for(i in 1:tlist$m)
    Ulist[i] <- sum(tlist$tvec[[i]])-length(tlist$tvec[[i]])*tlist$tauvec[i]/2
  GL <- sum(Ulist)/sqrt(sum(Ulist^2))
  return(GL)
}
  
# Run the test 
GLtest_multi <- function(tlist){
  GL <- GLtestobs_multi(tlist)
  pGL <- 2*pnorm(-abs(GL))
  cat("Generalized Laplace test for trend in multiple processes.\n")
  cat(paste("Test statistic: GL =",round(GL,digits=3),"\n"))
  cat(paste("p-value =",round(pGL,digits=5)),"\n")
  cat(paste("m =", tlist$m ,"processes\n"))
}



