library(goftest)
source("Trendtests.R")
source("SimAsympTestDistr.R") # Takes some time, consider modifying the 
                              # settings in the call to simAELRdists


# Failure times for USS Halfbeak, section 6.1
USSH_ftimes <- c(1.382,   2.990,   4.124,   6.827,   7.472,   7.567,   8.845,   
                9.450,   9.794,   10.848, 11.993,  12.300,   15.413,   16.497,   
                17.352,   17.632,   18.122,   19.067,   19.172,   19.299,
                19.360,  19.686,   19.940,   19.944) 


USSH_tau <- 20
USSH_N <- length(USSH_ftimes)

USSH_xtimes <- diff(c(0,USSH_ftimes))
mean(USSH_xtimes)
sd(USSH_xtimes)
cv <- sd(USSH_xtimes)/mean(USSH_xtimes)
cv


# USSH
plot(USSH_ftimes,1:USSH_N,xlab="Time (hours)",ylab="Cumulative number of failures",
     pch=19, col="maroon", main=" ")
lines(c(0,USSH_ftimes[USSH_N]),c(0,USSH_N),lty=3)


LR <- LRtid(USSH_ftimes,USSH_tau,sigma="s")
LR
pLRtid <- 2*pnorm(-abs(LR))
pLRtid

pKS <- function(x,kmax=1000) 
  1-(sqrt(2*pi)/x)*sum(exp(-(2*(1:kmax)-1)^2*pi^2/(8*x^2)))
KS <- KStid(USSH_ftimes,USSH_tau,sigma="s")
KS
pKStid <- pKS(KS)
pKStid


CVm <- CVtid(USSH_ftimes,USSH_tau,sigma="s")
CVm
pCVmtid <- 1-pCvM(q=CVm)
pCVmtid

AD <- ADtid(USSH_ftimes,USSH_tau,sigma="s")
AD
pADtid <- 1-pAD(q=AD)
pADtid

ILR1 <- INTLR1tid(USSH_ftimes,USSH_tau,sigma="s")
ILR1
pILR1tid <- 2*pnorm(-abs(ILR1))
pILR1tid

ILR2 <- INTLR2tid(USSH_ftimes,USSH_tau,sigma="s")
ILR2
pILR2tid <- 2*pnorm(-abs(ILR2))
pILR2tid

ICVKS <- INTCvMKS1tid(USSH_ftimes,USSH_tau, sigma="s", agrid=seq(0,1,by=0.001))
ICVKS
ICvM <- ICVKS$ICvM
pICvM <-  mean(INT1CvMdist>ICvM)
pICvM
IKS <- ICVKS$IKS
pIKS <-  mean(INT1KSdist>IKS)
pIKS


AELRs <- AELR1tid(USSH_ftimes,USSH_tau, sigma="s", agrid=seq(0,1,by=0.001))
AELRs  
IELR1 <- AELRs$AELR1int1
#pIELR1 <- 2*pnorm(-abs(IELR1),sd=sqrt(0.17493)) # Two-sided
pIELR1 <- 1-pnorm(IELR1,sd=sqrt(0.17493))
pIELR1

SELR1 <- AELRs$AELR1sup1
SELR1
pSELR1 <-  mean(AELR1sup1dist>SELR1)
pSELR1

ELR50 <- ELRtid1(USSH_ftimes,USSH_tau, sigma="s",a=0.5)
ELR50
pELR50tid <- pnorm(-abs(ELR50))
pELR50tid






