library(goftest)
source("Trendtests.R")
source("SimAsympTestDistr.R") # Takes some time, consider modifying the 
                              # settings in the call to simAELRdists

# Failure times for load-haul-dump data, section 6.2
LHD_ftimes <- c(16, 39, 71, 95, 98, 110, 114, 226, 294,
                344, 555, 599, 757, 822, 963, 1077, 1167, 1202,
                1257, 1317, 1345, 1372, 1402, 1536, 1625, 1643, 1675,
                1726, 1736, 1772, 1796, 1799, 1814, 1868, 1894, 1970)

LHD_tau <- 2000
LHD_N <- length(LHD_ftimes)

LHD_xtimes <- diff(c(0,LHD_ftimes))
mean(LHD_xtimes)
sd(LHD_xtimes)
cv <- sd(LHD_xtimes)/mean(LHD_xtimes)
cv


plot(LHD_ftimes,1:LHD_N,xlab="Time (hours)",ylab="Cumulative number of failures",
     pch=19, col="maroon", main=" ")
lines(c(0,LHD_ftimes[LHD_N]),c(0,LHD_N),lty=3)


LR <- LRtid(LHD_ftimes,LHD_tau,sigma="s")
LR
pLRtid <- 2*pnorm(-abs(LR))
pLRtid


pKS <- function(x,kmax=1000) 
  1-(sqrt(2*pi)/x)*sum(exp(-(2*(1:kmax)-1)^2*pi^2/(8*x^2)))

KS <- KStid(LHD_ftimes,LHD_tau,sigma="s")
KS
pKStid <- pKS(KS)
pKStid

CVm <- CVtid(LHD_ftimes,LHD_tau,sigma="s")
CVm
pCVmtid <- 1-pCvM(q=CVm)
pCVmtid

AD <- ADtid(LHD_ftimes,LHD_tau,sigma="s")
AD
pADtid <- 1-pAD(q=AD)
pADtid

ILR1 <- INTLR1tid(LHD_ftimes,LHD_tau,sigma="s")
ILR1
pILR1tid <- 2*pnorm(-abs(ILR1))
pILR1tid

ILR2 <- INTLR2tid(LHD_ftimes,LHD_tau,sigma="s")
ILR2
pILR2tid <- 2*pnorm(-abs(ILR2))
pILR2tid

ICVKS <- INTCvMKS1tid(LHD_ftimes,LHD_tau, sigma="s", agrid=seq(0,1,by=0.001))
ICVKS
ICvM <- ICVKS$ICvM
pICvM <-  mean(INT1CvMdist>ICvM)
pICvM
IKS <- ICVKS$IKS
pIKS <-  mean(INT1KSdist>IKS)
pIKS

AELRs <- AELR1tid(LHD_ftimes,LHD_tau, sigma="s", agrid=seq(0,1,by=0.001))
AELRs  
IELR1 <- AELRs$AELR1int1
#pIELR1 <- 2*pnorm(-abs(IELR1),sd=sqrt(0.17493)) # Two-sided
pIELR1 <- 1-pnorm(IELR1,sd=sqrt(0.17493))
pIELR1

SELR1 <- AELRs$AELR1sup1
SELR1
pSELR1 <-  mean(AELR1sup1dist>SELR1)
pSELR1

ELR50 <- ELRtid1(LHD_ftimes,LHD_tau, sigma="s",a=0.5)
ELR50
pELR50tid <- pnorm(-abs(ELR50))
pELR50tid

  
