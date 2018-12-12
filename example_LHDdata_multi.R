### Examples of how to run the multiple process tests

# Clear the workspace
rm(list=ls())

# Load the tests
source("trendtests_multi.R")

## Read in data
# The following data are from:
#  Kumar and Klefsjø (1992), Reliability analysis of hydraulic systems of LHD machines
#  using the power law process model, Reliability Engineering and System Safety 35: 217-224.
# The data are reported as times between events for machine no 1, 3, 9, 11, 17 and 20.

LHD1_xtimes <- c(327, 125, 7, 6, 107, 277, 54, 332, 510, 110, 10, 9, 85, 
                 27, 59, 16, 8, 34, 21, 152, 158, 44, 18) # Times between events
LHD1_ftimes <- cumsum(LHD1_xtimes) # Times to event
LHD1_tau <- tail(LHD1_ftimes,1)   # Use the last event time as censoring time
LHD1_ftimes <- LHD1_ftimes[LHD1_ftimes<LHD1_tau] # Remove the censoring time from the event times

LHD3_xtimes <- c(637, 40, 397, 36, 54, 53, 97, 63, 216, 118, 125, 25, 4,
                 101, 184, 167, 81, 46, 18, 32, 219, 405, 20, 248, 140)
LHD3_ftimes <- cumsum(LHD3_xtimes)
LHD3_tau <- tail(LHD3_ftimes,1)
LHD3_ftimes <- LHD3_ftimes[LHD3_ftimes<LHD3_tau] 

LHD9_xtimes <- c(278, 261, 990, 191, 107, 32, 51, 10, 132, 176, 247, 165, 454, 142,
                 38, 249, 212, 204, 182, 116, 30, 24, 32, 38, 10, 311, 61)
LHD9_ftimes <- cumsum(LHD9_xtimes)
LHD9_tau <- tail(LHD9_ftimes,1)
LHD9_ftimes <- LHD9_ftimes[LHD9_ftimes<LHD9_tau] 

LHD11_xtimes <- c(353, 96, 49, 211, 82, 175, 79, 117, 26, 4, 5, 60, 39, 35, 258, 
                  97, 59, 3, 37, 8, 245, 79, 49, 31, 259, 283, 150, 24)
LHD11_ftimes <- cumsum(LHD11_xtimes)
LHD11_tau <- tail(LHD11_ftimes,1)
LHD11_ftimes <- LHD11_ftimes[LHD11_ftimes<LHD11_tau] 

LHD17_xtimes <- c(401, 36, 18, 159, 341, 171, 24, 350, 72, 303, 34, 45, 324, 2, 70, 57, 
                  103, 11, 5, 3, 144, 80, 53, 84, 218, 122)
LHD17_ftimes <- cumsum(LHD17_xtimes)
LHD17_tau <- tail(LHD17_ftimes,1)
LHD17_ftimes <- LHD17_ftimes[LHD17_ftimes<LHD17_tau] 

LHD20_xtimes <- c(231, 20, 361, 260, 176, 16, 101, 293, 5, 119, 9, 80, 112, 10, 
                  162, 90, 176, 370, 90, 15, 315, 32, 266)
LHD20_ftimes <- cumsum(LHD20_xtimes)
LHD20_tau <- tail(LHD20_ftimes,1)
LHD20_ftimes <- LHD20_ftimes[LHD20_ftimes<LHD20_tau] 

# Number of processes considered
npros <- 6
# Vector of all censoring times 
ttauvec <- c(LHD1_tau,LHD3_tau,LHD9_tau,LHD11_tau,LHD17_tau,LHD20_tau)

# List of the data set up in the required format
SBdatalist <- list(tvec=list(LHD1_ftimes,LHD3_ftimes,LHD9_ftimes,LHD11_ftimes,LHD17_ftimes,LHD20_ftimes),
                   tauvec=ttauvec,m=npros)


## Nelson-Aalen plot of the data 
nrisk <- npros:1
taus <- sort(ttauvec)
etimes <- sort(c(LHD1_ftimes,LHD3_ftimes,LHD9_ftimes,LHD11_ftimes,LHD17_ftimes,LHD20_ftimes))
nrisks <- numeric(length(etimes))
for(i in 1:length(etimes))
  nrisks[i] <- head(nrisk[taus>etimes[i]],1)
NAest <- cumsum(1/nrisks)
plot(etimes,NAest,type="s",xlab="Time (hours)",ylab="Estimated mean cumulative number",
     lwd=2,col="blue", main="Nelson-Aalen plot")
points(etimes,NAest,col="blue")
lines(c(head(etimes,1),tail(etimes,1)),c(0,tail(NAest,1)),lty=3,lwd=1,col="red")


## Run the trend tests
LRtest_multi(SBdatalist,weights="sqrtNCV",sigma="s")
maxtau <- max(ttauvec) # Used to define a in each process in ELR below to correspond to maxtau/2
ELRtest_multi(SBdatalist,weights="sqrtNCV",sigma="s",avec=(1-0.5*ttauvec/maxtau))
CvMtest_multi1(SBdatalist,weights="tau",sigma="s",Npsim=10000)
CvMtest_multi2(SBdatalist,weights="tau",sigma="s")
ADtest_multi1(SBdatalist,weights="tau",sigma="s",Npsim=10000)
ADtest_multi2(SBdatalist,weights="tau",sigma="s")
LinRanktest_multi(SBdatalist)
GLtest_multi(SBdatalist)



