### Examples of how to run the single process tests

# Clear the workspace
rm(list=ls())

# Load the tests
source("trendtests.R")


# A vector of failure time data:
LHD_ftimes <- c(16, 39, 71, 95, 98, 110, 114, 226, 294,
                344, 555, 599, 757, 822, 963, 1077, 1167, 1202,
                1257, 1317, 1345, 1372, 1402, 1536, 1625, 1643, 1675,
                1726, 1736, 1772, 1796, 1799, 1814, 1868, 1894, 1970)
# The data are from:
#  Kumar, Klefsjø and Granholm (1989), Reliability investigation for a fleet of load haul dump
#  machines in a Swedish mine, Reliability Engineering and System Safety 24: 341-361.

# The censoring time
LHD_tau <- 2000


# Plot the data
LHD_N <- length(LHD_ftimes) # Number of failures
plot(LHD_ftimes,1:LHD_N,xlab="Time (hours)",ylab="Cumulative number of failures")
lines(c(0,LHD_ftimes[LHD_N]),c(0,LHD_N),lty=2)


# Run the trend tests
LRtest(LHD_ftimes,LHD_tau,sigma="s")
KStest(LHD_ftimes,LHD_tau,sigma="s")
CvMtest(LHD_ftimes,LHD_tau,sigma="s")
ADtest(LHD_ftimes,LHD_tau,sigma="s")
ADtest(LHD_ftimes,LHD_tau,sigma="l")
ADtest(LHD_ftimes,LHD_tau,sigma="c")
ELRtest(LHD_ftimes,LHD_tau,sigma="s", a=0.5)
ELRtest(LHD_ftimes,LHD_tau,sigma="s", a=1/3)
ELRtest(LHD_ftimes,LHD_tau,sigma="s", a=2/3)
ELRtest(LHD_ftimes,LHD_tau,sigma="s", a=0)
ELRtest(LHD_ftimes,LHD_tau,sigma="s", a=1)
LinRanktest(LHD_ftimes)


# Using permutation rather than asymptotic distributions
Permtests(LHD_ftimes,LHD_tau,B=10000,sigma="s", a=0.5)
Permtests(LHD_ftimes,LHD_tau,B=10000,sigma="l", a=0.5)
Permtests(LHD_ftimes,LHD_tau,B=10000,sigma="c", a=0.5)

# More permutations - takes some time
Permtests(LHD_ftimes,LHD_tau,B=100000,sigma="s", a=0.5)











