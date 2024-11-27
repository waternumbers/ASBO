rm(list=ls())
devtools::load_all()
set.seed(0)
m <- runif(100)
s <- pmax(1e-4,runif(100))
x <- rnorm(100,m,s) ## example data
gM <- gaussMean$new(x,m,s) ## anomalies are changes in mean
gM$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
gM$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
gM$collectiveCost(90,95,0,3) ## collective anomaly cost for x[90:95] with 0 penalty and at least 3 observation
gV <- gaussVar$new(x,m,s) ## anomalies are changes in variance
gV$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
gV$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
gV$collectiveCost(90,95,57,3) ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
gMV <- gaussMeanVar$new(x,m,s) ## anomalies are changes in mean and variance
gMV$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
gMV$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
gMV$collectiveCost(90,95,57,3) ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation

X <- matrix(rnorm(500,m,s),100,5)

gM <- gaussMean$new(X,m,s) ## anomalies are changes in mean
gM$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
gM$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
gM$collectiveCost(90,95,0,3) ## collective anomaly cost for x[90:95] with 0 penalty and at least 3 observation

devtools::load_all()
set.seed(0)
m <- c(1:4)/sum(1:4)
X <- t(rmultinom(100, 1, m))
p <- categoricalCost$new(X,m)
p$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
p$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
p$collectiveCost(90,95,57,3) ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
