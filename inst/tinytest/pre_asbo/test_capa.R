## tests for the capa class objects
## these tests are based on matching the output of anomaly v4.0.2

## NOTES:
## 1. expect_equivalent used for testing collective anomalies since row.names differ
## 2. multivariate robust mean checks cancelled out since take to long to run - TODO fix this

data(simulated)

## Adapt data to contain a clear point anomaly
X <- sim.data
X[100,1] <- 45

## adapt X to match the robust scaling applied in v4.0.2
X <- apply(X,2,function(x){ (x-median(x))/mad(x) })

## read in the results from v4.0.2
out <- readRDS("capa_results402_v2.rds") ## TODO need changing to remove path

## #################################
## univariate tests
x <- X[,1]
beta <- 4*log(length(x))
betaP <- 3*log(length(x))

expect_silent({
    fCost <- gaussMeanVar$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})
expect_equal( point_anomalies(res)$location, out$single_meanvar$point$location )
expect_equivalent( collective_anomalies(res)[,c("start","end")],
                  out$single_meanvar$collective[,c("start","end")] )
## expect_silent({ summary(res) })
## expect_silent({ show(res) })
## expect_silent({ plot(res,variate_names=TRUE) })

expect_silent({
    fCost <- gaussMean$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})
expect_equal( point_anomalies(res)$location, out$single_mean$point$location )
expect_equivalent( collective_anomalies(res)[,c("start","end")],
                  out$single_mean$collective[,c("start","end")] )
## expect_silent({ summary(res) })
## expect_silent({ show(res) })
## expect_silent({ plot(res,variate_names=TRUE) })



##expect_silent({ res<-capa(X[,1,drop=F],type="robustmean",min_seg_len=2) })
## expect_equal( point_anomalies(res), out$single_robustmean$point ) ## fails in v4.0.2
##expect_equivalent( collective_anomalies(res), out$single_robustmean$collective )
## expect_silent({ summary(res) }) ## fails in v4.0.2
## expect_silent({ show(res) }) ## fails in v4.0.2
## expect_silent({ plot(res,variate_names=TRUE) }) ## fails in v4.0.2


