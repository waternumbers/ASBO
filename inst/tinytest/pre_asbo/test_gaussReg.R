## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
## computed using th BIC not the MBIC
## changepoint returns the last timestep before a change
##rm(list=ls())
##devtools::load_all()

## simple test
set.seed(10)
x <- list()
n <- 120
for(ii in 1:48){
    
    if(ii < 10){ theta = c(1,0); sigma <- 0.1 }
    if(ii >= 10 & ii <12){ theta <- c(10,0); sigma <- 2}
    if(ii >= 12 & ii < 44){ theta <- c(5,1); sigma <- 2}
    if(ii >= 44 ){ theta <- c(1,0); sigma <- 0.1}
    
    X <- cbind(rep(1,n),runif(n,ii-1,ii))
    y <- rnorm(n, X%*%theta, sigma)
    x[[ii]] <- list(y=y,X=X)
}

expect_silent({
    fCost <- gaussRegMeanVar$new(x)
    p <- partition(4*log(sum(sapply(y,length))),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(9,11,43,48) )


## test with singular design matrix
xx <- x
for(ii in c(1:10,44:48)){
    xx[[ii]]$X[,2] <- 0
}
expect_silent({
    fCost <- gaussRegMeanVar$new(xx)
    p <- partition(4*log(sum(sapply(y,length))),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(9,11,43,48) )

## #############################################################
## replicate univariate tests in capa
data(simulated)

## Adapt data to contain a clear point anomaly
X <- sim.data
X[100,1] <- 45

## adapt X to match the robust scaling applied in v4.0.2
X <- apply(X,2,function(x){ (x-median(x))/mad(x) })

## read in the results from v4.0.2
out <- readRDS("capa_results402_v2.rds")

x <- rep(list(list(y=NA,X=matrix(1,1,1))),nrow(X))
for(ii in 1:length(x)){
    x[[ii]]$y <- X[ii,1]
}

beta <- 4*log(length(x))
betaP <- 3*log(length(x))

expect_silent({
    fCost <- gaussRegMeanVar$new(x)
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
    fCost <- gaussRegMean$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})
expect_equal( point_anomalies(res)$location, out$single_mean$point$location )
expect_equivalent( collective_anomalies(res)[,c("start","end")],
                  out$single_mean$collective[,c("start","end")] )
## expect_silent({ summary(res) })
## expect_silent({ show(res) })
## expect_silent({ plot(res,variate_names=TRUE) })
