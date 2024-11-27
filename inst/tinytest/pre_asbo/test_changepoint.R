## test the pelt and op algorithms with example data
## test results from the changepoint package unless stated
## computed using th BIC not the MBIC
## changepoint returns the last timestep before a change

set.seed(10)
x <- c(rnorm(100, 0, 1), rnorm(100, 1, 1), rnorm(100, 0, 1), rnorm(100, 0.2, 1))
expect_silent({
    fCost <- gaussMean$new(x)
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(97,192,400) )

expect_silent({
    p <- partition(1.5*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(97,192,273,400) )

## check adding missing data has no impact
expect_silent({
    xx <- c(x[1:200],NA,NA,NA,x[201:400])
    fCost <- gaussMean$new(xx)
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(97,192,403) )


## check change in mu has no impact
expect_silent({
    x[50:150] <- x[50:150] + 5
    mu <- rep(0,length(x))
    mu[50:150] <- mu[50:150] + 5
    fCost <- gaussMean$new(x,mu)
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(97,192,400) )
          


expect_silent({
    data("Lai2005fig4")#, package = "changepoint")
    x <- Lai2005fig4[, 5]
    fCost <- gaussMean$new(x)
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(81,85,89,96,123,133,193) )

    
expect_silent({
    data("wind") #, package = "gstat")
    x <- diff(wind[, 11])
    fCost <- gaussVar$new(x)
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end,
             c(3409, 3496, 5054, 5184, 5203, 5373, 5583, 5678, 5728, 6235, 6241, 6542, 6573) )

expect_silent({
    data("discoveries", package = "datasets")
    x <- discoveries
    fCost <- poisCost$new(as.numeric(x))
    p <- partition(2*log(length(x)),NA,2)
    res <- pelt(p,fCost)
})
expect_equal( collective_anomalies(res)$end, c(24, 29, 73, 100) )


