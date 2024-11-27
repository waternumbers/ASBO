## tests for the poisson distribution cost functions
## test results from the changepoint and CAPA packages unless stated
## changepoint returns the last timestep before a change

#library(tinytest)
#library(ASBO)

## ##############################################
## Test we can call the various forms of cost function
set.seed(1)
x <- rpois(100,lambda=1)
lhat <- mean(x)

## base cost
expect_silent({ fC <- poisBase$new(x,gamma = 2) })
expect_equal( fC$k(), -2 )
expect_equal( fC$min_length(), 1 )
expect_true( is.na( fC$param(1,100)) )
expect_equal( fC$length(), 100 )
expect_equal( fC$cost(1,100), -2*sum(dpois(x,lambda=1,log=TRUE)) + 2 )

## proportional cost
expect_silent({ fC <- poisProp$new(x,gamma = 2) })
expect_equal( fC$k(), -2 )
expect_equal( fC$min_length(), 1 )
expect_equal( fC$param(1,100), lhat )
expect_equal( fC$length(), 100 )
expect_equal( fC$cost(1,100), -2*sum(dpois(x,lambda=lhat,log=TRUE)) + 2 )

## proportional cost and penalty
expect_silent({ fC <- poisPropG$new(x,gamma = 2) })
expect_equal( fC$k(), 0 )
expect_equal( fC$min_length(), 1 )
expect_equal( fC$param(1,100), lhat )
expect_equal( fC$length(), 100 )
expect_equal( fC$cost(1,100), -2*sum(dpois(x,lambda=lhat,log=TRUE)) + 100*2 )

## ###############################################
## example from changepoint package
expect_silent({
    data("discoveries", package = "datasets")
    x <- discoveries
    fCost <- poisProp$new(as.numeric(x),gamma=0)
    p <- partition(2*log(length(x)))
    res <- pelt(p,list("pop"= fCost),min_length=2)
})
expect_equal( summary(res)$end, c(24, 29, 73, 100) )


