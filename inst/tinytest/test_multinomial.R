## tests for the multinomial distribution cost functions

#library(tinytest)
#devtools::load_all()
#library(ASBO)

## ##############################################
## Test we can call the various forms of cost function
set.seed(1)
m <- c(1:4)/sum(1:4)
X <- t(rmultinom(100, 144, m))
mhat <- colSums(X)/sum(X)

## base cost
expect_silent({ fC <- multinomialBase$new(X,m=m,gamma = 2) })
expect_equal( fC$k(), -2 )
expect_equal( fC$min_length(), 1 )
expect_true( is.na( fC$param(1,100)) )
expect_equal( fC$length(), 100 )
expect_equal( fC$cost(1,100), -2*sum(apply(X,1,function(x){dmultinom(x,prob=m,log=TRUE)})) + 2 )

## proportional cost
expect_silent({ fC <- multinomialReplace$new(X,gamma = 2) })
expect_equal( fC$k(), -2 )
expect_equal( fC$min_length(), 1 )
expect_equal( fC$param(1,100), mhat )
expect_equal( fC$length(), 100 )
expect_equal( fC$cost(1,100), -2*sum(apply(X,1,function(x){dmultinom(x,prob=mhat,log=TRUE)})) + 2 )

