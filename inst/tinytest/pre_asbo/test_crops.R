## tests for the CROPS implimnetation
set.seed(10)
x <- c(rnorm(100, 0, 1), rnorm(100, 1, 1), rnorm(100, 0, 1), rnorm(100, 0.2, 1))

fCost <- gaussMean$new(x)
expect_silent({ p <- crops(2*log(length(x)),50*log(length(x)),fCost) })

expect_silent({p <- crops(2*log(length(x)),50*log(length(x)),fCost,alg=capa,betaP=3*log(length(x)))})
