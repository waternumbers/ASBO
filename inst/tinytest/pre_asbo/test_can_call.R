## test that the penalty functions can be called and return silently
#rm(list=ls())
#library(tinytest)
#devtools::load_all()


## ###########################
## Univariate continuous
set.seed(10)
x <- c(rnorm(100),rnorm(100,10,1))

expect_silent({
    beta <- 4*log(length(x))
    betaP <- 3*log(length(x))
    fCost <- gaussMeanVar$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 4*log(length(x))
    betaP <- 3*log(length(x))
    fCost <- gaussMean$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 4*log(length(x))
    betaP <- 3*log(length(x))
    fCost <- gaussVar$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 4*log(length(x))
    betaP <- 3*log(length(x))
    fCost <- ladCost$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

## #################################
## univariate count
set.seed(10)
x <- c(rpois(100,0.1),rpois(100,3))
expect_silent({
    beta <- 3*log(length(x))
    betaP <- 2*log(length(x))
    fCost <- poisCost$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})




## #################################
## multivariate count
set.seed(10)
x <- matrix(sample(10,200*4,replace=TRUE),200,4)
expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- multinomialCost$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

set.seed(10)
x <- matrix(rnorm(800),200,4)
x <- (x==apply(x,1,max))*as.integer(1)
expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- categoricalCost$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

## #################################
## multivariate
set.seed(10)
x <- rbind( matrix(rnorm(400),100,4), matrix(rnorm(400,10,1),100,4) )

expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- rankCost$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- gaussMean$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- gaussMeanVar$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})

expect_silent({
    beta <- 6*log(length(x))
    betaP <- 5*log(length(x))
    fCost <- gaussVar$new(x)
    p <- partition(beta,betaP,2)
    res <- capa(p,fCost)
})



## ###################################
## regressions of various forms

## gauss_reg
set.seed(10)
n <- 120
x <- list()
for(ii in 1:48){
    
    if(ii < 10){ theta = c(1,0); sigma <- 0.1 }
    if(ii >= 10 & ii <12){ theta <- c(10,0); sigma <- 2}
    if(ii >= 12 & ii < 44){ theta <- c(5,1); sigma <- 2}
    if(ii >= 44 ){ theta <- c(1,0); sigma <- 0.1}
    
    X <- cbind(rep(1,n),runif(n,ii-1,ii))
    y <- rnorm(n, X%*%theta, sigma)
    x[[ii]] <- list(y=y,X=X)
}
n <- sum(sapply(x,function(xx){length(xx$y)}))
expect_silent({
    fCost <- gaussRegMeanVar$new(x)
    p <- partition(4*log(n),4*log(n) ,2)
    res <- capa(p,fCost)
})

## local reg
xx <- lapply(x,function(z){z$X <- z$X[,2,drop=FALSE]; z})
expect_silent({
    fCost <- localRegCost$new(xx)
    p <- partition(4*log(n),4*log(n) ,2)
    res <- capa(p,fCost)
})

## bspline
#set.seed(10)
#x <- rbind( matrix(rnorm(400),100,4), matrix(rnorm(400,10,1),100,4) )
#D <- matrix(

