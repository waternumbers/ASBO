## need to check not all values in a column of X are equal - in which case loess fails!
#' R6 class for local regression Cost Functions
#'
#' @description VERY experimental - do not use
#'
#' @export
localRegCost <- R6Class(
    "localRegCost",
    private = list(
        summaryStats = NULL,
        maxT = 0,
        theta0 = NA,
        sigma0 = NA,
        family = NA,
        non_neg = FALSE
    ),
    public=list(
        #' @description Get the length of time series
        length = function(){ private$maxT },
        #' @description Initialise the cost function
        #' @param x observations as for gauss_reg
        #' @param family for fitting
        initialize = function(x,family = c("gaussian", "symmetric")){
            private$family <- match.arg(family)
            private$summaryStats <- list(y = list(), X=list())
            for(ii in 1:length(x)){
                private$summaryStats$y[[ii]] <- x[[ii]]$y
                private$summaryStats$X[[ii]] <- x[[ii]]$X
            }
            private$maxT <- length(x)
            invisible(self)
        },
        #' @description Compute the non-anomalous cost of a segment
        #' @param a start of period
        #' @param b end of period
        #' @param pen penalty cost
        baseCost = function(a,b,pen=0){Inf},
        #' @description Compute the point anomaly cost of a time step
        #' @param b time step
        #' @param pen penalty cost
        pointCost = function(b,pen){
            if( length( private$summaryStats$y[[b]] ) == 1 ){ return(Inf) }
            
            self$collectiveCost(b,b,pen,1)
        },
        #' @description Compute the non-anomalous cost of a segment
        #' @param a start of period
        #' @param b end of period
        #' @param pen penalty cost
        #' @param len minimum number of observations
        collectiveCost = function(a,b,pen,len){
            if( (b-a+1) < len ){ return(NA) } ## check length
            D <- do.call(rbind,private$summaryStats$X[a:b])
            D <- data.frame( y = unlist( private$summaryStats$y[a:b] ),
                            D)
            
            yhat <- predict( loess(y~.,data=D,degree=1,family=private$family) )
            
            ## sigma <- sqrt( mean( (D$y-yhat)^2) )
            sigma <- mad( D$y-yhat )
            -2*sum(dnorm(D$y,yhat,sigma,log=TRUE)) + pen
        },
        #' @description Compute parameters of a segment if anomalous
        #' @param a start of period
        #' @param b end of period
        param = function(a,b){
            D <- do.call(rbind,private$summaryStats$X[a:b])
            D <- data.frame( y = unlist( private$summaryStats$y[a:b] ),
                            D)
            yhat <- predict( loess(y~.,data=D,degree=1, family=private$family) )
            sigma <- mad( D$y-yhat )
            ##                             sigma <- sqrt( mean( (D$y-yhat)^2) )
            n <- c(0,sapply(private$summaryStats$y[a:b],length))
            n <- cumsum(n)
            ##browser()
            out <- list()
            for(ii in 2:length(n)){
                idx <- (n[ii-1]+1):n[ii]
                out[[ii-1]] <- data.frame(mu = yhat[idx], sigma = sigma)
            }
            return(out)
        }        
    )
)
