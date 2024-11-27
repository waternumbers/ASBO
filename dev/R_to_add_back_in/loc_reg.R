## need to check not all values in a column of X are equal - in which case loess fails!

#' @export
localRegCost <- R6Class("localRegCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         theta0 = NA,
                         sigma0 = NA,
                         family = NA,
                         non_neg = FALSE,
                         initialize = function(x,family = c("gaussian", "symmetric")){

                             self$family <- match.arg(family)
                             
                             self$summaryStats <- list(y = list(), X=list())
                             
                             for(ii in 1:length(x)){
                                 self$summaryStats$y[[ii]] <- x[[ii]]$y
                                 self$summaryStats$X[[ii]] <- x[[ii]]$X
                             }
                             
                             self$maxT <- length(x)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){Inf},
                         pointCost = function(b,pen){
                             if( length( self$summaryStats$y[[b]] ) == 1 ){ return(Inf) }

                             self$collectiveCost(b,b,pen,1)
                         },
                         collectiveCost = function(a,b,pen,len){
                             if( (b-a+1) < len ){ return(NA) } ## check length
                             ##browser()
                             ## X <- do.call(rbind,self$summaryStats$X[a:b])
                             ## y <- unlist( self$summaryStats$y[a:b] )
                             ## yhat <- quantdr::llqr(X,y)$ll_est
                             
                             ## chck = function(x){x*(self$tau - (x<0))}
                             ## 2*( sum(chck(y-yhat)) + length(y)*log(0.5*(1-0.5)) ) + pen
                             
                             D <- do.call(rbind,self$summaryStats$X[a:b])
                             D <- data.frame( y = unlist( self$summaryStats$y[a:b] ),
                                             D)
                             
                             yhat <- predict( loess(y~.,data=D,degree=1,family=self$family) )
                             
                             ## sigma <- sqrt( mean( (D$y-yhat)^2) )
                             sigma <- mad( D$y-yhat )
                             -2*sum(dnorm(D$y,yhat,sigma,log=TRUE)) + pen
                         },
                         param = function(a,b){
                             D <- do.call(rbind,self$summaryStats$X[a:b])
                             D <- data.frame( y = unlist( self$summaryStats$y[a:b] ),
                                             D)
                             yhat <- predict( loess(y~.,data=D,degree=1, family=self$family) )
                             sigma <- mad( D$y-yhat )
##                             sigma <- sqrt( mean( (D$y-yhat)^2) )
                             n <- c(0,sapply(self$summaryStats$y[a:b],length))
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
