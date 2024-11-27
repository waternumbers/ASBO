#' @export
gaussRegRobCost <- R6Class("gaussRegCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         theta0 = NA,
                         sigma0 = NA,
                         initialize = function(x,theta0=NULL,sigma0=1){
                             self$summaryStats <- list(y=lapply(x,function(x){x$y}),
                                                       X=lapply(x,function(x){x$X}))
                             self$maxT <- length(x)
                             if(is.null(theta0)){ theta0 <- rep(0,ncol(x[[1]]$X)) }
                             if(length(theta0) != ncol(x[[1]]$X)){ stop("theta0 incorrect size") }
                             self$theta0 <- theta0
                             self$sigma0 <- sigma0
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             X <- do.call(rbind,self$summaryStats$X[a:b])
                             y <- unlist( self$summaryStats$y[a:b] )
                             sigma <- self$sigma0
                             theta <- self$theta0
                             -2*sum(dnorm(y,X%*%theta,sqrt(sigma),log=TRUE)) + pen
                         },
                         pointCost = function(a,pen){
                             self$collectiveCost(a,a,pen)
                         },
                         collectiveCost = function(a,b,pen){
                             p <- self$param(a,b)
                             p$cst + pen
                         },
                         param = function(a,b){
                             X <- do.call(rbind,self$summaryStats$X[a:b])
                             y <- unlist( self$summaryStats$y[a:b] )

                             ## Check for non singular regression
                             QR <- qr(crossprod(X))                 # Get the QR decomposition
                             vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
                             XX <- X[,vars]
                             ##browser()
                             
                             mdl <- robustbase::lmrob(y~XX-1,
                                                      control=robustbase::lmrob.control(fast.s.large.n = Inf))
                             ## mdl <- robustbase::ltsReg(XX,y) ## to slow
                             ## solve- set other values to 0## trim variables
                             ##mdl <- rq(y~XX-1, tau=0.5)
                             theta <- rep(0,ncol(X))
                             theta[vars] <- coef( mdl )
                             sigma <- mdl$scale^2 ## mad( residuals(mdl),center=0 )
                             cst <- -2*sum(dnorm(residuals(mdl),0,sqrt(sigma),log=TRUE))
                             list(theta=theta, sigma=sigma, cst=cst)
                         }
                     )
                     )

