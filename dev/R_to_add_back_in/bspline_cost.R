#' @export
splineCost <- R6Class("splineCost",
                      public=list(
                          summaryStats = NULL,
                          maxT = 0,
                          
                          non_neg = FALSE,
                          initialize = function(x,m=NULL,S=NULL,D ){
                              if( nrow(D) != nrow(x) ){ stop("incorrect dimensions") }
                              if( ncol(D) != length(m) ){ stop("incorrect model dimensions") }
                              self$summaryStats <- list(Y = x,
                                                        ybase = as.numeric(D%*%m),
                                                        D = D)
                              self$maxT <- ncol(x)
                              invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             if(a<1){
                                 Y <- self$summaryStats$Y[,1:b,drop=F]
                             }else{
                                 Y <- self$summaryStats$Y[,a:b,drop=F]
                             }
                             Y <- Y - self$summaryStats$ybase

                             n <- sum(is.finite(Y))
                             s <- sum(Y^2,na.rm=TRUE)
                             sigma2 <- (1/n)*s
                             
                             n*log(2*pi) + n*log( sigma2) + s + pen
                         },
                         pointCost = function(b,pen){
                             self$collectiveCost(b,b,pen)
                         },
                         collectiveCost = function(a,b,pen){
                             if(a<1){
                                 Y <- self$summaryStats$Y[,1:b,drop=F]
                             }else{
                                 Y <- self$summaryStats$Y[,a:b,drop=F]
                             }
                             nsamp <- ncol(Y) ##number of repeat observations
                             X <- self$summaryStats$D
                             XtX <- matrix(0,ncol(X),ncol(X))
                             Xty <- matrix(0,ncol(X),1)
                             
                             n <- 0
                             for(ii in 1:nsamp){
                                 idx <- is.finite(Y[,ii])
                                 n <- n + sum(idx)
                                 XtX <- XtX + t(X[idx,]) %*% X[idx,]
                                 Xty <- Xty + t(X[idx,])%*%Y[idx,ii]
                             }
                             
                             allCoeff <- as.vector(solve( XtX, Xty )  ) ##t(X)%*%X, t(X)%*%y))

                             Y <- Y - as.vector(X%*%allCoeff)  # residuals
                             s <- sum(Y^2,na.rm=TRUE)
                             sigma2 <- (1/n)*s
                             
                             n*log(2*pi) + n*log( sigma2) + s + pen
                             
                         }
                         
                      )
                      )

