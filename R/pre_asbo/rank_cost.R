#' R6 class for rank based anomalies
#'
#' @description This is VERY developmental - do not use
#' 
#' @examples
#' set.seed(0)
#' m <- runif(100)
#' s <- pmax(1e-4,runif(100))
#' x <- rnorm(100,m,s) ## example data
#' @export
rankCost <- R6Class("rankCost",
                    private = list(
                        summaryStats = NULL,
                        nStep = NA,
                        Sr = NA,
                        maxT = 0
                    ),
                    public=list(
                        #' @description Get the length of time series
                        length = function(){ private$maxT },
                        #' @description Initialise the cost function
                        #' @param x numeric matrix of observations
                        #' @param m numeric vector or matrix of mean values
                        initialize = function(x,m=0){
                            x <- as.matrix(x)
                            x <- x-m
                            private$maxT <- nrow(x)
                            S <- apply(x,2,rank,na.last="keep") - (nrow(x)+1)/2
                            private$summaryStats <- apply(S,2,cumsumNA)
                            private$nStep <- apply(is.finite(S),2,cumsumNA)
                            private$Sr <- var(S,na.rm=T)
                            invisible(self)
                        },
                        #' @description Compute the non-anomalous cost of a segment
                        #' @param a start of period
                        #' @param b end of period
                        #' @param pen penalty cost
                        baseCost = function(a,b,pen=0){
                            self$collectiveCost(a,b,pen,1)
                        },
                        #' @description Compute the point anomaly cost of a time step
                        #' @param b time step
                        #' @param pen penalty cost
                        pointCost = function(b,pen){
                            self$collectiveCost(b,b,pen,1)
                        },
                        #' @description Compute the non-anomalous cost of a segment
                        #' @param a start of period
                        #' @param b end of period
                        #' @param pen penalty cost
                        #' @param len minimum number of observations
                        collectiveCost = function(a,b,pen,len){
                            a <- a-1
                            if(a<1){
                                sumStat <- self$summaryStats[b,]
                                nS <- self$nStep[b,]
                            }else{
                                sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                nS <- self$nStep[b,] - self$nStep[a,]
                            }
                            if( any(is.na(nS) | nS<len) ){ return(NA) } ## check length and if NA
                            
                            rbar <- sumStat/nS
                            idx <- is.finite(rbar)
                            if(!any(idx)){ return(NA) }
                            
                            -(b-a-1)* rbar %*% solve(self$Sr,rbar) + pen
                        },
                        #' @description Compute parameters of a segment if anomalous
                        #' @param a start of period
                        #' @param b end of period
                        param = function(a,b){ return(NULL) }
                    ))

                     
