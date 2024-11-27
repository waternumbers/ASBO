#' R6 class for Univariate Poisson Cost Functions
#'
#' @description Cost functions for the univariate Poission distribution
#'
#' @details Collective anomalies are represented as multiplicative changes in rate
#' 
#' @examples
#' set.seed(0)
#' r <- 8 + runif(100)*2
#' x <- rpois(100,lambda = r)
#'
#' p <- poisCost$new(x,r)
#' p$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' p$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' p$collectiveCost(90,95,57,3) 
#' @export
poisCost <- R6Class("poisCost",
                    private = list(
                        summaryStats = NULL,
                        maxT = 0
                    ),
                    public=list(
                        #' @description Get the length of time series
                        length = function(){ private$maxT },
                        #' @description Initialise the cost function
                        #' @param x numeric vector of observations
                        #' @param rate numeric vector of rate parameters
                        initialize = function(x,rate=1){
                            x <- as.numeric(x)
                            rate <- as.numeric(rate)
                            idx <- is.finite(rate)
                            stopifnot(
                                "zero or negative rates are not allowed" = all(rate[idx]>0),
                                "negative x values rates are not allowed" = all(x>=0)
                            )
                            S <- cbind(x,lfactorial(x),rate,x*log(rate),1)
                            S[is.na(x)|!idx,] <- NA
                            private$summaryStats <- apply(S,2,cumsumNA)
                            private$maxT <- length(x)
                            invisible(self)
                        },
                        #' @description Compute the non-anomalous cost of a segment
                        #' @param a start of period
                        #' @param b end of period
                        #' @param pen penalty cost
                        baseCost = function(a,b,pen=0){
                            a <- a-1
                            if(a<1){
                                sumStat <- private$summaryStats[b,]
                            }else{
                                sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                            }
                            2*sumStat[3] - 2*sumStat[4] + 2*sumStat[2] + pen
                        },
                        #' @description Compute the point anomaly cost of a time step
                        #' @param b time step
                        #' @param pen penalty cost
                        pointCost = function(b,pen){    
                            a <- b-1
                            if(a<1){
                                sumStat <- private$summaryStats[b,]
                            }else{
                                sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                            }
                            if( is.na(sumStat[1]) | (sumStat[1] == 0) ){
                                xlgx <- 0
                            }else{
                                xlgx <- sumStat[1]*log(sumStat[1])
                            }
                            2*(sumStat[1] - xlgx  + sumStat[2]) + pen
                        },
                        #' @description Compute the anomalous cost of a segment
                        #' @param a start of period
                        #' @param b end of period
                        #' @param pen penalty cost
                        #' @param len minimum number of observations
                        collectiveCost = function(a,b,pen,len){
                            a <- a-1
                            if(a<1){
                                sumStat <- private$summaryStats[b,]
                            }else{
                                sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                            }
                            if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                            
                            rhat <- sumStat[1] / sumStat[3]
                            rhat <- max(rhat, .Machine$double.eps)
                            
                            2*rhat*sumStat[3] - 2*sumStat[1]*log(rhat) - 2*sumStat[4] + 2*sumStat[2] + pen
                        },
                        #' @description Compute parameters of a segment if anomalous
                        #' @param a start of period
                        #' @param b end of period
                        param = function(a,b){
                            a <- a-1
                            if(a<1){
                                sumStat <- private$summaryStats[b,]
                            }else{
                                sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                            }
                            
                            rhat <- sumStat[1] / sumStat[3]
                            max(rhat, .Machine$double.eps)
                        }
                    )
                    )
