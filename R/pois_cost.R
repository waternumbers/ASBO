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
#' p <- poisProp$new(x,r,0)
#' p$cost(90,95) ## cost of anomalous distribution for x[90:95]
#' @name poisCost

#' @rdname poisCost
poisCost <- R6Class("poisCost",
                    private = list(
                        summaryStats = NULL,
                        maxT = 0,
                        gamma = NULL
                    ),
                    public=list(
                        #' @description Get the length of time series
                        length = function(){ private$maxT },
                        #' @description Initialise the cost function
                        #' @param x numeric vector of observations
                        #' @param rate numeric vector of rate parameters
                        #' @param gamma penalty parameter
                        initialize = function(x,rate=1,gamma){
                            x <- as.numeric(x)
                            rate <- as.numeric(rate)
                            gamma <- as.numeric(gamma)
                            idx <- is.finite(rate)
                            stopifnot(
                                "zero or negative rates are not allowed" = all(rate[idx]>0),
                                "negative x values rates are not allowed" = all(x>=0)
                            )
                            S <- cbind(x,lfactorial(x),rate,x*log(rate),1)
                            S[is.na(x)|!idx,] <- NA
                            colnames(S) <- NULL ## else pollute the output
                            private$summaryStats <- apply(S,2,cumsumNA)
                            private$maxT <- length(x)
                            private$gamma <- gamma
                            invisible(self)
                        },
                        #' @description constant required in PELT algorithm
                        k = function(){ -private$gamma },
                        #' @description minimum length required for evaluation
                        min_length = function(){ 1 },
                        #' @description Compute parameters of a segment
                        #' @param a start of period
                        #' @param b end of period
                        param = function(a,b){ NA }
                    )
                    )



#' @description Compute the cost of a non-anomalous segment
#' @rdname poisCost
#' @export
poisBase = R6Class("poisBase",
                   inherit = poisCost,
                   public = list(
                       #' @description Compute the segment cost
                       #' @param a start of period
                       #' @param b end of period
                       cost = function(a,b){
                           a <- a-1
                           if(a<1){
                               sumStat <- private$summaryStats[b,]
                           }else{
                               sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                           }
                           2*sumStat[3] - 2*sumStat[4] + 2*sumStat[2] + private$gamma
                       })
                   )
                   ##     #' @description constant required in PELT algorithm
                   ##     k = function(){ -private$gamma },
                   ##     min_length = function(){ 1 },
                   ##     #' @description Compute parameters of a segment
                   ##     #' @param a start of period
                   ##     #' @param b end of period
                   ##     param = function(a,b){ NA }
                   ## )
                   ## )

#' @description Compute the cost of a anomalous segment with rate scaling parameter and constant cost
#' @rdname poisCost
#' @export
poisProp = R6Class("poisProp",
                   inherit = poisCost,
                   public = list(
                       #' @description Compute the segment cost
                       #' @param a start of period
                       #' @param b end of period
                       #' @param pen penalty cost
                       cost = function(a,b){
                           a <- a-1
                           if(a<1){
                               sumStat <- private$summaryStats[b,]
                           }else{
                               sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                           }
                                        #if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                           
                           rhat <- sumStat[1] / sumStat[3]
                           rhat <- max(rhat, .Machine$double.eps)
                           
                           2*rhat*sumStat[3] - 2*sumStat[1]*log(rhat) - 2*sumStat[4] + 2*sumStat[2] + private$gamma
                       },
                       ## #' @description constant required in PELT algorithm
                       ## k = function(){ -private$gamma },
                       ## min_length = function(){ 1 },
                       #' @description Compute parameters of a segment
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

#' @description Compute the cost of a anomalous segment with rate scaling parameter and cost proportional to expectation
#' @rdname poisCost
#' @export
poisPropG = R6Class("poisPropG",
                   inherit = poisCost,
                   public = list(
                       #' @description Compute the segment cost
                       #' @param a start of period
                       #' @param b end of period
                       #' @param pen penalty cost
                       cost = function(a,b){
                           a <- a-1
                           if(a<1){
                               sumStat <- private$summaryStats[b,]
                           }else{
                               sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                           }
                                        #if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                           
                           rhat <- sumStat[1] / sumStat[3]
                           rhat <- max(rhat, .Machine$double.eps)
                           
                           2*rhat*sumStat[3] - 2*sumStat[1]*log(rhat) - 2*sumStat[4] + 2*sumStat[2] + sumStat[3]*private$gamma
                       },
                       #' @description constant required in PELT algorithm
                       k = function(){ 0 },                       
                       ## min_length = function(){ 1 },
                       #' @description Compute parameters of a segment
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

