#' R6 class for the multinomial distribution
#'
#' @description Cost functions for the multinomial distribution
#'
#' @details Collective anomalies are represented as chnages to the expected proportions.
#' Time varying expected proportions are currently not handled.
#' 
#' @examples
#' set.seed(0)
#' m <- c(1:4)/sum(1:4)
#' X <- t(rmultinom(100, 144, m))
#'
#' p <- multinomialCost$new(X,m)
#' p$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' p$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' p$collectiveCost(90,95,57,3)
#' @export
multinomialCost <- R6Class("multinomialCost",
                           private = list(
                               m = NA, ## expected proportions
                               summaryStats = NA,
                               maxT = 0
                           ),
                           public=list(
                               #' @description Get the length of time series
                               length = function(){ private$maxT },
                               #' @description Initialise the cost function
                               #' @param x integer matrix of observations
                               #' @param m numeric vector of expected proportions
                               initialize = function(x,m=rep(1/ncol(x),ncol(x))){
                                   m <- as.numeric(m)
                                   stopifnot("x should be a matrix" = inherits(x,"matrix"),
                                             "x should be integer" = inherits(x[1],"integer"),
                                             "negative x values are not allowed" = all(x>=0),
                                             "dimension of m does not match size of x" = length(m)==ncol(x)
                                             )
                                   private$maxT <- nrow(x)
                                   S <- cbind(x,
                                              lfactorial(rowSums(x)),
                                              rowSums(lfactorial(x)),
                                              1 )
                                   S[is.na(rowSums(S)),] <- NA
                                   private$summaryStats <- apply(S,2,cumsumNA)
                                   private$m <- m 
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
                                   nm <- length(private$m)
                                   -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[1:nm]*log(private$m)) ) + pen
                               },
                               #' @description Compute the point anomaly cost of a time step
                               #' @param b time step
                               #' @param pen penalty cost
                               pointCost = function(b,pen){
                                   self$collectiveCost(b,b,pen,1)
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
                                   nm <- length(private$m)
                                   if( is.na(sumStat[nm+3]) | sumStat[nm+3]<len ){ return(NA) } ## check length and if NA
                                   m <- sumStat[1:nm]/sum(sumStat[1:nm])
                                   idx <- (1:nm)[m>0]
                                   -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[idx]*log(m[idx])) ) + pen
                                   
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
                                   nm <- length(private$m)
                                   sumStat[1:nm]/sum(sumStat[1:nm])
                               }
                           ))

                     
