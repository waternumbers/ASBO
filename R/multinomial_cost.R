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
#' p <- multinomialBase$new(X,m,2)
#' p$cost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' @name multinomialCost

multinomialCost <- R6Class("multinomialCost",
                           private = list(
                               m = NA, ## expected proportions
                               summaryStats = NA,
                               maxT = 0,
                               gamma = NULL
                           ),
                           public=list(
                               #' @description Get the length of time series
                               length = function(){ private$maxT },
                               #' @description Initialise the cost function
                               #' @param x integer matrix of observations
                               #' @param m numeric vector of expected proportions
                               #' @param gamma penalty parameter
                               initialize = function(x,m=rep(1/ncol(x),ncol(x)),gamma){
                                   m <- as.numeric(m)
                                   gamma <- as.numeric(gamma)
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
                                   private$gamma <- gamma
                                   invisible(self)
                               },
                               #' @description constant required in PELT algorithm
                               k = function(){ -private$gamma },
                               #' @description minimum length of sample required
                               min_length = function(){ 1 },
                               #' @description Compute parameters of a segment
                               #' @param a start of period
                               #' @param b end of period
                               param = function(a,b){ NA }
                           )
                           )


#' @description Compute the non-anomalous cost of a segment
#' @rdname multinomialCost
#' @export
multinomialBase <- R6Class("multinomialBase",
                           inherit = multinomialCost,
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
                                   nm <- length(private$m)
                                   -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[1:nm]*log(private$m)) ) + private$gamma
                               }
                           )
                           )


#' @description Compute the segment cost for an anomalous segment with replaced rates
#' @rdname multinomialCost
#' @export
multinomialReplace <- R6Class("multinomialReplace",
                              inherit = multinomialCost,
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
                                      nm <- length(private$m)
                                      ##                                      if( is.na(sumStat[nm+3]) | sumStat[nm+3]<len ){ return(NA) } ## check length and if NA
                                      m <- sumStat[1:nm]/sum(sumStat[1:nm])
                                      idx <- (1:nm)[m>0]
                                      -2*( sumStat[nm+1] - sumStat[nm+2] + sum(sumStat[idx]*log(m[idx])) ) + private$gamma
                                  },
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
                                      nm <- length(private$m)
                                      sumStat[1:nm]/sum(sumStat[1:nm])
                               }
                           ))
                     
