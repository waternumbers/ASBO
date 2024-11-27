#' R6 class for the categorical distribution
#'
#' @description Cost functions for the multinomial distribution
#'
#' @details Collective anomalies are represented as chnages to the expected proportions.
#' Time varying expected proportions are currently not handled.
#' 
#' @examples
#' set.seed(0)
#' m <- c(1:4)/sum(1:4)
#' X <- t(rmultinom(100, 1, m))
#'
#' p <- categoricalCost$new(X,m)
#' p$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' p$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' p$collectiveCost(90,95,57,3)
#' @export
categoricalCost <- R6Class("categoricalCost",
                           private = list(
                               m = NA, ## expected proportions
                               summaryStats = NA,
                               maxT = 0
                           ),
                           public=list(
                               #' @description Get the length of time series
                               length = function(){ private$maxT },
                               #' @description Initialise the cost function
                               #' @param x integer matrix of observations of 0,1
                               #' @param m numeric vector of expected proportions
                               initialize = function(x,m=rep(1/ncol(x),ncol(x))){
                                   m <- as.numeric(m)
                                   stopifnot("x should be a matrix" = inherits(x,"matrix"),
                                             "x should be integer" = inherits(x[1],"integer"),
                                             "x should have row sums of 1" = all(rowSums(x)==1, na.rm=T),
                                             "negative x values rates are not allowed" = all(x>=0),
                                             "dimension of m does not match size of x" = length(m)==ncol(x)
                                             )
                                   #' ## need to check row sums are 1 - possibly just use multinomial??
                                   private$maxT <- nrow(x)
                                   S <- cbind(1,x)
                                   S[is.na(rowSums(S)),] <- NA
                                   private$summaryStats <- apply(S,2,cumsumNA)
                                   private$m <- log(m)
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
                                   -2*sum( sumStat[-1]*private$m ) + pen
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
                                   if( is.na(sumStat[1]) | sumStat[1]<len ){ return(NA) } ## check length and if NA
                                   
                                   n <- sum(sumStat[-1])
                                   lambda <- sumStat[-1]/n
                                   lambda[lambda==0] <- 1e-10 ## this is to get finite logs which are multiplied by 0
                                   -2*sum( sumStat[-1]*log(lambda) ) + pen
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
                                   n <- sum(sumStat[-1])
                                   sumStat/n
                               }
                           ))

                     
