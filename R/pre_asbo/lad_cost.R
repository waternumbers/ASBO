#' R6 class for LAD Cost Functions
#'
#' @description Cost functions for the Least Absolution Deviation from a given quantile
#'
#' @details this is a very niaive and slow implimentation
#' 
#' @examples
#' set.seed(0)
#' m <- runif(100)
#' x <- rnorm(100,m)
#'
#' p <- ladCost$new(x,m,0.5)
#' p$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' p$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' p$collectiveCost(90,95,57,3)
#' @export
ladCost <- R6Class("ladCost",
                   private = list(
                       summaryStats = NULL,
                       tau = NA,
                       maxT = 0,
                       chck = function(x){x*(private$tau - (x<0))}
                   ),
                   public=list(
                       #' @description Get the length of time series
                       length = function(){ private$maxT },
                       #' @description Initialise the cost function
                       #' @param x numeric vector of observations
                       #' @param m expected value of x
                       #' @param tau the quantile
                       initialize = function(x,m=0,tau=0.5){
                           x <- as.numeric(x)
                           m <- as.numeric(m)
                           tau <- as.numeric(tau)[1]
                           stopifnot(
                               "tau should be between 0 & 1" = is.finite(tau) && tau>0 && tau <1
                           )
                           private$summaryStats <- x-m
                           private$maxT <- length(x)
                           private$tau <- tau
                           invisible(self)
                       },
                       #' @description Compute the non-anomalous cost of a segment
                       #' @param a start of period
                       #' @param b end of period
                       #' @param pen penalty cost
                       baseCost = function(a,b,pen=0){
                           sumStat <- private$summaryStats[a:b]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           2*( sum(private$chck(sumStat)) + length(sumStat)*log(private$tau*(1-private$tau)) ) + pen
                       },
                       #' @description Compute the point anomaly cost of a time step
                       #' @param a time step
                       #' @param pen penalty cost
                       pointCost = function(a,pen){ 
                           sumStat <- private$summaryStats[a]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           theta <- sumStat
                           2*( sum(private$chck(sumStat-theta)) + length(sumStat)*log(private$tau*(1-private$tau)) ) + pen
                           
                       },
                       #' @description Compute the anomalous cost of a segment
                       #' @param a start of period
                       #' @param b end of period
                       #' @param pen penalty cost
                       #' @param len minimum number of observations
                       collectiveCost = function(a,b,pen,len){
                           if( (b-1+1)<len ){ return(NA) } ## check length and if NA
                           sumStat <- private$summaryStats[a:b]
                           sumStat <- sumStat[ is.finite(sumStat) ]
                           if( length(sumStat) == 0 ){ return(Inf) }
                           if( private$tau == 0.5 ){ theta <- fast_med(sumStat) } ##median(sumStat) }
                           else{ theta <- quantile(sumStat,private$tau) }
                           2*( sum(private$chck(sumStat-theta)) + length(sumStat)*log(private$tau*(1-private$tau)) ) + pen
                       },
                       #' @description Compute parameters of a segment if anomalous
                       #' @param a start of period
                       #' @param b end of period
                       param = function(a,b){
                           quantile( private$summaryStats[a:b],tau )
                       }
                       
                   )
                   )

