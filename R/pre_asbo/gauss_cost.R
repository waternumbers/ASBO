#' R6 class for Univariate Gaussian Cost Functions
#'
#' @description Cost functions for differening types of univariate gaussion anomalies
#'
#' @param x numeric vector or matrix of observations (see details)
#' @param m numeric vector of mean values
#' @param s numeric vector of standard deviation values
#' @param point_type representation of point anomalies as either a change in mean or variance
#' @param a start of period
#' @param b end of period
#' @param pen penalty cost
#' @param len minimum number of obseervations
#' 
#' @details Collective anomalies are represented either as changes in mean (\code{gaussMean}),
#' variance (\code{gaussVar}) or mean and variance (\code{gaussMeanvar}). See vignettes for details.
#' If x is a matrix then the values in each row are treated as IID replicate observations.
#' 
#' @examples
#' set.seed(0)
#' m <- runif(100)
#' s <- pmax(1e-4,runif(100))
#' x <- rnorm(100,m,s) ## example data
#'
#' gM <- gaussMean$new(x,m,s) ## anomalies are changes in mean
#' gM$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' gM$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' gM$collectiveCost(90,95,57,3)
#' 
#' gV <- gaussVar$new(x,m,s) ## anomalies are changes in variance
#' gV$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' gV$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' gV$collectiveCost(90,95,57,3)
#'
#' gMV <- gaussMeanVar$new(x,m,s) ## anomalies are changes in mean and variance
#' gMV$baseCost(90,95) ## cost of non-anomalous distribution for x[90:95]
#' gMV$pointCost(90,0) ## point anomaly cost for x[90] with 0 penalty
#' ## collective anomaly cost for x[90:95] with penalty of 57 and at least 3 observation
#' gMV$collectiveCost(90,95,57,3)
#' @name gaussCost
NULL

#' @rdname gaussCost
gaussCost <- R6Class("gaussCost",
                     public=list(
                         #' @description Get the length of time series
                         length = function(){ private$maxT },
                         #' @description Initialise the cost function
                         #' @param x numeric vector of observations
                         #' @param m numeric vector of mean values
                         #' @param s numeric vector of standard deviation values
                         #' @param point_type representation of point anomalies as either a change in mean or variance
                         initialize = function(x,m=0,s=1,point_type=c("var","mean")){
                             if(is.matrix(x)){
                                 stopifnot(
                                     "x should be numeric" = inherits(x[1],"numeric")
                                 )
                             }else{
                                 x <- matrix(as.numeric(x))
                             }
                             m <- as.numeric(m)
                             s <- as.numeric(s)
                             
                             n <- rowSums(is.finite(x))

                             private$point_type <- match.arg(point_type)
                             private$maxT <- length(n)
                             
                             S <- cbind( n/s,n*log(s),rowSums( (x-m)/s ,na.rm=T ),
                                        rowSums( ((x-m)^2)/s ,na.rm=T),n )
                             S[n==0,] <- NA
                             private$summaryStats <- apply(S,2,cumsumNA)
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
                             sumStat[5]*log(2*pi) + sumStat[2] + sumStat[4] + pen
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
                             gamma <- max(.Machine$double.xmin, exp(-(1+pen)))
                             mhat <- sumStat[3] / sumStat[1]
                             n <- 1
                             cst <- switch(
                                 private$point_type,
                                 "var" = log(2*pi) + sumStat[2] + log(gamma + sumStat[4]) + 1,
                                 "mean" = n*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1],
                                 stop("Unknown point anomaly type")
                             )
                             cst + pen
                         }
                     ),
                     private=list(
                         summaryStats = NULL,
                         point_type = "var",
                         maxT = 0,
                         meanChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             
                             mhat <- sumStat[3] / sumStat[1]
                             sumStat[5]*log(2*pi) + sumStat[2] + sumStat[4] - (mhat^2)*sumStat[1] + pen
                         },
                         varChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             shat <- sumStat[4] / sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[5]*log(2*pi*shat) + sumStat[2] + sumStat[5] + pen
                         },
                         meanVarChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             if( is.na(sumStat[5]) | sumStat[5]<len ){ return(NA) } ## check length and if NA
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[5]*log(2*pi*shat) + sumStat[2] + sumStat[5] + pen
                         },
                         paramMean = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             sumStat[3] / sumStat[1]
                         },
                         paramVar = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             shat <- sumStat[4] / sumStat[5]
                             max(shat,.Machine$double.xmin)
                         },
                         paramMeanVar = function(a,b){
                             a <- a-1
                             if(a<1){
                                 sumStat <- private$summaryStats[b,]
                             }else{
                                 sumStat <- private$summaryStats[b,] - private$summaryStats[a,]
                             }
                             mhat <- sumStat[3] / sumStat[1]
                             shat <- (sumStat[4] - (mhat^2)*sumStat[1])/sumStat[5]
                             shat <- max(shat,.Machine$double.xmin)
                             c(mhat,shat)
                         }

             

                     )
                     )

#' @rdname gaussCost                     
#' @export
gaussMean <- R6Class("gaussMean",
                     inherit = gaussCost,
                     public = list(
                         #' @description Compute the anomalous cost of a segment
                         #' @param a start of period
                         #' @param b end of period
                         #' @param pen penalty cost
                         #' @param len minimum number of observations
                         collectiveCost = function(a,b,pen,len){ private$meanChange(a,b,pen,len) },
                         #' @description Compute parameters of a segment if anomalous
                         #' @param a start of period
                         #' @param b end of period
                         param = function(a,b){ private$paramMean(a,b) }
                     )
                     )
#' @rdname gaussCost
#' @export
gaussVar <- R6Class("gaussVar",
                    inherit = gaussCost,
                    public = list(
                        #' @description Compute the anomalous cost of a segment
                        #' @param a start of period
                        #' @param b end of period
                        #' @param pen penalty cost
                        #' @param len minimum number of observations
                        collectiveCost = function(a,b,pen,len){ private$varChange(a,b,pen,len) },
                        #' @description Compute parameters of a segment if anomalous
                        #' @param a start of period
                        #' @param b end of period
                        param = function(a,b){ private$paramVar(a,b) }
                    )
                    )

#' @rdname gaussCost
#' @export
gaussMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussCost,
                        public = list(
                            #' @description Compute the non-anomalous cost of a segment
                            #' @param a start of period
                            #' @param b end of period
                            #' @param pen penalty cost
                            #' @param len minimum number of observations
                            collectiveCost = function(a,b,pen,len){ private$meanVarChange(a,b,pen,len) },
                            #' @description Compute parameters of a segment if anomalous
                            #' @param a start of period
                            #' @param b end of period
                            param = function(a,b){ private$paramMeanVar(a,b) }
                        )
                        )



