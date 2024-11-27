#' @export
gammaRepCost <- R6Class("gaussRepCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         initialize = function(x,m=1e-6,s=1){
                             if(is.list(x)){
                                 tmp <- matrix(NA,length(x),9)
                                 m <- rep(m,length(x))
                                 s <- rep(s,length(x))
                                 #browser()
                                 for(ii in 1:length(x)){
                                     tmp[ii,] <- c(length(x[[ii]]),
                                                   length(x[[ii]])/s[ii],
                                                   length(x[[ii]])*log(s[ii]),
                                                   sum(x[[ii]]-m[ii])/s[ii],
                                                   sum( (x[[ii]]-m[ii])^2 )/s[ii],
                                                   sum( x[[ii]] ),
                                                   sum( log(x[[ii]]) ),
                                                   m[ii],
                                                   s[ii]
                                                   )
                                 }
                                 self$summaryStats <- tmp
                             }else{
                                 self$summaryStats <- cbind(1,1/s,log(s),(x-m)/s,((x-m)^2)/s,x,log(x),m,s)
                             }
                             self$maxT <- length(x)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             n <- self$summaryStats[a:b,1]
                             mn <- self$summaryStats[a:b,8]
                             vr <- self$summaryStats[a:b,9]
                             Sx <- self$summaryStats[a:b,6]
                             Slx <- self$summaryStats[a:b,7]
                             scale <- vr/mn
                             shape <- (mn^2)/vr
                             
                             tmp <- (shape-1)*Slx - (Sx/scale) - n*lgamma(shape) - n*shape*log(scale)
                             -2*sum(tmp) + pen
                         },
                         pointCost = function(a,pen){
                             private$varChange(a,a,pen) 
                              
                         }
                     ),
                     private=list(
                         paramMeanChange = function(a,b){
                             mhat <- sum(self$summaryStats[a:b,4]) / sum(self$summaryStats[a:b,2])
                             mn <- self$summaryStats[a:b,8] + mhat
                             vr <- self$summaryStats[a:b,9]
                             vr <- max(vr, .Machine$double.xmin)
                             mn <- max(mn, .Machine$double.xmin)
                             return(list(m=mn,v=vr))
                         },
                         meanChange = function(a,b,pen){
                             p <- private$paramMeanChange(a,b)
                             n <- self$summaryStats[a:b,1]
                             Sx <- self$summaryStats[a:b,6]
                             Slx <- self$summaryStats[a:b,7]
                             scale <- p$v/p$m
                             shape <- (p$m^2)/p$v
                             
                             tmp <- (shape-1)*Slx - (Sx/scale) - n*lgamma(shape) - n*shape*log(scale)
                             -2*sum(tmp) + pen
                         },
                         paramVarChange = function(a,b){
                             shat <- sum(self$summaryStats[a:b,5]) / sum(self$summaryStats[a:b,1])
                             shat <- max(shat,.Machine$double.xmin)
                             mn <- self$summaryStats[a:b,8]
                             vr <- self$summaryStats[a:b,9] * shat
                             vr <- max(vr, .Machine$double.xmin)
                             mn <- max(mn, .Machine$double.xmin)
                             return(list(m=mn,v=vr))
                         },
                         varChange = function(a,b,pen){
                             p <- private$paramVarChange(a,b)
                             n <- self$summaryStats[a:b,1]
                             Sx <- self$summaryStats[a:b,6]
                             Slx <- self$summaryStats[a:b,7]
                             scale <- p$v/p$m
                             shape <- (p$m^2)/p$v
                             if( any(scale<0) ){ browser() }
                             tmp <- (shape-1)*Slx - (Sx/scale) - n*lgamma(shape) - n*shape*log(scale)
                             -2*sum(tmp) + pen
                         },
                         paramMeanVarChange = function(a,b){
                             mhat <- sum(self$summaryStats[a:b,4]) / sum(self$summaryStats[a:b,2])
                             
                             shat <- sum(self$summaryStats[a:b,5] - (mhat^2)*self$summaryStats[a:b,2] ) / sum(self$summaryStats[a:b,1])
                             shat <- max(shat,.Machine$double.xmin)
                             
                             mn <- self$summaryStats[a:b,8] + mhat
                             vr <- self$summaryStats[a:b,9] * shat
                             mn <- max(mn, .Machine$double.xmin)
                             vr <- max(vr, .Machine$double.xmin)
                             return(list(m=mn,v=vr))
                         },
                         meanVarChange = function(a,b,pen){
                             p <- private$paramMeanVarChange(a,b)
                             n <- self$summaryStats[a:b,1]
                             Sx <- self$summaryStats[a:b,6]
                             Slx <- self$summaryStats[a:b,7]
                             scale <- p$v/p$m
                             shape <- (p$m^2)/p$v
                             if( any(scale<0) ){ browser() }
                             tmp <- (shape-1)*Slx - (Sx/scale) - n*lgamma(shape) - n*shape*log(scale)
                             -2*sum(tmp) + pen
                         }
                     )
                     )
                     
#' @export
gammaRepMean <- R6Class("gaussMean",
                     inherit = gammaRepCost,
                     public = list(
                         collectiveCost = function(a,b,pen){ private$meanChange(a,b,pen) },
                         param = function(a,b){private$paramMeanChange(a,b) }
                     )
                     )

#' @export
gammaRepVar <- R6Class("gaussVar",
                    inherit = gammaRepCost,
                    public = list(
                        collectiveCost = function(a,b,pen){ private$varChange(a,b,pen) },
                        param = function(a,b){private$paramVarChange(a,b) }
                    )
                    )

#' @export
gammaRepMeanVar <- R6Class("gaussMeanVar",
                        inherit = gammaRepCost,
                        public = list(
                            collectiveCost = function(a,b,pen){ private$meanVarChange(a,b,pen) },
                            param = function(a,b){private$paramMeanVarChange(a,b) }
                        )
                        )



