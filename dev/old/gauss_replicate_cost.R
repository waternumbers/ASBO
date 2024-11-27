#' @export
gaussRepCost <- R6Class("gaussRepCost",
                     public=list(
                         summaryStats = NULL,
                         nStep = NULL,
                         maxT = 0,
                         initialize = function(x,m=0,s=1){
                             if(is.list(x)){
                                 S <- matrix(NA,length(x),5)
                                 m <- rep(m,length(x))
                                 s <- rep(s,length(x))
                                 for(ii in 1:length(x)){
                                     n <- length(x[[ii]])
                                     if( n>0 ){
                                         S[ii,] <- c(n,
                                                     n/s[ii],
                                                     n*log(s[ii]),
                                                     sum(x[[ii]]-m[ii])/s[ii],
                                                     sum( (x[[ii]]-m[ii])^2 )/s[ii]
                                                     )
                                     }
                                 }
                             }else{
                                 S <- cbind(1,1/s,log(s),(x-m)/s,((x-m)^2)/s)
                             }
                             S[is.na(rowSums(S)),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
                             nS <- rep(1,nrow(S))
                             nS[is.na(S[,1])] <- NA
                             self$nStep <- cumsumNA(nS)
                             self$maxT <- length(x)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             
                             sumStat[1]*log(2*pi) + sumStat[3] + sumStat[5] + pen
                         },
                         pointCost = function(a,pen){
                             private$meanVarChange(a,a,pen,1)
                         }
                     ),
                     private=list(
                         meanChange = function(a,b,pen,len){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                                 nS <- self$nStep[b]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                 nS <- self$nStep[b] - self$nStep[a] + 1
                             }
                             if( is.na(nS) | nS<len ){ return(NA) } ## check length and if NA
                             
                             mhat <- sumStat[4] / sumStat[2]
                             sumStat[1]*log(2*pi) + sumStat[3] + sumStat[5] - (mhat^2)*sumStat[2] + pen
                         },
                         varChange = function(a,b,pen,len){
                             a <- a-1
                             n <- b-a
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                                 nS <- self$nStep[b]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                 nS <- self$nStep[b] - self$nStep[a] + 1
                             }
                             if( is.na(nS) | nS<len ){ return(NA) } ## check length and if NA
                             
                             shat <- sumStat[5] / sumStat[1]
                             shat <- max(shat,.Machine$double.xmin)
                             sumStat[1]*log(2*pi*shat) + sumStat[3] + sumStat[1] + pen
                         },
                         meanVarChange = function(a,b,pen,len){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                                 nS <- self$nStep[b]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                                 nS <- self$nStep[b] - self$nStep[a] + 1
                             }
                             if( is.na(nS) | nS<len ){ return(NA) } ## check length and if NA
                             
                             mhat <- sumStat[4] / sumStat[2]
                             shat <- (sumStat[5] - (mhat^2)*sumStat[2])/sumStat[1]
                             shat <- max(shat,.Machine$double.xmin)
                             
                             sumStat[1]*log(2*pi*shat) + sumStat[3] + sumStat[1] + pen
                         }
                     )
                     )
                     
#' @export
gaussRepMean <- R6Class("gaussMean",
                     inherit = gaussRepCost,
                     public = list(
                         collectiveCost = function(a,b,pen,len){ private$meanChange(a,b,pen,len) }
                     )
                     )

#' @export
gaussRepVar <- R6Class("gaussVar",
                    inherit = gaussRepCost,
                    public = list(
                        collectiveCost = function(a,b,pen,len){ private$varChange(a,b,pen,len) }
                    )
                    )

#' @export
gaussRepMeanVar <- R6Class("gaussMeanVar",
                        inherit = gaussRepCost,
                        public = list(
                            collectiveCost = function(a,b,pen,len){ private$meanVarChange(a,b,pen,len) }
                        )
                        )



