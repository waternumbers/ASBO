#' @export
poisCost <- R6Class("poisCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         initialize = function(x,rate=1){
                             S <- cbind(x,lfactorial(x),rate,x*log(rate),1)
                             S[is.na(x),] <- NA
                             self$summaryStats <- apply(S,2,cumsumNA)
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
                             2*sumStat[3] - 2*sumStat[4] + 2*sumStat[2] + pen
                         },
                         pointCost = function(b,pen){
                             a <- b-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             if( sumStat[1] == 0){
                                 xlgx <- 0
                             }else{
                                 xlgx <- sumStat[1]*log(sumStat[1])
                             }
                             2*(sumStat[1] - xlgx  + sumStat[2]) + pen
                         },
                         collectiveCost = function(a,b,pen){
                             a <- a-1
                             if(a<1){
                                 sumStat <- self$summaryStats[b,]
                             }else{
                                 sumStat <- self$summaryStats[b,] - self$summaryStats[a,]
                             }
                             rhat <- sumStat[1] / sumStat[3]
                             rhat <- max(rhat, .Machine$double.eps)
                             
                             2*rhat*sumStat[3] - 2*sumStat[1]*log(rhat) - 2*sumStat[4] + 2*sumStat[2] + pen
                         }
                     )
                    )
