#' Test computations of the cost
#'
#' Direct computation of the cost of a  partition for IID gaussian variables, with a background distribution of N(0,1)
#' 
#' @param x univariate time series
#' @param ca data frame of collective anomalies, must contain columns called start and end
#' @param pa data frame of point anomalies, must contain column called location
#' @param beta penalty for introducing a collective anomaly
#' @param betaP penalty for introducing a point anomaly
#' @param type type of gaussian anomaly
#'
#' @return total cost
#'
#' @details The point anomaly evaluation does not include the correction term gamma
#' 
#' @export
actual_cost <- function(x,ca,pa,beta,betaP,type=c("meanvar","mean")){
    ## hanlde point anomalies
    pdx <- NULL
    cost_point <- 0
    if(nrow(pa) > 0){
        pdx <- pa$location
        cost_point <- -2*dnorm(x[pdx],0,abs(x[pdx]),log=TRUE) + betaP
    }    
    ## handle collective anomalies
    cdx <- list()
    cost_col <- 0
    if(nrow(ca) > 0){
        cdx <- lapply(1:nrow(ca),function(ii){ca$start[ii]:ca$end[ii]})
        cdx <- lapply(cdx,function(idx){setdiff(idx,pdx)}) ## remove point anomalies
        fc <- switch(type,
                     mean = function(idx){sum(-2*dnorm(x[idx],mean(x[idx]),1,log=TRUE)) + beta},
                     meanvar = function(idx){
                         m <- mean(x[idx]); v <- mean( (x-m)^2 )
                         sum(-2*dnorm(x[idx],m,sqrt(v),log=TRUE)) + beta
                     },
                     stop("unknown type")
                     )
      cost_col <- sapply(cdx,fc)
    }
    ## base times
    bdx <- 1:length(x)
    bdx <- setdiff(bdx,pdx)
    bdx <- setdiff(bdx,unlist(cdx))
    cost_base <- -2*sum(dnorm(x[bdx],0,1,log=TRUE))
    return( sum(c(cost_base,cost_point,cost_col))
           ##c(cost_base,sum(cost_point),sum(cost_col))
           )
}
