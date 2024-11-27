#' Get the parameters for a partitioning result
#' @param res the result of a partitioning algorithm
#' @param fCost the cost function
#'
#' @return list of parameters
#' 
#' @details Not yet implimented for all cost functions
#' @export
param <- function(res,fCost){

    if( !("param" %in%  names(fCost)) ){ stop("Parameter method not available for cost function") }
    out <- summary(res)
    tmp <- list()
    for(ii in 1:nrow(out)){
        tmp[[ii]] <- fCost$param(out$start[ii], out$end[ii])
    }
    tmp <- do.call(rbind,tmp)
    tmp[out$type=="background",] <- NA
    cbind(out,tmp)
}
