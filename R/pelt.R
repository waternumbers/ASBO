#' An R implimentation of the segmented search algorithmpelt algorithm
#' 
#' @param part the starting partition
#' @param fCost the cost function
#' @param min_length minimum segment length
#' @param prune logical, should pruning be used
#' @param verbose logical, print out progress
#'
#' @return the optimal partition
#' 
#' @details Basic R implimentation of pelt - not efficent
#' @export
pelt <- function(part,fCost,min_length=1,prune = TRUE,verbose=FALSE){

    k <- 0
    min_len <- min_length
    maxT <- 0
    for(fC in fCost){
        k <- max(k, fC$k())
        min_len <- max(min_len, fC$min_length())
        maxT <- max(maxT,fC$length())
    }
    
    
    ## initialise the changes to part
    part$endPoint <- rep(NA,maxT)
    part$cost <- rep(NA,maxT)
    part$type <- rep(NA_integer_,maxT)

    ## cnst and for use in pruning
    cnst <- part$beta + k

    endPoints <- c(0) ## end points to search over

    for(tt in 1:maxT){ ##fCost$validTimes){ ##maxT){
        if(verbose && (tt %% 100==0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }

        endPointCosts <- rep(NA,length(endPoints)) ## costs at those end Points
        endPointTypes <- rep(NA,length(endPoints)) ## type of segment from those end Points
        for(ii in seq_along(endPoints)){
            jj <- endPoints[ii]
            if(jj + min_len > tt){ next } ## apply min length

            if(jj == 0){ jjCost <- 0 }else{ jjCost <- part$cost[jj] }
            
            segCost <- sapply(fCost,function(fC){ fC$cost(jj+1,tt) }) ## segment costs
            endPointCosts[ii] <- min(segCost) + jjCost + part$beta
            endPointTypes[ii] <- which.min(segCost)
        }

        if( all(is.na(endPointCosts)) ){ next } ## can't evaluate at tt

        ## find minimum
        idx <- which.min(endPointCosts) ## ignores NA
        part$endPoint[tt] <- endPoints[idx]
        part$cost[tt] <- endPointCosts[idx]
        part$type[tt] <- endPointTypes[idx]

        if(prune){
            idx <- is.na(endPointCosts) | (endPointCosts <= part$cost[tt] + cnst) ## change <= to <
            endPoints <- endPoints[ idx ]
        }

        endPoints <- c(endPoints, tt)
        
    }
    
    return(part)
}

