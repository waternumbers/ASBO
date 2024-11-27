#' An R implimentation of the segmented search algorithmpelt algorithm
#' 
#' @param part the starting partition
#' @param fCost the cost function
#' @param prune logical, should pruning be used
#' @param verbose logical, print out progress
#'
#' @return the optimal partition
#' 
#' @details Basic R implimentation of pelt - not efficent
#' @export
pelt <- function(part,fCost,prune = TRUE,verbose=FALSE){

    ## initialise the changes to part
    part$endPoint <- rep(NA,fCost$length())
    part$cost <- rep(NA,fCost$length())
    part$type <- rep(NA_character_,fCost$length())


    cnst <- part$beta

    endPoints <- c(0) ## end points to search over

    for(tt in 1:fCost$length()){ ##fCost$validTimes){ ##maxT){
        if(verbose && (tt %% 100==0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
        }

        endPointCosts <- rep(NA,length(endPoints)) ## costs at those end Points
        for(ii in seq_along(endPoints)){
            jj <- endPoints[ii]
            if(jj==0){ endPointCosts[ii] <- fCost$collectiveCost(jj+1,tt,part$beta,part$min_length) }
            else{ endPointCosts[ii] <- part$cost[ jj ] + fCost$collectiveCost(jj+1,tt,part$beta,part$min_length) }
        }

        if( all(is.na(endPointCosts)) ){ next } ## can't evaluate at tt

        ## find minimum
        idx <- which.min(endPointCosts) ## ignores NA
        part$endPoint[tt] <- endPoints[idx]
        part$cost[tt] <- endPointCosts[idx]
        part$type[tt] <- "collective"

        if(prune){
            idx <- is.na(endPointCosts) | (endPointCosts <= part$cost[tt] + cnst) ## change <= to <
            endPoints <- endPoints[ idx ]
        }

        endPoints <- c(endPoints, tt)
        
    }
    
    return(part)
}

