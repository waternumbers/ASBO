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
capa <- function(part,fCost,prune = TRUE,verbose=FALSE){

    ## initialise the changes to part
    part$endPoint <- rep(NA,fCost$length())
    part$cost <- rep(NA,fCost$length())
    part$type <- rep(NA_character_,fCost$length())


    
    cnst <- max(part$beta, part$betaP,na.rm=T)
    
    endPoints <- c(0) ## end points to search over

    step <- ceiling( fCost$length() / 100 )
    for(tt in 1:fCost$length()){ ##fCost$validTimes){ ##maxT){
        
        if(verbose && (tt %% step == 0)) {
            ## Print on the screen some message
            cat(paste0("time step: ", tt, "\n"))
            cat(paste0("Number of endPoints: ", length(endPoints), "\n"))
            
        }
       
        endPointCosts <- rep(NA,length(endPoints)) ## costs at those end Points
        endPointType <- rep(NA,length(endPoints)) ## what type is it?

        ## can't evaluate so skip to next time step
        ##if( is.na( fCost$baseCost(tt,tt,0) ) ){ next }
        
        for(ii in seq_along(endPoints)){
            jj <- endPoints[ii]
            if(jj == 0){ jjCost <- 0 }else{ jjCost <- part$cost[jj] }            
            if( jj == tt-1 ){
                ## test adding as base or point anomaly
                ## compute C2 from paper
                endPointCosts[ii] <- jjCost + fCost$baseCost(tt,tt,0)
                endPointType[ii] <- "background"
                ## compute C3 from paper
                tmp <- jjCost + fCost$pointCost(tt,part$betaP)
                
                if(is.finite(tmp) & is.finite(endPointCosts[ii]) & tmp < endPointCosts[ii]){
                    endPointCosts[ii] <- tmp
                    endPointType[ii] <- "point"
                }
            }else{
                endPointCosts[ii] <- jjCost + fCost$collectiveCost(jj+1,tt,part$beta,part$min_length) 
                endPointType[ii] <- "collective"
            }
        }
        
        
        ## old catch for NA values
        if( all(is.na(endPointCosts)) ){ next } ## can't evaluate at tt

        ## find minimum
        idx <- which.min(endPointCosts) ## ignores NA
        if(length(idx)==0){ browser() }
        part$endPoint[tt] <- endPoints[idx]
        part$cost[tt] <- endPointCosts[idx]
        part$type[tt] <- endPointType[idx]
        
        ## apply merging of background periods - could do this in summary?    
        if( ( part$endPoint[tt] > 0 ) &&
            ( part$type[tt] == "background" ) &&
            ( part$type[ part$endPoint[tt] ] == "background" )
           ){
            part$endPoint[tt] <- part$endPoint[ part$endPoint[tt] ]
        }
        
        if(prune){
            idx <- is.na(endPointCosts) | (endPointCosts < part$cost[tt] + cnst) ## removed equals
            endPoints <- endPoints[ idx ]
        }

        endPoints <- c(endPoints, tt)
        
    }
    
    return(part)
}

