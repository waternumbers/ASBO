#' An implimentation of the CROPS algorithm in 1D
#' @param betaMin lower bound of penalisation window
#' @param betaMax upper bound of penalisation window
#' @param fCost the cost function
#' @param alg algorithm either capa of pelt
#' @param betaP penalty for adding a point anomaly - only for use with capa
#' @param min_length minimum number of values in a collective anomaly
#' @param prune logical, should pruning be used
#' @param verbose logical, print out progress
#' @param maxIter maximum number of algorithm evaluations to perform
#' @return something...
#'
#' @details This will only work for cost functions where the beta is additive!!!
#' @export
crops <- function(betaMin,betaMax,fCost,alg=pelt,betaP=Inf,min_length=2,prune=TRUE,verbose=FALSE,maxIter=100){
    
    betaRec <- c(betaMin,betaMax)
    qRec <- rep(NA,2)
    mRec <- rep(NA,2)
    outRec <- list()

    for(ii in 1:length(betaRec)){
        tmp <- alg(partition(betaRec[ii],betaP,min_length),fCost)
        mRec[ii] <- max(0, nrow(collective_anomalies(tmp)))
        qRec[ii] <- max(tmp$cost,na.rm=TRUE) - mRec[ii]*betaRec[ii]
        outRec[[paste(mRec[ii])]] <- tmp
    }
    
    intervals <- list( c(1,2) )
    if( all(mRec==mRec[1]) ){
        warning("Starting conditions produce the same solution")
        return(list(betaRec=betaRec,
                    qRec=qRec,
                    mRec=mRec,
                    outRec=outRec))
    }

    
    cnt <- 0
    while(length(intervals)>0 & cnt < maxIter){
        ii <- intervals[[1]][1]
        jj <- intervals[[1]][2]
        if( mRec[ii] > mRec[jj]+1 ){
            bint <- (qRec[jj] - qRec[ii]) / (mRec[ii] - mRec[jj])
            if(!is.finite(bint)){ browser() }
            if( bint %in% betaRec ){ break } ## seems to be needed but not sure why..
            tmp <- alg(partition(bint,betaP,min_length),fCost)
            kk <- nrow(collective_anomalies(tmp))
            if(kk!=mRec[jj]){
                betaRec <- c(betaRec,bint)
                mRec <- c(mRec,kk)
                qRec <- c(qRec,max(tmp$cost,na.rm=TRUE) - kk*bint)
                
                if(!(paste(kk) %in% names(outRec))){
                    outRec[[paste(kk)]] <- tmp
                }
                intervals <- c(intervals,
                               list(c(ii,length(betaRec))),
                               list(c(length(betaRec),jj)))
            }
            cnt <- cnt + 1
        }
        intervals <- intervals[-1]
        ##print(intervals)
    }
    if(cnt==maxIter){ warning("Maximum numbers of iterations reached") }
    ##return(outRec)
    list(betaRec=betaRec,
                 qRec=qRec,
                 mRec=mRec,
                 outRec=outRec)
}

    
