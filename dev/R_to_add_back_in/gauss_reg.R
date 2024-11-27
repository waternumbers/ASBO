#' @export
gaussRegCost <- R6Class("gaussRegCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         theta0 = NA,
                         sigma0 = NA,
                         non_neg = FALSE,
                         initialize = function(x,m=NULL,S=NULL, non_neg=FALSE){
                             self$summaryStats <- list(XtSX=list(),XtSy=list(),
                                                       ytSy=rep(NA,length(x)),
                                                       n=rep(NA,length(x)),
                                                       K=rep(NA,length(x)))
                             nx <- ncol(x[[1]]$X)
                             XtSX <- matrix(0,nx,nx)
                             XtSy <- matrix(0,nx,1)
                             ytSy <- 0
                             n <- 0
                             K <- 0

                             if(!is.null(m)){
                                 ## TODO check inputs
                             }
                             if(!is.null(S)){
                                 ## TODO check inputs
                             }
                             
                             for(ii in 1:length(x)){
                                 if(is.null(m)){
                                     yhat <- x[[ii]]$y
                                 }else{
                                     yhat <- x[[ii]]$y - x[[ii]]$X %*% m[[ii]]
                                 }
                                 
                                 nt <- length(yhat)
                                 if(is.null(S)){
                                     iS <- diag(nt)
                                     detS <- 1
                                 }else{
                                     iS <- solve(S[[ii]])
                                     detS <- det(S[[ii]])
                                 }

                                 ##print(ii)
                                 XtSX <- XtSX + ( t(x[[ii]]$X) %*% iS %*% x[[ii]]$X )
                                 XtSy <- XtSy + ( t(x[[ii]]$X) %*% iS %*% yhat )
                                 ytSy <- ytSy + ( t(yhat) %*% iS %*% yhat )
                                 K <- K + nt*log(2*pi) + log(detS) ##nt*log(2*pi) + detS
                                 n <- n + nt
                                 
                                 self$summaryStats$XtSX[[ii]] <- XtSX
                                 self$summaryStats$XtSy[[ii]] <- XtSy
                                 self$summaryStats$n[ii] <- n
                                 self$summaryStats$ytSy[ii] <- ytSy
                                 self$summaryStats$K[ii] <- K
                             }
                             self$non_neg <- non_neg
                             self$maxT <- length(x)
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             a <- a-1
                             if(a<1){
                                 K <- self$summaryStats$K[b]
                                 ytSy <- self$summaryStats$ytSy[b]
                             }else{
                                 ytSy <- self$summaryStats$ytSy[b] - self$summaryStats$ytSy[a]
                                 K <- self$summaryStats$K[b] - self$summaryStats$K[a]
                             }

                             as.numeric( K + ytSy + pen)
                         },
                         pointCost = function(b,pen){
                             a <- b-1
                             if(a<1){
                                 K <- self$summaryStats$K[b]
                                 ytSy <- self$summaryStats$ytSy[b]
                                 nk <- self$summaryStats$n[b]
                             }else{
                                 ytSy <- self$summaryStats$ytSy[b] - self$summaryStats$ytSy[a]
                                 K <- self$summaryStats$K[b] - self$summaryStats$K[a]
                                 nk <- self$summaryStats$n[b] - self$summaryStats$n[a]
                             }
                             sigma <- max(1, ytSy/nk)
                             
                             as.numeric( K + nk*log(sigma) + (ytSy/sigma) + pen )
                         }
                     ),
                     private = list(
                         meanChange = function(a,b,pen,len,param=FALSE){
                             if( (b-1+1) < len ){ return(NA) }
                             a <- a-1
                             if(a<1){
                                 XtSX <- self$summaryStats$XtSX[[b]]
                                 XtSy <- self$summaryStats$XtSy[[b]]
                                 nk <- self$summaryStats$n[b]
                                 ytSy <- self$summaryStats$ytSy[b]
                                 K <- self$summaryStats$K[b]
                             }else{
                                 XtSX <- self$summaryStats$XtSX[[b]] - self$summaryStats$XtSX[[a]]
                                 XtSy <- self$summaryStats$XtSy[[b]] - self$summaryStats$XtSy[[a]]
                                 nk <- self$summaryStats$n[b] - self$summaryStats$n[[a]]
                                 ytSy <- self$summaryStats$ytSy[b] - self$summaryStats$ytSy[[a]]
                                 K <- self$summaryStats$K[b] - self$summaryStats$K[[a]]
                             }
                             ## Check for non singular regression
                             QR <- qr(XtSX)                 # Get the QR decomposition
                             vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
                             ## solve- set other values to 0
                             theta <- rep(0,ncol(XtSX))
                             if( self$non_neg ){
                                 theta[vars] <- quadprog::solve.QP(D=XtSX,d=XtSy,A=diag(ncol(XtSX)),b=rep(0,ncol(XtSX)),meq=0)$solution
                             }else{
                                 theta[vars] <- solve(XtSX[vars,vars],XtSy[vars])
                             }


                             
                             if(param){
                                 return( list(theta = theta,
                                              sigma = 1) )
                             }

                             as.numeric(K + ytSy - (t(XtSy) %*% theta) +pen)
                             
                         },
                         varChange = function(a,b,pen,len,param=FALSE){
                             if( (b-1+1) < len ){ return(NA) }
                             a <- a-1
                             if(a<1){
                                 K <- self$summaryStats$K[b]
                                 ytSy <- self$summaryStats$ytSy[b]
                                 nk <- self$summaryStats$n[b]
                             }else{
                                 ytSy <- self$summaryStats$ytSy[b] - self$summaryStats$ytSy[a]
                                 K <- self$summaryStats$K[b] - self$summaryStats$K[a]
                                 nk <- self$summaryStats$n[b] - self$summaryStats$n[a]
                             }
                             
                             if(param){
                                 return( list(theta = matrix(0,ncol(XtSX),1),
                                              sigma = ytSy/nk) )
                             }
                             
                             as.numeric( K + nk*log(sigma) + nk + pen )
                         },
                         meanVarChange = function(a,b,pen,len,param=FALSE){
                             if( (b-1+1) < len ){ return(NA) }
                             a <- a-1
                             if(a<1){
                                 XtSX <- self$summaryStats$XtSX[[b]]
                                 XtSy <- self$summaryStats$XtSy[[b]]
                                 nk <- self$summaryStats$n[b]
                                 ytSy <- self$summaryStats$ytSy[b]
                                 K <- self$summaryStats$K[b]
                             }else{
                                 XtSX <- self$summaryStats$XtSX[[b]] - self$summaryStats$XtSX[[a]]
                                 XtSy <- self$summaryStats$XtSy[[b]] - self$summaryStats$XtSy[[a]]
                                 nk <- self$summaryStats$n[b] - self$summaryStats$n[[a]]
                                 ytSy <- self$summaryStats$ytSy[b] - self$summaryStats$ytSy[[a]]
                                 K <- self$summaryStats$K[b] - self$summaryStats$K[[a]]
                             }
                             ## Check for non singular regression
                             QR <- qr(XtSX)                 # Get the QR decomposition
                             vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
                             ## solve- set other values to 0
                             theta <- rep(0,ncol(XtSX))
                             if( self$non_neg ){
                                 ##browser()
                                 theta[vars] <- quadprog::solve.QP(D=XtSX[vars,vars,drop=FALSE],d=XtSy[vars],A=diag(length(vars)),b=rep(0,length(vars)))$solution
                             }else{
                                 theta[vars] <- solve(XtSX[vars,vars],XtSy[vars])
                             }
##                             theta[vars] <- solve(XtSX[vars,vars],XtSy[vars])
                             
                             sigma <- (ytSy - t(XtSy) %*% theta)/nk

                             ##if(sigma<=-1e-6){ browser() } ## TODO hit this to often...
                             sigma <- max(0,sigma)
                             if(param){
                                 return( list(theta = theta,
                                              sigma = sigma) )
                             }

                             as.numeric( K + nk*log(sigma) + nk + pen )
                         }
                     )
                     )

#' @export
gaussRegMean <- R6Class("gaussRegMean",
                     inherit = gaussRegCost,
                     public = list(
                         collectiveCost = function(a,b,pen,len){ private$meanChange(a,b,pen,len) },
                         param = function(a,b){ private$meanChange(a,b,NA,0,TRUE) }

                     )
                     )

#' @export
gaussRegVar <- R6Class("gaussRegVar",
                    inherit = gaussRegCost,
                    public = list(
                        collectiveCost = function(a,b,pen,len){ private$varChange(a,b,pen,len) },
                        param = function(a,b){ private$varChange(a,b,NA,0,TRUE) }
                    )
                    )

#' @export
gaussRegMeanVar <- R6Class("gaussRegMeanVar",
                        inherit = gaussRegCost,
                        public = list(
                            collectiveCost = function(a,b,pen,len){ private$meanVarChange(a,b,pen,len) },
                            param = function(a,b){ private$meanVarChange(a,b,NA,0,TRUE) }

                        )
                        )

