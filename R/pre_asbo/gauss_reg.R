#' R6 class for Univariate Gaussian Cost Functions
#'
#' @description Cost functions for differening types of univariate gaussion anomalies
#' 
#' @details Collective anomalies are represented either as changes in parameters describing the mean (\code{gaussRegMean}),
#' variance (\code{gaussRegVar}) or mean and variance (\code{gaussRegMeanVar}). Changes in variance are represented as a scaling parameter, not changes to the covariance. See vignettes for details.
#'
#' Each element of the input \code{x} should be a list containing a vector of observations \code{y} and corresponding design matrix \code{X}. Optionally in can also include a vector of parameter \code{m} and covariance matrix \code{S}.
#' 
#' @examples
#' ## simple test
#' set.seed(10)
#' x <- list()
#' n <- 120
#' for(ii in 1:48){
#'     
#'     if(ii < 10){ theta = c(1,0); sigma <- 0.1 }
#'     if(ii >= 10 & ii <12){ theta <- c(10,0); sigma <- 2}
#'     if(ii >= 12 & ii < 44){ theta <- c(5,1); sigma <- 2}
#'     if(ii >= 44 ){ theta <- c(1,0); sigma <- 0.1}
#'     
#'     X <- cbind(rep(1,n),runif(n,ii-1,ii))
#'     y <- rnorm(n, X%*%theta, sigma)
#'     x[[ii]] <- list(y=y,X=X)
#' }
#'
#' fCost <- gaussRegMeanVar$new(x)
#' p <- partition(4*log(sum(sapply(y,length))),NA,2)
#' res <- pelt(p,fCost)
#' 
#' @name gaussRegCost
NULL

#' @rdname gaussRegCost
gaussRegCost <- R6Class(
    "gaussRegCost",
    public = list(
        #' @description Get the length of time series
        length = function(){ private$maxT },
        #' @description Initialise the cost function
        #' @param x a list of regressions (see details)
        #' @param non_neg should only non-negative paraemter solutions be considered
        initialize = function(x, non_neg=FALSE){
            stopifnot("x should be a list" = is.list(x))
            ## function to check a list
            fcheck <- function(z){
                stopifnot("y and X required" = all(c("y","X") %in% names(z)),
                          "y should be a numeric vector" = is.numeric(z$y),
                          "m should be a numeric vector" = !("m"%in% names(z)) || is.numeric(z$m),
                          "S should be a numeric matrix" = !("S"%in% names(z)) || (inherits(z$S,"matrix") &&  inherits(z$S,"numeric")),
                          "X should be a numeric matrix" = inherits(z$X,"matrix") &&  inherits(z$X[1],"numeric"),
                          "y and X should have the same number of observation" = length(z$y)==nrow(z$X),
                          
                          "m and X should be for the same number of parameters" = !("m"%in%names(z)) || length(z$m) == ncol(z$X),
                          "S is for the wrong observation dimension" = !("S"%in%names(z)) || all( dim(z$S)[1:2] == length(z$y))
                          )
                return(ncol(z$X))
            }
            #print("running fcheck")
            tmp <- sapply(x,fcheck)
            nt <- length(tmp)
            stopifnot("all X should have the same number of parameters" = length(unique(tmp))==1,
                      "ther must be more then one observation set" = nt>1)

            fprocess <- function(z){
                if("m" %in% names(z)){ z$y <- z$y - z$X %*% z$m }
                idx <- is.finite(z$y)
                if("S" %in% names(z)){
                    iS <- solve(z$S[idx,idx]); detS <- det(z$S[idx,idx])
                }else{ iS <- diag(sum(idx)); detS <- 1 }
                list(XtSX = t(z$X[idx,]) %*% iS %*% z$X[idx,],
                     XtSy = t(z$X[idx,]) %*% iS %*% z$y[idx],
                     ytSy = t(z$y[idx]) %*% iS %*% z$y[idx],
                     K = sum(idx)*log(2*pi) + log(detS),
                     n = sum(idx))
            }
            
            #print("f process")
            x <- lapply(x,fprocess)
            for(ii in 2:nt){
                for(jj in c("XtSX","XtSy","ytSy","K","n")){
                    x[[ii]][[jj]] <- x[[ii]][[jj]] + x[[ii-1]][[jj]]
                }
            }
            private$non_neg <- non_neg
            private$maxT <- nt
            private$summaryStats <- x
            invisible(self)
        },
        #' @description Compute the non-anomalous cost of a segment
        #' @param a start of period
        #' @param b end of period
        #' @param pen penalty cost
        baseCost = function(a,b,pen=0){
            a <- a-1
            if(a<1){
                K <- private$summaryStats[[b]]$K
                ytSy <- private$summaryStats[[b]]$ytSy
            }else{
                ytSy <- private$summaryStats[[b]]$ytSy - private$summaryStats[[a]]$ytSy
                K <- private$summaryStats[[b]]$K - private$summaryStats[[a]]$K
            }
            
            as.numeric( K + ytSy + pen)
        },
        #' @description Compute the point anomaly cost of a time step
        #' @param b time step
        #' @param pen penalty cost
        pointCost = function(b,pen){
            a <- b-1
            if(a<1){
                K <- private$summaryStats[[b]]$K
                ytSy <- private$summaryStats[[b]]$ytSy
                nk <- private$summaryStats[[b]]$n
            }else{
                ytSy <- private$summaryStats[[b]]$ytSy - private$summaryStats[[a]]$ytSy
                K <- private$summaryStats[[b]]$K - private$summaryStats[[a]]$K
                nk <- private$summaryStats[[b]]$n - private$summaryStats[[a]]$n
            }
            sigma <- max(1e-6, ytSy/nk)
            
            as.numeric( K + nk*log(sigma) + (ytSy/sigma) + pen )
        }
    ),
    private = list(
        summaryStats = NULL,
        maxT = 0,
        sigma0 = NA,
        non_neg = FALSE,
        meanChange = function(a,b,pen,len,param=FALSE){
            if( (b-1+1) < len ){ return(NA) }
            a <- a-1
            if(a<1){
                XtSX <- private$summaryStats[[b]]$XtSX
                XtSy <- private$summaryStats[[b]]$XtSy
                nk <- private$summaryStats[[b]]$n
                ytSy <- private$summaryStats[[b]]$ytSy
                K <- private$summaryStats[[b]]$K
            }else{
                XtSX <- private$summaryStats[[b]]$XtSX - private$summaryStats[[a]]$XtSX
                XtSy <- private$summaryStats[[b]]$XtSy - private$summaryStats[[a]]$XtSy
                nk <- private$summaryStats[[b]]$n - private$summaryStats[[a]]$n
                ytSy <- private$summaryStats[[b]]$ytSy - private$summaryStats[[a]]$ytSy
                K <- private$summaryStats[[b]]$K - private$summaryStats[[a]]$K
            }
            ## Check for non singular regression
            QR <- qr(XtSX)                 # Get the QR decomposition
            vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
            ## solve- set other values to 0
            theta <- rep(0,ncol(XtSX))
            if( private$non_neg ){
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
                K <- private$summaryStats[[b]]$K
                ytSy <- private$summaryStats[[b]]$ytSy
                nk <- private$summaryStats[[b]]$n
            }else{
                ytSy <- private$summaryStats[[b]]$ytSy - private$summaryStats[[a]]$ytSy
                K <- private$summaryStats[[b]]$K - private$summaryStats[[a]]$K
                nk <- private$summaryStats[[b]]$n - private$summaryStats[[a]]$n
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
                XtSX <- private$summaryStats[[b]]$XtSX
                XtSy <- private$summaryStats[[b]]$XtSy
                nk <- private$summaryStats[[b]]$n
                ytSy <- private$summaryStats[[b]]$ytSy
                K <- private$summaryStats[[b]]$K
            }else{
                XtSX <- private$summaryStats[[b]]$XtSX - private$summaryStats[[a]]$XtSX
                XtSy <- private$summaryStats[[b]]$XtSy - private$summaryStats[[a]]$XtSy
                nk <- private$summaryStats[[b]]$n - private$summaryStats[[a]]$n
                ytSy <- private$summaryStats[[b]]$ytSy - private$summaryStats[[a]]$ytSy
                K <- private$summaryStats[[b]]$K - private$summaryStats[[a]]$K
            }
            ## Check for non singular regression
            QR <- qr(XtSX)                 # Get the QR decomposition
            vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
            ## solve- set other values to 0
            theta <- rep(0,ncol(XtSX))
            if( private$non_neg ){
                theta[vars] <- quadprog::solve.QP(D=XtSX[vars,vars,drop=FALSE],d=XtSy[vars],A=diag(length(vars)),b=rep(0,length(vars)))$solution
            }else{
                theta[vars] <- solve(XtSX[vars,vars],XtSy[vars])
            }
            ##                             theta[vars] <- solve(XtSX[vars,vars],XtSy[vars])
            
            sigma <- (ytSy - t(XtSy) %*% theta)/nk
            
            ##if(sigma<=-1e-6){ browser() } ## TODO hit this to often...
            sigma <- max(1e-6,sigma)
            if(param){
                return( list(theta = theta,
                             sigma = sigma) )
            }
            
            as.numeric( K + nk*log(sigma) + nk + pen )
        }
    )
)

#' @rdname gaussRegCost      
#' @export
gaussRegMean <- R6Class("gaussRegMean",
                        inherit = gaussRegCost,
                        public = list(
                            #' @description Compute the anomalous cost of a segment
                            #' @param a start of period
                            #' @param b end of period
                            #' @param pen penalty cost
                            #' @param len minimum number of observations
                            collectiveCost = function(a,b,pen,len){ private$meanChange(a,b,pen,len) },
                            #' @description Compute parameters of a segment if anomalous
                            #' @param a start of period
                            #' @param b end of period
                            param = function(a,b){ private$meanChange(a,b,NA,0,TRUE) }
                            
                     )
                     )

#' @rdname gaussRegCost 
#' @export
gaussRegVar <- R6Class("gaussRegVar",
                       inherit = gaussRegCost,
                       public = list(
                           #' @description Compute the anomalous cost of a segment
                           #' @param a start of period
                           #' @param b end of period
                           #' @param pen penalty cost
                           #' @param len minimum number of observations
                           collectiveCost = function(a,b,pen,len){ private$varChange(a,b,pen,len) },
                           #' @description Compute parameters of a segment if anomalous
                           #' @param a start of period
                           #' @param b end of period
                           param = function(a,b){ private$varChange(a,b,NA,0,TRUE) }
                       )
                       )

#' @rdname gaussRegCost 
#' @export
gaussRegMeanVar <- R6Class("gaussRegMeanVar",
                           inherit = gaussRegCost,                           
                           public = list(
                               #' @description Compute the anomalous cost of a segment
                               #' @param a start of period
                               #' @param b end of period
                               #' @param pen penalty cost
                               #' @param len minimum number of observations
                               collectiveCost = function(a,b,pen,len){ private$meanVarChange(a,b,pen,len) },
                               #' @description Compute parameters of a segment if anomalous
                               #' @param a start of period
                               #' @param b end of period
                               param = function(a,b){ private$meanVarChange(a,b,NA,0,TRUE) }   
                           )
                           )

