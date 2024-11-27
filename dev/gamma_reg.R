#' @export
gammaRegCost <- R6Class("gammaRegCost",
                     public=list(
                         summaryStats = NULL,
                         maxT = 0,
                         theta0 = NA,
                         sigma0 = NA,
                         initialize = function(x,theta0=NULL,sigma0=1){
                             self$summaryStats <- list(y=lapply(x,function(x){x$y}),
                                                       X=lapply(x,function(x){x$X}))
                             self$maxT <- length(x)
                             nx <- ncol(x[[1]]$X)
                             if(is.null(theta0)){ theta0 <- rep(0,nx) }
                             if(length(theta0) != nx){ stop("theta0 incorrect size") }
                             self$theta0 <- theta0
                             self$sigma0 <- sigma0
                             invisible(self)
                         },
                         baseCost = function(a,b,pen=0){
                             X <- do.call(rbind,self$summaryStats$X[a:b])
                             y <- unlist( self$summaryStats$y[a:b] )
                             m <- X %*% self$theta0
                             v <- self$sigma0
                             sum( dgamma(y, shape = (m^2)/v, scale = v/m, log=TRUE) ) + pen
                         },
                         pointCost = function(a,pen){
                             self$collectiveCost(a,a,pen)
                         },
                         collectiveCost = function(a,b,pen){
                             p <- self$param(a,b)
                             y <- unlist( self$summaryStats$y[a:b] )
                             sum( dgamma(y, shape = (p$m^2)/p$v, scale = p$v/p$m, log=TRUE) ) + pen
                         },
                         param = function(a,b){
                             X <- do.call(rbind,self$summaryStats$X[a:b])
                             y <- unlist( self$summaryStats$y[a:b] )

                             ## Check for non singular regressio
                             QR <- qr(crossprod(X))                 # Get the QR decomposition
                             vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers that are OK to include
                             XX <- X[,vars,drop=FALSE]
                             strt <- solve(crossprod(XX),crossprod(XX,y))
                             ##browser()
                             strt <- rep(1,ncol(XX))
                             mdl <- glm(y~XX-1,family = Gamma(link = "identity"),start = as.numeric(strt))
                             ##browser()
                             theta <- rep(0,ncol(X))
                             theta[vars] <- coef( mdl )
                             mn = predict(mdl, type = "response")
                             vr = sum(residuals(mdl, type = "response")^2)/mdl$df.residual
                             list(theta=theta, v = vr, m=mn)
                         }
                     )
                     )

