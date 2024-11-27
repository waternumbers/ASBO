## This function is a (simplified) R implementation of the bs()
## function in the splines library and illustrates how the Cox-de Boor
## recursion formula is used to construct B-splines.

# code taken from Xueheng folder's

# last edit: 8/Dec/2022

Bi <- function(x, degree, i=1, knots) {
  # computes one basis function
  if(degree == 0){
    B <- ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0)
  } else {
    if((knots[degree+i] - knots[i]) == 0) {
      alpha1 <- 0
    } else {
      alpha1 <- (x - knots[i])/(knots[degree+i] - knots[i])
    }
    if((knots[i+degree+1] - knots[i+1]) == 0) {
      alpha2 <- 0
    } else {
      alpha2 <- (knots[i+degree+1] - x)/(knots[i+degree+1] - knots[i+1])
    }
    B <- alpha1*Bi(x, (degree-1), i, knots) + alpha2*Bi(x, (degree-1), (i+1), knots)
  }
  return(B)
}

bs.raw <- function(x, degree=3, knots, domain=NULL){
  # x: points to be evaluated
  # knots: a vector representing knots locations
  # domain: the domain of the function
  if(missing(x)) stop("You must provide x")
  if(degree < 1) stop("The spline degree must be at least 1")
  if(is.null(domain)) domain <- range(x)
  knots <- sort(knots)
  K <- length(knots) - (degree + 1) 
  B.mat <- matrix(0,length(x),K)
  for(j in 1:K) B.mat[,j] <- Bi(x, degree, j, knots)
  if(any(x == knots[length(knots)])) B.mat[x == knots[length(knots)], K] <- 1
  return(B.mat)
}

#' @export
bs <- function(x, degree=3, interior.knots=NULL, Boundary.knots = c(0,1)) {
  if(missing(x)) stop("You must provide x")
  if(degree < 1) stop("The spline degree must be at least 1")
  Boundary.knots <- sort(Boundary.knots)
  interior.knots.sorted <- NULL
  if(!is.null(interior.knots)) interior.knots.sorted <- sort(interior.knots)
  knots <- c(rep(Boundary.knots[1], (degree+1)), interior.knots.sorted, rep(Boundary.knots[2], (degree+1)))
  K <- length(interior.knots) + degree + 1 
  B.mat <- matrix(0,length(x),K)
  for(j in 1:K) B.mat[,j] <- Bi(x, degree, j, knots)
  if(any(x == Boundary.knots[2])) B.mat[x == Boundary.knots[2], K] <- 1
  return(B.mat)
}

# Example ---------------------------------------------------------------------
# spline curve 
if(FALSE){
  # define knots to simulate spline with fda package
  library(fda)
  xlim <- c(0,18)
  int.knots <- c(2,4,6,9,11,15)
  knots0 = c(xlim[1],2,4,6,9,11,15,xlim[2])
  
  #define basis functions
  basis0 <- create.bspline.basis(rangeval = xlim,
                                 norder = 4,
                                 breaks = knots0)
  # define coefficients to simulate spline 
  plot(basis0)
  
  knots.test <- c(0,2,4,6,9)
  s <- seq(knots.test[1],knots.test[5], length.out=50)
  lines(s,Bi(x=s,degree = 3, knots = knots.test), type = 'l',col=2, lty=3, lwd=3)
  knots.test <- c(0,0,2,4,6)
  s <- seq(knots.test[1],knots.test[5], length.out=50)
  lines(s,Bi(x=s,degree = 3,knots = knots.test), type = 'l',col=3, lty=2, lwd=3)
}


