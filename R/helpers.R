cumsumNA <- function(x){
    idx <- is.na(x)
    x[idx] <- 0
    out <- cumsum(x)
    out[idx] <- NA
    return(out)
}
