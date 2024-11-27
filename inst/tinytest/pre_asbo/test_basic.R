## basic tests of the correspondance to the theory in the univariate case

## testing point anomalies
expect_silent({
    x <- seq(0,5,length=1000)
    ## manual evaluation of the solution
    idx <- (x^2) > log( exp(-(1+2)) + x^2 ) + 1 + 2
    ans <- which(idx)
    ## various capa forma
    mv <- capa(partition(1e300,2,10),gaussMeanVar$new(x))
    m <- capa(partition(1e300,2,10),gaussMean$new(x))
})
expect_equal( point_anomalies(mv)$location, ans )
expect_equal( point_anomalies(m)$location, ans )
expect_equal( nrow(collective_anomalies(mv)), integer(1) )
expect_equal( nrow(collective_anomalies(m)), integer(1) )

## test for detection of collective anomalies with mean
expect_silent({
    n <- 20
    beta <- 2
    x <- rep(0,n) ## should be the same length as min_seg_len
    z <- c(x,x+sqrt(beta/n))
    res <- capa(partition(beta,100,10),gaussMean$new(z))
})
expect_equal( nrow(collective_anomalies(res)), integer(1) )
expect_equal( nrow(point_anomalies(res)), integer(1) )

expect_silent({
    z <- c(x,x+sqrt((beta+1e-10)/n))
    res <- capa(partition(beta,1e300,10),gaussMean$new(z))
    ca <- collective_anomalies(res)
})
expect_equal( c(ca$start,ca$end), c(21,40) )
expect_equal( nrow(point_anomalies(res)), integer(1) )

## TODO test for collective anomaly detection for meanvar
