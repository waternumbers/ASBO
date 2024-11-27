## TODO needs actual value tests
Data <- readRDS("example_days.rds")

s <- 1:24 - .5 # hrs at where data is observed
int.knots <- seq(2,22, by=3)
rangeval= c(0,24)

X <- bs(x=s,interior.knots=int.knots, Boundary.knots = rangeval)

fCost <- splineCost$new(Data,m=rep(0,ncol(X)),D=X)

expect_silent({ tmp <- fCost$collectiveCost(1,1,0) })
expect_silent({ tmp <- fCost$pointCost(1,0) })
expect_silent({ tmp <- fCost$baseCost(1,1) })
expect_silent({ tmp <- fCost$baseCost(1,4) })
expect_silent({ tmp <- fCost$collectiveCost(1,5,0) })
