## this code is adapted from that of tha anomaly package.
## Major changes are that:
## - ported from testhat to tinytest
## - multivariate tests are not used...
## - we try to second guess the defaults of the capa package
## - we have to allow for the differnt package output formats

target <- readRDS("trimmed_anomaly_test_objects.rds")

## ########################################################
## test_that("Example 1, 2, and 2a from vignettes",
set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)

## Example 1
expect_silent({
    fCost <- gaussMeanVar$new(x)
    p <- partition(4*log(length(x)), 3*log(length(x)), 10)
    res <- capa(p,fCost)
})
expect_equivalent( target$anomaly_paper_example_1_collective_anomalies[,c("start","end")],
                  collective_anomalies(res)[,c("start","end")] )
expect_equal( target$anomaly_paper_example_1_point_anomalies$location,
             point_anomalies(res)$location )

	
## Example 2
expect_silent({
    fCost <- gaussMean$new(x, point_type="mean")
    p <- partition(3*log(length(x)), 3*log(length(x)), 10)
    res <- capa(p,fCost)
})
expect_equivalent( target$anomaly_paper_example_2_collective_anomalies[,c("start","end")],
                  collective_anomalies(res)[,c("start","end")] )
expect_equal( target$anomaly_paper_example_2_point_anomalies$location,
             point_anomalies(res)$location ) # this fails

## Example 2a
expect_silent({
    fCost <- gaussMean$new(1 + 2*x, point_type="mean")
    p <- partition(3*log(length(x)), 3*log(length(x)), 10)
    res <- capa(p,fCost)
})
expect_equivalent( target$anomaly_paper_example_2_a_collective_anomalies[,c("start","end")],
                  collective_anomalies(res)[,c("start","end")] )
expect_equal( target$anomaly_paper_example_2_a_point_anomalies$location,
             point_anomalies(res)$location )

## ###############################
## ## "Example 4 from vignettes",
## ## this is quite slow so commented out 
## data("machinetemp")
## x <- (machinetemp$temperature - median(machinetemp$temperature)) / mad(machinetemp$temperature)

## ## Part 1
## expect_silent({
##     fCost <- gaussMean$new(x)
##     p <- partition(4*log(length(x)), 3*log(length(x)), 10)
##     res <- capa(p,fCost)
## })
## expect_equivalent( target$anomaly_paper_example_4_1_collective_anomalies[,c("start","end")],
##                   collective_anomalies(res)[,c("start","end")] )
## expect_equal( target$anomaly_paper_example_4_1_point_anomalies$location,
##              point_anomalies(res)$location )

## ## Part 2
## library(robustbase)
## n <- length(x)
## x.lagged <- matrix(c(x[1:(n - 1)], x[2:n]), n - 1, 2)
## phi <- robustbase::covMcd(x.lagged, cor = TRUE)$cor[1,2]
## inflated_penalty <- 3 * (1 + phi) / (1 - phi) * log(n)


## expect_silent({
##     fCost <- gaussMean$new(x)
##     p <- partition(inflated_penalty, inflated_penalty, 10)
##     res <- capa(p,fCost)
## })
## expect_equal( target$anomaly_paper_example_4_2_collective_anomalies[, c("start","end")],
##                   collective_anomalies(res)[,c("start","end")] )
## expect_equal( target$anomaly_paper_example_4_2_point_anomalies$location,
##              point_anomalies(res)$location )
             

## test_that("Example 5 from vignettes",
## {
## 	data(simulated)

## 	# Part 2
## 	{
## 		res <- capa(sim.data, type = "mean", min_seg_len = 2)

## 		expect_equal(anomaly_paper_example_5_1_result,res,ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_point_plot,plot(res),ignore_attr=TRUE)
##     	}

## 	# Part 2
## 	{
## 		beta <- 2 * log(ncol(sim.data):1)
## 		beta[1] <- beta[1] + 3 * log(nrow(sim.data))
## 		res<-capa(sim.data, type= "mean", min_seg_len = 2,beta = beta)

## 		expect_equal(anomaly_paper_example_5_2_result,res)
## 		expect_equal(anomaly_paper_example_5_2_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_2_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_2_point_plot,plot(res),ignore_attr=TRUE)

## 	}
## })


## test_that("Example 6 from vignettes",
## {
## 	set.seed(0)
## 	x1 <- rnorm(500)
## 	x2 <- rnorm(500)
## 	x3 <- rnorm(500)
## 	x4 <- rnorm(500)
## 	x1[151:200] <- x1[151:200] + 2
## 	x2[171:200] <- x2[171:200] + 2
## 	x3[161:190] <- x3[161:190] - 3
## 	x1[351:390] <- x1[371:390] + 2
## 	x3[351:400] <- x3[351:400] - 3
## 	x4[371:400] <- x4[371:400] + 2
## 	x4[451] <- x4[451] * max(1, abs(1 / x4[451])) * 6
## 	x4[100] <- x4[100] * max(1, abs(1 / x4[100])) * 6
## 	x2[050] <- x2[050] * max(1, abs(1 / x2[050])) * 6
## 	x1 <- (x1 - median(x1))/mad(x1)
## 	x2 <- (x2 - median(x2))/mad(x2)
## 	x3 <- (x3 - median(x3))/mad(x3)
## 	x4 <- (x4 - median(x4))/mad(x4)
## 	x <- cbind(x1, x2, x3, x4)
## 	res <- capa(x, max_lag = 20, type = "mean")

## 	expect_equal(anomaly_paper_example_6_result,res,ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_point_plot,plot(res),ignore_attr=TRUE)
## })








