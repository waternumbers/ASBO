rm(list=ls())
devtools::load_all()
set.seed(123)
x <- rnorm(100,0,0.1)
x[1:10] <- rnorm(10,0.5,0.1)
x[30:40] <- rnorm(11,-0.5,0.1)
fCost <- ladCost$new(x)
tmp <- crops(1,100,fCost,alg=capa,betaP=3)
plot(tmp$mRec,tmp$qRec)
tmp$outRec[["2"]]



chck <- function(u,tau){u*(tau-(u<0))}
u <- rnorm(100)
md <- median(u) + seq(-0.1,0.1,length=1000)
sc <- matrix(NA,length(md),2)
for(ii in 1:length(md)){
    sc[ii,1] <- sum(chck(u-md[ii],0.5))
    sc[ii,2] <- mad(u,center=md[ii],constant=1)
}
sc <- apply(sc,2,function(x){ (x-min(x))/(max(x)-min(x)) })
matplot(md,sc,type="l")
