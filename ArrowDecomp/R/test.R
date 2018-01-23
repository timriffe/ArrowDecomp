
# Author: tim
###############################################################################

# rates and probabilities testing.
Dat <- local(get(load("/home/tim/git/DistributionTTD/DistributionTTD/Data/HMDltper.Rdata")))

mx2qx <- function(mx){
	mx / (1 + .5 * mx)
}

# cumprod of px
qx2lx <- function(qx){
	cumprod(1-c(0,qx))[1:length(qx)]
}

#
lx2dx <- function(lx){
	-diff(c(lx,0))
}

lx2Lx <- function(lx,closeout=0){
	(lx + c(lx[-1],closeout)) / 2
	
}

lx2ex <- function(lx){
	lx <- lx / lx[1]
	Lx <- lx2Lx(lx)
	sum(Lx)
}

mx2ex <- function(mx){
	qx <- mx2qx(mx)
	lx <- qx2lx(qx)
	lx2ex(lx)
}

qx2ex <- function(qx){
	lx <- qx2lx(qx)
	lx2ex(lx)
}

px2ex <- function(px){
	qx <- 1-px
	qx2ex(qx)
}

head(Dat)
mx1 <- Dat$mx[Dat$CNTRY == "USA" & Dat$Year == 2000 & Dat$Sex == "f"]
qx1 <- mx2qx(mx1)
px1 <- 1 - qx1

mx2 <- Dat$mx[Dat$CNTRY == "USA" & Dat$Year == 2000 & Dat$Sex == "m"]
qx2 <- mx2qx(mx2)
px2 <- 1 - qx2

library(DecompHoriuchi)

decmx <- DecompContinuousOrig(mx2ex,mx1,mx2, N=100)
decqx <- DecompContinuousOrig(qx2ex,qx1,qx2, N=100)
decpx <- DecompContinuousOrig(px2ex,px1,px2, N=100)

plot(decmx)
lines(decqx)
lines(decpx,col="red",lty=2)

# makes no difference, so exit-only decomp is a matter of interpretation.
#-------------------------------------------




