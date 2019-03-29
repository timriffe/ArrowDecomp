# Author: tim
###############################################################################


# There is an essential question on what to do when an element
# is perturbed so as to hold constraint that all arrows originating
# in a state sum to 1. Can hold one constant, making the other(s) a
# residual. So far the residual is the thing left out of the decomp,
# namely self arrows. However, tests using different arrow types as
# residuals (below) make this potentially probematic. Need to do another
# round of testing where one item is perturbed and all remaining items
# constrained rather than just a single residual. But that would need
# to somehow happen inside horiuchi? There may be a functional solution
# but the problem at hand is mostly one of programming strategy, and will
# have to be woven into the guts of horiuchi (or maybe kitagawa better?).

source("/home/tim/git/HLEDecomp/HLEDecomp/Code/R/Functions.R")



f_dec_out <- function(outvec){
	datout  <- v2m(outvec, ntrans = 2)
	# remove death rates and replace with self-arrows, needed to make
	# transition matrices
	datself <- out2self(datout, ntrans = 2)
	sum(U2N(data_2_U(datself,2),1)[,c(1,6)] %*% c(.5,.5))
}


f_dec_self <- function(selfvec){
	datself  <- v2m(selfvec, ntrans = 2)
	colnames(datself) <- c("m11", "m12", "m21", "m22")
	# remove death rates and replace with self-arrows, needed to make
	sum(U2N(data_2_U(datself,2),1)[,c(1,6)] %*% c(.5,.5))
}

other2self <- function(datother,ntrans=2){
	# assumes datother has column names
	trans.self      <- getcols(ntrans = ntrans, self = TRUE)
	trans.other     <- colnames(datother)
	other.order     <- as.data.frame(matrix(trans.other, ntrans, byrow = TRUE),stringsAsFactors=FALSE)
	
	btwn           <- apply(other.order, 1, function(x, datother){
				1 - rowSums(datother[, x])
			}, datother = datother)
	colnames(btwn) <- paste0("m", 1:ntrans, ntrans:1)
	
	out             <- cbind(datother, btwn)[, trans.self]
	out
}
getcolsall <- function(ntrans,dead=ntrans+1){
	sort(paste0("m",outer(1:ntrans,c(1:ntrans,dead),paste0)))
}

f_dec_other <- function(othervec, names = c("m11", "m13","m22","m23")){
	datother           <- v2m(othervec, ntrans = 2)
	colnames(datother) <- names
	datself            <-other2self(datother,2)
	
	# remove death rates and replace with self-arrows, needed to make
	sum(U2N(data_2_U(datself,2),1)[,c(1,6)] %*% c(.5,.5))
}

#f_dec_self(selfvec)
#f_dec_out(outvec)
#f_dec_other(othervec)
# now t1
m12.1 <- c(.1,.2,.3,0); m21.1 <- c(.4,.2,.1,0)
m13.1 <- c(.1,.2,.4,1); m23.1 <- c(.3,.4,.8,1)
m11.1 <- 1 - c(m12.1 + m13.1);m22.1 <- 1 - c(m21.1 + m23.1)
# t2
m12.2 <- c(.1,.2,.3,0)/2; m21.2 <- c(.4,.2,.1,0)/2
m13.2 <- c(.08,.18,.3,1); m23.2 <- c(.4,.5,.8,1)
m11.2 <- 1 - c(m12.2 + m13.2);m22.2 <- 1 - c(m21.2 + m23.2)

outvec.1  <- c(m12.1,m13.1,m21.1,m23.1)
selfvec.1 <- c(m11.1,m12.1,m21.1,m22.1)

outvec.2  <- c(m12.2,m13.2,m21.2,m23.2)
selfvec.2 <- c(m11.2,m12.2,m21.2,m22.2)

othervec.1 <- c(m11.1, m13.1,m22.1,m23.1)
othervec.2 <- c(m11.2, m13.2,m22.2,m23.2)
# self vs out vs other identical
f_dec_self(selfvec.1);f_dec_out(outvec.1);f_dec_other(othervec.1)
# decrease in e0
f_dec_self(selfvec.2);f_dec_out(outvec.2);f_dec_other(othervec.2)

f_dec_self(selfvec.2)-f_dec_self(selfvec.1)


library(DemoDecomp)
s1 <- v2m(horiuchi(f_dec_self,selfvec.1,selfvec.2,200),2)
o1 <- v2m(horiuchi(f_dec_out,outvec.1,outvec.2,200),2)
oth1 <- v2m(horiuchi(f_dec_other,othervec.1,othervec.2,200),2)
colnames(s1) <- c("m11","m12","m21","m22")
colnames(o1) <- c("m12","m13","m21","m23")
colnames(oth1) <- c("m11","m13","m22","m23")

compareCols <- function(M,ntrans=2){
	allcols <- getcolsall(ntrans)
	MM <- matrix(nrow=nrow(M),ncol=length(allcols),dimnames=list(NULL,allcols))
	MM[,colnames(M)]<- M
	MM
}
colSums(compareCols(s1))
colSums(compareCols(o1))
colSums(compareCols(oth1))
