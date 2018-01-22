
# Author: tim
###############################################################################

# redo POC to remove self-arrows and include death transitions in decomp:

setwd("/home/tim/workspace/ArrowDecomp")



getsub <- function(Mat){
	Mat[row(Mat) == (col(Mat) + 1)]
}
extractfrom <- function(A){
	todead    <- A[nrow(A), ]
	
	# get transitions to death
	W2D       <- todead[1:49]
	I2D       <- todead[51:99]
	R2D       <- todead[101:149]
		
	A         <- A[1:150, 1:150]
	# extract positive vectors of px:
	rowBlocks <- ((row(A)-1) - (row(A)-1) %% 50) / 50 + 1
	colBlocks <- ((col(A)-1) - (col(A)-1) %% 50) / 50 + 1
	
	pivecs <- matrix(NA, ncol = 9, nrow = 49)
	

	W2I                	<- A[rowBlocks == 2 & colBlocks == 1]
	W2R                	<- A[rowBlocks == 3 & colBlocks == 1]
	I2W                	<- A[rowBlocks == 1 & colBlocks == 2]
	I2R                	<- A[rowBlocks == 3 & colBlocks == 2]
	R2W                	<- A[rowBlocks == 1 & colBlocks == 3]
	R2I                	<- A[rowBlocks == 2 & colBlocks == 3]
	
	# dims to extract diag
	dim(W2I) 			<- c(50,50)
	dim(W2R) 			<- c(50,50)
	dim(I2W) 			<- c(50,50)
	dim(I2R) 			<- c(50,50)
	dim(R2W) 			<- c(50,50)
	dim(R2I) 			<- c(50,50)
	
	W2I              	<- getsub(W2I)
	W2R              	<- getsub(W2R)
	I2W              	<- getsub(I2W)
	I2R              	<- getsub(I2R)
	R2W              	<- getsub(R2W)
	R2I              	<- getsub(R2I)
	# remix:
	Wfrom <- cbind(W2I, W2R, W2D) # to IRD
	Ifrom <- cbind(I2W, I2R, I2D) # to WRD
	Rfrom <- cbind(R2W, R2I, R2D) # to WID
	
	c(cbind(Wfrom, Ifrom, Rfrom))
}

piout2U <- function(pilongout){
	
	dim(pilongout) <- c(49,9)
	
	# extract parts, then get complement and reorder
	W2IRD <- pilongout[,1:3]
	I2WRD <- pilongout[,4:6]
	R2WID <- pilongout[,7:9]
	
	W2WIR <- cbind(1-rowSums(W2IRD), W2IRD[,1:2])
	I2WIR <- cbind(I2WRD[,1], 1-rowSums(I2WRD), I2WRD[,2])
	R2WIR <- cbind(R2WID[,1:2], 1-rowSums(R2WID))
	
	# make block mat to return
	cbind(rbind(pi2u(W2WIR[,1]),pi2u(W2WIR[,2]),pi2u(W2WIR[,3])),
		  rbind(pi2u(I2WIR[,1]),pi2u(I2WIR[,2]),pi2u(I2WIR[,3])),
		  rbind(pi2u(R2WIR[,1]),pi2u(R2WIR[,2]),pi2u(R2WIR[,3])))
	
}
# identical:
#max(abs(piout2U(extractfrom(A1)) -A[1:150,1:150]))

U2WLE50 <- function(U,comp = c(.7,.25,.05)){
	I   <- diag(nrow(U))
	Nsx <- solve(I - U)
	sum(colSums(Nsx[1:50, c(1,51,101)]) * comp)
}

piout2WLE50 <- function(pilongout,comp= c(.7,.25,.05)){
	U <- piout2U(pilongout)
	U2WLE50(U, comp)
}


decompDudel2 <- function(A1,A2,N=10,comp= c(.7,.25,.05)){
	
	pi1          <- extractfrom(A1)
	pi2          <- extractfrom(A2)
	
	# arrow decomposition:
	deltaxs      <- DecompHoriuchi::DecompContinuousOrig(
			          func = piout2WLE50, 
			          rates1 = pi1, 
			          rates2 = pi2, 
			          N = N, 
			          comp = comp)
	dim(deltaxs) <- c(length(deltaxs) / 9,9)
	deltaxs      <- colSums(deltaxs)
# from in columns, to in rows, in this case most arrows similar,
	wirnot  <- matrix(deltaxs, 3, 3, byrow = TRUE, dimnames = list(c("W", "I", "R"), c("W", "I", "R")))
	
	todeath <- diag(wirnot)
	diag(wirnot) <- NA
	rbind(wirnot, D = todeath)
}

decompbars  <- function(decmat,add = FALSE,cols=c("red", "blue")){
	
	x <- col(decmat) - 1
	y <- abs(row(decmat) - nrow(decmat))
	
	if(! add){
		plot(NULL, 
				xlim = c(0, ncol(decmat)),
				ylim =c (0, nrow(decmat)),
				type = "n", 
				asp = 1, 
				xlab = "", 
				ylab = "", 
				axes = FALSE)
	}
	x       <- c(x)
	y       <- c(y)
	values  <- c(decmat)
	lim     <- max(abs(values),na.rm=TRUE)
	valuesb <- values / (2*lim)
	rect(0:2,3:1,1:3,4:2,border=NA,col = gray(.9))
	rect(x, y + .5, x + 1, y + .5 + valuesb,
			col = ifelse(sign(valuesb) == 1, cols[1], cols[2]))	
	
	rect(x,y,x+1,y+1,border = gray(.8),lwd = 2)
	
	text(1:ncol(decmat)-.5,nrow(decmat),c("W","I","R"),font=2,pos=3,cex=1.5,xpd=TRUE)
	text(0,1:nrow(decmat)-.5,c("D","R","I","W"),font=2,pos=2,cex=1.5,xpd=TRUE)
	
}

decompbars2 <- function(
		decmat,
		add = FALSE,
		cols = c("red", "blue"),
		lim  =  max(abs(decmat), na.rm = TRUE)){
	
	x <- col(decmat) - 1
	y <- abs(row(decmat) - nrow(decmat))
	
	if(! add){
		plot(NULL, 
				xlim = c(0, ncol(decmat)),
				ylim =c (0, nrow(decmat)),
				type = "n", 
				asp = 1, 
				xlab = "", 
				ylab = "", 
				axes = FALSE)
	}
	x       <- c(x)
	y       <- c(y)
	values  <- c(decmat)
	valuesb <- values / lim
	rect(0:2,3:1,1:3,4:2,border=NA,col = gray(.9))
	rect(x, y , x + 1, y + abs(valuesb),
			col = ifelse(sign(valuesb) == 1, cols[1], cols[2]))	
	
	rect(x,y,x+1,y+1,border = gray(.8),lwd = 2)
	
	text(1:ncol(decmat)-.5,nrow(decmat),c("W","I","R"),font=2,pos=3,cex=1.5,xpd=TRUE)
	text(0,1:nrow(decmat)-.5,c("D","R","I","W"),font=2,pos=2,cex=1.5,xpd=TRUE)
	
}

A1 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_2004.csv"))
A2 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_2009.csv"))
#devtools::install_github("timriffe/DecompHoriuchi/DecompHoriuchi")
library(DecompHoriuchi)

decompA <- decompDudel2(A1,A2)
sum(decompA, na.rm=TRUE)
# works:
piout2WLE50(extractfrom(A2)) - piout2WLE50(extractfrom(A1))

decompbars2(decompA)
# problem: how to visualize contrib with directionality?
B1 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_edu0_2009.csv"))
B2 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_edu2_2009.csv"))

decompB  <- decompDudel2(B1,B2)
decompbars2(decompB)
sum(decompB, na.rm = TRUE)
piout2WLE50(extractfrom(B2)) - piout2WLE50(extractfrom(B1))


library(xtable)
(thetaA1 <- piout2WLE50(extractfrom(A1)))
(thetaA2 <- piout2WLE50(extractfrom(A2)))

(thetaB1 <- piout2WLE50(extractfrom(B1)))
(thetaB2 <- piout2WLE50(extractfrom(B2)))

thetaA1 - thetaA2
xtable(decompA)
pdf("/home/tim/git/ArrowDecomp/ArrowDecomp/Figures/decA.pdf",width=5,height=6)
par(mai=c(.1,.5,.5,.1))
decompbars2(decompA)
dev.off()

lim <- max(abs(c(decompA,decompB)),na.rm=TRUE)


pdf("/home/tim/git/ArrowDecomp/ArrowDecomp/Figures/decA2.pdf", width = 5, height = 6)
par(mai=c(.1,.5,.5,.1))
decompbars2(decompA,lim=lim)
dev.off()

pdf("/home/tim/git/ArrowDecomp/ArrowDecomp/Figures/decB.pdf", width = 5, height = 6)
par(mai=c(.1,.5,.5,.1))
decompbars2(decompB, lim = lim)
dev.off()
thetaB2 - thetaB1
