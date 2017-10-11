
# Author: tim
###############################################################################

# redo POC to remove self-arrows and handle death transitions properly:

A <- A1
#  X   2   3
#  4   X   6
#  7   8   X
# 10  11  12


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
pilongout <- extractfrom(A1)
diag(A1[2:50,1:49])
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



decomptest <- decompDudel2(A1,A2)
sum(decomptest, na.rm=TRUE)
# works:
piout2WLE50(extractfrom(A2)) - piout2WLE50(extractfrom(A1))
