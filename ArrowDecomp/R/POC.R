
# Proof of concept: arrow decomposition.

setwd("/home/tim/workspace/ArrowDecomp")


# 1 2 3
# 4 5 6
# 7 8 9
extractpi <- function(A){
	A         <- A[1:150, 1:150]
	# extract positive vectors of px:
	rowBlocks <- ((row(A)-1) - (row(A)-1) %% 50) / 50 + 1
	colBlocks <- ((col(A)-1) - (col(A)-1) %% 50) / 50 + 1
	
	pivecs <- matrix(NA, ncol = 9, nrow = 49)
	for (j in 1:3){ # columns
		for (i in 1:3){ # rows
			blockji <- A[rowBlocks == i & colBlocks == j]
			dim(blockji) <- c(50,50)
			pivec <- blockji[row(blockji) == (col(blockji) + 1)]
			pivecs[,(i-1)*3+j] <- pivec
		}
	}
	pivecs
}

pi2u <- function(pivec){
	cbind(rbind(0,diag(pivec)),0)
}
pi2U <- function(pivecs){
	u1 <- pi2u(pivecs[, 1])
	u2 <- pi2u(pivecs[, 2])
	u3 <- pi2u(pivecs[, 3])
	u4 <- pi2u(pivecs[, 4])
	u5 <- pi2u(pivecs[, 5])
	u6 <- pi2u(pivecs[, 6])
	u7 <- pi2u(pivecs[, 7])
	u8 <- pi2u(pivecs[, 8])
	u9 <- pi2u(pivecs[, 9])	
	rbind(cbind(u1,u2,u3),
	cbind(u4,u5,u6),
	cbind(u7,u8,u9))
}

U2WLE50 <- function(U,comp = c(.7,.25,.05)){
	I   <- diag(nrow(U))
	Nsx <- solve(I - U)
	sum(colSums(Nsx[1:50, c(1,51,101)]) * comp)
}


pilong2U <- function(pilong){
	dim(pilong) <- c(49,9)
	pi2U(pilong)
}

pilong2WLE50 <- function(pilong,comp= c(.7,.25,.05)){
	U <- pilong2U(pilong)
	U2WLE50(U, comp)
}


# this is somehow ignoring mortality.

decompDudel <- function(A1,A2,N=10,comp= c(.7,.25,.05)){
	
	pi1          <- c(extractpi(A1))
	pi2          <- c(extractpi(A2))
	
    # arrow decomposition:
	deltaxs      <- DecompHoriuchi::DecompContinuousOrig(
			          func = pilong2WLE50, 
			          rates1 = pi1, 
			          rates2 = pi2, 
			          N = N, 
			          comp = comp)
	dim(deltaxs) <- c(length(deltaxs) / 9,9)
	deltaxs      <- colSums(deltaxs)
# from in columns, to in rows, in this case most arrows similar,
	matrix(deltaxs, 3, 3, byrow = TRUE, dimnames = list(c("W", "I", "R"), c("W", "I", "R")))
}


A1 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_2004.csv"))
A2 <- as.matrix(read.csv("Data/Dudel/Transition matrices/Pmat_b_f_2009.csv"))


#devtools::install_github("timriffe/DecompHoriuchi/DecompHoriuchi")
library(DecompHoriuchi)

this.decomp <- decompDudel(A1,A2) # takes ca 1 min on Tim's old laptop

sum(this.decomp) # check sums

pilong2WLE50(c(extractpi(A2))) - pilong2WLE50(c(extractpi(A1))) 

# note, this for now ignores mortality differentials.
# Typically I'd decompose wrt event-exposure rates and not
# transition probabilities, but the difference probably isn't that hard.
# however, in a continuous matrix framework (Q matrix), I think
# one would get a separate estimate for the mortality differential.


