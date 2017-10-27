# Author: tim
###############################################################################

##################################################
# This work moved to separate project 'Spells'
##################################################
# below code will not be developed here.
# ----------------------------------------
do.this <- FALSE
if(do.this){
library(markovchain)
library(spatstat)
# here I just read in the first matrix given in the supplementary material to:
# Christian Dudel & Mikko Myrskylä (2017)
# Working Life Expectancy at Age 50 in the United States and the Impact of the Great Recession
# Demography https://link.springer.com/article/10.1007/s13524-017-0619-6


# get a transition matrix
TM <- as.matrix(
		read.csv("/home/tim/git/TransientSymmetry/TransientSymmetry/Data/Pmat_b_f_1994.csv",
				check.names=FALSE)
)


# transpose for standard MC stuff.
# demographers do stuff backwards I guess
TM <- t(TM)

# make s4 transition matrix from markovchain package
mcEmpl <- new("markovchain", states = rownames(TM),
		byrow = TRUE, transitionMatrix = TM,
		name = "Empl")

# how many sequences should we generate?
N      <- 1e4

# each assuming a start in employment at age 50.
RTraj  <- replicate(N,
		rmarkovchain(n = 50, object = mcEmpl, t0 = "50::Employed", parallel = TRUE)
) 
RTraj                 <- rbind(rep("50::Employed",N), RTraj)
RTraj_clean           <- gsub(".*:","",RTraj)
rownames(RTraj_clean) <- 50:100
x <- RTraj_clean[,1]
rep(1:3,1:3)
spell_durAge <- function(x, state = "Inactive", not_first = FALSE){
	
	# starting left x position for each observation
	x_age             <- 1:length(x) - 1
	names(x_age)      <- x

    sec               <- rle(x)
	spells            <- sec$values
	durs              <- sec$lengths
	durs2             <- durs
	durs2[spells != state] <- NA
    if (not_first){
		if (state %in% x){
			durs2[which(spells == state)[1]] <- NA
		}
		
	}
	
    #n_spells          <- sum(spells == state)
	spell_age         <- rep(durs2, durs)
	#spell_age[x != state] <- NA
	
	spell_age
}
EmplSpells <- apply(RTraj_clean, 2, spell_durAge, state = "Employed")
EmplSpellsH <- apply(RTraj_clean, 2, spell_durAge, state = "Employed",TRUE)


avg_dur_given_empl              <- rowMeans(EmplSpells,na.rm=TRUE)
employed_given_empl50           <- rowMeans(!is.na(EmplSpells))
avg_dur_given_emplH             <- rowMeans(EmplSpellsH,na.rm=TRUE)
employedH                       <- rowMeans(!is.na(EmplSpellsH))



plot(50:100, avg_dur_given_empl, type = 'l')
plot(50:100, employed_given_empl50, type = 'l')

plot(50:100, avg_dur_given_emplH, type = 'l')
plot(50:100, employedH, type = 'l')

plot(50:100, evg_dur_given_empl * employed_given_empl50)
plot(50:100, employedH * avg_dur_given_emplH)



InactSpells  <- apply(RTraj_clean, 2, spell_durAge, state = "Inactive")
InactSpellsH <- apply(RTraj_clean, 2, spell_durAge, state = "Inactive",TRUE)

avg_dur_InactSpells              <- rowMeans(InactSpells,na.rm=TRUE)
Inact_given_empl50               <- rowMeans(!is.na(InactSpells))
avg_dur_given_InactH             <- rowMeans(InactSpellsH,na.rm=TRUE)
InactH                           <- rowMeans(!is.na(InactSpellsH))


plot(50:100, avg_dur_InactSpells, type = 'l')
plot(50:100, Inact_given_empl50, type = 'l')

plot(50:100, avg_dur_given_InactH, type = 'l')
plot(50:100, InactH, type = 'l')

qdens <- function(X,p=.5){
	apply(X,1,function(x,p){
				if (sum(!is.na(x))>2){
					dx <- density(x,na.rm=TRUE)
					return(quantile(density(x,na.rm=TRUE),prob=p))
				} else {
					return(NA)
				}
			},p=p)
}

intervalpoly <- function(x,u,l,...){
	una <- is.na(u)
	lna <- is.na(l)
	keep <- !(una | lna)
	polygon(c(x[keep],rev(x[keep])),
			c(u[keep],rev(l[keep])),...)
}

intervalfan <- function(X,x=50:100,probs = seq(.1,.9,by=.1),...){
	for (i in probs){
		alpha <- i / 2
		ua    <- 1- alpha
		la    <- alpha
		ql <- qdens(X,la)
		qu <- qdens(X,ua)
		intervalpoly(x,ql,qu,...)
	}
}

sdInact <- sqrt(rowMeans((InactSpells - avg_dur_InactSpells) ^ 2, na.rm=TRUE))
png("/home/tim/git/ArrowDecomp/ArrowDecomp/Figures/spelldurationagepattern.png")
plot(50:100, avg_dur_InactSpells, type = 'l',ylim=c(0,18),
		main = "Conditional average 'inactive' spell duration\ngiven spell membership (not spell entry or exit)",
		sub = "quantile contours in 2.5% intervals",
		xlab = "Age", ylab = "spell duration",col="red",
		xlim = c(50,80))
intervalfan(InactSpells,50:100,border=NA,col = "#00000010",probs=seq(.05,.95,by=.05))
lines(50:100,qdens(InactSpells,.5),col = "black",lwd=2)
lines(50:100,avg_dur_InactSpells,col = "red",lwd=2)
text(55,5,"median",cex=1.5,font=)
text(55,8,"mean",cex=1.5,col = "red",font=2)
dev.off()

x <- RTraj_clean[,3]

spell_dur_before <- function(x, state = "Inactive", not_first = FALSE){
	
	# starting left x position for each observation
	x_age             <- 1:length(x) - 1
	names(x_age)      <- x
	
	sec               <- rle(x)
	spells            <- sec$values
	durs              <- sec$lengths
	durs2             <- durs
	durs2[spells != state] <- NA
	
	n          <- length(durs)
	x_lefts    <- cumsum(c(0,durs[-n]))
	x_lefts[spells != state] <- NA          
	if (not_first){
		if (state %in% x){
			x_lefts[which(spells == state)[1]] <- NA
		}
	}
	
	#n_spells          <- sum(spells == state)
	spell_starts         <- rep(x_lefts, durs)
	time_spent           <- x_age - spell_starts
	#spell_age[x != state] <- NA
	names(time_spent)    <- x_age
	time_spent + .5
}
spell_dur_after <- function(x, state = "Inactive", not_first = FALSE){
	
	# starting left x position for each observation
	x_age             <- 1:length(x) - 1
	names(x_age)      <- x
	
	sec               <- rle(x)
	spells            <- sec$values
	durs              <- sec$lengths
	durs2             <- durs
	durs2[spells != state] <- NA
	
	n          <- length(durs)
	x_lefts    <- cumsum(c(0,durs[-n]))
	x_rights   <- x_lefts + durs
	x_rights[spells != state] <- NA          
	if (not_first){
		if (state %in% x){
			x_rights[which(spells == state)[1]] <- NA
		}
	}
	
	#n_spells          <- sum(spells == state)
	spell_ends           <- rep(x_rights, durs)
	time_left            <- spell_ends - x_age
	#spell_age[x != state] <- NA
	names(time_left)    <- x_age
	time_left - .5
}

InactLeft  <- apply(RTraj_clean, 2, spell_dur_after, state = "Inactive")
InactSpent <- apply(RTraj_clean, 2, spell_dur_before, state = "Inactive")

InactLeftMean             <- rowMeans(InactLeft,na.rm=TRUE)
InactSpentMean             <- rowMeans(InactSpent,na.rm=TRUE)

plot(50:100, InactLeftMean, type = 'l',ylim=c(0,18),
		main = "Conditional average 'inactive' time left in spell\ngiven spell membership (not spell entry or exit)",
		sub = "quantile contours in 2.5% intervals",
		xlab = "Age", ylab = "spell duration",col="red",
		xlim = c(50,80))
intervalfan(InactLeft,50:100,border=NA,col = "#00000010",probs=seq(.05,.95,by=.05))
lines(50:100,qdens(InactLeft,.5),col = "black",lwd=2)
lines(50:100,InactLeftMean,col = "red",lwd=2)


plot(50:100, InactSpentMean, type = 'l',ylim=c(0,18),
		main = "Conditional average 'inactive' time spent in spell\ngiven spell membership (not spell entry or exit)",
		sub = "quantile contours in 2.5% intervals",
		xlab = "Age", ylab = "spell duration",col="red",
		xlim = c(50,80))
intervalfan(InactSpent,50:100,border=NA,col = "#00000010",probs=seq(.05,.95,by=.05))
lines(50:100,qdens(InactSpent,.5),col = "black",lwd=2)
lines(50:100,InactSpentMean,col = "red",lwd=2)
}

