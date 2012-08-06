random.effects <- function(object)
{
	if(!inherits(object, "Bayesthresh"))
		stop("Use an object of class Bayesthresh")
	tc <- (dim(object$compVar)[1])-1
	if (tc == 1){
		aleat <- data.frame(object$EfRandom, row.names=object$NomesZ)
		class(aleat) <- c("random.effects","Bayesthresh", "data.frame")
		return(aleat)
	}
	else {
		aleat <- NULL
		il <-  1
		ifc <- object$fl
		ic <- object$fl[1]
		for(i in 2:tc){
			ic[i] <- ifc[i]+ic[i-1]
			il[i] <- ic[i]-ifc[i]+1
		}
		for(i in 1:tc){
			aleat[[i]] <- data.frame(object$EfRandom[il[i]:ic[i],], row.names=object$NomesZ[il[i]:ic[i]])
		}
		names(aleat) <- rownames(object$compVar[1:(nrow(object$compVar)-1),])
		class(aleat) <- c("random.effects","Bayesthresh","list")
		return(aleat)
	}
}


