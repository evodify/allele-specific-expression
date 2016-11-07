##
## Routines for doing MCMC
##
##
opaque <- function(color, opacity, max=255) {
       # return a color that is the same as the color
       # argument of the function, but with the
       # specified opacity
       return(rgb(t(col2rgb(color)), alpha=opacity, max=max))
}


do.mcmc.safe <- function(n.iter, paramList, n.scaling.iter=4000, thin=1,
    verbose=TRUE, max.rounds.of.scaling=NULL, plot.scaling.progress=FALSE) {
	# first ensure that scaling is satisfactory
	if(n.scaling.iter > 0) {
		if(verbose) cat("starting scaling run ...\n")
		scalingList <- try(adapt.scale(n.scaling.iter, paramList, thin, 
            verbose, max.rounds.of.scaling, 
            plot.progress=plot.scaling.progress), silent=TRUE)
		if(inherits(scalingList, "try-error"))
			return(scalingList)
		else
			paramList <- scalingList
	}
	
	# now that scaling is good, do the MCMC
	if(n.iter > 0) {
		safeResult <- try(do.mcmc(n.iter, paramList, thin, verbose), silent=TRUE)
		if(inherits(safeResult, "try-error")) {
			return(safeResult)
		} else {
			attr(safeResult, "class") <- "MCMC"
			return(safeResult)
		}
	} else {
		return(paramList)
	}
}

do.mcmc <- function(n.iter, paramList, thin, verbose) {
	# initialize output data structure
	MCMCresult <<- vector("list")
	for(name in names(paramList))
		if(paramList[[name]]$infer) {		# only save results for parameters we care about for this run
			n <- length(paramList[[name]]$value)
			MCMCresult[[name]] <<- list(sims=matrix(nrow=floor(n.iter/thin), ncol=n), accept.rate=numeric(n), scale=paramList[[name]]$scale)
			dimnames(MCMCresult[[name]]$sims) <<- list(NULL, ifelse(rep(n, n) == 1, name, paste(name, '_', 1:n, sep='')))
		}
	
	iterGroups <- floor(n.iter/5)
	paramGroups <- names(paramList)
	for(i in 1:n.iter) {
		if(verbose && i %% iterGroups == 0) cat(sprintf("On iteration %d\n", i))
		for(j in 1:length(paramList)) {
			if(paramList[[j]]$infer) {
				if(paramList[[j]]$gibbs) {
					gibbsCall <<- paramList[[j]]$updateFunc(paramList)
					acceptVec <- gibbsCall$accept    # will typically be vector of 1's, unless gibbs has been hijacked to infer parameters that are not independent (e.g. p_i, e_i, component in ASE model)
					paramList[[j]]$value <- gibbsCall$value
				} else {
					# get current and proposal log-unnormalized posterior density
					currentVals <<- paramList[[j]]$value
					proposalVals <<- rnorm(length(currentVals), mean=currentVals, sd=paramList[[j]]$scale)
					currentLUD <<- paramList[[j]]$updateFunc(currentVals, paramList)
					proposalLUD <<- paramList[[j]]$updateFunc(proposalVals, paramList)
					if(any(proposalLUD == Inf)) {
						# changed from is.finite(): if the proposal is -Inf, we will always reject!
						stop(sprintf("While working on %s, at least one of your proposal values gave an infinite posterior value! Your proposalVals were %s", paramGroups[j], toString(proposalVals)))
					}
					
					# choose which proposals to accept and which to reject
					acceptVec <- MH.algorithm(currentLUD, proposalLUD)
					paramList[[j]]$value[acceptVec] <- proposalVals[acceptVec]
				}
				
				# save results
				MCMCresult[[paramGroups[j]]]$accept.rate <<- MCMCresult[[paramGroups[j]]]$accept.rate + acceptVec/n.iter
				if(i %% thin == 0)
					MCMCresult[[paramGroups[j]]]$sims[i/thin,] <<- paramList[[j]]$value
			}
		}
	}
	return(MCMCresult)
}

# adaptively scale proposal jumps to try to achieve close to target.acceptance.rate
adapt.scale <- function(n.scaling.iter, paramList, thin, verbose,
    max.rounds.of.scaling=NULL, target.acceptance.rate=0.3, expansion=20,
    plot.progress=FALSE) {
	if(interactive() && plot.progress) { 
        dev.new(width=8, height=4); par(mfrow=c(1, 2)); 
    }
	scaleResult <<- do.mcmc.safe(n.iter=n.scaling.iter, paramList=paramList, n.scaling.iter=0, thin, verbose=FALSE)
	if(inherits(scaleResult, "try-error")) {
		return(scaleResult)
	}
	
	DONE <- FALSE
	scaling.round <- 1
	while(! DONE) {
		DONE <- TRUE
		accept.allparams <<- numeric()
		for(name in names(scaleResult)) {
			if(!paramList[[name]]$infer || paramList[[name]]$gibbs) next
			accept.rate <- scaleResult[[name]]$accept.rate
			accept.allparams <<- c(accept.allparams, accept.rate)
			if(any(accept.rate < 0.25) || any(accept.rate > 0.35))
				DONE <- FALSE
			too.low <- which(accept.rate < target.acceptance.rate)
			too.high <- which(accept.rate > target.acceptance.rate)
			paramList[[name]]$scale[too.low] <- paramList[[name]]$scale[too.low]/(expansion*abs(accept.rate[too.low] - target.acceptance.rate)^1.5 + 1)
			paramList[[name]]$scale[too.high] <- paramList[[name]]$scale[too.high]*(expansion*abs(accept.rate[too.high] - target.acceptance.rate)^1.5 + 1)
			paramList[[name]]$value <- as.numeric(tail(scaleResult[[name]]$sims, 1))
		}
		if(interactive() && plot.progress) {
			hist(accept.allparams, col='lightblue', xlim=c(0, 1), ylim=c(0, length(accept.allparams)), xlab='accept rates', main='Histogram of accept rates')
			plot(accept.allparams, pch=19, cex=0.8, col='lightblue', ylim=c(0, 1), xlab='parameter index', ylab='accept rate')
			abline(h=target.acceptance.rate, col=opaque('red', 120), lty=2, lwd=2)
		}
		if(verbose) 
			cat(sprintf('after last %d iterations, accept rate was: mean=%.3f, range=(%.2f, %.2f)\n', n.scaling.iter, mean(accept.allparams), range(accept.allparams)[1], range(accept.allparams)[2]))
	
		if(! is.null(max.rounds.of.scaling) && scaling.round == max.rounds.of.scaling) break
		if(! DONE) {
			newParamList <<- paramList
			scaleResult <<- do.mcmc.safe(n.iter=n.scaling.iter, paramList=paramList, n.scaling.iter=0, thin, verbose=FALSE)
			scaling.round <- scaling.round + 1
			if(inherits(scaleResult, "try-error"))
				return(scaleResult)
		}
	}
	if(verbose)
		cat("Scaling is acceptable (Or hit max rounds of scaling). Continuing on to inference ...\n")
	return(paramList)
}

# Metropolis-Hastings algorithm
# takes log unnormalized posterior densities as input
MH.algorithm <- function(currentPostVec, proposalPostVec) {
	ratioVec <- proposalPostVec - currentPostVec
	if(any(is.na(ratioVec)))
		stop("Error in Metropolis-Hastings algorithm: one of proposalPostVec - currentPostVec is NA!")
	randDraws <- runif(n=length(ratioVec))
	return(exp(ratioVec) >= randDraws)
}

# set up a list corresponding to a single "parameter group"
setup <- function(initVec, updateFunc, scaleVec=NULL, gibbs=FALSE, infer=FALSE) {
	if(is.null(scaleVec)) {
		scaleVec <- rep(1, length(initVec))
	}
	myList <- list(value=initVec, scale=scaleVec, updateFunc=updateFunc, gibbs=gibbs, infer=infer)
	return(myList)
}

print.MCMC <- function(x, ...) {
	cat(sprintf("List of MCMC results with %d elements:\n", length(x)))
	cat(sprintf("%24s  %7s   %9s   %9s\n", "", "MEAN", "MEDIAN", "MEAN"))
	cat(sprintf("%-17s  %4s   %-7s   %9s   %9s\n", "PARAMETER", "N", "%ACCEPT", "VALUE", "SCALE"))
	options(warn=-17)   # turn off warnings
	for(name in names(x)) {
		if(is.list(x[[name]])) {
			acc <- x[[name]]$accept.rate
			meanVec <- c(mean(acc, na.rm=TRUE), median(x[[name]]$sims), mean(x[[name]]$scale))
			meanVec[abs(meanVec) < 1e-20] <- 0
			cat(sprintf("%-17s  %4d   %7.3f   ", name, ncol(x[[name]]$sims), meanVec[1]))
			cat(sprintf("%9.3g   %9.3g\n", meanVec[2], meanVec[3]))
		}
	}
	options(warn=0)
}
