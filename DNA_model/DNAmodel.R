#!/usr/bin/env Rscript

seed <- round(runif(1, 1, 10000))
set.seed(seed)
cat(sprintf("Running with seed %d for reproducibility. PID was %d\n", seed, Sys.getpid()))

# model for genomic DNA data that allows for within-gene variability
# params <- c("LOG_A", "LOG_D", "LOGIT_P_I", "LOGIT_E_I")
args <- commandArgs()
scriptname <- ifelse(any(grepl("--file", args)), strsplit(args[grep("--file", args)], '=', fixed=TRUE)[[1]][2], "stdin")
trailingArgs <- commandArgs(trailingOnly=TRUE)
if(interactive()) {
	infile <- file.choose()
	outfile <- file.choose()
	n.iter <- 10000
	thin <- 10
	n.scaling.iter <- 0
} else {
	stopifnot(length(trailingArgs) == 5)
	infile <- trailingArgs[1]
	outfile <- trailingArgs[2]
	n.iter <- as.integer(trailingArgs[3])
	thin <- as.integer(trailingArgs[4])
	n.scaling.iter <- as.integer(trailingArgs[5])
}
cat(sprintf("writing output to %s\n", outfile))

# read in data
cat("loading data\n", file=stderr())
dd <- read.table(infile, header=TRUE, stringsAsFactors=FALSE)
dd <- dd[dd[,3] + dd[,4] > 0,]     # don't consider SNPs with no data
dd <- dd[order(dd[,1]),]
features <- unique(dd[,1])
y <- dd[,3]
N <- dd[,3] + dd[,4]
geneIndices <- tapply(N, dd[,1])    # keeps track of which gene each y/N value belongs to
cat("loaded data successfully\n", file=stderr())

logit <- function(val) log(val/(1-val))
inv.logit <- function(val) exp(val)/(1+exp(val))
dbetabin <- function(y, n, alpha, beta, log=TRUE) {
	# works with matrix y, N and vector alpha, beta (gene-specific)
	toReturn <- lchoose(n, y) + lgamma(alpha+beta) - lgamma(n+alpha+beta) + lgamma(y+alpha) - lgamma(alpha) + lgamma(n-y+beta) - lgamma(beta)
	if(isTRUE(any(toReturn > 0))) {
		cat("Error: one of toReturn in dbetabin > 0: range(toReturn) was\n", file=stderr())
        	print(range(toReturn), file=stderr())
        	stop()
	}
	toReturn
}
dbetabin.slow <- function(y, n, alpha, beta, log=TRUE) {
 	# A (much) slower version of dbetabin that should be more numerically stable
	# works for vector y, n, alpha, beta (even for very, very large alpha, beta)
	stopifnot(length(y) == length(alpha))
	result <- lchoose(n, y)
	for(i in 1:length(y)) {
		if(is.na(result[i])) next
		if(n[i] == 0) next      # everything in dbetabin cancels out to be nothing
		nRange <- 0:(n[i]-1)
		yRange <- 0:(y[i]-1)
		nyRange <- 0:(n[i]-y[i]-1)
		if(y[i] == 0) {				# don't use yRange
			result[i] <- result[i] - sum(log(alpha[i]+beta[i]+nRange)) + sum(log(beta[i] + nyRange))
		} else if(n[i] - y[i] == 0) {    	# don't use nyRange
			result[i] <- result[i] - sum(log(alpha[i]+beta[i]+nRange)) + sum(log(alpha[i] + yRange))
		} else {				# use all ranges
			result[i] <- result[i] - sum(log(alpha[i]+beta[i]+nRange)) + sum(log(alpha[i] + yRange)) + sum(log(beta[i] + nyRange))
		}
	}
	result
}
# Metropolis-Hastings algorithm
# takes log unnormalized posterior densities as input (can be a vector)
MH.algorithm <- function(currentPosteriorVec, proposalPosteriorVec, params) {
	ratioVec <- proposalPosteriorVec - currentPosteriorVec
	if(any(proposalPosteriorVec == Inf)) {
		save(params, file=sprintf("%s.debug.infPosterior.Rsaved", scriptname))
	    stop("Error in Metropolis-Hastings algorithm: at least one of your proposal values gave an infinite posterior value!")
	}
	if(any(is.na(ratioVec))) {
		save(params, file=sprintf("%s.debug.NAinMH.Rsaved", scriptname))
		stop("Error in Metropolis-Hastings algorithm: one of proposalPostVec - currentPostVec is NA!")
	}
	randDraws <- runif(n=length(ratioVec))
	return(exp(ratioVec) >= randDraws)
}

Y_i.density <- function(y, N, p, e) {
	# robust calculation of the Y_i density given p, e
	# dbetabin starts having trouble once alpha or beta > approx 1e9
	stopifnot(length(p) == n.genes)
	alpha <- p*(1-e)/e
	beta <- (1-p)*(1-e)/e
	## next line is not necessary
	##stopifnot(all(alpha > 1e-9) & all(beta > 1e-9))   # alpha/beta should really never be super, super small

	alphaTooBig <- alpha > 1e9
	betaTooBig <- beta > 1e9
	genesTooBig <- sum(alphaTooBig + betaTooBig) 
	newAlpha <- alpha
	newBeta <- beta
	if(genesTooBig > 0) {
		newAlpha[alphaTooBig | betaTooBig] <- NA
		newBeta[alphaTooBig | betaTooBig] <- NA
	}
	betabinResult <- dbetabin(y, N, newAlpha[geneIndices], newBeta[geneIndices])
	byGeneResult <- tapply(betabinResult, geneIndices, sum)
	if(genesTooBig > 0) {
		# numerical issues: some alphas or betas were very big
		stopifnot(all(is.na(byGeneResult[alphaTooBig | betaTooBig])))
		y[! alphaTooBig[geneIndices] & ! betaTooBig[geneIndices]] <- NA
		N[! alphaTooBig[geneIndices] & ! betaTooBig[geneIndices]] <- NA
		p[! alphaTooBig & ! betaTooBig] <- NA
		
		# can approximate using dbinom if alpha == beta
		alphaAndBetaInfinite <- ! is.finite(alpha) & ! is.finite(beta)
		alphaAndBetaEqual <- alphaAndBetaInfinite | (alphaTooBig & betaTooBig & abs(1 - alpha/beta) < 0.001)
		if(sum(alphaAndBetaEqual) > 0) {
			binomResult <- dbinom(y, N, p[geneIndices], log=TRUE)
			binomByGene <- tapply(binomResult, geneIndices, sum)
			byGeneResult[alphaAndBetaEqual] <- binomByGene[alphaAndBetaEqual]
			y[alphaAndBetaEqual[geneIndices]] <- NA
			N[alphaAndBetaEqual[geneIndices]] <- NA
		}
		
		# genes with alpha != beta and at least one big:
		oneOfAlphaOrBetaBig <- (alphaTooBig | betaTooBig) & ! alphaAndBetaEqual
		if(sum(oneOfAlphaOrBetaBig) > 0) {
			betabinSlowResult <- dbetabin.slow(y, N, alpha[geneIndices], beta[geneIndices])
			betabinSlowByGene <- tapply(betabinSlowResult, geneIndices, sum)
			byGeneResult[oneOfAlphaOrBetaBig] <- betabinSlowByGene[oneOfAlphaOrBetaBig]
		}
	}
	if(any(byGeneResult > 0)) {
		# error somewhere
		save(y, N, p, e, file=sprintf("%s.%d.debugging.Rout", scriptname, Sys.getpid()))
		stop()
	}
	byGeneResult
}

update.logA <- function(logA.old, logD, logit_p_i, logit_e_i, scale) {
	d <- exp(logD[1])
	logA.new <- rnorm(1, logA.old[1], scale)
	a.old <- exp(logA.old[1])
	a.new <- exp(logA.new)
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	loglik.old <- sum(dbeta(p_i, a.old, a.old, log=TRUE) + dbeta(e_i, 1, d, log=TRUE)) + dnorm(logA.old[1], mean=4.3, sd=1.8, log=TRUE)
	loglik.new <- sum(dbeta(p_i, a.new, a.new, log=TRUE) + dbeta(e_i, 1, d, log=TRUE)) + dnorm(logA.new, mean=4.3, sd=1.8, log=TRUE)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logA.old=logA.old, logA.new=logA.new, logD=logD, logit_p_i=logit_p_i[,1], logit_e_i=logit_e_i[,1]))
	if(accept) {
		logA.old[1] <- logA.new
		logA.old[2] <- logA.old[2] + 1		# increment the number of acceptances
	}
	logA.old
}
update.logD <- function(logA, logD.old, logit_p_i, logit_e_i, scale) {
	a <- exp(logA[1])
	logD.new <- rnorm(1, logD.old[1], scale)
	d.old <- exp(logD.old[1])
	d.new <- exp(logD.new)
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	loglik.old <- sum(dbeta(p_i, a, a, log=TRUE) + dbeta(e_i, 1, d.old, log=TRUE)) + dexp(d.old, 0.0001, log=TRUE) + logD.old[1]
	loglik.new <- sum(dbeta(p_i, a, a, log=TRUE) + dbeta(e_i, 1, d.new, log=TRUE)) + dexp(d.new, 0.0001, log=TRUE) + logD.new
	accept <- MH.algorithm(loglik.old, loglik.new, list(logA=logA, logD.old=logD.old, logD.new=logD.new, logit_p_i=logit_p_i[,1], logit_e_i=logit_e_i[,1]))
	if(accept) {
		logD.old[1] <- logD.new
		logD.old[2] <- logD.old[2] + 1		# increment the number of acceptances
	}
	logD.old
}
update.logit_p_i <- function(y, N, logA, logD, logit_p_i, logit_e_i, scale) {
	a <- exp(logA[1])
	d <- exp(logD[1])
	e_i <- inv.logit(logit_e_i[,1])
	logit_p_i.old <- logit_p_i[,1]
	logit_p_i.new <- rnorm(length(logit_p_i.old), logit_p_i.old, scale)
	p_i.old <- inv.logit(logit_p_i.old)
	p_i.new <- inv.logit(logit_p_i.new)
	loglik.old <- Y_i.density(y, N, p_i.old, e_i) + dbeta(p_i.old, a, a, log=TRUE) + log(p_i.old) + log(1 - p_i.old)
	loglik.new <- Y_i.density(y, N, p_i.new, e_i) + dbeta(p_i.new, a, a, log=TRUE) + log(p_i.new) + log(1 - p_i.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logA=logA, logD=logD, logit_p_i.old=logit_p_i.old, logit_p_i.new=logit_p_i.new, logit_e_i=logit_e_i[,1]))
	logit_p_i[accept, 1] <- logit_p_i.new[accept]
	logit_p_i[accept, 2] <- logit_p_i[accept, 2] + 1
	stopifnot(all(is.finite(logit_p_i[,2])))
	logit_p_i
}
update.logit_e_i <- function(y, N, logA, logD, logit_p_i, logit_e_i, scale) {
	a <- exp(logA[1])
	d <- exp(logD[1])
	p_i <- inv.logit(logit_p_i[,1])
	logit_e_i.old <- logit_e_i[,1]
	logit_e_i.new <- rnorm(length(logit_e_i.old), logit_e_i.old, scale)
	e_i.old <- inv.logit(logit_e_i.old)
	e_i.new <- inv.logit(logit_e_i.new)
	loglik.old <- Y_i.density(y, N, p_i, e_i.old) + dbeta(e_i.old, 1, d, log=TRUE) + log(e_i.old) + log(1 - e_i.old)
	loglik.new <- Y_i.density(y, N, p_i, e_i.new) + dbeta(e_i.new, 1, d, log=TRUE) + log(e_i.new) + log(1 - e_i.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logA=logA, logD=logD, logit_p_i=logit_p_i[,1], logit_e_i.old=logit_e_i.old, logit_e_i.new=logit_e_i.new))
	logit_e_i[accept, 1] <- logit_e_i.new[accept]
	logit_e_i[accept, 2] <- logit_e_i[accept, 2] + 1
	stopifnot(all(is.finite(logit_e_i[,2])))
	logit_e_i
}

do.mcmc <- function(n.iter, thin, y, N, logA, logD, logit_p_i, logit_e_i, scaleList, verbose=FALSE, outfile=NULL) {
    if(is.null(outfile)) write.output <- FALSE
    else                 write.output <- TRUE
    for(i in 1:n.iter) {
        logA <- update.logA(logA, logD, logit_p_i, logit_e_i, scaleList$logA)
	logD <- update.logD(logA, logD, logit_p_i, logit_e_i, scaleList$logD)
	logit_p_i <- update.logit_p_i(y, N, logA, logD, logit_p_i, logit_e_i, scaleList$logit_p_i)
	logit_e_i <- update.logit_e_i(y, N, logA, logD, logit_p_i, logit_e_i, scaleList$logit_e_i)
	if(i %% thin == 0) {
	    if(verbose && i %% floor(n.iter/5) == 0) cat(sprintf("On iteration %d\n", i), file=stderr())
	    if(write.output) cat(logA[1], logD[1], logit_p_i[,1], logit_e_i[,1], "\n", file=outfile)
	}
    }
    list(logA=logA, logD=logD, logit_p_i=logit_p_i, logit_e_i=logit_e_i)
}

adjust.scale <- function(scale, accept.rate, target.acceptance.rate=0.30) {
    stopifnot(length(scale) == length(accept.rate))
    if(any(accept.rate < target.acceptance.rate-0.05) || any(accept.rate > target.acceptance.rate+0.05)) {
		too.low <- which(accept.rate < target.acceptance.rate & accept.rate != 0)
		too.high <- which(accept.rate > target.acceptance.rate & accept.rate != 1)
		scale[too.low] <- scale[too.low]/(20*abs(accept.rate[too.low] - target.acceptance.rate)^1.5 + 1)
		scale[too.high] <- scale[too.high]*(20*abs(accept.rate[too.high] - target.acceptance.rate)^1.5 + 1)
		scale[accept.rate == 0] <- scale[accept.rate == 0]/2
		scale[accept.rate == 1] <- scale[accept.rate == 1]*2
		list(done=FALSE, scale=scale)
    } else {
		list(done=TRUE, scale=scale)
    }
}

print.scaling.message <- function(paramName, oldScale, acceptRate, newScale) {
	cat(sprintf("for %s, old scale was %f and acceptance rate was %f so new scale is now %f\n", paramName, oldScale, acceptRate, newScale))
}
mySummary <- function(vals) { 
	  result <- summary(vals)[c(1,4,6)]   # just min/max/mean
	  class(result) <- c("summaryDefault", "table")
	  result
}


# Inference: a, d, p_i, e_i
# parameters below are matrices with column1 == value, column2 == number of acceptances
n.genes <- length(features); cat(sprintf("inferring for %d genes\n", n.genes))
LOG_A <- c(rnorm(1, 4.3, 1.8), 0)
LOG_D <- c(log(rexp(1, 0.0001)), 0)
LOGIT_P_I <- matrix(c(logit(rbeta(n.genes, exp(LOG_A[1]), exp(LOG_A[1]))), rep(0, n.genes)), ncol=2)
LOGIT_E_I <- matrix(c(logit(rbeta(n.genes, 1, exp(LOG_D[1]))), rep(0, n.genes)), ncol=2)
start.time <- proc.time()[3] 

# scaling ...
datasetPrefix <- strsplit(infile, ".", fixed=TRUE)[[1]][1]
cat("working on scaling to ensure that parameter acceptance rates for MCMC are reasonable\n", file=stderr())
if(file.exists(sprintf("%s.scaling.Rout", datasetPrefix))) {
        cat("loading scaling values from file\n", file=stderr())
        load(sprintf("%s.scaling.Rout", datasetPrefix))
} else if(file.exists("scaling.Rout")) {
        cat("loading scaling values from file\n", file=stderr())
        load("scaling.Rout")
} else if(file.exists("scaling.out")) {
        cat("loading scaling values from file\n", file=stderr())
	load("scaling.out")
} else {
        cat("randomly guessing at scaling values\n", file=stderr())
	cat("In the future, you may wish to do a scaling run, and save the scaling parameters for future runs\n", file=stderr())
	cat("To do this, run this program with n.iter=0, saving the output as scaling.Rout or scaling.out\n", file=stderr())
	scaleList <- list(logA=0.5, logD=0.5, logit_p_i=rep(0.25, n.genes), logit_e_i=rep(0.25, n.genes))
}
###  for running this script on mini datasets:
if(length(scaleList$logit_p_i) > n.genes) scaleList$logit_p_i <- scaleList$logit_p_i[1:n.genes]
if(length(scaleList$logit_e_i) > n.genes) scaleList$logit_e_i <- scaleList$logit_e_i[1:n.genes]
###
acceptRates <- list(logA=0, logD=0, logit_p_i=rep(0, n.genes), logit_e_i=rep(0, n.genes))
target.acceptance.rate <- 0.30
DONE <- ifelse(n.scaling.iter > 0, FALSE, TRUE)
k <- 1
while(! DONE) {
    cat("working on scaling the proposal distributions ...\n")
    result <- do.mcmc(n.scaling.iter, 1, y, N, LOG_A, LOG_D, LOGIT_P_I, LOGIT_E_I, scaleList, verbose=TRUE)
    acceptRates <- list(logA=result$logA[2]/n.scaling.iter, logD=result$logD[2]/n.scaling.iter, logit_p_i=result$logit_p_i[,2]/n.scaling.iter, logit_e_i=result$logit_e_i[,2]/n.scaling.iter)
    cat("acceptRates:\n")
    print(sapply(acceptRates, mySummary))
    logA.scaling <- adjust.scale(scaleList$logA, acceptRates$logA)
    logD.scaling <- adjust.scale(scaleList$logD, acceptRates$logD)
    logit_p_i.scaling <- adjust.scale(scaleList$logit_p_i, acceptRates$logit_p_i)
    logit_e_i.scaling <- adjust.scale(scaleList$logit_e_i, acceptRates$logit_e_i)

    scaleList <- list(logA=logA.scaling$scale, logD=logD.scaling$scale, logit_p_i=logit_p_i.scaling$scale, logit_e_i=logit_e_i.scaling$scale)
    if(all(logA.scaling$done, logD.scaling$done, logit_p_i.scaling$done, logit_e_i.scaling$done)) DONE <- TRUE
    # reset acceptance counts
    LOG_A <- result$logA; LOG_A[2] <- 0
    LOG_D <- result$logD; LOG_D[2] <- 0
    LOGIT_P_I <- result$logit_p_i; LOGIT_P_I[,2] <- 0
    LOGIT_E_I <- result$logit_e_i; LOGIT_E_I[,2] <- 0
    save(scaleList, acceptRates, file=sprintf("%s.%s.scaling%d", scriptname, Sys.getpid(), k))
    k <- k+1
}
if(n.iter == 0) {
	# a scaling run only
	save(scaleList, acceptRates, file=outfile)
	quit(status=0)
}

# prepare the output file
outfile.gz <- gzfile(outfile, open="wb")
cat("#seed=", seed, "\n", sep="", file=outfile.gz)
cat("#scriptname=", scriptname, "\n", sep="", file=outfile.gz)
cat("#infile=", infile, "\n", sep="", file=outfile.gz)
cat("#n.iter=", n.iter, "\n", sep="", file=outfile.gz)
cat("#thin=", thin, "\n", sep="", file=outfile.gz)
# can easily read dataset below back in with scan()
cat("#dataset nrow=", n.iter/thin, "\n", sep="", file=outfile.gz)
cat("#dataset ncol=", n.genes*2+2, "\n", sep="", file=outfile.gz)
cat("#y=", file=outfile.gz); cat(y, sep=" ", file=outfile.gz)
cat("\n#N=", file=outfile.gz); cat(N, sep=" ", file=outfile.gz)
cat("\n#geneIndices=", file=outfile.gz); cat(geneIndices, sep=" ", file=outfile.gz)
cat("\n#features=", file=outfile.gz); cat(features, sep=" ", file=outfile.gz)
cat("\n#names=logA logD", paste("logit_p_", 1:n.genes, sep=""), paste("logit_e_", 1:n.genes, sep=""), "\n", file=outfile.gz)

# do the real MCMC run
result <- do.mcmc(n.iter, thin, y, N, LOG_A, LOG_D, LOGIT_P_I, LOGIT_E_I, scaleList, verbose=TRUE, outfile=outfile.gz)

# print final acceptance rates (logA[2], logD[2], logit_p_i[,2], logit_e_i[,2])
# print final time to finish (in seconds)
cat("#acceptanceRates=", file=outfile.gz); cat(result$logA[2]/n.iter, result$logD[2]/n.iter, result$logit_p_i[,2]/n.iter, result$logit_e_i[,2]/n.iter, sep=" ", file=outfile.gz)
cat("\n#scaling=", file=outfile.gz); cat(scaleList$logA, scaleList$logD, scaleList$logit_p_i, scaleList$logit_e_i, file=outfile.gz)
cat("\n#time=", proc.time()[3] - start.time, "\n", sep="", file=outfile.gz)
close(outfile.gz)
