#!/usr/bin/env Rscript
options(warn=2)

# model for RNA data that allows for within-gene variability
# params <- c("LOGIT_PI0", "LOG_F", "LOG_G", "LOG_H", "LOGIT_P_I", "LOGIT_E_I")
suppressPackageStartupMessages(library(optparse))
optionList <- list(make_option(c("-n", "--n.iter"), type="integer", default=500, help="Number of iterations [default %default]"),
	make_option(c("-t", "--thin"), type="integer", help="save every THIN'th iteration of MCMC [default %default]", default=10),
	make_option(c("-s", "--n.scaling.iter"), type="integer", help="number of scaling iterations [default %default]", default=0),
	make_option(c("-a", "--a.hat"), type="double", help="estimate of a (from genomic DNA model)"),
	make_option(c("-d", "--d.hat"), type="double", help="estimate of d (from genomic DNA model)"),
	make_option(c("-m", "--max.rounds.of.scaling"), type="integer", help="maximum number of scaling rounds to undertake"),
	make_option("--seed", type="integer", help="seed for random number generator"))
parser <- OptionParser(usage="%prog [options] dataset1 [dataset2 dataset3 ...] outfile", option_list=optionList)
# a hack to fix a bug in optparse that won't let you use positional args if you also have non-boolean optional args:
getOptionStrings <- function(parserObj) {
	optionStrings <- character()
	for(item in parserObj@options) {
		optionStrings <- append(optionStrings, c(item@short_flag, item@long_flag))
	}
	optionStrings
}
optStrings <- getOptionStrings(parser)
arguments <- parse_args(parser, positional_arguments=TRUE)
args <- Filter(function(arg) !(strsplit(arg, "=")[[1]][1] %in% optStrings), arguments$args)
options <- Filter(function(arg) strsplit(arg, "=")[[1]][1] %in% optStrings, arguments$args)
opts <- arguments$options
for(opt in options) {
	optsplit <- strsplit(opt, "=")
	opts[[sub("--", "", optsplit[[1]][1])]] <- as.integer(optsplit[[1]][2])
}
scriptname <- ifelse(any(grepl("--file", commandArgs())), strsplit(commandArgs()[grep("--file", commandArgs())], '=', fixed=TRUE)[[1]][2], "stdin")
start.time <- proc.time()[3]

n.iter <- opts$n.iter
thin <- opts$thin
n.scaling.iter <- opts$n.scaling.iter
a.hat <- opts$a.hat
d.hat <- opts$d.hat

dataset <- list(y=list(), N=list(), geneIndices=list())
datasetArgs <- Filter(function(val) !val %in% options, arguments$args)
if(length(datasetArgs) < 2) {
	cat("Error: You must specify at least one dataset and outfile!\n")
	print_help(parser)
	quit(status=1)
}
readData <- function(filename) {
	dd <- read.table(filename, header=TRUE, stringsAsFactors=FALSE)
	##dd <- dd[dd[,3] + dd[,4] > 0,]     # don't consider SNPs with no data
	dd <- dd[order(dd[,1]),]
	features <- unique(dd[,1])
	y <- dd[,3]
	N <- dd[,3] + dd[,4]
	geneIndices <- tapply(N, dd[,1])    # keeps track of which gene each y/N value belongs to
	list(y=y, N=N, geneIndices=geneIndices, features=features)
}

if(is.null(opts$seed)) {
	seed <- round(runif(1, 1, 5000))
} else {
	seed <- opts$seed
}
set.seed(seed)
cat(sprintf("Running with seed %d for reproducibility. PID was %d\n", seed, Sys.getpid()))
n.datasets <- length(datasetArgs) - 1
outfile <- datasetArgs[n.datasets + 1]
cat(sprintf("writing output to %s\n", outfile))
cat("reading in data ...\n")
for(i in 1:n.datasets) {
	z <- readData(datasetArgs[i])
	if(! exists("features")) {
		features <- z$features
	} else {
		if(! all.equal(features, z$features)) {
			cat("Error: features in this file are not the same as those in the last file!\n")
			quit(status=1)
		}
	}
	dataset$y[[i]] <- z$y
	dataset$N[[i]] <- z$N
	if(! exists("geneIndices")) {
		geneIndices <- z$geneIndices
	} else {
		if(! all.equal(geneIndices, z$geneIndices)) {
			cat("Error: difference in geneIndices!\n")
			quit(status=1)
		}
	}
}
dataset$geneIndices <- geneIndices
cat("loaded data successfully\n", file=stderr())

logit <- function(val) log(val/(1-val))
inv.logit <- function(val) {
	result <- exp(val)/(1+exp(val))
	if(any(is.nan(result))) {
		stopifnot(all(! is.nan(result[val < 700])))
		result[val >= 700] <- 1
	}
	result
}
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
piei.dens <- function(p_i, e_i, f, g, h, pi0, log=TRUE) {
	log(pi0*dbeta(p_i, a.hat, a.hat)*dbeta(e_i, 1, d.hat) + (1-pi0)*dbeta(p_i, f, g)*dbeta(e_i, 1, h))
}
Y_i.density <- function(yList, nList, p, e) {	
	# robust calculation of the Y_i density given p, e
	# dbetabin starts having trouble once alpha or beta > approx 1e9
	stopifnot(length(p) == n.genes)
	alpha <- p*(1-e)/e
	beta <- (1-p)*(1-e)/e
	##stopifnot(all(alpha > 1e-9) & all(beta > 1e-9))   # alpha/beta should really never be super, super small

	alphaTooBig <- alpha > 1e9
	betaTooBig <- beta > 1e9
	alphaTooSmall <- alpha < 1e-300
	betaTooSmall <- beta < 1e-300
	genesTooBigOrSmall <- sum(alphaTooBig + betaTooBig + alphaTooSmall + betaTooSmall)
	newAlpha <- alpha
	newBeta <- beta
	if(genesTooBigOrSmall > 0) {
		tooBigOrSmall <- alphaTooBig | betaTooBig | alphaTooSmall | betaTooSmall
		newAlpha[tooBigOrSmall] <- NA
		newBeta[tooBigOrSmall] <- NA
	}
	betabinResult <- 0
	for(i in 1:n.datasets) {
		betabinResult <- betabinResult + dbetabin(yList[[i]], nList[[i]], newAlpha[geneIndices], newBeta[geneIndices])
	}
	byGeneResult <- tapply(betabinResult, geneIndices, sum)
	
	if(genesTooBigOrSmall > 0) {
		# numerical issues: some alphas or betas were very big or very small
		stopifnot(all(is.na(byGeneResult[tooBigOrSmall])))
		for(i in 1:n.datasets) {
			yList[[i]][! tooBigOrSmall[geneIndices]] <- NA
			nList[[i]][! tooBigOrSmall[geneIndices]] <- NA
		}
		p[! tooBigOrSmall] <- NA
		
		# can approximate using dbinom if alpha == beta
		alphaInfinite <- ! is.finite(alpha)
		betaInfinite <- ! is.finite(beta)
		alphaAndBetaInfinite <- alphaInfinite & betaInfinite
		alphaAndBetaEqual <- alphaAndBetaInfinite | (alphaTooBig & betaTooBig & abs(1 - alpha/beta) < 0.001)
		if(sum(alphaAndBetaEqual) > 0) {
			binomResult <- 0
			for(i in 1:n.datasets) {
				binomResult <- binomResult + dbinom(yList[[i]], nList[[i]], p[geneIndices], log=TRUE)
			}
			binomByGene <- tapply(binomResult, geneIndices, sum)
			byGeneResult[alphaAndBetaEqual] <- binomByGene[alphaAndBetaEqual]
			
			for(i in 1:n.datasets) {
				yList[[i]][alphaAndBetaEqual[geneIndices]] <- NA
				nList[[i]][alphaAndBetaEqual[geneIndices]] <- NA
			}
		}
		
		# if only one of alpha, beta is infinite, reject
		oneOfAlphaOrBetaInfinite <- xor(alphaInfinite, betaInfinite)
		if(sum(oneOfAlphaOrBetaInfinite) > 0) {
			byGeneResult[oneOfAlphaOrBetaInfinite] <- -Inf
			for(i in 1:n.datasets) {
				yList[[i]][oneOfAlphaOrBetaInfinite[geneIndices]] <- NA
				nList[[i]][oneOfAlphaOrBetaInfinite[geneIndices]] <- NA
			}
		}
		
		# genes with alpha or beta small
		alphaOrBetaSmall <- alphaTooSmall | betaTooSmall
		if(sum(alphaOrBetaSmall) > 0) {
		    byGeneResult[alphaOrBetaSmall] <- -Inf    # can't calculate betabinomial density for alpha/beta --> 0
		    for(i in 1:n.datasets) {
				yList[[i]][alphaOrBetaSmall[geneIndices]] <- NA
				nList[[i]][alphaOrBetaSmall[geneIndices]] <- NA
		    }
		}

		# genes with alpha != beta and at least one big:
		oneOfAlphaOrBetaBig <- (alphaTooBig | betaTooBig) & ! alphaAndBetaEqual & ! oneOfAlphaOrBetaInfinite
		if(sum(oneOfAlphaOrBetaBig) > 0) {
			betabinSlowResult <- 0
			for(i in 1:n.datasets) {
				betabinSlowResult <- betabinSlowResult + dbetabin.slow(yList[[i]], nList[[i]], alpha[geneIndices], beta[geneIndices])
			}
			betabinSlowByGene <- tapply(betabinSlowResult, geneIndices, sum)
			byGeneResult[oneOfAlphaOrBetaBig] <- betabinSlowByGene[oneOfAlphaOrBetaBig]
		}
	}
	if(any(byGeneResult > 0)) {
		# error somewhere
		cat("Error: one of the byGeneResult values was > 0! Check out scriptname.PID.debugging.Rout\n", file=stderr())
		save(yList, NList, p, e, file=sprintf("%s.%d.debugging.Rout", scriptname, Sys.getpid()))
		stop()
	}
	byGeneResult
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

update.logitPi0 <- function(logitPi0.old, logF, logG, logH, logit_p_i, logit_e_i, scale) {
	logitPi0.new <- rnorm(1, logitPi0.old[1], scale)
	pi0.old <- inv.logit(logitPi0.old[1])
	pi0.new <- inv.logit(logitPi0.new)
	f <- exp(logF[1])
	g <- exp(logG[1])
	h <- exp(logH[1])
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	loglik.old <- sum(piei.dens(p_i, e_i, f, g, h, pi0.old)) + log(pi0.old) + log(1 - pi0.old)
	loglik.new <- sum(piei.dens(p_i, e_i, f, g, h, pi0.new)) + log(pi0.new) + log(1 - pi0.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0.old=logitPi0.old[1], logitPi0.new=logitPi0.new, logF=logF[1], logG=logG[1], logH=logH[1], logit_p_i=logit_p_i[,1], logit_e_i=logit_e_i[,1]))
	if(accept) {
		logitPi0.old[1] <- logitPi0.new
		logitPi0.old[2] <- logitPi0.old[2] + 1    # increment the number of acceptances
	}
	logitPi0.old
}
update.logF <- function(logitPi0, logF.old, logG, logH, logit_p_i, logit_e_i, scale) {
	pi0 <- inv.logit(logitPi0[1])
	logF.new <- rnorm(1, logF.old[1], scale)
	f.old <- exp(logF.old[1])
	f.new <- exp(logF.new)
	g <- exp(logG[1])
	h <- exp(logH[1])
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	
	q.old <- f.old/(f.old + g)
	r.old <- 1/(1 + f.old + g)
	q.new <- f.new/(f.new + g)
	r.new <- 1/(1 + f.new + g)
	loglik.old <- sum(piei.dens(p_i, e_i, f.old, g, h, pi0)) + dbeta(q.old, 100, 100, log=TRUE) + dbeta(r.old, 1, 20, log=TRUE) + logF.old[1] + logG[1] - log(f.old + g) - 2*log(1 + f.old + g)
	loglik.new <- sum(piei.dens(p_i, e_i, f.new, g, h, pi0)) + dbeta(q.new, 100, 100, log=TRUE) + dbeta(r.new, 1, 20, log=TRUE) + logF.new + logG[1] - log(f.new + g) - 2*log(1 + f.new + g)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0=logitPi0[1], logF.old=logF.old[1], logF.new=logF.new, logG=logG[1], logH=logH[1], logit_p_i=logit_p_i[,1], logit_e_i=logit_e_i[,1]))
	if(accept) {
		logF.old[1] <- logF.new
		logF.old[2] <- logF.old[2] + 1		# increment the number of acceptances
	}
	logF.old
}
update.logG <- function(logitPi0, logF, logG.old, logH, logit_p_i, logit_e_i, scale) {
	pi0 <- inv.logit(logitPi0[1])
	f <- exp(logF[1])
	logG.new <- rnorm(1, logG.old[1], scale)
	g.old <- exp(logG.old[1])
	g.new <- exp(logG.new)
	h <- exp(logH[1])
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	
	q.old <- f/(f + g.old)
	r.old <- 1/(1 + f + g.old)
	q.new <- f/(f + g.new)
	r.new <- 1/(1 + f + g.new)
	loglik.old <- sum(piei.dens(p_i, e_i, f, g.old, h, pi0)) + dbeta(q.old, 100, 100, log=TRUE) + dbeta(r.old, 1, 20, log=TRUE) + logF[1] + logG.old[1] - log(f + g.old) - 2*log(1 + f + g.old)
	loglik.new <- sum(piei.dens(p_i, e_i, f, g.new, h, pi0)) + dbeta(q.new, 100, 100, log=TRUE) + dbeta(r.new, 1, 20, log=TRUE) + logF[1] + logG.new - log(f + g.new) - 2*log(1 + f + g.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0=logitPi0[1], logF=logF[1], logG.old=logG.old[1], logG.new=logG.new, logH=logH[1], logit_p_i=logit_p_i[,1], logit_e_i=logit_e_i[,1]))
	if(accept) {
		logG.old[1] <- logG.new
		logG.old[2] <- logG.old[2] + 1		# increment the number of acceptances
	}
	logG.old
}
update.logH <- function(logitPi0, logF, logG, logH.old, logit_p_i, logit_e_i, scale) {
	pi0 <- inv.logit(logitPi0[1])
	f <- exp(logF[1])
	g <- exp(logG[1])
	logH.new <- rnorm(1, logH.old, scale)
	h.old <- exp(logH.old[1])
	h.new <- exp(logH.new)
	p_i <- inv.logit(logit_p_i[,1])
	e_i <- inv.logit(logit_e_i[,1])
	loglik.old <- sum(piei.dens(p_i, e_i, f, g, h.old, pi0)) + dexp(h.old, 0.03, log=TRUE) + logH.old[1]
	loglik.new <- sum(piei.dens(p_i, e_i, f, g, h.new, pi0)) + dexp(h.new, 0.03, log=TRUE) + logH.new
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0=logitPi0[1], logF=logF[1], logG=logG[1], logH.old=logH.old[1], logH.new=logH.new, logit_p_i=logit_p_i, logit_e_i=logit_e_i[,1]))
	if(accept) {
		logH.old[1] <- logH.new
		logH.old[2] <- logH.old[2] + 1    # increment the number of acceptances
	}
	logH.old
}
update.logit_p_i <- function(y, N, logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scale) {
	pi0 <- inv.logit(logitPi0[1])
	f <- exp(logF[1])
	g <- exp(logG[1])
	h <- exp(logH[1])
	e_i <- inv.logit(logit_e_i[,1])
	logit_p_i.old <- logit_p_i[,1]
	logit_p_i.new <- rnorm(length(logit_p_i.old), logit_p_i.old, scale)
	p_i.old <- inv.logit(logit_p_i.old)
	p_i.new <- inv.logit(logit_p_i.new)
	loglik.old <- Y_i.density(y, N, p_i.old, e_i) + piei.dens(p_i.old, e_i, f, g, h, pi0) + log(p_i.old) + log(1 - p_i.old)
	loglik.new <- Y_i.density(y, N, p_i.new, e_i) + piei.dens(p_i.new, e_i, f, g, h, pi0) + log(p_i.new) + log(1 - p_i.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0=logitPi0[1], logF=logF[1], logG=logG[1], logH=logH[1], logit_p_i.old=logit_p_i.old, logit_p_i.new=logit_p_i.new, logit_e_i=logit_e_i[,1]))
	logit_p_i[accept, 1] <- logit_p_i.new[accept]
	logit_p_i[accept, 2] <- logit_p_i[accept, 2] + 1
	stopifnot(all(is.finite(logit_p_i[,2])))
	logit_p_i
}
update.logit_e_i <- function(y, N, logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scale) {
	pi0 <- inv.logit(logitPi0[1])
	f <- exp(logF[1])
	g <- exp(logG[1])
	h <- exp(logH[1])
	p_i <- inv.logit(logit_p_i[,1])
	logit_e_i.old <- logit_e_i[,1]
	logit_e_i.new <- rnorm(length(logit_e_i.old), logit_e_i.old, scale)
	e_i.old <- inv.logit(logit_e_i.old)
	e_i.new <- inv.logit(logit_e_i.new)	
	loglik.old <- Y_i.density(y, N, p_i, e_i.old) + piei.dens(p_i, e_i.old, f, g, h, pi0) + log(e_i.old) + log(1 - e_i.old)
	loglik.new <- Y_i.density(y, N, p_i, e_i.new) + piei.dens(p_i, e_i.new, f, g, h, pi0) + log(e_i.new) + log(1 - e_i.new)
	accept <- MH.algorithm(loglik.old, loglik.new, list(logitPi0=logitPi0[1], logF=logF[1], logG=logG[1], logH=logH[1], logit_p_i=logit_p_i[,1], logit_e_i.old=logit_e_i.old, logit_e_i.new=logit_e_i.new))
	logit_e_i[accept, 1] <- logit_e_i.new[accept]
	logit_e_i[accept, 2] <- logit_e_i[accept, 2] + 1
	stopifnot(all(is.finite(logit_e_i[,2])))
	logit_e_i
}

do.mcmc <- function(n.iter, thin, y, N, logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList, verbose=FALSE, outfile=NULL) {
    if(is.null(outfile)) {
    	write.output <- FALSE
    } else {
    	write.output <- TRUE
    }
    for(i in 1:n.iter) {
    	logitPi0 <- update.logitPi0(logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$pi0)
    	logF <- update.logF(logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$f)
    	logG <- update.logG(logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$g)
    	logH <- update.logH(logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$h)
    	logit_p_i <- update.logit_p_i(y, N, logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$p_i)
    	logit_e_i <- update.logit_e_i(y, N, logitPi0, logF, logG, logH, logit_p_i, logit_e_i, scaleList$e_i)
    	if(i %% thin == 0) {
	    	if(verbose && i %% floor(n.iter/5) == 0) cat(sprintf("On iteration %d\n", i), file=stderr())
	    	if(write.output) cat(logitPi0[1], logF[1], logG[1], logH[1], logit_p_i[,1], logit_e_i[,1], "\n", file=outfile)
		}
    }
    list(logitPi0=logitPi0, logF=logF, logG=logG, logH=logH, logit_p_i=logit_p_i, logit_e_i=logit_e_i)
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
# params <- c("LOGIT_PI0", "LOG_F", "LOG_G", "LOG_H", "LOGIT_P_I", "LOGIT_E_I")
# parameters below are matrices with column1 == value, column2 == number of acceptances
n.genes <- length(features); cat(sprintf("inferring for %d genes\n", n.genes))
LOGIT_PI0 <- c(logit(runif(1)), 0)
LOG_F <- c(log(runif(1, 1, 30)), 0)
LOG_G <- c(log(runif(1, 1, 30)), 0)
LOG_H <- c(log(runif(1, 1, 30)), 0)
LOGIT_P_I <- matrix(c(logit(rep(0.5, n.genes)), rep(0, n.genes)), ncol=2)
LOGIT_E_I <- matrix(c(logit(rep(0.01, n.genes)), rep(0, n.genes)), ncol=2)
cat("Initial values: logitPi0=", LOGIT_PI0[1], " LOG_F=", LOG_F[1], " LOG_G=", LOG_G[1], " LOG_H=", LOG_H[1], "\n", file=stderr())
start.time <- proc.time()[3] 

# scaling ...
datasetPrefix <- substr(datasetArgs[1], 1, 6)
if(file.exists(sprintf("%s.scaling.Rout", datasetPrefix))) {
	cat("loading scaling values from file\n", file=stderr())
	load(sprintf("%s.scaling.Rout", datasetPrefix))
} else if(file.exists("scaling.Rout")) {
	cat("loading scaling values from file\n", file=stderr())
	load("scaling.Rout")
} else {
	cat("randomly guessing at scaling values\n", file=stderr())
	scaleList <- list(pi0=3, f=0.7, g=0.7, h=0.7, p_i=rep(0.5, n.genes), e_i=rep(1.2, n.genes))
}
cat("Initial values: logitPi0=", LOGIT_PI0[1], " LOG_F=", LOG_F[1], " LOG_G=", LOG_G[1], " LOG_H=", LOG_H[1], "\n", file=stderr())
cat("Initial scales: pi0=", scaleList$pi0, " f=", scaleList$f, " g=", scaleList$g, " h=", scaleList$h, "\n", file=stderr())


acceptRates <- list(logitPi0=0, logF=0, logG=0, logH=0, logit_p_i=rep(0, n.genes), logit_e_i=rep(0, n.genes))
target.acceptance.rate <- 0.30
DONE <- ifelse(n.scaling.iter > 0, FALSE, TRUE)
k <- 1
while(! DONE) {
    cat("working on scaling the proposal distributions ...\n")
    result <- do.mcmc(n.scaling.iter, 1, dataset$y, dataset$N, LOGIT_PI0, LOG_F, LOG_G, LOG_H, LOGIT_P_I, LOGIT_E_I, scaleList, verbose=FALSE)
    acceptRates <- list(logitPi0=result$logitPi0[2]/n.scaling.iter, logF=result$logF[2]/n.scaling.iter, logG=result$logG[2]/n.scaling.iter, logH=result$logH[2]/n.scaling.iter, logit_p_i=result$logit_p_i[,2]/n.scaling.iter, logit_e_i=result$logit_e_i[,2]/n.scaling.iter)
    cat("acceptRates:\n")
    print(sapply(acceptRates, mySummary))
    logitPi0.scaling <- adjust.scale(scaleList$pi0, acceptRates$logitPi0)
    logF.scaling <- adjust.scale(scaleList$f, acceptRates$logF)
    logG.scaling <- adjust.scale(scaleList$g, acceptRates$logG)
    logH.scaling <- adjust.scale(scaleList$h, acceptRates$logH)
    logit_p_i.scaling <- adjust.scale(scaleList$p_i, acceptRates$logit_p_i)
    logit_e_i.scaling <- adjust.scale(scaleList$e_i, acceptRates$logit_e_i)

    scaleList <- list(pi0=logitPi0.scaling$scale, f=logF.scaling$scale, g=logG.scaling$scale, h=logH.scaling$scale, p_i=logit_p_i.scaling$scale, e_i=logit_e_i.scaling$scale)
    if(all(logitPi0.scaling$done, logF.scaling$done, logG.scaling$done, logH.scaling$done, logit_p_i.scaling$done, logit_e_i.scaling$done)) DONE <- TRUE
    # reset acceptance counts
    LOGIT_PI0 <- result$logitPi0; LOGIT_PI0[2] <- 0
    LOG_F <- result$logF; LOG_F[2] <- 0
    LOG_G <- result$logG; LOG_G[2] <- 0
    LOG_H <- result$logH; LOG_H[2] <- 0
    LOGIT_P_I <- result$logit_p_i; LOGIT_P_I[,2] <- 0
    LOGIT_E_I <- result$logit_e_i; LOGIT_E_I[,2] <- 0
    #save(scaleList, acceptRates, file=sprintf("%s.%s.scaling%d", scriptname, Sys.getpid(), k))
    k <- k+1
    if(! is.null(opts$max.rounds.of.scaling) & k > opts$max.rounds.of.scaling) {
		cat("finishing scaling: hit max.rounds.of.scaling specified ...\n", file=stderr())
    	break
    }
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
cat("#infiles=", datasetArgs[1:n.datasets], "\n", sep=" ", file=outfile.gz)
cat("#n.iter=", n.iter, "\n", sep="", file=outfile.gz)
cat("#thin=", thin, "\n", sep="", file=outfile.gz)
cat("#a.hat=", a.hat, "\n", sep="", file=outfile.gz)
cat("#d.hat=", d.hat, "\n", sep="", file=outfile.gz)
# can easily read dataset below back in with scan()
cat("#dataset nrow=", n.iter/thin, "\n", sep="", file=outfile.gz)
cat("#dataset ncol=", n.genes*2+4, "\n", sep="", file=outfile.gz)
cat("#n.datasets=", n.datasets, "\n", sep="", file=outfile.gz)
for(i in 1:n.datasets) {
	cat(sprintf("#y%d=", i), file=outfile.gz); cat(dataset$y[[i]], "\n", sep=" ", file=outfile.gz)
	cat(sprintf("#N%d=", i), file=outfile.gz); cat(dataset$N[[i]], "\n", sep=" ", file=outfile.gz)
}
cat("#geneIndices=", file=outfile.gz); cat(geneIndices, sep=" ", file=outfile.gz)
cat("\n#features=", file=outfile.gz); cat(features, sep=" ", file=outfile.gz)
cat("\n#names=logitPi0 logF logG logH", paste("logit_p_", 1:n.genes, sep=""), paste("logit_e_", 1:n.genes, sep=""), "\n", file=outfile.gz)

# do the real MCMC run
result <- do.mcmc(n.iter, thin, dataset$y, dataset$N, LOGIT_PI0, LOG_F, LOG_G, LOG_H, LOGIT_P_I, LOGIT_E_I, scaleList, verbose=TRUE, outfile=outfile.gz)

# print final acceptance rates (logA[2], logD[2], logit_p_i[,2], logit_e_i[,2])
# print final time to finish (in seconds)
cat("#acceptanceRates=", file=outfile.gz); cat(result$logitPi0[2]/n.iter, result$logF[2]/n.iter, result$logG[2]/n.iter, result$logH[2]/n.iter, result$logit_p_i[,2]/n.iter, result$logit_e_i[,2]/n.iter, sep=" ", file=outfile.gz)
cat("\n#scaling=", file=outfile.gz); cat(scaleList$pi0, scaleList$f, scaleList$g, scaleList$h, scaleList$p_i, scaleList$e_i, file=outfile.gz)
cat("\n#time=", proc.time()[3] - start.time, "\n", sep="", file=outfile.gz)
close(outfile.gz)
