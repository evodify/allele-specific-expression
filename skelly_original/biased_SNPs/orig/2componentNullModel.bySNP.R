#!/usr/bin/env Rscript

# model for null ASE data (i.e. allelic read counts obtained from genomic DNA)
# read in the data
suppressPackageStartupMessages(library(optparse))
optionList <- list(make_option(c("-f", "--filename"), type="character", help="file of count data"),
	make_option(c("-n", "--n.iter"), type="integer", default=50000, help="Number of iterations [default %default]"),
	make_option(c("-t", "--thin"), type="integer", help="save every THIN'th iteration of MCMC [default %default]", default=10))
parser <- OptionParser(option_list=optionList)
opts <- parse_args(parser)
start.time <- proc.time()[3]

outfile <- sprintf("2componentNullModel.R.%s.out", Sys.getpid())
if(interactive()) {   # do not parse args
	if (! exists("filename")) {
		cat("Choose a file of read counts to analyze.\n"); flush.console()
		filename <- file.choose()
	}
	if (! exists("n.iter")) {
		opts$n.iter <- as.integer(readline("Enter the number of iterations: "))
	}
	if (! exists("thin")) {
		opts$thin <- as.integer(readline("Enter a value for thin: "))
	}
} else {
	if(is.null(opts$filename)) {
		cat("Incorrect number of required args!\n")
		print_help(parser)
		quit(status=1)
	}
	filename <- opts$filename
	cat(sprintf("saving output to %s\n", outfile))
}
dd <- read.table(filename, header=TRUE, stringsAsFactors=FALSE)
dd$total <- dd[, 3] + dd[, 4]
dd <- dd[dd$total > 0 & dd[, 1] != "intergenic",]
y <- dd[, 4]
N <- dd$total
cat("done reading in data ...\n")

dbetabin <- function(y, n, alpha, beta, log=TRUE) {
	lchoose(n, y) + lgamma(alpha+beta) + lgamma(y+alpha) + 
		lgamma(n-y+beta) - lgamma(alpha) - lgamma(beta) - lgamma(n+alpha+beta)
}
dYi <- function(y, n, pi0, alpha, delta, epsilon, log=TRUE) {
	result.vec <- pi0*exp(dbetabin(y, n, alpha, alpha, log=TRUE)) + (1-pi0)*exp(dbetabin(y, n, delta, epsilon, log=TRUE))
	sum(log(result.vec))
}
logit <- function(val) log(val/(1-val))
inv.logit <- function(val) exp(val)/(1+exp(val))
convertParams <- function(paramList) {
	pi0 <- inv.logit(paramList$logitPi0$value)
	alpha <- exp(paramList$logAlpha$value)
	delta <- inv.logit(paramList$logitDelta$value)
	epsilon <- inv.logit(paramList$logitEpsilon$value)
	list(pi0=pi0, alpha=alpha, delta=delta, epsilon=epsilon)
}

update.pi0 <- function(logitPi0, paramList) {
	params <- convertParams(paramList)
	pi0 <- inv.logit(logitPi0)
	dYi(y, N, pi0, params$alpha, params$delta, params$epsilon) + log(pi0) + log(1 - pi0) # last part is Jacobian for change of variable
}
update.alpha <- function(logAlpha, paramList) {
	params <- convertParams(paramList)
	alpha <- exp(logAlpha)
	dYi(y, N, params$pi0, alpha, params$delta, params$epsilon) + dnorm(logAlpha, mean=4.3, sd=1.8, log=TRUE)
}
update.delta <- function(logitDelta, paramList) {
	params <- convertParams(paramList)
	delta <- inv.logit(logitDelta)
	dYi(y, N, params$pi0, params$alpha, delta, params$epsilon) + log(delta) + log(1 - delta)
}
update.epsilon <- function(logitEpsilon, paramList) {
	params <- convertParams(paramList)
	epsilon <- inv.logit(logitEpsilon)
	dYi(y, N, params$pi0, params$alpha, params$delta, epsilon) + log(epsilon) + log(1 - epsilon)
}

# Inference
source('generic_MCMC.R')
params <- list(
	logitPi0=setup(init=logit(0.5), updateFunc=update.pi0, scale=0.2),
	logAlpha=setup(init=log(10), updateFunc=update.alpha, scale=1),
	logitDelta=setup(init=logit(0.5), updateFunc=update.delta, scale=4.4),
	logitEpsilon=setup(init=logit(0.5), updateFunc=update.epsilon, scale=4.4)
)
to.infer <- c("logitPi0", "logAlpha", "logitDelta", "logitEpsilon")
for(name in to.infer) params[[name]]$infer <- TRUE

simsToFinal <- function(simsList) {
	toReturn <- list()
	toReturn$pi0 <- list(accept.rate=simsList$logitPi0$accept.rate, scale=simsList$logitPi0$scale, sims=inv.logit(simsList$logitPi0$sims))
	toReturn$alpha <- list(accept.rate=simsList$logAlpha$accept.rate, scale=simsList$logAlpha$scale, sims=exp(simsList$logAlpha$sims))
	toReturn$delta <- list(accept.rate=simsList$logitDelta$accept.rate, scale=simsList$logitDelta$scale, sims=inv.logit(simsList$logitDelta$sims))
	toReturn$epsilon <- list(accept.rate=simsList$logitEpsilon$accept.rate, scale=simsList$logitEpsilon$scale, sims=inv.logit(simsList$logitEpsilon$sims))
	attr(toReturn, "class") <- "MCMC"
	toReturn
}

# do the MCMC
n.iter <- opts$n.iter
thin <- opts$thin
mcmc <- do.mcmc.safe(n.iter=n.iter, paramList=params, n.scaling.iter=1000, 
	thin=thin, verbose=TRUE, max.rounds.of.scaling=20,
	plot.scaling.progress=TRUE)
if(inherits(mcmc, "try-error"))
	stop(mcmc)
result <- list(scriptname="2componentNullModel.bySNP.R", n.iter=n.iter,
	thin=thin, dataset=list(y=y, N=N), mcmc=mcmc)
result$command <- sprintf("%s --filename=%s --n.iter=%d --thin=%d", 		
	result$scriptname, filename, n.iter, thin)
final <- try(simsToFinal(mcmc), silent=TRUE)
if(inherits(final, "try-error")) {
	result$final <- NULL
} else {
	result$final <- final
}
result$time <- proc.time()[3] - start.time
cat(sprintf("finished in %.1f seconds == %.2f minutes == %.2f hours\n", result$time, result$time/60, result$time/3600))

if(! interactive()) {
	save(result, logit, inv.logit, print.MCMC, file=outfile)
}
