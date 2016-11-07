# Given estimates of pi0, alpha, delta, and epsilon from the 
# biased SNPs model, calculate the posterior probability
# that each SNP is biased.
#
# The estimates should be named pi0.hat, alpha.hat, delta.hat, and
# epsilon.hat. In the paper we used posterior medians for this estimate.

dbetabin <- function(y, n, alpha, beta, log=TRUE) {
	lchoose(n, y) + lgamma(alpha + beta) + lgamma(y + alpha) + 
		lgamma(n - y + beta) - lgamma(alpha) - lgamma(beta) -
		lgamma(n + alpha + beta)
}

posteriorProbBiased <- function(y, n) {
	comp1 <- pi0.hat*exp(dbetabin(y, n, alpha.hat, alpha.hat))
	comp2 <- (1 - pi0.hat)*exp(dbetabin(y, n, delta.hat, epsilon.hat))
	comp2/(comp1 + comp2)
}

# Change Y and N to reflect real data here if necessary.
# This code is for simulated data:
Y <- dat$Y
N <- dat$N
probs <- posteriorProbBiased(y=Y, n=N)
biased <- probs > 0.5
biased[N == 0] <- TRUE      # we'll want to filter out any SNPs with no data

# plot posterior probability biased:
hist(probs, col='grey80', main='Posterior Probability of ASE', 
    xlab='probability', cex.axis=1.5, cex.main=1.5, cex.lab=1.5)
# Most SNPs should have P(biased) close to 0 or 1

# Now we should filter out any SNPs with P(biased) > 0.5
