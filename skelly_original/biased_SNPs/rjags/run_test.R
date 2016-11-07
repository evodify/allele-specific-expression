# Run a test of the model using simulated data
source('../simulate_data.R')

dat <- simulate.data()   # can specify n.snps=100 for a very quick run
n.snps <- length(dat$Y)

mat <- data.frame(Y=dat$Y, N=dat$N)
mat$p <- apply(mat, 1, function(row) binom.test(row[1], row[2])$p.value)

# set cluster1 for the SNP with least ASE and cluster2 the most, 
# according to the binomial exact test. This is as per instructions in
# http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of\
# -normal-distributions.html, for implementing a mixture model in JAGS
isNull <- rep(NA, n.snps)
isNull[which.min(mat$p)] <- FALSE
isNull[which.max(mat$p)] <- TRUE
datlist <- list(nsnps=n.snps, Y=dat$Y, N=dat$N)

source('biased_snps_rjags.R')
