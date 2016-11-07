# This is a simple test to see how accurately we can infer
# the shape parameter alpha in a beta-binomial model
# where the beta prior is symmetric.
library(rjags)
rjags::load.module("mix")       # for dbetabin in JAGS

simulate.data <- function(seed, n.snps=1000, avg.coverage=30) {
    set.seed(seed)
    alpha <- rlnorm(1, 4.3, 1.8)
    cat(sprintf("True value of alpha was %.3g\n", alpha))
    N <- rpois(n.snps, lambda=avg.coverage)
    Y <- rbinom(n.snps, N, rbeta(n.snps, alpha, alpha))
    list(nsnps=n.snps, Y=Y, N=N)
}

model_string <- "model {
    alpha ~ dlnorm(4.3, 0.31)
    
    for (j in 1:nsnps) {
        Y[j] ~ dbetabin(alpha, alpha, N[j])
    }
}"
run.model <- function(datlist, quiet=TRUE) {
    jags <- jags.model(textConnection(model_string), data=datlist, 
        n.chains=4, quiet=quiet)
    progress.bar <- "text"
    if (quiet) progress.bar <- "none"
    update(jags, 2000, progress.bar=progress.bar)
    samp <- coda.samples(jags, c("alpha"), n.iter=10000, thin=5)
    cat("Posterior CIs for alpha from each chain:\n")
    print(sapply(samp, quantile, c(0.025, 0.25, 0.75, 0.975)))
    cat("\n")
    invisible(samp)
}

datlist <- simulate.data(seed=17311)
run.model(datlist)

datlist <- simulate.data(seed=14)
run.model(datlist)

datlist <- simulate.data(seed=1140)
run.model(datlist)
