# Simulate some data that we can use to test whether the model
# is giving us reasonable results.
# Note that we are simulating data according to the exact 
# statistical model we implemented, so this should be an easy task.
simulate.data <- function(seed=12345, n.snps=10000) {
    set.seed(seed)
    pi0 <- 0.95
    alpha <- rlnorm(1, 4.3, 1.8)
    delta <- 0.1
    epsilon <- 0.2
    avg.coverage <- 50
    N <- rpois(n.snps, lambda=avg.coverage)
    # generate beta-binomially distributed values:
    isNull <- runif(n.snps) < pi0
    Y <- isNull*rbinom(n.snps, N, rbeta(n.snps, alpha, alpha)) + 
        (1 - isNull)*rbinom(n.snps, N, rbeta(n.snps, delta, epsilon))
    list(Y=Y, N=N, isNull=isNull, alpha=alpha, delta=delta, epsilon=epsilon,
        pi0=pi0, n.snps=n.snps)
}
