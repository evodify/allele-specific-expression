#!/usr/bin/env Rscript
library(rjags)
library(runjags)
rjags::load.module("mix")       # for dbetabin in JAGS

### Run the model and examine results ###
mod <- runjags::run.jags("biased_snps.bug", data=datlist, n.chains=4,
    monitor=c("pi0", "shape1", "shape2"))

# print medians of posterior distributions:
posterior.medians <- apply(as.mcmc(mod), 2, median)
print(posterior.medians)

pi0.hat <- posterior.medians["pi0"]
alpha.hat <- posterior.medians["shape1[2]"]    # == shape2[2]
delta.hat <- posterior.medians["shape1[1]"]
epsilon.hat <- posterior.medians["shape2[1]"]
