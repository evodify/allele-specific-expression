# Run a test of the model using simulated data
source('../simulate_data.R')

dat <- simulate.data()   # can specify n.snps=100 for a very quick run
n.snps <- length(dat$Y)
# columns 1 and 2 are currently ignored.
df <- data.frame(gene=rep("gene", n.snps), SNP=1:n.snps, dat$Y, dat$N-dat$Y)
filename <- "test_SNP_counts.txt"
write.table(df, file=filename, col.names=F, row.names=F, sep="\t", quote=F)

n.iter <- 100000
thin <- 100
source('2componentNullModel.bySNP.R')
unlink(filename)

posterior.medians <- sapply(result$final, function(x) median(x$sims))
print(posterior.medians)

pi0.hat <- posterior.medians["pi0"]
alpha.hat <- posterior.medians["alpha"]
delta.hat <- posterior.medians["delta"]
epsilon.hat <- posterior.medians["epsilon"]
