# modification by Dmytro Kryvokhyzha:
#
# multiple file input, filter of biased SNPs on the fly
#
# command: 
#
# Rscript RNAmodel_exctract_results.R input1.gz input2.gz input3.gz 

cutoff <- 0.99 # significance level

# the posterior probability of ASE
posteriorProbASE <- function(p, e, f, g, a, d, p0) {
  postN <- dbeta(p,f,g)*dbeta(e,1,h)*(1-p0)
  postD <- dbeta(p,a,a)*dbeta(e,1,d)*p0
  postN / (postN + postD)
}

inv.logit <- function(val) exp(val)/(1 + exp(val))

ff <- commandArgs(trailingOnly = T) # a list of files to process
source('readGzippedMcmcOutput.R')

# loop through all the files
for (file in ff) {
  # load the results
  result <- read.mcmc(file)
  message(paste(file, 'loaded\n', sep=' '))
  # character vector of gene names
  a.hat <- result$a.hat
  d.hat <- result$d.hat
  # below we take the posterior median of MCMC samples that are vectors of length equal to n.iter
  pi0 <- median(inv.logit(result$mcmc$logitPi0))
  f <- median(exp(result$mcmc$logF))
  g <- median(exp(result$mcmc$logG))
  h <- median(exp(result$mcmc$logH))
  # below we take the posterior median across columns of matrices with nrow==n.iter and ncol==length(genes)
  p <- apply(inv.logit(result$logit_p), 2, median)
  e <- apply(inv.logit(result$logit_e), 2, median)
  # calculate posterior prob of ASE
  postProbAse <- rep(NA, length(p))
  for (i in c(1:length(p))) {
    postProbAse[i] <- posteriorProbASE(p[i], e[i], f, g, a.hat, d.hat, pi0)
  }
  # select only signiciant genes
  genes <- result$features
  calledSignificant <- which(postProbAse >= cutoff)
  SignASE <- cbind(genes[calledSignificant],postProbAse[calledSignificant])
  # output the results
  ASE <- cbind(genes,postProbAse)
  write.table(SignASE, paste(file, 'signASE.csv', sep='_'), quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(ASE, paste(file, 'allASE.csv', sep='_'), quote = F, sep = "\t", row.names = F, col.names = F)
  message(paste(file, 'processed\n', sep=' '))
}
