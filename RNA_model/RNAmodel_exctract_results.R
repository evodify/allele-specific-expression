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
  # plotestimates
  jpeg(paste(file, 'jpeg', sep='.'), width=900, height = 600)
  par(mar=c(5, 5, 3, 2), mfcol = c(2,3))
  par(cex=1)
  hist(inv.logit(result$mcmc$logitPi0), main='', xlab='pi0', ylab = 'frequency', breaks = 50)
  hist(exp(result$mcmc$logF), main='', xlab='f', ylab = 'frequency', breaks = 50)
  hist(exp(result$mcmc$logG), main='', xlab='g', ylab = '', breaks = 50)
  hist(exp(result$mcmc$logH), main='', xlab='h', ylab = '', breaks = 50)
  hist(p, main='', xlab='p', ylab = '', breaks = 50)
  hist(e, main='', xlab='e', ylab = '', breaks = 50)
  dev.off()
  # calculate posterior prob of ASE
  postProbAse <- rep(NA, length(p))
  for (i in c(1:length(p))) {
    postProbAse[i] <- posteriorProbASE(p[i], e[i], f, g, a.hat, d.hat, pi0)
  }
  # get the FDR
  postProbSig <- which(postProbAse >= cutoff)
  fdr <- mean((1 - postProbAse)[postProbSig])
  # select only FDR signiciant genes
  genes <- result$features
  FDRsign <- which((1-postProbAse) <= fdr)
  SignASE <- cbind(genes[FDRsign],(1-postProbAse)[FDRsign])
  SignASEprint <- rbind(SignASE,c('overallFRD', fdr))
  # output the results
  ASE <- cbind(genes,(1-postProbAse))
  write.table(SignASEprint, paste(file, 'signASE.csv', sep='_'), quote = F, sep = "\t", row.names = F, col.names = F)
  write.table(ASE, paste(file, 'allASE.csv', sep='_'), quote = F, sep = "\t", row.names = F, col.names = F)
  message(paste(file, 'processed\n', sep=' '))
}
