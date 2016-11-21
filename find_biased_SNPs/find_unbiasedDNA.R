# modification of the original file calculate_posterior_prob_biased.R
#
# original instactions:
# 
# Given estimates of pi0, alpha, delta, and epsilon from the 
# biased SNPs model, calculate the posterior probability
# that each SNP is biased.
# The estimates should be named pi0.hat, alpha.hat, delta.hat, and
# epsilon.hat. In the paper we used posterior medians for this estimate.
#
# modification by Dmytro Kryvokhyzha:
#
# multiple file input, filter of biased SNPs on the fly
#
# command: 
#
# Rscript unbiasedSNPs_sampleList.R input1.csv input2.csv input3.csv 


ff <- commandArgs(trailingOnly = T) # a list of files to process

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

# loop through all the files
for (file in ff) {

  dnaD <- na.omit(read.table(file, header = T, na.strings = 'NA'))
  message('\n')
  message(paste(file, 'loaded\n', sep=' '))
  dnaD$total <- apply(dnaD[, c(3,4)], 1, sum)
  
  n.snps <- length(dnaD[,1])
  mat <- data.frame(Y=dnaD[,3], N=dnaD$total)
  mat$p <- apply(mat, 1, function(row) binom.test(row[1], row[2])$p.value)
  

  isNull <- rep(NA, n.snps)
  isNull[which.min(mat$p)] <- FALSE
  isNull[which.max(mat$p)] <- TRUE
  datlist <- list(nsnps=n.snps, Y=dnaD[,3], N=dnaD$total)
  
  source('biased_snps_rjags.R')
  
  Y <- dnaD[,3]
  N <- dnaD$total
  probs <- posteriorProbBiased(y=Y, n=N)

  # plot posterior probability biased:
  jpeg(paste(file, 'jpeg', sep='.'))
  par(mar=c(5, 5, 3, 2), cex =1.5)
  hist(probs, col='grey80', main='Posterior Probability of ASE', 
       xlab='probability')
  dev.off()
  # Most SNPs should have P(biased) close to 0 or 1
  
  # filter out biased SNPs:
  biased <- probs > 0.5  # filter out any SNPs with P(biased) > 0.5
  biased[N == 0] <- TRUE  # filter out any SNPs with no data
  datUnbiased <- dnaD[!biased,]
  write.table(datUnbiased, paste(file, 'unbiased.csv', sep='_'), row.names =F)
  
  # plot allelic ratio of all and unbiased SNPs, to visualize the results of filtering.
  jpeg(paste(file, 'unbiased.jpeg', sep='_'), width=960, height = 480)
  par(mar=c(5, 5, 3, 2), cex =1.5, mfcol = c(1,2))
  hist(dnaD[,3]/(dnaD[,3] + dnaD[,4]), col='grey80', main='Allelic ratio (all SNPs)', ylab = "Frequency",
       xlab='Proportion of the homeologue A')
  hist(datUnbiased[,3]/(datUnbiased[,3] + datUnbiased[,4]), col='grey80', main='Allelic ratio (unbiased SNPs)', ylab = "Frequency",
       xlab='Proportion of the homeologue A')
  dev.off()
  message('\n')
  message(paste(file, 'processed\n', sep=' '))
}
