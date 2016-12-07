# command: 
#
# Rscript DNAmodel_exctract_results.R input1.gz input2.gz input3.gz

source('readGzippedMcmcOutput.R')

ff <- commandArgs(trailingOnly = T) # a list of files to process

outputList <- data.frame(matrix(nrow = c(1:length(ff)), ncol = 3))
colnames(outputList) <-	c('name', 'a.hat', 'd.hat')

# loop through all the files
for (i in c(1:length(ff))) {
  result <- read.mcmc(ff[i])
  message(paste(ff[i], 'loaded\n', sep=' '))
  n.iter <- result$n.iter/result$thin
  burnin <- 0.1*n.iter
  a.hat <- median(exp(result$mcmc$logA[burnin:n.iter]))
  d.hat <- median(exp(result$mcmc$logD[burnin:n.iter]))
  outputList[i,] <- c(ff[i], a.hat, d.hat)
  # plot
  jpeg(paste(ff[i], 'jpeg', sep='_'), width=960, height = 480)
  par(mar=c(5, 5, 3, 2), mfcol = c(1,2))
  par(cex=1)
  hist(exp(result$mcmc$logA[burnin:n.iter]), breaks = 50,
       main='Distribution of a', xlab='a')
  abline(v = a.hat, col ='red', lwd = 2)
  hist(exp(result$mcmc$logD[burnin:n.iter]), breaks = 50,
       col='grey80', main='Distribution of d', xlab='d')
  abline(v = d.hat, col ='red', lwd = 2)
  dev.off()
  message(paste(ff[i], 'processed\n', sep=' '))
}

write.table(outputList, "a.hat_d.hat.txt", row.names = F)
