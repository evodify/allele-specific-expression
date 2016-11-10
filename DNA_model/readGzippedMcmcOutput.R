# read in the output of an MCMC run, which is gzipped
read.mcmc <- function(file=file.choose()) {
	toReturn <- list()
        toReturn$nameOfFileStoringMcmcResults <- file
	# first get comments
	file.gz <- gzfile(file, open="rb")
	comments <- grep("^#", readLines(file.gz), value=TRUE)
	close(file.gz)
	
	toReturn$scriptname <- strsplit(grep("#scriptname", comments, value=TRUE), "=")[[1]][2]
	toReturn$infile <- strsplit(grep("#infile", comments, value=TRUE), "=")[[1]][2]
	toReturn$n.iter <- as.integer(strsplit(grep("#n.iter", comments, value=TRUE), "=")[[1]][2])
	toReturn$thin <- as.integer(strsplit(grep("#thin", comments, value=TRUE), "=")[[1]][2])
	toReturn$features <- strsplit(strsplit(grep("#features", comments, value=TRUE), "=")[[1]][2], " ")[[1]]
	toReturn$n.datasets <- ifelse(any(grepl("#n.datasets", comments)), as.integer(strsplit(grep("#n.datasets", comments, value=TRUE), "=")[[1]][2]), 1)
	if(any(grepl("#a.hat", comments))) {
            toReturn$a.hat <- as.integer(strsplit(grep("#a.hat", comments, value=TRUE), "=")[[1]][2])
        }
        if(any(grepl("#d.hat", comments))) {
            toReturn$d.hat <- as.integer(strsplit(grep("#d.hat", comments, value=TRUE), "=")[[1]][2])
        }
	n.genes <- length(toReturn$features)
	nDataCol <- as.integer(strsplit(comments[grep("ncol", comments)], "=")[[1]][2])
	nDataRow <- as.integer(strsplit(comments[grep("nrow", comments)], "=")[[1]][2])
	stopifnot(nDataRow == toReturn$n.iter/toReturn$thin)
	toReturn <- readInNumericComment(toReturn, comments, "scaling")
	toReturn <- readInNumericComment(toReturn, comments, "acceptanceRates")
	if(toReturn$n.datasets == 1) {
		toReturn <- readInNumericComment(toReturn, comments, "y")
		toReturn <- readInNumericComment(toReturn, comments, "N")
	} else {
		for(i in 1:toReturn$n.datasets) {
			toReturn <- readInNumericComment(toReturn, comments, sprintf("y%d", i))
			toReturn <- readInNumericComment(toReturn, comments, sprintf("N%d", i))
		}
	}
	toReturn <- readInNumericComment(toReturn, comments, "geneIndices")
	timeComment <- grep("#time", comments, value=TRUE)
	if(length(timeComment) > 0) {
		toReturn$time <- as.double(strsplit(timeComment, "=")[[1]][2])
	}
	
	# now get the actual data
	paramNames <- strsplit(strsplit(grep("#names", comments, value=TRUE), "=")[[1]][2], " ")[[1]]
	file.gz <- gzfile(file, open="rb")
	mat <- matrix(scan(file.gz, comment.char="#"), ncol=nDataCol, byrow=TRUE)
	colnames(mat) <- paramNames
	logit_p_names <- grepl("logit_p", paramNames)
	logit_e_names <- grepl("logit_e", paramNames)
	nullComponentNames <- grepl("nullComponent", paramNames)
	otherNames <- !logit_p_names & !logit_e_names & !nullComponentNames
	toReturn$logit_p <- data.frame(mat[,logit_p_names])
	toReturn$logit_e <- data.frame(mat[,logit_e_names])
	nullComp <- mat[,nullComponentNames]
	####toReturn$nullComponent <- data.frame(apply(nullComp, c(1, 2), as.logical))
	toReturn$nullComponent <- as.logical(nullComp)
	dim(toReturn$nullComponent) <- dim(nullComp)
	toReturn$mcmc <- data.frame(mat[,otherNames])
	close(file.gz)
	toReturn
}

readInNumericComment <- function(resultsList, comments, fieldName) {
	matchingComment <- grep(sprintf("#%s", fieldName), comments, value=TRUE)
	if(length(matchingComment) > 0) {
		resultsList[[fieldName]] <- as.numeric(strsplit(strsplit(matchingComment, "=")[[1]][2], " ")[[1]])
	}
	resultsList
}

getMatrixStats <- function(mat, burnin, n.iter) {
	med <- apply(mat[burnin:n.iter,], 2, median)
	upper <- apply(mat[burnin:n.iter,], 2, quantile, probs=0.75)
	lower <- apply(mat[burnin:n.iter,], 2, quantile, probs=0.25)
	list(med=med, upper=upper, lower=lower)
}

logit <- function(val) log(val/(1-val))
inv.logit <- function(val) exp(val)/(1 + exp(val))
