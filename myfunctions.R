
col2hex <- function(col){
	rgbval <- col2rgb(col)
	return(rgb(rgbval[1,], rgbval[2,], rgbval[3,], maxColorVal = 255))
}

argumentChange <- function(defaultArgs, userArgs){
		userArgs <- list(...)
		defaultArgs[names(userArgs)] <- userArgs
		return(defaultArgs)
}

manhattan <- function(pvals, cumpos = FALSE, alpha = 0.05, transparent = "4D", cols = 3, ...){
	if(is.null(pvals$pval) & !names(pvals)[length(pvals)] %in% c('chrom', 'pos', 'cumpos')) pvalName <- tail(names(pvals), 1) else if(!is.null(pvals$pval)) pvalName <- "pval" else stop("must provide a column labeled 'pval' with pvalues")
	if(!all(c('chrom', 'pos') %in% names(pvals))) stop("dataframe must have columns labeled 'chrom', 'pos' and 'pval'")
	if(cumpos){
		pos <- "pos"
	} else {
		pvals <- cbind(pvals, cumpos = NA)
		maxpos <- 0
		buffer <- 0
		for(i in unique(pvals[,"chrom"])){
		    pvals[pvals[, "chrom"] == i, "cumpos"] <- pvals[pvals[, "chrom"] == i, "pos"] + maxpos + buffer
		    if(buffer == 0) buffer <-  1e6
		    maxpos <- max(pvals[pvals[,"chrom"] == i, "cumpos"])
		}
		pos <- "cumpos"
	}
	if(is.null(pvals$col)){
		if(is.numeric(cols) | is.integer(cols)) cols <- c("dodgerblue4", "skyblue3", "slategray4") 
		chroms <- unique(pvals$chrom)
		nrep <- ceiling(length(chroms) / length(cols))
		snpCols <- data.frame(chrom = chroms, col = col2hex(rep(cols[1:length(cols)], times = nrep))[1:length(chroms)])
		pvals <- merge(pvals, snpCols, by = "chrom")
	} else {
		if(!grepl("^\\#", pvals$col[1])) pvals$col <- col2hex(pvals$col)
	}
	l <- list(...)
	plotArgs <- argumentChange(list(col = pvals$col, bg = paste0(pvals$col, transparent), pch = 21, cex = 1, xaxt = "n", xlab = "Chromosome position", ylab = expression(-log[10](pvalue))), l)	
	do.call(plot, c(list(x = pvals[[pos]], y = pvals[[pvalName]]), plotArgs))
	# plot(pvals[[pos]], pvals$pval, col = pvals$col, bg = paste0(pvals$col, transparent), pch = 21, cex = 1, xaxt = "n")
	bonf <- alpha / nrow(pvals)
	abline(h = -log10(bonf))
	axis(1, at = labPos, labels = names(labPos))
}