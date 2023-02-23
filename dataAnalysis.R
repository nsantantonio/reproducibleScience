load("/home/nsant/Dropbox/reproducibleScience/SRWobs_2022_WARVA_forPMgwas.RData")

X <- model.matrix(~Block, data = obsWR)
Z <- model.matrix(~Line -1, data = obsWR)
colnames(Z) <- gsub("^Line", "", colnames(Z))

if (all(colnames(Z) == rownames(Mwr))) {
	ZM <- Z %*% Mwr
} else {
	stop("colnames of Z dont match row names of M!")
}

################ 
# for plotting #
################


pos <- snpInfo
pos <- cbind(pos, cumpos = NA)
maxpos <- 0
buffer <- 0
for(i in unique(pos[,"chrom"])){
    pos[pos[, "chrom"] == i, "cumpos"] <- pos[pos[, "chrom"] == i, "pos"] + maxpos + buffer
    if(buffer == 0) buffer <-  1e6
    maxpos <- max(pos[pos[,"chrom"] == i, "cumpos"])
}

labPos <- tapply(pos$cumpos, pos$chrom, function(x) (max(x) - min(x)) / 2 + min(x))


###############
# niave model #
###############

trait <- "PowderyMildew"

obsWRzm <- cbind(obsWR, as.data.frame(ZM))
pval = NULL
eff = NULL
for(i in colnames(ZM)){
	# i <- colnames(ZM)[1]
	fiti <- lm(formula(paste0(trait, "~ Block + ", i)), data = obsWRzm)
	pval[i] <- -log10(coef(summary(fiti))[i, "Pr(>|t|)"])
}



# make function to make manhattan plot
pvals <- cbind(snpInfo, pval = pval)

col2hex <- function(col){
	rgbval <- col2rgb(col)
	return(rgb(rgbval[1,], rgbval[2,], rgbval[3,], maxColorVal = 255))
}

manhattan <- function(pvals, cumpos = FALSE, alpha = 0.05, transparent = "4D", cols = 3, ...){
	argumentChange <- function(defaultArgs, userArgs){
		userArgs <- list(...)
		defaultArgs[names(userArgs)] <- userArgs
		return(defaultArgs)
	}
	col2hex <- function(col){
		rgbval <- col2rgb(col)
		return(rgb(rgbval[1,], rgbval[2,], rgbval[3,], maxColorVal = 255))
	}
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
# manhattan(pvals) # default colors

# # if you are feeling particularly patriotic
# manhattan(pvals, cols = c("firebrick", "royalblue4"))

# # just in case you want to be extra garish
# manhattan(pvals, cols = c("#861F41", "#E87722", "#75787b"))

# make it a pdf
pdf("naiveGWAS.pdf", width = 14)
manhattan(pvals)
dev.off()





###############
# GWAS model #
###############
library(rrBLUP)
ZMt <- as.data.frame(t(ZM)) # to make 'geno' for GWAS()
colnames(ZMt) <- obsWRzm$Line 

ZKZt <- Z %*% Kwr %*% t(Z) # expand K onto data structure
rownames(ZKZt) <- colnames(ZKZt) <- obsWRzm$Line

# note this is super repetative and generally unneccessary. We were just rushing to analyze multiple traits before lab meeting. 
# really, it should be looped over traits, and the pvals saved, etc. 
traits <- c("PowderyMildew", "GrainYield")#, "TestWeight", "HeadingDate", "Height")
# triats <- traits[1]
resultList <- list()
for (i in traits){
	pheno <- obsWRzm[c("Line", "Block", i)]
	geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
	resultList[[i]] <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt, plot = FALSE)
}

pdf("p3dGWASmanhattan.pdf")
manhattan(resultList[[1]])
manhattan(resultList[[2]])
# manhattan(resultList[[3]])
# manhattan(resultList[[4]])
dev.off()

# pheno <- obsWRzm[c("Line", "Block", "GrainYield")]
# ZMt <- as.data.frame(t(ZM))
# colnames(ZMt) <- pheno$Line
# geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
# geno[1:10, 1:20]
# ZKZt <- Z %*% Kwr %*% t(Z) 

# rownames(ZKZt) <- colnames(ZKZt) <- pheno$Line
# # all(pheno$Line %in% colnames(ZKZt))
# rrgwas <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt)


# pheno <- obsWRzm[c("Line", "Block", "TestWeight")]
# ZMt <- as.data.frame(t(ZM))
# colnames(ZMt) <- pheno$Line
# geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
# geno[1:10, 1:20]
# ZKZt <- Z %*% Kwr %*% t(Z) 

# rownames(ZKZt) <- colnames(ZKZt) <- pheno$Line
# # all(pheno$Line %in% colnames(ZKZt))
# rrgwas <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt)

# pheno <- obsWRzm[c("Line", "Block", "HeadingDate")]
# ZMt <- as.data.frame(t(ZM))
# colnames(ZMt) <- pheno$Line
# geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
# geno[1:10, 1:20]
# ZKZt <- Z %*% Kwr %*% t(Z) 

# rownames(ZKZt) <- colnames(ZKZt) <- pheno$Line
# # all(pheno$Line %in% colnames(ZKZt))
# rrgwas <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt)

# pheno <- obsWRzm[c("Line", "Block", "Height")]
# ZMt <- as.data.frame(t(ZM))
# colnames(ZMt) <- pheno$Line
# geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
# geno[1:10, 1:20]
# ZKZt <- Z %*% Kwr %*% t(Z) 

# rownames(ZKZt) <- colnames(ZKZt) <- pheno$Line
# # all(pheno$Line %in% colnames(ZKZt))
# rrgwas <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt)