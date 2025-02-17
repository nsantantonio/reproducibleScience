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

source("myfunctions.R")

manhattan(pvals) # use user defined colors in data.frame 'pvals' 
manhattan(pvals[!names(pvals) %in% "col"]) # default colors from function

# # if you are feeling particularly patriotic
manhattan(pvals[!names(pvals) %in% "col"], cols = c("firebrick", "royalblue4"))

# # just in case you want to be extra garish
manhattan(pvals[!names(pvals) %in% "col"], cols = c("#861F41", "#E87722", "#75787b"))

# make it a pdf
pdf("naiveGWAS.pdf", width = 14)
manhattan(pvals)
dev.off()


# this is what is different!




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