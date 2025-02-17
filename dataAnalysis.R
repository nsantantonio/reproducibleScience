# load rrblup library
# install.packages("rrBLUP") # install if not done, only needs to be run one time!
library(rrBLUP)

# source our functions from a script (even better make it a package!)
source("myfunctions.R")

# get path to data. 
local <- FALSE
if(local){
	path <- "/home/nsant/Dropbox/reproducibleScience/" 
} else {
	path <- "./" 
}
load(paste0(path, "SRWobs_2022_WARVA_forPMgwas.RData"))

# make design and incidence matrices for model fit
X <- model.matrix(~Block, data = obsWR)
Z <- model.matrix(~Line -1, data = obsWR)
colnames(Z) <- gsub("^Line", "", colnames(Z))

if (all(colnames(Z) == rownames(Mwr))) {
	ZM <- Z %*% Mwr
} else {
	stop("colnames of Z dont match row names of M!")
}

####################################################
# get cumulative positions of markers for plotting #
####################################################

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


#############################
# niave GWAS model with OLS #
#############################

trait <- "PowderyMildew"

obsWRzm <- cbind(obsWR, as.data.frame(ZM))
pval = NULL
eff = NULL
for(i in colnames(ZM)){
	# i <- colnames(ZM)[1]
	fiti <- lm(formula(paste0(trait, "~ Block + ", i)), data = obsWRzm)
	pval[i] <- -log10(coef(summary(fiti))[i, "Pr(>|t|)"])
}


# combine snp positions with pvalues
pvals <- cbind(snpInfo, pval = pval)

# user defined function to make a manhattan plot
manhattan(pvals) # use user defined colors in data.frame 'pvals' 
manhattan(pvals[!names(pvals) %in% "col"]) # default colors from function

# # just in case you want to be extra garish
manhattan(pvals[!names(pvals) %in% "col"], cols = c("#861F41", "#E87722", "#75787b"))

# make it a pdf
pdf("naiveGWAS.pdf", width = 14)
manhattan(pvals[!names(pvals) %in% "col"])
dev.off()



##################################################
# GWAS model correcting for population structure #
##################################################

ZMt <- as.data.frame(t(ZM)) # to make 'geno' for GWAS()
colnames(ZMt) <- obsWRzm$Line 

ZKZt <- Z %*% Kwr %*% t(Z) # expand K onto data structure
rownames(ZKZt) <- colnames(ZKZt) <- obsWRzm$Line

# loop over traits, pvals saved
traits <- c("PowderyMildew", "GrainYield")#, "TestWeight", "HeadingDate", "Height")
# triats <- traits[1]
resultList <- list()
for (i in traits){
	pheno <- obsWRzm[c("Line", "Block", i)]
	geno <- cbind(snpInfo[c("rs", "chrom", "pos")], ZMt)
	resultList[[i]] <- GWAS(pheno = pheno, geno = geno, fixed = c("Block"), K = ZKZt, plot = FALSE)
}



# plot results
pdf("p3dGWASmanhattan.pdf")
	manhattan(resultList[[1]], main = traits[1])
	manhattan(resultList[[2]], main = traits[2])
dev.off()

# add function to make qq plot in class. 

# print table of markers with significant effects for trait 1
alpha <- 0.05
bonferroni <- alpha / ncol(Mwr)
QTL <- resultList[[1]][resultList[[1]][,4] >= -log10(bonferroni),]
names(QTL) <- c("SNP", "Chromosome", "Position", "-log(P-value)")
print(QTL)


# below will make a latex table for a latex document

# library(devtools) # if you dont have devtools installed already.
# install_github("summaryTools",username="nsantantonio")
# library(summaryTools)
# xtable2(QTL, con = "qtlTable.tex", include.rownames = FALSE)