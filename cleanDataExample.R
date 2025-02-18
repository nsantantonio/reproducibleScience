local = TRUE
if(local){
	genoDataPath <- "/home/nsant/Dropbox/genotypingGBS2022/"
	phenoDataPath <- "/home/nsant/Dropbox/VTsmallgrainsAnalysis2022/analysis/compiledDataSets/"
} else {
	genoDataPath <- "dataFiles/"
	phenoDataPath <- "dataFiles/"
}

# read in K (additive covariance) and M (marker score) matrices
# K <- read.table(paste0(genoDataPath, "SRWaddCov2022.txt"), check.names = FALSE)
M <- read.table(paste0(genoDataPath, "SRWmarkersConsensus2022.txt"), check.names = FALSE)
M <- as.matrix(M)

vanRaden1 <- function(M){
	Z <- scale(M, scale = FALSE)
	p <- attributes(Z)[["scaled:center"]] / 2
	ZZt <- tcrossprod(Z)
	ZZt / (2 * crossprod(p, 1-p)[[1]])
}

K <- vanRaden1(K)

# get snp information
snps <- colnames(M)
head(snps)
snpInfo <- data.frame(rs = snps, chrom = gsub("^S|_.*", "", snps), pos = as.integer(gsub(".*_", "", snps)))
head(snpInfo)

# read in pheno data
obs <- read.csv(paste0(phenoDataPath, "SRWObs_2022_withURBIL.csv"))
# get rid of units on 
names(obs) <- gsub("\\..*", "", names(obs))
head(obs)

# keep observations from Warsaw location
obsWR <- droplevels(obs[obs$Location %in% c("WARVA"), ])

# subset matrices 
Kwr <- K[rownames(K) %in% obsWR$Line, colnames(K) %in% obsWR$Line]
Mwr <- M[rownames(M) %in% obsWR$Line, ]

# remove markers with maf < 0.01
p <- colMeans(Mwr) / 2 
# hist(p)
tooLowMaf <- p < 0.01
sum(tooLowMaf)
Mwr <- Mwr[, !tooLowMaf]
snpInfo <- snpInfo[snpInfo$rs %in% colnames(Mwr),]

# check that all lines in pheno have been genotyped
all(colnames(Kwr) %in% obsWR$Line)
all(obsWR$Line %in% colnames(Kwr))
notGeno <- unique(obsWR$Line[!obsWR$Line %in% colnames(Kwr)])
length(notGeno)
# drop ungenotyped lines from pheno data
obsWR <- droplevels(obsWR[!obsWR$Line %in% notGeno, ])
all(obsWR$Line %in% colnames(Kwr))

# check that rows of marker matrix matches relationship matrix
all(colnames(Kwr) == rownames(Mwr))

# make Line, Block, Location factors (location not necessary since there is only one location for this example)
# set levels equal to rows/columns of relationship matrix / rows of marker matrix
obsWR$Line <- factor(obsWR$Line, levels = colnames(Kwr))
obsWR$Block <- factor(obsWR$Block, levels = unique(obsWR$Block))
obsWR$Location <- factor(obsWR$Location, levels = unique(obsWR$Location))

# check that levels of Lines in pheno data same 
all(levels(obsWR$Line) == colnames(Kwr))

# subset marker information table
snpInfo <- snpInfo[snpInfo$rs %in% colnames(Mwr),]
snpInfo$chrom <- gsub("^S", "", snpInfo$chrom)
snpCols <- data.frame(chrom = unique(snpInfo$chrom), col = rep(c("darkblue", "blue", "lightblue"), times = 7))
snpInfo <- merge(snpInfo, snpCols, by = "chrom")

# save relavent objects in an RData file for use in analysis
save(obsWR, Kwr, Mwr, snpInfo, file = "SRWobs_2022_WARVA_forPMgwas.RData")
