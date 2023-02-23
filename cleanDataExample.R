K <- read.table("/home/nsant/Dropbox/genotypingGBS2022/SRWaddCov2022.txt", check.names = FALSE)
M <- read.table("/home/nsant/Dropbox/genotypingGBS2022/SRWmarkersConsensus2022.txt", check.names = FALSE)
K <- as.matrix(K)
M <- as.matrix(M)

snps <- colnames(M)
snpInfo <- data.frame(rs = snps, chrom = gsub("^S|_.*", "", snps), pos = as.integer(gsub(".*_", "", snps)))
head(snpInfo)

obs <- read.csv("/home/nsant/Dropbox/VTsmallgrainsAnalysis2022/analysis/compiledDataSets/SRWObs_2022_withURBIL.csv")
names(obs) <- gsub("\\..*", "", names(obs))
head(obs)


obsWR <- droplevels(obs[obs$Location %in% c("WARVA"), ])

Kwr <- K[rownames(K) %in% obsWR$Line, colnames(K) %in% obsWR$Line]
Mwr <- M[rownames(M) %in% obsWR$Line, ]

p <- colMeans(Mwr) / 2 
# hist(p)
tooLowMaf <- p < 0.01
sum(tooLowMaf)
Mwr <- Mwr[, !tooLowMaf]
snpInfo <- snpInfo[snpInfo$rs %in% colnames(Mwr),]

all(colnames(Kwr) %in% obsWR$Line)
all(obsWR$Line %in% colnames(Kwr))
notGeno <- unique(obsWR$Line[!obsWR$Line %in% colnames(Kwr)])
length(notGeno)
obsWR <- droplevels(obsWR[!obsWR$Line %in% notGeno, ])
all(obsWR$Line %in% colnames(Kwr))

all(colnames(Kwr) == rownames(Mwr))

obsWR$Line <- factor(obsWR$Line, levels = colnames(Kwr))
obsWR$Block <- factor(obsWR$Block, levels = unique(obsWR$Block))
obsWR$Location <- factor(obsWR$Location, levels = unique(obsWR$Location))

all(levels(obsWR$Line) == colnames(Kwr))

snpInfo <- snpInfo[snpInfo$rs %in% colnames(Mwr),]
snpInfo$chrom <- gsub("^S", "", snpInfo$chrom)
snpCols <- data.frame(chrom = unique(snpInfo$chrom), col = rep(c("darkblue", "blue", "lightblue"), times = 7))
snpInfo <- merge(snpInfo, snpCols, by = "chrom")


save(obsWR, Kwr, Mwr, snpInfo, file = "SRWobs_2022_WARVA_forPMgwas.RData")
