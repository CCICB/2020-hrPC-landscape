#####################################
## Preprocess gc19_pc.promCore.bed ##
#####################################


setwd("/Users/ewilkie/Documents/Landscape_paper/PCAWG/NBR/NBR_EEW2")

##############################################
## step 1 : convert PromCore bed12 to bed 6 ##
##############################################

library(data.table)
library(GenomicRanges)

## this file if from the original NBR_TO_SHARE
bed <- fread("Regions/gc19_pc.promCore.bed")

bed$V11 <- gsub(",$", "", bed$V11)
bed$V12 <- gsub(",$", "", bed$V12)

bed <- as.data.frame(bed)

regiondf = vector(mode="list", nrow(bed))
for (j in 1:nrow(bed)) {
  start = bed[j,2]
  segment_sizes = as.numeric(strsplit(bed[j,11],split=",")[[1]])
  segment_starts = as.numeric(strsplit(bed[j,12],split=",")[[1]])
  name <- gsub("gc19_pc.promCore::gencode::", "", bed[j,4])
  regiondf[[j]] = data.frame(chr=gsub("chr","",bed[j,1]), start=start+segment_starts, end=start+segment_starts+segment_sizes, name=name, score=1000, strand=bed[j,6], stringsAsFactors=F)
}

names(regiondf) = sapply(regiondf, function(x) x[1,4])

## add transcript name
regionls <- list()
for(i in 1:length(regiondf)){
  gene <- regiondf[[i]]
  gene$name <- paste(gene$name, "::TS_", seq(1:nrow(gene)), sep="")
  regionls[[i]] <- gene
}

regiondf_collapse <- do.call(rbind, regionls)

##########################################
## Check all candidates are there - yep ##
##########################################

transdf <- read.table("Candidates/Data/Transcript_Gene_TopExpress_promoter_regions_for_Mustafa.txt", stringsAsFactors = F, header=T)

setdiff(transdf$NBR_name,regiondf_collapse$name)
#character(0)

##############################################
### Check and merge overlapping transcripts ##
##############################################

gr <- makeGRangesFromDataFrame(regiondf_collapse,keep.extra.columns=T)
fo <- findOverlaps(gr, type = "any", select = "all", ignore.strand=F)
fodf <- as.data.frame(fo, stringsAsFactors=F)
fodfu <- fodf[which(fodf[,1] != fodf[,2]),]

QH <- regiondf_collapse[fodfu[,"queryHits"],c(1:4)]
colnames(QH)[1:4] <- paste("QH_", colnames(QH)[1:4], sep="")
QH_sp <- strsplit(QH$QH_name, "::")
QH_name_df <- do.call(rbind, QH_sp)
QH <- cbind(QH, QH_name_df)
colnames(QH)[5:7] <- c("QH_gene", "QH_ENS", "QH_TS")
head(QH)

SH <- regiondf_collapse[fodfu[,"subjectHits"],c(1:4)]
colnames(SH)[1:4] <- paste("SH_", colnames(SH)[1:4], sep="")
SH_sp <- strsplit(SH$SH_name, "::")
SH_name_df <- do.call(rbind, SH_sp)
SH <- cbind(SH, SH_name_df)
colnames(SH)[5:7] <- c("SH_gene", "SH_ENS", "SH_TS")
head(SH)

AH <- cbind(QH,SH)
head(AH)
dim(AH)
## 1524 overlaps

length(unique(c(QH$QH_name, SH$SH_name))) 
## 1311

## using the entire string there are no matches
AH[which(AH$QH_name == AH$SH_name),]

## are any within the same genes - based on ENS - no
AH[which(AH$QH_ENS == AH$SH_ENS),]

##this is based on gene name - seems like multiple gene names with differen ENS IDS
AH[which(AH$QH_gene == AH$SH_gene),]

overlaps <- list()
overlaps$df <- AH
overlaps$impacted_TS <- unique(c(QH$QH_name, SH$SH_name))

#  overlaps in the transcripts of interest
## no overlaps with candidates
which(transdf$NBR_name %in% overlaps$impacted_TS)
# integer(0)

## transcripts indicate that each line is a transcript promoter
write.table(regiondf_collapse, "Regions/gc19_pc.promCore_transcripts.bed", row.names=F, col.names=F, quote = F, sep="\t")
