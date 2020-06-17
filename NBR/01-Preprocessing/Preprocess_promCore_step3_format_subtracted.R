##############################################################
## Step 3: Investigate discontinous regions due to subtract ##
##############################################################

bedtools_subtract <- fread("gc19_pc.promCore_transcripts.sorted.merged2.sorted.tier3_subtract.bed")
colnames(bedtools_subtract) <- c("chr", "start","end","name", "score", "strand")
bedtools_subtract$width <- bedtools_subtract$end - bedtools_subtract$start

setdiff(transdf$NBR_name,bedtools_subtract$name)
## "CCNE1::ENSG00000105173.9::TS_1"  "JAK1::ENSG00000162434.7::TS_2"   "JAK2::ENSG00000096968.8::TS_2"   "KIT::ENSG00000157404.11::TS_1"  
## "KRAS::ENSG00000133703.7::TS_3"   "MDM2::ENSG00000135679.17::TS_1"  "VEGFA::ENSG00000112715.16::TS_4"

## get all duplicated
dup <- unique(as.vector(bedtools_subtract$name[duplicated(bedtools_subtract$name) == TRUE]))
length(dup)
# 5625

dups <- bedtools_subtract[which(bedtools_subtract$name %in% dup),]

sp <- split(dups, f=as.vector(dups$name))

spls <- list()
for(i in 1:length(sp)){
  widthsum <- sum(sp[[i]]$width)
  spls[[i]] <- c(as.vector(unique(sp[[i]]$name)),widthsum)
}

dupsdf <- do.call(rbind, spls)
dupsdf <- cbind(dupsdf, "Discontinous")
colnames(dupsdf) <- c("name", "width","type")
dupsdf <- as.data.frame(dupsdf)
dupsdf <- sapply(dupsdf, as.character )

notdups <- cbind(bedtools_subtract[-which(bedtools_subtract$name %in% dup),c("name", "width")], "Continuous")
colnames(notdups) <- c("name", "width","type")
notdups <- sapply( notdups, as.character )

all <- rbind(dupsdf,notdups)
alldf <- as.data.frame(all, stringsAsFactors = F)
alldf$width <- as.numeric(alldf$width)

table(alldf$type)
#Continuous Discontinous 
#47867         5625 
##5625/(47867+5625) * 100 - 10.49% discontinous

###############################################################
### Step4: Create bed12 to work with buildrefnbr_oldoutput.R ##
###############################################################

#################
## Transcripts ##
#################

## remove those with width < 2 since can't get trinucleotides from those
width <- bedtools_subtract$end - bedtools_subtract$start
length(which(width < 3))
# 660

df <- bedtools_subtract[-which(width < 3),]

## combine discontinous regions for transcripts 
sp <- split(df, f=as.vector(df$name))

library(parallel)
bed12 <- mclapply(sp, function(x) data.frame(chr=unique(x$chr),start=min(x$start),end=max(x$end), name=unique(x$name), score=1000,strand=unique(x$strand),thickStart=min(x$start),thickEnd=max(x$end),itemRgb=0,blockCount=nrow(x),blockSizes=paste(x$end-x$start, collapse=",", sep=","), blockStarts=paste(x$start - min(x$start), collapse=",", sep=","), stringsAsFactors = F ),mc.cores = 4)

bed12df <- do.call(rbind, bed12)
rownames(bed12df) <- NULL
head(bed12df)
str(bed12df)

write.table(bed12df, "Regions/gc19_pc.promCore_transcripts.sorted.merged2.sorted.tier3_subtract_bed12.bed", quote=F, row.names=F, col.names=F, sep="\t")
