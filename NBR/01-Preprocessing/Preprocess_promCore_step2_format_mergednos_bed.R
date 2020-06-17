##############################################
## check hybrid bed for overlaps and format ##
##############################################

setwd("/Users/ewilkie/Documents/Landscape_paper/PCAWG/NBR/NBR_EEW2/")
hybrid_bed <- fread("Regions/gc19_pc.promCore_transcripts.sorted.merged.nos.bed")
hybrid_bed <- as.data.frame(hybrid_bed)
colnames(hybrid_bed) <- c("chr", "start", "end","name", "strand")
head(hybrid_bed)

## check candidates
transdf <- read.table("Candidates/Data/Transcript_Gene_TopExpress_promoter_regions_for_Mustafa.txt", stringsAsFactors = F, header=T)

hybrid_bed[grep(paste(transdf$NBR_name, collapse="|"), hybrid_bed$name),]
hybrid_bed_candid <- hybrid_bed[which(hybrid_bed$name %in% transdf$NBR_name),]

hybrid_bed_names <- unique(unlist(strsplit(hybrid_bed$name , ",")))
hybrid_bed_names[which(hybrid_bed_names %in% transdf$NBR_name)]

## no difference, they're all there
setdiff(transdf$NBR_name,hybrid_bed_candid$name )
## character(0)

## ignore strand since NBR and mustations are strand agnostic
## do formatting to remove these 
unique(hybrid_bed$strand)
#[1] "+"           "-"           "-,+"         "+,-"         "+,+"         "+,-,+"       "-,+,-"       "-,-,-"       "-,-"         "-,+,+,+"     "+,+,+"       "-,+,-,+"    
#[13] "-,+,+"       "-,-,+"       "+,-,-"       "+,+,-"       "-,-,-,+"     "-,-,-,+,-,-" "-,-,+,-"     "-,+,+,-"  

hybrid_bed_form <- hybrid_bed[,1:4]

gr <- makeGRangesFromDataFrame(hybrid_bed_form,keep.extra.columns=T)
fo <- findOverlaps(gr, type = "any", select = "all", ignore.strand=F)
fodf <- as.data.frame(fo, stringsAsFactors=F)
fodfu <- fodf[which(fodf[,1] != fodf[,2]),]

print(fodfu)
## 0 overlaps

merge_out <- cbind(hybrid_bed_form[,1:4], 1000, ".")
write.table(merge_out, "Regions/gc19_pc.promCore_transcripts.sorted.merged2nos.bed", quote=F, row.names=F, col.names=F, sep="\t")
