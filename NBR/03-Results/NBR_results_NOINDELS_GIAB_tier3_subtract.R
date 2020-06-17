##############################
## results without INDELS ####
##############################

setwd("/Users/ewilkie/Documents/Work/CCI/Landscape_paper/PCAWG/Github_submission_code/NBR")

library(data.table)

minhyper_results <- fread("03-Results/gc19_pc.promCore_transcripts.sorted.merged2.sorted.GRCh37_alldifficultregions.tier3_subtract_bed12.bed-Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.contigrm-Selection_output_landscape.txt")
minhyper_results_df <- as.data.frame(minhyper_results)
#854 significant results based on Pvalue
length(which(minhyper_results_df$pval_both < 0.05))
## 0 based in qvals
length(which(minhyper_results_df$qval_subs < 0.05))

######################
## Get patient info ##
######################

source("NBR_functions.R")

## intermediate file from NBR to see whether samples have been filtered 
infile_mutinter <- "03-Results/gc19_pc.promCore_transcripts.sorted.merged2.sorted.GRCh37_alldifficultregions.tier3_subtract_bed12.bed-Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.contigrm-Intermediate_mutations.txt"

infile_minhyper_results <- "03-Results/gc19_pc.promCore_transcripts.sorted.merged2.sorted.GRCh37_alldifficultregions.tier3_subtract_bed12.bed-Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.contigrm-Selection_output_landscape.txt"

top10_res <- get_patient_info(infile_mutinter, infile_minhyper_results,nreg=10,glist=NULL,pvaltype="pval_subs_CV")

transdf <- read.table("03-Results/Transcript_Gene_TopExpress_promoter_regions.txt", stringsAsFactors = F, header=T)
candid_res <- get_patient_info(infile_mutinter, infile_minhyper_results,nreg=NULL,glist=transdf$NBR_name, pvaltype="pval_subs_CV")

#################################
## calculate p.adjust manually ##
#################################

cu <- unique(candid_res[,c("region","pval_subs_CV")])

pval_subs_CV_BH <- round(p.adjust(cu$pval_subs_CV,method="BH"),7)

cubu <- cbind(cu$region,pval_subs_CV_BH)
colnames(cubu)[1] <- "region"

candida_bh <- merge(candid_res, cubu, by="region")
candida_bho <- candida_bh[order(candida_bh$pval_subs_CV),]

write.table(candida_bho,"03-Results/gc19_pc.promCore_transcripts.sorted.merged2.sorted.GRCh37_alldifficultregions.tier3_subtract_bed12.bed-Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.contigrm_candidates_results.txt", col.names= T, row.names=F, sep="\t")

