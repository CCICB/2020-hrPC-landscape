##################################
## Prepare regions to be remove ##
##################################

## giab genome file
ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v2.0/GRCh37/union/GRCh37_alldifficultregions.bed.gz

## combine files since bedtools merge only works on a single file
cat GRCh37_alldifficultregions.bed tier3.bed > GRCh37_alldifficultregions.tier3.bed

## sort files
sort -k1,1 -k2,2n GRCh37_alldifficultregions.tier3.bed > GRCh37_alldifficultregions.tier3.sorted.bed

## merge files to combine overlapping regions
bedtools merge -i GRCh37_alldifficultregions.tier3.sorted.bed > GRCh37_alldifficultregions.tier3.sorted.merged.bed

## sort again (probably not necessary)
sort -k1,1 -k2,2n GRCh37_alldifficultregions.tier3.sorted.merged.bed > GRCh37_alldifficultregions.tier3.sorted.merged.sorted.bed

#########################################
## Preparing promCore transcript files ##
#########################################

## get promCore regions from:
https://dcc.icgc.org/releases/PCAWG/drivers/metadata/genomic_intervals_lists
file is: gc19_pc.promCore.bed 

## run script with above file as input
Preprocess_promCore_step1_convert_to_transcripts.R

## sort output from above script
sort -k1,1 -k2,2n gc19_pc.promCore_transcripts.bed > gc19_pc.promCore_transcripts.sorted.bed

# merge transcripts
bedtools merge -c 4,6 -o collapse -i gc19_pc.promCore_transcripts.sorted.bed > gc19_pc.promCore_transcripts.sorted.merged.nos.bed

## runs script with merged output to format file
Preprocess_promCore_step2_format_mergednos_bed.R

## sort above output 
sort -k1,1 -k2,2n gc19_pc.promCore_transcripts.sorted.merged2.bed > gc19_pc.promCore_transcripts.sorted.merged2.sorted.bed 

## subtract difficult to sequence regions from above output
bedtools subtract -a gc19_pc.promCore_transcripts.sorted.merged2nos.sorted.bed -b GRCh37_alldifficultregions.tier3.sorted.merged.sorted.bed > gc19_pc.promCore_transcripts.sorted.merged2.sorted.GRCh37_alldifficultregions.tier3_subtract.bed

## format above 
Preprocess_promCore_step3_format_subtracted.R

## Run buildrefnbr_oldoutput.R on NCI to generate background regions

## on NCI
## module load  R/4.0.0
source("buildrefnbr_oldoutput.R")
## genome.fa contains MT chr therefore use onlychrs subset
chrs <- c(seq(1:22), "X", "Y")
buildrefnbr("gc19_pc.promCore_transcripts.sorted.merged2.sorted.tier3_subtract_bed12.bed", "genomes/hs37d5/genome.fa", onlychrs =chrs)