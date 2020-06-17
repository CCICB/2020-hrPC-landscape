###############################
## Remove contiguous regions ##
###############################

mutations_minhyper_min_indel <- fread("Samples/Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.txt")

mutations <- as.data.frame(mutations_minhyper_min_indel)
head(mutations)
dim(mutations)

## remove mutation column burder
mutations <- mutations[,-6]
mutations = mutations[order(mutations$sampleID,mutations$chr,mutations$pos),]

## finds all the positions where the difference is 1
ind = which(diff(mutations$pos)==1)

contiguous <- mutations[unique(sort(c(ind,ind+1))),]
dim(contiguous)
discontiguous <- mutations[-unique(sort(c(ind,ind+1))),]
dim(discontiguous)
head(discontiguous)

write.table(discontiguous, "Landscape_somatic_variants.mut_burden_mb.filter10.without_indels.contigrm.txt", sep="\t", quote=F, row.names=F)
