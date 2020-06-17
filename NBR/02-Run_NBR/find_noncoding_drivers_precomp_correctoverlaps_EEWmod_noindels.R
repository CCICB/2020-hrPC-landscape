# Inigo Martincorena - 2014
# Non-coding driver detection with precomputed trinucleotide composition and using the
# local density of putatively neutral mutations as covariates ("t") in the framework of a
# negative binomial regression.
#
#
##########################################################################################
# Instructions for Mark Cowley et al:
# Required R packages: "GenomicRanges", "Rsamtools", "MASS"
# Mutations file must be a tab-separated matrix, with columns: "sampleID","chr","pos","ref","mut"
# Example:
#  sampleID	chr	pos	ref	mut
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1230448	G	A
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1609723	C	T
#  0009b464-b376-4fbc-8a56-da538269a02f	1	1903276	C	T
#  0009b464-b376-4fbc-8a56-da538269a02f	1	2574999	C	T
#  ...
#
# Use chromosome names without "chr", and use build 37 (hg19) coordinates
#
# Modify PATH_TO_DATA_AND_GENOME and genomeFile to match your paths and file names
# PATH_TO_DATA_AND_GENOME: path to GRanges_driver_regions.RData, 
#                                  Trinucfreqs_within_100kb_bins.txt 
#                                  Neutral_regions_within_100kb_bins.txt
#                                  Regions/
#
# Examples:
# miRNA promoters:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" "Regions/mirna.prom.bed" 
# miRNA precursors:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" Regions/mirna.pre.bed
# miRNA mature:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" Regions/mirna.mat.bed
# Protein-coding genes promoters:
#    Rscript find_noncoding_drivers_precomp_correctoverlaps.R "mutations.txt" "Regions/gc19_pc.promCore.bed" 
#
# Input files:
# - mutations.txt: Table of mutations in the 5 column format (sampleID\tchr\tpos\tref\tmut)
# - trinucfreqs_prefix: Region_name\tIntervals_string\tTrinucleotide composition
#
# Test:
#     Rscript find_noncoding_drivers_precomp_correctoverlaps.R  example_Thy-AdenoCa.muts.5cols Regions/gc19_pc.promCore.bed
#   Main output:
#     Regions/gc19_pc.promCore.bed-Selection_output.txt
#       region	chr	start	end	exp_subs	exp_indels	obs_subs	obs_indels	local_t	local_t_indels	pval_subs	pval_indels	pval_both	qval_subs	qval_indels	qval_both	exclude_forfit	cv_predicted_subs	cv_predicted_indels	pval_subs_CV	pval_indels_CV	pval_both_CV	qval_subs_CV	qval_indels_CV	qval_both_CV	qval_both_CV_tiered	obsexp_subs_mle	obsexp_subs_low	obsexp_subs_high	obsexp_indels_mle	obsexp_indels_low	obsexp_indels_high
#       gc19_pc.promCore::gencode::TERT::ENSG00000164362.14	5	1295105	1295362	0.009328714	0.000690961	12	0	1.360598307	0.838612166	7.35E-21	1	0	1.48E-16	1	0	TRUE	0.008465574	0.000666505	6.23E-26	1	0	1.26E-21	1	0	0	1417.505812	721.8452856	2562.228503	0	0	2935.603388
#       gc19_pc.promCore::gencode::PLEKHS1::ENSG00000148735.10	10	115511013	115534874	0.043607871	0.004911716	4	0	0.932850568	1.141215427	2.48E-05	1	0.000287769	0.249859112	1	1	FALSE	0.042695237	0.005474887	2.45E-06	1	3.41E-05	0.024713939	1	0.343976004	0.680680582	93.68726473	28.34270653	229.041039	0	0	357.5185204
#       gc19_pc.promCore::gencode::OR4F5::ENSG00000186092.4	1	68891	69090	0.005602044	0.000535629	0	0	0.721785491	0.552476592	1	1	1	1	1	1	FALSE	0.005694223	0.00045065	1	1	1	1	1	1	1	0	0	343.7822825	0	0	4341.780179
#       gc19_pc.promCore::gencode::AL627309.1::ENSG00000237683.5	1	139310	139579	0.008023746	0.000723099	0	0	0.703778928	0.537202648	1	1	1	1	1	1	FALSE	0.008181887	0.000603953	1	1	1	1	1	1	1	0	0	239.1610734	0	0	3239.682298
#       gc19_pc.promCore::gencode::OR4F29::ENSG00000235249.1	1	367440	367658	0.004864302	0.000586513	0	0	0.727080947	0.509529225	1	1	1	1	1	1	FALSE	0.004939695	0.000483438	1	1	1	1	1	1	1	0	0	396.0928871	0	0	4047.297819
#
#   Columns of interest: mainly "qval_subs_CV","qval_indels_CV","qval_both_CV","qval_both_CV_tiered"
#   The corresponding pvalues can be of interest if you want to do restricted hypotheses testing on a group of
#    genes or regions.
# 
#   If the script is going to be run for different sets of mutations and the same regions, 
#    consider moving the output files to avoid overwritting. Or modify the "write.table" lines
#    to choose your own output naming choices.
#
##########################################################################################


# OPTIONS

PATH_TO_DATA_AND_GENOME = "/g/data/tx70/private/ew8662/PCAWC/NBR_EEW2" # Set this path accordingly.
#genomeFile              = paste(PATH_TO_DATA_AND_GENOME,"/genome.fa",sep="") # Set genome file accordingly
genomeFile              = "/g/data/tx70/genomes/hs37d5/genome.fa" # Set genome file accordingly

max_num_muts_perRegion_perSample = 2
unique_indelsites_FLAG = 1
load(paste(PATH_TO_DATA_AND_GENOME,"/Background/GRanges_driver_regions.RData",sep="")) # Regions to be excluded for the background model fit (GRanges object called "gr_drivers")
trinucfreq_neutralbins_file = paste(PATH_TO_DATA_AND_GENOME,"/Background/Trinucfreqs_within_100kb_bins.txt",sep="")
regions_neutralbins_file = paste(PATH_TO_DATA_AND_GENOME,"/Background/Neutral_regions_within_100kb_bins.txt",sep="")


# Environment
library("GenomicRanges")
library("Rsamtools")
library("MASS")
args = commandArgs(TRUE)
mutations_file = args[1]
trinucfreqs_prefix = args[2]

## put results in seperate folder
on1 <- gsub(".*/", "", trinucfreqs_prefix)
on2 <- gsub(".*/", "", mutations_file)
on3 <- gsub("\\.txt", "", on2)
output_prefix = paste(PATH_TO_DATA_AND_GENOME, "/Results/", on1,"-",on3, sep="")

## 1. Loading the mutations and extracting the trinucleotide context

cat("Loading mutations...\n")
chr_list = c(1:22,"X","Y")
mutations = read.table(mutations_file, header=0, sep = "\t", stringsAsFactors=F) # Loading the file
mutations = mutations[,1:5]
colnames(mutations) = c("sampleID","chr","pos","ref","mut")
mutations = mutations[as.vector(mutations$chr) %in% chr_list,]
mutations$pos = as.numeric(as.vector(mutations$pos))


# Extracting the trinucleotide context of each substitution

cat("Indels...\n")
indels_pos = as.vector(mutations$ref)=="-" | as.vector(mutations$mut)=="-" | nchar(as.vector(mutations$ref))!=1 | nchar(as.vector(mutations$mut))!=1
indels = mutations[indels_pos,]
subs = mutations[!indels_pos,]

cat("Trinucleotides...\n")
nt = c("A","C","G","T"); base1 = rep(nt,each=16,times=1); base2 = rep(nt,each=4,times=4); base3 = rep(nt,each=1,times=16)
trinuc_list = paste(base1,base2,base3, sep="")

seqs = scanFa(genomeFile, GRanges(as.vector(subs$chr), IRanges(subs$pos-1, subs$pos+1)))
muts_trinuc = as.vector(seqs)
subs$trinuc_ref = muts_trinuc

if (nrow(indels)>0) { indels$trinuc_ref = NA }

# Annotating unique indels: indels with the same coordinates and same ref>mut event are flagged as duplicates
indels$unique_indels = rownames(indels) %in% rownames(unique(indels[,2:5]))
subs$unique_indels = TRUE

mutations = rbind(subs,indels)
mutations = mutations[order(mutations$sampleID, mutations$chr, mutations$pos),]

# Loading the precomputed region files
cat("Loading target regions...\n")

target_regions = read.table(sprintf("%s.regions",trinucfreqs_prefix), header=1, sep="\t")
neutral_regions = read.table(regions_neutralbins_file, header=1, sep="\t")
trinucfreq_table = read.table(sprintf("%s.txt",trinucfreqs_prefix), header=1, sep="\t")
trinucfreq_table_neutral = read.table(trinucfreq_neutralbins_file, header=1, sep="\t")

# Creating the 2 global L_matrix
cat("Defining L matrix...\n")

L_matrix_ref = array(0,c(64,4)); rownames(L_matrix_ref) = trinuc_list; colnames(L_matrix_ref) = nt
for (j in 1:64) { L_matrix_ref[j,base2[j]] = NA }

trin_freqs = colSums(trinucfreq_table[,trinuc_list])
L_matrix_global = L_matrix_ref
L_matrix_global[names(trin_freqs),] = L_matrix_global[names(trin_freqs),] + array(rep(trin_freqs,4), dim=c(length(trin_freqs),4))

trin_freqs = colSums(trinucfreq_table_neutral[,trinuc_list])
L_matrix_global_neutral = L_matrix_ref
L_matrix_global_neutral[names(trin_freqs),] = L_matrix_global_neutral[names(trin_freqs),] + array(rep(trin_freqs,4), dim=c(length(trin_freqs),4))



## 2. Mapping the mutations to the regions and to neutral bins

# Subfunction intersecting the mutations with the regions of interest
cat("Defining functions...\n")

map_mutations = function(mutations, intervals) {

    rangesREGS = GRanges(intervals[,2], IRanges(as.numeric(intervals[,3]), as.numeric(intervals[,4])))
    mut_start = as.numeric(as.vector(mutations$pos))
    mut_end = mut_start + nchar(as.vector(mutations$ref)) - 1
    rangesMUTS = GRanges(as.vector(mutations$chr), IRanges(mut_start, mut_end))
    olmatrix = as.matrix(findOverlaps(rangesMUTS, rangesREGS, type="any", select="all"))
    return(olmatrix)
}

# Subfunction that subsamples mutations to a maximum of nmax mutations per sample per region of interest

cap_mutations = function(m = sample_and_region, nmax = max_num_muts_perRegion_perSample) {

    nmut = dim(m)[1]
    mapped_muts = which(m[,2]!="0")
    m = m[mapped_muts,]

    rows = paste(m[,1],m[,2],sep=",")
    freqs = table(rows)
    freqs = freqs[freqs>nmax]
    duplrows = strsplit(names(freqs),split=",")
    rmrows = array(0,dim(m)[1])
    
    if (length(freqs)>0) {
        for (j in 1:length(freqs)) {
            vals = duplrows[[j]]
            if (vals[2]!="0") { # Only for mapped mutations we apply the capping
                pos = which(m[,1]==vals[1] & m[,2]==vals[2])
                rmrows[sample(pos, length(pos)-nmax)] = 1
            }
        }
    }
    rmrows_allmuts = array(0,nmut)
    rmrows_allmuts[mapped_muts[rmrows==1]] = 1
    return(rmrows_allmuts)
}



# a. Mapping the mutations to the regions of interest
cat("Mapping mutations...\n")

olm = map_mutations(mutations, target_regions)
mutations$target_region = 0
m1 = mutations[olm[,1],] # Duplicating subs if they hit more than one region
m1$target_region = target_regions[olm[,2],1] # Annotating hit region
m2 = mutations[-unique(olm[,1]),] # Mutations not mapping to any element
mutations = rbind(m1,m2)

# b. Mapping the remaining mutations to the neutral regions
cat("Mapping mutations (2)...\n")

olm = map_mutations(mutations, neutral_regions)
mutations$neutral_region = 0
m1 = mutations[olm[,1],] # Duplicating subs if they hit more than one region
m1$neutral_region = neutral_regions[olm[,2],1] # Annotating hit region
m2 = mutations[-unique(olm[,1]),] # Mutations not mapping to any element
mutations = rbind(m1,m2)

# c. Masking out duplicate (non-unique) indels (if desired)

if (unique_indelsites_FLAG==1) { # We mask out duplicate indels
    mutations$target_region[mutations$unique_indels==F] = 0
    mutations$neutral_region[mutations$unique_indels==F] = 0
}


# c. Subsampling mutations to a maximum of nmax mutations per sample per region of interest

sample_and_region = cbind(as.vector(mutations$sampleID), as.character(mutations$target_region))
maskmuts = cap_mutations(sample_and_region, max_num_muts_perRegion_perSample)
mutations$masked_muts_bynmax = maskmuts
mutations$target_region[maskmuts==1] = 0


### 3. Calculating the global trinucleotide rates and the local density of mutations

mutations$trinuc_ref[!(as.vector(mutations$trinuc_ref) %in% trinuc_list)] = NA

n_matrix_global = L_matrix_ref
n_matrix_global_neutral = L_matrix_ref

subsin = (!is.na(mutations$trinuc_ref)) & (mutations$target_region!=0)
aux = cbind(as.vector(mutations$trinuc_ref[subsin]), as.vector(mutations$mut[subsin]))
for (j in 1:dim(aux)[1]) {
    n_matrix_global[aux[j,1],aux[j,2]] = n_matrix_global[aux[j,1],aux[j,2]] + 1
}
numindels_target = sum((is.na(mutations$trinuc_ref)) & (mutations$target_region!=0))

subsin = (!is.na(mutations$trinuc_ref)) & (mutations$neutral_region!=0) & (mutations$target_region==0)
aux = cbind(as.vector(mutations$trinuc_ref[subsin]), as.vector(mutations$mut[subsin]))
for (j in 1:dim(aux)[1]) {
    n_matrix_global_neutral[aux[j,1],aux[j,2]] = n_matrix_global_neutral[aux[j,1],aux[j,2]] + 1
}
numindels_neutral = sum((is.na(mutations$trinuc_ref)) & (mutations$neutral_region!=0))

rate_names = array(NA, dim=c(64,4))
for (j in 1:dim(L_matrix_ref)[1]) {
    for (h in 1:dim(L_matrix_ref)[2]) {
        rate_names[j,h] = sprintf("%s>%s%s%s", trinuc_list[j], substr(trinuc_list[j],1,1), nt[h], substr(trinuc_list[j],3,3))
    }
}


## a. Rates
rates_target = c(n_matrix_global/L_matrix_global); names(rates_target) = c(rate_names)
rates_neutral = c(n_matrix_global_neutral/L_matrix_global_neutral); names(rates_neutral) = c(rate_names)

indelsrate_target = numindels_target/sum(L_matrix_global,na.rm=T)*3
indelsrate_neutral = numindels_neutral/sum(L_matrix_global_neutral,na.rm=T)*3


## b. Local density of mutations

# Expected number of subs and indels

targetregions_df = data.frame(region=as.vector(trinucfreq_table[,1]))
aux = strsplit(as.vector(trinucfreq_table[,2]), split=":")
targetregions_df$chr = sapply(aux, function(x) x[1])
#aux2 = sapply(aux, function(x) min(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))

targetregions_df$start = sapply(aux, function(x) min(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))
targetregions_df$end = sapply(aux, function(x) max(suppressWarnings(as.numeric(unlist(strsplit(x[-1], split="-")))), na.rm=T))

#targetregions_df$start = sapply(aux2, function(x) min(x,na.rm=T))
#targetregions_df$end = sapply(aux2, function(x) max(x,na.rm=T))
# Shouldn't the previous two lines be replaced by the following two instead? At least on two occasions I have needed to run it like this...
#targetregions_df$start = apply(aux2,2,min)
#targetregions_df$end = apply(aux2,2,max)


tf = as.matrix(trinucfreq_table[,trinuc_list])
targetregions_df$exp_subs = apply(tf, 1, function(x) sum(rep(x,4)*rates_target, na.rm=T) )
targetregions_df$exp_indels = apply(tf, 1, function(x) sum(x)*indelsrate_target )

neutralregions_df = data.frame(region=as.vector(trinucfreq_table_neutral[,1:3]))
tf = as.matrix(trinucfreq_table_neutral[,trinuc_list])
neutralregions_df$exp_subs = apply(tf, 1, function(x) sum(rep(x,4)*rates_neutral, na.rm=T) )
neutralregions_df$exp_indels = apply(tf, 1, function(x) sum(x)*indelsrate_neutral )

# Observed number of subs and indels

targetregions_df$obs_subs = 0
targetregions_df$obs_indels = 0
indel_pos = is.na(mutations$trinuc_ref)
numsubs = table( mutations$target_region[mutations$target_region > 0 & !indel_pos] )
targetregions_df$obs_subs[as.numeric(names(numsubs))] = numsubs
numinds = table( mutations$target_region[mutations$target_region > 0 & indel_pos] )
targetregions_df$obs_indels[as.numeric(names(numinds))] = numinds

neutralregions_df$obs_subs = 0
neutralregions_df$obs_indels = 0
indel_pos = is.na(mutations$trinuc_ref)
numsubs = table( mutations$neutral_region[mutations$neutral_region > 0 & !indel_pos] )
neutralregions_df$obs_subs[as.numeric(names(numsubs))] = numsubs
numinds = table( mutations$neutral_region[mutations$neutral_region > 0 & indel_pos] )
neutralregions_df$obs_indels[as.numeric(names(numinds))] = numinds


# Estimating the neighbourhood "t" for every target region
# We choose as neighbouring neutral regions of a given target region all those that are 
# within "neighbourhood_localrate" distance of the target region. We consider any overlap
# between the segments as valid, which means that the neighbourhood of a target region
# will always be contained within the interval "neighbourhood_localrate + 100kb" around
# the target region (for a 100kb binning of the genome in the neutral reference)
# For example, if neighbourhood_localrate=1e5, only regions less or equal to 200kb away 
# from the ends of the target region are considered in the calculation of the local rate.

neighbourhood_localrate = ceiling(0.001/dim(mutations)[1]*3e9)*100000

rangesTARGET = GRanges(target_regions[,2], IRanges(as.numeric(target_regions[,3])-neighbourhood_localrate, as.numeric(target_regions[,4])+neighbourhood_localrate))
rangesNEUTRAL = GRanges(neutralregions_df[,1], IRanges(as.numeric(neutralregions_df[,2]), as.numeric(neutralregions_df[,3])))
ol = findOverlaps(rangesTARGET, rangesNEUTRAL, type="any", select="all")
olmatrix = as.matrix(ol)
neighbours = unique(cbind(target_regions[olmatrix[,1],1],olmatrix[,2]))

targetregions_df$local_t = NA
targetregions_df$local_t_indels = NA

for (j in 1:dim(targetregions_df)[1]) {
    neutral_bins = neighbours[neighbours[,1]==j,2]
    targetregions_df$local_t[j] = sum(neutralregions_df$obs_subs[neutral_bins])/sum(neutralregions_df$exp_subs[neutral_bins])
    targetregions_df$local_t_indels[j] = sum(neutralregions_df$obs_indels[neutral_bins])/sum(neutralregions_df$exp_indels[neutral_bins])
    if ((j/1000)==round(j/1000)) { print(j/dim(targetregions_df)[1]) }
}



### 4. Negative binomial regression with local density of mutations and covariates

# Masking out regions with ZERO expected values

targetregions_df$exp_subs[targetregions_df$exp_subs==0] = NA
targetregions_df$exp_indels[targetregions_df$exp_indels==0] = NA

## 4a. MODEL 1: No use of local mutation rates or covariates

# Negative binomial regression
model_subs = glm.nb(formula = targetregions_df$obs_subs ~ offset(log(targetregions_df$exp_subs)) -1 )
nb_size_subs = model_subs$theta
if (numindels_target>0) {
    model_indels = glm.nb(formula = targetregions_df$obs_indels ~ offset(log(targetregions_df$exp_indels)) -1 )
    nb_size_indels = model_indels$theta
}

# P-values (Neg Binom)

targetregions_df$pval_subs = NA
targetregions_df$pval_indels = NA
targetregions_df$pval_both = NA

for (j in 1:dim(targetregions_df)[1]) {
    targetregions_df$pval_subs[j] = pnbinom(q=targetregions_df$obs_subs[j]-0.1, mu=targetregions_df$exp_subs[j], size=nb_size_subs, lower.tail=F)
    
    if (numindels_target>0) {
        targetregions_df$pval_indels[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1, mu=targetregions_df$exp_indels[j], size=nb_size_indels, lower.tail=F)
        # Fisher combined p-value
        p_vec = c(targetregions_df$pval_subs[j], targetregions_df$pval_indels[j])
        targetregions_df$pval_both[j] = 1-pchisq(-2*sum(log(p_vec)),length(p_vec)*2)
    } else {
        # We use only subs
        targetregions_df$pval_both[j] = targetregions_df$pval_subs[j]
    }
    if (round(j/1000)==(j/1000)) { print(j/dim(targetregions_df)[1]) }
}
targetregions_df$qval_subs = p.adjust(targetregions_df$pval_subs, method="BH") # Adjusted q-value
targetregions_df$qval_indels = p.adjust(targetregions_df$pval_indels, method="BH") # Adjusted q-value
targetregions_df$qval_both = p.adjust(targetregions_df$pval_both, method="BH") # Adjusted q-value


## 4b. MODEL 2: Using local mutation rates and covariates

# Excluding elements overlapping "driver" genomic regions from the negbin background fit
gr_elements = GRanges(targetregions_df$chr, IRanges(targetregions_df$start, targetregions_df$end))
olmatrix = as.matrix(findOverlaps(gr_elements, gr_drivers, type="any", select="all"))
exclude_forfit = (1:nrow(targetregions_df)) %in% unique(olmatrix[,1])
targetregions_df$exclude_forfit = exclude_forfit

# Negative binomial regression
model_subs = glm.nb(formula = obs_subs ~ offset(log(exp_subs)) + local_t , data=targetregions_df[!exclude_forfit,])
nb_size_subs_cv = model_subs$theta
targetregions_df$cv_predicted_subs = exp(predict(model_subs, newdata=targetregions_df))

if (numindels_target>0) {
    model_indels = glm.nb(formula = obs_indels ~ offset(log(exp_indels)) + local_t_indels , data=targetregions_df[!exclude_forfit,]) ## FOR INDELS THE LOCAL RATE IS A COVARIATE
    nb_size_indels_cv = model_indels$theta
    targetregions_df$cv_predicted_indels = exp(predict(model_indels, newdata=targetregions_df))
}

# P-values (Neg Binom)

targetregions_df$pval_subs_CV = NA
targetregions_df$pval_indels_CV = NA
targetregions_df$pval_both_CV = NA

for (j in 1:dim(targetregions_df)[1]) {
    targetregions_df$pval_subs_CV[j] = pnbinom(q=targetregions_df$obs_subs[j]-0.1, mu=targetregions_df$cv_predicted_subs[j], size=nb_size_subs_cv, lower.tail=F)
    #targetregions_df$pval_indels_CV[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1, mu=targetregions_df$cv_predicted_indels[j], size=nb_size_indels_cv, lower.tail=F)
    
    if (numindels_target>0) {
        targetregions_df$pval_indels_CV[j] = pnbinom(q=targetregions_df$obs_indels[j]-0.1, mu=targetregions_df$cv_predicted_indels[j], size=nb_size_indels_cv, lower.tail=F)
        # Fisher combined p-value
        p_vec = c(targetregions_df$pval_subs_CV[j], targetregions_df$pval_indels_CV[j])
        targetregions_df$pval_both_CV[j] = 1-pchisq(-2*sum(log(p_vec)),length(p_vec)*2)
    } else {
        # We use only subs
        targetregions_df$pval_both[j] = targetregions_df$pval_subs[j]
    }
    if (round(j/1000)==(j/1000)) { print(j/dim(targetregions_df)[1]) }
}
targetregions_df$qval_subs_CV = p.adjust(targetregions_df$pval_subs_CV, method="BH") # Adjusted q-value
targetregions_df$qval_indels_CV = p.adjust(targetregions_df$pval_indels_CV, method="BH") # Adjusted q-value
targetregions_df$qval_both_CV = p.adjust(targetregions_df$pval_both_CV, method="BH") # Adjusted q-value


# Tiered FDR correction

targetregions_df$qval_both_CV_tiered = NA
inds = targetregions_df$exclude_forfit
targetregions_df$qval_both_CV_tiered[inds] = p.adjust(targetregions_df$pval_both_CV[inds], method="BH") # Adjusted q-value
targetregions_df$qval_both_CV_tiered[!inds] = p.adjust(targetregions_df$pval_both_CV[!inds], method="BH") # Adjusted q-value
targetregions_df = targetregions_df[order(targetregions_df$qval_both_CV_tiered),]

write.table(targetregions_df, file=paste(output_prefix,"-Selection_output.txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)


# Plot of the underlying gamma distributions of the 4 models (subs and indels with and without covariates)

pdf(paste(output_prefix,"-Underlying_gamma_distributions_landscape.pdf",sep=""),height=4,width=4)
xvec = seq(0,4,by=0.0001)
plot(xvec,0.0001*dgamma(x=xvec,shape=nb_size_subs_cv,rate=nb_size_subs_cv),type="l",lty=1,xlab="Relative mutation rate",ylab="",col="cadetblue",las=1,main="PDFs of the underlying Gamma distributions",cex.axis=0.8,cex.main=0.85)
lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_subs,rate=nb_size_subs),lty=2,col="cadetblue")
#lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_indels,rate=nb_size_subs),lty=2,col="chocolate")
#lines(xvec,0.0001*dgamma(x=xvec,shape=nb_size_indels_cv,rate=nb_size_subs_cv),lty=1,col="chocolate")
#legend("topright", lty=c(2,1,2,1), col=c("cadetblue","cadetblue","chocolate","chocolate"), title="Size parameters", legend=c(sprintf("subs=%0.3g",nb_size_subs),sprintf("subs_cv=%0.3g",nb_size_subs_cv),sprintf("indels=%0.3g",nb_size_indels),sprintf("indels_cv=%0.3g",nb_size_indels_cv)), box.col=NA, cex=0.8)
legend("topright", lty=c(2,1), col=c("cadetblue","cadetblue"), title="Size parameters", legend=c(sprintf("subs=%0.3g",nb_size_subs),sprintf("subs_cv=%0.3g",nb_size_subs_cv)), box.col=NA, cex=0.8)
dev.off()


## Confidence intervals for the obs/exp ratios in NBR (23.12.2016)

calculate_CI95_flag = 1

# Subfunction to calculate CI95% for the obs/exp ratios of subs and indels in NBR

ci95nbr = function(n_obs,n_exp,nb_size) {

    wmax = 100000; iter = 6; grid_size = 10; cutoff = qchisq(p=0.95,df=1) # Default params
    wmle = n_obs/n_exp # MLE for w
    ml = dnbinom(x=n_obs, mu=n_exp*wmle, size=nb_size, log=T) # LogLik under MLE

    if (!is.na(n_exp)) {
      if (wmle<wmax) {

        # 1. Iterative search of lower bound CI95%
        if (wmle>0) {
            search_range = c(1e-9, wmle)
            for (it in 1:iter) {
                wvec = seq(search_range[1], search_range[2], length.out=grid_size)
                ll = dnbinom(x=n_obs, mu=n_exp*wvec, size=nb_size, log=T)
                lr = 2*(ml-ll) > cutoff
                ind = max(which(wvec<=wmle & lr))
                search_range = c(wvec[ind], wvec[ind+1])
            }
            w_low = wvec[ind]
        } else {
            w_low = 0
        }
        
        # 2. Iterative search of higher bound CI95%
        search_range = c(wmle, wmax)
        llhighbound = dnbinom(x=n_obs, mu=n_exp*wmax, size=nb_size, log=T)
        outofboundaries = !(2*(ml-llhighbound) > cutoff)
        if (!outofboundaries) {
            for (it in 1:iter) {
                wvec = seq(search_range[1], search_range[2],length.out=grid_size)
                ll = dnbinom(x=n_obs, mu=n_exp*wvec, size=nb_size, log=T)
                lr = 2*(ml-ll) > cutoff
                ind = min(which(wvec>=wmle & lr))
                search_range = c(wvec[ind-1], wvec[ind])
            }
            w_high = wvec[ind]
        } else {
            w_high = wmax
        }
      } else {
        wmle = w_low = w_high = wmax # Out of bounds
      }
    } else {
      wmle = w_low = w_high = NA # invalid
    }

    return(c(wmle,w_low,w_high))
}

if (calculate_CI95_flag == 1) {
    ci95_subs = t(apply(as.matrix(targetregions_df[,c("obs_subs","cv_predicted_subs")]), 1, function(x) tryCatch(ci95nbr(x[1], x[2], nb_size_indels_cv), error=function(err) rep(NA,3))))
    targetregions_df$obsexp_subs_mle = ci95_subs[,1]
    targetregions_df$obsexp_subs_low = ci95_subs[,2]
    targetregions_df$obsexp_subs_high = ci95_subs[,3]
    
    #ci95_indels = t(apply(as.matrix(targetregions_df[,c("obs_indels","cv_predicted_indels")]), 1, function(x) tryCatch(ci95nbr(x[1], x[2], nb_size_indels_cv), error=function(err) rep(NA,3))))
    #targetregions_df$obsexp_indels_mle = ci95_indels[,1]
    #targetregions_df$obsexp_indels_low = ci95_indels[,2]
    #targetregions_df$obsexp_indels_high = ci95_indels[,3]

    write.table(targetregions_df, file=paste(output_prefix,"-Selection_output_landscape.txt",sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)
    write.table(mutations, file=paste(output_prefix,"-Intermediate_mutations.txt", sep=""), sep = "\t", row.names = FALSE, col.names = TRUE, append=FALSE, quote=FALSE)
}

