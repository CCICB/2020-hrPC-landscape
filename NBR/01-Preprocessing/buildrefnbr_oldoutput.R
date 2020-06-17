#' buildrefnbr
#' 
#' Function to generate reference files for NBR from a reference genome and interval or bed-like files.
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Rheinbay, Muhlig Nielsen, Abascal, et al. (2020) Analyses of non-coding somatic drivers in 2,658 cancer whole genomes. Nature. 578, 102–111.
#' 
#' @param bedfile Path to a bed file with the regions to be analysed. Files must be 4-column or 12-column standard bed files with names. 
#' @param genomefile Path to the indexed reference genome file.
#' @param onlychrs Vector of valid chromosome names (default: all chromosomes will be included)
#' 
#' @export

library(Rsamtools)
library(GenomicRanges)

buildrefnbr = function(bedfile, genomeFile, onlychrs = NULL) {
    
    ## 1. Processing the input bedfile file
    message("[1/3] Processing the input bedfile file...")
    
    # Loading and processing the input bed file
    bed = read.table(bedfile, header=0, sep="\t", stringsAsFactors=F)
    
    if (ncol(bed)==4) { # bed4
        
        bed[,3] = bed[,3] + 1 # Adding a 1bp flank (the start field in a bed file is already 0-based)
        colnames(bed) = c("chr","start","end","name")
        if (!is.null(onlychrs)) {
            bed = bed[bed$chr %in% onlychrs, ]
        }
        if (nrow(bed)==0) {
            stop("Empty bed file: if you are using the onlychrs argument, ensure that chromosome names match those in the bed file.")
        }
        regiondf = split(bed, f=bed[,4]) # List the regions for each element
        for (j in 1:length(regiondf)) {
            regiondf[[j]]$bin = j
        }

    } else if (ncol(bed)==12) { # bed12
        
        if (!is.null(onlychrs)) {
            bed = bed[bed[,1] %in% onlychrs, ] # Restricting the analysis to chrs in the onlychr argument
        }
        if (nrow(bed)==0) {
            stop("Empty bed file: if you are using the onlychrs argument, ensure that chromosome names match those in the bed file.")
        }
        regiondf = vector(mode="list", nrow(bed))
        for (j in 1:nrow(bed)) {
            start = bed[j,2] + 1
            segment_sizes = as.numeric(strsplit(bed[j,11],split=",")[[1]])
            segment_starts = as.numeric(strsplit(bed[j,12],split=",")[[1]])
            regiondf[[j]] = data.frame(chr=bed[j,1], start=start+segment_starts-1, end=start+segment_starts+segment_sizes, name=bed[j,4], bin=j, stringsAsFactors=F) # Regions with 1bp flanks
        }
        names(regiondf) = sapply(regiondf, function(x) x[1,4])
        
    } else { # Invalid: bed file must have 4 or 12 columns
        stop("Invalid input bed file: please use a bed file with 4 columns, with the 4th column being the element name, or a bed12 file with 12 columns.")
    }
    
    ## 2. Calculating the trinucleotide frequencies for each element
    message("[2/3] Calculating the trinucleotide frequencies for each element...")
    
    regionstrin = array(0, dim=c(length(regiondf),64)) # Trinucleotide composition
    regionslist = data.frame(region=names(regiondf), segments=NA, stringsAsFactors=F) # Name and coordinates of each element
    regionspos = data.table::rbindlist(regiondf) # All regions as a single data.frame
        
    for (j in 1:length(regiondf)) {
        
        gr = reduce(GenomicRanges::GRanges(regiondf[[j]][,1], IRanges::IRanges(regiondf[[j]][,2], regiondf[[j]][,3]))) # Sequences with 1bp flanks
        seqs = scanFa(genomeFile, gr)
        regionstrin[j, ] = colSums(Biostrings::trinucleotideFrequency(seqs)) # Trinucleotide frequencies
        regionslist[j,2] = paste(regiondf[[j]]$chr, paste(regiondf[[j]]$start+1,regiondf[[j]]$end-1,sep="-"), sep=":", collapse=",") # String with all segments (1-based)
        
        if ((j/1000)==round(j/1000)) { print(round(j/nrow(bed),2)) }
    }
    
    nt = c("A","C","G","T")
    base1 = rep(nt,each=16,times=1)
    base2 = rep(nt,each=4,times=4)
    base3 = rep(nt,each=1,times=16)
    trinuc_list = paste(base1,base2,base3, sep="")
    
    colnames(regionstrin) = trinuc_list
    out1 = cbind(regionslist, regionstrin)
    out2 = regionspos[,c("bin","chr","start","end")]
    out2$start = out2$start + 1 # Removing the trinucleotide flank
    out2$end = out2$end - 1 # Removing the trinucleotide flank
    colnames(out2) = c("bin","seqnames","start","end")
    
    # Writing output files
    write.table(out1, file=sprintf("%s.txt",bedfile), col.names=T, row.names=F, sep="\t", quote=F)
    write.table(out2, file=sprintf("%s.regions",bedfile), col.names=T, row.names=F, sep="\t", quote=F)
}
