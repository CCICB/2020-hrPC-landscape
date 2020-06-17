###############################################################
## General functions used on the process of NBR calculations ##
###############################################################

#################################
## Collapse contiguous regions ##
#################################

Collapse_contiguous_landscape <- function(mutations){
  mutations <- as.data.frame(mutations)
  head(mutations)
  
  ## remove mutation column burder
  mutations <- mutations[,-6]
  mutations = mutations[order(mutations$sampleID,mutations$chr,mutations$pos),]
  ## finds all the positions where the difference is 1
  ind = which(diff(mutations$pos)==1)
  contiguous <- mutations[unique(sort(c(ind,ind+1))),]
  dim(contiguous)
  
  ## split each contig into a list
  dat <- unique(sort(c(ind,ind+1)))
  indsp <- split(dat, cumsum(c(1, diff(dat) != 1)))
  
  
  colp <- list()
  for(i in 1:length(indsp)){
    mut <- mutations[indsp[[i]],]
    if(length(unique(mut$sampleID))!= 1){
      print("error not the same sample")
    }else{
      out <- data.frame(sampleID = unique(mut$sampleID),chromosome=unique(mut$chromosome), position=mut$position[1], ref=paste(mut$ref, collapse =""), alt = paste(mut$alt, collapse =""),  stringsAsFactors = F)
      colp[[i]] <- out
    }
  }
  
  coldf <- do.call(rbind, colp)
  notcontig <- mutations[-unique(sort(c(ind,ind+1))),]
  
  all_collapsed <- rbind(coldf, notcontig)
  
  allo <- all_collapsed[order(all_collapsed$sampleID,all_collapsed$chr,all_collapsed$pos),]
  return(allo)
}


#################################################
## Manually calculate number of BPs subtracted ##
#################################################

manual_subtract_calculations <- function(infile_WAO, infile_subtract, subgenes, infile_sample_mutations.ls){
  
  ## intersected transcript regions
  bed_intersect <- fread(infile_WAO)
  bed_intersect_df <- as.data.frame(bed_intersect)
  bed_intersect_df <- bed_intersect_df[c(1:4,6:10)]
  colnames(bed_intersect_df) <- c("Input_chr","Input_start","Input_end","Input_name","Input_strand","Blacklist_chr","Blacklist_start","Blacklist_end","bp_ol")
  
  ### subtracted transcript regions
  bed_subtract <- fread(infile_subtract)
  bed_subtract_df <- as.data.frame(bed_subtract)
  bed_subtract_df <- bed_subtract_df[c(1:4,6)]
  colnames(bed_subtract_df) <- c("Subtract_chr","Subtracted_start","Subtracted_end","Subtracted_name","Subtracted_strand")
  
  
  bed_intersect_candid <- bed_intersect_df[grep(paste(subgenes, collapse="|"),bed_intersect_df$Input_name),]
  bed_subtract_candid <- bed_subtract_df[grep(paste(subgenes, collapse="|"),bed_subtract_df$Subtracted_name),]
  
  ## merge data
  m <- merge(bed_intersect_candid, bed_subtract_candid, by.x="Input_name", by.y="Subtracted_name", all=T)
  
  ## if region entirely contained in tier3.bed will not be included in output, hence add back in to make it easier downstream 
  sna <- which(is.na(m$Subtract_chr))
  m$Subtract_chr[sna] <- m$Input_chr[sna]
  m$Subtracted_start[sna] <- m$Input_start[sna]
  m$Subtracted_end[sna] <- m$Input_end[sna]
  m$Subtracted_strand[sna] <- m$Input_strand[sna]
  
  sp <- split(m, f=as.vector(m$Input_name))
  
  mutsfiles <- list()
  for(file in 1:length(infile_sample_mutations.ls)){
    hyper <- fread(infile_sample_mutations.ls[[file]])
    hyperdf <- as.data.frame(hyper)
    mutsfiles[[file]] <- hyperdf
  }
  
  df.out <- list()
  ## in this loop collapse discontinous transcripts
  for(i in 1:length(sp)){
    ## can have the scenario of mutiple tier3, but end up with one subtracted region - hence make unique
    subr <- unique(sp[[i]][,10:13])

    ## count number of mutations in each region of the subtracted transcript
    mutsfiles.count <- list()
    for (file in 1:length(mutsfiles)){
        sn_hyper <- 0
        uniqmuts <- 0
        for(j in 1: nrow(subr)){
          ## get index of samples in region
          w2 <- which(mutsfiles[[file]]$chromosome  == subr[j,"Subtract_chr"] & mutsfiles[[file]]$position > subr[j,"Subtracted_start"] & mutsfiles[[file]]$position <= subr[j,"Subtracted_end"])
          
          ## info about unqique samples
          samp2 <- unique(mutsfiles[[file]][w2,"sampleID"])
          sn_hyper <- sn_hyper + length(samp2)
          mutsfiles.count[[file]] <- sn_hyper
          # 
          
        }
    }
    ## format so names are correct
    names(mutsfiles.count) <- lapply(infile_sample_mutations.ls, function(x) gsub(".txt", "", gsub(".*/", "", x)))
    numtsdf <- do.call(cbind, mutsfiles.count)  
      
    ## get the size of all subtracted regions combined 
    Subtracted_size <- 0
    for(j in 1: nrow(subr)){  
      sz <- subr[j, "Subtracted_end"] - subr[j, "Subtracted_start"]
      Subtracted_size  <- Subtracted_size + sz
    }
    
    subra <- as.data.frame(sp[[i]])
    Input_size <- unique(subra[,"Input_end"] - subra[,"Input_start"])
    
    ## distinguish between 1 tier3 bed or multiple tier3 bed
    if(nrow(unique(subra[,6:8])) == 1){
      BP_OL <- unique(subra[,"bp_ol"])
    }else{
      BP_OL <- sum(subra[,"bp_ol"])
    }
    
    ## set subtracted to 0 if its the same as input
    if(Subtracted_size == Input_size){
      remaining_BP = 0
      Subtracted_size = 0 
      OL_plus_Subtracted = Input_size
    }else{
      remaining_BP = Subtracted_size
      OL_plus_Subtracted <- Subtracted_size + BP_OL
    }
    ## collapse tier3 regions so there is only one line
    if(unique(subra$Blacklist_start) == -1){
      Blacklist_region = NA
    }else{
      a <- apply(subra,1, function(x) paste(as.vector(x["Blacklist_start"]), as.vector(x["Blacklist_end"]), sep="-"))
      Blacklist_region <- paste(a,collapse=",")
    }
    
    df <- data.frame(Input_name= unique(sp[[i]][,"Input_name"]),Input_chr= unique(sp[[i]][,"Input_chr"]), Input_start= unique(sp[[i]][,"Input_start"]), Input_end= unique(sp[[i]][,"Input_end"]), Input_strand = unique(sp[[i]][,"Input_strand"]), Blacklist_region = Blacklist_region, remaining_BP = remaining_BP, BP_OL = BP_OL, Input_size=Input_size, OL_plus_Remaining = OL_plus_Subtracted,  nmuts = numtsdf)
    df.out[[i]] <- df
  }
  
  df.all <- do.call(rbind, df.out)
  
  Blacklist_regionbin <- rep(NA, nrow(df.all))
  Blacklist_regionbin[which(is.na(df.all$Blacklist_region))] <- "N"
  Blacklist_regionbin[which(!is.na(df.all$Blacklist_region))] <- "Y"
  df.fin <- cbind(df.all[,1:5],Blacklist_regionbin,df.all[,6:ncol(df.all)])
  df.fin[, c(1,10:11)]
  df.fin
  return(df.fin)
}

########################################################
## Extracgt patient name and mutations within regions ##
########################################################

get_patient_info <- function(infile_mutinter, infile_minhyper_results,nreg=NULL,glist=NULL,pvaltype){
  
  mutinter <- fread(infile_mutinter)
  mutinterdf <- as.data.frame(mutinter)
  colnames(mutinterdf)[c(2:3,5)] <- c("chromosome","position", "alt")
  mutinterdf$mutID <- paste(mutinterdf$chromosome,":", mutinterdf$position, "-", mutinterdf$ref, "|", mutinterdf$alt, sep="")
  
  minhyper_results <- fread(infile_minhyper_results)
  minhyper_results_df <- as.data.frame(minhyper_results)
  
  if(!is.null(nreg)){
    #minhyper_top10 <- minhyper_results_dfo[1:10,c(1:8,11:13)]
    minhyper_results_dfo <- minhyper_results_df[order(minhyper_results_df[,grep(pvaltype,colnames(minhyper_results_df))], decreasing=F),]
    minhyper_top10 <- minhyper_results_dfo[1:nreg,]
  }else if(!is.null(glist)){
    minhyper_results_df_glist <- minhyper_results_df[grep(paste(glist, collapse="|"), minhyper_results_df$region),]
    minhyper_top10 <- minhyper_results_df_glist[order(minhyper_results_df_glist$pval_both, decreasing=F),]
  }
  
  ## which patients have mutations in the top significant regions?? 
  minhyper_top10_plus_samples.ls <- list()
  for(i in 1:nrow(minhyper_top10)){
    ## get patient and mutation info
    w <- which(mutinterdf$chromosome == minhyper_top10[i,"chr"] & mutinterdf$position > minhyper_top10[i,"start"] & mutinterdf$position <= minhyper_top10[i,"end"])
    
    ## if zero patient info 
    if(length(w) == 0){
      sub.df <- data.frame(mutID=NA, nsamps=NA, sample_names=NA)
      out <- cbind(i, minhyper_top10[i,c("region", "chr", "start", "end")], sub.df, minhyper_top10[i,c("exp_subs","exp_indels","obs_subs","obs_indels", "pval_subs_CV","pval_indels_CV","pval_both_CV","qval_subs_CV","qval_indels_CV","qval_both_CV")])
      colnames(out)[1] <- "Rank"
      minhyper_top10_plus_samples.ls[[i]] <- out
      ## if there is sample information  
    }else{
      samples_mutations <- mutinterdf[w,]
      ## split on each mutation
      smds <- split(samples_mutations, f=samples_mutations$mutID)
      
      ## loop through each mutation to get 
      sub.ls <- list()
      print(i)
      for(j in 1:length(smds)){
        ## remove those that don't have target region as 0 - means they're filtered
        wk <- which(smds[[j]]$target_region != 0)
        if(length(wk) != 0 ){
          rel <- smds[[j]][wk,]
          usamps <- unique(rel$sampleID)
          nsamps <- length(usamps)
          sample_names <- paste(usamps, collapse=";")
          mutID <- unique(rel$mutID)
          subdf <- cbind(mutID,nsamps,sample_names)
          sub.ls[[j]] <- subdf
        }
      }
      
      sub.df <- do.call(rbind, sub.ls)
      out <- cbind(i, minhyper_top10[i,c("region", "chr", "start", "end")], sub.df, minhyper_top10[i,c("exp_subs","exp_indels","obs_subs","obs_indels", "pval_subs_CV","pval_indels_CV","pval_both_CV","qval_subs_CV","qval_indels_CV","qval_both_CV")])
      colnames(out)[1] <- "Rank"
      minhyper_top10_plus_samples.ls[[i]] <- out
    }
  }
  
  ## combine results - end up with some duplication
  outdf <- do.call(rbind, minhyper_top10_plus_samples.ls)
}


