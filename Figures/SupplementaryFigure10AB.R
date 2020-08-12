###########################################################
### Supplementary Figure 10 code - PLatform concordance ###
###########################################################

library(stringr)
library(ggsci)
library(dplyr)
library(ggplot2)
library(reshape)

##############
## New Data ##
##############

## file location: https://docs.google.com/spreadsheets/d/1eWpk3iPvRIZwWDUlQ7fwo0URox8y94Zvszv7IWFMJ_w/edit#gid=900627208
mutcounts <- read.csv("Data/landscape_platform_mutcounts - 15June2020_data.csv", sep=",", header=T)
str(mutcounts)

## looking at the data
str(mutcounts)
table(mutcounts$platforms)
table(mutcounts$muttype)

## add new type collumn
mutcounts$type_mod <- mutcounts$muttype

#####################
## Modify plotform ##
#####################

## #Remove the P from everything, W, R, WR
mutcounts$platform_mod <- sub("P", "", mutcounts$platforms)

## add germline
mutcounts$platform_mod[which(mutcounts$muttype == "germline")] <- "G"

#########################################
## TMB patients - remove and summarize ##
#########################################

tmb_pat <- c("LKCGP-P000202-248737-01-01-01-R1", "LKCGP-P001202-255259-01-01-01-R1", "LKCGP-P002504-271859-01-04-01-R1", "LKCGP-P004104-280143-01-01-01-D1", "LKCGP-T005601-241252-01-01-01-D1", "LKCP-T002401-031582-01-01-01-D1", "LKCP-T003001-168602-01-01-01-D1")

tmbdf <- as.data.frame(tmb_pat)
tmbdf$platforms <- "W"
tmbdf$muttype <- "tmb"
tmbdf$mutcount <- 1
tmbdf$type_mod <- "tmb"
tmbdf$platform_mod <- "W"
colnames(tmbdf)[1] <- "sample"

## remove tmb patients
mutcounts_rm <- mutcounts[-which(mutcounts$sample %in% tmbdf$sample),]

## add tmb patients with different platform
mutcounts <- rbind(mutcounts_rm, tmbdf)
str(mutcounts)

#####################
## add methylation ##
#####################

final3 <- read.table("Data/Methylation_Match.txt", sep="\t", header=T, stringsAsFactors = F)

str(final3)
table(final3$diag_vec)

meth_sample <- final3[which(final3$diag_vec %in% c("Same", "Changed")),"sample_id"]

names_meth <- as.data.frame(meth_sample)
names_meth$platforms <- "M"
names_meth$muttype <- "meth"
names_meth$mutcount <- 1
names_meth$type_mod <- "meth"
names_meth$platform_mod <- "M"
colnames(names_meth)[1] <- "sample"

head(names_meth)

## add tmb patients with different platform
mutcounts <- rbind(mutcounts, names_meth)
str(mutcounts)

table(mutcounts$platform_mod)
table(mutcounts$type_mod)

## change column name so I don't have to change code below
colnames(mutcounts)[1] <- "sample_id"
head(mutcounts)

## checking - good
table(mutcounts[,c(6,5)])

####################
## Modify Muttype ##
####################

## change cnv_rna to cnv
mutcounts[which(mutcounts$type_mod == "cnv_rna"),"type_mod"] <- gsub("_rna$", "", mutcounts[which(mutcounts$muttype == "cnv_rna"),"muttype"])

## change naming
mutcounts$type_mod <- toupper(mutcounts$type_mod)
mutcounts[which(mutcounts$type_mod == "GERMLINE"),"type_mod"] <- "GL"
mutcounts[which(mutcounts$type_mod == "RNA"),"type_mod"] <- "EXP"

#############################
## Add zero count patients ## 
#############################

marie <- read.csv("Data/landscape_samples_mutations - marie_concat_prism_summ_v3.csv" , sep=",", header=T)
zero <- as.vector(marie$sample_id[-which(marie$sample_id %in% unique(mutcounts$sample))])
z <- zero[-which(zero == "")]
pz <- paste0('"', paste(z, collapse='", "'), '"')
print(pz, quote=F)

zs_df <- as.data.frame(z)
zs_df$type_mod <- NA
zs_df$mutcount <- 0
zs_df$platforms <- NA
zs_df$muttype <- NA
zs_df$platform_mod <- NA
colnames(zs_df)[1] <- "sample_id"

mutcounts <- rbind(mutcounts, zs_df)

##########################
## get total per sample ##
##########################

ps <-split(mutcounts, f=mutcounts$sample_id)
psl <- lapply(ps, function(x) sum(x$mutcount))
psm <- do.call(rbind, psl)
colnames(psm) <- "total"
fo <- order(psm[,1])

plot_all <- merge(mutcounts, psm, by.x="sample_id", by.y="row.names")
head(plot_all)

###################################
###################################
## Plot A -  individual patients ##
###################################
###################################

## third most abundant platfomr is WR
plattot <- table(plot_all$platform_mod)
ord <- names(plattot)[order(plattot, decreasing =T )]

plattot[order(plattot, decreasing =T )]
#WR   W   R   M   G 
#241 169 124  55  34 

out_df <- plot_all
for(cat in 1:length(ord)){
  print(cat)
  wpr <- out_df[which(out_df$platform_mod == ord[cat]),]
  wprs <- split(wpr,f=as.vector(wpr$sample_id))
  wprsl <- lapply(wprs, function(x) sum(x$mutcount))
  wprsm <- do.call(rbind, wprsl)
  out_df <- merge(out_df, wprsm, by.x="sample_id", by.y="row.names", all=T)
  colnames(out_df)[ncol(out_df)] <- ord[cat]
  out_df[,ncol(out_df)][is.na(out_df[,ncol(out_df)])] <- 0
}


plot_mcs <- out_df[with(out_df, order(total,WR,W,R,M,G,decreasing=T)), ] 

## so I don't have to change code below
plot_mcs <- plot_mcs[,c(1,4:ncol(plot_mcs))]
colnames(plot_mcs)[4] <- "platform"

plot_mcs$sample_id<- factor(plot_mcs$sample_id, ordered = FALSE )
plot_mcs$sample_id <- ordered(plot_mcs$sample_id, levels = unique(plot_mcs$sample_id) )

plot_mcs$platform <- factor(plot_mcs$platform, ordered = FALSE )
plot_mcs$platform <- ordered(plot_mcs$platform, levels=rev(ord))

## change colour
col_val <- c(pal_npg("nrc")(10)[c(6,8,5,7,2)])
col_nam <- c("W","WR","R","G", "M")
names(col_val) <- col_nam


g6 <- ggplot(plot_mcs, aes(x=sample_id, y=mutcount, fill=platform), ordered=T) + geom_bar(position="stack", stat="identity") + theme_bw() + scale_fill_manual(breaks=c("G", "M", "R", "W", "WR"), labels=c("W" = "WGS only", "WR" = "WGS + RNA", "R" ="RNA", "G" = "Germline WGS", "M" = "Methylation only"), values=col_val) + labs(fill="Platform",y="Number of Reportable Findings", x="Patients") + scale_y_continuous(expand=c(0,0),breaks =c(5,10,15),limits=c(0,14)) + theme(axis.text.x = element_blank(),panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) + geom_hline(yintercept=1, linetype="dashed", color = "black") +   theme(
  axis.text.x = element_blank(),
  axis.title.y = element_text(angle = 90, vjust = 0.5, size=18),
  axis.title.x = element_text(hjust = 0.5, size=18),
  plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
  panel.spacing = unit(0, "lines"),
  panel.grid.major.x=element_blank(), 
  panel.grid.minor.x=element_blank(), 
  text = element_text(size = 20),
  legend.title=element_blank(),
  legend.justification=c(1,0), 
  legend.position=c(0.95, 0.75),  
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.box.margin=margin(c(0,0,0,0)),
  legend.text=element_text(size=20)
)

ggsave("Plots_output/Supplementary_figure10A.pdf", g6, device="pdf", units="in", width=19, height=14, limitsize=FALSE, dpi=330)

######################
######################
## Plot B - Summary ##
######################
######################

## this I reinvented the wheel with this - could have just used table
count.f <- function(x) {
  ## to padd out results for easy merge later
  null_vec <- c(0,0,0,0,0,0,0) 
  names(null_vec) <- c("SNV","METH","SV","CNV","GL","TMB","EXP") 
  null_df <- as.data.frame(null_vec)
  
  ## get counts per type
  y <- split(x,f=x$type_mod) 
  yz <- lapply(y, function(z) sum(z$mutcount))
  yzu <- as.data.frame(unlist(yz))
  colnames(yzu) <- "count"
  
  ## merge with 0 to get padded results
  m <- merge(null_df, yzu, by="row.names", all=T)
  rownames(m) <- m[,1]
  m <- m[,-1]
  out <- rowSums(m, na.rm=T)
  return(out)
  
}

## total reportable per category
sub <- plot_all[, c("platform_mod" ,"type_mod", "mutcount")]
subs <- split(sub, f=as.factor(as.vector(sub$platform_mod)))
subsl <- lapply(subs, count.f )

## ploting df
subsm <- as.data.frame(do.call(rbind, subsl))
subsm$plt <- rownames(subsm)
subsmlt <- melt(subsm, id=c("plt"))
subsmlt$plt <- factor(subsmlt$plt, levels = c("WR","W","R","M","G"))

## remove 0 values 
subsmlt <- subsmlt[-which(subsmlt$value == 0),]

## Set levels
subsmlt$variable <- as.vector(subsmlt$variable)
subsmlt$variable <- factor(subsmlt$variable, levels = c("EXP", "SNV", "CNV", "SV", "GL", "METH", "TMB"))

## change colour
col_val <- c(pal_npg("nrc")(10)[c(6,8,5,7,2)])
col_nam <- c("W","WR","R","G", "M")
names(col_val) <- col_nam

## labels for x-axis
labs <- c("WR" = "WGS\n +\n RNA", "W" = "WGS", "R" = "RNA", "G" = "WGS", "M" = "METH")

##create plot
g9e <- ggplot(subsmlt, aes(x=plt, y=value, fill=plt)) + geom_col(position = position_dodge2(width = 0.9, preserve = "single")) + facet_grid(~variable, scales="free", space="free",labeller=label_wrap_gen(width=15)) + scale_fill_manual(values=col_val) + theme_bw() + theme(axis.title.x= element_blank(), legend.position="none",panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(), text = element_text(size = 12),panel.spacing.x = unit(0.25, "lines")) + labs(y="Number of Reportable Findings") + scale_x_discrete(labels= labs)

ggsave("Plots_output/Supplementary_figure10B.pdf", g9e, device="pdf", width=6, height=4, units="in")
