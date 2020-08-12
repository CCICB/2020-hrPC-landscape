##############
## Figure 1 ##
##############

library(ggplot2)
library(reshape)
library(grid)

###############################
## Plot A - Patient age plot ##
###############################

marie <- read.csv("Data/landscape_samples_mutations - marie_concat_prism_summ_v3.csv" , sep=",", header=T, stringsAsFactors = F)
marie2 <- marie[-which(marie$patient_id == ""),]

sub <- marie2[,c("event","age_at_sample")]
sub$event[which(sub$event == "Refractory")] <- "Progression"

## specify colours
Diagnosis <- "#87CEFF"
Refractory <- "#AB82FF"
Relapse <- "#EEAEEE"
Secondary <- "#9AFF9A"

## colour values
col_val <- c(Diagnosis,Refractory,Relapse,Secondary)
names(col_val) <- c("Diagnosis", "Progression", "Relapse", "Secondary")

## bin all age > 25 - change this on x-axis
sub$age_at_sample <- as.numeric(sub$age_at_sample)
sub[which(sub$age_at_sample > 24), "age_at_sample"] <- 25

subdff <- as.data.frame(table(sub), stringsAsFactors = F)
r <- c("Diagnosis", 24, 0)
r2 <- c("Progression", 24, 0)
r3 <- c("Relapse", 24, 0)
r4 <- c("Secondary", 24, 0)
subdff <- rbind(subdff,r,r2,r3,r4)

## order
subdff$event <- ordered(subdff$event,levels=unique(subdff$event))
subdff$age_at_sample <- ordered(subdff$age_at_sample,levels=sort(as.numeric(unique(subdff$age_at_sample))))
subdff$Freq <- as.numeric(subdff$Freq)

## plot
g <- ggplot(subdff, aes(x=age_at_sample,y=Freq, fill=event))  + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=col_val) + labs(fill="Stage",y="Number of Patients", x="Age (years)") + scale_y_continuous(expand=c(0,0), limits=c(0,21)) + theme_bw() + scale_x_discrete(breaks=seq(0,25),labels=c(seq(0,24), ">25")) + theme(legend.position="none") + theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

ggsave("Final_plots/Plots_output/Figure1A.pdf", g)

###############################
## Plot C - Patient age plot ##
###############################

manifest <- read.csv("Final_plots/Data/ZCC_Landscape_ManifestInfo - ZCC_Landscape_ManifestInfo_updated_v4.csv", sep=",", header=T)

sub <- manifest[,c("Patient.ID","Paper.category", "Diagnosis")]

sum <- table(sub[,2:3])
summ <- melt(sum)

## ordering
subs <- split(summ, f=summ$Paper.category)

names <- vector()
for(x in 1:length(subs)){
  y <- subs[[x]][which(subs[[x]]$value != 0),]
  o <-  y[order(y$value, decreasing=T),]
  n <- as.vector(o$Diagnosis)
  names <- c(names, n)
}

## remove 0 categories 
summ <- summ[-which(summ$value == 0),]

sub$Diagnosis <- factor(sub$Diagnosis, levels = rev(names))
summ$Diagnosis <- factor(summ$Diagnosis, levels = rev(names))

###########################
## colour by cancer type ##
###########################

## import colours
cancertype_colours <- read.table("Final_plots/Data/Figure1C_cancertype_colours.txt", header=T, sep="\t", comment.char = "")

cancersub_colours <- read.table("Final_plots/Data/Figure1C_PaperCategoryHexColours.txt", header=T, sep="\t", comment.char = "")

## create colour vector manually 
##https://gist.github.com/Jfortin1/72ef064469d1703c6b30
lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

col_val <- as.vector(cancersub_colours$HexCol)

colc <- list()
nam <- list()
for(i in 1:length(col_val)){
  l <- length(which(sum[i,] != 0))
  alc <- vector()
  nam[[i]] <- names(which(sum[i,] != 0))
  for(n in 1:l){
    x <- lighten(col_val[i], factor=(n-0.5)/l)
    alc <-c(alc,rev(x))
  }
  colc[[i]] <- alc
}

col_val <- unlist(colc)
names(col_val) <- names

################################################
## Make bars narrower - flip test, formatting ##
################################################

p2 <- ggplot(summ,aes(x=factor(Paper.category), y=value, fill=Diagnosis, label=Diagnosis)) + facet_wrap(~Paper.category, strip.position = "top", nrow=1, scales="free_x")  + geom_bar( position="stack", stat="identity") + scale_fill_manual(values=col_val,guide = guide_legend(ncol=1, reverse = TRUE)) + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.border = element_blank(), axis.title.x= element_blank(),  panel.background = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", aspect.ratio = 6/1) + geom_text(size = 3, position = position_stack(vjust = 0.5), stat = 'identity', check_overlap = TRUE) + scale_x_discrete(expand = c(0,0)) + scale_y_reverse() + labs(y="Patient Counts")


## from https://htykuut.blogspot.com/2019/01/r-ggplot2-change-colour-of-font-and.html#
g <- grid.force(ggplotGrob(p2))

# Get the names of grobs and their gPaths into a data.frame structure
grobs_df <- do.call(cbind.data.frame, grid.ls(g, print = FALSE))

# Build optimal gPaths that will be later used to identify grobs and edit them
grobs_df$gPath_full <- paste(grobs_df$gPath, grobs_df$name, sep = "::")

grobs_df$gPath_full <- gsub(pattern = "layout::", replacement = "", x = grobs_df$gPath_full, fixed = TRUE)

# Get the gPaths of the strip background grobs
strip_bg_gpath <- grobs_df$gPath_full[grepl(pattern = ".*strip\\.background.*", x = grobs_df$gPath_full)]

# example of a gPath for strip background 
strip_bg_gpath[1] 

# Get the gPaths of the strip titles
strip_txt_gpath <- grobs_df$gPath_full[grepl(pattern = "strip.*titleGrob.*text.*", x = grobs_df$gPath_full)]

# example of a gPath for strip title
strip_txt_gpath[1] 


# Generate some color
n_cols <- length(cancersub_colours$HexCol)
fills <- cancersub_colours$HexCol
txt_colors <- rep("#FFFFFF", n_cols)

# Edit the grobs
for (i in 1:length(strip_bg_gpath)){
  g <- editGrob(grob = g, gPath = strip_bg_gpath[i], gp = gpar(fill = fills[i]))
  g <- editGrob(grob = g, gPath = strip_txt_gpath[i], gp = gpar(col = txt_colors[i]))
}

grid.newpage()
grid.draw(g)

ggsave("Final_plots/Plots_output/Figure1C.pdf", g)

