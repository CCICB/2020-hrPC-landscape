library(tidyr);library(dplyr);library(devtools);library(sunburstR);library(htmltools);library(readxl)
library(d3r);library(treemap);library(mapview);library(ggsci)


##load clinical
file<-"SupplementaryTable1_revision.xlsx"
clinFile<-read_excel(file,"Samples")
clin<-as.data.frame(clinFile,stringsAsFactors=F)

###get assays per model
clin$assay <- ifelse(clin$Tumour.WGS != "NULL" & clin$RNAseq == "FALSE" & clin$Methylation == "FALSE", "WGS",
                     ifelse(clin$Tumour.WGS != "NULL" & clin$RNAseq == "TRUE" & clin$Methylation == "FALSE", "WGS,RNA",
                            ifelse(clin$Tumour.WGS == "NULL" & clin$RNAseq == "TRUE" & clin$Methylation == "FALSE", "RNA",
                                   ifelse( clin$Tumour.WGS != "NULL" & clin$RNAseq == "TRUE" & clin$Methylation == "TRUE", "WGS,RNA,Methylation",
                                           ifelse( clin$Tumour.WGS != "NULL" & clin$RNAseq == "FALSE" & clin$Methylation == "TRUE", "WGS,Methylation",
                                                   ifelse( clin$Tumour.WGS == "NULL" & clin$RNAseq == "TRUE" & clin$Methylation == "TRUE", "RNA,Methylation",
                                                           ifelse( clin$Tumour.WGS == "NULL" & clin$RNAseq == "FALSE" & clin$Methylation == "TRUE", "Methylation",
                                                                   "")))))))
clin$RNAseq <- gsub("FALSE", "", clin$RNAseq)
clin$RNAseq <- gsub("TRUE", "RNA", clin$RNAseq)
clin$Tumour.WGS <- gsub("FALSE", "", clin$Tumour.WGS)
clin$Tumour.WGS <- gsub("TRUE", "WGS", clin$Tumour.WGS)
clin$Methylation <- gsub("FALSE", "", clin$Methylation)
clin$Methylation <- gsub("TRUE", "Methylation", clin$Methylation)


###subset clinical
clin.sub <- clin[,c("sample", "cancer.type", "stage", "Tumour.WGS","RNAseq", "Methylation")]
clin.sub$nodes <- paste0(clin.sub[[1]], clin.sub[[2]], clin.sub[[3]], clin.sub[[4]], clin.sub[[5]])

###make node patterns
clin.sub[,1]<- NULL
names(clin.sub) <- c("level1", "level2", "level3", "level4", "level5", "nodes")
counts <- clin.sub %>% 
  group_by(nodes) %>%
  dplyr::summarise(size = n())

final <- merge(clin.sub, counts, all.x = T)
final$nodes <- NULL
write.table(final, "output/counts.csv", col.names = T, row.names = F, quote = F, sep = ",")

###create treemap
sunburst(
  d3_nest(dat, value_cols="size"),
  count = TRUE, percent = F,
  legend = T, withD3 = T)
tm <- treemap(dat,
              index = c("level1", "level2", "level3", "level4", "level5"),
              vSize = "size",
              draw = FALSE
)$tm

###subset for colors
level4 <- subset(tm, level == 4)
level5 <- subset(tm, level == 5)
level3 <- subset(tm, level == 3)
level2 <- subset(tm, level == 2)
level1 <- subset(tm, level == 1)

###color dx/relapse/broad hist/assay
level2$color <- ifelse(level2$level2 == "Diagnosis", "#87CEFF", 
                       ifelse(level2$level2 == "Relapse", "#EEAEEE",
                              ifelse(level2$level2 == "Refractory", "#AB82FF", "#9AFF9A")))
level1$color <- ifelse(level1$level1 == "Other", "#F39B7FFF",
                       ifelse(level1$level1 == "Leukaemia", "#4DBBD5B2", 
                              ifelse(level1$level1 == "Lymphoma","#4DBBD5B2",
                                     ifelse(level1$level1 == "NBL", "#00A087FF",
                                            ifelse(level1$level1 == "Sarcoma", "#3C5488FF", "#E64B35FF")))))
level4$color <- ifelse(level4$level4 == "RNAseq", "#CCCCCC", "#FFFFFF")
level3$color <- ifelse(level3$level3 == "Tumour.WGS", "#969696", "#FFFFFF")
level5$color <- ifelse(level5$level5 == "Methylation", "#E0E0E0", "#FFFFFF")


new.col <- bind_rows(list(level1, level2, level3, level4, level5))
##plot pie with new colors
p<- sund2b(
  d3_nest(new.col, value_cols = colnames(new.col)[-(1:5)]),
  colors = htmlwidgets::JS(
    "function(name, d){return d.color || '#ccc';}"
  ),
  valueField = "vSize"
)
 
mapshot(p, url="output/sampleOverview_pie.html")
