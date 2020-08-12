library(ComplexHeatmap);library(data.table);library(ggsci);library(readxl)

reportTable<-read.delim("./ZCC-Reportable.txt",sep="\t",header=T,stringsAsFactors=F)
oncomat<-data.table::dcast(data=reportTable[,.(gene,var_type,sample)],formula=gene ~ sample,
                           fun.aggregate=function(x){
                             x=unique(as.character(x))
                             return(x)
                           }, value.var='var_type',fill=';',drop=FALSE)

write.table(oncomat,"Oncoprint_matrix.txt",sep="\t",row.names=F,quote=F)

#order oncomat by pathways - this is according to the order in the below file
genePath<-read.delim("./Gene_Pathway.txt",sep="\t",header=T)
reOrd<-match(genePath$gene,oncomat$gene)
newMat<-oncomat[reOrd,]
oncomat<-newMat

mat<-as.matrix(oncomat)
rownames(mat)<-mat[,1]
mat<-mat[,-1]
mat[is.na(mat)]<-""

col=c("HOMDEL"="blue","DELETION"="lightskyblue1","GAIN"="indianred1","AMP"="red2","SNV"="mediumseagreen","INDEL"="yellow","Fusion"="purple","TS"="orchid4","ACT"="orchid","GERM"="black","UPREG"="red4","DOWNREG"="navy")
mutLabels=c("Biallelic Deletion","Deletion","Gain","Amplification","SNV","Indel","Fusion","Damaging SV","Activating SV","Germline","Over-Expressed","Under-Expressed")
alter_fun=function(x,y,w,h,v){
  if(v["HOMDEL"]) grid.rect(x,y.w*0.9,h*0.9,gp=gpar(fill=col["HOMDEL"],col=NA))
  if(v["DEL"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["DEL"],col=NA))
  if(v["GAIN"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["GAIN"],col=NA))
  if(v["AMP"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["AMP"],col=NA))
  if(v["SNV"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["SNV"],col=NA))
  if(v["INDEL"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["INDEL"],col=NA))
  if(v["Fusion"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["Fusion"],col=NA))
  if(v["TS"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["TS"],col=NA))
  if(v["ACT"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["ACT"],col=NA))
  if(v["GERM"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["GERM"],col=NA))
  if(v["UPREG"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["UPREG"],col=NA))
  if(v["DOWNREG"]) grid.rect(x,y.w*0.9,h*0.4,gp=gpar(fill=col["DOWNREG"],col=NA))
}

#INSERT ANNOTATION INFORMATION
file<-"SupplementaryTable1_revision.xlsx"
histFile<-read_excel(file,"Samples")
histData<-as.data.frame(clinFile,stringsAsFactors=F)

cncCol=c("CNS"="#E64B35FF","HM"="#4DBBD5FF","NBL"="#00A087FF","Sarcoma"="#3C5488FF","Solid tumour"="#F39B7FFF")
sexCol=c("Female"="#FFB6C1","Male"="#BCD2EE")
eventCol<-c("Diagnosis"="#87CEFF","Refractory"="#AB82FF","Relapse"="#EEAEEE","Secondary"="#9AFF9A")

Stage<-histData$stage
CancerType<-histData$cancer.type
Gender<-histData$sex

dfHist<-cbind.data.frame(CancerType,Stage,Gender)
heatAnno=HeatmapAnnotation(df=dfHist,col=list(CancerType=cncCol,Stage=eventCol,Gender=sexCol),
                           annotation_height=1, na_col="snow1")


oncoPlot<-oncoPrint(mat,get_type=function(x) strsplit(x,";")[[1]],
                    alter_fun=function(x,y,w,h,v){
                      n=sum(v)
                      h=h*0.75
                      grid.rect(x,y,w,h,gp=gpar(fill="snow1",col=NA))
                      if(n) grid.rect(x,y-h*0.5+1:n/n*h,w*0.75,1/n*h,
                                      gp=gpar(fill=col[names(which(v))],col="NA"),just="top")
                      
                    },
                    col=col,
                    heatmap_legend_param=list(title="Alteration",at=names(col),labels=mutLabels),
                    show_column_names=F,
                    remove_empty_columns = T,
                    row_order=1:nrow(mat),
                    split=genePath$Pathway,
                    row_names_gp=gpar(fontsize=6),
                    pct_gp=gpar(fontsize=6),
                    bottom_annotation=heatAnno
                    
)
pdf("Oncoprint_SupFig2.pdf",height=13,width=9.22)
draw(oncoPlot)
dev.off()

