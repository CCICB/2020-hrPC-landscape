library(ggplot2)

#code needs to be run after Oncoplot_SupplementaryFigure2.R as pulls from the matrix
oncomat<-read.delim("./Oncoprint_matrix.txt",sep="\t",header=T,stringsAsFactors=F)
colnames(oncomat)<-gsub(pattern="\\.",replacement="-",colnames(oncomat))

genePath<-read.delim("Gene_Pathway_muts.txt",sep="\t",header=T,stringsAsFactors=F)
reOrd<-match(genePath$gene,oncomat$gene)
newMat<-oncomat[reOrd,]
oncomat<-newMat

mutDat<-data.frame(gene=as.character(),mutType=as.character(),pathway=as.character(),stringsAsFactors=F)
k=1
for(i in 1:nrow(oncomat)){
  geneName<-oncomat[i,1]
  index<-which(genePath$gene==geneName)
  pathName<-genePath[index,2]
  for(j in 2:ncol(oncomat)){
    if(oncomat[i,j]=="" || is.na(oncomat[i,j])){
      
    }else{
      entry<-strsplit(oncomat[i,j],";")
      for(a in 1:length(entry[[1]])){
        mutDat[k,1]<-geneName
        mutDat[k,2]<-entry[[1]][a]
        mutDat[k,3]<-pathName
        k<-k+1
      }
    }
  }
}

mutDat$gene<-factor(mutDat$gene,levels=genePath$gene)

col=c("HOMDEL"="blue","DELETION"="lightskyblue1","GAIN"="indianred1","AMP"="red2","SNV"="mediumseagreen","INDEL"="yellow","Fusion"="purple","TS"="orchid4","ACT"="orchid","GERM"="black","UPREG"="red4","DOWNREG"="navy")
mutLabels=c("Activating SV","Amplification","Deletion","Under-Expressed","Fusion","Gain","Germline","Biallelic Deletion","Indel","SNV","Damaging SV","Over-Expressed")

png("GeneMutations_v3.png",res=100,height=600,width=1200)
p<-ggplot(mutDat,aes(x=gene,fill=mutType),faceting=T,facetingVarNames=pathway)+
  geom_bar(stat="count",width=0.8)+
  scale_fill_manual(values=col,labels=mutLabels,name="Alteration")+
  scale_y_continuous(expand=c(0,0),lim=c(0,65))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=7),axis.text.y=element_text(size=10),axis.title=element_text(size=14,face="bold"),
        panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  labs(x="Gene",y="Count")
p+coord_flip()
ggsave("GeneMutations_SupFig3.pdf",p,useDingbats=F)
dev.off()


genePath<-read.delim("Gene_Pathway_muts_modify.txt",sep="\t",header=T,stringsAsFactors=F)
reOrd<-match(genePath$gene,oncomat$gene)
newMat<-oncomat[reOrd,]
oncomat<-newMat

mutDat<-data.frame(gene=as.character(),mutType=as.character(),pathway=as.character(),stringsAsFactors=F)
k=1
for(i in 1:nrow(oncomat)){
  geneName<-oncomat[i,1]
  index<-which(genePath$gene==geneName)
  pathName<-genePath[index,2]
  for(j in 2:ncol(oncomat)){
    if(oncomat[i,j]=="" || is.na(oncomat[i,j])){
      
    }else{
      entry<-strsplit(oncomat[i,j],";")
      for(a in 1:length(entry[[1]])){
        mutDat[k,1]<-geneName
        mutDat[k,2]<-entry[[1]][a]
        mutDat[k,3]<-pathName
        k<-k+1
      }
    }
  }
}

mutDat$gene<-factor(mutDat$gene,levels=genePath$gene)
uniqPath<-unique(genePath$Pathway)
mutDat$pathway<-factor(mutDat$pathway,levels=uniqPath)

col=c("HOMDEL"="blue","DELETION"="lightskyblue1","GAIN"="indianred1","AMP"="red2","SNV"="mediumseagreen","INDEL"="yellow","Fusion"="purple","TS"="orchid4","ACT"="orchid","GERM"="black","UPREG"="red4","DOWNREG"="navy")
mutLabels=c("Activating SV","Amplification","Deletion","Under-Expressed","Fusion","Gain","Germline","Biallelic Deletion","Indel","SNV","Damaging SV","Over-Expressed")

p<-ggplot(mutDat,aes(x=gene,fill=mutType),faceting=T,facetingVarNames=pathway)+
  geom_bar(stat="count",width=0.8)+
  scale_fill_manual(values=col,labels=mutLabels,name="Alteration")+
  scale_y_continuous(expand=c(0,0),lim=c(0,65))+
  facet_wrap(~pathway,nrow=1,scales="free_x",labeller=label_wrap_gen(width=10))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10),axis.title=element_text(size=14,face="bold"),
        strip.text.x=element_text(size=12),panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
  labs(x="Gene",y="Count")
p
ggsave("GeneMutations_truncated_Fig2.pdf",p,useDingbats=F)


