library(ggplot2);library(ggsci);library(ggbeeswarm);library(readxl)

###Load RNA TPM data
file<-"SupplementaryTable1_revision.xlsx"
rnaFile<-read_excel(file,"RNA")
rnaTpm<-as.data.frame(rnaFile,stringsAsFactors=F)
reportTable<-read.delim("./ZCC-Reportable.txt",sep="\t",header=T,stringsAsFactors=F)


#subset to pull out only top genes of interest as previously identified
goi<-c("PDGFRA","CDK4","MYC","KRAS","KIT","MDM2","CCNE1","ERBB2","NRAS","STAT5B","VEGFA","JAK3","PHOX2B","BCL2","JAK1","JAK2","MTOR","PIK3CD","CDKN2A","TP53","SMARCB1")
gn<-rnaTpm$gene
index<-which(gn %in% goi)
subRnaTpm<-rnaTpm[index,]

#merge RNA data with reportable
rnaRep<-merge(subRnaTpm,reportTable,by=c("sample"),all=F)

#order genes for plotting
subRnaTpm$gene<-factor(subRnaTpm$gene,levels=c("PDGFRA","CDK4","MYC","KRAS","KIT","MDM2","CCNE1","ERBB2","NRAS","STAT5B","VEGFA","JAK3","PHOX2B","BCL2","JAK1","JAK2","MTOR","PIK3CD","CDKN2A","TP53","SMARCB1"))

colsHit<-c(RNA="#DC0000FF",CNV="#4DBBD5FF",SNV="#7E6148FF",FUS="#3C5488FF",Multi="#F39B7FFF",AN="grey87")


p<-ggplot(tGoiTpm)+
  geom_hline(aes(yintercept=-2))+
  geom_hline(aes(yintercept=2))+
  geom_jitter(aes(gene,zscore,colour=typeOfHit,size=reportable),width=0.27,height=0.5)+
  scale_colour_manual(values=colsHit,labels=c("Not Reported","CNV","Fusion","Multi","Expression Only","SNV"),name="Alterations")+ 
  scale_y_continuous(breaks=c(5,2,0,-2,-5))+
  scale_size_manual(values=c(0.4,1.4),guide=FALSE)+
  theme_bw()+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=10),axis.title=element_text(size=12,face="bold"))+
  labs(x="Gene",y="Z-score")
p+coord_flip()
ggsave("Fig2_RnaExpression.pdf",p,useDingbats=F)
dev.off()

