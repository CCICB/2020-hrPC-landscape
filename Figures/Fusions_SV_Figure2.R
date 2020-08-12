library(ggplot2);library(ggsci);library(forcats)

sv<-read.delim("./ZCC-Reportable.txt",sep="\t",header=T,stringsAsFactors=F)
subSv<-sv[,c("sample","fus_sv_name","fus_sv_class","cancer.type")]

index<-which(subSv$fus_sv_class == "Fusion")
fus<-subSv[index,]
othSv<-subSv[-index,]

colsCnc<-c(CNS = "#E64B35FF", HM = "#4DBBD5B2", NBL = "#00A087FF", Sarcoma = "#3C5488FF", `Solid tumour` = "#F39B7FFF")

p<-ggplot(fus,aes(fct_infreq(fus_sv_name)))+
  geom_bar(stat="count",width=0.7,aes(fill=Paper.category))+
  scale_y_continuous(expand=c(0,0),limits=c(0,12))+
  scale_fill_manual(values=colsCnc,name="")+
  theme_classic()+
  theme(legend.position=c(0.85,0.9),legend.key.size=unit(0.8,"lines"),axis.text.x=element_text(hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=7),axis.title=element_blank())
p+coord_flip()
ggsave("FusionCounts_Figure2.pdf",p,width=95,height=85,units="mm",useDingbats=F)


colsType<-c(TS = "#8491B4FF", ACT = "#B09C85FF")

p<-ggplot(othSv,aes(fct_infreq(fus_sv_name)))+
  geom_bar(stat="count",width=0.7,aes(fill=fus_sv_class))+
  scale_y_continuous(expand=c(0,0),limits=c(0,7))+
  scale_fill_manual(values=colsType,name="",labels=c("Activating","Damaging"))+
  theme_classic()+
  theme(legend.position=c(0.8,0.95),legend.key.size=unit(0.8,"lines"),axis.text.x=element_text(hjust=1,vjust=0.5,size=10),axis.text.y=element_text(size=8),axis.title=element_blank())+
  labs(x="Rearrangement",y="Count")
p+coord_flip()
ggsave("SVCounts_Figure2.pdf",width=95,height=85,units="mm",useDingbats=F)



