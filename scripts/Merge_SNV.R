
patients<-c("...") #insert comma separated list of patients to process

for(i in 1:length(patients)){
  pid<-strsplit(patients[i],'-')
  id<-pid[[1]][1]
  dirLoc<-paste("./SNP_analysis/",id,"/",sep="")
  setwd(dirLoc)
  patientLabel<-patients[i]
  filename=paste(patientLabel,".vaf.txt",sep="")
  rawdata<-read.delim(filename,sep="\t",header=T)
  pos<-rawdata$Chrom.Pos
  
  filename=paste(patientLabel,"_snp.txt",sep="")
  if(file.exists(filename)){
    snp<-read.delim(filename,sep="\t",header=T)
    gl<-snp$Genomic.Location
    gn<-snp$Gene
    
    index<-which(pos %in% gl)
    data<-rawdata[index,]
    
    filename<-paste("./ExpressionAnalysis/",id,"/",patientLabel,"-FC.txt",sep="")
    exp<-read.delim(filename,sep="\t",header=T)
    subexp<-exp[,c(1,4)]
    gene<-subexp$gene_id
    index<-which(gene %in% gn)
    tpm<-subexp[index,]
    
    snp$id<-rownames(snp)
    snpExp<-merge(snp,tpm,by.x=c("Gene"),by.y=c("gene_id"),all.x=T,sort=F)
    ordSnpExp<-snpExp[order(as.numeric(snpExp$id)),]
    
    snpExpData<-merge(ordSnpExp,data,by.x=c("Genomic.Location"),by.y=c("Chrom.Pos"),all.x=T,sort=F)
    combData<-snpExpData[order(as.numeric(snpExpData$id)),c(1,2,7,4,9)]
    filename=paste(patientLabel,"_snp_rna.txt",sep="")
    write.table(combData,filename,sep="\t",row.names=F)
  }
}


