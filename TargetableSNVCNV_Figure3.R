library(ggplot2);library(ggsci);library(forcats)

CNV_target<-read.delim("CNV_target.txt",sep="\t",header=T,stringsAsFactors=F)

df<-CNV_target %>% separate_rows(CNV_reprt_2,sep=",",convert=FALSE)
df2<-df %>% separate(CNV_reprt_2,c("CNV_reprt_2_Gene","CNV_reprt_2_CNVtype"),sep="_",convert=FALSE)

cols<-c(CNS="#E64B35FF",HM="#4DBBD5FF",NBL="#00A087FF",Sarcoma="#3C5488FF",`Solid tumour`="#F39B7FFF")

p<-ggplot(df2, aes(fct_rev(fct_infreq(CNV_reprt_2_Gene)), ordered = TRUE, fill=cancer_category)) + 
  geom_bar(width = 0.5) + 
  labs(x = "Targetable CNV", y = "frequency") + 
  scale_x_discrete(expand = c(0.02, 0.02)) + 
  scale_y_continuous(expand = c(0, 0.5)) + 
  scale_fill_manual(name = "cancer type", values = cols) +
  theme_bw() +
  theme(
    legend.title=element_blank(),
    legend.position=c(0.9, 0.02),
    legend.justification = c("right", "bottom"),
    legend.box.margin=margin(c(0,0,0,0)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    strip.text.x = element_text(size=10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
  ) +
  coord_flip() + 
  labs(title = "") +
  guides(fill=guide_legend(ncol=1))
p
ggsave("CNV_reportable_Fig3.pdf",p,width=100,height=200,units="mm",useDingbats=F)


####SNV
SNV_target<-read.delim("SNV_action.txt",sep="\t",header=T,stringsAsFactors=F)

df<-SNV_target %>% separate_rows(SNV_Gene,sep=",",convert=FALSE)

cols<-c(CNS="#E64B35FF",HM="#4DBBD5FF",NBL="#00A087FF",Sarcoma="#3C5488FF",`Solid tumour`="#F39B7FFF")

p <- ggplot(df, aes(fct_rev(fct_infreq(SNV_Gene)), ordered = TRUE, fill=cancer_category)) + 
  geom_bar(width = 0.5) + 
  labs(x = "Targetable SNV", y = "frequency") + 
  scale_x_discrete(expand = c(0.02, 0.02)) + 
  scale_y_continuous(expand = c(0, 0.5),breaks=c(0,2,4,6,8,10)) + 
  scale_fill_manual(name = "cancer type", values = cols) +
  theme_bw() +
  theme(
    legend.title=element_blank(),
    legend.position=c(0.9, 0.02),
    legend.justification = c("right", "bottom"),
    legend.box.margin=margin(c(0,0,0,0)),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x=element_blank(),
    axis.title.y = element_text(angle = 90, vjust = 0.5),
    strip.text.x = element_text(size=10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
  ) +
  coord_flip() + 
  labs(title = "") +
  guides(fill=guide_legend(ncol=1))
p
ggsave("SNV_reportable_Fig3.pdf", p, width=100, height=200, units="mm",useDingbats=F)