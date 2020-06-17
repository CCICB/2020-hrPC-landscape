##############################################################
### Supplementary Figure 6 code - methylation vs diagnosis ###
##############################################################

library(ggplot2)
library(ggsci)

#############
## updated ##
#############

final3 <- read.table("Data/Methylation_Match.txt", sep="\t", header=T, stringsAsFactors = F)

diag_plot_data <- final3[,c("sample_id","meth_vec","diag_vec")]
diag_plot_data$diag_vec[which(diag_plot_data$diag_vec == "Low_Meth_Score")] <- "No Match"

diag_plot_data$meth_vec <- factor(diag_plot_data$meth_vec, levels=c("NOMatch", "Match_Weak","Match_Strong"))

fils <- c(pal_npg("nrc")(3), "#D3D3D3")
cols <- c("Same" = fils[3], "Different" = fils[2], "Changed" = fils[1], "No Match" = fils[4])

mc <- ggplot(diag_plot_data, aes(factor(meth_vec), fill=diag_vec)) +  geom_bar() +  scale_fill_manual(breaks=c("Same","Different","Changed","No Match"), values=cols, guide_legend(title = "Concordance type", reverse=TRUE)) + theme_bw() + coord_flip() + scale_x_discrete(labels = c('No Match [0.0-0.5)','Weak Match [0.5 -0.9)','Strong Match [0.9-1.0)')) + labs(x="Cancer Subtype classification", y="Number of Patients") + theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank())

ggsave("Plots_output/SupplementaryFigure6.pdf", plot=mc, width=6.4, height=3.8, units= "in")
