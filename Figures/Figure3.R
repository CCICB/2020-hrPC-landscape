setwd("C:/Users/mpinese/Downloads/gl_summary")

# Wanted:
# 1. Frequency ~ gene, coloured by cancer type
# 2. Frequency ~ cancer type, coloured by cancer type
# 3. Estimate of how often WGS-detected GL vars were seen in panel seq

library(openxlsx)

gl_vars = read.xlsx("query_result_combined_annotated_201908141543_NF_MPlong2.xlsx")

table(gl_vars$paper_cat[gl_vars$include_in_gl_stats == 1], gl_vars$germline_reported_or_reportable_now[gl_vars$include_in_gl_stats == 1])
plot(table(gl_vars$paper_cat[gl_vars$include_in_gl_stats == 1], gl_vars$germline_reported_or_reportable_now[gl_vars$include_in_gl_stats == 1]))

library(plyr)
collapsed = ddply(gl_vars[gl_vars$include_in_gl_stats == 1,], .(patient_id), function(d) {
    data.frame(paper_cat = d$paper_cat[[1]], germline_reported_or_reportable_now = any(d$germline_reported_or_reportable_now==1), nrep = sum(d$germline_reported_or_reportable_now==1))})
collapsed$paper_cat = factor(collapsed$paper_cat, levels = c("CNS", "HM", "NBL", "Sarcoma", "Solid tumour"))
table(collapsed$paper_cat, collapsed$germline_reported_or_reportable_now)
plot(table(collapsed$paper_cat, collapsed$germline_reported_or_reportable_now))


gl_vars = gl_vars[gl_vars$include_in_gl_stats == 1 & gl_vars$germline_reported_or_reportable_now==1,]

var_counts = table(gl_vars$germline_presented_to_mdt[gl_vars$germline_reported_or_reportable_now==1])
var_counts = sort(var_counts, decreasing = TRUE)
gl_vars$plot_idx = sapply(gl_vars$germline_presented_to_mdt, function(rep) which(names(var_counts) == rep))

gl_vars = gl_vars[order(gl_vars$plot_idx),]
gl_vars$germline_presented_to_mdt_gene = gsub("_.*", "", gl_vars$germline_presented_to_mdt)
gl_vars$germline_presented_to_mdt = ordered(gl_vars$germline_presented_to_mdt, levels = names(var_counts))
gl_vars$germline_presented_to_mdt_gene = ordered(gl_vars$germline_presented_to_mdt_gene, levels = gsub("_.*", "", names(var_counts)))

library(ggplot2)
library(ggsci)

ggsave("gl_fig2a.svg", ggplot(gl_vars[gl_vars$germline_reported_or_reportable_now==1,], aes(x = germline_presented_to_mdt_gene, y = ..count.., fill = paper_cat)) + geom_bar() + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 45)) + labs(x = "Reported germline gene", y = "Frequency", fill = "Cancer group"))


ggsave("gl_fig2b.svg", ggplot(gl_vars[gl_vars$germline_reported_or_reportable_now==1,], aes(x = paper_cat, y = ..count.., fill = paper_cat)) + geom_bar() + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 45)) + labs(x = "Cancer group", y = "Frequency of reported germline findings", fill = "Cancer group"))

ggsave("gl_fig2b_pct.svg", ggplot(gl_vars[gl_vars$germline_reported_or_reportable_now==1,], aes(x = paper_cat, y = ..count../sum(gl_vars$germline_reported_or_reportable_now==1)*100, fill = paper_cat)) + geom_bar() + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 45)) + labs(x = "Cancer group", y = "Frequency of reported germline findings", fill = "Cancer group"))


ggplot(collapsed, aes(x = paper_cat, alpha = germline_reported_or_reportable_now, fill = paper_cat, col = paper_cat)) + geom_bar() + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 45)) + labs(x = "Cancer group", y = "Frequency of patients with reported germline findings", fill = "Cancer group")

ggsave("gl_fig2b_collapsed_pct.svg", ggplot(collapsed, aes(x = paper_cat, alpha = germline_reported_or_reportable_now, fill = paper_cat, col = paper_cat)) + geom_bar(position = position_fill()) + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 45)) + labs(x = "Cancer group", y = "Frequency of patients with reported germline findings", fill = "Cancer group", col = "Cancer group", alpha = "Finding made") + ylim(0, 0.25))
