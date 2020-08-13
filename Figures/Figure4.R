library(openxlsx)

# Load data and mark non-reported vars
gl_vars = read.xlsx("glvars_redacted.xlsx")
gl_vars$cancer_category = factor(gl_vars$cancer_category, levels = c("CNS", "HM", "NBL", "Sarcoma", "Solid tumour"))
gl_vars$reported = !is.na(gl_vars$germline_allele) & !grepl("^ +", gl_vars$germline_allele)

# Summarise the fraction of patients with germline findings by cancer type
library(plyr)
patient_summary = ddply(gl_vars, .(cancer_category), function(d) {
	patients_without_finding = length(unique(d$sample_id[d$reported == FALSE]))
	patients_with_finding = length(unique(d$sample_id[d$reported == TRUE]))
	data.frame(cancer_category = d$cancer_category[[1]], finding_rate = patients_with_finding / (patients_with_finding + patients_without_finding))})

# Code to order plots by decreasing gene allele count
gl_vars$gene = gsub("_.*", "", gl_vars$reported_germline)
gene_allele_counts = table(gl_vars$gene[gl_vars$reported])
gene_allele_counts = sort(gene_allele_counts, decreasing = TRUE)
gl_vars$plot_idx = sapply(gl_vars$gene, function(gene) if (is.na(gene)) { NA } else { which(names(gene_allele_counts) == gene) })
gl_vars$plot_idx[is.na(gl_vars$plot_idx)] = max(gl_vars$plot_idx, na.rm = TRUE) + 1

gl_vars = gl_vars[order(gl_vars$plot_idx),]
gl_vars$gene = ordered(gl_vars$gene, levels = names(gene_allele_counts))

library(ggplot2)
library(ggsci)

# Figure 4b
ggplot(patient_summary, aes(x = cancer_category, fill = cancer_category, col = cancer_category, y = finding_rate)) + geom_bar(stat = "identity") + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 90)) + labs(x = "Cancer group", y = "Frequency of patients with reported germline findings", fill = "Cancer group", col = "Cancer group", alpha = "Finding made") + ylim(0, 0.25)

# Figure 4c
ggplot(gl_vars[gl_vars$reported,], aes(x = gene, y = ..count.., fill = cancer_category)) + geom_bar() + 
    theme_bw() + scale_color_npg() + scale_fill_npg() + theme(axis.text.x =  element_text(angle = 90)) + labs(x = "Reported germline gene", y = "Frequency", fill = "Cancer group")
