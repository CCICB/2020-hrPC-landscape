# Whole genome, transcriptome and methylome profiling enhances actionable target discovery in high-risk paediatric cancer 
Wong, Mayoh, Lau et al, 2020, Nature Medicine

This repository contains scripts and figure generation code used for this publication. They are shared as is, with no warranty.

# Genomic Methods
Very extensive genomic methods are available in the online publication.

# Data
The raw and processed WGS, RNAseq and methylation data generated by this study are available from the European Genome Archive, accession number EGAS00001004572. As of 12/08/2020, data are still being migrated to EGA. ETA is 4-6 weeks from now.

# Analysis Pipeline
The analytical pipeline and docker containers will soon be made available through the AusBioCommons GitHub repository. To be updated when we have a link. ETA is by 31/12/2020.

# Scripts
Additional scripts mentioned in the publication are in the scripts folder here, or:

* https://bitbucket.org/cciacb/cci-strelka-filter/src/master/

# Software versions used
We are so grateful to all these open source software developers for developing and maintaining these critical software.

* ANNOVAR (v20190929)
* Arriba (v1.1.0)
* bedtools (v2.28.0)
* Branchpointer (v1.3.1)
* BWA-MEM (v0.17.10-r789)
* CADD (v1.3)
* dbNSFP (v2.9)
* dbscSNV (v1.1)
* deconstructSigs (v1.8.0)
* FastQC (v0.11.5)
* FATHMM (via dbNSFP v2.9)
* GATK GenotypeVCFs (v3.3) 
* GATK HaplotypeCaller (v3.3 for WGS and v3.6 for RNAseq) 
* GATK Indel Realignment (v3.3) 
* GATK ReassignOneMappingQuality (v3.6)
* GATK SplitNCigarReads (v3.6) 
* GATK VQSR (v3.3) 
* GEMINI (v0.11.0)
* ggplot2 (v3.3.2)
* ggsashimi (v0.4.0) 
* GRIDSS (v2.7.2)
* IGV (v2.6.2)
* Introme (v0.5.1)
* JAFFA (v1.09)
* LINX (v1.7)
* MetaLR (v1.0)93 
* MetaSVM (v1.0)
* MMSplice (v2.1.0)
* MNP Classifier (https://www.molecularneuropathology.org/mnp, versions may have changed over time)
* NBR (see NBR folder in this repo)
* Novosort (v1.03.01)
* Polyphen2 (v2.2.2)
* PROVEAN (v1.1)
* PURPLE (v2.39)
* R (v3.5.3)
* Refynr (v1.17.8)
* RNA VAF estimator (see scripts folder in this repo)
* RSEM (v1.2.31)
* RStudio (v1.2.1335)
* SAMTools (v1.3.1)
* Seave (https://seave.bio; updated throughout the project)
* SIFT (v5.0.2)
* SnpEff (v4_3t)
* SPIDEX (v1.0)
* SpliceAI (v1.3.1)
* STAR (v2.5)
* STAR-Fusion (v1.3.1)
* Strelka (v2.0.17)
* Strelka filter (https://bitbucket.org/cciacb/cci-strelka-filter/src/master/)
* Variant Effect Predictor, VEP (v87)

# Databases used
* 1000 genomes (https://www.internationalgenome.org/data)
* COSMIC (https://cancer.sanger.ac.uk/cosmic)
* Cancer Gene Census (https://cancer.sanger.ac.uk/census)
* ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/)
* dbNSFP (https://sites.google.com/site/jpopgen/dbNSFP)
* dbscSNV (http://www.liulab.science/dbscsnv.html)
* ESP (https://evs.gs.washington.edu/EVS/)
* ExAC (http://exac.broadinstitute.org/)
* GIAB (https://jimb.stanford.edu/giab-resources)
* gnomAD (https://gnomad.broadinstitute.org/)
* MGRB (https://sgc.garvan.org.au/)
* Pecan (https://pecan.stjude.cloud/)
* Platinum Genomes (https://github.com/Illumina/PlatinumGenomes)
