# Prepare to GSEA

setwd("C:/Users/François/Dropbox/FunctionalGenomics project")

rnaseq_backup <- read.table('HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep='\t', fill=TRUE, header = TRUE)
rnaseq <- rnaseq_backup

row.names (rnaseq) <- rnaseq[,1] # 1st col = gene ID
rnaseq <- rnaseq[,-1] # removes 1st col (now rownames)
rnaseq <- rnaseq[-1,] # first row is just "normalized counts"

# Only keep patients for which I have DNAm tumor age
rna_patients <- colnames (rnaseq)
# !! there are normal tissues in there
tissue <- substring (rna_patients, 14, 15) # 01 = primary tumor; 06 = metastase; 11 = normal tissue

meta_inds <- which (tissue == '06') # only 2 metastases (indexes 489; 546)
normal_inds <- which (tissue == '11') # 44 normal tissues rnaseq
tumor_inds <- which (tissue == '01') # 520 tumor tissues rnaseq

# sum of the lengths of the indexes should be total number of samples
length (meta_inds) + length (normal_inds) + length (tumor_inds) == length (rna_patients) # TRUE. OK!

# separate rnaseq from normal and tumor tissues
rnaseq_normal <- rnaseq[,normal_inds] # 44 cols, ok.
rnaseq_tumor <- rnaseq[,tumor_inds] # 520 cols, ok.

# I'm going to work on tumor tissue. Not normal.
# but normal can be interested. Just not the question.
tumor_rna_patients <- colnames (rnaseq_tumor)
tumor_rna_patients <- substring (tumor_rna_patients, 9, 12)
tumor_rna_patients <- tolower (tumor_rna_patients) # patients for who I have RNAseq

age_patients <- Clinical$patient
age_patients <- as.character(age_patients) # patients for who I have tumor DNAm age

colnames (rnaseq_tumor) <- tumor_rna_patients # Each col of rnaseq = one patient

# duplicates were actually normal tissues.
# no more duplicates

# remove all patients not in clinical.
# If I don't have tumor DNAm age, I cannot use it.
length (which (tumor_rna_patients %in% age_patients) == TRUE) # 342 (max = 521) rna patients are in age_patients
length (which (age_patients %in% tumor_rna_patients) == TRUE) # 342 (max = 348) age patients are in rna_patients
# so there are some patients too much in both
length(intersect (age_patients, tumor_rna_patients)) # 342 intersect = those who are in both

patients_keep <- intersect (age_patients, tumor_rna_patients) # these are the patients I should have in both at the end

rna_keep_inds <- which (tumor_rna_patients %in% patients_keep == TRUE) # indexes of the patients I want from rnaseq
rnaseq_tumor <- rnaseq_tumor [,rna_keep_inds] # take only these cols = 342 cols now

# keep the same in clinical
clinical_keep_inds <- which (age_patients %in% patients_keep == TRUE) # indexes (rows) of the patients in clinical who are in patients_keep
Clinical_rna <- Clinical [clinical_keep_inds,] # keep only these rows # 342 rows now

# Check this was done properly.
# They should be exactly the same patients (not necessarily in the same order for now)
Clinical_rna$patient %in% colnames (rnaseq_tumor) # all TRUE
length (intersect (Clinical_rna$patient, colnames (rnaseq_tumor))) # and intersection is 342. OK.

# + They need to be in the same order too
clinical_patients <- as.character(Clinical_rna$patient)
tumor_rna_patients <- colnames (rnaseq_tumor)

clinical_patients == tumor_rna_patients # All FALSE = they are not in the right order

rnaseq_tumor <- rnaseq_tumor[,order(names(rnaseq_tumor))] # Order rnaseq alphabetically according to col names

tumor_rna_patients <- colnames (rnaseq_tumor)

clinical_patients == tumor_rna_patients # All TRUE = they are in the same order. Good!
# same patients in same order in both

# Check that no gene is duplicated, all unique
genes <- rownames(rnaseq_tumor)
which (duplicated(genes)) == TRUE # no duplicated gene, OK!

# Use only one from ENTREZID / GENESYMBOL
# I'll use ENTREZID
# because I'm scared first ones which have ? as ENTREZID might cause problems later

genes <- rownames (rnaseq_tumor)
genes <- gsub ("\\|.*$", "", genes) # keep only the ENTREZID part (uses regexp)
which (genes == '?') # > first 29 genes are ?

rnaseq_tumor <- rnaseq_tumor [-c(1:29),] # Remove first 29 rows in rnaseq
genes <- genes [-(which(genes=='?'))] # remove same in genes
# 20502 rows in rnaseq & 20502 genes in genes. And first ones are the same. Looks good.

# Actually one duplicate in genes, remove it from genes and rnaseq.
# Not spotted because it has two different GENESYMBOL
# Row 16273
genes <- genes [-16273]
rnaseq_tumor <- rnaseq_tumor [-16273,]
# rnaseq has 20501 rows; length genes is 20501. Good!

rownames(rnaseq_tumor) <- genes # change the rownames by the genes (strings treated)

# I need to extract the DNAm tumor ages from Clinical
# they will match the patients (as it is same patients in the same order)
tumor_age_rna <- as.numeric (as.character(Clinical_rna$cancer_age))

# Export DNAm tumor age for patients
#  write.table (tumor_age_rna, "C:/Users/François/Dropbox/FunctionalGenomics project/tumor_age.txt", sep = '\t')

# Export RNAseq data
#  write.table (rnaseq_tumor, "C:/Users/François/Dropbox/FunctionalGenomics project/rnaseq_tumor.txt", sep = '\t')

# do the same with tumour acceleration 2
# just need to export tumour acceleration 2

tumor_acc2 <- as.numeric (as.character(Clinical_rna$tumoracc2))
# write.table (tumor_acc2, "C:/Users/François/Dropbox/FunctionalGenomics project/tumoracc2.txt", sep = '\t')

########################################################################################################
# RUN SAMR

# source("https://bioconductor.org/biocLite.R")
# biocLite("impute")
# install.packages ('samr')

library (samr)
# maybe convert into matrix
rnaseq_tumor_mat <- as.matrix(rnaseq_tumor)
# columns are characters, weird. They should be numerics
rnaseq_tumor_mat <- apply (rnaseq_tumor_mat, 2, as.numeric)

# I think problem because counts are not integers
rnaseq_tumor_mat <- round (rnaseq_tumor_mat)

# try again
# samfit_age <- SAMseq (rnaseq_tumor_mat, tumor_age_rna, resp.type = 'Quantitative')

samfit_acc <- SAMseq (rnaseq_tumor_mat, tumor_acc2, resp.type = 'Quantitative')

# 
# # Import non-normalized data
# #  rnaseq_raw <- read.table('HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt', sep='\t', fill=TRUE)
# # takes like 6 min. Don't do it again.
# 
# # rnaseq_raw_backup <- rnaseq_raw
# 
# # recover like this
# rnaseq_raw <- rnaseq_raw_backup
# 
# # 20533 rows vs. 1699 cols
# # first column is gene
# ncol (rnaseq_raw) / 3 # 566.333, matches number of patients in rnaseq
# 
# # first column are genes, put that as row name
# rownames(rnaseq_raw) <- rnaseq_raw[,1]
# # and delete first column then
# rnaseq_raw <- rnaseq_raw [,-1]
# 
# 
# # only keep raw count columns = first columns for each sample = one column each 3
# # columns to keep are # 1, 4, 
# 
# # I THINK IT IS WRONG THAT I NEED NON_NORMALIZED
# 
# 
