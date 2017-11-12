# number of DNA copy number breakpoints

setwd("C:/Users/François/Dropbox/FunctionalGenomics project/")

genome <- read.table('genome_wide.txt', sep='\t', fill=TRUE, header = TRUE)
options (scipen = 10000)
patient1 <- head (genome, 280)

plot (Segment_Mean ~ Start, data = patient1, xlab = 'Genome position (start)', ylab = 'Log2(Copy Number/2)',
       pch = 20, col = '#343335')
 abline (h = 0, lty = 'dashed', lwd = 2, col = 'red')

codes <- genome$Sample
tissue <- substring (codes, 14, 15)
blood_inds <- which (tissue == '10') # 30927 # blood normal
tumour_inds <- which (tissue == '01') # 74970 # primary tumour
normal_inds <- which (tissue == '11') # 4011 # normal tissue matched
meta_inds <- which (tissue == '06') # 381 # metastasis
# genome is 110289. Total here is 110289. # OK

genome_blood <- genome [blood_inds,]

Genome <- genome [tumour_inds,]
rownames (Genome) <- NULL

# patients
patient <- substring (codes, 9, 12)

# plot (Segment_Mean ~ Start, data = Genome, xlab = 'Genome position (start)', ylab = 'Log2(Copy Number/2)',
#       pch = 20, col = '#343335')
# abline (h = 0, lty = 'dashed', lwd = 2, col = 'red')

# grep ('AA8J', codes, value = FALSE)


takeDataPerPatient <- function (patients, genome) {
  for (p in patients) {
    print (p)
    inds <- grep (p, Genome$Sample, value = FALSE)
    PerPatient [[length(PerPatient) + 1]] <<- Genome [inds,]
  }
}

# you can take only patients that are in Clinical
patient_clinical <- toupper (Clinical$patient)
patient <- substring (Genome$Sample, 9, 12)
common <- intersect (patient, patient_clinical)
rowstokeep <- which (patient %in% common) # 345 (all from Clinical are there)
Genome <- Genome [rowstokeep,]
rownames (Genome) <- NULL

# order genome according to patient code
patient <- substring (Genome$Sample, 9, 12)
Genome <- cbind (Genome, patient)
Genome <- Genome[order(Genome$patient),]
common <- common[order(common)]


PerPatient <- list()
takeDataPerPatient(common, Genome)

# extract the CNV and compute the mean for each

cnvs <- c()
for (PP in PerPatient) {
  cnv <- mean(PP[,6], na.rm = TRUE)
  print (cnv)
  cnvs <- c(cnvs, cnv)
}

# means are in cnvs
# for patient in clinical order

# now take the data from Clinical and you'll have your datapoint

patient_clinical <- toupper (Clinical$patient)
patient_genome <- unique(Genome$patient)
patient_genome %in% patient_clinical # all TRUE
inds <- patient_clinical %in% patient_genome # 3 F as expected. 3 too much in clinical
which (inds == FALSE) #13, 40, 84
tumourage <- Clinical$cancer_age
tumourage <- tumourage [-c(13,40,84)] # 345 now!
# same in both and they are patient matched alphabetically
#merge them

tumouracc <- Clinical$tumoracc2
tumouracc <- tumouracc [-c(13,40,84)]

cnv_age <- as.data.frame (cbind (tumourage,tumouracc,cnvs))

# sort by age and you're done
cnv_age <- cnv_age [order(cnv_age$tumourage),]


plot (cnvs ~ tumourage, data = cnv_age, type = 'p', pch = 20,
      xlab = 'Tumour DNAm age (years)', ylab = 'Average CNV', col = '#871F07', xaxs = 'i', xlim = c(0, 90))
abline (h = 0, lty = 'dashed', lwd = 2, col = '#343335')
cor.test (cnv_age$tumourage, cnv_age$cnvs, use = 'pairwise.complete.obs', method = 'spearman')

plot (cnvs ~ tumouracc, data = cnv_age, type = 'p', pch = 20,
      xlab = 'Tumour age acceleration', ylab = 'Average CNV', col = '#025D8C', xaxs = 'i', xlim = c(0, 2.0))
abline (h = 0, lty = 'dashed', lwd = 2, col = '#343335')
cor.test (cnv_age$tumouracc, cnv_age$cnvs, use = 'pairwise.complete.obs', method = 'spearman')
