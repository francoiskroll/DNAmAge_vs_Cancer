# FUNCTIONAL GENOMICS PROJECT

setwd("C:/Users/François/Dropbox/FunctionalGenomics project")

load ('HNSC-8.Rda')
meth <- as.data.frame (dm.age.subset)

# Import RNAseq data
# rnaseq_backup <- read.table('HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt', sep='\t', fill=TRUE, header = TRUE)
# rnaseq <- rnaseq_raw
# 
# row.names (rnaseq) <- rnaseq[,1]
# rnaseq <- rnaseq[,-1]
# rnaseq <- rnaseq[-1,]

# Import clinical data
clinical1 <- read.table('All_CDEs.txt', sep='\t', fill=TRUE, header = TRUE)
clinical2 <- read.table('HNSC.clin.merged.picked.txt', sep='\t', fill=TRUE, header = TRUE)
codenames <- read.table('stage3_params_clin_selection_HNSC.tsv', sep='\t', fill=TRUE, header = TRUE)

# Input are: rnaseq, clinical1, clinical2, meth. All dataframes

# From meth, separate normal and cancer tissues
meth2 <- meth

clinical2 <- t(clinical2)
colnames (clinical2) <- clinical2[1,]
clinical2 <- clinical2[-1,]
clinical2 <- as.data.frame (clinical2)
patient <- row.names(clinical2)
patient <- substring (patient, 9) # takes only patient code. clinical2 is only clinical annotations, no tissue code
clinical2 <- cbind(patient, clinical2)

patient <- tolower (rownames(meth2))
patient <- substring (patient, 9, 12) # takes only patient code. !! Tissue code here

tissue <- substring (rownames(meth2), 14, 15) # takes only tissue code. 
tissue <- as.numeric (tissue)

meth2 <- cbind (meth2, tissue, patient) # 3 types: 1, 6, 11. 1 = primary tumor, 6 = metastase; 11 = normal tissue
# you can just cbind them to patient, order won't have changed

# Patient a4ca has no age in clinical2. Remove it.
meth2 <- meth2 [- which (meth2$patient == 'a4ca'),]
clinical2 <- clinical2 [- which (clinical2$patient == 'a4ca'),]

unique_patients <- unique (meth2$patient)
clinical2 <- clinical2 [clinical2$patient %in% unique_patients,]

cancer_inds <- which (meth2$tissue == 1)# One metastase, did not pick it
normal_inds <- which (meth2$tissue == 11)

(length (cancer_inds) + length(normal_inds)) + 1 == nrow (meth2) # correct

cancer <- meth2 [cancer_inds,]
cancer <- cancer [,-2]
colnames (cancer) [1] <- c('cancer_age')

normal <- meth2 [normal_inds,]
normal <- normal [,-2]
colnames (normal) [1] <- c('normal_age')

Clinical <- merge (clinical2, normal, by = 'patient', all = TRUE) 
Clinical <- merge (Clinical, cancer, by = 'patient', all = TRUE)
# !!! conversion of a factor. Be extra careful.
Clinical$years_to_birth <- as.numeric (as.character(Clinical$years_to_birth))
Clinical$number_pack_years_smoked <- as.numeric (as.character(Clinical$number_pack_years_smoked))
Clinical$number_of_lymph_nodes <- as.numeric (as.character(Clinical$number_of_lymph_nodes))
Clinical$days_to_last_followup <- as.numeric (as.character(Clinical$days_to_last_followup))
Clinical$days_to_death <- as.numeric (as.character(Clinical$days_to_death))

# clean up clinical1 and add a couple of interesting variables to Clinical
clinical1 <- t(clinical1)
colnames (clinical1) <- clinical1[1,]
clinical1 <- clinical1[-1,]
clinical1 <- as.data.frame (clinical1)

patient1 <- row.names(clinical1)
patient1 <- substring (patient1, 9)
clinical1 <- cbind(patient1, clinical1)
colnames (clinical1) [1] <- 'patient'
clinical1 <- clinical1 [- which (clinical1$patient == 'a4ca'),]

# Adding some data from clinical1 to Clinical
Clinical <- merge (Clinical, clinical1 [,c('patient', 'alcohol_history_documented')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'amount_of_alcohol_consumption_per_day')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'clinical_m')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'clinical_n')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'clinical_stage')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'clinical_t')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'days_to_completion_of_curative_tx')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'disease_after_curative_tx')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'followup_treatment_success')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'frequency_of_alcohol_consumption')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'history_of_neoadjuvant_treatment')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'hpv_call_1')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'hpv_status')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'hpv_status_by_ish_testing')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'hpv_status_by_p16_testing')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'icd_10')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'icd_o_3_histology')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'icd_o_3_site')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'laterality')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'lymph_node_examined_count')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'lymphovascular_invasion_present')], by = 'patient')
Clinical <- merge (Clinical, clinical1 [,c('patient', 'margin_status')], by = 'patient')

# Conversion of factor
Clinical$amount_of_alcohol_consumption_per_day <- as.numeric (as.character(Clinical$amount_of_alcohol_consumption_per_day))
Clinical$days_to_completion_of_curative_tx <- as.numeric (as.character(Clinical$days_to_completion_of_curative_tx))
Clinical$frequency_of_alcohol_consumption <- as.numeric (as.character(Clinical$frequency_of_alcohol_consumption))
Clinical$lymph_node_examined_count <- as.numeric (as.character(Clinical$lymph_node_examined_count))
Clinical$date_of_initial_pathologic_diagnosis <- as.numeric (as.character(Clinical$date_of_initial_pathologic_diagnosis))
Clinical$year_of_tobacco_smoking_onset <- as.numeric (as.character(Clinical$year_of_tobacco_smoking_onset))

########################################################################################################
# QUESTION 1

# Spearman checks for monotonic correlation (Pearson is only linear).
# Also, distributions do not have to be normal.
par (bg = NA)
# Q1.1: Chronological age vs Tumor DNAm age
plot (Clinical$years_to_birth, Clinical$cancer_age, type = 'p', pch = 20, 
      xlab = 'Chronological age (years)', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 90), bg = NA)

cor (Clinical$years_to_birth, Clinical$cancer_age, use = 'pairwise.complete.obs', method = 'spearman') # 0.21 - low correlation
cor (Clinical$years_to_birth, Clinical$cancer_age, use = 'pairwise.complete.obs', method = 'pearson')
cor.test (Clinical$years_to_birth, Clinical$cancer_age, use = 'pairwise.complete.obs', method = 'spearman', exact = FALSE) # very low p-value

# Q1.2: Chronological age vs Normal tissue DNAm age
plot (Clinical$years_to_birth, Clinical$normal_age, type = 'p', pch = 20, 
      xlab = 'Chronological age (years)', ylab = 'Normal tissue DNAm age (years)', col = '#043264',
      yaxs = 'i', ylim = c(0, 90), bg = NA)
cor (Clinical$years_to_birth, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'spearman') # 0.67 - correlate
cor (Clinical$years_to_birth, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'pearson')
cor.test (Clinical$years_to_birth, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'spearman', exact = FALSE) # low p-value
cor.test (Clinical$years_to_birth, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'spearman')
# Q1.3: Tumor DNAm age vs Normal tissue DNAm age
plot (Clinical$normal_age, Clinical$cancer_age, type = 'p', pch = 20, 
      xlab = 'Normal tissue DNAm age (years)', ylab = 'Tumour DNAm age (years)', col = '#506620',
      yaxs = 'i', ylim = c(0, 90), bg = NA)
cor (Clinical$cancer_age, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'spearman') # 0.23 - do not correlate
cor (Clinical$cancer_age, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'pearson')
cor.test (Clinical$cancer_age, Clinical$normal_age, use = 'pairwise.complete.obs', method = 'spearman') # 0.1 p-value, do not trust.

# Conclusions are the same if Pearson used.

# Q1.4: Two definitions tumor age acceleration.
# Not really acceleration but:
  # def1: ratio between tumor age and normal tissue age. Eg. tumor is twice younger than normal tissue.
  # def2: ratio between tumor age and chronological age. Eg. tumor is twice younger than the patient.
tumoracc1 <- Clinical$cancer_age / Clinical$normal_age
tumoracc2 <- Clinical$cancer_age / Clinical$years_to_birth
Clinical <- cbind(Clinical, tumoracc1)
Clinical <- cbind(Clinical, tumoracc2)

plot (Clinical$tumoracc1, Clinical$tumoracc2, type = 'p', pch = 20,
      xlab = 'Tumour acceleration (1)', ylab = 'Tumour acceleration (2)', col = '#FF4040',
      yaxs = 'i', xaxs = 'i', ylim = c(0, 2), xlim = c(0, 2.2), bg = NA)

cor (Clinical$tumoracc1, Clinical$tumoracc2, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$tumoracc1, Clinical$tumoracc2, use = 'pairwise.complete.obs', method = 'spearman', exact = TRUE)

# what about their distributions?
plot(density(Clinical$tumoracc1, na.rm = TRUE), 
     type = 'l', lwd = 3, col = '#506620', main = '',
     yaxs = 'i', xlab = 'Tumour age acceleration', ylab = 'Density',
     ylim = c(0.0, 2.0), bty = 'l')
lines (density(Clinical$tumoracc2, na.rm = TRUE), type = 'l', lwd = 3, col = '#025D8C')

ks.test (Clinical$tumoracc1[which (! is.na(Clinical$tumoracc1))], Clinical$tumoracc2) # low pvalue normally means not likely from the same distribution

########################################################################################################

# QUESTION 2

# use ...._picked, not All_CDEs
# You can add some from clinical1, we'll see after.

# Before using tests, check distribution normality

# Normality of tumor age?
plot(density(Clinical$cancer_age), type = 'l', lwd = 3, col = '#871F07', main = '',
     yaxs = 'i', xlab = 'Tumour DNAm age', ylab = 'Density', bty = 'l')

plot (density(Clinical$tumoracc2), type = 'l', lwd = 3, col = '#025D8C', main = '',
       yaxs = 'i', xlab = 'Tumour age acceleration', ylab = 'Density', bty = 'l')

qqnorm(Clinical$cancer_age, main = '', col = '#871F07', bty = 'l',
       xlab = 'Theoretical quantiles', ylab = 'Sample quantiles', pch = 21)
qqline(Clinical$cancer_age, col = '#6C686E', lwd = 2)

qqnorm(Clinical$tumoracc2, main = '', col = '#025D8C', bty = 'l',
       xlab = 'Theoretical quantiles', ylab = 'Sample quantiles', pch = 21)
qqline(Clinical$tumoracc2, col = '#6C686E', lwd = 2)

shapiro.test(Clinical$cancer_age) # do not trust that too much. Sample is big, expected to be low.
# low means safe to reject H0 = distribution is normal
# -> Basically it is normal enough so you can run tests.

# Normality of tumor acceleration 1?
plot(density(na.omit(Clinical$tumoracc1)))
qqnorm(Clinical$tumoracc1);qqline(Clinical$tumoracc1, col = 2)
shapiro.test(Clinical$tumoracc1) # same here but QQ plot looks fine

# Normality of tumor acceleration 2? # USE THIS ONE
plot(density(na.omit(Clinical$tumoracc2)))
qqnorm(Clinical$tumoracc1);qqline(Clinical$tumoracc2, col = 2) # maybe not be ok here, careful
shapiro.test(Clinical$tumoracc2)  # and p value extremely low!

# Normality of normal tissue age?
plot(density(na.omit(Clinical$normal_age)))
qqnorm(Clinical$normal_age);qqline(Clinical$normal_age, col = 2)
shapiro.test (Clinical$normal_age)

# Normality of chrono age?
plot(density(na.omit(Clinical$years_to_birth)))
qqnorm(Clinical$years_to_birth);qqline(Clinical$years_to_birth, col = 2)
shapiro.test (Clinical$years_to_birth)

# BINARY VARIABLES

# vs. vital status
boxplot (cancer_age ~ vital_status, data = Clinical, 
         names = c('Alive', 'Dead'), xlab = 'Vital status', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

boxplot (tumoracc2 ~ vital_status, data = Clinical, 
         names = c('Alive', 'Dead'), xlab = 'Vital status', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0.0, 2.0))

t.test (cancer_age ~ vital_status, data = Clinical, alternative = 'two.sided')
wilcox.test(cancer_age ~ vital_status, data = Clinical, alternative = 'two.sided')

t.test (tumoracc1 ~ vital_status, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc1 ~ vital_status, data = Clinical, alternative = 'two.sided')

t.test (tumoracc2 ~ vital_status, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc2 ~ vital_status, data = Clinical, alternative = 'two.sided')

t.test (normal_age ~ vital_status, data = Clinical, alternative = 'two.sided')
wilcox.test (normal_age ~ vital_status, data = Clinical, alternative = 'two.sided')

# tumor age vs gender
boxplot (cancer_age ~ gender, data = Clinical, 
         names = c('Female', 'Male'), xlab = 'Gender', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

boxplot (tumoracc2 ~ gender, data = Clinical, 
         names = c('Female', 'Male'), xlab = 'Gender', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0.0, 2.0))

boxplot (normal_age ~ gender, data = Clinical, 
         names = c('Female', 'Male'), xlab = 'Gender', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

t.test (cancer_age ~ gender, data = Clinical, alternative = 'two.sided') # not significant
wilcox.test (cancer_age ~ gender, data = Clinical, alternative = 'two.sided')

t.test (tumoracc1 ~ gender, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc1 ~ gender, data = Clinical, alternative = 'two.sided')

t.test (tumoracc2 ~ gender, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc2 ~ gender, data = Clinical, alternative = 'two.sided')

wilcox.test (normal_age ~ gender, data = Clinical, alternative = 'two.sided')
wilcox.test (normal_age ~ gender, data = Clinical, alternative = 'two.sided')

# tumor age vs radiation therapy
boxplot (cancer_age ~ radiation_therapy, data = Clinical, 
         names = c('No', 'Yes'), xlab = 'Radiation therapy', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

boxplot (tumoracc2 ~ radiation_therapy, data = Clinical, 
         names = c('No', 'Yes'), xlab = 'Radiation therapy', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0, 2))

t.test (cancer_age ~ radiation_therapy, data = Clinical, alternative = 'two.sided')
wilcox.test (cancer_age ~ radiation_therapy, data = Clinical, alternative = 'two.sided')

t.test (tumoracc1 ~ radiation_therapy, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc1 ~ radiation_therapy, data = Clinical, alternative = 'two.sided')

t.test (tumoracc2 ~ radiation_therapy, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc2 ~ radiation_therapy, data = Clinical, alternative = 'two.sided')

t.test (normal_age ~ radiation_therapy, data = Clinical, alternative = 'two.sided')
wilcox.test (normal_age ~ radiation_therapy, data = Clinical, alternative = 'two.sided')

# tumor age vs hipanic or latino
boxplot (cancer_age ~ ethnicity, data = Clinical)

t.test (cancer_age ~ ethnicity, data = Clinical, alternative = 'two.sided')
wilcox.test (cancer_age ~ ethnicity, data = Clinical, alternative = 'two.sided')

t.test (tumoracc1 ~ ethnicity, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc1 ~ ethnicity, data = Clinical, alternative = 'two.sided')

t.test (tumoracc2 ~ ethnicity, data = Clinical, alternative = 'two.sided')
wilcox.test (tumoracc2 ~ ethnicity, data = Clinical, alternative = 'two.sided')

t.test (normal_age ~ ethnicity, data = Clinical, alternative = 'two.sided')
wilcox.test (normal_age ~ ethnicity, data = Clinical, alternative = 'two.sided')

# CATEGORICAL VARIABLES
# more than 2 categories

# vs pathologic_stage
boxplot (cancer_age ~ pathologic_stage, data = Clinical, xlab = 'Pathologic stage', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100)) 
# careful, only 1 stage ivc
# replace the only IVc with NA. It is row 344, patient name a8z8
Clinical$pathologic_stage [344] <- NA
#  Clinical$pathologic_stage <- droplevels (Clinical$pathologic_stage) # drops the IVc level which is not there anymore

# bartlett.test (cancer_age ~ vital_status + , data = Clinical)
anova_pathologic_stage <- aov (cancer_age ~ pathologic_stage, data = Clinical)
summary (anova_pathologic_stage)
kruskal.test (cancer_age ~ pathologic_stage, data = Clinical)

boxplot (tumoracc1 ~ pathologic_stage, data = Clinical)
boxplot (tumoracc2 ~ pathologic_stage, data = Clinical, xlab = 'Pathologic stage', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0, 2), bg = NA) # nothing interesting

anova_pathologic_stage_acc <- aov (tumoracc2 ~ pathologic_stage, data = Clinical)
summary (anova_pathologic_stage_acc)

boxplot (normal_age ~ pathologic_stage, data = Clinical)

# vs pathology_T_stage
# Clinical$pathology_T_stage <- droplevels (Clinical$pathology_T_stage)
boxplot (cancer_age ~ pathology_T_stage, data = Clinical, xlab = 'T pathologic stage', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100)) 
anova_pathology_T_stage <- aov (cancer_age ~ pathology_T_stage, data = Clinical)
summary (anova_pathology_T_stage)

boxplot (tumoracc2 ~ pathology_T_stage, data = Clinical, xlab = 'T pathologic stage', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0, 2), bg = NA) 

anova_pathology_T_stage_acc <- aov (tumoracc2 ~ pathology_T_stage, data = Clinical)
summary (anova_pathology_T_stage_acc)

# vs pathology_N_stage
boxplot (cancer_age ~ pathology_N_stage, data = Clinical, xlab = 'N pathologic stage', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

anova_pathology_N_stage <- aov (cancer_age ~ pathology_N_stage, data = Clinical)
summary (anova_pathology_N_stage)

boxplot (tumoracc2 ~ pathology_N_stage, data = Clinical, xlab = 'N pathologic stage', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0, 2), bg = NA) 
anova_pathology_N_stage_acc <- aov (tumoracc2 ~ pathology_N_stage, data = Clinical)
summary (anova_pathology_N_stage_acc)

# vs pathology_M_stage
boxplot (cancer_age ~ pathology_M_stage, data = Clinical)

# vs histological_type
boxplot (cancer_age ~ histological_type, data = Clinical)

# vs race
boxplot (cancer_age ~ race, data = Clinical, xlab = 'Race', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100), names = c('native', 'asian', 'black', 'white'))

anova_race <- aov (cancer_age ~ race, data = Clinical)
summary (anova_race)

kruskal.test (cancer_age ~ race, data = Clinical)

boxplot (tumoracc1 ~ race, data = Clinical)

boxplot (tumoracc2 ~ race, data = Clinical, xlab = 'Race', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0, 2), bg = NA, names = c('native', 'asian', 'black', 'white')) 

anova_acc_race <- aov (tumoracc2 ~ race, data = Clinical)
summary (anova_acc_race)

boxplot (normal_age ~ race, data = Clinical)

# CONTINUOUS VARIABLES

# number of years the patient smoked
years_smoking <- Clinical$date_of_initial_pathologic_diagnosis - Clinical$year_of_tobacco_smoking_onset
Clinical <- cbind (years_smoking, Clinical)

plot (Clinical$years_smoking, Clinical$cancer_age, type = 'p', pch = 20, 
      xlab = 'Number of years smoking', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 90), xlim = c(0, 70), bg = NA, xaxs = 'i')

plot (Clinical$years_smoking, Clinical$tumoracc2, type = 'p', pch = 20, 
      xlab = 'Number of years smoking', ylab = 'Tumour age acceleration', col = '#025D8C',
      yaxs = 'i', ylim = c(0, 2.0), xlim = c(0, 70), bg = NA, xaxs = 'i')

plot (Clinical$years_smoking, Clinical$normal_age, type = 'p', pch = 20, 
      xlab = 'Number of years smoking', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 90), xlim = c(0, 70), bg = NA, xaxs = 'i')

plot (Clinical$normal_age, Clinical$years_smoking, type = 'p', pch = 20, 
      xlab = 'Number of years smoking', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 100), xlim = c(0, 100), bg = NA, xaxs = 'i')

plot (Clinical$years_smoking, Clinical$years_to_birth, type = 'p', pch = 20, 
      xlab = 'Number of years smoking', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 90), xlim = c(0, 70), bg = NA, xaxs = 'i')

cor.test (Clinical$years_to_birth, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'spearman')

cor.test (Clinical$normal_age, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'spearman')

cor.test (Clinical$cancer_age, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'spearman')
cor (Clinical$cancer_age, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'spearman')
cor (Clinical$cancer_age, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'pearson')

cor.test (Clinical$tumoracc2, Clinical$years_smoking, use = 'pairwise.complete.obs', method = 'spearman')
# smoker vs non smoker

smoker_status <- Clinical$number_pack_years_smoked
smoker_status [is.na(smoker_status)] <- 'no'
smoker_status [which (smoker_status != 'no')] <- 'yes'

Clinical <- cbind (smoker_status, Clinical)

# vs. number_pack_years_smoked
plot (Clinical$number_pack_years_smoked, Clinical$cancer_age, type = 'p', pch = 20, 
      xlab = 'Packs smoked/year', ylab = 'Tumour DNAm age (years)', col = '#871F07',
      yaxs = 'i', ylim = c(0, 90), xlim = c(0, 200), bg = NA, xaxs = 'i')

plot (Clinical$number_pack_years_smoked, Clinical$tumoracc2, type = 'p', pch = 20, 
      xlab = 'Packs smoked/year', ylab = 'Tumour age acceleration', col = '#025D8C',
      yaxs = 'i', ylim = c(0, 2.0), xlim = c(0, 200), bg = NA, xaxs = 'i')


cor (Clinical$cancer_age, Clinical$number_pack_years_smoked, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$cancer_age, Clinical$number_pack_years_smoked, use = 'pairwise.complete.obs', method = 'spearman')

cor.test (Clinical$tumoracc2, Clinical$number_pack_years_smoked, use = 'pairwise.complete.obs', method = 'spearman')

# vs. smoking status

boxplot (cancer_age ~ smoker_status, data = Clinical, 
         names = c('No', 'Yes'), xlab = 'Smoker status', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

boxplot (tumoracc2 ~ smoker_status, data = Clinical, 
         names = c('No', 'Yes'), xlab = 'Smoker status', ylab = 'Tumour age acceleration',
         bty = 'l', yaxs = 'i', ylim = c(0.0, 2.0))

t.test (cancer_age ~ smoker_status, data = Clinical, alternative = 'two.sided')
t.test (tumoracc2 ~ smoker_status, data = Clinical, alternative = 'two.sided')

boxplot (normal_age ~ smoker_status, data = Clinical, 
         names = c('No', 'Yes'), xlab = 'Smoker status', ylab = 'Tumour DNAm age (years)',
         bty = 'l', yaxs = 'i', ylim = c(0, 100))

# vs. number of lymph_nodes
plot (Clinical$cancer_age, Clinical$number_of_lymph_nodes, type = 'p', pch = 20, 
      xlab = 'Tumour DNAm age (years)', ylab = 'Number of lymph nodes affected', col = '#871F07',
      yaxs = 'i', ylim = c(-1, 20), xlim = c(0, 100), bg = NA, xaxs = 'i')

cor.test (Clinical$cancer_age, Clinical$number_of_lymph_nodes, use = 'pairwise.complete.obs', method = 'spearman')

plot (Clinical$tumoracc2, Clinical$number_of_lymph_nodes, type = 'p', pch = 20, 
      xlab = 'Tumour age acceleration', ylab = 'Number of lymph nodes affected', col = '#025D8C',
      yaxs = 'i', ylim = c(-1, 20), xlim = c(0, 2.0), bg = NA, xaxs = 'i')
cor.test (Clinical$number_of_lymph_nodes, Clinical$tumoracc2, use = 'pairwise.complete.obs', method = 'spearman')
# ADDITIONAL DATA FROM CLINICAL 1
 
# vs. alcohol_history_documented
boxplot (cancer_age ~ alcohol_history_documented, data = Clinical)
t.test (cancer_age ~ alcohol_history_documented, data = Clinical, alternative = 'two.sided')
wilcox.test (cancer_age ~ alcohol_history_documented, data = Clinical, alternative = 'two.sided')
boxplot (tumoracc1 ~ alcohol_history_documented, data = Clinical)
boxplot (tumoracc2 ~ alcohol_history_documented, data = Clinical)
boxplot (normal_age ~ alcohol_history_documented, data = Clinical)

# vs. amount_of_alcohol_consumption_per_day (continuous)
plot (Clinical$cancer_age, Clinical$amount_of_alcohol_consumption_per_day, type = 'p', pch = 20)
cor.test (Clinical$cancer_age, Clinical$amount_of_alcohol_consumption_per_day, use = 'pairwise.complete.obs', method = 'spearman')
plot (Clinical$tumoracc1, Clinical$amount_of_alcohol_consumption_per_day, type = 'p', pch = 20)
plot (Clinical$tumoracc2, Clinical$amount_of_alcohol_consumption_per_day, type = 'p', pch = 20)
cor.test (Clinical$tumoracc2, Clinical$amount_of_alcohol_consumption_per_day, use = 'pairwise.complete.obs', method = 'spearman')
plot (Clinical$normal_age, Clinical$amount_of_alcohol_consumption_per_day, type = 'p', pch = 20)
cor.test (Clinical$normal_age, Clinical$amount_of_alcohol_consumption_per_day, use = 'pairwise.complete.obs', method = 'spearman')

# vs. clinical_m
boxplot (cancer_age ~ clinical_m, data = Clinical)
anova_clinical_m <- aov (cancer_age ~ clinical_m, data = Clinical)
summary (anova_clinical_m)
kruskal.test (cancer_age ~ clinical_m, data = Clinical)
boxplot (tumoracc1 ~ clinical_m, data = Clinical)
boxplot (tumoracc2 ~ clinical_m, data = Clinical)

# vs. clinical_n
boxplot (cancer_age ~ clinical_n, data = Clinical)

# vs. clinical_stage
boxplot (cancer_age ~ clinical_stage, data = Clinical)

# vs. clinical_t
boxplot (cancer_age ~ clinical_t, data = Clinical)

# vs. days_to_completion_of_curative_tx
plot (Clinical$cancer_age, Clinical$days_to_completion_of_curative_tx, type = 'p', pch = 20)
cor (Clinical$cancer_age, Clinical$days_to_completion_of_curative_tx, use = 'pairwise.complete.obs', method = 'spearman')

# vs. disease_after_curative_tx
boxplot (cancer_age ~ disease_after_curative_tx, data = Clinical)
t.test (cancer_age ~ disease_after_curative_tx, data = Clinical)

boxplot (tumoracc2 ~ disease_after_curative_tx, data = Clinical)
t.test (tumoracc2 ~ disease_after_curative_tx, data = Clinical)

# vs. followup_treatment_success
boxplot (cancer_age ~ followup_treatment_success, data = Clinical)
boxplot (tumoracc1 ~ followup_treatment_success, data = Clinical)
boxplot (tumoracc2 ~ followup_treatment_success, data = Clinical)
kruskal.test (cancer_age ~ followup_treatment_success, data = Clinical)
aov_followup <- aov (cancer_age ~ followup_treatment_success, data = Clinical)
summary (aov_followup)

aov_followup_acc <- aov (tumoracc2 ~ followup_treatment_success, data = Clinical)
summary (aov_followup_acc)

kruskal.test(tumoracc2 ~ followup_treatment_success, data = Clinical)
# vs. frequency_of_alcohol_consumption
plot (Clinical$cancer_age, Clinical$frequency_of_alcohol_consumption, type = 'p', pch = 20)
cor (Clinical$cancer_age, Clinical$frequency_of_alcohol_consumption, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$cancer_age, Clinical$frequency_of_alcohol_consumption, use = 'pairwise.complete.obs', method = 'spearman')

cor.test (Clinical$tumoracc2, Clinical$frequency_of_alcohol_consumption, use = 'pairwise.complete.obs', method = 'spearman')

# vs. history_of_neoadjuvant_treatment
boxplot (cancer_age ~ history_of_neoadjuvant_treatment, data = Clinical)
t.test (cancer_age ~ history_of_neoadjuvant_treatment, data = Clinical, alternative = 'two.sided') # SIGNIFICANT
wilcox.test(cancer_age ~ history_of_neoadjuvant_treatment, data = Clinical, alternative = 'two.sided') # SIGNIFICANT
# but only 5 Yes, unlikely to be interesting

# vs. hpv_call_1
boxplot (cancer_age ~ hpv_call_1, data = Clinical)
kruskal.test (cancer_age ~ hpv_call_1, data = Clinical)
anova_hpv_call_1 <- aov (cancer_age ~ hpv_call_1, data = Clinical)
summary (anova_hpv_call_1) # !!! Interesting

# vs. hpv_status
boxplot (cancer_age ~ hpv_status, data = Clinical)
kruskal.test (cancer_age ~ hpv_status, data = Clinical)
t.test (cancer_age ~ hpv_status, data = Clinical)
# vs. laterality
boxplot (cancer_age ~ laterality, data = Clinical)

# vs. lymph_node_examined_count
plot (Clinical$cancer_age, Clinical$lymph_node_examined_count, type = 'p', pch = 20)
cor (Clinical$cancer_age, Clinical$lymph_node_examined_count, use = 'pairwise.complete.obs', method = 'spearman')

# vs. lymphovascular_invasion_present
boxplot (cancer_age ~ lymphovascular_invasion_present, data = Clinical)

# vs. margin_status
boxplot (cancer_age ~ margin_status, data = Clinical)
kruskal.test (cancer_age ~ margin_status, data = Clinical)

#######################################################################################################

# Cox/Logrank analysis survival data

library(survival)

days <- ifelse(is.na(Clinical$days_to_death),Clinical$days_to_last_followup,Clinical$days_to_death)
# merging the two columns. It means: if days_to_death is NA, pick days_to_last_followup
# else (days_to_death is not NA), pick it
Clinical <- cbind (Clinical, days)
# days is nbr of days between diagnosis and death (if vital_status == 1 = dead)
# or nbr of days between diagnosis and last followup (if vital_status == 0 = alive)

Clinical$survival <- with(Clinical, Surv(days, vital_status == 1))
# days are noted with a + if right-censored (ie. patient is alive)
# if no +, then patient is dead
# it is a 'Surv' object

plot (survfit (survival ~ 1, data = Clinical), mark.time = TRUE, mark = 'line')

# Survival vs. gender
gender_surv <- survfit (survival ~ gender, data = Clinical)
plot (gender_surv, col = c('red', 'blue'), mark.time = TRUE, mark = 'line', yaxs = 'i',
      xlab = 'Time (days)', ylab = 'Proportion alive', conf = TRUE) # red is female

cox_gender <- coxph (survival ~ gender, data = Clinical)
summary (cox_gender) # significant
survdiff (survival ~ gender, data = Clinical)

# Survival vs. radiation therapy
radiation_surv <- survfit (survival ~ radiation_therapy, data = Clinical)
plot (radiation_surv, col = c('red', 'blue'), mark.time = TRUE, mark = 'line', yaxs = 'i', 
      xlab = 'Time (days)', ylab = 'Proportion alive') # red is no
cox_radiation <- coxph (survival ~ radiation_therapy, data = Clinical)
summary (cox_radiation) # highly significant
survdiff (survival ~ radiation_therapy, data = Clinical) # actually it is included in Cox

# Survival vs. cancer_age
plot (Clinical$cancer_age, Clinical$days, type = 'p', pch = 20)
cox_cancer_age <- coxph (survival ~ cancer_age, data = Clinical) 
summary (cox_cancer_age) # not significant

# Survival vs. tumor acceleration 1
plot (Clinical$tumoracc1, Clinical$days, type = 'p', pch = 20)
cor (Clinical$tumoracc1, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$tumoracc1, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')

cox_tumoracc1 <- coxph (survival ~ tumoracc1, data = Clinical)
summary (cox_tumoracc1) # SIGNIFICANT

cox_normalage <- coxph (survival ~ normal_age, data = Clinical)
summary (cox_normalage) #NS

cox_chrono <- coxph (survival ~ years_to_birth, data = Clinical)
summary (cox_chrono) #S

# Survival vs. tumor acceleration 2
plot (Clinical$tumoracc2, Clinical$days, type = 'p', pch = 20)
cor (Clinical$tumoracc2, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$tumoracc2, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')

cox_tumoracc2 <- coxph (survival ~ tumoracc2, data = Clinical)
summary (cox_tumoracc2) # Not significant

# Survival vs normal_age
plot (Clinical$normal_age, Clinical$days, type = 'p', pch = 20)
cor (Clinical$normal_age, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')
cor.test (Clinical$normal_age, Clinical$days, use = 'pairwise.complete.obs', method = 'spearman')

cox_normal_age <- coxph (survival ~ normal_age, data = Clinical)
summary (cox_normal_age) # Not significant

cox_sum <- coxph (survival ~ cancer_age + tumoracc2, data = Clinical)
summary (cox_sum)

cox_sum2 <- coxph (survival ~ cancer_age + years_to_birth, data = Clinical)
summary (cox_sum2)