setwd("~/Dropbox/Research/2022/UPMC_MS/Code/Phenotype")
library(readr)
library(dplyr)
library(readxl)
library(ggplot2)
source("Phenotyping_eval.R")
library(Matrix)
library(pROC)
library(KOMAP)

#### Data preprocessing ####
data_path <- "~/Dropbox/Research/2022/UPMC_MS/Data/"

# Since we need to use both NLP and PheCode, we use truncated version from 2011 - 2021
ehr <- read_csv(paste(data_path, "Latest/UPMC_MS_2011_to_2021_Codified+NLP_data_2023-05-24.csv", sep=""),
                show_col_types = FALSE)
ehr.main <- ehr %>% group_by(patient_num) %>% 
  dplyr::summarize(utl = n(),
                   ICD = sum(feature_id == "PheCode:335"),
                   NLP = sum(feature_id == "C0026769"),.groups = "drop")
ehr.main$ICD.NLP <- ehr.main$ICD + ehr.main$NLP
ehr.main <- ehr.main[!is.na(ehr.main$patient_num),]
# 19 Relevant codified features from ONCE
codified.feat <- read_csv(paste0(data_path, "Latest/ONCE_Multiple sclerosis_Codified.csv"),
                          show_col_types = FALSE)
feat.set.thresholded <- codified.feat %>% filter(phenotyping_features == 1)
feat.set <- setdiff(feat.set.thresholded$Variable, "PheCode:335")
# 30 relevant CUI features from ONCE
nlp.feat <- read_csv(paste0(data_path, "Latest/ONCE_Multiple sclerosis_CUI.csv"),
                     show_col_types = FALSE)
feat.set <- c(feat.set, nlp.feat$cui[nlp.feat$phenotyping_features == TRUE])
feat.set <- setdiff(feat.set, "C0026769")
write.csv(feat.set, paste0(data_path, "Processed/Once_Feature_List_20231103.csv"), row.names=FALSE)

ehr.set <- ehr %>% filter(feature_id %in% feat.set)
ehr.set <- ehr.set %>% group_by(patient_num, feature_id) %>% summarize(count=n(),.groups = "drop")
ehr.set = tidyr::spread(ehr.set, key = feature_id, value = count, fill = 0)
ehr.dat <- left_join(ehr.main, ehr.set, by="patient_num")
ehr.dat[is.na(ehr.dat)] <- 0

#### Add labels from MS chart review or registry data ####

# labels from MS chart review
MS_Annotation <- read_excel(paste0(data_path,"/Latest/MS_Annotation_499_20230626.xlsx"))
colnames(MS_Annotation) = c("PATIENT_NUM", "SEX", "ENCOUNTER", "MS", "DATE", "ALIVE")

pt.id.missingEHR = setdiff(MS_Annotation$PATIENT_NUM, ehr.dat$patient_num)
write.csv(pt.id.missingEHR, paste0(data_path, "Processed/UPMC_MS_Phenotype_PatientList_MissingEHR_20231103.csv"), row.names=FALSE)


pos.id.ehr = as.matrix(MS_Annotation %>% filter(MS == 1) %>% select(PATIENT_NUM)) 
neg.id.ehr = as.matrix(MS_Annotation %>% filter(MS %in% c(2,3)) %>% select(PATIENT_NUM))
ehr.dat$is.MS = NA
ehr.dat$is.MS[ehr.dat$patient_num %in% pos.id.ehr] = 1 
ehr.dat$is.MS[ehr.dat$patient_num %in% neg.id.ehr] = 0

# labels from MS registry data
linkage_file <- paste(data_path, "ID Linkage/R3_2468_Xia_PROMOTE_20220223_studyid_NoPHI.xls", sep="")
linkage <- read_excel(linkage_file) # 2445 pts in registry
demo <- read_csv(paste0(data_path, "Latest/ClinicalDemographics_230425.csv"))
registry_ehr <- merge(linkage, demo[,c("id_participant","enrollment_diagnosis")], 
                      by="id_participant")
ehr_id = unique(ehr$patient_num)
registry_ehr <- registry_ehr[registry_ehr$PATIENT_NUM %in% ehr_id,]
ehr.dat$is.Registry <- ehr.dat$patient_num %in% registry_ehr$PATIENT_NUM
pos.id.registry = unique(as.matrix(registry_ehr %>% 
                                     filter(enrollment_diagnosis %in% c(1:6,9,10)) %>% 
                                     select(PATIENT_NUM)))
neg.id.registry = unique(as.matrix(registry_ehr %>% 
                                     filter(enrollment_diagnosis %in% c(7, 8, 11, 12)) %>% 
                                     select(PATIENT_NUM)))
ehr.dat$is.MS[ehr.dat$patient_num %in% pos.id.registry] = 1 
ehr.dat$is.MS[ehr.dat$patient_num %in% neg.id.registry] = 0

ehr.dat$is.Labeled <- 1- is.na(ehr.dat$is.MS)
save(ehr.dat, file=paste0(data_path,"/Processed/UPMC_MS_Phenotyping_20231103.RData"))
write.csv(ehr.dat, paste0(data_path, "Processed/UPMC_MS_Phenotyping_20231103.csv"), row.names=FALSE)

#### MAP Phenotyping ####
load(file=paste0(data_path,"/Processed/UPMC_MS_Phenotyping_20231103.RData")) ##-> ehr.dat

ehr.map <- ehr.dat %>% filter(ICD > 0) # only keep pts with main PheCode
cor(log(1+ehr.map[,"ICD"]), log(1+ehr.map[,"NLP"]))
mat <- Matrix(data=as.matrix(ehr.map[,c("ICD","NLP")]),sparse = TRUE)
note <- Matrix(ehr.map$utl,ncol=1,sparse = TRUE)
## predict risk score ##
library("MAP")
res <- MAP(mat = mat,  note=note)
phenotype.res <- ehr.map[,c("patient_num", "is.MS", "is.Registry", "is.Labeled")]
phenotype.res$MAP <- as.numeric(res$scores)
# phenotype.res$MAP.pred.cut <- as.numeric(res$scores > res$cut.MAP)

## evaluation using only cases ##
# s <- as.numeric(ehr.sum$patient_num %in% linkage$PATIENT_NUM)
# phenotype_res <- Phenotype_eval(x=log(pred+0.1), s=s, cutoff= log(res$cut.MAP+0.1), bootstrap = FALSE)
# rocobj <- roc(phenotype_res$y, pred)
# phenotype_res_CI <- Phenotype_eval(x=log(pred+0.1), s=s, cutoff= log(res$cut.MAP+0.1),bootstrap = TRUE)
# phenotype_res_CI$results
# phenotype_res_CI$results_CI

#Evaluation metrics
# tmp <- matrix(data.frame(phenotype_res_CI$results), ncol=2, byrow=TRUE)
# rownames(tmp) <- c("auc","fpr","sens","ppv","f_score")
# colnames(tmp) <- c("mean", "sd")
# tmp

#### KOMAP Phenotyping ####
set.seed(2023)
ehr_logcount_wide <- ehr.map[,c(1:4,6:53)]
ehr_logcount_wide[,-1] <- log(ehr_logcount_wide[, -1] + 1)
main_code = 'ICD'
ehr_corrupt = as.data.frame(ehr_logcount_wide[,-1])
ehr_corrupt$corrupt_mainICD = ehr_corrupt[, main_code]
ehr_corrupt$corrupt_mainICD[sample(1:nrow(ehr_corrupt), round(nrow(ehr_corrupt) * 0.2), replace = FALSE)] = mean(as.matrix(ehr_corrupt[, main_code]))
main_code = 'NLP'
ehr_corrupt$corrupt_mainNLP = ehr_corrupt[, main_code]
ehr_corrupt$corrupt_mainNLP[sample(1:nrow(ehr_corrupt), round(nrow(ehr_corrupt) * 0.2), replace = FALSE)] = mean(as.matrix(ehr_corrupt[, main_code]))

id.train = sample(1:nrow(ehr_corrupt), round(nrow(ehr_corrupt) / 3 * 2))
id.valid = setdiff(1:nrow(ehr_corrupt), id.train)
dat.cov.train = cov(ehr_corrupt[id.train, ])
dat.cov.valid = cov(ehr_corrupt[id.valid, ])
out <- KOMAP_corrupt(dat.cov.train, dat.cov.valid, is.wide = TRUE, 'ICD', 'NLP', 
                       nm.disease = 'MS', 'utl', nm.multi = NULL, nm.corrupt.code = 'corrupt_mainICD', 
                       nm.corrupt.cui = 'corrupt_mainNLP', 
                       pred = TRUE, eval.real = FALSE, eval.sim = FALSE,
                       dat.part = ehr_logcount_wide, nm.id = 'patient_num')
colnames(out$pred_prob$pred.score)
phenotype.res$KOMAP.ICD.EHR = out$pred_prob$pred.score[,2]
phenotype.res$KOMAP.ICDNLP.EHR = out$pred_prob$pred.score[,3]

#### PheNorm Phenotyping ####
ehr.phenorm <- ehr.map
ehr.phenorm$patient_num <- NULL
ehr.phenorm$is.MS <- NULL
ehr.phenorm$is.Registry <- NULL
ehr.phenorm$is.Labeled <- NULL
ehr.phenorm <- log(ehr.phenorm + 1)

library("PheNorm")
set.seed(1234) 
phenorm.icd=PheNorm.Prob("ICD", "utl", ehr.phenorm, nm.X = NULL, 
                         corrupt.rate=0.3, train.size=nrow(ehr.phenorm)) 
phenorm.icd.cov=PheNorm.Prob("ICD", "utl", ehr.phenorm, 
                             nm.X = setdiff(colnames(ehr.phenorm), c("ICD", "ICD.NLP","utl")), 
                             corrupt.rate=0.3, train.size=nrow(ehr.phenorm)) 
phenorm.icdnlp=PheNorm.Prob(c("ICD", "NLP","ICD.NLP"), "utl", ehr.phenorm, 
                            nm.X = NULL, 
                            corrupt.rate=0.3, train.size=nrow(ehr.phenorm)) 
phenorm.icdnlp.cov=PheNorm.Prob(c("ICD", "NLP","ICD.NLP"), "utl", ehr.phenorm, 
                                nm.X = setdiff(colnames(ehr.phenorm), c("ICD", "NLP","ICD.NLP","utl")), 
                                corrupt.rate=0.3, train.size=nrow(ehr.phenorm)) 

phenotype.res$PheNorm.ICD = phenorm.icd$probs
phenotype.res$PheNorm.ICD.EHR = phenorm.icd.cov$probs
phenotype.res$PheNorm.ICDNLP = phenorm.icdnlp$probs
phenotype.res$PheNorm.ICDNLP.EHR = phenorm.icdnlp.cov$probs

## evaluation using only cases ##
# phe.icd <- Phenotype_eval(x=log(phenorm.icd$probs+0.1), s=s, bootstrap = FALSE)
# phe.icd.cov <- Phenotype_eval(x=log(phenorm.icd.cov$probs+0.1), s=s, bootstrap = FALSE)
# phe.nlp <- Phenotype_eval(x=log(phenorm.nlp$probs+0.1), s=s, bootstrap = FALSE)
# phe.nlp.cov <- Phenotype_eval(x=log(phenorm.nlp.cov$probs+0.1), s=s, bootstrap = FALSE)
# phe.icdnlp <- Phenotype_eval(x=log(phenorm.icdnlp$probs+0.1), s=s, bootstrap = FALSE)
# phe.icdnlp.cov <- Phenotype_eval(x=log(phenorm.icdnlp.cov$probs+0.1), s=s, bootstrap = FALSE)
# roc_list = list(ICD=roc(phe.icd$y, phenorm.icd$probs),
#                 ICD.EHR.Feature=roc(phe.icd.cov$y, phenorm.icd.cov$probs),
#                 NLP=roc(phe.nlp$y, phenorm.nlp$probs),
#                 NLP.EHR.Feature=roc(phe.nlp.cov$y, phenorm.nlp.cov$probs),
#                 ICD.NLP=roc(phe.icdnlp$y, phenorm.icdnlp$probs),
#                 ICD.NLP.EHR.Feature=roc(phe.icdnlp.cov$y, phenorm.icdnlp.cov$probs))
# pp <- ggrocs(roc_list, legendTitel="")
# pp

# pred <- phenorm.icdnlp.cov$probs
# rocobj <- roc(phe.icdnlp.cov$y, pred)
# phenotype_res_CI <- Phenotype_eval(x=log(pred+0.1), s=s, bootstrap = TRUE)

#Evaluation metrics
# tmp <- matrix(data.frame(phenotype_res_CI$results), ncol=2, byrow=TRUE)
# rownames(tmp) <- c("auc","fpr","sens","ppv","f_score")
# colnames(tmp) <- c("mean", "sd")
# tmp

#### Baselines ####
phenotype.res$Log.ICD = ehr.phenorm$ICD
phenotype.res$Log.NLP = ehr.phenorm$NLP

save(phenotype.res, file=paste0(data_path,"/Processed/UPMC_MS_Phenotyping_res_20231103.RData"))

#### Evaluation ####
load(file=paste0(data_path,"/Processed/UPMC_MS_Phenotyping_res_20231103.RData")) ##-> phenotype.res

library(htmlTable)
# Registry patients
registry.res <- phenotype.res %>% filter(is.Labeled == 1,
                                         is.Registry == 1)
tmp2 <- as.list(registry.res[,c(5:ncol(registry.res))])
res.97 <- binary_eval(registry.res$is.MS, tmp2, cutoff = 0.97)
res.95 <- binary_eval(registry.res$is.MS, tmp2, cutoff = 0.95)
res.all <- merge(res.97, res.95, by=c("Method", "AUROC", "AUPRC"), 
                 suffixes = c("(at 97% Specificity)", "(at 95% Specificity)"))
t(res.all) %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable

res.90 <- binary_eval(registry.res$is.MS, tmp2, cutoff = 0.90)
res.80 <- binary_eval(registry.res$is.MS, tmp2, cutoff = 0.80)
res.all <- merge(res.90, res.80, by=c("Method", "AUROC", "AUPRC"), 
                 suffixes = c("(at 90% Specificity)", "(at 80% Specificity)"))
t(res.all) %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable

# Non-Registry patients
ehr.res <- phenotype.res %>% filter(is.Labeled == 1,
                                    is.Registry == 0)
tmp3 <- as.list(ehr.res[,c(5:ncol(ehr.res))])
res.97 <- binary_eval(ehr.res$is.MS, tmp3, cutoff = 0.97)
res.95 <- binary_eval(ehr.res$is.MS, tmp3, cutoff = 0.95)
ehr.res.all <- merge(res.97, res.95, by=c("Method", "AUROC", "AUPRC"), 
                 suffixes = c("(at 97% Specificity)", "(at 95% Specificity)"))
t(ehr.res.all) %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable

res.90 <- binary_eval(ehr.res$is.MS, tmp3, cutoff = 0.90)
res.80 <- binary_eval(ehr.res$is.MS, tmp3, cutoff = 0.80)
ehr.res.all <- merge(res.90, res.80, by=c("Method", "AUROC", "AUPRC"), 
                     suffixes = c("(at 90% Specificity)", "(at 80% Specificity)"))
t(ehr.res.all) %>% 
  addHtmlTableStyle(col.rgroup = c("none", "#F7F7F7")) %>%
  htmlTable

#### Patients with MS phenotpye ####
# The best model on non-registry patients is PheNorm.ICD.EHR
# To make it consistent with MGB algorithm, we save results for KOMAP
ehr.res <- phenotype.res %>% filter(is.Labeled == 1,
                                    is.Registry == 0)
tmp <- roc(ehr.res$is.MS, ehr.res$KOMAP.ICDNLP.EHR)
cut90 <- tmp$thresholds[which.min(abs(tmp$specificities - 0.90))]
cut95 <- tmp$thresholds[which.min(abs(tmp$specificities - 0.95))]

MS.out <- phenotype.res[,c("patient_num", "is.Labeled", "is.MS")]
MS.out$pred.MS.90specificity <- phenotype.res$KOMAP.ICDNLP.EHR > cut90
MS.out$pred.MS.95specificity <- phenotype.res$KOMAP.ICDNLP.EHR > cut95
MS.out = left_join(ehr.dat[,c("patient_num", "is.Labeled", "is.MS", "is.Registry")], 
                   MS.out, 
                   by=c("patient_num", "is.Labeled", "is.MS"))
summary(MS.out)
write.csv(MS.out, paste0(data_path, "Processed/UPMC_KOMAP_MS_Phenotype_PatientList_20231103.csv"), row.names=FALSE)

KOMAP.coef <- out$est$long_df
dat1 <- KOMAP.coef %>% filter(method == "mainICD + codify") %>% select(feat, coeff)
dat2 <- KOMAP.coef %>% filter(method == "mainICDNLP + codify & NLP") %>% select(feat, coeff)
KOMAP.coef = merge(dat1, dat2, by=c("feat"), 
                   suffixes = c(".KOMAP.ICD.EHR", ".KOMAP. ICDNLP.EHR"))
write.csv(KOMAP.coef, paste0(data_path, "Processed/UPMC_KOMAP_MS_Coeff_20231103.csv"), row.names=FALSE)


# MS.out$is.MS[is.na(MS.out$is.MS)] <- 99
# MS.KG$MS <- MS.KG$is.Labeled * MS.KG$is.MS + 
#   (1-MS.KG$is.Labeled) * MS.KG$pred.MS
# write.csv(MS.KG$patient_num[MS.KG$MS == 1], paste0(data_path, "Processed/UPMC_MS_Phenotype_PatientList_20230501.csv"), row.names=FALSE)


## ROC plots
# pos_id = pos.id.registry
# neg_id = neg.id.registry
# labeled_id = intersect(unique(c(pos_id, neg_id)), ehr.sum$patient_num)
# y = as.numeric(ehr.sum$patient_num[ehr.sum$patient_num %in% labeled_id] %in% pos_id)
# roc_list = list(ICD=roc(y, phenorm.icd$probs[ehr.sum$patient_num %in% labeled_id]),
#                 ICD.EHRFeature=roc(y, phenorm.icd.cov$probs[ehr.sum$patient_num %in% labeled_id]),
#                 # NLP=roc(y, phenorm.nlp$probs[ehr.sum$patient_num %in% labeled_id]),
#                 # NLP.EHRFeature=roc(y, phenorm.nlp.cov$probs[ehr.sum$patient_num %in% labeled_id]),
#                 ICD.NLP=roc(y, phenorm.icdnlp$probs[ehr.sum$patient_num %in% labeled_id]),
#                 ICD.NLP.EHRFeature=roc(y, phenorm.icdnlp.cov$probs[ehr.sum$patient_num %in% labeled_id]))
# pp <- ggrocs(roc_list, legendTitel="")
# pp

