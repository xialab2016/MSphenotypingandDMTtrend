library(readr)
library(dplyr)
library(parallel)
library(readxl)

## patient list ####
# Load phenotyped MS patients from UPMC and MGB (Mass General Brigham) cohorts.
# Retain only patients predicted as MS with >= 90% specificity.
pe_UPMC <- read_csv("UPMC_KOMAP_MS_Phenotype_PatientList_20241104.csv")
pe_UPMC_list <- unique(pe_UPMC$patient_num[pe_UPMC$pred.MS.90specificity==1])
pe_MGB <- read_csv("MS_algorithm_diagnosis.csv")
pe_MGB_list <- unique(pe_MGB$patient_num[pe_MGB$MS.spec.90==1])

## handle UPMC data ####
# Load UPMC longitudinal EHR data, restrict to MS patients, and parse the year from visit date.
# Normalize procedure codes: "CCS-PCS:" prefix is collapsed to "CCS:" for consistency with MGB.
# Keep only the four coded feature types used downstream (RXNORM, PheCode, LOINC, CCS).
x_UPMC <- read_csv("UPMC_MS_2004_to_2022_Codified_processed_data_2023-08-25.csv")
x_UPMC <- x_UPMC %>%
  filter(patient_num %in% pe_UPMC_list) %>%
  mutate(year = format(start_date, "%Y"))
x_UPMC$feature_id <- gsub("CCS-PCS:", "CCS:", x_UPMC$feature_id)
x_UPMC$group <- sapply(x_UPMC$feature_id, function(x) strsplit(x,":")[[1]][1])
x_UPMC <- x_UPMC %>%
  filter(group %in% c("RXNORM","PheCode","LOINC","CCS"))
dict <- x_UPMC %>%
  group_by(feature_id) %>%
  summarize(n_UPMC = n())

## handle MGB data ####
# MGB data arrives as monthly aggregate counts rather than individual visit dates.
# We reconstruct approximate event-level records by distributing each month's count
# evenly across the 30-day window, anchored to the patient's MS start date.
load("MGB_MS_patient_StartDate.Rdata")
MGB.pe <- MGB.pe %>%
  filter(PatientNum %in% pe_MGB_list)
colnames(MGB.pe) <- c("patient_num", "start_date")

# Each file contains one month's worth of aggregate (feature, count) data across all patients.
files <- list.files(path = "MGB/MS-AD/monthly_codified/MS", pattern = "^daily_patient_", full.names = TRUE)
x_MGB <- data.frame(paitent_num = numeric(),
                    start_time = numeric(),
                    feature_id = character(),
                    year = character(),
                    group = character())
for(file in files){
  load(file)
  colnames(data.grouped) <- c("patient_num", "feature_id", "month", "n")
  data.grouped <- data.grouped %>%
    filter(grepl(":", feature_id)) %>%          # drop features without a coded prefix
    filter(patient_num %in% pe_MGB_list) %>%
    mutate(group = sapply(feature_id, function(x) strsplit(x,":")[[1]][1])) %>%
    filter(group %in% c("RXNORM","PheCode","LOINC","CCS")) %>%
    merge(MGB.pe)                               # attach each patient's MS start date

  # Compute each month's offset from the patient's first MS code (PheCode:335),
  # so that time is expressed relative to disease onset rather than calendar date.
  pe.tmp <- data.grouped %>%
    filter(feature_id=="PheCode:335") %>%
    group_by(patient_num) %>%
    summarize(start_month = min(month))
  data.grouped <- data.grouped %>%
    merge(pe.tmp) %>%
    mutate(n_month = month - start_month)

  # Expand aggregate monthly counts into individual pseudo-event rows.
  # Each of the n occurrences in a month is spread evenly across the 30-day window.
  tmp <- data.grouped$n
  data.grouped <- data.grouped[rep(1:length(tmp), tmp),]
  data.grouped$tadd <- unlist(sapply(tmp, function(x) 1:x)) - 1
  data.grouped$time <- data.grouped$n_month * 30 + round(data.grouped$tadd / data.grouped$n * 30)
  data.grouped$start_date <- data.grouped$time + data.grouped$start_date
  data.grouped$year <- format(data.grouped$start_date, "%Y")
  x_MGB <- rbind(x_MGB, data.grouped[,c('patient_num','start_date','feature_id','year','group')])
}
rm(MGB.pe)

save(x_MGB, x_UPMC, file = "SPPMI_rawdata_filtered_pe.Rdata")

## select features ####
# Merge both sites into a single long-format table and build the feature vocabulary.
# Include a feature if it appears > 100 times OR belongs to the MS drug list / KOMAP MS feature set.
rm(list = ls())
load("SPPMI_rawdata_filtered_pe.Rdata")
x <- rbind(x_MGB, x_UPMC)
rm(x_MGB, x_UPMC)

# Load curated MS drug list and map RxNorm IDs to "RXNORM:XXXXX" format.
Drug <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
Drug <- na.omit(unique(Drug$`RxNorm Ingredient id`))
Drug <- paste("RXNORM", Drug, sep=":")

dict <- x %>%
  group_by(feature_id) %>%
  summarise(n = n())

summary(dict$n[dict$feature_id%in%Drug])
# ONCE output: cosine-similarity-based MS feature selection (threshold 0.165)
mydict <- read_csv("ONCE_multiple sclerosis_PheCode335_cos0.165.csv")
summary(dict$n[dict$feature_id%in%mydict$Variable])

# Retain features that are either in the curated drug/ONCE lists or have sufficient global frequency.
dict <- dict %>%
  filter((feature_id %in% c(Drug,mydict$Variable)) | n > 100)

dict$id <- 1:nrow(dict)
x <- x %>%
  filter(feature_id %in% dict$feature_id) %>%
  mutate(id = dict$id[match(feature_id, dict$feature_id)])

save(x, dict, file = "SPPMI_rawdata_filtered.Rdata")

## Adding features ####
# For two specific drugs (RXNORM:117055 = glatiramer, RXNORM:44157 = interferon beta-1a),
# a single prescription record is propagated forward for 24 months to reflect ongoing treatment,
# since these are continuous-use disease-modifying therapies.
rm(list = setdiff(ls(),c("x")))
library(tidyr)
library(lubridate)

tmpx <- x %>% filter(feature_id %in% c("RXNORM:117055","RXNORM:44157"))
generate_dates <- function(date, months = 24) {
  seq(from = date, by = "month", length.out = months)
}
expanded_tmpx <- tmpx %>%
  rowwise() %>%
  mutate(
    additional_dates = list(generate_dates(start_date, 24))
  ) %>%
  unnest(additional_dates) %>%
  mutate(
    start_date = additional_dates
  ) %>%
  select(-additional_dates)

# Deduplicate within patient-feature-year-month to avoid inflating co-occurrence counts.
expanded_tmpx <- expanded_tmpx %>%
  mutate(
    year = format(start_date, "%Y"),
    month = format(start_date, "%m")
  ) %>%
  group_by(patient_num, feature_id, year, month) %>%
  distinct(patient_num, feature_id, year, month, .keep_all = TRUE) %>%
  ungroup() %>%
  arrange(patient_num, start_date)
expanded_tmpx <- expanded_tmpx[,match(colnames(x),colnames(expanded_tmpx))]
x <- x %>% filter(!feature_id %in% c("RXNORM:117055","RXNORM:44157"))
x <- rbind(x, expanded_tmpx) %>%
  arrange(patient_num, start_date)
rownames(x) <- NULL
dict <- x %>%
  group_by(feature_id, id) %>%
  summarise(n = n())
save(x, dict, file = "SPPMI_rawdata_filtered.Rdata")

# Add drug treatment-class features: collapse individual drugs to their class label
# so that class-level embeddings can be computed alongside drug-level ones.
load("SPPMI_rawdata_filtered.Rdata")
Drug <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
yearlist = 2000:2023
Drug = Drug[!is.na(Drug$`RxNorm Ingredient id`),]
Drug$feature_id = paste("RXNORM", Drug$`RxNorm Ingredient id`,sep=":")
y = x %>%
  filter(feature_id %in% Drug$feature_id) %>%
  mutate(feature_id = Drug$Treatment_Class[match(feature_id,Drug$feature_id)])
x = rbind(x, y)
x <- x %>%
  arrange(patient_num, start_date)
rownames(x) <- NULL
dict <- x %>%
  group_by(feature_id) %>%
  summarise(n = n())
dict$id <- 1:nrow(dict)
x$id <- dict$id[match(x$feature_id,dict$feature_id)]
save(x, dict, file = "SPPMI_rawdata_filtered.Rdata")

x <- x[,c('patient_num','start_date','id')]
write_csv(x, file = "SPPMI_rawdata_filtered.csv")
rm(list = ls())

## get feature frequency ####
# Compute per-year patient-level prevalence (proportion of active patients who had the feature)
# and total event count for each feature, used later to normalize cosine similarity plots.
library(dplyr)
load("SPPMI_rawdata_filtered.Rdata")
data = x %>%
  select(patient_num, year, feature_id) %>%
  unique() %>%                            # one record per patient-year-feature
  group_by(year, feature_id) %>%
  summarise(n = n())                      # number of patients with the feature that year
data2 = x %>%
  select(patient_num, year) %>%
  unique() %>%
  group_by(year) %>%
  summarise(npe = n())                    # total active patients per year (denominator)
data = data %>%
  merge(data2) %>%
  mutate(freq = n/npe) %>%
  filter(year >= 2000, year <= 2023)
data2 = x %>%
  group_by(year, feature_id) %>%
  summarise(totalcount = n()) %>%        # total event count (not patient-deduplicated)
  filter(year >= 2000, year <= 2023)
data = data %>%
  merge(data2)

save(data, file = "rawdata_frequency_dict.Rdata")


## compute SVD SPPMI ####
# For each calendar year, build a co-occurrence matrix from precomputed pair counts,
# then apply Shifted Positive PMI (SPPMI) and SVD to obtain low-dimensional embeddings.
# Features with row-sum <= 10 are excluded from a given year's embedding to avoid
# instability from sparse co-occurrence data.
library(dplyr)
library(readr)
library(readxl)
load("SPPMI_rawdata_filtered.Rdata")

# Exclude patients whose records fall outside 2000-2022 to keep embeddings comparable.
pe = unique(x$patient_num[which(x$year<2000|x$year>2023)])
countdict = x %>%
  filter(!patient_num%in%pe) %>%
  group_by(feature_id, year) %>%
  summarize(n = n())
rm(list = setdiff(ls(), c("dict","countdict")))

# Load pre-aggregated co-occurrence counts (computed externally, one CSV per year).
cooc = list()
file = "result/"
for(year in 2000:2022){
  tmp = read_csv(paste(file, "MS_cooc_", year, ".csv", sep = ""))
  cooc[[as.character(year)]] = tmp
}
n = length(cooc)
d = max(dict$id)
namelist = dict$feature_id[match(1:d, dict$id)]
embed = list()
Drug <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
Drug = Drug[!is.na(Drug$`RxNorm Ingredient id`),]
Drug$feature_id = paste("RXNORM", Drug$`RxNorm Ingredient id`,sep=":")

for(i in 1:n){
  freq = countdict %>% filter(year == i+1999)
  freq = freq$n[match(namelist, freq$feature_id)]
  names(freq) = namelist

  # Reconstruct the symmetric co-occurrence matrix from the stored (i, j, count) triplets.
  tmp = cooc[[i]]
  co = matrix(0, nrow = d, ncol = d)
  co[cbind(tmp$i, tmp$j)] = tmp$count
  co = co + t(co)
  diag(co) = diag(co)/2   # diagonal was doubled by the symmetrization above
  rownames(co) = colnames(co) = namelist

  co = co + t(co)
  diag(co) = diag(co)/2

  cooc[[i]] = co

  # SPPMI transformation: normalize co-occurrence by marginal frequencies to get PMI,
  # then shift by log(k)=0 (i.e., standard PPMI) and zero out negative values.
  tmp = co
  rs = rowSums(tmp)

  idx = which(rs>10)   # drop low-frequency features before SVD
  tmp = tmp[idx,idx]
  rs = rs[idx]
  tmp = (t(tmp/rs))/rs * sum(rs)
  tmp = log(tmp)
  tmp[tmp<0] = 0       # positive PMI: clip negatives to 0

  # SVD of the SPPMI matrix; keep components where u and v agree in sign
  # (i.e., semantically coherent directions), up to 100 components.
  fit = svd(tmp)
  idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
  idx = idx[1:min(100, length(idx))]
  U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))   # symmetric factorization: W = U * sqrt(D)
  rownames(U) = rownames(tmp)
  embed[[i]] = U
}
save(embed, file = "raw_embedding_SVD_SPPMI.Rdata")
save(cooc, file = "raw_threshold_cooc_SVD_SPPMI.Rdata")


## estimate variance ####
# Estimate per-feature embedding variance by computing within-patient cosine similarities
# to the MS anchor (PheCode:335) and measuring how much they vary across patients.
# This variance is used as the basis for empirical confidence intervals in the trend plots.
library(dplyr)
options(dplyr.summarise.inform = FALSE)

# Compute a co-occurrence matrix for a single patient's visit sequence within a time window.
get_cooc <- function(xx, windowsize){
  xx$feature <- match(xx$feature_id, dict$feature_id)
  xx$time <- as.numeric(difftime(as.character(xx$start_date), "2000-01-01", units="days"))
  xx <- xx[,c('feature','time')]

  # Build all (j, k) day pairs where day_j and day_k are both present and |j-k| <= windowsize.
  unique_time <- unique(xx$time)
  time_table <- data.frame(time_j = rep(unique_time, windowsize),
                           time_k = rep(unique_time, windowsize) -
                             rep(1:windowsize, each = length(unique_time))) %>%
    filter(time_j %in% unique_time & time_k %in% unique_time)
  colnames(xx) <- c('feature_j', 'time_j')
  x1 <- merge(xx, time_table, by = "time_j")
  xx <- xx[,c('feature_j', 'time_j')]
  colnames(xx) <- c('feature_k', 'time_k')
  xx <- inner_join(x1, xx, by = "time_k", relationship = "many-to-many")
  xx <- xx %>% group_by(feature_j, feature_k) %>% summarise(n = n())
  return(xx)
}

# Compute the SPPMI cosine similarity between a target feature and all others
# for a single patient's event sequence.
get_cos <- function(xx, windowsize, target){
  xx <- get_cooc(xx = xx, windowsize = windowsize)
  co <- matrix(0, nrow = nrow(dict), ncol = nrow(dict))
  co[cbind(xx$feature_j, xx$feature_k)] <- xx$n
  co <- co + t(co)
  diag(co) <- diag(co)/2
  rownames(co) <- dict$feature_id
  rs <- rowSums(co)
  idx <- which(rs > 5)
  if(length(idx) < 5 | !(target%in%rownames(co)[idx])) return(rep(NA, nrow(dict)))
  co <- co[idx,idx]
  rs <- rowSums(co)
  co <- t(co/rs)/rs * sum(rs)
  co <- log(co)
  co[co<0] <- 0
  fit <- svd(co)
  idx <- which(sign(fit$u[1,])==sign(fit$v[1,]))
  if(length(idx)<3) return(rep(NA, nrow(dict)))
  idx <- idx[1:min(100, length(idx))]
  U <- fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
  U <- U/apply(U, 1, norm, '2')   # L2-normalize rows before dot-product similarity
  U <- U %*% t(U)                 # cosine similarity matrix
  return(U[match(target,rownames(co)),match(dict$feature_id,rownames(co))])
}

# For a given year's data, iterate over patients with at least 20 records,
# compute their individual cosine similarity vectors, and return the column-wise variance.
get_cos_summary <- function(data, windowsize, target){
  patient <- data %>% group_by(patient_num) %>% summarise(n = n()) %>% filter(n > 20)
  patient <- unique(patient$patient_num)
  cosine <- matrix(nrow = 0, ncol = nrow(dict))
  for(i in 1:length(patient)){
    pe <- patient[i]
    cosine <- rbind(cosine, get_cos(data%>%filter(patient_num==pe),
                                    windowsize = windowsize,
                                    target = target))
    if(i%%10==0) cat(i,",")
    if(i%%100==0){
      save(cosine, file = "result/part_cos_var.Rdata")
      cat("\n")
    }
  }
  rownames(cosine) <- patient
  return(list(emph_var = apply(cosine, 2, var, na.rm=TRUE),
              ntotal = nrow(cosine),
              nvalid = apply(!is.na(cosine), 2, sum)))
}
load("SPPMI_rawdata_filtered.Rdata")
for(curryear in 2000:2022){
  mycosine[[as.character(curryear)]] = get_cos_summary(x %>% filter(year==curryear), windowsize = 30, target = "PheCode:335")
  save(mycosine, file = "summary_cos_var.Rdata")
}

# Aggregate variance estimates from chunked result files (computed in parallel batches).
load("dict.Rdata")
estedvar = list()
for(theyear in 0:8){
  load(paste("result/summary_cos_var_",theyear,".Rdata",sep=""))
  estedvar = c(estedvar, mycosine)
}
yearlist = names(estedvar)

# Scale variance by effective sample size: features with fewer than 50 effective observations
# are set to NA so they are excluded from CI construction in the plots.
estedvar = lapply(estedvar, function(x){
  tmp <- x$emph_var/sqrt(x$noriginal/x$ntotal*x$nvalid)
  idx <- which(x$noriginal/x$ntotal*x$nvalid <= 50)
  tmp[idx] <- NA
  return(tmp)
})
estedvar = do.call("rbind", estedvar)
rownames(estedvar) = yearlist
colnames(estedvar) = dict$feature_id
save(estedvar, file = "summary_est_var.Rdata")

## plots with empirical CI ####
# For each MS drug (and each drug class), plot the cosine similarity to PheCode:335 over time
# with 95% empirical confidence intervals derived from the patient-level variance estimates.
# Smoothing uses loess when >= 4 valid years are available, otherwise falls back to lm.

library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(patchwork)
load("raw_embedding_SVD_SPPMI.Rdata")
load("rawdata_frequency_dict.Rdata")
load("summary_est_var.Rdata")
Drug <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
yearlist = 2000:2021
Drug = Drug[!is.na(Drug$`RxNorm Ingredient id`),]
Drug$feature_id = paste("RXNORM", Drug$`RxNorm Ingredient id`,sep=":")
Drug = Drug %>%
  filter(!Treatment_GenericName%in%c('inebilizumab-cdon'))  # exclude biosimilar with incomplete data
druglist <- unique(Drug$feature_id)
myplot <- list()
cosinetable <- data.frame(cos = numeric(),
                          freq = numeric(),
                          year = numeric(),
                          drug = character(),
                          sd = numeric())
names(embed) <- as.character(2000:2022)
pdf("drug_trend.pdf", width = 6.5, height = 4.2)

# --- Individual drug plots ---
for(drug in druglist){
  # Extract cosine similarity between the drug and MS (PheCode:335) in each year's embedding.
  mycosine <- lapply(embed[as.character(yearlist)], function(x){
    if(drug%in%rownames(x) & "PheCode:335"%in%rownames(x)){
      x = x/apply(x, 1, norm, '2')   # L2-normalize before computing dot product
      cos = sum(x[which(rownames(x)==drug),]*x[which(rownames(x)=="PheCode:335"),])
      return(cos)
    }else{
      return(NA)
    }
  })
  freq = data %>% filter(feature_id==drug)
  freq = freq$freq[match(yearlist, freq$year)]
  mysd <- estedvar[match(yearlist, rownames(estedvar)),match(drug,colnames(estedvar))]
  mycosine <- data.frame(cos = unlist(mycosine),
                         freq = freq,
                         year = yearlist,
                         sd = mysd) %>%
    mutate(high_se = pmin(cos + sd*1.96,0.9),   # cap CI at [−0.1, 0.9] for readability
           low_se = pmax(cos - sd*1.96,-0.1))
  idx <- which(is.na(mycosine$cos))
  if(length(idx)>0){
    mycosine$high_se[idx] = NA
    mycosine$low_se[idx] = NA
  }
  cosinetable <- rbind(cosinetable, cbind(mycosine, drug = drug))
  nvalid = sum(!is.na(mycosine$cos))
  if(nvalid==0) next
  currmethod = ifelse(nvalid >= 4, "loess", "lm")
  curve_plot <- ggplot(mycosine, aes(x = year, y = cos)) +
    geom_point() +
    geom_smooth(method = currmethod, se = FALSE) +
    geom_errorbar(aes(ymin = low_se, ymax = high_se), color = "blue", alpha = 0.75) +
    scale_x_continuous(breaks=seq(2000, 2023, 2), limits = c(1999.5, 2023.5)) +
    scale_y_continuous(breaks=seq(-0.1, 0.8, 0.1), limits = c(-0.1, 0.87)) +
    ggtitle(paste(Drug$Treatment_GenericName[match(drug, Drug$feature_id)], ": ", drug, sep = ""))
  histogram_plot <- ggplot(mycosine, aes(x = year, y = freq)) +
    geom_bar(stat = "identity", fill = "gray", width = 0.5) +
    labs(x = "year", y = "Frequency") +
    scale_x_continuous(breaks=seq(2000, 2023, 2), limits = c(1999.5, 2023.5)) +
    theme_minimal()
  combined_plot <- curve_plot / histogram_plot +
    plot_layout(heights = c(2, 1))
  plot(combined_plot)
}

cosinetable$GName <- Drug$Treatment_GenericName[match(cosinetable$drug,Drug$feature_id)]

# --- Drug-class-level plots ---
# Each panel overlays all individual drugs in the class (colored) with the class-level
# embedding trend (black) and its empirical CI error bars.
for(drug in unique(Drug$Treatment_Class)){
  drugname = drug
  mycosine <- lapply(embed[as.character(yearlist)], function(x){
    if(drug%in%rownames(x) & "PheCode:335"%in%rownames(x)){
      x = x/apply(x, 1, norm, '2')
      cos = sum(x[which(rownames(x)==drug),]*x[which(rownames(x)=="PheCode:335"),])
      return(cos)
    }else{
      return(NA)
    }
  })
  freq = data %>% filter(feature_id==drug)
  freq = freq$freq[match(yearlist, freq$year)]
  mycosine <- data.frame(cos = unlist(mycosine),
                         freq = freq,
                         year = yearlist)
  mysd <- estedvar[match(yearlist, rownames(estedvar)),match(drug,colnames(estedvar))]
  cosinetable <- rbind(cosinetable,
                       cbind(mycosine,
                             sd = mysd,
                             drug = drugname,
                             GName = drugname) %>%
                         mutate(high_se = pmin(cos + sd*1.96,0.9),
                                low_se = pmax(cos - sd*1.96,-0.1)))
  idx <- which(is.na(cosinetable$cos))
  if(length(idx)>0){
    cosinetable$high_se[idx] = NA
    cosinetable$low_se[idx] = NA
  }
  # Combine the class-level entry with all its member drugs for the overlay plot.
  curr_druglist = c(drugname, sort(unique(Drug$feature_id[which(Drug$Treatment_Class==drug)])))
  mycosine <- cosinetable %>%
    filter(drug %in% curr_druglist)
  if(length(unique(mycosine$drug))<=2) next   # skip classes with only one drug
  unique_gnames <- sort(unique(mycosine$GName[mycosine$GName != drugname]))
  palette <- setNames(scales::hue_pal()(length(unique_gnames)), unique_gnames)
  palette[drugname] <- "black"   # class-level trend always plotted in black
  mycosine$alpha <- ifelse(mycosine$GName == drugname, 1, 0.7)
  mycosine$GName <- factor(mycosine$GName, levels = c(drugname, unique_gnames))

  filtered_mycosine <- mycosine %>% filter(!is.na(cos))

  # Split into loess vs lm groups based on number of valid years per drug.
  loess_data <- filtered_mycosine %>%
    group_by(GName) %>%
    filter(sum(!is.na(cos)) > 4)

  lm_data <- filtered_mycosine %>%
    group_by(GName) %>%
    filter(sum(!is.na(cos)) <= 4)

  curve_plot <- ggplot() +
    geom_point(data = filtered_mycosine %>% filter(GName != drugname),
               aes(x = year, y = cos, color = GName, alpha = alpha, shape = GName)) +
    geom_smooth(data = loess_data %>% filter(GName != drugname),
                aes(x = year, y = cos, color = GName),
                method = "loess", formula = 'y ~ x', se = FALSE) +
    geom_smooth(data = lm_data %>% filter(GName != drugname),
                aes(x = year, y = cos, color = GName),
                method = "lm", formula = 'y ~ x', se = FALSE) +
    geom_point(data = filtered_mycosine %>% filter(GName == drugname),
               aes(x = year, y = cos, alpha = alpha),
               color = "black", shape = 16) +
    geom_smooth(data = loess_data %>% filter(GName == drugname),
                aes(x = year, y = cos),
                method = "loess", formula = 'y ~ x', se = FALSE, color = "black") +
    geom_errorbar(data = loess_data %>% filter(GName == drugname),
                  aes(x = year, ymin = low_se, ymax = high_se), alpha = 0.75, color = "black") +
    scale_color_manual(values = palette) +
    scale_alpha_identity() +
    scale_shape_manual(values = 1:length(unique(filtered_mycosine$GName))) +
    scale_x_continuous(breaks = seq(2000, 2023, 2), limits = c(1999.5, 2023.5)) +
    scale_y_continuous(breaks = seq(-0.1, 0.8, 0.1), limits = c(-0.1, 0.87)) +
    theme(legend.position = "none") +
    ggtitle(drugname)
  histogram_plot <- ggplot(mycosine, aes(x = year, y = freq, fill = GName)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    scale_x_continuous(breaks=seq(2000, 2023, 2), limits = c(1999.5, 2023.5)) +
    scale_fill_manual(values = palette) +
    labs(x = "year", y = "Frequency") +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())
  if(length(unique(mycosine$GName))>=4){
    histogram_plot <- histogram_plot +
      guides(color = guide_legend(nrow = 2))
  }
  combined_plot <- curve_plot / histogram_plot +
    plot_layout(heights = c(2, 1))
  plot(combined_plot)
}
dev.off()

rm(embed)
rm(list = ls())

## get summary statistics ####
# Rebuild the cosine similarity table (without CIs) over the full 2000-2022 range,
# then classify each drug's temporal trend via linear regression and segmented regression.

library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(patchwork)
load("raw_embedding_SVD_SPPMI.Rdata")
load("rawdata_frequency_dict.Rdata")
Drug <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
yearlist = 2000:2022
Drug = Drug[!is.na(Drug$`RxNorm Ingredient id`),]
Drug$feature_id = paste("RXNORM", Drug$`RxNorm Ingredient id`,sep=":")
Drug = Drug %>%
  filter(!Treatment_GenericName%in%c('inebilizumab-cdon'))
druglist <- unique(Drug$feature_id)
myplot <- list()
cosinetable <- data.frame(cos = numeric(),
                          freq = numeric(),
                          year = numeric(),
                          drug = character())
for(drug in druglist){
  mycosine <- lapply(embed, function(x){
    if(drug%in%rownames(x) & "PheCode:335"%in%rownames(x)){
      x = x/apply(x, 1, norm, '2')
      cos = sum(x[which(rownames(x)==drug),]*x[which(rownames(x)=="PheCode:335"),])
      return(cos)
    }else{
      return(NA)
    }
  })
  freq = data %>% filter(feature_id==drug)
  freq = freq$freq[match(yearlist, freq$year)]
  mycosine <- data.frame(cos = unlist(mycosine),
                         freq = freq,
                         year = yearlist)
  cosinetable <- rbind(cosinetable, cbind(mycosine, drug = drug))
  nvalid = sum(!is.na(mycosine$cos))
  if(nvalid==0) next
}

cosinetable$GName <- Drug$Treatment_GenericName[match(cosinetable$drug,Drug$feature_id)]

for(drug in unique(Drug$Treatment_Class)){
  drugname = drug
  mycosine <- lapply(embed, function(x){
    if(drug%in%rownames(x) & "PheCode:335"%in%rownames(x)){
      x = x/apply(x, 1, norm, '2')
      cos = sum(x[which(rownames(x)==drug),]*x[which(rownames(x)=="PheCode:335"),])
      return(cos)
    }else{
      return(NA)
    }
  })
  freq = data %>% filter(feature_id==drug)
  freq = freq$freq[match(yearlist, freq$year)]
  mycosine <- data.frame(cos = unlist(mycosine),
                         freq = freq,
                         year = yearlist)
  cosinetable <- rbind(cosinetable, cbind(mycosine, drug = drugname, GName = drugname))
}

# Simple linear trend test: fit cos ~ year and record the slope, t-statistic, and p-value.
summarytable <- data.frame(t = numeric(),
                           p = numeric(),
                           beta = numeric(),
                           se_beta = numeric(),
                           drug = character(),
                           n = numeric())
for(currdrug in unique(cosinetable$drug)){
  mycosine <- cosinetable %>%
    filter(drug == currdrug) %>%
    select(cos, year) %>%
    na.omit()
  if(nrow(mycosine) > 2){
    lmmod <- lm(cos~year, mycosine)
    slmmod <- summary(lmmod)
    summarytable <- rbind(summarytable,
                          data.frame(t = slmmod$coefficients['year','t value'],
                                     p = slmmod$coefficients['year','Pr(>|t|)'],
                                     beta = slmmod$coefficients['year','Estimate'],
                                     se_beta = slmmod$coefficients['year','Std. Error'],
                                     drug = currdrug,
                                     n = nrow(mycosine)))
  }else{
    summarytable <- rbind(summarytable,
                          data.frame(t = NA,
                                     p = NA,
                                     beta = NA,
                                     se_beta = NA,
                                     drug = currdrug,
                                     n = nrow(mycosine)))
  }
}

# Classify drugs by trend direction based on p < 0.05 from the linear model.
increase_drug = summarytable$drug[which(summarytable$p < 0.05 & summarytable$beta > 0)]
decrease_drug = summarytable$drug[which(summarytable$p < 0.05 & summarytable$beta < 0)]
random_drug = setdiff(summarytable$drug, c(increase_drug,decrease_drug))

summarytable$trend = "None"
summarytable$trend[which(summarytable$drug%in%increase_drug)] = "Increase"
summarytable$trend[which(summarytable$drug%in%decrease_drug)] = "Decrease"

# Plot normalized cosine similarity (z-scored within each drug) grouped by trend type,
# followed by individual drug plots for each trend group.
pdf("drug.pdf")
for(idx in 1:3){
  if(idx==1){
    currdrug = increase_drug
    mytitle = "increase"
  }else if(idx==2){
    currdrug = decrease_drug
    mytitle = "decrease"
  }else{
    currdrug = random_drug
    mytitle = "random"
  }
  currtable <- cosinetable %>%
    filter(drug %in% currdrug) %>%
    na.omit() %>%
    group_by(drug) %>%
    mutate(mean_cos = mean(cos),
           sd_cos = sd(cos)) %>%
    ungroup() %>%
    mutate(norm_cos = (cos-mean_cos)/sd_cos)   # z-score to put drugs on a common scale
  plot(ggplot(currtable, aes(x = year, y = norm_cos, group = drug, color = drug)) +
         geom_point() +
         geom_smooth(se = FALSE) +
         ggtitle(mytitle))
  for(mydrug in currdrug){
    currtable <- cosinetable %>%
      filter(drug == mydrug) %>%
      na.omit()
    if(nrow(currtable) <= 3) next
    plot(ggplot(currtable, aes(x = year, y = cos, group = drug, color = drug)) +
           geom_point() +
           geom_smooth(se = FALSE) +
           ggtitle(paste(mytitle, ":", mydrug)))
  }
}
dev.off()

library(segmented)

# analyze_trend: classify a drug's cosine-similarity time series into one of five trend types:
#   "Increase", "Decrease", "None", "I-D" (increase then decrease), "D-I" (decrease then increase).
# Strategy:
#   1. Fit a simple linear model (lm).
#   2. If >= 6 data points, also attempt a one-breakpoint segmented regression.
#   3. If segmented fit succeeds and both segments are significant, classify by their sign combination.
#   4. Otherwise fall back to the global linear trend.
# Returns slopes (beta, beta1, beta2), p-values, half-widths (w, w1, w2), and the change point year.
# NOTE: The first definition below is superseded by the second definition; the second one is used.
analyze_trend <- function(df) {
  if(nrow(df) <= 3) return(list(type = "None", change_point = NA,
                                beta1 = NA, beta2 = NA, p1 = NA, p2 = NA,
                                w1 = NA, w2 = NA))

  lm_fit <- lm(cos ~ year, data = df)

  seg_fit <- tryCatch(
    segmented(lm_fit, seg.Z = ~ year, npsi = 1),
    error = function(e) NULL
  )

  if (is.null(seg_fit) | nrow(df) < 6) {
    summary_lm <- summary(lm_fit)
    slope <- coef(summary_lm)["year", "Estimate"]
    beta1 <- slope
    beta2 <- NA
    pval <- coef(summary_lm)["year", "Pr(>|t|)"]
    w1 <- coef(summary_lm)["year", "Std. Error"] * qt(1-alpha/2, summary_lm$df[2])
    w2 <- NA

    if (pval >= 0.05) {
      return(list(type = "None", change_point = NA,
                  beta1 = beta1, beta2 = beta2, p1 = pval, p2 = NA,
                  w1 = w1, w2 = w2))
    } else if (slope > 0) {
      return(list(type = "Increase", change_point = NA,
                  beta1 = beta1, beta2 = beta2, p1 = pval, p2 = NA,
                  w1 = w1, w2 = w2))
    } else {
      return(list(type = "Decrease", change_point = NA,
                  beta1 = beta1, beta2 = beta2, p1 = pval, p2 = NA,
                  w1 = w1, w2 = w2))
    }
  } else {
    slopes <- slope(seg_fit)$year[, 1]
    slope_se <- slope(seg_fit)$year[, 2]
    slope_se[is.nan(slope_se)] <- Inf
    slope_t <- slopes / slope_se
    slope_p <- 2 * pt(-abs(slope_t), df = df.residual(seg_fit))
    w1 = slope_se[1] * qt(1-alpha/2, df = df.residual(seg_fit)[1])
    w2 = slope_se[2] * qt(1-alpha/2, df = df.residual(seg_fit)[2])

    cp <- seg_fit$psi[2]

    beta1 <- slopes[1]
    beta2 <- slopes[2]
    p1 <- slope_p[1]
    p2 <- slope_p[2]

    sig1 <- slope_p[1] < 0.05
    sig2 <- slope_p[2] < 0.05

    if (!sig1 && !sig2) {
      summary_lm <- summary(lm_fit)
      slope <- coef(summary_lm)["year", "Estimate"]
      pval <- coef(summary_lm)["year", "Pr(>|t|)"]
      w1 <- coef(summary_lm)["year", "Std. Error"] * qt(1-alpha/2, summary_lm$df[2])
      w2 <- NA
      if (pval >= 0.05) {
        return(list(type = "None", change_point = NA,
                    beta1 = slope, beta2 = NA, p1 = pval, p2 = NA,
                    w1 = w1, w2 = w2))
      } else if (slope > 0) {
        return(list(type = "Increase", change_point = NA,
                    beta1 = slope, beta2 = NA, p1 = pval, p2 = NA,
                    w1 = w1, w2 = w2))
      } else {
        return(list(type = "Decrease", change_point = NA,
                    beta1 = slope, beta2 = NA, p1 = pval, p2 = NA,
                    w1 = w1, w2 = w2))
      }
    } else if (sig1 && sig2) {
      if (slopes[1] > 0 && slopes[2] < 0) {
        type <- "I-D"
      } else if (slopes[1] < 0 && slopes[2] > 0) {
        type <- "D-I"
      } else if (slopes[1] > 0 && slopes[2] > 0) {
        type <- "Increase"
      } else if (slopes[1] < 0 && slopes[2] < 0) {
        type <- "Decrease"
      } else {
        type <- "None"
      }
    } else if (sig1) {
      type <- ifelse(slopes[1] > 0, "Increase", "Decrease")
    } else if (sig2) {
      type <- ifelse(slopes[2] > 0, "Increase", "Decrease")
    } else {
      type <- "None"
    }

    if(type %in% c("Decrease", "Increase", "None")){
      summary_lm <- summary(lm_fit)
      beta1 <- coef(summary_lm)["year", "Estimate"]
      beta2 <- NA
      p1 <- coef(summary_lm)["year", "Pr(>|t|)"]
      p2 <- NA
      w1 <- coef(summary_lm)["year", "Std. Error"] * qt(1-alpha/2, summary_lm$df[2])
      w2 <- NA
    }

    return(list(type = type, change_point = cp,
                beta1 = beta1, beta2 = beta2,
                p1 = p1, p2 = p2, w1 = w1, w2 = w2))
  }
}

# Revised version of analyze_trend: also returns global (unsegmented) beta/p/w
# alongside segment-specific values, so downstream tables have a single beta column.
analyze_trend <- function(df) {
  if(nrow(df) <= 3) return(list(type = "None", change_point = NA,
                                beta = NA, beta1 = NA, beta2 = NA,
                                p = NA, p1 = NA, p2 = NA,
                                w = NA, w1 = NA, w2 = NA))

  lm_fit <- lm(cos ~ year, data = df)

  seg_fit <- tryCatch(
    segmented(lm_fit, seg.Z = ~ year, npsi = 1),
    error = function(e) NULL
  )

  # Always compute global linear slope as the primary summary statistic.
  summary_lm <- summary(lm_fit)
  slope <- coef(summary_lm)["year", "Estimate"]
  beta <- slope
  p <- coef(summary_lm)["year", "Pr(>|t|)"]
  w <- coef(summary_lm)["year", "Std. Error"] * qt(1-alpha/2, summary_lm$df[2])

  if ((!is.null(seg_fit)) & nrow(df) >= 6) {
    slopes <- slope(seg_fit)$year[, 1]
    slope_se <- slope(seg_fit)$year[, 2]
    slope_se[is.nan(slope_se)] <- Inf
    slope_t <- slopes / slope_se
    slope_p <- 2 * pt(-abs(slope_t), df = df.residual(seg_fit))
    w1 = slope_se[1] * qt(1-alpha/2, df = df.residual(seg_fit))
    w2 = slope_se[2] * qt(1-alpha/2, df = df.residual(seg_fit))

    cp <- seg_fit$psi[2]   # estimated change-point year

    beta1 <- slopes[1]
    beta2 <- slopes[2]
    p1 <- slope_p[1]
    p2 <- slope_p[2]

    sig1 <- slope_p[1] < alpha
    sig2 <- slope_p[2] < alpha

    if (!sig1 && !sig2){
      # Neither segment significant: fall back to global linear classification.
      if (p > alpha) {
        type = "None"
      } else if (slope > 0) {
        type = "Increase"
      } else {
        type = "Decrease"
      }
    } else if (sig1 && sig2) {
      if (slopes[1] > 0 && slopes[2] < 0) {
        type <- "I-D"
      } else if (slopes[1] < 0 && slopes[2] > 0) {
        type <- "D-I"
      } else if (slopes[1] > 0 && slopes[2] > 0) {
        type <- "Increase"
      } else if (slopes[1] < 0 && slopes[2] < 0) {
        type <- "Decrease"
      } else {
        type <- "None"
      }
    } else if (sig1) {
      type <- ifelse(slopes[1] > 0, "Increase", "Decrease")
    } else if (sig2) {
      type <- ifelse(slopes[2] > 0, "Increase", "Decrease")
    } else {
      type <- "None"
    }
  }else{
    # Segmented fit failed or too few points: use global linear model only.
    if (p > alpha) {
      type = "None"
    } else if (slope > 0) {
      type = "Increase"
    } else {
      type = "Decrease"
    }
    cp = p1 = p2 = w1 = w2 = beta1 = beta2 = NA
  }

  return(list(type = type, change_point = cp,
              beta = beta, beta1 = beta1, beta2 = beta2,
              p = p, p1 = p1, p2 = p2,
              w = w, w1 = w1, w2 = w2))

}

alpha = 0.05
results <- data.frame(drug = character(),
                      type = character(),
                      change_point = numeric(),
                      n = numeric(),
                      start = numeric(),
                      end = numeric(),
                      beta = numeric(),
                      beta1 = numeric(),
                      beta2 = numeric(),
                      p = numeric(),
                      p1 = numeric(),
                      p2 = numeric(),
                      w = numeric(),
                      w1 = numeric(),
                      w2 = numeric())
for(drug_name in unique(cosinetable$drug)){
  df <- cosinetable %>% filter(drug == drug_name) %>% na.omit()
  res <- analyze_trend(df)
  results <- rbind(results,
                   data.frame(drug = drug_name,
                              type = res$type,
                              change_point = res$change_point,
                              n = nrow(df),
                              start = min(df$year),
                              end = max(df$year),
                              beta = res$beta,
                              beta1 = res$beta1,
                              beta2 = res$beta2,
                              p = res$p,
                              p1 = res$p1,
                              p2 = res$p2,
                              w = res$w,
                              w1 = res$w1,
                              w2 = res$w2))
}


fullsummarytable <- merge(summarytable, results)

# Plot normalized cosine similarity by trend type, with a vertical change-point line
# for non-monotone trends (I-D or D-I).
pdf("drug1.pdf", width = 8, height = 7)
for(currtype in unique(results$type)){
  currdrug <- results$drug[which(results$type==currtype)]
  currtable <- cosinetable %>%
    filter(drug %in% currdrug) %>%
    na.omit() %>%
    group_by(drug) %>%
    mutate(mean_cos = mean(cos),
           sd_cos = sd(cos)) %>%
    ungroup() %>%
    mutate(norm_cos = (cos-mean_cos)/sd_cos)
  plot(ggplot(currtable, aes(x = year, y = norm_cos, group = drug, color = drug)) +
         geom_point() +
         geom_smooth(se = FALSE) +
         ggtitle(paste(currtype, ":", nrow(currtype))))
  for(mydrug in currdrug){
    if(currtype%in%c("D-I","I-D")){
      cp <- results$change_point[which(results$drug==mydrug)]
    }
    currtable <- cosinetable %>%
      filter(drug == mydrug) %>%
      na.omit()
    if(nrow(currtable) <= 3) next

    p <- ggplot(currtable, aes(x = year, y = cos, group = drug, color = drug)) +
      geom_point() +
      ggtitle(paste(currtype, ":", mydrug))
    if(nrow(currtable) < 6){
      p <- p +
        geom_smooth(method = "lm", se = FALSE)
    }else{
      p <- p +
        geom_smooth(method = "loess", se = FALSE)
    }
    if(currtype%in%c("D-I","I-D")){
      p <- p + geom_vline(xintercept = cp, linetype = "dashed", color = "red", size = 1)
    }
    plot(p)
  }
}
dev.off()

# Merge trend classification results back into the drug metadata table and export.
Drugtable <- read_excel("MS Drugs Mechanism Citations Efficacy 20240927.xlsx")
Drugtable$drug <- paste("RXNORM", Drugtable$`RxNorm Ingredient id`, sep = ":")
Drugtable <- merge(Drugtable, results, all = TRUE)
write.csv(Drugtable, file = "MS Drugs Mechanism Citations Efficacy with Trend group and Beta.csv")
