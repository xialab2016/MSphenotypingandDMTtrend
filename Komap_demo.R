# Title: plot
# Output target in Rmd (html_document) is removed for .R use
# Date: 2025-04-30

# =========================
# Setup
# =========================

# Global chunk options removed (knitr). Echo is naturally controlled in .R via console.

# --- Packages ---
library(data.table)
library(stringr)
library(lubridate)   # Date handling
library(readr)
library(tableone)
require(openxlsx)
library(naniar)
library(writexl)
library(dplyr)
library(ggplot2)
library(caTools)
library(corrplot)
library(patchwork)
library(gplots)
library(gridExtra)
library(tidyr)
library(readxl)
library(conflicted)

# Resolve common function name conflicts in favor of dplyr
conflict_prefer("select", "dplyr")
conflict_prefer("first", "dplyr")
conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::slice)

# =========================
# Data load (initial demo)
# =========================

# Using here() was in the Rmd, but files are referenced via explicit paths or CWD.
# If needed, uncomment the next line and configure a project root.
# library(here)

demo_MGB <- read_csv("demo_MGB.csv")
MGB_MS_algorithm_diagnosis_20240201 <- read_csv("MGB_MS_algorithm_diagnosis 20240201.csv")

# EHR mapping can be read from CSV or XLSX; both lines were in the Rmd.
# Keep both; the latter overwrites the former if both exist in your environment.
EHR_mapping <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/20240202/EHR_mapping.csv")
EHR_mapping <- read_excel("EHR_mapping.xlsx")

# Quick check for missing PatientNum
EHR_mapping %>% filter(is.na(PatientNum))

# Build MGB cohort (MS specificity 0.90), join mappings and demographics, derive NHW flag, keep records with non-missing SEX
MGB <- MGB_MS_algorithm_diagnosis_20240201 %>%
  filter(MS.spec.90 == 1) %>%
  select(patient_num, is.registry) %>%
  rename(PatientNum = patient_num) %>%
  left_join(EHR_mapping, by = "PatientNum") %>%
  left_join(demo_MGB, by = c("MAPPING_INDEX" = "PATIENT_NUM")) %>%
  mutate(NHW = ifelse(RACE == 7 & ETHNICITY == 2, 1, 0)) %>%
  filter(!is.na(SEX))

# =========================
# TableOne example (note: object `summary` is not defined in the source Rmd)
# =========================
# The original Rmd used `summary` here, which is not defined.
# This code is kept as-is from the Rmd; ensure `summary` exists in your environment
# or replace `summary` with the intended data frame.

var <- c("AGE_NOW","SEX","Disease_Subtype", "THNICITY_DESC","RACE_DESC")
catVar <- c("SEX","THNICITY_DESC","RACE_DESC","Disease_Subtype")

# CreateTableOne on `summary` (as in the original Rmd)
summary_tbl  <- CreateTableOne(vars = var, data = summary, factorVars = catVar)
tab3Mat  <- print(summary_tbl)

DemotableUniqueN <- summary %>%
  group_by(PATIENT_NUM) %>%
  arrange(PATIENT_NUM, date) %>%
  slice(1) %>%
  arrange(PATIENT_NUM)

summary_tbl2 <- CreateTableOne(vars = var, data = DemotableUniqueN, factorVars = catVar)
tab3Mat2 <- as.data.frame(print(summary_tbl2, exact = "stage", smd = FALSE))

lis <- list("MGB unique patient" = tab3Mat2, "MGB PDDS" = tab3Mat)

# =========================
# chat - demo
# =========================

library(dplyr)
library(lubridate)
library(tableone)
library(openxlsx)

# Step 1: Construct MGB data
MGB <- MGB_MS_algorithm_diagnosis_20240201 %>%
  filter(MS.spec.90 == 1) %>%
  select(patient_num, is.registry) %>%
  rename(PatientNum = patient_num) %>%
  left_join(EHR_mapping, by = "PatientNum") %>%
  left_join(demo_MGB, by = c("MAPPING_INDEX" = "PATIENT_NUM")) %>%
  mutate(
    NHW = ifelse(RACE == 7 & ETHNICITY == 2, 1, 0),
    # Step 2: Build diagnosis date fields (fallback to first symptom if diagnosis year is NA)
    DIAG_YEAR  = ifelse(!is.na(FIRSTDIAG_YEAR),  FIRSTDIAG_YEAR,  FIRSTSYMP_YEAR),
    DIAG_MONTH = ifelse(!is.na(FIRSTDIAG_YEAR),  FIRSTDIAG_MONTH, 7),
    DIAG_DAY   = ifelse(!is.na(FIRSTDIAG_YEAR),  FIRSTDIAG_DAY,   15),
    # If year is NA, set DIAG_DATE to NA to avoid ymd() errors
    DIAG_DATE = case_when(
      is.na(DIAG_YEAR) ~ as.Date(NA),
      TRUE ~ ymd(sprintf("%04d-%02d-%02d",
                         DIAG_YEAR, coalesce(DIAG_MONTH, 7), coalesce(DIAG_DAY, 15)))
    ),
    DISEASE_DURATION = as.numeric(difftime(ymd("2024-02-01"), DIAG_DATE, units = "days")) / 365.25
  ) %>%
  filter(!is.na(SEX))

# Step 3: Select variables and generate TableOne
var <- c("AGE_NOW", "SEX", "Disease_Subtype", "THNICITY_DESC", "RACE_DESC", "DISEASE_DURATION")
catVar <- c("SEX", "THNICITY_DESC", "RACE_DESC", "Disease_Subtype")

summary_tbl <- CreateTableOne(vars = var, data = MGB, factorVars = catVar)
tab3Mat <- print(summary_tbl, printToggle = FALSE)

# Step 4: Export to Excel
write.xlsx(tab3Mat, file = "MGB_summary_table.xlsx", rowNames = TRUE)
getwd()

# =========================
# MGB demo
# =========================

setwd("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh")

# 18,868
MS_algorithm_diagnosis <- read_csv("Gennetech/E4A/MGB/Diagnosis Phenotype/20240201/MS_algorithm_diagnosis.csv") %>%
  filter(MS.spec.90 == 1) %>%
  select(patient_num)

EHR_mapping <- read_csv("20240202/EHR_mapping.csv")
MS_demographics <- read_csv("20240131/MS_demographics.csv")
MGB_demographics_imputed_race_ethnicity_separate <- read_csv("Raw Data/Data Queries/NIWreport/MGB demographics imputed race ethnicity separate.csv")

table(MS_demographics$Race1)

mapping <- EHR_mapping %>% filter(EHR_mapping$PatientNum %in% MS_algorithm_diagnosis$patient_num)

a <- MGB_MS_algorithm_diagnosis_20240201 %>%
  rename(PatientNum = patient_num) %>%
  left_join(MS_demographics, by = "PatientNum") %>%
  select(PatientNum, MS.spec.90,
         Gender_Legal_Sex, Race_Group, Ethnic_Group, Date_of_Birth) %>%
  mutate(
    Ethnic_Group = case_when(
      is.na(Ethnic_Group) | Ethnic_Group %in% c("DECLINED") ~ "Unknown/Missing",
      TRUE ~ Ethnic_Group
    ),
    Race_Group = case_when(
      is.na(Race_Group) | Race_Group %in% c("Other", "Declined") ~ "Unknown/Missing",
      TRUE ~ Race_Group
    )
  ) %>%
  left_join(
    select(MGB_demographics_imputed_race_ethnicity_separate, PatientNum, race, ethnicity),
    by = "PatientNum", suffix = c("", "_imputed")
  ) %>%
  mutate(
    Ethnic_Group = case_when(
      is.na(Ethnic_Group) | Ethnic_Group %in% c("Unknown", "Unknown/Missing") ~ ethnicity,
      TRUE ~ Ethnic_Group
    ),
    Race_Group = case_when(
      is.na(Race_Group) | Race_Group %in% c("Unknown", "Unknown/Missing") ~ race,
      TRUE ~ Race_Group
    ),
    Gender_Legal_Sex = case_when(Gender_Legal_Sex == "O" ~ NA, TRUE ~ Gender_Legal_Sex),
    Race_Group = case_when(Race_Group == "Other" ~ NA, TRUE ~ Race_Group)
  ) %>%
  filter(MS.spec.90 == 1)

# =========================
# Read EHR data for: age at 1st PheCode:335
# =========================

setwd("/Users/kay/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Raw Data")
target <- read_excel("MS Drugs FDA app.xlsx") %>% mutate(feature_id = paste("RXNORM:", `RxNorm Ingredient id`, sep = ""))
rx <- target %>% select(Treatment_BrandName, Treatment_Class, feature_id)

# Initialize container for daily data merge
all_data <- data.frame()

# Generate file names
file_names <- c()
start <- seq(1, 43001, by = 1000)
end   <- seq(1000, 44000, by = 1000)

setwd("/Users/kay/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Gennetech/E4A/MGB/20240613/daily")

for (i in 1:length(start)) {
  file_name <- sprintf("daily_patient_%d-%d.RData", start[i], end[i])
  file_names <- c(file_names, file_name)
}

for (file_name in file_names) {
  if (file.exists(file_name)) {
    print(file_name)
    load(file_name)
    # Expecting `ms.daily` to be present after load()
    ms.daily <- ms.daily %>% filter(feature_id %in% rx$feature_id | feature_id == "PheCode:335")
    all_data <- rbind(all_data, ms.daily)
    rm(ms.daily)
  } else {
    print(paste("File not found:", file_name))
  }
}

# =========================
# First diagnosis date and DMT flags
# =========================

mgb_date <- all_data %>%
  group_by(PatientNum, feature_id) %>%
  arrange(Date) %>%
  slice(1)

dig <- mgb_date %>%
  filter(feature_id == "PheCode:335") %>%
  select(PatientNum, Date, feature_id) %>%
  mutate(DISEASE_DURATION = as.numeric(difftime(ymd("2023-01-31"), Date, units = "days")) / 365.25)

library(dplyr)
library(stringr)

# Step 1: Mark which records are DMT (RXNORM)
mgb_date <- mgb_date %>%
  mutate(is_rxnorm = str_detect(feature_id, "RXNORM"))

# Step 2: Build ever_dmt per patient
ever_dmt_df <- mgb_date %>%
  group_by(PatientNum) %>%
  summarise(ever_dmt = ifelse(any(is_rxnorm), "yes", "no"), .groups = "drop")

# Combine with demographics (object `a` from above), compute AGE_AT_DIAG
mgbSummary <- a %>%
  left_join(dig, by = "PatientNum") %>%
  mutate(
    Date_of_Birth = as.Date(Date_of_Birth),
    Date = as.Date(Date),
    AGE_AT_DIAG = as.numeric(difftime(Date, Date_of_Birth, units = "days")) / 365.25
  ) %>%
  select(-feature_id) %>%
  arrange(PatientNum) %>%
  mutate(RACE_TITLE = toupper(Race_Group)) %>%
  mutate(
    RACE_TITLE = case_when(
      Race_Group %in% c("ALASKA NATIVE", "AMERICAN INDIAN/ALASKA NATIVE") ~ "American Indian/Alask  a Native",
      Race_Group %in% c("CHINESE", "ASIAN INDIAN", "FILIPINO", "JAPANESE", "KOREAN", "OTHER ASIAN") ~ "Asian",
      Race_Group %in% c("NATIVE HAWAIIAN", "OTHER PACIFIC ISLANDER") ~ "Native Hawaiian or Other Pacific Islander",
      is.na(Race_Group) | Race_Group %in% c("OTHER", "UNREPORTED", "CHOSE NOT TO DISCLOSE RACE", "Unknown or Not Reported",
                                            "UNREPORTED,CHOSE NOT TO DISCLOSE RACE", "NA") ~ "Unknown or Not Reported",
      TRUE ~ Race_Group
    )
  )
# Optionally join `ever_dmt_df` if desired:
# %>% left_join(ever_dmt_df, by = "PatientNum")

# =========================
# TableOne for MGB summary
# =========================

dput(names(mgbSummary))
library(tableone)
library(openxlsx)

# 1) Variables
vars <- c("AGE_AT_DIAG","Gender_Legal_Sex", "race", "ethnicity", "ever_dmt","DISEASE_DURATION")

# 2) Categorical variables
catVars <- c("race","Gender_Legal_Sex", "ethnicity", "ever_dmt")

# 3) Create TableOne
summary_tbl <- CreateTableOne(vars = vars, data = mgbSummary, factorVars = catVars)

# 4) Matrix for Excel
tabMat <- print(summary_tbl, printToggle = FALSE)

# 5) Save Excel
write.xlsx(tabMat, file = "mgb_summary_table_full.xlsx", rowNames = TRUE)
summary(mgbSummary)

# =========================
# Pitt demo
# =========================

library(readr)

UPMC_MS_note3_patient_status_90_specificity <- read_csv("UPMC_MS_note3_patient_status_90_specificity.csv")
# ClinicalDemographics (older path commented in original)
ClinicalDemographics2025_01_27 <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Raw Data/ClinicalDemographics2025-01-27.csv") %>%
  select(date_msdx, id_participant)

R3_2468_Xia_PROMOTE <- read_excel("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Gennetech/E4A/ICD/R3_2468_Xia_PROMOTE_20220223_studyid_NoPHI.xls") %>%
  select(-LINE_NUM) %>%
  mutate(PATIENT_NUM = as.numeric(PATIENT_NUM))

dx <- ClinicalDemographics2025_01_27 %>%
  left_join(R3_2468_Xia_PROMOTE, by = c("id_participant" = "PATIENT_NUM")) %>%
  rename(PATIENT_STUDY_ID = PATIENT_NUM)

ms_demo <- read.delim(
  "~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Raw Data/MS_Demographics.dsv",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

pitt <- UPMC_MS_note3_patient_status_90_specificity %>%
  select(patient_num, is.Registry, pred_y) %>%
  rename(PATIENT_STUDY_ID = patient_num) %>%
  left_join(ms_demo, by = "PATIENT_STUDY_ID") %>%
  filter(pred_y == 1) %>%
  left_join(dx, by = "PATIENT_STUDY_ID")

dput(names(pitt))

# =========================
# Pitt codified data and cleanup
# =========================

setwd("/Users/kay/Library/CloudStorage/OneDrive-UniversityofPittsburgh")

UPMC_MS_2011_to_2021_Codified_rolled_up_data_2023_02_04 <- read_csv(
  "Gennetech/E4A/ICD/UPMC_MS_2011_to_2021_Codified_rolled_up_data_2023-02-04.csv"
) %>%
  mutate(PATIENT_NUM = as.character(patient_num)) %>%
  select(-patient_num) %>%
  filter(feature_id %in% rx$feature_id | feature_id == "PheCode:335")

Pitt <- pitt %>%
  mutate(
    RACE_TITLE = case_when(
      RACE_TITLE %in% c("NOT SPECIFIED", "UNREPORTED", "CHOSE NOT TO DISCLOSE RACE") ~ "Unknown or Not Reported",
      TRUE ~ RACE_TITLE
    ),
    RACE_TITLE = case_when(
      RACE_TITLE %in% c("ALASKA NATIVE", "AMERICAN INDIAN/ALASKA NATIVE") ~ "American Indian/Alask  a Native",
      RACE_TITLE %in% c("CHINESE", "ASIAN INDIAN", "FILIPINO", "JAPANESE", "KOREAN", "OTHER ASIAN") ~ "Asian",
      RACE_TITLE %in% c("NATIVE HAWAIIAN", "OTHER PACIFIC ISLANDER") ~ "Native Hawaiian or Other Pacific Islander",
      is.na(RACE_TITLE) | RACE_TITLE %in% c("OTHER", "UNREPORTED", "CHOSE NOT TO DISCLOSE RACE", "Unknown or Not Reported",
                                            "UNREPORTED,CHOSE NOT TO DISCLOSE RACE", "NA") ~ "Unknown or Not Reported",
      TRUE ~ RACE_TITLE
    ),
    ETHNIC_TITLE = case_when(
      is.na(ETHNIC_TITLE) | ETHNIC_TITLE %in% c("NOT SPECIFIED", "UNREPORTED/CHOSE NOT TO DISCLOSE", "UNREPORTED", "CHOSE NOT TO DISCLOSE") ~ "unknown",
      TRUE ~ ETHNIC_TITLE
    )
  ) %>%
  select(PATIENT_STUDY_ID, BIRTH_DATE, GENDER_TITLE, RACE_TITLE, ETHNIC_TITLE)

dput(names(Pitt))

# =========================
# Disease duration (Pitt EHR chunked files)
# =========================

setwd("/Users/kay/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Raw Data/EHR/20250129 EHR Data 2025 Export/MSObsStructured200101282025")

file_names <- c(
  "MSObsStructured200101282025.csv",
  "MSObsStructured400101282025.csv",
  "MSObsStructured600101282025.csv",
  "MSObsStructured800101282025.csv",
  "MSObsStructured1000101282025.csv",
  "MSObsStructured1200101282025.csv",
  "MSObsStructured1400101282025.csv",
  "MSObsStructured1600101282025.csv",
  "MSObsStructured1800101282025.csv",
  "MSObsStructured2000001282025.csv"
)

phe_rx_df <- c()
for (file in file_names) {
  # Read each CSV; chunked reading not implemented here—adjust if files are very large.
  getwd()
  df <- read.csv(file)
  rxnorm <- df %>% filter(feature_id == "PheCode:335")
  phe_rx_df <- bind_rows(phe_rx_df, rxnorm)
}

phe_rx_df <- phe_rx_df %>% rename(PATIENT_NUM = patient_num)

# =========================
# First PheCode:335 and ever DMT for Pitt
# =========================

library(dplyr)
library(lubridate)

# Step 1: Sort and cast dates
test335 <- phe_rx_df %>%
  arrange(PATIENT_NUM, start_date) %>%
  filter(feature_id == "PheCode:335") %>%
  mutate(start_date = as.Date(start_date))

summary(test335$start_date)

# Step 2: Extract first PheCode:335 date per patient
phe335_date <- phe_rx_df %>%
  group_by(PATIENT_NUM) %>%
  summarise(first_phe335_date = min(as.Date(start_date)), .groups = "drop")

# Step 3: Ever DMT flag (any RXNORM MS drug)
ever_dmt_df <- phe_rx_df %>%
  filter(feature_id %in% rx$feature_id) %>%
  distinct(PATIENT_NUM) %>%
  mutate(ever_dmt = "yes")

patients <- unique(phe_rx_df$PATIENT_NUM)
ever_dmt_dfP <- tibble(PATIENT_NUM = patients) %>%
  left_join(ever_dmt_df, by = "PATIENT_NUM") %>%
  mutate(ever_dmt = ifelse(is.na(ever_dmt), "no", ever_dmt))

# Step 4: Combine first Phe335 and ever_dmt
summary_df <- full_join(phe335_date, ever_dmt_dfP, by = "PATIENT_NUM")

# Step 5: Merge demographics and compute ages/durations
summary_df <- summary_df %>%
  mutate(PATIENT_NUM = as.numeric(PATIENT_NUM)) %>%
  right_join(Pitt %>% rename(PATIENT_NUM = PATIENT_STUDY_ID), by = "PATIENT_NUM") %>%
  mutate(
    BIRTH_DATE = dmy(BIRTH_DATE),
    AGE_AT_DIAG = as.numeric(difftime(first_phe335_date, BIRTH_DATE, units = "days")) / 365.25,
    DISEASE_DURATION = as.numeric(difftime(ymd("2022-11-30"), first_phe335_date, units = "days")) / 365.25
  )

vars <- c("AGE_AT_DIAG", "DISEASE_DURATION", "GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE", "ever_dmt")
catVars <- c("GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE", "ever_dmt")

# Step 6: TableOne for Pitt
summary_tbl <- CreateTableOne(vars = vars, data = summary_df, factorVars = catVars)
tabMat <- print(summary_tbl, printToggle = FALSE)

# Quick overview
summary(summary_df)

# Optional: write Excel
# write.xlsx(tabMat, file = "pitt_summary.xlsx", rowNames = TRUE)

# =========================
# Combine MGB and Pitt for analysis
# =========================

a1 <- summary_df %>%
  rename(PatientNum = PATIENT_NUM) %>%
  select(PatientNum, DISEASE_DURATION, AGE_AT_DIAG) %>%
  mutate(group = "Pitt")

a2 <- mgbSummary %>%
  select(PatientNum, DISEASE_DURATION, AGE_AT_DIAG) %>%
  mutate(group = "MGB")

b <- rbind(a1, a2)

vars <- c("AGE_AT_DIAG", "DISEASE_DURATION", "GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE", "ever_dmt")
catVars <- c("GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE", "ever_dmt")

summary_tbl <- CreateTableOne(vars = vars, data = b, factorVars = catVars)
tabMat <- print(summary_tbl, printToggle = FALSE)
summary(b)

# For reference
names(summary_df)
names(mgbSummary)

# =========================
# Pitt: age now and disease duration (alt calc)
# =========================

library(dplyr)
library(lubridate)
library(tableone)
library(openxlsx)

# Step 1: Compute current age and duration relative to 2024-02-01
pitt <- pitt %>%
  mutate(
    BIRTH_DATE = dmy(BIRTH_DATE),
    AGE_NOW = as.numeric(difftime(ymd("2024-02-01"), BIRTH_DATE, units = "days")) / 365.25,
    DISEASE_DURATION = as.numeric(difftime(ymd("2024-02-01"), ymd(date_msdx), units = "days")) / 365.25
  )

# =========================
# Pitt: summary table for demographics (AGE_NOW)
# =========================

library(dplyr)
library(lubridate)
library(tableone)
library(openxlsx)

# Step 1: Parse birth date and compute AGE_NOW
pitt <- pitt %>%
  mutate(
    BIRTH_DATE = dmy(BIRTH_DATE),
    AGE_NOW = as.numeric(difftime(ymd("2024-02-01"), BIRTH_DATE, units = "days")) / 365.25
  )

# Step 2: Variables for TableOne
var <- c("AGE_NOW", "GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE")
catVar <- c("GENDER_TITLE", "RACE_TITLE", "ETHNIC_TITLE")

# Step 3: Create TableOne
summary_tbl <- CreateTableOne(vars = var, data = pitt, factorVars = catVar)
tabMat <- print(summary_tbl, printToggle = FALSE)

# Step 4: Export
write.xlsx(tabMat, file = "pitt_summary_table.xlsx", rowNames = TRUE)

# =========================
# AUROC plots: circular and side-by-side bar
# =========================

library(tidyverse)
library(ggplot2)

# Data frame based on extracted AUROC values
df <- tribble(
  ~Method,             ~Site,    ~AUROC,
  "KOMAP Codified",    "UPMC",   0.912,
  "KOMAP Cod+NLP",     "UPMC",   0.922,
  "PheNorm Codified",  "UPMC",   0.913,
  "PheNorm Cod+NLP",   "UPMC",   0.909,
  "KOMAP Codified",    "MGB",    0.974,
  "KOMAP Cod+NLP",     "MGB",    0.994,
  "PheNorm Codified",  "MGB",    0.959,
  "PheNorm Cod+NLP",   "MGB",    0.983
)

# Order for display and compute polar text angles
df <- df %>%
  arrange(AUROC) %>%
  mutate(
    Method = factor(Method, levels = Method),
    angle = 90 - 360 * (row_number() - 0.5) / n(),
    label = paste(Method, Site, sep = " - ")
  )

# Circular bar chart (polar)
ggplot(df, aes(x = Method, y = AUROC, fill = Method)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  geom_text(aes(y = AUROC + 0.03, label = Method, angle = angle),
            color = "black", size = 3.5, hjust = 0) +
  ggtitle("Mean AUC (UPMC)")

# Note: The next block in Rmd referenced objects (id, value_trans, radial_breaks, radial_labels)
# that were not defined. We keep the plotting call as in source, but it will not run without those.
# If needed, define those variables before running.

# Placeholder for the referenced plot object `p` (kept as-is from Rmd; variables must exist):
# p <- ggplot(df, aes(x = factor(id), y = value_trans, fill = Method)) + ...

# Save (expects `p` to exist)
# ggsave("auroc_circular_plot.svg", plot = p, width = 10, height = 10, dpi = 600)
# p

# Side-by-side bar plot by site
sidebyside <- ggplot(df, aes(x = Method, y = AUROC, fill = Site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = round(AUROC, 3)),
            position = position_dodge(width = 0.8),
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "AUROC Comparison by Method and Site",
       x = "Method",
       y = "AUROC") +
  theme_minimal()

ggsave("auroc_comparison.svg", plot = sidebyside, width = 8, height = 6, dpi = 600)

# Save circular plot again if `p` exists; kept to mirror original Rmd
# ggsave("auroc_circular_0.9_to_1.svg", plot = p, width = 10, height = 10, dpi = 600)
# p
