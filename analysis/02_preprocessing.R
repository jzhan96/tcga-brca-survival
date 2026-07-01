# ==============================================================================
# 02_preprocessing.R
# Merge TCGA-BRCA clinical data with PAM50 molecular subtypes.
# Define survival endpoints, clean covariates, and save analytic dataset.
# Output: data/brca_analytic.rds
# ==============================================================================

library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# Load raw data
# ------------------------------------------------------------------------------
clinical_raw <- readRDS("data/clinical_raw.rds")
subtype_raw  <- readRDS("data/subtype_raw.rds")

# ------------------------------------------------------------------------------
# Clean clinical data
# Retain columns needed for survival analysis
# ------------------------------------------------------------------------------
clinical <- clinical_raw %>%
  select(
    submitter_id,
    days_to_death,
    days_to_last_follow_up,
    vital_status,
    age_at_index,
    tumor_stage = ajcc_pathologic_stage,
    gender = sex_at_birth
  ) %>%
  # Define overall survival time (months) and event indicator
  mutate(
    os_time = case_when(
      vital_status == "Dead"  ~ as.numeric(days_to_death) / 30.44,
      vital_status == "Alive" ~ as.numeric(days_to_last_follow_up) / 30.44,
      TRUE ~ NA_real_
    ),
    os_event = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  # Remove records with missing or non-positive survival time
  filter(!is.na(os_time), os_time > 0) %>%
  # Clean stage: collapse to broad groups
  mutate(
    stage_group = case_when(
      grepl("Stage I",  tumor_stage, ignore.case = TRUE) &
        !grepl("Stage II|III|IV", tumor_stage, ignore.case = TRUE) ~ "Stage I",
      grepl("Stage II", tumor_stage, ignore.case = TRUE) &
        !grepl("Stage III|IV", tumor_stage, ignore.case = TRUE)    ~ "Stage II",
      grepl("Stage III", tumor_stage, ignore.case = TRUE) &
        !grepl("Stage IV", tumor_stage, ignore.case = TRUE)        ~ "Stage III",
      grepl("Stage IV", tumor_stage, ignore.case = TRUE)           ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    # Binary stage grouping for forest plot subgroup
    stage_binary = case_when(
      stage_group %in% c("Stage I", "Stage II")   ~ "Stage I-II",
      stage_group %in% c("Stage III", "Stage IV") ~ "Stage III-IV",
      TRUE ~ NA_character_
    ),
    # Age group
    age_group = ifelse(age_at_index < 60, "< 60", ">= 60"),
    age_group = factor(age_group, levels = c("< 60", ">= 60"))
  ) %>%
  rename(patient_id = submitter_id)

# ------------------------------------------------------------------------------
# Clean subtype data
# PAM50 subtypes: Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like
# ------------------------------------------------------------------------------
subtype <- subtype_raw %>%
  select(patient_id = patient, subtype = BRCA_Subtype_PAM50) %>%
  mutate(
    subtype = case_when(
      subtype == "LumA"   ~ "Luminal A",
      subtype == "LumB"   ~ "Luminal B",
      subtype == "Her2"   ~ "HER2-enriched",
      subtype == "Basal"  ~ "Basal-like",
      subtype == "Normal" ~ "Normal-like",
      TRUE ~ NA_character_
    ),
    # Binary subtype for forest plot
    subtype_lumA = ifelse(subtype == "Luminal A", "Luminal A", "Other")
  ) %>%
  filter(!is.na(subtype))

# ------------------------------------------------------------------------------
# Merge clinical and subtype data
# ------------------------------------------------------------------------------
brca <- inner_join(clinical, subtype, by = "patient_id") %>%
  # Set reference levels for Cox model
  mutate(
    subtype     = factor(subtype,
                         levels = c("Luminal A", "Luminal B",
                                    "HER2-enriched", "Basal-like",
                                    "Normal-like")),
    stage_group  = factor(stage_group,
                          levels = c("Stage I", "Stage II",
                                     "Stage III", "Stage IV")),
    stage_binary = factor(stage_binary,
                          levels = c("Stage I-II", "Stage III-IV")),
    subtype_lumA = factor(subtype_lumA,
                          levels = c("Other", "Luminal A"))
  )

cat("Final analytic sample:", nrow(brca), "subjects\n")
cat("Subtype distribution:\n")
print(table(brca$subtype))
cat("Stage distribution:\n")
print(table(brca$stage_group))
cat("Events (deaths):", sum(brca$os_event), "\n")
cat("Median follow-up (months):", round(median(brca$os_time), 1), "\n")

saveRDS(brca, "data/brca_analytic.rds")
