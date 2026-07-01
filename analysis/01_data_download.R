# ==============================================================================
# 01_data_download.R
# Download TCGA-BRCA clinical and molecular subtype data via TCGAbiolinks
# Output: data/clinical_raw.rds
# ==============================================================================

# install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(dplyr)

dir.create("data", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Download TCGA-BRCA clinical data
# ------------------------------------------------------------------------------
clinical_raw <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# ------------------------------------------------------------------------------
# Download molecular subtype annotations (PAM50)
# Available via TCGAbiolinks::TCGAquery_subtype()
# ------------------------------------------------------------------------------
subtype_raw <- TCGAquery_subtype(tumor = "brca")

# Save raw downloads for reproducibility
saveRDS(clinical_raw, "data/clinical_raw.rds")
saveRDS(subtype_raw,  "data/subtype_raw.rds")

cat("Clinical records downloaded:", nrow(clinical_raw), "\n")
cat("Subtype records downloaded: ", nrow(subtype_raw), "\n")
