# ==============================================================================
# 04_cox_regression.R
# Multivariable Cox regression for overall survival in TCGA-BRCA.
# Model 1: Molecular subtype (PAM50) adjusted for age and stage
# Model 2: Clinical stage adjusted for age and subtype
# Output: results/cox_model1.csv, results/cox_model2.csv
#         plots/cox_forest_model1.png, plots/cox_forest_model2.png
# ==============================================================================

library(survival)
library(broom)
library(dplyr)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("plots",   showWarnings = FALSE)

brca <- readRDS("data/brca_analytic.rds")

# Use only subjects with complete covariates
brca_complete <- brca %>%
  filter(!is.na(stage_group), !is.na(subtype), !is.na(age_at_index))

cat("Subjects in Cox analysis:", nrow(brca_complete), "\n")
cat("Events:", sum(brca_complete$os_event), "\n")

# ------------------------------------------------------------------------------
# Helper: tidy Cox output with HR and 95% CI
# Reference levels: Luminal A (subtype), Stage I (stage)
# ------------------------------------------------------------------------------
tidy_cox <- function(fit) {
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, HR = estimate, ci_lower = conf.low, ci_upper = conf.high, p.value) %>%
    mutate(across(c(HR, ci_lower, ci_upper), ~ round(., 3)),
           p.value = round(p.value, 4))
}

# ------------------------------------------------------------------------------
# Model 1: Subtype as main exposure, adjusted for age + stage
# Reference: Luminal A
# ------------------------------------------------------------------------------
cox1 <- coxph(
  Surv(os_time, os_event) ~ subtype + age_at_index + stage_group,
  data = brca_complete
)

cat("\n--- Model 1: Subtype (ref = Luminal A) ---\n")
print(summary(cox1))

results_cox1 <- tidy_cox(cox1)
write.csv(results_cox1, "results/cox_model1.csv", row.names = FALSE)
cat("Saved: results/cox_model1.csv\n")

# Proportional hazards check
ph_test1 <- cox.zph(cox1)
cat("\nPH test (p < 0.05 suggests violation):\n")
print(ph_test1)

# ------------------------------------------------------------------------------
# Model 2: Stage as main exposure, adjusted for age + subtype
# Reference: Stage I
# ------------------------------------------------------------------------------
cox2 <- coxph(
  Surv(os_time, os_event) ~ stage_group + age_at_index + subtype,
  data = brca_complete
)

cat("\n--- Model 2: Stage (ref = Stage I) ---\n")
print(summary(cox2))

results_cox2 <- tidy_cox(cox2)
write.csv(results_cox2, "results/cox_model2.csv", row.names = FALSE)
cat("Saved: results/cox_model2.csv\n")

ph_test2 <- cox.zph(cox2)
cat("\nPH test:\n")
print(ph_test2)

# ------------------------------------------------------------------------------
# Forest plot helper: clean term labels
# ------------------------------------------------------------------------------
label_terms <- function(df) {
  df %>%
    mutate(term = case_when(
      term == "subtypeLuminal B"    ~ "Luminal B vs Luminal A",
      term == "subtypeHER2-enriched"~ "HER2-enriched vs Luminal A",
      term == "subtypeBasal-like"   ~ "Basal-like vs Luminal A",
      term == "subtypeNormal-like"  ~ "Normal-like vs Luminal A",
      term == "stage_groupStage II" ~ "Stage II vs Stage I",
      term == "stage_groupStage III"~ "Stage III vs Stage I",
      term == "stage_groupStage IV" ~ "Stage IV vs Stage I",
      term == "age_at_index"        ~ "Age (per year)",
      TRUE ~ term
    ))
}

plot_cox_forest <- function(results_df, title_str, filename) {
  df <- label_terms(results_df) %>%
    filter(!grepl("subtype|stage_group", term) |
             grepl("vs", term)) %>%
    mutate(sig = p.value < 0.05,
           term = factor(term, levels = rev(term)))

  p <- ggplot(df, aes(x = HR, y = term, colour = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                   height = 0.25, linewidth = 0.8) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                        labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
                        name   = NULL) +
    labs(
      title = title_str,
      x     = "Hazard Ratio (95% CI)",
      y     = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y      = element_text(size = 12),
      axis.text.x      = element_text(size = 12),
      axis.title.x     = element_text(size = 13),
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 14, face = "bold"),
      plot.margin      = margin(20, 30, 20, 20)
    )

  ggsave(filename, p, width = 7, height = 5, dpi = 300, bg = "white")
  cat("Saved:", filename, "\n")
}

plot_cox_forest(results_cox1,
                "Multivariable Cox: Molecular Subtype\n(ref = Luminal A)",
                "plots/cox_forest_subtype.png")

plot_cox_forest(results_cox2,
                "Multivariable Cox: Clinical Stage\n(ref = Stage I)",
                "plots/cox_forest_stage.png")
