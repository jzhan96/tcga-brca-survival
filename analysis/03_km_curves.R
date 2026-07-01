# ==============================================================================
# 03_km_curves.R
# Kaplan-Meier survival curves stratified by PAM50 molecular subtype and stage.
# Log-rank test p-values displayed on each plot.
# Output: plots/km_subtype.png, plots/km_stage.png
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)

brca <- readRDS("data/brca_analytic.rds")

# ------------------------------------------------------------------------------
# Shared theme for publication-quality figures
# ------------------------------------------------------------------------------
km_theme <- theme_minimal(base_size = 13) +
  theme(
    axis.title.x     = element_text(size = 14),
    axis.title.y     = element_text(size = 14),
    axis.text        = element_text(size = 12),
    legend.text      = element_text(size = 12),
    legend.title     = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(20, 20, 20, 20)
  )

# ------------------------------------------------------------------------------
# Figure 1: KM curves by PAM50 molecular subtype (5 groups)
# ------------------------------------------------------------------------------
fit_subtype <- survfit(Surv(os_time, os_event) ~ subtype, data = brca)

km_subtype <- ggsurvplot(
  fit_subtype,
  data          = brca,
  pval          = TRUE,
  pval.size     = 5,
  conf.int      = FALSE,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  xlab          = "Time since diagnosis (months)",
  ylab          = "Overall Survival Probability",
  legend.title  = "",
  legend.labs   = levels(brca$subtype),
  palette       = c("#2196F3", "#4CAF50", "#FF9800", "#F44336", "#9C27B0"),
  ggtheme       = km_theme,
  fontsize      = 4
)

ggsave("plots/km_subtype.png",
       print(km_subtype),
       width = 9, height = 7.5, dpi = 300, bg = "white")

cat("Saved: plots/km_subtype.png\n")

# ------------------------------------------------------------------------------
# Figure 2: KM curves by clinical stage (I–IV)
# ------------------------------------------------------------------------------
brca_stage <- brca %>% filter(!is.na(stage_group))

fit_stage <- survfit(Surv(os_time, os_event) ~ stage_group, data = brca_stage)

km_stage <- ggsurvplot(
  fit_stage,
  data          = brca_stage,
  pval          = TRUE,
  pval.size     = 5,
  conf.int      = FALSE,
  risk.table    = TRUE,
  risk.table.height = 0.28,
  xlab          = "Time since diagnosis (months)",
  ylab          = "Overall Survival Probability",
  legend.title  = "",
  legend.labs   = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  palette       = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"),
  ggtheme       = km_theme,
  fontsize      = 4
)

ggsave("plots/km_stage.png",
       print(km_stage),
       width = 9, height = 7.5, dpi = 300, bg = "white")

cat("Saved: plots/km_stage.png\n")

# ------------------------------------------------------------------------------
# Log-rank test results (for reporting)
# ------------------------------------------------------------------------------
cat("\nLog-rank test — subtype:\n")
print(survdiff(Surv(os_time, os_event) ~ subtype, data = brca))

cat("\nLog-rank test — stage:\n")
print(survdiff(Surv(os_time, os_event) ~ stage_group, data = brca_stage))
