# ==============================================================================
# 05_subgroup_forest.R
# Subgroup analysis: effect of molecular subtype (Luminal A vs Other)
# stratified across clinically relevant subgroups.
# Produces interaction-test p-values and a publication-ready forest plot.
# Output: results/subgroup_results.csv, plots/forest_subgroup.png
# ==============================================================================

library(survival)
library(broom)
library(dplyr)
library(ggplot2)
library(tidyr)

brca <- readRDS("data/brca_analytic.rds") %>%
  filter(!is.na(stage_group), !is.na(subtype), !is.na(age_at_index)) %>%
  mutate(subtype_lumA = factor(subtype_lumA, levels = c("Other", "Luminal A")))

# ------------------------------------------------------------------------------
# Define subgroups
# ------------------------------------------------------------------------------
subgroups <- list(
  list(var = "age_group",    label = "Age",
       levels = c("< 60", ">= 60")),
  list(var = "stage_binary", label = "Stage",
       levels = c("Stage I-II", "Stage III-IV")),
  list(var = "gender",       label = "Sex",
       levels = c("female", "male"))
)

# ------------------------------------------------------------------------------
# Helper: fit Cox within one stratum
# ------------------------------------------------------------------------------
fit_stratum <- function(data, subgroup_var, stratum_label) {
  na_row <- data.frame(
    subgroup = subgroup_var, stratum = stratum_label,
    HR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    n = nrow(data), events = sum(data$os_event), p.value = NA_real_
  )
  if (sum(data$os_event) < 5)                return(na_row)
  if (length(unique(data$subtype_lumA)) < 2) return(na_row)
  tryCatch({
    fit <- coxph(Surv(os_time, os_event) ~ subtype_lumA + age_at_index,
                 data = data)
    res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "subtype_lumALuminal A")
    if (nrow(res) == 0) return(na_row)
    data.frame(
      subgroup = subgroup_var, stratum = stratum_label,
      HR       = round(res$estimate, 3),
      ci_lower = round(res$conf.low, 3),
      ci_upper = round(res$conf.high, 3),
      n        = nrow(data),
      events   = sum(data$os_event),
      p.value  = round(res$p.value, 4)
    )
  }, error = function(e) na_row)
}

# ------------------------------------------------------------------------------
# Helper: interaction p-value
# ------------------------------------------------------------------------------
interaction_pval <- function(data, subgroup_var) {
  if (length(unique(data$subtype_lumA)) < 2)    return(NA_real_)
  if (length(unique(data[[subgroup_var]])) < 2) return(NA_real_)
  tryCatch({
    fit <- coxph(
      as.formula(paste0(
        "Surv(os_time, os_event) ~ subtype_lumA * ", subgroup_var, " + age_at_index"
      )),
      data = data
    )
    res <- tidy(fit) %>% filter(grepl(":", term)) %>% slice(1)
    if (nrow(res) == 0) return(NA_real_)
    round(res$p.value, 4)
  }, error = function(e) NA_real_)
}

# ------------------------------------------------------------------------------
# Run subgroup analysis
# ------------------------------------------------------------------------------
results_list     <- list()
interaction_list <- list()

for (sg in subgroups) {
  var   <- sg$var
  label <- sg$label
  lvls  <- sg$levels

  data_sg <- brca %>% filter(!is.na(.data[[var]]))

  for (lv in lvls) {
    stratum_data <- data_sg %>% filter(.data[[var]] == lv)
    results_list[[paste(var, lv)]] <- fit_stratum(stratum_data, label, lv)
  }

  p_int <- interaction_pval(data_sg, var)
  interaction_list[[label]] <- data.frame(subgroup = label, p_interact = p_int)
}

subgroup_results <- bind_rows(results_list)
interaction_df   <- bind_rows(interaction_list)

subgroup_results <- subgroup_results %>%
  left_join(interaction_df, by = "subgroup")

write.csv(subgroup_results, "results/subgroup_results.csv", row.names = FALSE)
cat("Saved: results/subgroup_results.csv\n")
print(subgroup_results)

# ------------------------------------------------------------------------------
# Overall effect
# ------------------------------------------------------------------------------
fit_overall <- coxph(
  Surv(os_time, os_event) ~ subtype_lumA + age_at_index + stage_group,
  data = brca
)
overall_res <- tidy(fit_overall, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "subtype_lumALuminal A") %>%
  transmute(
    subgroup   = "Overall", stratum = "Overall",
    HR         = round(estimate, 3),
    ci_lower   = round(conf.low, 3),
    ci_upper   = round(conf.high, 3),
    n          = nrow(brca),
    events     = sum(brca$os_event),
    p.value    = round(p.value, 4),
    p_interact = NA_real_
  )

# ------------------------------------------------------------------------------
# Build plot_data
# y_label embeds p-interaction for header rows; data rows are indented
# Row order (top to bottom): Overall, Age header, <60, >=60,
#                             Stage header, I-II, III-IV, Sex header, female, male
# ------------------------------------------------------------------------------

# Merge p_interact into subgroup_results for label building
all_results <- bind_rows(overall_res, subgroup_results)

make_header_label <- function(group_name, p_int) {
  if (is.na(p_int)) {
    group_name
  } else {
    paste0(group_name, "  (p-interaction = ", p_int, ")")
  }
}

# Header rows: one per subgroup with p-interaction in label
header_rows <- interaction_df %>%
  rowwise() %>%
  mutate(
    y_label  = make_header_label(subgroup, p_interact),
    stratum  = subgroup,
    HR       = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    n        = NA_integer_, events = NA_integer_,
    p.value  = NA_real_,
    is_header = TRUE,
    ci_text  = ""
  ) %>%
  ungroup()

# Data rows: indented label, HR text
data_rows <- all_results %>%
  mutate(
    y_label   = paste0("  ", stratum),
    is_header = FALSE,
    ci_text   = ifelse(is.na(HR), "",
                       sprintf("%.2f (%.2f\u2013%.2f)", HR, ci_lower, ci_upper))
  )

# Define y-axis order (factor levels, reversed so top = first)
row_order_labels <- c(
  "  Overall",
  make_header_label("Age",   interaction_df$p_interact[interaction_df$subgroup == "Age"]),
  "  < 60", "  >= 60",
  make_header_label("Stage", interaction_df$p_interact[interaction_df$subgroup == "Stage"]),
  "  Stage I-II", "  Stage III-IV",
  make_header_label("Sex",   interaction_df$p_interact[interaction_df$subgroup == "Sex"]),
  "  female", "  male"
)

plot_data <- bind_rows(data_rows, header_rows) %>%
  mutate(y_label = factor(y_label, levels = rev(row_order_labels)))

# ------------------------------------------------------------------------------
# Forest plot
# ------------------------------------------------------------------------------
x_lim <- c(0.2, 3.0)

forest_plot <- ggplot(plot_data,
                      aes(x = HR, y = y_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.3, linewidth = 0.8,
                 colour = "#2C7BB6", na.rm = TRUE) +
  geom_point(aes(size = events), shape = 15,
             colour = "#2C7BB6", na.rm = TRUE) +
  # HR (95% CI) text on the right
  geom_text(aes(x = x_lim[2] * 1.05, label = ci_text),
            hjust = 0, size = 3.5) +
  scale_x_continuous(
    limits = c(x_lim[1], x_lim[2] * 1.6),
    breaks = c(0.25, 0.5, 1, 2, 3),
    trans  = "log",
    labels = c("0.25", "0.5", "1.0", "2.0", "3.0")
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(
    x     = "Hazard Ratio for Luminal A vs Other (95% CI)",
    y     = NULL,
    title = "Subgroup Analysis: Overall Survival\nLuminal A vs Other Subtypes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y        = element_text(size = 11),
    axis.text.x        = element_text(size = 11),
    axis.title.x       = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(size = 13, face = "bold"),
    plot.margin        = margin(20, 80, 20, 20)
  )

ggsave("plots/forest_subgroup.png",
       forest_plot, width = 10, height = 6, dpi = 300, bg = "white")

cat("Saved: plots/forest_subgroup.png\n")
