# Statistical analysis of HBB E7V editing rates across treatment groups
library(dplyr)
library(ggplot2)
library(broom)

# Create the dataset
e7v_data <- data.frame(
  Sample = c(
    "HBB E7V UV #4", "HBB E7V UV #3", "HBB E7V UV #2", "HBB E7V UV #1",
    "HBB E7V Unmodified #4", "HBB E7V Unmodified #3", "HBB E7V Unmodified #1",
    "HBB E7V EtOH #4", "HBB E7V EtOH #3", "HBB E7V EtOH #2", "HBB E7V EtOH #1",
    "HBB E7V 3X #4", "HBB E7V 3X #3", "HBB E7V 3X #2", "HBB E7V 3X #1",
    "HBB E7V 30X #4", "HBB E7V 30X #3", "HBB E7V 30X #2", "HBB E7V 30X #1",
    "HBB E7V 1X #4", "HBB E7V 1X #3", "HBB E7V 1X #2", "HBB E7V 1X #1",
    "HBB E7V 10X #4", "HBB E7V 10X #3", "HBB E7V 10X #2", "HBB E7V 10X #1",
    "HBB E7V 0.3X #4", "HBB E7V 0.3X #3", "HBB E7V 0.3X #2", "HBB E7V 0.3X #1",
    "HBB E7V 0.1X #4", "HBB E7V 0.1X #3", "HBB E7V 0.1X #2", "HBB E7V 0.1X #1",
    "HBB E7V 0.03X #3", "HBB E7V 0.03X #2", "HBB E7V 0.03X #1",
    "HBB E7V 0.01X #4", "HBB E7V 0.01X #3", "HBB E7V 0.01X #2", "HBB E7V 0.01X #1"
  ),
  Value = c(
    0.80270388, 0.21311475, 0.50772505, 0.45266078,
    0.01178011, 0.03391089, 0.00609244,
    0.7919341, 0.6160719, 0.87703765, 0.8436154,
    0.90578881, 0.80345756, 0.84785294, 0.86636796,
    0.58864452, 0.56678787, 0.22049689, 0.22458257,
    0.88888889, 0.73877289, 0.82085562, 0.81520189,
    0.70445533, 0.59811709, 0.55953843, 0.53599277,
    0.30592073, 0.62, 0.79148665, 0.76306588,
    0.88230382, 0.54952488, 0.80016269, 0.44827586,
    0.64881742, 0.85260969, 0.86035948,
    0.48459384, 0.71401114, 0.81892174, 0.8359375
  )
)
# Extract treatment groups
e7v_data$Group <- gsub(" #\\d+$", "", e7v_data$Sample)
e7v_data$Group <- gsub("HBB E7V ", "", e7v_data$Group)

# Create proper factor levels in logical order
group_order <- c("Unmodified", "UV", "EtOH", "0.01X", "0.03X", "0.1X", "0.3X", 
                 "1X", "3X", "10X", "30X")
e7v_data$Group <- factor(e7v_data$Group, levels = group_order)

# Print summary statistics by group
cat("=== Summary Statistics by Treatment Group ===\n")
summary_stats <- e7v_data %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    mean = mean(Value),
    sd = sd(Value),
    se = sd/sqrt(n),
    min = min(Value),
    max = max(Value),
    .groups = "drop"
  )
print(summary_stats)

# Test specific comparisons of interest
cat("\n=== Specific Comparisons of Interest ===\n")

# Compare UV vs all other treatment groups
uv_mean <- mean(e7v_data$Value[e7v_data$Group == "UV"])
cat(sprintf("UV group mean: %.4f\n", uv_mean))

treatment_groups <- setdiff(levels(e7v_data$Group), "UV")
for (group in treatment_groups) {
  group_data <- e7v_data$Value[e7v_data$Group == group]
  group_mean <- mean(group_data)
  
  # Perform t-test
  t_result <- t.test(e7v_data$Value[e7v_data$Group == "UV"],
                     group_data)
  
  cat(sprintf("%s vs UV: %.4f vs %.4f, p = %.6f\n",
              group, group_mean, uv_mean, t_result$p.value))
}

# Apply Bonferroni correction to these specific tests
cat("\n=== Bonferroni Corrected p-values for comparisons with UV ===\n")
n_comparisons <- length(treatment_groups)
for (group in treatment_groups) {
  group_data <- e7v_data$Value[e7v_data$Group == group]
  t_result <- t.test(e7v_data$Value[e7v_data$Group == "UV"],
                     group_data)
  
  bonferroni_p <- min(t_result$p.value * n_comparisons, 1.0)
  significance <- ifelse(bonferroni_p < 0.05, "*", "")
  
  cat(sprintf("%s vs UV: p.bonferroni = %.6f %s\n",
              group, bonferroni_p, significance))
}

# Check for dose-response relationship in crosslinking concentrations
cat("\n=== Dose-Response Analysis ===\n")
crosslink_groups <- c("0.01X", "0.03X", "0.1X", "0.3X", "1X", "3X", "10X", "30X")
crosslink_data <- e7v_data[e7v_data$Group %in% crosslink_groups, ]

# Convert to numeric doses for correlation
dose_mapping <- c("0.01X" = 0.01, "0.03X" = 0.03, "0.1X" = 0.1, "0.3X" = 0.3, 
                  "1X" = 1, "3X" = 3, "10X" = 10, "30X" = 30)
crosslink_data$Dose <- dose_mapping[as.character(crosslink_data$Group)]

# Log-transform dose for better correlation
crosslink_data$LogDose <- log10(crosslink_data$Dose)

# Correlation analysis
cor_result <- cor.test(crosslink_data$LogDose, crosslink_data$Value)
cat(sprintf("Correlation between log10(dose) and editing rate: r = %.4f, p = %.6f\n", 
            cor_result$estimate, cor_result$p.value))

# Linear regression
lm_result <- lm(Value ~ LogDose, data = crosslink_data)
lm_summary <- summary(lm_result)
cat(sprintf("Linear regression R² = %.4f, p = %.6f\n", 
            lm_summary$r.squared, lm_summary$coefficients[2,4]))

cat("\n=== Analysis Complete ===\n")
