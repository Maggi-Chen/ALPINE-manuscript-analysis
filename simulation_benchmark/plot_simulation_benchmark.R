#!/usr/bin/env Rscript
setwd('ALPINE_Manuscript/revision/')
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Read the simulation benchmark data
data <- read_csv("simulation_benchmark.csv")

# Clean up column names and remove empty columns
data <- data %>%
  select(Site, Template, Class, TP, FP, FN, Precision, Recall, F1, Tool, Data) %>%
  filter(!is.na(Site))

# Define integration and variant categories
integration_categories <- c("HDR", "ITR", "Truncated HDR", "Cross-contamination-HDR", "Cross-contamination-ITR")
variant_categories <- c("Wild type", "SNP", "large DEL", "small DEL", "large INS", "small INS")

# Create category groupings
data <- data %>%
  mutate(
    Category_Type = case_when(
      Template %in% integration_categories ~ "Integration",
      Template %in% variant_categories ~ "Variants",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Category_Type != "Other")

# Create a combined identifier for tool and platform
data <- data %>%
  mutate(
    Tool_Platform = paste(Tool, Data, sep = "-"),
    Panel = paste(Site, Category_Type, sep = " ")
  )

# Set factor levels for Tool_Platform to control ordering
data$Tool_Platform <- factor(data$Tool_Platform, levels = c(
  "ALPINE-PacBio HiFi",
  "ALPINE-ONT",
  "knock-knock-PacBio HiFi",
  "knock-knock-ONT"
))

# Create multi-line labels for better readability
data <- data %>%
  mutate(
    Template_Label = case_when(
      Template == "Cross-contamination-HDR" ~ "Cross-target\nHDR",
      Template == "Cross-contamination-ITR" ~ "Cross-target\nITR",
      Template == "Truncated HDR" ~ "Truncated\nHDR",
      Template == "Wild type" ~ "Wild-type",
      Template == "large DEL" ~ "Large\nDEL",
      Template == "small DEL" ~ "Small\nDEL",
      Template == "large INS" ~ "Large\nINS",
      Template == "small INS" ~ "Small\nINS",
      TRUE ~ as.character(Template)
    )
  )

# Reorder template categories for better visualization
data$Template_Label <- factor(data$Template_Label, levels = c(
  # Integration categories
  "HDR", "ITR", "Truncated\nHDR", "Cross-target\nHDR", "Cross-target\nITR",
  # Variant categories - large variants first
  "Large\nDEL", "Large\nINS", "Small\nDEL", "Small\nINS", "SNP", "Wild-type"
))

# Reorder panels
data$Panel <- factor(data$Panel, levels = c(
  "TRAC Integration", "TRAC Variants", "TRBC Integration", "TRBC Variants"
))

# Create color palette for tools and platforms - similar colors for each tool family
colors <- c(
  "ALPINE-PacBio HiFi" = "#cb181d",     # Dark red
  "ALPINE-ONT" = "#fb6a4a",             # Light red
  "knock-knock-PacBio HiFi" = "#08519c", # Dark blue
  "knock-knock-ONT" = "#6baed6"         # Light blue
)

# Create the plot
p <- ggplot(data, aes(x = Template_Label, y = F1, fill = Tool_Platform)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8, preserve = "total"), width = 0.7) +
  facet_wrap(~ Panel, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = colors, name = "Tool-Platform") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(
    x = "Category",
    y = "F1-Score (%)",
    fill = "Tool-Platform"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 20)),
    axis.text.x = element_text(hjust = 0.5, size = 9, lineheight = 0.8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Print the plot
print(p)

# Save the plot
ggsave("simulation_benchmark_f1_comparison.pdf", plot = p, 
       width = 10, height = 8, units = "in", dpi = 300)


# Print summary statistics
cat("\nSummary Statistics:\n")
cat("==================\n")

summary_stats <- data %>%
  group_by(Panel, Tool) %>%
  summarise(
    Mean_F1 = round(mean(F1, na.rm = TRUE), 2),
    Median_F1 = round(median(F1, na.rm = TRUE), 2),
    Min_F1 = round(min(F1, na.rm = TRUE), 2),
    Max_F1 = round(max(F1, na.rm = TRUE), 2),
    .groups = 'drop'
  )

print(summary_stats)

# Print platform comparison
cat("\nPlatform Comparison (Mean F1-scores):\n")
cat("=====================================\n")

platform_stats <- data %>%
  group_by(Tool, Data) %>%
  summarise(
    Mean_F1 = round(mean(F1, na.rm = TRUE), 2),
    N_Categories = n(),
    .groups = 'drop'
  ) %>%
  arrange(Tool, Data)

print(platform_stats)

cat("\nPlot saved as:\n")
cat("- simulation_benchmark_f1_comparison.pdf\n")
cat("- simulation_benchmark_f1_comparison.png\n")
