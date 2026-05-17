# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(stringr)
library(purrr)  # For map_dfr function

setwd('ALPINE_Manuscript/revision/')
# Read the detailed data with full group information
data <- read.delim("combined_all_full_group_update_knock.txt", stringsAsFactors = FALSE)

# Check the structure
cat("Data structure:\n")
print(head(data))
print(str(data))

# Rename columns for consistency and convert to long format
colnames(data)[colnames(data) == "DEL"] <- "Deletion"
colnames(data)[colnames(data) == "INS"] <- "Insertion"  
colnames(data)[colnames(data) == "WT"] <- "Wild-type"

# Convert to long format for plotting
data_long <- data %>%
  pivot_longer(
    cols = c("DEL.Small", "DEL.Large", "Deletion", "Insertion", "Complex.Indel", "HDR", "Wild-type", "Other"),
    names_to = "Outcome",
    values_to = "Proportion"
  ) %>%
  mutate(
    # Clean up outcome names to use consistent format
    Outcome_Clean = case_when(
      Outcome == "DEL.Small" ~ "DEL-Small",
      Outcome == "DEL.Large" ~ "DEL-Large", 
      Outcome == "Deletion" ~ "Deletion",
      Outcome == "Insertion" ~ "Insertion",
      Outcome == "Complex.Indel" ~ "Complex-Indel",
      Outcome == "HDR" ~ "Integration",
      Outcome == "Wild-type" ~ "Wild-type",
      Outcome == "Other" ~ "Other",
      TRUE ~ Outcome
    )
  ) %>%
  # Keep all rows including zero proportions to preserve all samples
  # filter(Proportion > 0) %>%
  mutate(
    # Ensure proper factor order for stacking
    Outcome_Clean = factor(Outcome_Clean, levels = c(
      "Other", "Wild-type", "Integration", "Complex-Indel", "Insertion", "Deletion", "DEL-Large", "DEL-Small"
    )),
    # Extract treatment and time information from detailed groups
    Treatment = case_when(
      grepl("Unmodified", Group) ~ "Unmodified",
      grepl("Untreated", Group) ~ "Untreated",
      grepl("200M", Group) ~ "200M",  # Check 200M first before 10M
      grepl("10M", Group) ~ "10M",
      grepl("4M", Group) ~ "4M",
      grepl("UV", Group) ~ "UV",
      grepl("E7V", Group) ~ "CRISPResso2_Groups",
      TRUE ~ "Other"
    ),
    # Extract time point from groups with "hours"
    Time_Hours = case_when(
      grepl("\\d+ hours", Group) ~ str_extract(Group, "(\\d+) hours", group = 1),
      TRUE ~ NA_character_
    ),
    # Create detailed group name combining treatment and time
    Group_Simple = case_when(
      is.na(Time_Hours) ~ Treatment,  # E7V groups don't have time
      TRUE ~ paste0(Treatment, "_", Time_Hours, "h")  # Add time info
    ),
    # Create modification status for faceting
    Modification_Status = case_when(
      grepl("Unmodified", Group) ~ "Unmodified",
      TRUE ~ "Modified"
    ),
    Tool = factor(Tool, levels = c("ALPINE", "knock-knock", "crispresso2"))
  )

# Define colors for each outcome type (same as original)
outcome_colors <- c(
  "DEL-Small" = "#FF6347",      # Tomato (red-orange family) 
  "DEL-Large" = "#DC143C",      # Crimson red
  "Deletion" = "#FF4500",       # Orange red 
  "Insertion" = "#FFB6C1",      # Light pink 
  "Complex-Indel" = "#8B0000",          # Fire brick red 
  "Integration" = "#9370DB",            # Medium purple
  "Wild-type" = "#4169E1",      # Royal blue
  "Other" = "#505050"           # Dark grey
)

# Create summary data for tool comparison
tool_summary <- data_long %>%
  group_by(Tool, Outcome_Clean) %>%
  summarise(
    Mean_Proportion = mean(Proportion, na.rm = TRUE),
    SD_Proportion = sd(Proportion, na.rm = TRUE),
    N_Samples = n(),
    .groups = "drop"
  )

# Calculate sample counts by modification status for facet labels
sample_counts_by_mod <- data_long %>%
  group_by(Modification_Status) %>%
  summarise(N_Samples = n_distinct(Sample), .groups = "drop") %>%
  mutate(Mod_Status_Label = paste0(Modification_Status, " (n=", N_Samples, ")"))

# Create summary data for tool comparison by modification status
tool_summary_by_mod <- data_long %>%
  group_by(Tool, Modification_Status, Outcome_Clean) %>%
  summarise(
    Mean_Proportion = mean(Proportion, na.rm = TRUE),
    SD_Proportion = sd(Proportion, na.rm = TRUE),
    N_Samples = n(),
    .groups = "drop"
  ) %>%
  # Add the sample count labels
  left_join(sample_counts_by_mod, by = "Modification_Status") %>%
  mutate(Modification_Status = factor(Mod_Status_Label))

# Plot 1: Overall tool comparison (mean proportions)
p1 <- ggplot(tool_summary, aes(x = Tool, y = Mean_Proportion, fill = Outcome_Clean)) +
  geom_col(position = "stack", color = "white", size = 0.2) +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  labs(
    title = "Comparison of CRISPR Outcome Classification Tools",
    subtitle = "Mean proportions across all samples",
    x = "Tool",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
p1

# Plot 1.5: Tool comparison separated by modification status
p1.5 <- ggplot(tool_summary_by_mod, aes(x = Tool, y = Mean_Proportion, fill = Outcome_Clean)) +
  geom_col(position = "stack", color = "white", size = 0.2) +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  facet_wrap(~ Modification_Status, scales = "free_y") +
  labs(
    title = "CRISPR Tool Comparison: Modified vs Unmodified Samples",
    subtitle = "Mean proportions by modification status",
    x = "Tool",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )
p1.5

# Plot 2: By-sample comparison (showing individual samples)
# Filter to include all time variants of the main groups and E7V groups
main_treatments <- c("Unmodified", "Untreated", "UV", "4M", "10M", "200M")
e7v_treatments <- c("E7V_Unmodified", "E7V_EtOH", "E7V_UV", "E7V_0.01X", "E7V_0.03X", "E7V_0.1X", "E7V_0.3X", "E7V_1X", "E7V_3X", "E7V_10X", "E7V_30X")

data_filtered <- data_long %>%
  filter(Treatment %in% c(main_treatments, e7v_treatments)) %>%
  mutate(
    Tool_Type = case_when(
      Tool %in% c("ALPINE", "knock-knock") ~ "Long-read tools (PacBio)",
      Tool == "crispresso2" ~ "Short-read tool (Illumina)",
      TRUE ~ Tool
    ),
    # Use Group_Simple as is (already includes time info like "10M_48h")
    Group_Clean = Group_Simple,
    # Create ordering - group by treatment first, then by time
    Treatment_Order = case_when(
      Treatment == "Unmodified" ~ 1,
      Treatment == "E7V_Unmodified" ~ 2,
      Treatment == "Untreated" ~ 3,
      Treatment == "UV" ~ 4,
      Treatment == "E7V_UV" ~ 5,
      Treatment == "E7V_EtOH" ~ 6,
      Treatment == "E7V_0.01X" ~ 7,
      Treatment == "E7V_0.03X" ~ 8,
      Treatment == "E7V_0.1X" ~ 9,
      Treatment == "E7V_0.3X" ~ 10,
      Treatment == "E7V_1X" ~ 11,
      Treatment == "E7V_3X" ~ 12,
      Treatment == "E7V_10X" ~ 13,
      Treatment == "E7V_30X" ~ 14,
      Treatment == "4M" ~ 15,
      Treatment == "10M" ~ 16,
      Treatment == "200M" ~ 17,
      TRUE ~ 18
    ),
    Time_Order = case_when(
      is.na(Time_Hours) ~ 0,  # E7V groups come first
      Time_Hours == "48" ~ 1,
      Time_Hours == "96" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  # Order by treatment first, then time
  arrange(Treatment_Order, Time_Order, Sample) %>%
  mutate(
    Group_Clean = factor(Group_Clean, levels = unique(Group_Clean))
  )

# Update the ordering logic to create individual sample labels
data_filtered <- data_filtered %>%
  mutate(
    # Extract actual replicate numbers from original Group column
    Replicate_Num = case_when(
      # E7V format: "HBB E7V 0.01X #1" -> extract number after #
      str_detect(Group, "#\\d+") ~ str_extract(Group, "#(\\d+)", group = 1),
      # Regular format: "HBB 10M 2 48 hours" -> extract number before " \\d+ hours"
      str_detect(Group, "\\d+ hours") ~ str_extract(Group, "\\s(\\d+)\\s+\\d+ hours", group = 1),
      TRUE ~ "1"  # default
    ),
    # Create clean treatment label without time suffix for E7V groups
    Treatment_Clean = case_when(
      str_detect(Treatment, "E7V") ~ Treatment,  # E7V groups don't have time
      TRUE ~ str_replace(Group_Simple, "_\\d+h$", "")  # Remove time suffix for others
    ),
    # Create clean sample labels: Treatment_Rep#
    Sample_Label = paste0(Treatment_Clean, "_", Replicate_Num)
  ) %>%
# Order samples by treatment, time, then replicate
arrange(Treatment_Order, Time_Order, as.numeric(Replicate_Num)) %>%
mutate(Sample_Label = factor(Sample_Label, levels = unique(Sample_Label)))

p2 <- ggplot(data_filtered, aes(x = Sample_Label, y = Proportion, fill = Outcome_Clean)) +
  geom_col(position = "stack", color = "white", size = 0.1) +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  facet_wrap(~ Tool, scales = "free_x", ncol = 1) +
  labs(
    title = "Sample-level CRISPR Outcome Classifications",
    subtitle = "Individual samples grouped by experimental condition",
    x = "Experimental Group", 
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.ticks.x = element_line(size = 0.2),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )
p2

# Plot 3: Group-level comparison for ALPINE vs knock-knock only
pacbio_data <- data_long %>%
  filter(Tool %in% c("ALPINE", "knock-knock"), Treatment %in% main_treatments) %>%
  mutate(
    # Apply same group ordering as p2 with time information
    Treatment_Order = case_when(
      Treatment == "Unmodified" ~ 1,
      Treatment == "Untreated" ~ 2,
      Treatment == "UV" ~ 3,
      Treatment == "4M" ~ 4,
      Treatment == "10M" ~ 5,
      Treatment == "200M" ~ 6,
      TRUE ~ 7
    ),
    Time_Order = case_when(
      is.na(Time_Hours) ~ 0,
      Time_Hours == "48" ~ 1,
      Time_Hours == "96" ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(Treatment_Order, Time_Order) %>%
  mutate(
    Group_Simple = factor(Group_Simple, levels = unique(Group_Simple))
  ) %>%
  group_by(Tool, Group_Simple, Outcome_Clean) %>%
  summarise(
    Mean_Proportion = mean(Proportion, na.rm = TRUE),
    .groups = "drop"
  )

p3 <- ggplot(pacbio_data, aes(x = Group_Simple, y = Mean_Proportion, fill = Outcome_Clean)) +
  geom_col(position = "stack", color = "white", size = 0.2) +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  facet_wrap(~ Tool, ncol = 2) +
  labs(
    title = "ALPINE vs knock-knock: PacBio Long-read Comparison",
    subtitle = "Mean proportions by experimental group",
    x = "Experimental Group",
    y = "Proportion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )

# Plot 4: HDR-only comparison per sample
hdr_data <- data_filtered %>%
  filter(Outcome_Clean == "HDR") %>%
  # Ensure we have the same ordering as p2
  arrange(Group_Order, Sample)

p4 <- ggplot(hdr_data, aes(x = Sample_Label, y = Proportion, fill = Tool)) +
  geom_col(position = "dodge", color = "white", size = 0.1) +
  scale_fill_manual(values = c("ALPINE" = "#2E86C1", "knock-knock" = "#E74C3C", "crispresso2" = "#28B463"), 
                   name = "Tool") +
  facet_wrap(~ Tool, scales = "free_x", ncol = 1) +
  labs(
    title = "HDR Detection Rates by Tool and Sample",
    subtitle = "Homology-directed repair proportions across individual samples",
    x = "Experimental Group", 
    y = "HDR Proportion"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.ticks.x = element_line(size = 0.2),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )

# Plot 5: Simple per-sample plot using Group names directly
# Order groups logically
data_p5 <- data_long %>%
  mutate(
    # Remove "HBB " prefix from Group names for cleaner display
    Group_Clean = str_replace(Group, "^HBB ", ""),
    # Create ordering based on raw Group names
    Group_Order = case_when(
      # Proper control hierarchy: least to most intervention
      # 1. Unmodified groups first (no CRISPR)
      str_detect(Group, "Unmodified") ~ 1,
      str_detect(Group, "E7V.*Unmodified") ~ 2,
      # 2. Untreated groups (CRISPR only)
      str_detect(Group, "Untreated") ~ 3,
      # 3. UV groups (CRISPR + UV)
      str_detect(Group, "UV") ~ 4,
      str_detect(Group, "E7V.*UV") ~ 5,
      # 4. EtOH groups (CRISPR + solvent)
      str_detect(Group, "E7V.*EtOH") ~ 6,
      # 5. E7V dose groups in ascending order
      str_detect(Group, "E7V.*0\\.01X") ~ 7,
      str_detect(Group, "E7V.*0\\.03X") ~ 8,
      str_detect(Group, "E7V.*0\\.1X") ~ 9,
      str_detect(Group, "E7V.*0\\.3X") ~ 10,
      str_detect(Group, "E7V.*1X") ~ 11,
      str_detect(Group, "E7V.*3X") ~ 12,
      str_detect(Group, "E7V.*10X") ~ 13,
      str_detect(Group, "E7V.*30X") ~ 14,
      # 6. PacBio dose groups (absolute concentrations)
      str_detect(Group, "4M") ~ 15,
      str_detect(Group, "10M") ~ 16,
      str_detect(Group, "200M") ~ 17,
      TRUE ~ 18
    ),
    # Extract time and replicate for secondary sorting
    Time_Hours = case_when(
      str_detect(Group, "(\\d+) hours") ~ str_extract(Group, "(\\d+) hours", group = 1),
      TRUE ~ "0"
    ),
    Replicate = case_when(
      str_detect(Group, "#\\d+") ~ str_extract(Group, "#(\\d+)", group = 1),
      str_detect(Group, "\\s(\\d+)\\s+\\d+ hours") ~ str_extract(Group, "\\s(\\d+)\\s+\\d+ hours", group = 1),
      TRUE ~ "1"
    )
  ) %>%
  # Order by group, then time, then replicate
  arrange(Group_Order, as.numeric(Time_Hours), as.numeric(Replicate)) %>%
  mutate(Group_Clean = factor(Group_Clean, levels = unique(Group_Clean)))

p5 <- ggplot(data_p5, aes(x = Group_Clean, y = Proportion, fill = Outcome_Clean)) +
  geom_col(position = "stack", color = "white", size = 0.1) +
  scale_fill_manual(values = outcome_colors, name = "Outcome") +
  facet_wrap(~ Tool, scales = "free_x", ncol = 1) +
  labs(
    x = "Sample",
    y = "Proportion of outcome class"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.ticks.x = element_line(size = 0.2),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold")
  )
p5

# Plot 5.5: HDR-only comparison using same sample order as p5 (ALPINE and knock-knock only)
data_p55 <- data_p5 %>%
  filter(Outcome_Clean == "Integration", Tool %in% c("ALPINE", "knock-knock"))

# Colorblind-friendly palette with blue gradient for doses
simple_group_colors <- c(
  "Unmodified" = "#E31A1C",        # Bright red (control)
  "Untreated" = "#FF7F00",         # Orange (control)
  "UV" = "#9370DB",                # Light purple (physical treatment)
  "4M" = "#87CEEB",                # Light blue (low dose - very visible for tiny bars)
  "10M" = "#4682B4",               # Medium blue (medium dose)
  "200M" = "#1E3A8A"               # Dark blue (high dose)
)

# Statistical testing: Compare each treatment vs Untreated for HDR rates
# Separate tests for ALPINE and knock-knock
stat_results <- data_p55 %>%
  filter(Treatment != "Unmodified") %>%  # Exclude unmodified from comparisons
  group_by(Tool) %>%
  do({
    tool_data <- .
    untreated_data <- tool_data %>% filter(Treatment == "Untreated")
    
    # Compare each treatment group vs Untreated
    treatments <- c("UV", "4M", "10M", "200M")
    
    map_dfr(treatments, function(treat) {
      treat_data <- tool_data %>% filter(Treatment == treat)
      
      if(nrow(treat_data) > 0 && nrow(untreated_data) > 0) {
        # Use Wilcoxon test (non-parametric, good for proportions)
        test_result <- wilcox.test(treat_data$Proportion, untreated_data$Proportion,
                                   alternative = "two.sided")
        
        data.frame(
          Treatment = treat,
          p_value = test_result$p.value,
          significance = case_when(
            test_result$p.value < 0.01 ~ "***",
            test_result$p.value < 0.05 ~ "**",
            test_result$p.value < 0.1 ~ "*",
            TRUE ~ "ns"
          )
        )
      } else {
        data.frame(Treatment = treat, p_value = NA, significance = "")
      }
    })
  }) %>%
  ungroup()

# Statistical testing: Compare each treatment vs UV for HDR rates
stat_results_vs_uv <- data_p55 %>%
  filter(Treatment != "Unmodified") %>%  # Exclude unmodified from comparisons
  group_by(Tool) %>%
  do({
    tool_data <- .
    uv_data <- tool_data %>% filter(Treatment == "UV")
    
    # Compare each treatment group vs UV
    treatments <- c("Untreated", "4M", "10M", "200M")
    
    map_dfr(treatments, function(treat) {
      treat_data <- tool_data %>% filter(Treatment == treat)
      
      if(nrow(treat_data) > 0 && nrow(uv_data) > 0) {
        # Use Wilcoxon test (non-parametric, good for proportions)
        test_result <- wilcox.test(treat_data$Proportion, uv_data$Proportion,
                                   alternative = "two.sided")
        
        data.frame(
          Treatment = treat,
          p_value = test_result$p.value,
          significance = case_when(
            test_result$p.value < 0.01 ~ "***",
            test_result$p.value < 0.05 ~ "**",
            test_result$p.value < 0.1 ~ "*",
            TRUE ~ "ns"
          )
        )
      } else {
        data.frame(Treatment = treat, p_value = NA, significance = "")
      }
    })
  }) %>%
  ungroup()

# Print statistical results
cat("\n=== HDR RATE STATISTICAL COMPARISONS (vs Untreated) ===\n")
print(stat_results)

#cat("\n=== HDR RATE STATISTICAL COMPARISONS (vs UV) ===\n")
#print(stat_results_vs_uv)

# Reorder Treatment factor for logical legend ordering
data_p55_ordered <- data_p55 %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M"))
  )

p5.5 <- ggplot(data_p55_ordered, aes(x = Group_Clean, y = Proportion, fill = Treatment)) +
  geom_col(position = "identity", color = "white", size = 0.1) +
  scale_fill_manual(values = simple_group_colors, name = "Treatment",
                    breaks = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M")) +
  facet_wrap(~ Tool, scales = "free_x", ncol = 1) +
  labs(
    x = "Group",
    y = "Integration Proportion"
  ) +
  ylim(NA, 0.02) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
    axis.title.y = element_text(margin = margin(r = 20)), # Add margin to y-axis title
    axis.ticks.x = element_line(size = 0.2),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(l = 10, r = 10, t = 10, b = 10) # Add overall plot margins
  )
p5.5

# Plot 6: Indel correlation between ALPINE and knock-knock (all samples)
# Filter data to only ALPINE and knock-knock
data_indel <- data[data$Tool %in% c("ALPINE", "knock-knock"),]

# Add total indel column
data_indel$Total_Indel <- data_indel$DEL.Small + data_indel$DEL.Large +
                          data_indel$Deletion + data_indel$Insertion +
                          data_indel$Complex.Indel

# Reshape to wide format for correlation
library(reshape2)
indel_wide <- dcast(data_indel, Sample + Group ~ Tool, value.var = "Total_Indel")

# Extract treatment information for coloring
indel_wide$Treatment <- case_when(
  grepl("Unmodified", indel_wide$Group) ~ "Unmodified",
  grepl("Untreated", indel_wide$Group) ~ "Untreated",
  grepl("200M", indel_wide$Group) ~ "200M",  # Check 200M first before 10M
  grepl("10M", indel_wide$Group) ~ "10M",
  grepl("4M", indel_wide$Group) ~ "4M",
  grepl("UV", indel_wide$Group) ~ "UV",
  TRUE ~ "Other"
)

# Order treatment factor same as p5.5
indel_wide$Treatment <- factor(indel_wide$Treatment, levels = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M"))

# Use same color scheme as p5.5
simple_group_colors <- c(
  "Unmodified" = "#E31A1C",        # Bright red (control)
  "Untreated" = "#FF7F00",         # Orange (control)
  "UV" = "#9370DB",                # Light purple (physical treatment)
  "4M" = "#87CEEB",                # Light blue (low dose - very visible for tiny bars)
  "10M" = "#4682B4",               # Medium blue (medium dose)
  "200M" = "#1E3A8A"               # Dark blue (high dose)
)

# Calculate correlation
cor_result <- cor(indel_wide$ALPINE, indel_wide$`knock-knock`, use = "complete.obs")

# Plot 6: Correlation plot (all samples)
p6 <- ggplot(indel_wide, aes(x = `knock-knock` * 100, y = ALPINE * 100, color = Treatment)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue",size=0.5) +
  scale_color_manual(values = simple_group_colors, name = "Treatment",
                     breaks = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M")) +
  labs(
    title = "All samples",
    subtitle = paste("Pearson r =", round(cor_result, 3)),
    x = "knock-knock Indel Proportion (%)",
    y = "ALPINE Indel Proportion (%)"
  ) +
  theme_minimal() + xlim(0,100) + ylim(0,100)

print(p6)

# Plot 6.5: Indel correlation excluding unmodified samples
# Filter out unmodified samples
data_indel_modified <- data_indel[!grepl("Unmodified", data_indel$Group),]
indel_wide_modified <- dcast(data_indel_modified, Sample + Group ~ Tool, value.var = "Total_Indel")

# Extract treatment information for coloring (modified samples only)
indel_wide_modified$Treatment <- case_when(
  grepl("Untreated", indel_wide_modified$Group) ~ "Untreated",
  grepl("200M", indel_wide_modified$Group) ~ "200M",  # Check 200M first before 10M
  grepl("10M", indel_wide_modified$Group) ~ "10M",
  grepl("4M", indel_wide_modified$Group) ~ "4M",
  grepl("UV", indel_wide_modified$Group) ~ "UV",
  TRUE ~ "Other"
)

# Order treatment factor same as p5.5 (excluding Unmodified)
indel_wide_modified$Treatment <- factor(indel_wide_modified$Treatment, levels = c("Untreated", "UV", "4M", "10M", "200M"))

# Calculate correlation for modified samples only
cor_result_modified <- cor(indel_wide_modified$ALPINE, indel_wide_modified$`knock-knock`, use = "complete.obs")

# Plot 6.5: Correlation plot (modified samples only)
p6.5 <- ggplot(indel_wide_modified, aes(x = `knock-knock` * 100, y = ALPINE * 100, color = Treatment)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue",size=0.5) +
  scale_color_manual(values = simple_group_colors, name = "Treatment",
                     breaks = c("Untreated", "UV", "4M", "10M", "200M")) +
  labs(
    title = "Modified Samples",
    subtitle = paste("Pearson r =", round(cor_result_modified, 3)),
    x = "knock-knock Indel Proportion (%)",
    y = "ALPINE Indel Proportion (%)"
  ) +
  theme_minimal()

print(p6.5)

# Plot 7: Individual outcome correlations (DEL.Small, DEL.Large, Insertion, Wild-type, HDR)
outcome_types <- c("DEL.Small", "DEL.Large", "Insertion", "Wild-type", "HDR")
outcome_labels <- c("Small Deletion", "Large Deletion", "Insertion", "Wild-type", "Integration")

# Create correlation plots for each outcome type (both all samples and modified samples)
correlation_plots_all <- list()
correlation_plots_modified <- list()

for(i in 1:length(outcome_types)) {
  outcome <- outcome_types[i]
  label <- outcome_labels[i]
  
  # === ALL SAMPLES PLOT ===
  # Reshape data for this outcome
  outcome_wide <- dcast(data_indel, Sample + Group ~ Tool, value.var = outcome)
  
  # Extract treatment information for coloring
  outcome_wide$Treatment <- case_when(
    grepl("Unmodified", outcome_wide$Group) ~ "Unmodified",
    grepl("Untreated", outcome_wide$Group) ~ "Untreated",
    grepl("200M", outcome_wide$Group) ~ "200M",
    grepl("10M", outcome_wide$Group) ~ "10M",
    grepl("4M", outcome_wide$Group) ~ "4M",
    grepl("UV", outcome_wide$Group) ~ "UV",
    TRUE ~ "Other"
  )
  
  # Order treatment factor
  outcome_wide$Treatment <- factor(outcome_wide$Treatment, levels = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M"))
  
  # Calculate correlation
  cor_result_outcome_all <- cor(outcome_wide$ALPINE, outcome_wide$`knock-knock`, use = "complete.obs")
  
  # Create plot for all samples
  p_outcome_all <- ggplot(outcome_wide, aes(x = `knock-knock` * 100, y = ALPINE * 100, color = Treatment)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", size = 0.5) +
    scale_color_manual(values = simple_group_colors, name = "Treatment",
                       breaks = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M")) +
    labs(
      title = paste(label, "- All Samples"),
      x = "knock-knock Proportion (%)",
      y = "ALPINE Proportion (%)"
    ) +
    annotate("text",
             x = Inf, y = -Inf,
             label = paste("r =", round(cor_result_outcome_all, 3)),
             hjust = 1.1, vjust = -0.5,
             size = 4, color = "black", fontface = "bold") +
    theme_minimal()
  
  # Store plot
  correlation_plots_all[[outcome]] <- p_outcome_all
  
  # === MODIFIED SAMPLES PLOT ===
  # Filter out unmodified samples
  data_indel_mod <- data_indel[!grepl("Unmodified", data_indel$Group),]
  outcome_wide_mod <- dcast(data_indel_mod, Sample + Group ~ Tool, value.var = outcome)
  
  # Extract treatment information for coloring (modified samples only)
  outcome_wide_mod$Treatment <- case_when(
    grepl("Untreated", outcome_wide_mod$Group) ~ "Untreated",
    grepl("200M", outcome_wide_mod$Group) ~ "200M",
    grepl("10M", outcome_wide_mod$Group) ~ "10M",
    grepl("4M", outcome_wide_mod$Group) ~ "4M",
    grepl("UV", outcome_wide_mod$Group) ~ "UV",
    TRUE ~ "Other"
  )
  
  # Order treatment factor (excluding Unmodified)
  outcome_wide_mod$Treatment <- factor(outcome_wide_mod$Treatment, levels = c("Untreated", "UV", "4M", "10M", "200M"))
  
  # Calculate correlation for modified samples
  cor_result_outcome_mod <- cor(outcome_wide_mod$ALPINE, outcome_wide_mod$`knock-knock`, use = "complete.obs")
  
  # Create plot for modified samples
  p_outcome_mod <- ggplot(outcome_wide_mod, aes(x = `knock-knock` * 100, y = ALPINE * 100, color = Treatment)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue", size = 0.5) +
    scale_color_manual(values = simple_group_colors, name = "Treatment",
                       breaks = c("Untreated", "UV", "4M", "10M", "200M")) +
    labs(
      title = paste(label, "- Modified Samples"),
      x = "knock-knock Proportion (%)",
      y = "ALPINE Proportion (%)"
    ) +
    annotate("text",
             x = Inf, y = -Inf,
             label = paste("r =", round(cor_result_outcome_mod, 3)),
             hjust = 1.1, vjust = -0.5,
             size = 4, color = "black", fontface = "bold") +
    theme_minimal()
  
  # Store plot
  correlation_plots_modified[[outcome]] <- p_outcome_mod
  
  # Print correlation info
  cat(paste("\n", label, "- All samples correlation:", round(cor_result_outcome_all, 3)))
  cat(paste("\n", label, "- Modified samples correlation:", round(cor_result_outcome_mod, 3)))
}

# Print all outcome correlation plots
cat("\n\n=== INDIVIDUAL OUTCOME CORRELATIONS ===\n")
for(i in 1:length(outcome_types)) {
  outcome <- outcome_types[i]
  print(correlation_plots_all[[outcome]])
  print(correlation_plots_modified[[outcome]])
}

# Save plots as PDFs
ggsave("tool_comparison_overall_detailed.pdf", p1, width = 10, height = 6)
ggsave("tool_comparison_by_modification_detailed.pdf", p1.5, width = 12, height = 6)
ggsave("tool_comparison_samples_detailed.pdf", p2, width = 16, height = 10)
ggsave("alpine_vs_knockknock_detailed.pdf", p3, width = 12, height = 6)
ggsave("hdr_comparison_samples_detailed.pdf", p4, width = 16, height = 10)
ggsave("raw_groups_comparison_detailed.pdf",
       p5, width = 10, height = 9)
ggsave("hdr_alpine_knockknock_comparison_detailed.pdf",
       p5.5, width = 10, height = 11)

# Save individual outcome correlation plots (both all samples and modified samples)
for(i in 1:length(outcome_types)) {
  outcome <- outcome_types[i]
  label <- outcome_labels[i]
  
  # Save all samples plot
  filename_all <- paste0("outcome_correlation_",
                        gsub("[^A-Za-z0-9]", "_", tolower(label)), "_all_samples.pdf")
  ggsave(filename_all, correlation_plots_all[[outcome]], width = 4.2, height = 2.8)
  
  # Save modified samples plot
  filename_mod <- paste0("outcome_correlation_",
                        gsub("[^A-Za-z0-9]", "_", tolower(label)), "_modified_samples.pdf")
  ggsave(filename_mod, correlation_plots_modified[[outcome]], width = 4.2, height = 2.8)
}

# Print summary statistics
cat("=== DETAILED DATA SUMMARY ===\n")
print(data_long %>% 
  group_by(Tool, Group_Simple) %>%
  summarise(N_Samples = n_distinct(Sample), .groups = "drop"))

# Display plots
print(p1)
print(p1.5)
print(p2)
print(p3)
print(p4)
print(p5)
print(p5.5)

cat("\nPlots saved as:\n")
cat("- tool_comparison_overall_detailed.pdf\n")
cat("- tool_comparison_by_modification_detailed.pdf\n")
cat("- tool_comparison_samples_detailed.pdf\n")
cat("- alpine_vs_knockknock_detailed.pdf\n")
cat("- hdr_comparison_samples_detailed.pdf\n")
cat("- raw_groups_comparison_detailed.pdf\n")
cat("- hdr_alpine_knockknock_comparison_detailed.pdf\n")

# Plot 8: 9bp deletion comparison between ALPINE and knock-knock
# Read the 9bp deletion data
del_9bp_data <- read.delim("Merged_count_9bp_del.txt",
                           header = FALSE,
                           col.names = c("Tool", "Sample", "Treatment", "Count_9bp", "Total_Reads", "Percentage_9bp"))

# Clean up the data
del_9bp_data <- del_9bp_data %>%
  mutate(
    # Extract replicate number from Sample ID
    Replicate = str_extract(Sample, "\\d+$"),
    
    # Create treatment ordering similar to previous plots
    Treatment_Order = case_when(
      Treatment == "Unmodified" ~ 1,
      Treatment == "Untreated" ~ 2,
      Treatment == "UV" ~ 3,
      Treatment == "4M" ~ 4,
      Treatment == "10M" ~ 5,
      Treatment == "200M" ~ 6,
      TRUE ~ 7
    ),
    
    # Create sample labels combining treatment and replicate
    Sample_Label = paste0(Treatment, "_", Replicate),
    
    # Factor the treatment for proper ordering
    Treatment = factor(Treatment, levels = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M"))
  ) %>%
  # Order samples by treatment and replicate
  arrange(Treatment_Order, as.numeric(Replicate)) %>%
  mutate(Sample_Label = factor(Sample_Label, levels = unique(Sample_Label)))

# Use same colors as previous plots
del_9bp_colors <- c(
  "Unmodified" = "#E31A1C",        # Bright red (control)
  "Untreated" = "#FF7F00",         # Orange (control)
  "UV" = "#9370DB",                # Light purple (physical treatment)
  "4M" = "#87CEEB",                # Light blue (low dose)
  "10M" = "#4682B4",               # Medium blue (medium dose)
  "200M" = "#1E3A8A"               # Dark blue (high dose)
)

# Calculate mean and standard error for each treatment group
del_9bp_summary_plot <- del_9bp_data %>%
  group_by(Tool, Treatment) %>%
  summarise(
    Mean_9bp_Percentage = mean(Percentage_9bp, na.rm = TRUE),
    SE_9bp_Percentage = sd(Percentage_9bp, na.rm = TRUE) / sqrt(n()),
    N_Samples = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Unmodified", "Untreated", "UV", "4M", "10M", "200M"))
  )

# Create a reference-style panel (without actual data) to add to our comparison
p8_reference <- ggplot() +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = 0, ymax = 0.34), fill = "gray80", color = "white") +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = 0, ymax = 0.33), fill = "gray80", color = "white") +
  geom_rect(aes(xmin = 2.5, xmax = 3.5, ymin = 0, ymax = 0.29), fill = "gray80", color = "white") +
  geom_rect(aes(xmin = 3.5, xmax = 4.5, ymin = 0, ymax = 0.28), fill = "gray80", color = "white") +
  geom_rect(aes(xmin = 4.5, xmax = 5.5, ymin = 0, ymax = 0.25), fill = "gray80", color = "white") +
  geom_rect(aes(xmin = 5.5, xmax = 6.5, ymin = 0, ymax = 0.25), fill = "gray80", color = "white") +
  # Add reference dots
  geom_point(aes(x = 1, y = c(0.002, 0.001)), color = "black", size = 2) +
  geom_point(aes(x = 2, y = c(0.35, 0.34, 0.33)), color = "black", size = 2) +
  geom_point(aes(x = 3, y = c(0.32, 0.31, 0.28)), color = "black", size = 2) +
  geom_point(aes(x = 4, y = c(0.32, 0.28, 0.26)), color = "black", size = 2) +
  geom_point(aes(x = 5, y = c(0.27, 0.26)), color = "black", size = 2) +
  geom_point(aes(x = 6, y = c(0.30, 0.25, 0.20)), color = "black", size = 2) +
  scale_x_continuous(breaks = 1:6, labels = c("Unedit.", "0μM", "UV", "4μM", "10μM", "200μM")) +
  scale_y_continuous(limits = c(0, 0.4)) +
  labs(
    x = "Sample",
    y = "frequency",
    title = "Indel size: -9bp"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Create the actual data plot with mean values and error bars
p8 <- ggplot() +
  # Add bars for means
  geom_col(data = del_9bp_summary_plot,
           aes(x = Treatment, y = Mean_9bp_Percentage, fill = Treatment),
           position = "dodge", color = "white", size = 0.3) +
  # Add error bars
  geom_errorbar(data = del_9bp_summary_plot,
                aes(x = Treatment,
                    ymin = Mean_9bp_Percentage - SE_9bp_Percentage,
                    ymax = Mean_9bp_Percentage + SE_9bp_Percentage),
                position = position_dodge(0.9), width = 0.25, size = 0.5) +
  # Add individual data points
  geom_point(data = del_9bp_data,
             aes(x = Treatment, y = Percentage_9bp),
             color = "black", size = 2, alpha = 0.7) +
  scale_fill_manual(values = del_9bp_colors, name = "Treatment") +
  facet_wrap(~ Tool, scales = "free_x", ncol = 2) +
  labs(
    x = "Treatment",
    y = "9bp Deletion Percentage (%)",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.ticks.x = element_line(size = 0.2),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(l = 10, r = 10, t = 10, b = 10)
  )

#print(p8_reference)
print(p8)

# Save plot 8 (both reference and actual data)
ggsave("9bp_deletion_reference_panel.pdf", p8_reference, width = 8, height = 6)
ggsave("9bp_deletion_comparison.pdf", p8, width = 12, height = 8)

# Create summary statistics for 9bp deletions
del_9bp_summary <- del_9bp_data %>%
  group_by(Tool, Treatment) %>%
  summarise(
    Mean_9bp_Percentage = mean(Percentage_9bp, na.rm = TRUE),
    SD_9bp_Percentage = sd(Percentage_9bp, na.rm = TRUE),
    N_Samples = n(),
    .groups = "drop"
  )

cat("\n=== 9BP DELETION SUMMARY STATISTICS ===\n")
print(del_9bp_summary)

# Statistical analysis: Compare modified groups within each tool
cat("\n=== STATISTICAL ANALYSIS: MODIFIED GROUPS COMPARISON ===\n")

# Filter to only modified groups (exclude Unmodified)
del_9bp_modified <- del_9bp_data %>%
  filter(Treatment != "Unmodified")

# ALPINE: Compare all modified groups
alpine_data <- del_9bp_modified %>% filter(Tool == "ALPINE")
knockknock_data <- del_9bp_modified %>% filter(Tool == "knock-knock")

cat("\n--- ALPINE: Pairwise comparisons of modified groups ---\n")
if(nrow(alpine_data) > 0) {
  # One-way ANOVA first
  alpine_anova <- aov(Percentage_9bp ~ Treatment, data = alpine_data)
  alpine_anova_summary <- summary(alpine_anova)
  cat("ANOVA F-test p-value:", alpine_anova_summary[[1]][["Pr(>F)"]][1], "\n")
  
  # Pairwise t-tests with Bonferroni correction
  alpine_pairwise <- pairwise.t.test(alpine_data$Percentage_9bp,
                                     alpine_data$Treatment,
                                     p.adjust.method = "bonferroni")
  print(alpine_pairwise)
}

cat("\n--- knock-knock: Pairwise comparisons of modified groups ---\n")
if(nrow(knockknock_data) > 0) {
  # One-way ANOVA first
  knockknock_anova <- aov(Percentage_9bp ~ Treatment, data = knockknock_data)
  knockknock_anova_summary <- summary(knockknock_anova)
  cat("ANOVA F-test p-value:", knockknock_anova_summary[[1]][["Pr(>F)"]][1], "\n")
  
  # Pairwise t-tests with Bonferroni correction
  knockknock_pairwise <- pairwise.t.test(knockknock_data$Percentage_9bp,
                                         knockknock_data$Treatment,
                                         p.adjust.method = "bonferroni")
  print(knockknock_pairwise)
}

# Summary of significance within each tool
cat("\n--- Summary: Are modified groups significantly different within each tool? ---\n")
cat("ALPINE modified groups (Untreated, UV, 4M, 10M, 200M):\n")
if(exists("alpine_anova_summary")) {
  if(alpine_anova_summary[[1]][["Pr(>F)"]][1] < 0.05) {
    cat("  ANOVA p-value =", round(alpine_anova_summary[[1]][["Pr(>F)"]][1], 4), "- SIGNIFICANT differences exist\n")
    cat("  See pairwise comparisons above for specific group differences\n")
  } else {
    cat("  ANOVA p-value =", round(alpine_anova_summary[[1]][["Pr(>F)"]][1], 4), "- NO significant differences\n")
  }
}

cat("\nknock-knock modified groups (Untreated, UV, 4M, 10M, 200M):\n")
if(exists("knockknock_anova_summary")) {
  if(knockknock_anova_summary[[1]][["Pr(>F)"]][1] < 0.05) {
    cat("  ANOVA p-value =", round(knockknock_anova_summary[[1]][["Pr(>F)"]][1], 4), "- SIGNIFICANT differences exist\n")
    cat("  See pairwise comparisons above for specific group differences\n")
  } else {
    cat("  ANOVA p-value =", round(knockknock_anova_summary[[1]][["Pr(>F)"]][1], 4), "- NO significant differences\n")
  }
}
cat("- 9bp_deletion_comparison.pdf\n")
