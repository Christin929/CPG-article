# ==============================================================================
# SECTION 1: OVERALL ANALYSIS - DOMAIN SCORES
# ==============================================================================

# Load necessary libraries
library(e1071)
library(effsize)
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(car)
library(boot)
library(rstatix)

# 1. Data Import
file_path <- "E:/Nature communications/Uploaded files/Source Data.xlsx"
data <- read_excel(file_path, sheet = "Figure 4a")

# 2. Extract 6 AGREE II domains
metrics <- c("Scope and purpose", "Stakeholder involvement", 
             "Rigor of development", "Clarity of presentation", 
             "Applicability", "Editorial independence")

domain_data <- data %>% select(all_of(metrics))

# 3. Reshape to long format
data_long <- domain_data %>%
  pivot_longer(
    cols = everything(),
    names_to = "Domain",
    values_to = "Score"
  ) %>%
  mutate(Domain = factor(Domain))

# 4. Normality Test (Shapiro-Wilk)
check_normality_simple <- function(x, name) {
  x <- na.omit(x)
  sw_test <- shapiro.test(x)
  return(data.frame(
    Domain = name,
    SW_p = sw_test$p.value,
    Is_Normal = sw_test$p.value > 0.05,
    stringsAsFactors = FALSE
  ))
}

normality_results <- map_dfr(metrics, ~check_normality_simple(domain_data[[.]], .))
cat("\n--- Normality Test Results ---\n")
print(normality_results)

# 5. Homogeneity of Variance Test (Levene's Test)
levene_result <- leveneTest(Score ~ Domain, data = data_long)
cat("\n--- Homogeneity of Variance (Levene's Test) ---\n")
print(levene_result)

# 6. Statistical Method Selection
use_parametric <- all(normality_results$Is_Normal) && (levene_result$`Pr(>F)`[1] > 0.05)

if(use_parametric) {
  cat("\nData meet assumptions. Using One-way ANOVA.\n")
  anova_result <- aov(Score ~ Domain, data = data_long)
  anova_summary <- summary(anova_result)
  print(anova_summary)
  
  # Effect Size: Eta Squared (η²)
  ss_total <- sum(anova_summary[[1]]$`Sum Sq`)
  ss_between <- anova_summary[[1]]$`Sum Sq`[1]
  eta_squared <- ss_between / ss_total
  cat(sprintf("\nEffect Size η² = %.3f\n", eta_squared))
  
  # Post-hoc Test: Tukey HSD
  if(anova_summary[[1]]$`Pr(>F)`[1] < 0.05) {
    tukey_result <- TukeyHSD(anova_result)
    tukey_df <- as.data.frame(tukey_result$Domain) %>%
      rownames_to_column("Comparison") %>%
      mutate(Significance = case_when(
        `p adj` < 0.001 ~ "***", `p adj` < 0.01 ~ "**", 
        `p adj` < 0.05 ~ "*", TRUE ~ "ns"
      ))
    write_csv(tukey_df, "tukey_hsd_results.csv")
  }
} else {
  cat("\nData violate assumptions. Using Kruskal-Wallis Test.\n")
  kruskal_result <- kruskal.test(Score ~ Domain, data = data_long)
  print(kruskal_result)
  
  # Effect Size: Epsilon Squared (ε²)
  H <- kruskal_result$statistic
  k <- length(metrics)
  n <- nrow(data_long)
  epsilon_squared <- (H - k + 1) / (n - k)
  cat(sprintf("\nEffect Size ε² = %.3f\n", epsilon_squared))
  
  # Post-hoc Test: Dunn's Test
  if(kruskal_result$p.value < 0.05) {
    dunn_result <- dunn_test(Score ~ Domain, data = data_long, p.adjust.method = "bonferroni")
    dunn_result <- dunn_result %>%
      mutate(Significance = case_when(
        p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**", 
        p.adj < 0.05 ~ "*", TRUE ~ "ns"
      ))
    write_csv(dunn_result, "dunn_test_results.csv")
  }
}

# Bootstrap for Effect Size CI (95%)
get_epsilon_sq <- function(data, indices) {
  d <- data[indices, ]
  k_res <- kruskal.test(Score ~ Domain, data = d)
  H <- k_res$statistic
  k <- length(unique(d$Domain))
  n <- nrow(d)
  return(as.numeric((H - k + 1) / (n - k)))
}

cat("\nCalculating 95% CI for effect size via Bootstrap (R=1000)...\n")
boot_results <- boot(data = data_long, statistic = get_epsilon_sq, R = 1000)
boot_ci <- boot.ci(boot_results, type = "bca")
cat(sprintf("Effect Size ε² = %.3f; 95%% CI: [%.3f, %.3f]\n", 
            boot_results$t0, boot_ci$bca[4], boot_ci$bca[5]))

# 7. Descriptive Statistics
descriptive_stats <- data_long %>%
  group_by(Domain) %>%
  summarise(
    N = n(), Mean = mean(Score, na.rm = TRUE), SD = sd(Score, na.rm = TRUE),
    Median = median(Score, na.rm = TRUE), Q1 = quantile(Score, 0.25, na.rm = TRUE),
    Q3 = quantile(Score, 0.75, na.rm = TRUE), IQR = IQR(Score, na.rm = TRUE),
    Skewness = skewness(Score), Kurtosis = kurtosis(Score) - 3
  ) %>%
  arrange(desc(Median))

# 8. Plotting: Main Boxplot (Nature Communications Style)
nc_text_size <- 5; nc_line_size <- 0.2

domain_colors <- c("#728EB9", "#B1273E", "#E3A85C", "#A5C7D9", "#C96C76", "#EABD82")

boxplot_main <- ggplot(data_long, aes(x = Domain, y = Score, fill = Domain)) +
  geom_boxplot(alpha = 0.8, color = "grey40", linewidth = nc_line_size,
               outlier.shape = 16, outlier.size = 0.3, width = 0.6) +
  scale_fill_manual(values = domain_colors) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "AGREE II domain scores (%)") +
  theme_void() +
  theme(text = element_text(size = nc_text_size, family = "sans", color = "black"),
        axis.line.x = element_line(color = "black", linewidth = nc_line_size),
        axis.line.y = element_line(color = "black", linewidth = nc_line_size),
        axis.text.x = element_text(angle = 45, hjust = 1, size = nc_text_size),
        axis.text.y = element_text(size = nc_text_size),
        axis.title.y = element_text(size = nc_text_size, angle = 90, margin = margin(r = 5)),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_AGREE_II_Overall.pdf", boxplot_main, width = 88, height = 50, units = "mm", device = cairo_pdf)


# ==============================================================================
# SECTION 2: TEMPORAL TRENDS (1995-2022)
# ==============================================================================

# 1. Grouping by Time Periods
data_clean <- data %>%
  select(all_of(metrics), Year) %>%
  mutate(Time_Period = case_when(
    Year <= 2000 ~ "1995-2000", Year <= 2005 ~ "2001-2005",
    Year <= 2010 ~ "2006-2010", Year <= 2015 ~ "2011-2015",
    Year <= 2020 ~ "2016-2020", Year <= 2022 ~ "2021-2022",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Time_Period)) %>%
  mutate(Time_Period = factor(Time_Period, levels = c("1995-2000", "2001-2005", "2006-2010", 
                                                      "2011-2015", "2016-2020", "2021-2022")))

# 2. Enhanced Trend Analysis per Domain
data_long_trend <- data_clean %>%
  pivot_longer(cols = all_of(metrics), names_to = "Domain", values_to = "Score")

trend_results <- data_long_trend %>%
  group_by(Domain) %>%
  do({
    kw_res <- kruskal_test(., Score ~ Time_Period)
    eff_res <- kruskal_effsize(., Score ~ Time_Period)
    data.frame(H = kw_res$statistic, p = kw_res$p, epsilon_sq = eff_res$effsize)
  })

write_csv(trend_results, "Trend_Analysis_Results.csv")

# 3. Trend Plotting (Grouped Boxplot)
domain_order_fixed <- c("Clarity of presentation", "Scope and purpose", "Editorial independence",
                        "Rigor of development", "Stakeholder involvement", "Applicability")

time_colors <- c("#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026")

boxplot_trend <- ggplot(data_long_trend %>% mutate(Domain = factor(Domain, levels = domain_order_fixed)), 
                        aes(x = Domain, y = Score, fill = Time_Period)) +
  geom_boxplot(position = position_dodge(width = 0.85), alpha = 0.9, color = "grey30",
               linewidth = nc_line_size, outlier.shape = 16, outlier.size = 0.3, width = 0.7) +
  scale_fill_manual(values = time_colors) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "AGREE II domain scores (%)", fill = "Period") +
  theme_void() +
  theme(text = element_text(size = nc_text_size, family = "sans", color = "black"),
        axis.line.x = element_line(color = "black", linewidth = nc_line_size),
        axis.line.y = element_line(color = "black", linewidth = nc_line_size),
        axis.text.x = element_text(angle = 45, hjust = 1, size = nc_text_size),
        axis.text.y = element_text(size = nc_text_size),
        axis.title.y = element_text(size = nc_text_size, angle = 90, margin = margin(r = 5)),
        legend.position = "right",
        legend.text = element_text(size = nc_text_size),
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_Temporal_Trends.pdf", boxplot_trend, width = 130, height = 60, units = "mm", device = cairo_pdf)


# ==============================================================================
# SECTION 3: PRE-VS POST-2011 COMPARISON
# ==============================================================================

# 1. Period Classification
data_period <- data %>%
  select(all_of(metrics), Year) %>%
  mutate(Period = factor(ifelse(Year < 2011, "Pre-2011", "Post-2011"), levels = c("Pre-2011", "Post-2011")))

# 2. Statistical Comparison (Wilcoxon Rank-Sum & Cliff's Delta)
period_comparison <- data_period %>%
  pivot_longer(cols = all_of(metrics), names_to = "Domain", values_to = "Score") %>%
  group_by(Domain) %>%
  summarise(
    W = wilcox.test(Score ~ Period)$statistic,
    p = wilcox.test(Score ~ Period)$p.value,
    cd_res = list(cliff.delta(Score ~ Period)),
    cliff_delta = unlist(lapply(cd_res, function(x) x$estimate)),
    CI_lower = unlist(lapply(cd_res, function(x) x$conf.int[1])),
    CI_upper = unlist(lapply(cd_res, function(x) x$conf.int[2]))
  ) %>%
  select(-cd_res)

write_csv(period_comparison, "Pre_Post_2011_Statistical_Results.csv")

# 3. Comparison Plotting
boxplot_comp <- ggplot(data_period %>% pivot_longer(cols = all_of(metrics), names_to = "Domain", values_to = "Score") %>%
                         mutate(Domain = factor(Domain, levels = domain_order_fixed)), 
                       aes(x = Domain, y = Score, fill = Period)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9, color = "grey30",
               linewidth = nc_line_size, outlier.shape = 16, outlier.size = 0.3, width = 0.7) +
  scale_fill_manual(values = c("Pre-2011" = "#00008B", "Post-2011" = "#8B0000")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "AGREE II domain scores (%)") +
  theme_void() +
  theme(text = element_text(size = nc_text_size, family = "sans", color = "black"),
        axis.line.x = element_line(color = "black", linewidth = nc_line_size),
        axis.line.y = element_line(color = "black", linewidth = nc_line_size),
        axis.text.x = element_text(angle = 45, hjust = 1, size = nc_text_size),
        axis.text.y = element_text(size = nc_text_size),
        axis.title.y = element_text(size = nc_text_size, angle = 90, margin = margin(r = 5)),
        legend.position = "right",
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_Pre_Post_Comparison.pdf", boxplot_comp, width = 88, height = 50, units = "mm", device = cairo_pdf)

cat("\nAnalysis complete. All figures and statistical tables have been generated successfully.\n")

# ==============================================================================
# SECTION 4: ANALYSIS OF AGREE II DOMAIN SCORES BY REGION
# ==============================================================================

# Load necessary libraries
library(tidyverse)
library(readxl)
library(rstatix)
library(ggpubr)
library(boot)

# 1. Data Import
file_path <- "E:/Nature communications/Uploaded files/Source Data.xlsx"
data <- read_excel(file_path, sheet = "Figure 4b")

# 2. Variable Selection
metrics <- c("Scope and purpose", "Stakeholder involvement",
             "Rigor of development", "Clarity of presentation",
             "Applicability", "Editorial independence")

# Filter data and handle missing regions
region_data <- data %>%
  select(all_of(metrics), region) %>%
  filter(!is.na(region))

cat("--- Regional Distribution ---\n")
print(region_data %>% count(region) %>% arrange(desc(n)))

# 3. Reshape to long format
data_long_region <- region_data %>%
  pivot_longer(cols = all_of(metrics), names_to = "Domain", values_to = "Score") %>%
  mutate(Domain = factor(Domain), region = factor(region))

# 4. Descriptive Statistics by Region
desc_stats_region <- data_long_region %>%
  group_by(region, Domain) %>%
  summarise(
    N = n(), Median = median(Score, na.rm = TRUE),
    Q1 = quantile(Score, 0.25, na.rm = TRUE), Q3 = quantile(Score, 0.75, na.rm = TRUE),
    Mean = mean(Score, na.rm = TRUE), SD = sd(Score, na.rm = TRUE), .groups = 'drop'
  )

# 5. Statistical Inference: Inter-regional Comparisons
get_kw_eps <- function(data, indices) {
  d <- data[indices, ]
  res <- kruskal.test(Score ~ region, data = d)
  H <- res$statistic
  k <- length(unique(d$region))
  n <- nrow(d)
  return(as.numeric((H - k + 1) / (n - k)))
}

region_test_results <- data_long_region %>%
  group_by(Domain) %>%
  group_modify(~ {
    # Assumptions Check
    levene_p <- tryCatch(levene_test(.x, Score ~ region)$p, error = function(e) 0)
    n_groups <- length(unique(.x$region))
    
    if(levene_p > 0.05 && n_groups > 1) {
      # ANOVA Branch
      res_anova <- anova_test(.x, Score ~ region, detailed = TRUE)
      out_main <- tibble(
        Test = "ANOVA", Stat = res_anova$F, DF = paste0(res_anova$DFn, ",", res_anova$DFd),
        P_Value = res_anova$p, Effect_Size = res_anova$ges, ES_Type = "ges"
      )
      posthoc <- if(res_anova$p < 0.05) tukey_hsd(.x, Score ~ region) else NULL
    } else {
      # Kruskal-Wallis Branch
      res_kw <- kruskal_test(.x, Score ~ region)
      res_eff <- kruskal_effsize(.x, Score ~ region)
      set.seed(123)
      boot_res <- boot(data = .x, statistic = get_kw_eps, R = 1000)
      boot_ci <- tryCatch(boot.ci(boot_res, type = "perc"), error = function(e) list(percent=c(NA,NA,NA,NA,NA)))
      
      out_main <- tibble(
        Test = "Kruskal-Wallis", Stat = res_kw$statistic, DF = as.character(res_kw$df),
        P_Value = res_kw$p, Effect_Size = res_eff$effsize, ES_Type = "epsilon_squared",
        ES_CI_Low = boot_ci$percent[4], ES_CI_High = boot_ci$percent[5]
      )
      posthoc <- if(res_kw$p < 0.05) dunn_test(.x, Score ~ region, p.adjust.method = "bonferroni") else NULL
    }
    out_main %>% mutate(Posthoc = list(posthoc))
  }) %>% ungroup()

# 6. Export Results
write.csv(region_test_results %>% select(-Posthoc), "Region_Main_Stats.csv", row.names = FALSE)
all_region_posthoc <- region_test_results %>% select(Domain, Posthoc) %>% unnest(Posthoc)
write.csv(all_region_posthoc, "Region_Pairwise_Comparisons.csv", row.names = FALSE)

# 7. Visualization (NC Style)
region_colors <- c("Europe"="#3D7EB5", "America"="#C56F93", "Asia"="#D4C36A", 
                   "Cross region"="#93B7D6", "Oceania"="#E4B1B5", "Africa"="#E4D89E")
domain_order <- c("Clarity of presentation", "Scope and purpose", "Editorial independence",
                  "Rigor of development", "Stakeholder involvement", "Applicability")

plot_region <- ggplot(data_long_region %>% mutate(Domain = factor(Domain, levels = domain_order), 
                                                  region = factor(region, levels = names(region_colors))), 
                      aes(x = Domain, y = Score, fill = region)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9, color = "grey30", 
               linewidth = 0.2, outlier.shape = 16, outlier.size = 0.4, width = 0.7) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "AGREE II domain scores (%)", fill = "Region") +
  theme_void() +
  theme(text = element_text(size = 5, family = "sans"),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 5, angle = 90, margin = margin(r = 8)),
        legend.position = "right", legend.key.size = unit(3, "mm"),
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_Region_Analysis.pdf", plot_region, width = 130, height = 85, units = "mm", device = cairo_pdf)


# ==============================================================================
# SECTION 5: ANALYSIS OF AGREE II DOMAIN SCORES BY DISEASE
# ==============================================================================

# 1. Filter and Prepare Disease Data
disease_data <- data %>%
  select(all_of(metrics), diseases) %>%
  filter(!is.na(diseases)) %>%
  mutate(diseases = as.character(diseases))

# Select diseases with sufficient sample size (Threshold = 10)
sample_threshold <- 10
disease_counts <- disease_data %>% count(diseases) %>% filter(n >= sample_threshold)
disease_data_filtered <- disease_data %>% filter(diseases %in% disease_counts$diseases)

# 2. Reshape to long format
data_long_disease <- disease_data_filtered %>%
  pivot_longer(cols = all_of(metrics), names_to = "Domain", values_to = "Score") %>%
  mutate(Domain = factor(Domain), diseases = factor(diseases))

# 3. Enhanced Statistical Comparison between Diseases
get_kw_eps_disease <- function(data, indices) {
  d <- data[indices, ]
  if(length(unique(d$diseases)) < 2) return(NA)
  res <- kruskal.test(Score ~ diseases, data = d)
  return(as.numeric((res$statistic - length(unique(d$diseases)) + 1) / (nrow(d) - length(unique(d$diseases)))))
}

disease_test_results <- data_long_disease %>%
  filter(diseases != "Other") %>%
  group_by(Domain) %>%
  group_modify(~ {
    n_diseases <- length(unique(.x$diseases))
    levene_p <- tryCatch(levene_test(.x, Score ~ diseases)$p, error = function(e) 0)
    
    if(levene_p > 0.05 && n_diseases > 1) {
      res_anova <- anova_test(.x, Score ~ diseases, detailed = TRUE)
      out_main <- tibble(Test = "ANOVA", Stat = res_anova$F, P_Value = res_anova$p, Effect_Size = res_anova$ges)
      posthoc <- if(res_anova$p < 0.05) tukey_hsd(.x, Score ~ diseases) else NULL
    } else {
      res_kw <- kruskal_test(.x, Score ~ diseases)
      res_eff <- kruskal_effsize(.x, Score ~ diseases)
      set.seed(123)
      boot_res <- boot(data = .x, statistic = get_kw_eps_disease, R = 1000)
      boot_ci <- tryCatch(boot.ci(boot_res, type = "perc"), error = function(e) list(percent=c(NA,NA,NA,NA,NA)))
      
      out_main <- tibble(
        Test = "Kruskal-Wallis", Stat = res_kw$statistic, P_Value = res_kw$p, 
        Effect_Size = res_eff$effsize, ES_CI_Low = boot_ci$percent[4], ES_CI_High = boot_ci$percent[5]
      )
      posthoc <- if(res_kw$p < 0.05) dunn_test(.x, Score ~ diseases, p.adjust.method = "bonferroni") else NULL
    }
    out_main %>% mutate(Posthoc = list(posthoc))
  }) %>% ungroup()

# 4. Export Disease Statistics
write.csv(disease_test_results %>% select(-Posthoc), "Disease_Main_Test_Results.csv", row.names = FALSE)
all_disease_posthoc <- disease_test_results %>% select(Domain, Posthoc) %>% filter(!map_lgl(Posthoc, is.null)) %>% unnest(Posthoc)
write.csv(all_disease_posthoc, "Disease_Pairwise_Comparisons.csv", row.names = FALSE)

# 5. Disease Visualization (NC Style)
target_disease_order <- c("Neoplasms", "Cardiovascular diseases", "Digestive diseases", 
                          "Musculoskeletal disorders", "Diabetes and kidney diseases", 
                          "Neurological disorders", "Chronic respiratory diseases", 
                          "Skin and subcutaneous diseases", "Mental disorders", "Sense organ diseases")

disease_colors <- c("Neoplasms"="#728EB9", "Cardiovascular diseases"="#B1273E", "Digestive diseases"="#E3A85C",
                    "Musculoskeletal disorders"="#A5C7D9", "Diabetes and kidney diseases"="#C96C76",
                    "Neurological disorders"="#EABD82", "Chronic respiratory diseases"="#C3E5EB",
                    "Skin and subcutaneous diseases"="#F5CCCD", "Mental disorders"="#F1D6B0", "Sense organ diseases"="#C0E3D9")

plot_disease <- ggplot(data_long_disease %>% filter(diseases %in% target_disease_order) %>% 
                         mutate(Domain = factor(Domain, levels = domain_order), 
                                diseases = factor(diseases, levels = target_disease_order)), 
                       aes(x = Domain, y = Score, fill = diseases)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9, color = "grey30", 
               linewidth = 0.2, outlier.shape = 16, outlier.size = 0.3, width = 0.7) +
  scale_fill_manual(values = disease_colors) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = NULL, y = "AGREE II domain scores (%)", fill = "Diseases") +
  theme_void() +
  theme(text = element_text(size = 5, family = "sans"),
        axis.line = element_line(color = "black", linewidth = 0.2),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.text.y = element_text(size = 5),
        axis.title.y = element_text(size = 5, angle = 90, margin = margin(r = 5)),
        legend.position = "right", legend.key.size = unit(2.5, "mm"),
        plot.margin = margin(5, 5, 5, 5))

ggsave("Figure_Disease_Analysis.pdf", plot_disease, width = 130, height = 85, units = "mm", device = cairo_pdf)

cat("\nAnalysis complete. Results exported to CSV and PDF.\n")