# ==========================================
# 1. Environment Preparation
# ==========================================
install.packages(c("readxl", "tidyverse"))
library(readxl)
library(tidyverse)

# Define the core function for EAPC calculation
# This function extracts EAPC, 95% CI (Z-distribution), 
# test statistics (t), degrees of freedom (df), and exact P-values.
calculate_eapc <- function(data, group_col_name) {
  data %>%
    group_by(!!sym(group_col_name)) %>%
    group_modify(~ {
      # Fit log-linear regression model: ln(y) = beta*x + epsilon
      fit <- lm(log(Count) ~ Year, data = .x)
      s_fit <- summary(fit)
      coef_table <- s_fit$coefficients
      
      # Check if the model successfully estimated the slope for 'Year'
      if (nrow(coef_table) >= 2) {
        beta1  <- coef_table[2, 1]   # Regression coefficient (Slope)
        se1    <- coef_table[2, 2]   # Standard Error (SE)
        t_stat <- coef_table[2, 3]   # Test statistic (t-value)
        p_val  <- coef_table[2, 4]   # Exact P-value
        df_val <- fit$df.residual    # Degrees of Freedom (n-k)
        
        # Calculate 95% CI using a critical value of 1.96 (Standard Normal Distribution)
        ci_low_beta <- beta1 - 1.96 * se1
        ci_up_beta  <- beta1 + 1.96 * se1
        
        tibble(
          EAPC = round((exp(beta1) - 1) * 100, 2),
          CI_Lower = round((exp(ci_low_beta) - 1) * 100, 2),
          CI_Upper = round((exp(ci_up_beta) - 1) * 100, 2),
          Test_Statistic_t = round(t_stat, 3),
          DF = df_val,
          P_Value = p_val
        )
      } else {
        # Return NA if the model cannot be fitted
        tibble(EAPC = NA, CI_Lower = NA, CI_Upper = NA, 
               Test_Statistic_t = NA, DF = NA, P_Value = NA)
      }
    }) %>%
    ungroup() %>%
    # Add significance markers based on exact P-values
    mutate(Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ "ns"
    ))
}

# Define the local file path for data import
file_path <- "E:/Nature communications/Uploaded files/Source Data.xlsx"

# ==========================================
# 2. Batch Processing: Dataset A (Region)
# ==========================================
res_region <- read_excel(file_path, sheet = "Figure 2a") %>%
  pivot_longer(cols = -1, names_to = "Region", values_to = "Count") %>%
  mutate(Count = as.numeric(Count), Year = as.numeric(Year)) %>%
  # Exclude groups containing any zero counts to ensure valid log-transformation
  group_by(Region) %>% filter(!any(Count == 0)) %>% ungroup() %>%
  calculate_eapc("Region")

# ==========================================
# 3. Batch Processing: Dataset B (Disease)
# ==========================================
res_disease <- read_excel(file_path, sheet = "Figure 2c") %>%
  pivot_longer(cols = -1, names_to = "Disease", values_to = "Count") %>%
  mutate(Count = as.numeric(Count), Year = as.numeric(Year)) %>%
  group_by(Disease) %>% filter(!any(Count == 0)) %>% ungroup() %>%
  calculate_eapc("Disease")

# ==========================================
# 4. Batch Processing: Dataset C (Theme)
# ==========================================
res_theme <- read_excel(file_path, sheet = "Figure 2d") %>%
  pivot_longer(cols = -1, names_to = "Theme", values_to = "Count") %>%
  mutate(Count = as.numeric(Count), Year = as.numeric(Year)) %>%
  group_by(Theme) %>% filter(!any(Count == 0)) %>% ungroup() %>%
  calculate_eapc("Theme")

# ==========================================
# 5. Exporting Results
# ==========================================
write.csv(res_region, "EAPC_Analysis_Region.csv", row.names = FALSE)
write.csv(res_disease, "EAPC_Analysis_Disease.csv", row.names = FALSE)
write.csv(res_theme, "EAPC_Analysis_Theme.csv", row.names = FALSE)

cat("Analysis complete. EAPC statistics for all datasets have been exported successfully.\n")