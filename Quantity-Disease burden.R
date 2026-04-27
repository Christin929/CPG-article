# ==============================================================================
# SENSITIVITY ANALYSIS: COMPLETE CASE VS. PCA IMPUTATION
# Script prepared for Nature Communications code availability
# ==============================================================================

# 1. Environment Setup
if(!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, geepack, broom, missMDA, boot)

# Path configuration
FILE_PATH  <- "E:/Nature communications/Uploaded files/Source Data.xlsx"
SHEET_NAME <- "Table 1"
M_COUNT    <- 20    # Number of imputations
SEED       <- 123   # For reproducibility
IV         <- "count"  # Independent Variable
COVARS     <- c("populationdensity", "urbanpopulationpercent", "GDPpercapita", 
                "HumanDevelopmentIndex", "year")
OUTCOMES   <- c("DALYs", "Deaths")
ALL_VARS   <- c(IV, COVARS, OUTCOMES)

set.seed(SEED)

# 2. Data Preprocessing
df_raw <- read_excel(FILE_PATH, sheet = SHEET_NAME)

df_clean <- df_raw %>%
  mutate(across(all_of(ALL_VARS), ~ as.numeric(gsub("[^0-9.-]", "", as.character(.))))) %>%
  filter(!is.na(country)) %>%
  mutate(country_id = as.numeric(as.factor(country)))

# Standardize time-series structure
df_std <- df_clean %>%
  group_by(country) %>%
  complete(year = min(df_clean$year, na.rm=TRUE):max(df_clean$year, na.rm=TRUE)) %>%
  fill(country_id, .direction = "downup") %>%
  ungroup() %>%
  arrange(country_id, year)

# 3. GEE Modeling Engine
run_gee_model <- function(data, outcome_name, lag_year) {
  if(!"country_id" %in% names(data)) return(NULL)
  
  dat_lag <- data %>%
    group_by(country_id) %>%
    arrange(year) %>%
    mutate(
      y_target    = lead(.data[[outcome_name]], n = lag_year),
      count_base  = .data[[IV]],
      across(all_of(COVARS), ~ identity(.x), .names = "{.col}_base")
    ) %>%
    ungroup() %>%
    filter(!is.na(y_target), y_target > 0) %>% 
    arrange(country_id, year) %>% 
    select(country_id, count_base, y_target, ends_with("_base")) %>%
    na.omit()
  
  if(nrow(dat_lag) < 25) return(NULL)
  
  formula_str <- paste0("y_target ~ count_base + ", 
                        paste(paste0(COVARS, "_base"), collapse = " + "))
  
  m <- tryCatch(
    geeglm(as.formula(formula_str), data = dat_lag, id = country_id, 
           family = Gamma(link = "log"), corstr = "independence",
           control = geese.control(maxit = 50)),
    error = function(e) NULL
  )
  
  if(is.null(m)) return(NULL)
  res <- summary(m)$coefficients
  
  if("count_base" %in% rownames(res)) {
    return(list(beta = res["count_base", 1], se = res["count_base", 2], n_obs = nrow(dat_lag)))
  }
  return(NULL)
}

# 4. Statistical Workflow (Rubin's Rules for Pooling)
analyze_workflow <- function(data_list, method_label) {
  final_res <- data.frame()
  for(out in OUTCOMES) {
    for(lag in 1:5) { 
      results_list <- lapply(data_list, function(d) run_gee_model(d, out, lag))
      valid_res <- Filter(Negate(is.null), results_list)
      
      if(length(valid_res) >= 1) {
        betas <- sapply(valid_res, function(x) x$beta)
        ses   <- sapply(valid_res, function(x) x$se)
        n_avg <- mean(sapply(valid_res, function(x) x$n_obs))
        
        # Pooling Estimates
        m <- length(betas)
        b_pool <- mean(betas)
        vw <- mean(ses^2) # Within-imputation variance
        vb <- if(m > 1) var(betas) else 0 # Between-imputation variance
        vt <- vw + vb + (vb/m) # Total variance
        se_pool <- sqrt(vt)
        
        final_res <- rbind(final_res, data.frame(
          Method  = method_label,
          Outcome = out,
          Lag     = lag,
          RR      = exp(b_pool),
          LL      = exp(b_pool - 1.96*se_pool),
          UL      = exp(b_pool + 1.96*se_pool),
          P_Value = 2 * (1 - pnorm(abs(b_pool/se_pool))),
          N_Obs   = round(n_avg)
        ))
      }
    }
  }
  return(final_res)
}

# 5. Execute Analysis
# Strategy 1: Complete Case Analysis (CCA)
res_cca <- analyze_workflow(list(df_std), "1_CCA")

# Strategy 2: Multiple Imputation with PCA (MIPCA)
df_pca_input <- as.data.frame(df_std[, ALL_VARS])
imp_pca_multi <- missMDA::MIPCA(df_pca_input, ncp = 2, nboot = M_COUNT) 

# Re-attach country IDs to imputed sets
pca_list <- lapply(imp_pca_multi$res.MI, function(d) {
  d$country_id <- df_std$country_id
  d$year <- df_std$year
  return(d)
})
res_pca <- analyze_workflow(pca_list, "2_MIPCA")

# 6. Summary and Export
sensitivity_summary <- bind_rows(res_cca, res_pca)
write.csv(sensitivity_summary, "Sensitivity_Analysis_Results.csv", row.names = FALSE)

cat("--- Analysis Complete ---\n")
print(sensitivity_summary)