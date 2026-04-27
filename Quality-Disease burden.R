# ==============================================================================
# SENSITIVITY ANALYSIS: COMPLETE CASE VS. PCA IMPUTATION (QUALITY COHORT)
# Prepared for Nature Communications Code Availability
# ==============================================================================

# 1. Environment and Package Loading
if(!require(pacman)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, geepack, broom, missMDA)

# Path and Parameter Configuration
FILE_PATH  <- "E:/Nature communications/Uploaded files/Source Data.xlsx" # Ensure the file is in the working directory
SHEET_NAME <- "Table 2-Scope and purpose"        # Standardized sheet name
M_COUNT    <- 20                         # Number of multiple imputations
SEED       <- 123                        # For statistical reproducibility

IV         <- "highpercent"              # Independent Variable (Exposure)
COVARS     <- c("populationdensity", "urbanpopulationpercent", "GDPpercapita", 
                "HumanDevelopmentIndex", "year")
OUTCOMES   <- c("DALYs", "Deaths")
ALL_VARS   <- c(IV, COVARS, OUTCOMES)

set.seed(SEED)

# 2. Data Cleaning and Standardization
df_raw <- read_excel(FILE_PATH, sheet = "All domains＞70%")

df_clean <- df_raw %>%
  mutate(
    highpercent = as.numeric(gsub("%", "", as.character(na_if(na_if(highpercent, "na"), "NA")))),
    across(all_of(setdiff(ALL_VARS, "highpercent")), ~ as.numeric(na_if(as.character(.), "na")))
  ) %>%
  filter(!is.na(country)) %>%
  mutate(country_id = as.numeric(as.factor(country)))

# Create continuous time-series for panel data
global_years <- min(df_clean$year, na.rm=TRUE):max(df_clean$year, na.rm=TRUE)
df_std <- df_clean %>%
  group_by(country) %>%
  complete(year = global_years) %>%
  fill(country_id, .direction = "downup") %>%
  ungroup() %>%
  arrange(country_id, year)

# 3. GEE Modeling Engine (Gamma distribution with Log link)
run_gee_model <- function(data, outcome_name, lag_year) {
  data <- as.data.frame(data)
  
  dat_lag <- data %>%
    group_by(country_id) %>%
    arrange(year) %>%
    mutate(
      y_target = lead(.data[[outcome_name]], n = lag_year),
      exposure_val = .data[[IV]],
      across(all_of(COVARS), ~ identity(.x), .names = "{.col}_base")
    ) %>%
    ungroup() %>%
    filter(!is.na(y_target), y_target > 0) %>%
    select(country_id, y_target, exposure_val, ends_with("_base")) %>%
    na.omit() 
  
  if(nrow(dat_lag) < 15) return(NULL) 
  
  # Construct formula: outcome ~ exposure + covariates
  formula_str <- paste0("y_target ~ exposure_val + ", 
                        paste(paste0(COVARS, "_base"), collapse = " + "))
  
  m <- tryCatch(
    geeglm(as.formula(formula_str), data = dat_lag, id = country_id, 
           family = Gamma(link = "log"), corstr = "independence"),
    error = function(e) NULL
  )
  
  if(is.null(m)) return(NULL)
  res <- summary(m)$coefficients
  
  if("exposure_val" %in% rownames(res)) {
    return(list(beta = res["exposure_val", 1], se = res["exposure_val", 2], n_obs = nrow(dat_lag)))
  }
  return(NULL)
}

# 4. Statistical Workflow (Applying Rubin's Rules for pooling)
analyze_workflow <- function(data_list, method_label) {
  final_res <- data.frame()
  for(out in OUTCOMES) {
    for(lag in 1:5) { # Analyzing 1 to 5-year lag effects
      results_list <- lapply(data_list, function(d) run_gee_model(d, out, lag))
      valid_res <- Filter(function(x) !is.null(x), results_list)
      
      if(length(valid_res) >= 1) {
        betas <- sapply(valid_res, function(x) x$beta)
        ses   <- sapply(valid_res, function(x) x$se)
        n_avg <- mean(sapply(valid_res, function(x) x$n_obs))
        
        # Rubin's Rules Pooling
        m <- length(betas)
        b_pool <- mean(betas)
        vw <- mean(ses^2) 
        vb <- if(m > 1) var(betas) else 0 
        vt <- vw + vb + (vb/m) 
        se_pool <- sqrt(vt)
        
        final_res <- rbind(final_res, data.frame(
          Method  = method_label,
          Outcome = out,
          Lag     = lag,
          RR      = exp(b_pool),
          LL      = exp(b_pool - 1.96 * se_pool),
          UL      = exp(b_pool + 1.96 * se_pool),
          P_Value = 2 * (1 - pnorm(abs(b_pool/se_pool))),
          N_Obs   = round(n_avg)
        ))
      }
    }
  }
  return(final_res)
}

# 5. Analysis Execution
# Strategy 1: Complete Case Analysis (CCA)
res_cca <- analyze_workflow(list(df_std), "1_CCA")

# Strategy 2: Multiple Imputation with PCA (MIPCA)
res_pca <- data.frame()
try({
  imp_pca_multi <- missMDA::MIPCA(df_std[, ALL_VARS], ncp = 2, nboot = M_COUNT)
  
  pca_list <- lapply(imp_pca_multi$res.MI, function(x) {
    d <- as.data.frame(x)
    d$country_id <- df_std$country_id
    d$year       <- df_std$year
    return(d)
  })
  
  res_pca <- analyze_workflow(pca_list, "2_MIPCA")
}, silent = FALSE)

# 6. Final Data Integration and Export
sensitivity_summary <- bind_rows(res_cca, res_pca)
write.csv(sensitivity_summary, "Sensitivity_Analysis_Quality_Results.csv", row.names = FALSE)

cat("--- Statistical Analysis for Quality Cohort Complete ---\n")
print(sensitivity_summary)