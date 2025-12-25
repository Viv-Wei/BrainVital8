# ==============================================================================
# Association Analysis of BrainVital8_aligned vs. dementia_status
# Analysis Pipeline: Data Loading → ID Standardization → Variable Imputation → 
#                   Covariate Scoring → Quartile Stratification → Logistic Regression
# ==============================================================================

# ==================== 1. Load Required Libraries ====================
library(readxl)        
library(dplyr)         
library(tidyr)         
library(stringr)       
library(broom)         
library(purrr)         
library(writexl)       
library(haven)         
library(magrittr)      
# ==================== 2. Set Output Path ====================
# Define output path for saving results
out_path <- file.path('./output/')
cat("output path:", out_path, "\n")
# ==================== 3. Data Loading ====================
cat("Loading datasets...\n")
# Specify ID column names for each dataset
brain_id_col <- "ID"          # ID column name in brain age dataset
stata_id_col <- "prim_key"    # ID column name in Stata dataset
lancet_id_col <- "ID1"        # ID column name in LANCET dataset
libra2_id_col <- "ID1"        # ID column name in LIBRA2 dataset
# Load datasets
brain_age_data <- read.csv(
  "./input/LASI_merged_final.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
lancet_data <- read.csv(
  "./input/lancet_total.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
libra2_data <- read.csv(
  "./input/libra2_total.csv",
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)
stata_data <- read_dta("./input/H_LASI_a3.dta")
# Validate existence of ID columns
if (!brain_id_col %in% colnames(brain_age_data)) stop(paste("Brain age dataset missing ID column:", brain_id_col))
if (!stata_id_col %in% colnames(stata_data)) stop(paste("Stata dataset missing ID column:", stata_id_col))
if (!lancet_id_col %in% colnames(lancet_data)) stop(paste("LANCET dataset missing ID column:", lancet_id_col))
if (!libra2_id_col %in% colnames(libra2_data)) stop(paste("LIBRA2 dataset missing ID column:", libra2_id_col))
# Define and validate exposure and outcome variables
exposure_var <- "BrainVital8_aligned"  # Primary exposure variable
outcome_var <- "dementia_status"       # Primary outcome variable (1=case, 0=control)
if (!exposure_var %in% colnames(brain_age_data)) stop(paste("Missing exposure variable:", exposure_var))
if (!outcome_var %in% colnames(brain_age_data)) stop(paste("Missing outcome variable:", outcome_var))
# Print dataset dimensions
cat("Brain age dataset dimensions:", dim(brain_age_data), "\n")
cat("LANCET dataset dimensions:", dim(lancet_data), "\n")
cat("LIBRA2 dataset dimensions:", dim(libra2_data), "\n")
cat("Stata dataset (H_LASI_a3.dta) dimensions:", dim(stata_data), "\n")
# ==================== 4. ID Standardization and Variable Supplementation ====================
cat("\nPerforming ID standardization and variable supplementation...\n")
# 4.1 Rename all ID columns to "ID" for consistency
brain_age_data <- rename(brain_age_data, ID = !!sym(brain_id_col))
stata_data <- rename(stata_data, ID = !!sym(stata_id_col))
lancet_data <- rename(lancet_data, ID = !!sym(lancet_id_col))
libra2_data <- rename(libra2_data, ID = !!sym(libra2_id_col))
# 4.2 Standardize ID format to character
brain_age_data$ID <- as.character(brain_age_data$ID)
stata_data$ID <- as.character(stata_data$ID)
lancet_data$ID <- as.character(lancet_data$ID)
libra2_data$ID <- as.character(libra2_data$ID)
# 4.3 Supplement missing variables from Stata dataset to brain age dataset
missing_vars <- setdiff(colnames(stata_data), colnames(brain_age_data))
if (length(missing_vars) > 0) {
  cat("Supplementing missing variables from Stata dataset:", paste(missing_vars, collapse = ", "), "\n")
  stata_supplement <- stata_data %>% select(ID, all_of(missing_vars))
  brain_age_data <- brain_age_data %>% left_join(stata_supplement, by = "ID")
  cat("Brain age dataset dimensions after supplementation:", dim(brain_age_data), "\n")
} else {
  cat("No missing variables to supplement from Stata dataset.\n")
}
# ==================== 5. Data Preprocessing and Merging ====================
cat("\nMerging datasets...\n")
# Function to extract total score column from dataset
extract_score_col <- function(data, score_type, data_name) {
  score_cols <- grep(paste0(score_type, "|score|total"), colnames(data), value = TRUE, ignore.case = TRUE)
  if (length(score_cols) == 0) score_cols <- setdiff(colnames(data), "ID")
  if (length(score_cols) == 0) stop(paste(data_name, "no score column found"))
  return(score_cols[1])
}
# Extract score columns for LANCET and LIBRA2 datasets
lancet_score_col <- extract_score_col(lancet_data, "LANCET", "LANCET dataset")
libra2_score_col <- extract_score_col(libra2_data, "LIBRA2", "LIBRA2 dataset")
cat("LANCET score column:", lancet_score_col, "\n")
cat("LIBRA2 score column:", libra2_score_col, "\n")
# Create subsets with standardized column names for merging
lancet_subset <- lancet_data %>% select(ID, LANCET_total = !!sym(lancet_score_col))
libra2_subset <- libra2_data %>% select(ID, LIBRA2_total = !!sym(libra2_score_col))
# Merge all datasets and filter out rows with missing outcome values
combined_data <- brain_age_data %>%
  left_join(lancet_subset, by = "ID") %>%
  left_join(libra2_subset, by = "ID") %>%
  filter(!is.na(!!sym(outcome_var)))
cat("Merged dataset dimensions:", dim(combined_data), "\n")
# ==================== 6. Covariate Missing Value Imputation ====================
cat("\nImputing missing values for covariates...\n")
# 6.1 Impute categorical variables with 9 (missing indicator)
categorical_vars <- c("ragender", "r1hibpe", "r1rxhibp")  # Categorical covariates: sex, hypertension-related
for (var in categorical_vars) {
  if (var %in% colnames(combined_data)) {
    missing_count <- sum(is.na(combined_data[[var]]))
    combined_data[[var]] <- replace(combined_data[[var]], is.na(combined_data[[var]]), 9)
    cat(sprintf("Categorical variable %s: %d missing values imputed with 9\n", var, missing_count))
  }
}
# 6.2 Impute continuous variables with median
continuous_vars <- c("r1agey", "r1mbmi", "r1drinkb", "r1systo", "r1diasto")  # Continuous covariates: age, BMI, alcohol consumption, blood pressure
for (var in continuous_vars) {
  if (var %in% colnames(combined_data)) {
    missing_count <- sum(is.na(combined_data[[var]]))
    median_val <- median(combined_data[[var]], na.rm = TRUE)
    combined_data[[var]] <- replace(combined_data[[var]], is.na(combined_data[[var]]), median_val)
    cat(sprintf("Continuous variable %s: %d missing values imputed with median (%.2f)\n", var, missing_count, median_val))
  }
}
# 6.3 Remove rows with missing exposure variable
exposure_missing <- sum(is.na(combined_data[[exposure_var]]))
combined_data <- filter(combined_data, !is.na(!!sym(exposure_var)))
cat(sprintf("Exposure variable %s: %d missing values removed, remaining sample size: %d\n", exposure_var, exposure_missing, nrow(combined_data)))
# Record final sample size
final_sample_size <- nrow(combined_data)
cat("Final analysis sample size:", final_sample_size, "\n")
# ==================== 7. Covariate Scoring and Quartile Stratification ====================
cat("\nCalculating covariate scores and exposure quartiles...\n")
# 7.1 Function to calculate hypertension status (0=non-hypertensive, 1=hypertensive)
calculate_hypertension_score <- function(r1hibpe, r1systo, r1diasto, r1rxhibp) {
  case_when(
    r1hibpe == 1 | r1systo > 140 | r1diasto > 90 | r1rxhibp == 1 ~ 1,  # Hypertensive
    TRUE ~ 0                                                           # Non-hypertensive (0-1 encoding for regression)
  )
}
# 7.2 Calculate hypertension score and stratify exposure into quartiles
combined_data <- combined_data %>%
  mutate(
    # Hypertension status (impute missing with 9)
    hypertension = replace(calculate_hypertension_score(r1hibpe, r1systo, r1diasto, r1rxhibp), 
                          is.na(calculate_hypertension_score(r1hibpe, r1systo, r1diasto, r1rxhibp)), 9),
    
    # Exposure quartiles (Q1: lowest 25%, Q4: highest 25%)
    exposure_q1 = quantile(!!sym(exposure_var), 0.25, na.rm = TRUE),
    exposure_q2 = quantile(!!sym(exposure_var), 0.5, na.rm = TRUE),
    exposure_q3 = quantile(!!sym(exposure_var), 0.75, na.rm = TRUE),
    exposure_quartile = case_when(
      !!sym(exposure_var) <= exposure_q1 ~ 1,
      !!sym(exposure_var) <= exposure_q2 ~ 2,
      !!sym(exposure_var) <= exposure_q3 ~ 3,
      TRUE ~ 4
    ),
    exposure_quartile_f = factor(exposure_quartile, levels = 1:4, labels = c("Q1", "Q2", "Q3", "Q4")),
    
    # Age stratification (<65 years vs. ≥65 years)
    age_group = factor(
      ifelse(r1agey < 65, 1, 2),
      levels = c(1, 2),
      labels = c("<65y", "≥65y")
    ),
    
    # Remove temporary quartile variables
    exposure_q1 = NULL, exposure_q2 = NULL, exposure_q3 = NULL
  )
# Create stratified datasets by gender and age
male_data <- combined_data %>% filter(ragender == 1)  # Male subset
female_data <- combined_data %>% filter(ragender == 2) # Female subset
age_lt65_data <- combined_data %>% filter(r1agey < 65) # <65 years subset
age_ge65_data <- combined_data %>% filter(r1agey >= 65) # ≥65 years subset
# 7.3 Calculate sample size and event count by exposure quartile
calculate_stratified_stats <- function(data, strata_name) {
  data %>%
    group_by(exposure_quartile_f) %>%
    summarise(
      sample_size = n(),
      event_count = sum(!!sym(outcome_var) == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(strata = strata_name)
}
calculate_total_stats <- function(data, strata_name) {
  data %>%
    summarise(
      sample_size = n(),
      event_count = sum(!!sym(outcome_var) == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      strata = strata_name,
      exposure_quartile_f = "Total"
    )
}
# Global descriptive statistics
global_stats <- calculate_stratified_stats(combined_data, "Global")
global_total <- calculate_total_stats(combined_data, "Global")
# Gender-stratified descriptive statistics
male_stats <- calculate_stratified_stats(male_data, "Male")
male_total <- calculate_total_stats(male_data, "Male")
female_stats <- calculate_stratified_stats(female_data, "Female")
female_total <- calculate_total_stats(female_data, "Female")
# Age-stratified descriptive statistics
age_lt65_stats <- calculate_stratified_stats(age_lt65_data, "<65y")
age_lt65_total <- calculate_total_stats(age_lt65_data, "<65y")
age_ge65_stats <- calculate_stratified_stats(age_ge65_data, "≥65y")
age_ge65_total <- calculate_total_stats(age_ge65_data, "≥65y")
# Combine all descriptive statistics
descriptive_stats <- bind_rows(
  global_total, global_stats,
  male_total, male_stats,
  female_total, female_stats,
  age_lt65_total, age_lt65_stats,
  age_ge65_total, age_ge65_stats
) %>%
  mutate(Case_N = sprintf("%d/%d", event_count, sample_size))

# 7.4 Create Case_N for continuous exposure (overall sample in each stratum)
create_continuous_case_n <- function(data_list, strata_names) {
  map2_dfr(data_list, strata_names, function(data, strata) {
    total_event <- sum(data[[outcome_var]] == 1, na.rm = TRUE)
    total_n <- nrow(data)
    tibble(
      strata = strata,
      exposure_quartile_f = "Continuous",
      Case_N = sprintf("%d/%d", total_event, total_n)
    )
  })
}
# Define stratum datasets and corresponding names
strata_datasets <- list(combined_data, male_data, female_data, age_lt65_data, age_ge65_data)
strata_names <- c("Global", "Male", "Female", "<65y", "≥65y")
# Generate Case_N for continuous exposure
continuous_case_n <- create_continuous_case_n(strata_datasets, strata_names)
# Merge continuous Case_N into descriptive statistics
descriptive_stats <- bind_rows(descriptive_stats, continuous_case_n)
# ==================== 8. Define Logistic Regression Function ====================
run_logistic_models <- function(data, formula_str_base, model_label, stratify_label) {
  tryCatch({
    if (nrow(data) == 0) stop("Empty dataset")
    
    # 8.1 Model 1: Continuous exposure variable
    formula_continuous <- paste(formula_str_base, exposure_var, sep = " + ")
    model_continuous <- glm(
      as.formula(formula_continuous),
      family = binomial(link = "logit"),
      data = data
    )
    
    # 8.2 Model 2: Categorical exposure (quartiles, Q1 as reference)
    formula_quartile <- paste(formula_str_base, "exposure_quartile_f", sep = " + ")
    model_quartile <- glm(
      as.formula(formula_quartile),
      family = binomial(link = "logit"),
      data = data
    )
    
    # 8.3 Tidy continuous exposure model results
    results_continuous <- tidy(model_continuous, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(str_detect(term, exposure_var)) %>%
      mutate(
        strata = stratify_label,
        model = model_label,
        exposure_category = "Continuous",
        OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
        P_value = case_when(
          is.na(p.value) ~ "NA",
          p.value < 0.001 ~ sprintf("%.2e", p.value),
          p.value < 0.05 ~ sprintf("%.4f", p.value),
          TRUE ~ sprintf("%.2f", p.value)
        )
      ) %>%
      select(strata, model, exposure_category, OR_95CI, P_value, p.value)
    
    # 8.4 Tidy categorical exposure (quartile) model results
    results_quartile <- tidy(model_quartile, conf.int = TRUE, exponentiate = TRUE) %>%
      filter(str_detect(term, "exposure_quartile_f")) %>%
      mutate(
        strata = stratify_label,
        model = model_label,
        OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
        P_value = case_when(
          is.na(p.value) ~ "NA",
          p.value < 0.001 ~ sprintf("%.2e", p.value),
          p.value < 0.05 ~ sprintf("%.4f", p.value),
          TRUE ~ sprintf("%.2f", p.value)
        ),
        exposure_category = case_when(
          str_detect(term, "Q2") ~ "Q2",
          str_detect(term, "Q3") ~ "Q3",
          str_detect(term, "Q4") ~ "Q4",
          TRUE ~ "Q1"
        )
      ) %>%
      select(strata, model, exposure_category, OR_95CI, P_value, p.value)
    
    # 8.5 Add reference group (Q1)
    reference_row <- data.frame(
      strata = stratify_label,
      model = model_label,
      exposure_category = "Q1",
      OR_95CI = "Reference",
      P_value = "-",
      p.value = NA,
      stringsAsFactors = FALSE
    )
    
    # 8.6 Combine continuous and quartile results
    results_combined <- bind_rows(results_continuous, reference_row, results_quartile) %>%
      filter(exposure_category %in% c("Continuous", "Q1", "Q2", "Q3", "Q4"))
    
    list(results = results_combined, success = TRUE)
  }, error = function(e) {
    cat("Model", model_label, "failed:", e$message, "\n")
    list(success = FALSE)
  })
}
# ==================== 9. Fit Target Logistic Regression Models ====================
cat("\nFitting target logistic regression models...\n")
# Define base covariates (age, sex, body mass index, alcohol consumption, hypertension)
# Variable mapping:
# age → r1agey (continuous age)
# sex → ragender (categorical sex)
# body mass index → r1mbmi (continuous BMI)
# alcohol consumption → r1drinkb (continuous alcohol consumption)
# hypertension → hypertension (binary hypertension status: 0=non-hypertensive, 1=hypertensive)
base_covariates <- "r1agey + ragender + r1mbmi + r1drinkb + hypertension"
# Initialize model list
model_list <- list()
# 9.1 Main model (Global + base covariates)
model_list[["main_global"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "Main Model",
  stratify_label = "Global"
)
# 9.2 LANCET-adjusted model (Global + base covariates + LANCET_total)
model_list[["lancet_adjusted"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~ LANCET_total +", base_covariates),
  model_label = "LANCET-adjusted",
  stratify_label = "Global"
)
# 9.3 LIBRA2-adjusted model (Global + base covariates + LIBRA2_total)
model_list[["libra2_adjusted"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~ LIBRA2_total +", base_covariates),
  model_label = "LIBRA2-adjusted",
  stratify_label = "Global"
)
# 9.4 Gender-stratified models (Male/Female + base covariates)
model_list[["male_stratified"]] <- run_logistic_models(
  data = male_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "Male-stratified",
  stratify_label = "Male"
)
model_list[["female_stratified"]] <- run_logistic_models(
  data = female_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "Female-stratified",
  stratify_label = "Female"
)
# 9.5 Age-stratified models (<65y/≥65y + base covariates)
model_list[["age_lt65_stratified"]] <- run_logistic_models(
  data = age_lt65_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "<65y-stratified",
  stratify_label = "<65y"
)
model_list[["age_ge65_stratified"]] <- run_logistic_models(
  data = age_ge65_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "≥65y-stratified",
  stratify_label = "≥65y"
)
# ==================== 10. Combine Results and Generate Output Table ====================
cat("\nCombining results and generating output table...\n")
# Filter successfully fitted models
successful_models <- model_list[sapply(model_list, function(x) x$success)]
if (length(successful_models) == 0) stop("No successful models fitted")
# Combine results from all successful models
model_results <- map_dfr(successful_models, function(x) x$results)
# Merge with descriptive statistics (Case_N)
final_results <- model_results %>%
  left_join(
    descriptive_stats %>% select(strata = strata, exposure_quartile_f, Case_N),
    by = c("strata", "exposure_category" = "exposure_quartile_f")
  ) %>%
  # Reorder columns to Nature Genetics style
  select(
    Model = model,
    Strata = strata,
    Exposure_Category = exposure_category,
    Case_N = Case_N,
    OR_95CI = OR_95CI,
    P_value = P_value
  ) %>%
  # Sort by model and exposure category
  arrange(
    factor(Model, levels = c("Main Model", "LANCET-adjusted", "LIBRA2-adjusted", 
                            "Male-stratified", "Female-stratified", "<65y-stratified", "≥65y-stratified")),
    factor(Exposure_Category, levels = c("Continuous", "Q1", "Q2", "Q3", "Q4"))
  )
# ==================== 11. Save Results ====================
output_filename <- paste0("Exposure_", exposure_var, "_Outcome_", outcome_var, "_Results.xlsx")
output_path <- file.path(out_path, output_filename)
write_xlsx(final_results, output_path)
cat("✅ Results saved to:", output_path, "\n")
# Preview first 20 rows of results
cat("\n=== Results Preview (First 20 Rows) ===\n")
print(head(final_results, 20), row.names = FALSE)
# ==================== 12. Generate Summary Statistics ====================
cat("\n=== Summary Statistics ===\n")
cat("Total sample size:", final_sample_size, "\n")
cat("Total cases:", sum(combined_data[[outcome_var]] == 1, na.rm = TRUE), "\n")
cat("Total controls:", sum(combined_data[[outcome_var]] == 0, na.rm = TRUE), "\n")
cat("Number of successfully fitted models:", length(successful_models), "\n")