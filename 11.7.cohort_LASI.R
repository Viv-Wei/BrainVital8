# ==============================================================================
# BrainVital8_aligned vs. dementia_status Association Analysis
# Analysis Pipeline: Data Loading → ID Standardization → Variable Imputation → 
#                   Covariate Scoring → Quartile Stratification → Logistic Regression
# ==============================================================================

# ==================== 1. Load Required Libraries ====================
library(readxl)        # Read Excel files
library(dplyr)         # Data manipulation
library(tidyr)         # Data tidying
library(stringr)       # String processing
library(broom)         # Model result tidying
library(purrr)         # Functional programming
library(writexl)       # Write Excel files
library(haven)         # Read Stata (.dta) files
library(magrittr)      # Enhanced pipe operations



# ==================== 3. Data Loading ====================
cat("Loading datasets...\n")

# Specify ID column names for each dataset
brain_id_col <- "ID"          # Brain age dataset ID column
stata_id_col <- "prim_key"    # Stata dataset ID column
lancet_id_col <- "ID1"        # LANCET dataset ID column
libra2_id_col <- "ID1"        # LIBRA2 dataset ID column

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

# Validate existence of exposure and outcome variables
exposure_var <- "BrainVital8_aligned"  # Primary exposure variable
outcome_var <- "dementia_status"       # Primary outcome variable (1=case, 0=control)
if (!exposure_var %in% colnames(brain_age_data)) stop(paste("Missing exposure variable:", exposure_var))
if (!outcome_var %in% colnames(brain_age_data)) stop(paste("Missing outcome variable:", outcome_var))

# Print dataset dimensions
cat("Brain age dataset dimensions:", dim(brain_age_data), "\n")
cat("LANCET dataset dimensions:", dim(lancet_data), "\n")
cat("LIBRA2 dataset dimensions:", dim(libra2_data), "\n")
cat("Stata dataset (H_LASI_a3.dta) dimensions:", dim(stata_data), "\n")

# ==================== 4. ID Standardization and Variable Supplement ====================
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

# 4.3 Supplement missing variables from Stata dataset
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

# Function to extract total score column
extract_score_col <- function(data, score_type, data_name) {
  score_cols <- grep(paste0(score_type, "|score|total"), colnames(data), value = TRUE, ignore.case = TRUE)
  if (length(score_cols) == 0) score_cols <- setdiff(colnames(data), "ID")
  if (length(score_cols) == 0) stop(paste(data_name, "no score column found"))
  return(score_cols[1])
}

# Extract score columns for LANCET and LIBRA2
lancet_score_col <- extract_score_col(lancet_data, "LANCET", "LANCET dataset")
libra2_score_col <- extract_score_col(libra2_data, "LIBRA2", "LIBRA2 dataset")
cat("LANCET score column:", lancet_score_col, "\n")
cat("LIBRA2 score column:", libra2_score_col, "\n")

# Create subsets with standardized column names
lancet_subset <- lancet_data %>% select(ID, LANCET_total = !!sym(lancet_score_col))
libra2_subset <- libra2_data %>% select(ID, LIBRA2_total = !!sym(libra2_score_col))

# Merge all datasets and filter out missing outcome values
combined_data <- brain_age_data %>%
  left_join(lancet_subset, by = "ID") %>%
  left_join(libra2_subset, by = "ID") %>%
  filter(!is.na(!!sym(outcome_var)))
cat("Merged dataset dimensions:", dim(combined_data), "\n")

# ==================== 6. Covariate Missing Value Imputation ====================
cat("\nImputing missing values for covariates...\n")

# 6.1 Impute categorical variables with 9 (missing indicator)
categorical_vars <- c("raeduc_l", "ragender", "r1smoken", "r1hibpe", "r1rxhibp", "r1diabe", "r1rxdiab")
for (var in categorical_vars) {
  if (var %in% colnames(combined_data)) {
    missing_count <- sum(is.na(combined_data[[var]]))
    combined_data[[var]] <- replace(combined_data[[var]], is.na(combined_data[[var]]), 9)
    cat(sprintf("Categorical variable %s: %d missing values imputed with 9\n", var, missing_count))
  }
}

# 6.2 Impute continuous variables with median
continuous_vars <- c("r1systo", "r1diasto", "r1mbmi", "r1agey", "r1cesd10", "hh1itot", "r1drinkb")
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

# 7.1 Functions to calculate covariate scores
calculate_education_score <- function(raeduc_l) {
  case_when(
    raeduc_l %in% c(0,1,2,3,4) ~ 4,  # Highest education level
    raeduc_l == 5 ~ 3,
    raeduc_l %in% c(6,7) ~ 2,
    raeduc_l %in% c(8,9) ~ 1,        # Lowest education level
    TRUE ~ NA_real_
  )
}

calculate_hypertension_score <- function(r1hibpe, r1systo, r1diasto, r1rxhibp) {
  case_when(
    r1hibpe == 1 | r1systo > 140 | r1diasto > 90 | r1rxhibp == 1 ~ 1,  # Hypertensive
    TRUE ~ 4                                                           # Non-hypertensive
  )
}

calculate_diabetes_score <- function(r1diabe, r1rxdiab) {
  case_when(
    r1diabe == 0 ~ 4,                                                  # No diabetes
    r1diabe == 1 & r1rxdiab == 0 ~ 1,                                  # Diabetes without medication
    r1diabe == 1 & r1rxdiab == 1 ~ 2.5,                                # Diabetes with medication
    TRUE ~ NA_real_
  )
}

calculate_smoking_score <- function(r1smoken) {
  case_when(
    r1smoken == 0 ~ 4,  # Non-smoker
    r1smoken == 1 ~ 1,  # Smoker
    TRUE ~ NA_real_
  )
}

calculate_depression_score <- function(r1cesd10) {
  qs <- quantile(r1cesd10, c(0.25, 0.5, 0.75), na.rm = TRUE)
  case_when(
    r1cesd10 <= qs[1] ~ 4,  # Lowest depression risk
    r1cesd10 <= qs[2] ~ 3,
    r1cesd10 <= qs[3] ~ 2,
    r1cesd10 > qs[3] ~ 1,   # Highest depression risk
    TRUE ~ NA_real_
  )
}

# 7.2 Calculate scores and stratify exposure into quartiles
combined_data <- combined_data %>%
  mutate(
    # Covariate scores (impute missing scores with 9)
    education_score = replace(calculate_education_score(raeduc_l), is.na(calculate_education_score(raeduc_l)), 9),
    hypertension_score = replace(calculate_hypertension_score(r1hibpe, r1systo, r1diasto, r1rxhibp), 
                                 is.na(calculate_hypertension_score(r1hibpe, r1systo, r1diasto, r1rxhibp)), 9),
    diabetes_score = replace(calculate_diabetes_score(r1diabe, r1rxdiab), 
                             is.na(calculate_diabetes_score(r1diabe, r1rxdiab)), 9),
    smoking_score = replace(calculate_smoking_score(r1smoken), is.na(calculate_smoking_score(r1smoken)), 9),
    depression_score = if ("r1cesd10" %in% colnames(.)) {
      replace(calculate_depression_score(r1cesd10), is.na(calculate_depression_score(r1cesd10)), 9)
    } else 9,
    
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
    
    # Remove temporary variables
    exposure_q1 = NULL, exposure_q2 = NULL, exposure_q3 = NULL
  )

# Create stratified datasets (gender and age)
male_data <- combined_data %>% filter(ragender == 1)  # Male
female_data <- combined_data %>% filter(ragender == 2) # Female
age_lt65_data <- combined_data %>% filter(r1agey < 65) # <65 years
age_ge65_data <- combined_data %>% filter(r1agey >= 65) # ≥65 years

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

# Global statistics
global_stats <- calculate_stratified_stats(combined_data, "Global")
global_total <- calculate_total_stats(combined_data, "Global")

# Gender-stratified statistics
male_stats <- calculate_stratified_stats(male_data, "Male")
male_total <- calculate_total_stats(male_data, "Male")
female_stats <- calculate_stratified_stats(female_data, "Female")
female_total <- calculate_total_stats(female_data, "Female")

# Age-stratified statistics
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
    
    # 8.3 Tidy continuous exposure results
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
    
    # 8.4 Tidy quartile exposure results
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
    
    # 8.6 Combine results
    results_combined <- bind_rows(results_continuous, reference_row, results_quartile) %>%
      filter(exposure_category %in% c("Continuous", "Q1", "Q2", "Q3", "Q4"))
    
    list(results = results_combined, success = TRUE)
  }, error = function(e) {
    cat("Model", model_label, "failed:", e$message, "\n")
    list(success = FALSE)
  })
}

# ==================== 9. Fit Target Models ====================
cat("\nFitting target logistic regression models...\n")

# Define base covariate formula
base_covariates <- "education_score + hypertension_score + diabetes_score + smoking_score + depression_score"

# Initialize model list
model_list <- list()

# 9.1 Main model (Global + covariates)
model_list[["main_global"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~", base_covariates),
  model_label = "Main Model",
  stratify_label = "Global"
)

# 9.2 LANCET-adjusted model (Global + covariates + LANCET_total)
model_list[["lancet_adjusted"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~ LANCET_total +", base_covariates),
  model_label = "LANCET-adjusted",
  stratify_label = "Global"
)

# 9.3 LIBRA2-adjusted model (Global + covariates + LIBRA2_total)
model_list[["libra2_adjusted"]] <- run_logistic_models(
  data = combined_data,
  formula_str_base = paste(outcome_var, "~ LIBRA2_total +", base_covariates),
  model_label = "LIBRA2-adjusted",
  stratify_label = "Global"
)

# 9.4 Gender-stratified models (Male/Female + covariates)
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

# 9.5 Age-stratified models (<65y/≥65y + covariates)
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

# Filter successful models
successful_models <- model_list[sapply(model_list, function(x) x$success)]
if (length(successful_models) == 0) stop("No successful models")

# Combine model results
model_results <- map_dfr(successful_models, function(x) x$results)

# Add descriptive statistics (Case/N)
final_results <- model_results %>%
  left_join(
    descriptive_stats %>% select(strata = strata, exposure_category = exposure_quartile_f, Case_N),
    by = c("strata", "exposure_category")
  ) %>%
  # Reorder columns (Nature Genetics style)
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
output_path <- file.path('./output/', output_filename)
write_xlsx(final_results, output_path)
cat("✅ Results saved to:", output_path, "\n")

# Preview results (first 20 rows)
cat("\n=== Results Preview (First 20 Rows) ===\n")
print(head(final_results, 20), row.names = FALSE)

# ==================== 12. Generate Summary Statistics ====================
cat("\n=== Summary Statistics ===\n")
cat("Total sample size:", final_sample_size, "\n")
cat("Total cases:", sum(combined_data[[outcome_var]] == 1, na.rm = TRUE), "\n")
cat("Total controls:", sum(combined_data[[outcome_var]] == 0, na.rm = TRUE), "\n")
cat("Number of models successfully fitted:", length(successful_models), "\n")