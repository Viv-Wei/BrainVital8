aa_0826_filled <- read.csv("./input/charls_data.csv")

library(dplyr)
library(tidyverse)
library(survival)

# Enhanced main function: Run Cox analysis and extract results
run_cox_analysis <- function(data, formula_str, analysis_name, is_quartile = FALSE) {
  formula <- as.formula(formula_str)
  model <- coxph(formula, data = data)
  
  # Extract model results
  summary_model <- summary(model)
  
  # Get HR and 95% CI
  results <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(summary_model$coefficients),
    HR = exp(summary_model$coefficients[, "coef"]),
    CI_lower = exp(summary_model$coefficients[, "coef"] - 1.96 * summary_model$coefficients[, "se(coef)"]),
    CI_upper = exp(summary_model$coefficients[, "coef"] + 1.96 * summary_model$coefficients[, "se(coef)"]),
    P_value = summary_model$coefficients[, "Pr(>|z|)"],
    N = model$n,
    Case = sum(data$status == 1),
    stringsAsFactors = FALSE
  )
  
  # Add group-specific N and case counts for quartile analysis
  if(is_quartile && grepl("quartile", formula_str)) {
    # Extract quartile variable name
    quartile_var <- gsub(".*~ (BrainVital8_quartile[^ ]*).*", "\\1", formula_str)
    if(quartile_var == formula_str) {
      quartile_var <- "BrainVital8_quartile"
    }
    
    # Calculate N and case counts per quartile group
    if(quartile_var %in% names(data)) {
      quartile_counts <- data %>%
        group_by(!!sym(quartile_var)) %>%
        summarise(
          N_group = n(),
          Case_group = sum(status == 1),
          .groups = 'drop'
        )
      
      # Initialize group columns with NA
      results$Group_N <- NA
      results$Group_Case <- NA
      
      # Map group counts to results
      for(i in 1:nrow(quartile_counts)) {
        quartile_level <- as.character(quartile_counts[[quartile_var]][i])
        var_name <- paste0(quartile_var, quartile_level)
        
        # Find corresponding row in results
        row_idx <- grep(var_name, results$Variable)
        if(length(row_idx) > 0) {
          results$Group_N[row_idx] <- quartile_counts$N_group[i]
          results$Group_Case[row_idx] <- quartile_counts$Case_group[i]
        }
      }
      
      # Add reference group (first quartile) results
      ref_group <- levels(data[[quartile_var]])[1]
      ref_counts <- quartile_counts[quartile_counts[[quartile_var]] == ref_group, ]
      ref_row <- data.frame(
        Analysis = analysis_name,
        Variable = paste0(quartile_var, ref_group, " (ref)"),
        HR = 1.00,
        CI_lower = NA,
        CI_upper = NA,
        P_value = NA,
        N = model$n,
        Case = sum(data$status == 1),
        Group_N = ref_counts$N_group,
        Group_Case = ref_counts$Case_group,
        stringsAsFactors = FALSE
      )
      
      # Merge reference group into results
      results <- bind_rows(ref_row, results)
    }
  }
  
  # Perform proportional hazards (PH) assumption test
  ph_test <- cox.zph(model)
  
  return(list(results = results, ph_test = ph_test, model = model))
}

# Helper function: Extract quartile group counts
extract_quartile_counts <- function(data, quartile_var_name) {
  counts <- data %>%
    filter(!is.na(!!sym(quartile_var_name))) %>%
    group_by(!!sym(quartile_var_name)) %>%
    summarise(
      N = n(),
      Case = sum(status == 1),
      .groups = 'drop'
    )
  return(counts)
}

# Prepare main analysis dataset (remove rows with missing time/status)
aa_0826_filled <- aa_0826_filled[!is.na(aa_0826_filled$time) & !is.na(aa_0826_filled$status), ]

# Recalculate quartiles (based on full cohort)
aa_0826_filled$BrainVital8_quartile <- cut(aa_0826_filled$BrainVital8, 
                                           breaks = quantile(aa_0826_filled$BrainVital8, 
                                                             probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                             na.rm = TRUE),
                                           labels = c("Q1", "Q2", "Q3", "Q4"),
                                           include.lowest = TRUE)

# Store all analysis results
all_results <- list()

## 1. Main analysis (continuous variable)
res_cont_main <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Main_Continuous"
)

## 2. Main analysis (quartile) - extract group counts first
quart_counts_main <- extract_quartile_counts(aa_0826_filled, "BrainVital8_quartile")
print("Main analysis quartile group counts:")
print(quart_counts_main)

res_quart_main <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Main_Quartile",
  is_quartile = TRUE
)

## 3. Stratified analysis (male/female separately)
# Male subgroup
data_male <- aa_0826_filled[aa_0826_filled$gender == 1, ]
data_male$BrainVital8_quartile_male <- cut(data_male$BrainVital8,
                                           breaks = quantile(data_male$BrainVital8, 
                                                             probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                             na.rm = TRUE),
                                           labels = c("Q1", "Q2", "Q3", "Q4"),
                                           include.lowest = TRUE)

# Extract male quartile group counts
quart_counts_male <- extract_quartile_counts(data_male, "BrainVital8_quartile_male")
print("Male subgroup quartile group counts:")
print(quart_counts_male)

# Female subgroup
data_female <- aa_0826_filled[aa_0826_filled$gender == 0, ]
data_female$BrainVital8_quartile_female <- cut(data_female$BrainVital8,
                                               breaks = quantile(data_female$BrainVital8, 
                                                                 probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                                 na.rm = TRUE),
                                               labels = c("Q1", "Q2", "Q3", "Q4"),
                                               include.lowest = TRUE)

# Extract female quartile group counts
quart_counts_female <- extract_quartile_counts(data_female, "BrainVital8_quartile_female")
print("Female subgroup quartile group counts:")
print(quart_counts_female)

# Run stratified analyses (keep original code, only modify function calls)
res_male_cont <- run_cox_analysis(
  data = data_male,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Male_Continuous"
)

res_male_quart <- run_cox_analysis(
  data = data_male,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_male + age + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Male_Quartile",
  is_quartile = TRUE
)

res_female_cont <- run_cox_analysis(
  data = data_female,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Female_Continuous"
)

res_female_quart <- run_cox_analysis(
  data = data_female,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_female + age + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Female_Quartile",
  is_quartile = TRUE
)

## 4. Sensitivity analysis: Exclude events within 2 years
data_exclude_2year <- aa_0826_filled[aa_0826_filled$time > 2 | (aa_0826_filled$time <= 2 & aa_0826_filled$status == 0), ]
data_exclude_2year$BrainVital8_quartile_ex2 <- cut(data_exclude_2year$BrainVital8,
                                                   breaks = quantile(data_exclude_2year$BrainVital8, 
                                                                     probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                                     na.rm = TRUE),
                                                   labels = c("Q1", "Q2", "Q3", "Q4"),
                                                   include.lowest = TRUE)

# Extract quartile counts for 2-year exclusion analysis
quart_counts_ex2 <- extract_quartile_counts(data_exclude_2year, "BrainVital8_quartile_ex2")
print("2-year exclusion quartile group counts:")
print(quart_counts_ex2)

res_ex2_cont <- run_cox_analysis(
  data = data_exclude_2year,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Exclude2yr_Continuous"
)

res_ex2_quart <- run_cox_analysis(
  data = data_exclude_2year,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_ex2 + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Exclude2yr_Quartile",
  is_quartile = TRUE
)

## 5. Adjusted analysis: Add lancet1 covariate
res_lancet_cont <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression + lancet1",
  analysis_name = "PlusLancet_Continuous"
)

res_lancet_quart <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression + lancet1",
  analysis_name = "PlusLancet_Quartile",
  is_quartile = TRUE
)

## 6. Adjusted analysis: Add LIBRA21 covariate
res_libra_cont <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression + LIBRA21",
  analysis_name = "PlusLIBRA_Continuous"
)

res_libra_quart <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression + LIBRA21",
  analysis_name = "PlusLIBRA_Quartile",
  is_quartile = TRUE
)

# Function to standardize result dataframe structure
standardize_results <- function(df) {
  # Check and add missing columns
  required_cols <- c("Analysis", "Variable", "HR", "CI_lower", "CI_upper", 
                     "P_value", "N", "Case", "Group_N", "Group_Case")
  
  for(col in required_cols) {
    if(!col %in% names(df)) {
      df[[col]] <- NA
    }
  }
  
  # Ensure consistent column order
  df <- df[, required_cols]
  return(df)
}

# Standardize all results
results_list <- list(
  res_cont_main$results,
  res_quart_main$results,
  res_male_cont$results,
  res_male_quart$results,
  res_female_cont$results,
  res_female_quart$results,
  res_ex2_cont$results,
  res_ex2_quart$results,
  res_lancet_cont$results,
  res_lancet_quart$results,
  res_libra_cont$results,
  res_libra_quart$results
)

# Standardize and merge all results
results_standardized <- lapply(results_list, standardize_results)
all_results_df <- do.call(rbind, results_standardized)

# Print merged results
print(all_results_df)

library(dplyr)
library(survival)

# 1. Data stratification: Filter populations <65 years and ≥65 years (ensure age column exists and is numeric)
data_age_under65 <- aa_0826_filled[aa_0826_filled$age < 65 & !is.na(aa_0826_filled$age), ]
data_age_over65 <- aa_0826_filled[aa_0826_filled$age >= 65 & !is.na(aa_0826_filled$age), ]

# 2. Recalculate quartiles for stratified populations (quartiles within each stratum)
# Quartiles for <65 years population
data_age_under65$BrainVital8_quartile_age65 <- cut(
  data_age_under65$BrainVital8, 
  breaks = quantile(data_age_under65$BrainVital8, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)

# Quartiles for ≥65 years population
data_age_over65$BrainVital8_quartile_age65 <- cut(
  data_age_over65$BrainVital8, 
  breaks = quantile(data_age_over65$BrainVital8, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)

# 3. Extract stratified quartile group information (optional: print for verification)
quart_counts_under65 <- extract_quartile_counts(data_age_under65, "BrainVital8_quartile_age65")
quart_counts_over65 <- extract_quartile_counts(data_age_over65, "BrainVital8_quartile_age65")
cat("==== Quartile group information for age <65 years ====\n")
print(quart_counts_under65)
cat("==== Quartile group information for age ≥65 years ====\n")
print(quart_counts_over65)

# 4. Run age-stratified Cox analysis (covariates: ragender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression)
## 4.1 <65 years - continuous independent variable
res_age_under65_cont <- run_cox_analysis(
  data = data_age_under65,
  formula_str = "Surv(time, status) ~ BrainVital8 + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Age_Under65_Continuous"
)

## 4.2 <65 years - quartile independent variable
res_age_under65_quart <- run_cox_analysis(
  data = data_age_under65,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_age65 + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Age_Under65_Quartile",
  is_quartile = TRUE
)

## 4.3 ≥65 years - continuous independent variable
res_age_over65_cont <- run_cox_analysis(
  data = data_age_over65,
  formula_str = "Surv(time, status) ~ BrainVital8 + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Age_Over65_Continuous"
)

## 4.4 ≥65 years - quartile independent variable
res_age_over65_quart <- run_cox_analysis(
  data = data_age_over65,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_age65 + gender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression",
  analysis_name = "Age_Over65_Quartile",
  is_quartile = TRUE
)

# 5. Merge age-stratified Cox results into the master dataframe
## 5.1 Standardize stratified results
age_results_list <- list(
  res_age_under65_cont$results,
  res_age_under65_quart$results,
  res_age_over65_cont$results,
  res_age_over65_quart$results
)
age_results_standardized <- lapply(age_results_list, standardize_results)
age_all_results_df <- do.call(rbind, age_results_standardized)

## 5.2 Merge into the original master results (optional: add age-stratified results to full results)
all_results_df_with_age <- rbind(all_results_df, age_all_results_df)