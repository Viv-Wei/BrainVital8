library(dplyr)
library(tidyverse)

aa_0826_filled <- read.csv("./input/charls.csv")

aa_0826_filled$sex <- as.factor(aa_0826_filled$sex)
aa_0826_filled$drink <- as.factor(aa_0826_filled$drink)
aa_0826_filled$hypertension <- as.factor(aa_0826_filled$hypertension)
table(aa_0826_filled$sex,useNA = "always")
table(aa_0826_filled$drink,useNA = "always")
table(aa_0826_filled$hypertension,useNA = "always")

##########################Sensitivity Analysis################################
library(survival)

# Enhanced main function: Run Cox analysis and extract results
run_cox_analysis <- function(data, formula_str, analysis_name, is_quartile = FALSE) {
  formula <- as.formula(formula_str)
  model <- coxph(formula, data = data)
  
  # Extract results
  summary_model <- summary(model)
  
  # Get HR and CI
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
  
  # If it is quartile analysis, add case and N for each group
  if(is_quartile && grepl("quartile", formula_str)) {
    # Extract quartile variable name
    quartile_var <- gsub(".*~ (BrainVital8_quartile[^ ]*).*", "\\1", formula_str)
    if(quartile_var == formula_str) {
      quartile_var <- "BrainVital8_quartile"
    }
    
    # Calculate case and N for each group
    if(quartile_var %in% names(data)) {
      quartile_counts <- data %>%
        group_by(!!sym(quartile_var)) %>%
        summarise(
          N_group = n(),
          Case_group = sum(status == 1),
          .groups = 'drop'
        )
      
      # Add group information to results
      results$Group_N <- NA
      results$Group_Case <- NA
      
      for(i in 1:nrow(quartile_counts)) {
        quartile_level <- as.character(quartile_counts[[quartile_var]][i])
        var_name <- paste0(quartile_var, quartile_level)
        
        # Find the corresponding row
        row_idx <- grep(var_name, results$Variable)
        if(length(row_idx) > 0) {
          results$Group_N[row_idx] <- quartile_counts$N_group[i]
          results$Group_Case[row_idx] <- quartile_counts$Case_group[i]
        }
      }
      
      # Add case and N for reference group
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
  
  # PH test
  ph_test <- cox.zph(model)
  
  return(list(results = results, ph_test = ph_test, model = model))
}

# Another function: Extract quartile group information specifically
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

# Prepare data for main analysis
aa_0826_filled <- aa_0826_filled[!is.na(aa_0826_filled$time) & !is.na(aa_0826_filled$status), ]

# Recalculate quartiles (based on the full population)
aa_0826_filled$BrainVital8_quartile <- cut(aa_0826_filled$BrainVital8, 
                                           breaks = quantile(aa_0826_filled$BrainVital8, 
                                                             probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                             na.rm = TRUE),
                                           labels = c("Q1", "Q2", "Q3", "Q4"),
                                           include.lowest = TRUE)

# Store all results
all_results <- list()

## 1. Main analysis (continuous)
res_cont_main <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + sex + bmi + drink + hypertension",
  analysis_name = "Main_Continuous"
)

## 2. Main analysis (quartile) - Extract group information first
quart_counts_main <- extract_quartile_counts(aa_0826_filled, "BrainVital8_quartile")
print("Main analysis quartile group information:")
print(quart_counts_main)

res_quart_main <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + drink + hypertension",
  analysis_name = "Main_Quartile",
  is_quartile = TRUE
)

## 3. Stratified analysis (separated by gender)
# Male
data_male <- aa_0826_filled[aa_0826_filled$sex == 1, ]
data_male$BrainVital8_quartile_male <- cut(data_male$BrainVital8,
                                           breaks = quantile(data_male$BrainVital8, 
                                                             probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                             na.rm = TRUE),
                                           labels = c("Q1", "Q2", "Q3", "Q4"),
                                           include.lowest = TRUE)

# Extract quartile information for males
quart_counts_male <- extract_quartile_counts(data_male, "BrainVital8_quartile_male")
print("Male quartile group information:")
print(quart_counts_male)

# Female
data_female <- aa_0826_filled[aa_0826_filled$sex == 0, ]
data_female$BrainVital8_quartile_female <- cut(data_female$BrainVital8,
                                               breaks = quantile(data_female$BrainVital8, 
                                                                 probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                                 na.rm = TRUE),
                                               labels = c("Q1", "Q2", "Q3", "Q4"),
                                               include.lowest = TRUE)

# Extract quartile information for females
quart_counts_female <- extract_quartile_counts(data_female, "BrainVital8_quartile_female")
print("Female quartile group information:")
print(quart_counts_female)

# Run each analysis (keep your original code, only modify function calls)
res_male_cont <- run_cox_analysis(
  data = data_male,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + bmi + drink + hypertension",
  analysis_name = "Male_Continuous"
)

res_male_quart <- run_cox_analysis(
  data = data_male,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_male + age + bmi + drink + hypertension",
  analysis_name = "Male_Quartile",
  is_quartile = TRUE
)

res_female_cont <- run_cox_analysis(
  data = data_female,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + bmi + drink + hypertension",
  analysis_name = "Female_Continuous"
)

res_female_quart <- run_cox_analysis(
  data = data_female,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_female + age  + bmi + drink + hypertension",
  analysis_name = "Female_Quartile",
  is_quartile = TRUE
)

## 4. Exclude incident cases within 2 years
data_exclude_2year <- aa_0826_filled[aa_0826_filled$time > 2 | (aa_0826_filled$time <= 2 & aa_0826_filled$status == 0), ]
data_exclude_2year$BrainVital8_quartile_ex2 <- cut(data_exclude_2year$BrainVital8,
                                                   breaks = quantile(data_exclude_2year$BrainVital8, 
                                                                     probs = c(0, 0.25, 0.50, 0.75, 1), 
                                                                     na.rm = TRUE),
                                                   labels = c("Q1", "Q2", "Q3", "Q4"),
                                                   include.lowest = TRUE)

# Extract quartile information after excluding 2-year incident cases
quart_counts_ex2 <- extract_quartile_counts(data_exclude_2year, "BrainVital8_quartile_ex2")
print("Quartile group information after excluding 2-year incident cases:")
print(quart_counts_ex2)

res_ex2_cont <- run_cox_analysis(
  data = data_exclude_2year,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + sex + bmi + drink + hypertension",
  analysis_name = "Exclude2yr_Continuous"
)

res_ex2_quart <- run_cox_analysis(
  data = data_exclude_2year,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_ex2 + age + sex + bmi + drink + hypertension",
  analysis_name = "Exclude2yr_Quartile",
  is_quartile = TRUE
)

## 5. Add lancet1
res_lancet_cont <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + sex + bmi + drink + hypertension + LANCET",
  analysis_name = "PlusLancet_Continuous"
)

res_lancet_quart <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + drink + hypertension + LANCET",
  analysis_name = "PlusLancet_Quartile",
  is_quartile = TRUE
)

## 6. Add LIBRA21
res_libra_cont <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8 + age + sex + bmi + drink + hypertension + LIBRA2",
  analysis_name = "PlusLIBRA_Continuous"
)

res_libra_quart <- run_cox_analysis(
  data = aa_0826_filled,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + drink + hypertension + LIBRA2",
  analysis_name = "PlusLIBRA_Quartile",
  is_quartile = TRUE
)

# Create a function to unify the data frame structure
standardize_results <- function(df) {
  # Check and add missing columns
  required_cols <- c("Analysis", "Variable", "HR", "CI_lower", "CI_upper", 
                     "P_value", "N", "Case", "Group_N", "Group_Case")
  
  for(col in required_cols) {
    if(!col %in% names(df)) {
      df[[col]] <- NA
    }
  }
  
  # Ensure the order of columns is consistent
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

# Standardize and merge
results_standardized <- lapply(results_list, standardize_results)
all_results_df <- do.call(rbind, results_standardized)

# Output merged results
print(all_results_df)


# Load required packages
library(dplyr)

# Create function: Extract and organize PH test results
extract_ph_results <- function(ph_test_obj, analysis_name) {
  # Extract core results of cox.zph
  ph_results <- as.data.frame(ph_test_obj$table) %>%
    tibble::rownames_to_column("Variable") %>%
    rename(
      chisq = `chisq`,
      df = `df`,
      p_value = `p`
    ) %>%
    mutate(
      Analysis = analysis_name,  # Add analysis name
      # Supplement column names for easy merging
      chisq = round(chisq, 4),
      df = as.integer(df),
      p_value = round(p_value, 4),
      # Add judgment of PH test results (violate PH assumption if p<0.05)
      ph_violation = ifelse(p_value < 0.05, "Yes", "No")
    ) %>%
    select(Analysis, Variable, chisq, df, p_value, ph_violation)
  
  # Supplement PH test results of the overall model (global test)
  global_chisq <- sum(ph_test_obj$table[, "chisq"])
  global_df <- sum(ph_test_obj$table[, "df"])
  global_p <- 1 - pchisq(global_chisq, global_df)
  
  global_row <- data.frame(
    Analysis = analysis_name,
    Variable = "Global (Overall)",
    chisq = round(global_chisq, 4),
    df = as.integer(global_df),
    p_value = round(global_p, 4),
    ph_violation = ifelse(global_p < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  # Merge variable-level and overall PH test results
  ph_results <- bind_rows(ph_results, global_row)
  
  return(ph_results)
}

# Extract PH test results of all analyses
ph_results_list <- list(
  extract_ph_results(res_cont_main$ph_test, "Main_Continuous"),
  extract_ph_results(res_quart_main$ph_test, "Main_Quartile"),
  extract_ph_results(res_male_cont$ph_test, "Male_Continuous"),
  extract_ph_results(res_male_quart$ph_test, "Male_Quartile"),
  extract_ph_results(res_female_cont$ph_test, "Female_Continuous"),
  extract_ph_results(res_female_quart$ph_test, "Female_Quartile"),
  extract_ph_results(res_ex2_cont$ph_test, "Exclude2yr_Continuous"),
  extract_ph_results(res_ex2_quart$ph_test, "Exclude2yr_Quartile"),
  extract_ph_results(res_lancet_cont$ph_test, "PlusLancet_Continuous"),
  extract_ph_results(res_lancet_quart$ph_test, "PlusLancet_Quartile"),
  extract_ph_results(res_libra_cont$ph_test, "PlusLIBRA_Continuous"),
  extract_ph_results(res_libra_quart$ph_test, "PlusLIBRA_Quartile")
)

# Merge all PH test results into one data frame
all_ph_results <- do.call(rbind, ph_results_list)





####################Supplement age stratification####################
# ===================== New: Age-stratified analysis (<65 years old / ≥65 years old) =====================
library(dplyr)
library(survival)

# 1. Data stratification: Filter populations under 65 and over 65 (ensure age column exists and is numeric)
data_age_under65 <- aa_0826_filled[aa_0826_filled$age < 65 & !is.na(aa_0826_filled$age), ]
data_age_over65 <- aa_0826_filled[aa_0826_filled$age >= 65 & !is.na(aa_0826_filled$age), ]

# 2. Recalculate quartiles for stratified populations (quartiles within each stratum)
# Quartiles for population under 65
data_age_under65$BrainVital8_quartile_age65 <- cut(
  data_age_under65$BrainVital8, 
  breaks = quantile(data_age_under65$BrainVital8, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)

# Quartiles for population over 65
data_age_over65$BrainVital8_quartile_age65 <- cut(
  data_age_over65$BrainVital8, 
  breaks = quantile(data_age_over65$BrainVital8, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE),
  labels = c("Q1", "Q2", "Q3", "Q4"),
  include.lowest = TRUE
)

# 3. Extract stratified quartile group information (optional: print for verification)
quart_counts_under65 <- extract_quartile_counts(data_age_under65, "BrainVital8_quartile_age65")
quart_counts_over65 <- extract_quartile_counts(data_age_over65, "BrainVital8_quartile_age65")
cat("==== Quartile group information for age <65 ====\n")
print(quart_counts_under65)
cat("==== Quartile group information for age ≥65 ====\n")
print(quart_counts_over65)

# 4. Run age-stratified Cox analysis (covariates: ragender + bmi + smoke + drink + raeduc_c + income_total + hyper + diabetes + depression)
## 4.1 Under 65 - Continuous independent variable
res_age_under65_cont <- run_cox_analysis(
  data = data_age_under65,
  formula_str = "Surv(time, status) ~ BrainVital8 + sex + bmi + drink + hypertension",
  analysis_name = "Age_Under65_Continuous"
)

## 4.2 Under 65 - Quartile independent variable
res_age_under65_quart <- run_cox_analysis(
  data = data_age_under65,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_age65 + sex + bmi + drink + hypertension",
  analysis_name = "Age_Under65_Quartile",
  is_quartile = TRUE
)

## 4.3 Over 65 - Continuous independent variable
res_age_over65_cont <- run_cox_analysis(
  data = data_age_over65,
  formula_str = "Surv(time, status) ~ BrainVital8 + sex + bmi + drink + hypertension",
  analysis_name = "Age_Over65_Continuous"
)

## 4.4 Over 65 - Quartile independent variable
res_age_over65_quart <- run_cox_analysis(
  data = data_age_over65,
  formula_str = "Surv(time, status) ~ BrainVital8_quartile_age65 + sex + bmi + drink + hypertension",
  analysis_name = "Age_Over65_Quartile",
  is_quartile = TRUE
)

# 5. Merge age-stratified Cox results into the total data frame
## 5.1 Standardize stratified results
age_results_list <- list(
  res_age_under65_cont$results,
  res_age_under65_quart$results,
  res_age_over65_cont$results,
  res_age_over65_quart$results
)
age_results_standardized <- lapply(age_results_list, standardize_results)
age_all_results_df <- do.call(rbind, age_results_standardized)

## 5.2 Merge into the original total results (optional: if need to add age-stratified results to full results)
all_results_df_with_age <- rbind(all_results_df, age_all_results_df)

# 6. Extract PH test results for age stratification (reuse the previous PH test extraction function)
## 6.1 Load PH test extraction function (run this section first if not defined before)
if (!exists("extract_ph_results")) {
  extract_ph_results <- function(ph_test_obj, analysis_name) {
    ph_results <- as.data.frame(ph_test_obj$table) %>%
      tibble::rownames_to_column("Variable") %>%
      rename(
        chisq = `chisq`,
        df = `df`,
        p_value = `p`
      ) %>%
      mutate(
        Analysis = analysis_name,
        chisq = round(chisq, 4),
        df = as.integer(df),
        p_value = round(p_value, 4),
        ph_violation = ifelse(p_value < 0.05, "Yes", "No")
      ) %>%
      select(Analysis, Variable, chisq, df, p_value, ph_violation)
    
    # Supplement overall PH test
    global_chisq <- sum(ph_test_obj$table[, "chisq"])
    global_df <- sum(ph_test_obj$table[, "df"])
    global_p <- 1 - pchisq(global_chisq, global_df)
    
    global_row <- data.frame(
      Analysis = analysis_name,
      Variable = "Global (Overall)",
      chisq = round(global_chisq, 4),
      df = as.integer(global_df),
      p_value = round(global_p, 4),
      ph_violation = ifelse(global_p < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    )
    
    ph_results <- bind_rows(ph_results, global_row)
    return(ph_results)
  }
}

## 6.2 Extract age-stratified PH test results
age_ph_results_list <- list(
  extract_ph_results(res_age_under65_cont$ph_test, "Age_Under65_Continuous"),
  extract_ph_results(res_age_under65_quart$ph_test, "Age_Under65_Quartile"),
  extract_ph_results(res_age_over65_cont$ph_test, "Age_Over65_Continuous"),
  extract_ph_results(res_age_over65_quart$ph_test, "Age_Over65_Quartile")
)
age_all_ph_results <- do.call(rbind, age_ph_results_list)

all_ph_results <- rbind(all_ph_results, age_all_ph_results)
