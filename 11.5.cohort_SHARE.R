############################################################
# End of Code —— Logic and Methodology Verified
############################################################
library(survival)
library(dplyr)
library(readr)

##----------------------------------------------------------
## 1. Set Working Directory & Read Data
##----------------------------------------------------------
setwd("./input/BrainVutal8_dementia")

df <- read_csv("./input/final_cox_merged_processed2.csv")

##----------------------------------------------------------
## 2. Data Preprocessing
##----------------------------------------------------------

# Exclude baseline dementia
df <- df %>%
  filter(dementia_status != "baseline")

# Construct survival time and outcome (Convert days to years)
df <- df %>%
  mutate(
    time_days   = dementia_to_jinzu_days,   # Survival time (days)
    time_years  = dementia_to_jinzu_days / 365.25,  # Survival time (years)
    status = combined_diagnosed             # 1=Dementia, 0=Censored
  )
# Create BrainVital8 quartile groups (Q1 lowest, Q4 highest)
df <- df %>%
  mutate(
    BrainVital8_quartile = cut(
      BrainVital8_aligned,
      breaks = quantile(BrainVital8_aligned, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
      labels = c("Q1", "Q2", "Q3", "Q4"),
      include.lowest = TRUE
    )
  )

##----------------------------------------------------------
## 3. Revised Cox Analysis Function (Add group stats, include Q1 ref)
##----------------------------------------------------------
run_cox_analysis_with_groups <- function(data, formula_str, analysis_name, 
                                         required_vars = NULL, model_type = "both") {
  
  # Filter samples with no missing values if required variables are specified
  if (!is.null(required_vars)) {
    data_filtered <- data
    for (var in required_vars) {
      data_filtered <- data_filtered %>% 
        filter(!is.na(!!sym(var)))
    }
  } else {
    data_filtered <- data
  }
  
  # Calculate total sample size and number of cases
  total_n <- nrow(data_filtered)
  total_cases <- sum(data_filtered$status == 1, na.rm = TRUE)
  
  # Calculate sample size and cases for each BrainVital8 quartile
  group_stats <- data_filtered %>%
    group_by(BrainVital8_quartile) %>%
    summarise(
      Group_N = n(),
      Group_Case = sum(status == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(BrainVital8_quartile)
  
  cat("\nAnalysis:", analysis_name)
  cat("\n  Total N:", total_n)
  cat("\n  Total Cases:", total_cases)
  cat("\n  Sample counts per group:")
  for (i in 1:nrow(group_stats)) {
    cat(paste0("\n    ", group_stats$BrainVital8_quartile[i], ": N=", 
               group_stats$Group_N[i], ", Case=", group_stats$Group_Case[i]))
  }
  cat("\n")
  
  # Check for sufficient sample size
  if (total_n < 10 || total_cases < 5) {
    warning(paste("Analysis", analysis_name, "insufficient sample size"))
    return(NULL)
  }
  
  # Run Cox model
  model <- coxph(as.formula(formula_str), data = data_filtered)
  sm <- summary(model)
  
  # Extract all variable results
  results_full <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(sm$coefficients),
    HR       = exp(sm$coefficients[, "coef"]),
    CI_lower = exp(sm$coefficients[, "coef"] - 1.96 * sm$coefficients[, "se(coef)"]),
    CI_upper = exp(sm$coefficients[, "coef"] + 1.96 * sm$coefficients[, "se(coef)"]),
    P_value  = sm$coefficients[, "Pr(>|z|)"],
    N        = total_n,           
    Case     = total_cases,       
    stringsAsFactors = FALSE
  )
  
  # Extract BrainVital8 related variables
  brainvital8_vars <- c("BrainVital8_aligned", "BrainVital8_quartileQ2", 
                        "BrainVital8_quartileQ3", "BrainVital8_quartileQ4")
  
  # Create a data frame containing all BrainVital8 variables
  brainvital8_results <- data.frame()
  
  # 1. Add continuous variable results
  continuous_row <- results_full %>% 
    filter(grepl("BrainVital8_aligned", Variable))
  
  if (nrow(continuous_row) > 0) {
    continuous_row$Group_N <- total_n
    continuous_row$Group_Case <- total_cases
    brainvital8_results <- bind_rows(brainvital8_results, continuous_row)
  }
  
  # 2. Add categorical variable results
  
  # Find statistics for Q1 first
  q1_stats <- group_stats %>% filter(BrainVital8_quartile == "Q1")
  q1_n <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_N[1], NA)
  q1_case <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_Case[1], NA)
  
  # Add Q1 as reference group
  q1_row <- data.frame(
    Analysis = analysis_name,
    Variable = "BrainVital8_quartileQ1 (ref)",
    HR = 1.00,
    CI_lower = NA,
    CI_upper = NA,
    P_value = NA,
    N = total_n,
    Case = total_cases,
    Group_N = q1_n,
    Group_Case = q1_case,
    stringsAsFactors = FALSE
  )
  
  brainvital8_results <- bind_rows(brainvital8_results, q1_row)
  
  # Add Q2, Q3, Q4
  for (quartile in c("Q2", "Q3", "Q4")) {
    var_name <- paste0("BrainVital8_quartile", quartile)
    
    # Extract results for this variable from Cox results
    quartile_row <- results_full %>% filter(Variable == var_name)
    
    # Extract statistics for this group from group stats
    quartile_stats <- group_stats %>% filter(BrainVital8_quartile == quartile)
    quartile_n <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_N[1], NA)
    quartile_case <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_Case[1], NA)
    
    if (nrow(quartile_row) > 0) {
      quartile_row$Group_N <- quartile_n
      quartile_row$Group_Case <- quartile_case
      brainvital8_results <- bind_rows(brainvital8_results, quartile_row)
    } else {
      # If the variable is not in the Cox model, still add group info
      missing_row <- data.frame(
        Analysis = analysis_name,
        Variable = var_name,
        HR = NA,
        CI_lower = NA,
        CI_upper = NA,
        P_value = NA,
        N = total_n,
        Case = total_cases,
        Group_N = quartile_n,
        Group_Case = quartile_case,
        stringsAsFactors = FALSE
      )
      brainvital8_results <- bind_rows(brainvital8_results, missing_row)
    }
  }
  
  # Reorder: Continuous variable, followed by Q1-Q4
  brainvital8_results <- brainvital8_results %>%
    mutate(order = case_when(
      Variable == "BrainVital8_aligned" ~ 1,
      Variable == "BrainVital8_quartileQ1 (ref)" ~ 2,
      Variable == "BrainVital8_quartileQ2" ~ 3,
      Variable == "BrainVital8_quartileQ3" ~ 4,
      Variable == "BrainVital8_quartileQ4" ~ 5,
      TRUE ~ 6
    )) %>%
    arrange(order) %>%
    select(-order)
  
  # PH test
  ph_test <- cox.zph(model)
  
  return(list(
    results_full = results_full,
    results = brainvital8_results,
    ph_test = ph_test,
    model   = model,
    data_used = data_filtered,
    group_stats = group_stats,
    model_type = model_type
  ))
}

##----------------------------------------------------------
## 4. Function to run both Cox models
##----------------------------------------------------------
run_both_models <- function(data, analysis_name, required_vars) {
  
  # Model 1: Using continuous variable... (Modified covariates)
  formula_continuous <- "Surv(time_years, status) ~ BrainVital8_aligned + baseline_age + gender + bmi + drinking + hypertension"
  
  res_continuous <- run_cox_analysis_with_groups(
    data = data,
    formula_str = formula_continuous,
    analysis_name = paste0(analysis_name, "_continuous"),
    required_vars = required_vars,
    model_type = "continuous"
  )
  
  # Model 2: Using categorical variable... (Modified covariates)
  formula_categorical <- "Surv(time_years, status) ~ BrainVital8_quartile + baseline_age + gender + bmi + drinking + hypertension"
  
  res_categorical <- run_cox_analysis_with_groups(
    data = data,
    formula_str = formula_categorical,
    analysis_name = paste0(analysis_name, "_categorical"),
    required_vars = required_vars,
    model_type = "categorical"
  )
  
  return(list(
    continuous = res_continuous,
    categorical = res_categorical
  ))
}

##----------------------------------------------------------
## 5. Main Analysis (Full Population) - Both Models
##----------------------------------------------------------
cat("=== Starting sample screening for all analyses ===\n")

# Variables required for Main Analysis
main_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
               "baseline_age", "gender", "bmi", "drinking", "hypertension")

res_main_both <- run_both_models(
  data = df,
  analysis_name = "Main",
  required_vars = main_vars
)

res_main_continuous <- res_main_both$continuous
res_main_categorical <- res_main_both$categorical

##----------------------------------------------------------
## 6. Sex-stratified Analysis - Both Models
##----------------------------------------------------------

# Variables for Male analysis
male_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
               "baseline_age", "bmi", "drinking", "hypertension")

# Screen male samples first
df_male <- df %>% 
  filter(gender == 1) %>%
  filter(!is.na(gender))

res_male_both <- run_both_models(
  data = df_male,
  analysis_name = "Male",
  required_vars = male_vars
)

res_male_continuous <- res_male_both$continuous
res_male_categorical <- res_male_both$categorical

# Female analysis
df_female <- df %>% 
  filter(gender == 2) %>%
  filter(!is.na(gender))

res_female_both <- run_both_models(
  data = df_female,
  analysis_name = "Female",
  required_vars = male_vars
)

res_female_continuous <- res_female_both$continuous
res_female_categorical <- res_female_both$categorical

##----------------------------------------------------------
## 6.5 Age-stratified Analysis - Both Models
##----------------------------------------------------------

# Age stratification (<65 years)
df_age_lt65 <- df %>% 
  filter(baseline_age < 65) %>%
  filter(!is.na(baseline_age))

# baseline_age is no longer included in age-stratified analysis
age_strat_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
                    "gender", "bmi", "drinking", "hypertension")

# Age stratification (<65 years)
res_age_lt65_both <- run_both_models(
  data = df_age_lt65,
  analysis_name = "Age_lt65",
  required_vars = age_strat_vars
)

res_age_lt65_continuous <- res_age_lt65_both$continuous
res_age_lt65_categorical <- res_age_lt65_both$categorical

# Age stratification (≥65 years)
df_age_ge65 <- df %>% 
  filter(baseline_age >= 65) %>%
  filter(!is.na(baseline_age))

res_age_ge65_both <- run_both_models(
  data = df_age_ge65,
  analysis_name = "Age_ge65",
  required_vars = age_strat_vars
)

res_age_ge65_continuous <- res_age_ge65_both$continuous
res_age_ge65_categorical <- res_age_ge65_both$categorical

##----------------------------------------------------------
## 7. Exclude cases within first 2 years - Both Models
##----------------------------------------------------------

if (!is.null(res_main_continuous)) {
  df_ex2yr <- res_main_continuous$data_used %>%
    filter(!(status == 1 & time_years < 2))
  
  res_ex2yr_both <- run_both_models(
    data = df_ex2yr,
    analysis_name = "Exclude_2yr",
    required_vars = main_vars
  )
  
  res_ex2yr_continuous <- res_ex2yr_both$continuous
  res_ex2yr_categorical <- res_ex2yr_both$categorical
}

##----------------------------------------------------------
## 8. Analysis including LIBRA2 - Both Models
##----------------------------------------------------------
libra2_vars <- c(main_vars, "LIBRA2")

res_libra2_both <- run_both_models(
  data = df,
  analysis_name = "LIBRA2",
  required_vars = libra2_vars
)

res_libra2_continuous <- res_libra2_both$continuous
res_libra2_categorical <- res_libra2_both$categorical

##----------------------------------------------------------
## 9. Analysis including Lancet - Both Models
##----------------------------------------------------------
lancet_vars <- c(main_vars, "lancet")

res_lancet_both <- run_both_models(
  data = df,
  analysis_name = "Lancet",
  required_vars = lancet_vars
)

res_lancet_continuous <- res_lancet_both$continuous
res_lancet_categorical <- res_lancet_both$categorical

##----------------------------------------------------------
## 10. Merge and Save Results
##----------------------------------------------------------
# Collect all valid results
all_results <- list()

# Add continuous/categorical model results
if (!is.null(res_main_continuous)) all_results[["Main_continuous"]] <- res_main_continuous
if (!is.null(res_male_continuous)) all_results[["Male_continuous"]] <- res_male_continuous
if (!is.null(res_female_continuous)) all_results[["Female_continuous"]] <- res_female_continuous
if (!is.null(res_age_lt65_continuous)) all_results[["Age_lt65_continuous"]] <- res_age_lt65_continuous
if (!is.null(res_age_ge65_continuous)) all_results[["Age_ge65_continuous"]] <- res_age_ge65_continuous
if (!is.null(res_ex2yr_continuous)) all_results[["Exclude_2yr_continuous"]] <- res_ex2yr_continuous
if (!is.null(res_libra2_continuous)) all_results[["LIBRA2_continuous"]] <- res_libra2_continuous
if (!is.null(res_lancet_continuous)) all_results[["Lancet_continuous"]] <- res_lancet_continuous
if (!is.null(res_main_categorical)) all_results[["Main_categorical"]] <- res_main_categorical
if (!is.null(res_male_categorical)) all_results[["Male_categorical"]] <- res_male_categorical
if (!is.null(res_female_categorical)) all_results[["Female_categorical"]] <- res_female_categorical
if (!is.null(res_age_lt65_categorical)) all_results[["Age_lt65_categorical"]] <- res_age_lt65_categorical
if (!is.null(res_age_ge65_categorical)) all_results[["Age_ge65_categorical"]] <- res_age_ge65_categorical
if (!is.null(res_ex2yr_categorical)) all_results[["Exclude_2yr_categorical"]] <- res_ex2yr_categorical
if (!is.null(res_libra2_categorical)) all_results[["LIBRA2_categorical"]] <- res_libra2_categorical
if (!is.null(res_lancet_categorical)) all_results[["Lancet_categorical"]] <- res_lancet_categorical

# Merge BrainVital8 related results
brainvital8_results <- data.frame()

for (analysis_name in names(all_results)) {
  res <- all_results[[analysis_name]]
  brainvital8_results <- bind_rows(brainvital8_results, res$results)
}

# Merge full results (excluding group statistics)
full_results <- data.frame()

for (analysis_name in names(all_results)) {
  res <- all_results[[analysis_name]]
  full_results <- bind_rows(full_results, res$results_full)
}

##----------------------------------------------------------
## 11. Extract and Save Global PH Test Results
##----------------------------------------------------------
extract_global_ph_test_results <- function(results_list) {
  ph_results <- data.frame()
  
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    
    if (!is.null(res$ph_test)) {
      # Extract PH test results
      ph_table <- res$ph_test$table
      
      # Determine model type
      model_type <- ifelse(grepl("_continuous$", analysis_name), "continuous", "categorical")
      model_label <- ifelse(model_type == "continuous", 
                            "BrainVital8_aligned (continuous)", 
                            "BrainVital8_quartile (categorical)")
      
      # Extract Global PH test results
      if ("GLOBAL" %in% rownames(ph_table)) {
        global_ph <- ph_table["GLOBAL", ]
        global_ph_df <- data.frame(
          Analysis = gsub("_(continuous|categorical)$", "", analysis_name),
          Model_Type = model_type,
          Model_Label = model_label,
          Global_chisq = global_ph["chisq"],
          Global_df = global_ph["df"],
          Global_p = global_ph["p"],
          stringsAsFactors = FALSE
        )
        ph_results <- bind_rows(ph_results, global_ph_df)
      }
    }
  }
  
  return(ph_results)
}

# Extract global PH test results for all analyses
global_ph_results <- extract_global_ph_test_results(all_results)

##----------------------------------------------------------
## 12. Extract and Save Variable-specific PH Test Results
##----------------------------------------------------------
extract_variable_specific_ph_results <- function(results_list) {
  ph_results <- data.frame()
  
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    
    if (!is.null(res$ph_test)) {
      # Extract PH test results
      ph_table <- res$ph_test$table
      
      # Determine model type
      model_type <- ifelse(grepl("_continuous$", analysis_name), "continuous", "categorical")
      
      # Extract PH test results for BrainVital8 related variables
      brainvital8_vars <- grep("BrainVital8", rownames(ph_table), value = TRUE)
      
      if (length(brainvital8_vars) > 0) {
        for (var in brainvital8_vars) {
          var_ph <- ph_table[var, ]
          var_ph_df <- data.frame(
            Analysis = gsub("_(continuous|categorical)$", "", analysis_name),
            Model_Type = model_type,
            Variable = var,
            chisq = var_ph["chisq"],
            df = var_ph["df"],
            p = var_ph["p"],
            stringsAsFactors = FALSE
          )
          ph_results <- bind_rows(ph_results, var_ph_df)
        }
      }
    }
  }
  
  return(ph_results)
}

# Extract variable-specific PH test results
variable_ph_results <- extract_variable_specific_ph_results(all_results)

##----------------------------------------------------------
## 13. Save Results
##----------------------------------------------------------

# Save global PH test results
write.csv(
  global_ph_results,
  "./output/BrainVital8_Global_PH_test_results.csv",
  row.names = FALSE
)

# Save variable-specific PH test results
write.csv(
  variable_ph_results,
  "./output/BrainVital8_Variable_PH_test_results.csv",
  row.names = FALSE
)

# Save Cox analysis results
write.csv(
  brainvital8_results,
  "./output/BrainVital8_dementia_Cox_results_with_group_stats_incl_Q1xxxxx.csv",
  row.names = FALSE
)

##----------------------------------------------------------
## 14. Create Summary Table
##----------------------------------------------------------
create_ph_summary_table <- function(global_results, variable_results) {
  # Combine global results and variable results
  summary_table <- data.frame()
  
  for (analysis in unique(global_results$Analysis)) {
    # Get all results for this analysis
    analysis_global <- global_results %>% filter(Analysis == analysis)
    analysis_variable <- variable_results %>% filter(Analysis == analysis)
    
    # Continuous model
    cont_global <- analysis_global %>% filter(Model_Type == "continuous")
    cont_variable <- analysis_variable %>% filter(Model_Type == "continuous" & Variable == "BrainVital8_aligned")
    
    if (nrow(cont_global) > 0) {
      cont_row <- data.frame(
        Analysis = analysis,
        Model = "Continuous (BrainVital8_aligned)",
        BrainVital8_Variable_PH = ifelse(nrow(cont_variable) > 0, 
                                         paste0("χ²=", round(cont_variable$chisq, 3), 
                                                ", df=", cont_variable$df, 
                                                ", p=", round(cont_variable$p, 4)),
                                         "N/A"),
        Global_PH = paste0("χ²=", round(cont_global$Global_chisq, 3), 
                           ", df=", cont_global$Global_df, 
                           ", p=", round(cont_global$Global_p, 4)),
        stringsAsFactors = FALSE
      )
      summary_table <- bind_rows(summary_table, cont_row)
    }
    
    # Categorical model
    cat_global <- analysis_global %>% filter(Model_Type == "categorical")
    cat_variables <- analysis_variable %>% filter(Model_Type == "categorical" & grepl("BrainVital8_quartile", Variable))
    
    if (nrow(cat_global) > 0) {
      # Create variable PH result string for categorical model
      var_ph_strings <- c()
      if (nrow(cat_variables) > 0) {
        for (i in 1:nrow(cat_variables)) {
          var_name <- gsub("BrainVital8_quartile", "", cat_variables$Variable[i])
          var_ph_strings <- c(var_ph_strings, 
                              paste0(var_name, ": χ²=", round(cat_variables$chisq[i], 3),
                                     ", p=", round(cat_variables$p[i], 4)))
        }
        var_ph_combined <- paste(var_ph_strings, collapse = "; ")
      } else {
        var_ph_combined <- "N/A"
      }
      
      cat_row <- data.frame(
        Analysis = analysis,
        Model = "Categorical (BrainVital8_quartile)",
        BrainVital8_Variable_PH = var_ph_combined,
        Global_PH = paste0("χ²=", round(cat_global$Global_chisq, 3), 
                           ", df=", cat_global$Global_df, 
                           ", p=", round(cat_global$Global_p, 4)),
        stringsAsFactors = FALSE
      )
      summary_table <- bind_rows(summary_table, cat_row)
    }
  }
  
  return(summary_table)
}

# Create summary table
ph_summary_table <- create_ph_summary_table(global_ph_results, variable_ph_results)

# Save summary table
write.csv(
  ph_summary_table,
  "./output/BrainVital8_PH_test_summary_table.csv",
  row.names = FALSE
)
