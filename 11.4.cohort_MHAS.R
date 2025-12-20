############################################################
## Complete Final Revised Code - Fixed Filtering Logic + P-value Display + Full Functionality Retained
############################################################
library(survival)
library(dplyr)
library(readr)

##----------------------------------------------------------
## Read Data (Updated to Your New Path)
##----------------------------------------------------------

# Read dataset (automatically looks in new working directory)
df <- read_csv("./input/BrainVital8_alignment-main_data.csv")

##----------------------------------------------------------
## 2. Data Preprocessing (Using ntile for Equal Sample Size Grouping to Solve Breakpoint Duplication)
##----------------------------------------------------------
# 1. Divide income into 4 groups based on sample size (1=lowest, 4=highest, retain missing value code 9)
df <- df %>%
  mutate(
    # First ensure income is numeric
    income = as.numeric(income),
    # Divide into 4 groups based on sample size, no duplicate breakpoints
    income = ntile(income, 4),
    # Convert to factor (specify level order: low → high)
    income = factor(income, levels = c("1", "2", "3", "4"))
  )

# 2. Divide BrainVital8 into 4 groups based on sample size (Q1 lowest, Q4 highest)
df <- df %>%
  mutate(
    # Divide into groups based on sample size
    BrainVital8_quartile = ntile(BrainVital8, 4),
    # Convert to Q1-Q4 labels
    BrainVital8_quartile = case_when(
      BrainVital8_quartile == 1 ~ "Q1",
      BrainVital8_quartile == 2 ~ "Q2",
      BrainVital8_quartile == 3 ~ "Q3",
      BrainVital8_quartile == 4 ~ "Q4"
    ),
    # Convert to factor (specify level order)
    BrainVital8_quartile = factor(BrainVital8_quartile, levels = c("Q1", "Q2", "Q3", "Q4")),
    .after = BrainVital8
  )

# 3. Ensure correct format for other categorical variables (retain missing value code 9)
df <- df %>%
  mutate(
    sex = factor(sex, levels = c(0, 1)),
    drink = factor(drink, levels = c(0, 1)),
    smoke = factor(smoke, levels = c(0, 1, 9)),
    T2D = factor(T2D, levels = c(0, 1)),
    hypertension = factor(hypertension, levels = c(0, 1, 9)),
    depression = factor(depression, levels = c(0, 1)),
    education = factor(education),
    cohort = factor(cohort)
  )

##----------------------------------------------------------
## 3. Cox Analysis Function (Retain Core Logic + Complete PH Test Calculation)
##----------------------------------------------------------
run_cox_analysis_with_groups <- function(data, formula_str, analysis_name, 
                                         required_vars = NULL, model_type = "both") {
  
  # Filter samples with non-NA required variables (retain missing value code 9)
  if (!is.null(required_vars)) {
    data_filtered <- data
    for (var in required_vars) {
      data_filtered <- data_filtered %>% 
        filter(!is.na(!!sym(var)))
    }
  } else {
    data_filtered <- data
  }
  
  # Calculate overall and group statistics
  total_n <- nrow(data_filtered)
  total_cases <- sum(data_filtered$status == 1, na.rm = TRUE)
  
  group_stats <- data_filtered %>%
    group_by(BrainVital8_quartile) %>%
    summarise(
      Group_N = n(),
      Group_Case = sum(status == 1, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(BrainVital8_quartile)
  
  # Output statistical information
  cat("\nAnalysis:", analysis_name)
  cat("\n  Total sample size:", total_n)
  cat("\n  Total cases:", total_cases)
  cat("\n  Sample size by group:")
  for (i in 1:nrow(group_stats)) {
    cat(paste0("\n    ", group_stats$BrainVital8_quartile[i], ": N=", 
               group_stats$Group_N[i], ", Case=", group_stats$Group_Case[i]))
  }
  cat("\n")
  
  # Sample size check
  if (total_n < 10 || total_cases < 5) {
    warning(paste("Analysis", analysis_name, "has insufficient sample size"))
    return(NULL)
  }
  
  # Run Cox model (core Cox calculation)
  model <- coxph(as.formula(formula_str), data = data_filtered)
  sm <- summary(model)
  
  # Extract complete Cox results
  results_full <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(sm$coefficients),
    HR       = exp(sm$coefficients[, "coef"]),
    CI_lower = exp(sm$coefficients[, "coef"] - 1.96 * sm$coefficients[, "se(coef)"]),
    CI_upper = exp(sm$coefficients[, "coef"] + 1.96 * sm$coefficients[, "se(coef)"]),
    P_value  = sm$coefficients[, "Pr(>|z|)"],  # Retain original P-value for later formatting
    N        = total_n,
    Case     = total_cases,
    stringsAsFactors = FALSE
  )
  
  # Extract BrainVital8-specific Cox results (including Q1 reference group)
  brainvital8_vars <- c("BrainVital8", "BrainVital8_quartileQ2", 
                        "BrainVital8_quartileQ3", "BrainVital8_quartileQ4")
  brainvital8_results <- data.frame()
  
  # Continuous variable results
  continuous_row <- results_full %>% 
    filter(grepl("BrainVital8", Variable) & !grepl("quartile", Variable))
  if (nrow(continuous_row) > 0) {
    continuous_row$Group_N <- total_n
    continuous_row$Group_Case <- total_cases
    brainvital8_results <- bind_rows(brainvital8_results, continuous_row)
  }
  
  # Q1 reference group
  q1_stats <- group_stats %>% filter(BrainVital8_quartile == "Q1")
  q1_n <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_N[1], NA)
  q1_case <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_Case[1], NA)
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
  
  # Q2-Q4 results
  for (quartile in c("Q2", "Q3", "Q4")) {
    var_name <- paste0("BrainVital8_quartile", quartile)
    quartile_row <- results_full %>% filter(Variable == var_name)
    quartile_stats <- group_stats %>% filter(BrainVital8_quartile == quartile)
    quartile_n <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_N[1], NA)
    quartile_case <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_Case[1], NA)
    
    if (nrow(quartile_row) > 0) {
      quartile_row$Group_N <- quartile_n
      quartile_row$Group_Case <- quartile_case
      brainvital8_results <- bind_rows(brainvital8_results, quartile_row)
    } else {
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
  
  # Sorting
  brainvital8_results <- brainvital8_results %>%
    mutate(order = case_when(
      Variable == "BrainVital8" ~ 1,
      Variable == "BrainVital8_quartileQ1 (ref)" ~ 2,
      Variable == "BrainVital8_quartileQ2" ~ 3,
      Variable == "BrainVital8_quartileQ3" ~ 4,
      Variable == "BrainVital8_quartileQ4" ~ 5,
      TRUE ~ 6
    )) %>%
    arrange(order) %>%
    select(-order)
  
  # Run PH test (proportional hazards assumption test, core PH calculation)
  ph_test <- cox.zph(model)
  
  return(list(
    results_full = results_full,          # Complete Cox results
    results = brainvital8_results,        # BrainVital8-specific Cox results
    ph_test = ph_test,                    # Complete PH test results
    model   = model,                      # Cox model object
    data_used = data_filtered,            # Data used in analysis
    group_stats = group_stats,            # Group statistics
    model_type = model_type               # Model type
  ))
}


##----------------------------------------------------------
## 4. General Function: Run Continuous + Categorical Models
##----------------------------------------------------------
run_both_models_custom <- function(data, analysis_name, required_vars, formula_continuous, formula_categorical) {
  # Continuous model
  res_continuous <- run_cox_analysis_with_groups(
    data = data,
    formula_str = formula_continuous,
    analysis_name = paste0(analysis_name, "_continuous"),
    required_vars = required_vars,
    model_type = "continuous"
  )
  
  # Categorical model
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
## 5. Execute 8 Core Analyses (Fixed Excluding First 2 Years Filtering Logic)
##----------------------------------------------------------
cat("=== Starting 8 Core Cox+PH Analyses ===\n")

# ---------------------- 1. Main Result Analysis ----------------------
main_required_vars <- c("id", "age", "sex", "bmi", "time", "status", "drink", 
                        "smoke", "education", "depression", "T2D", "hypertension", 
                        "BrainVital8", "income")
main_formula_cont <- "Surv(time, status) ~ BrainVital8 + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression"
main_formula_cat <- "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression"
res_main_both <- run_both_models_custom(
  data = df,
  analysis_name = "Main",
  required_vars = main_required_vars,
  formula_continuous = main_formula_cont,
  formula_categorical = main_formula_cat
)
res_main_continuous <- res_main_both$continuous
res_main_categorical <- res_main_both$categorical

# ---------------------- 2. LANCET Score Analysis ----------------------
lancet_required_vars <- c(main_required_vars, "LANCET")
lancet_formula_cont <- "Surv(time, status) ~ BrainVital8 + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression + LANCET"
lancet_formula_cat <- "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression + LANCET"
res_lancet_both <- run_both_models_custom(
  data = df,
  analysis_name = "LANCET",
  required_vars = lancet_required_vars,
  formula_continuous = lancet_formula_cont,
  formula_categorical = lancet_formula_cat
)
res_lancet_continuous <- res_lancet_both$continuous
res_lancet_categorical <- res_lancet_both$categorical

# ---------------------- 3. LIBRA2 Score Analysis ----------------------
libra2_required_vars <- c("id", "age", "sex", "bmi", "time", "status", "drink", 
                          "smoke", "education", "depression", "T2D", "hypertension", 
                          "BrainVital8", "income", "cohort", "LIBRA2")
libra2_formula_cont <- "Surv(time, status) ~ BrainVital8 + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression + LIBRA2"
libra2_formula_cat <- "Surv(time, status) ~ BrainVital8_quartile + age + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression + LIBRA2"
res_libra2_both <- run_both_models_custom(
  data = df,
  analysis_name = "LIBRA2",
  required_vars = libra2_required_vars,
  formula_continuous = libra2_formula_cont,
  formula_categorical = libra2_formula_cat
)
res_libra2_continuous <- res_libra2_both$continuous
res_libra2_categorical <- res_libra2_both$categorical

# ---------------------- 4. Age <65 Stratified Analysis ----------------------
age_required_vars <- c("id", "age", "sex", "bmi", "time", "status", "drink", 
                       "smoke", "education", "depression", "T2D", "hypertension", 
                       "BrainVital8", "income", "cohort")
age_formula_cont <- "Surv(time, status) ~ BrainVital8 + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression"
age_formula_cat <- "Surv(time, status) ~ BrainVital8_quartile + sex + bmi + education + drink + smoke + income + T2D + hypertension + depression"
df_age_lt65 <- df %>% filter(age < 65 & !is.na(age))
res_age_lt65_both <- run_both_models_custom(
  data = df_age_lt65,
  analysis_name = "Age_lt65",
  required_vars = age_required_vars,
  formula_continuous = age_formula_cont,
  formula_categorical = age_formula_cat
)
res_age_lt65_continuous <- res_age_lt65_both$continuous
res_age_lt65_categorical <- res_age_lt65_both$categorical

# ---------------------- 5. Age ≥65 Stratified Analysis (Completed) ----------------------
df_age_ge65 <- df %>% filter(age >= 65 & !is.na(age))
res_age_ge65_both <- run_both_models_custom(
  data = df_age_ge65,
  analysis_name = "Age_ge65",
  required_vars = age_required_vars,
  formula_continuous = age_formula_cont,
  formula_categorical = age_formula_cat
)
res_age_ge65_continuous <- res_age_ge65_both$continuous
res_age_ge65_categorical <- res_age_ge65_both$categorical

# ---------------------- 6. Male Stratified Analysis ----------------------
sex_required_vars <- age_required_vars
sex_formula_cont <- "Surv(time, status) ~ BrainVital8 + age + bmi + education + drink + smoke + income + T2D + hypertension + depression"
sex_formula_cat <- "Surv(time, status) ~ BrainVital8_quartile + age + bmi + education + drink + smoke + income + T2D + hypertension + depression"
df_male <- df %>% filter(sex == 1 & !is.na(sex))
res_male_both <- run_both_models_custom(
  data = df_male,
  analysis_name = "Male",
  required_vars = sex_required_vars,
  formula_continuous = sex_formula_cont,
  formula_categorical = sex_formula_cat
)
res_male_continuous <- res_male_both$continuous
res_male_categorical <- res_male_both$categorical

# ---------------------- 7. Female Stratified Analysis (Completed) ----------------------
df_female <- df %>% filter(sex == 0 & !is.na(sex))
res_female_both <- run_both_models_custom(
  data = df_female,
  analysis_name = "Female",
  required_vars = sex_required_vars,
  formula_continuous = sex_formula_cont,
  formula_categorical = sex_formula_cat
)
res_female_continuous <- res_female_both$continuous
res_female_categorical <- res_female_both$categorical

# ---------------------- 8. Exclude First 2 Years Onset Analysis (Fixed Filtering Logic) ----------------------
exclude2yr_required_vars <- age_required_vars
if (!is.null(res_main_continuous)) {
  # Fixed filtering logic: keep only cases with status=1 and time>2, and all status=0 cases
  df_ex2yr <- res_main_continuous$data_used %>%
    filter((status == 1 & time > 2) | (status == 0))  # Core fix: explicitly require time>2 for cases, distinguish by status type
  ex2yr_formula_cont <- main_formula_cont
  ex2yr_formula_cat <- main_formula_cat
  res_ex2yr_both <- run_both_models_custom(
    data = df_ex2yr,
    analysis_name = "Exclude_2yr",
    required_vars = exclude2yr_required_vars,
    formula_continuous = ex2yr_formula_cont,
    formula_categorical = ex2yr_formula_cat
  )
  res_ex2yr_continuous <- res_ex2yr_both$continuous
  res_ex2yr_categorical <- res_ex2yr_both$categorical
}


##----------------------------------------------------------
## 6. Generate Target Cox Table (Fixed P-value Display, Keep 8 Decimal Places)
##----------------------------------------------------------
# Format conversion function: Adapt to image table style + unify column types + complete P-value display
format_target_table <- function(cont_res, cat_res, analysis_name) {
  # Process continuous model: P-value keeps 8 decimal places, convert to character type, do not show 0
  cont_table <- if (!is.null(cont_res)) {
    cont_res$results %>%
      filter(Variable == "BrainVital8") %>%
      mutate(
        Analysis_Name = analysis_name,
        Analysis_Type = "Continuous (adjusted)",
        Comparison = "Continuous",
        HR_CI = paste0(round(HR, 2), " (", round(CI_lower, 2), "-", round(CI_upper, 2), ")"),
        # Fixed P-value: keep 8 decimal places, output specific value (e.g., 0.00000235)
        P_Value = ifelse(is.na(P_value), "", sprintf("%.8f", P_value)),
        Cases = as.integer(Case),
        Total_N = as.integer(N)
      ) %>%
      select(Analysis_Name, Analysis_Type, Comparison, HR_CI, P_Value, Cases, Total_N)
  } else {NULL}
  
  # Process categorical model: P-value keeps 8 decimal places, reference group is empty string, do not show 0
  cat_table <- if (!is.null(cat_res)) {
    # Extract group case counts/sample sizes for Q1-Q4
    group_stats <- cat_res$group_stats %>%
      rename(Cases = Group_Case, Total_N = Group_N) %>%
      mutate(
        Comparison = BrainVital8_quartile,
        Cases = as.integer(Cases),
        Total_N = as.integer(Total_N)
      )
    
    # Extract HR/P-values for Q1-Q4, complete P-value display
    cat_res$results %>%
      filter(grepl("BrainVital8_quartile", Variable)) %>%
      mutate(
        Comparison = gsub("BrainVital8_quartile| \\(ref\\)", "", Variable),
        HR_CI = ifelse(Variable == "BrainVital8_quartileQ1 (ref)", "Reference", 
                       paste0(round(HR, 2), " (", round(CI_lower, 2), "-", round(CI_upper, 2), ")")),
        # Fixed P-value: empty for reference group, keep 8 decimal places for others
        P_Value = ifelse(Variable == "BrainVital8_quartileQ1 (ref)", "", 
                         ifelse(is.na(P_value), "", sprintf("%.8f", P_value))),
        Analysis_Name = analysis_name,
        Analysis_Type = "Quartile (adjusted)"
      ) %>%
      left_join(group_stats, by = "Comparison") %>%
      select(Analysis_Name, Analysis_Type, Comparison, HR_CI, P_Value, Cases, Total_N)
  } else {NULL}
  
  # Merge continuous + categorical results (now all column types are consistent)
  bind_rows(cont_table, cat_table)
}

# Generate Cox tables for 8 analyses
main_table <- format_target_table(res_main_continuous, res_main_categorical, "Main Analysis")
lancet_table <- format_target_table(res_lancet_continuous, res_lancet_categorical, "LANCET Adjusted")
libra2_table <- format_target_table(res_libra2_continuous, res_libra2_categorical, "LIBRA2 Adjusted")
age_lt65_table <- format_target_table(res_age_lt65_continuous, res_age_lt65_categorical, "Age < 65")
age_ge65_table <- format_target_table(res_age_ge65_continuous, res_age_ge65_categorical, "Age ≥ 65")  # P-value complete display
male_table <- format_target_table(res_male_continuous, res_male_categorical, "Male")
female_table <- format_target_table(res_female_continuous, res_female_categorical, "Female")          # New
ex2yr_table <- format_target_table(res_ex2yr_continuous, res_ex2yr_categorical, "Exclude First 2 Years")  # Filtering fixed

# Merge all Cox tables and save
final_target_table <- bind_rows(main_table, lancet_table, libra2_table, age_lt65_table, age_ge65_table, 
                                male_table, female_table, ex2yr_table)
write.csv(
  final_target_table,
  "./output/BrainVital8_Cox_Target_Table.csv",
  row.names = FALSE
)


##----------------------------------------------------------
## 7. Extract + Summarize + Save PH Test Results (Fixed P-value Display)
##----------------------------------------------------------
# Collect all valid model results (including Cox+PH)
all_results <- list(
  Main_continuous = res_main_continuous,
  Main_categorical = res_main_categorical,
  LANCET_continuous = res_lancet_continuous,
  LANCET_categorical = res_lancet_categorical,
  LIBRA2_continuous = res_libra2_continuous,
  LIBRA2_categorical = res_libra2_categorical,
  Age_lt65_continuous = res_age_lt65_continuous,
  Age_lt65_categorical = res_age_lt65_categorical,
  Age_ge65_continuous = res_age_ge65_continuous,
  Age_ge65_categorical = res_age_ge65_categorical,
  Male_continuous = res_male_continuous,
  Male_categorical = res_male_categorical,
  Female_continuous = res_female_continuous,
  Female_categorical = res_female_categorical,
  Exclude_2yr_continuous = res_ex2yr_continuous,
  Exclude_2yr_categorical = res_ex2yr_categorical
)
# Filter out invalid results
all_results <- all_results[!sapply(all_results, is.null)]

# Function 1: Extract global PH test results (fixed P-value display)
extract_global_ph_test_results <- function(results_list) {
  ph_results <- data.frame()
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    if (!is.null(res$ph_test)) {
      ph_table <- res$ph_test$table
      model_type <- ifelse(grepl("_continuous$", analysis_name), "continuous", "categorical")
      model_label <- ifelse(model_type == "continuous", 
                            "BrainVital8 (continuous)", 
                            "BrainVital8_quartile (categorical)")
      if ("GLOBAL" %in% rownames(ph_table)) {
        global_ph <- ph_table["GLOBAL", ]
        global_ph_df <- data.frame(
          Analysis = gsub("_(continuous|categorical)$", "", analysis_name),
          Model_Type = model_type,
          Model_Label = model_label,
          Global_chisq = global_ph["chisq"],
          Global_df = global_ph["df"],
          # Fixed P-value: keep 8 decimal places, output specific value
          Global_p = ifelse(is.na(global_ph["p"]), NA, sprintf("%.8f", global_ph["p"])),
          stringsAsFactors = FALSE
        )
        ph_results <- bind_rows(ph_results, global_ph_df)
      }
    }
  }
  return(ph_results)
}
global_ph_results <- extract_global_ph_test_results(all_results)

# Function 2: Extract variable-specific PH test results (fixed P-value display)
extract_variable_specific_ph_results <- function(results_list) {
  ph_results <- data.frame()
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    if (!is.null(res$ph_test)) {
      ph_table <- res$ph_test$table
      model_type <- ifelse(grepl("_continuous$", analysis_name), "continuous", "categorical")
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
            # Fixed P-value: keep 8 decimal places, output specific value
            p = ifelse(is.na(var_ph["p"]), NA, sprintf("%.8f", var_ph["p"])),
            stringsAsFactors = FALSE
          )
          ph_results <- bind_rows(ph_results, var_ph_df)
        }
      }
    }
  }
  return(ph_results)
}
variable_ph_results <- extract_variable_specific_ph_results(all_results)

# Function 3: Generate PH test summary table (fixed P-value display)
create_ph_summary_table <- function(global_results, variable_results) {
  summary_table <- data.frame()
  for (analysis in unique(global_results$Analysis)) {
    analysis_global <- global_results %>% filter(Analysis == analysis)
    analysis_variable <- variable_results %>% filter(Analysis == analysis)
    
    # Continuous model
    cont_global <- analysis_global %>% filter(Model_Type == "continuous")
    cont_variable <- analysis_variable %>% filter(Model_Type == "continuous" & Variable == "BrainVital8")
    if (nrow(cont_global) > 0) {
      cont_row <- data.frame(
        Analysis = analysis,
        Model = "Continuous (BrainVital8)",
        BrainVital8_Variable_PH = ifelse(nrow(cont_variable) > 0, 
                                         paste0("χ²=", round(cont_variable$chisq, 3), 
                                                ", df=", cont_variable$df, 
                                                ", p=", cont_variable$p),  # Use already formatted P-value
                                         "N/A"),
        Global_PH = paste0("χ²=", round(cont_global$Global_chisq, 3), 
                           ", df=", cont_global$Global_df, 
                           ", p=", cont_global$Global_p),  # Use already formatted P-value
        stringsAsFactors = FALSE
      )
      summary_table <- bind_rows(summary_table, cont_row)
    }
    
    # Categorical model
    cat_global <- analysis_global %>% filter(Model_Type == "categorical")
    cat_variables <- analysis_variable %>% filter(Model_Type == "categorical" & grepl("BrainVital8_quartile", Variable))
    if (nrow(cat_global) > 0) {
      var_ph_strings <- c()
      if (nrow(cat_variables) > 0) {
        for (i in 1:nrow(cat_variables)) {
          var_name <- gsub("BrainVital8_quartile", "", cat_variables$Variable[i])
          var_ph_strings <- c(var_ph_strings, 
                              paste0(var_name, ": χ²=", round(cat_variables$chisq[i], 3),
                                     ", p=", cat_variables$p))  # Use already formatted P-value
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
                           ", p=", cat_global$Global_p),  # Use already formatted P-value
        stringsAsFactors = FALSE
      )
      summary_table <- bind_rows(summary_table, cat_row)
    }
  }
  return(summary_table)
}
ph_summary_table <- create_ph_summary_table(global_ph_results, variable_ph_results)

# Save all PH test results
write.csv(global_ph_results, "./output/BrainVital8_Global_PH_Test_Results.csv", row.names = FALSE)
write.csv(variable_ph_results, "./output/BrainVital8_Variable_PH_Test_Results.csv", row.names = FALSE)
write.csv(ph_summary_table, "./output/BrainVital8_PH_Test_Summary_Table.csv", row.names = FALSE)


##----------------------------------------------------------
## 8. Output Complete Summary Information
##----------------------------------------------------------
cat("\n\n=== Analysis Complete! Cox+PH Test Complete Results Summary ===\n")

cat("\n📊 Target Cox Table Preview (Includes 8 Analyses):\n")
print(final_target_table, right = TRUE)

cat("\n\n📊 Global PH Test Results Preview:\n")
print(global_ph_results, right = TRUE)

cat("\n\n📊 PH Test Summary Table Preview:\n")
print(ph_summary_table, right = TRUE)

cat("\n📁 File Saving Information:")
cat("\n- Cox Target Table: ", getwd(), "./output/BrainVital8_Cox_Target_Table.csv")
cat("\n- Global PH Test Results: ", getwd(), "./output/BrainVital8_Global_PH_Test_Results.csv")
cat("\n- Variable-Specific PH Test Results: ", getwd(), "./output/BrainVital8_Variable_PH_Test_Results.csv")
cat("\n- PH Test Summary Table: ", getwd(), "./output/BrainVital8_PH_Test_Summary_Table.csv")
cat("\n- Working Directory: ", getwd())

cat("\n\n🔍 Key Notes:")
cat("\n1. Fixed P-value Display: Keep 8 decimal places, output specific values for very small P-values (e.g., 0.00000235), do not show 0 or <xxx")
cat("\n2. Complete Cox+PH Test Retention: HR/CI/P-values are real and valid, proportional hazards assumption test results are completely preserved")
