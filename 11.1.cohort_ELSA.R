############################################
# Load required packages
############################################
library(dplyr)
library(survival)

############################################
# Read data
############################################
final_data_clean <- read.csv(
  "./input/ELSA.csv"
)

cat("=== Data summary ===\n")
cat("Total sample size:", nrow(final_data_clean), "\n")
cat("Total events:", sum(final_data_clean$status), "\n")
cat("Event rate:", round(mean(final_data_clean$status) * 100, 2), "%\n")

############################################
# Cox regression analysis function
############################################
run_analysis <- function(data, analysis_name, covariates_formula,
                         exposure_var = "BrainVital8") {
  
  covariate_vars <- all.vars(as.formula(paste("~", covariates_formula)))
  required_cols <- c(exposure_var, "time", "status", covariate_vars)
  
  data <- data %>%
    select(all_of(required_cols)) %>%
    na.omit() %>%
    mutate(
      total_score_group_q = cut(
        BrainVital8,
        breaks = quantile(
          BrainVital8,
          probs = c(0, 0.25, 0.5, 0.75, 1),
          na.rm = TRUE
        ),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE
      )
    )
  
  # Sample size and cases by quartile
  group_stats <- data %>%
    group_by(total_score_group_q) %>%
    summarise(N = n(), Cases = sum(status), .groups = "drop")
  
  total_n <- nrow(data)
  total_cases <- sum(data$status)
  
  # Continuous exposure model
  cox_cont <- coxph(
    as.formula(
      paste("Surv(time, status) ~", exposure_var, "+", covariates_formula)
    ),
    data = data
  )
  cont_sum <- summary(cox_cont)
  
  # Quartile-based model
  cox_q <- coxph(
    as.formula(
      paste("Surv(time, status) ~ total_score_group_q +", covariates_formula)
    ),
    data = data
  )
  q_sum <- summary(cox_q)
  
  res <- list()
  
  # Continuous results
  res[[1]] <- data.frame(
    Analysis = analysis_name,
    Model = "Continuous",
    Exposure = exposure_var,
    Comparison = "Per unit increase",
    N = total_n,
    Cases = total_cases,
    HR = exp(coef(cox_cont)[exposure_var]),
    HR_lower = exp(confint(cox_cont)[exposure_var, 1]),
    HR_upper = exp(confint(cox_cont)[exposure_var, 2]),
    p_value = cont_sum$coefficients[exposure_var, "Pr(>|z|)"]
  )
  
  # Q1 reference group
  q1 <- group_stats %>% filter(total_score_group_q == "Q1")
  res[[2]] <- data.frame(
    Analysis = analysis_name,
    Model = "Quartile",
    Exposure = exposure_var,
    Comparison = "Q1 (Reference)",
    N = q1$N,
    Cases = q1$Cases,
    HR = 1,
    HR_lower = NA,
    HR_upper = NA,
    p_value = NA
  )
  
  # Q2–Q4 vs Q1
  for (i in 2:4) {
    term <- paste0("total_score_group_qQ", i)
    qi <- group_stats %>%
      filter(total_score_group_q == paste0("Q", i))
    
    res[[i + 1]] <- data.frame(
      Analysis = analysis_name,
      Model = "Quartile",
      Exposure = exposure_var,
      Comparison = paste0("Q", i, " vs Q1"),
      N = qi$N,
      Cases = qi$Cases,
      HR = exp(coef(cox_q)[term]),
      HR_lower = exp(confint(cox_q)[term, 1]),
      HR_upper = exp(confint(cox_q)[term, 2]),
      p_value = q_sum$coefficients[term, "Pr(>|z|)"]
    )
  }
  
  bind_rows(res)
}

############################################
# Main analysis and subgroup analyses
############################################
main_results <- run_analysis(
  final_data_clean, "Main_Analysis",
  "age + sex + bmi + drink + hypertension"
)

male_data <- final_data_clean %>% filter(sex == 1)
male_results <- run_analysis(
  male_data, "Male_Subgroup",
  "age + bmi + drink + hypertension"
)

female_data <- final_data_clean %>% filter(sex == 0)
female_results <- run_analysis(
  female_data, "Female_Subgroup",
  "age + bmi + drink + hypertension"
)

############################################
# Exclude dementia cases within first 2 years
############################################
exclude_2year_data <- final_data_clean %>%
  filter(!(status == 1 & time <= 2))

exclude_2year_results <- run_analysis(
  exclude_2year_data, "Exclude_2year",
  "age + sex + bmi + drink + hypertension"
)

############################################
# Additional adjusted models
############################################
analysis6_results <- run_analysis(
  final_data_clean, "Analysis6_LANCET",
  "age + sex + bmi + drink + hypertension + LANCET"
)

analysis7_results <- run_analysis(
  final_data_clean, "Analysis7_LIBRA2",
  "age + sex + bmi + drink + hypertension + LIBRA2"
)

############################################
# Age-stratified analyses
############################################
age_ge65_data <- final_data_clean %>% filter(age >= 65)
age_ge65_results <- run_analysis(
  age_ge65_data, "Age_ge65",
  "sex + bmi + drink + hypertension"
)

age_lt65_data <- final_data_clean %>% filter(age < 65)
age_lt65_results <- run_analysis(
  age_lt65_data, "Age_lt65",
  "sex + bmi + drink + hypertension"
)

############################################
# Combine results and export
############################################
final_results <- bind_rows(
  main_results,
  male_results,
  female_results,
  exclude_2year_results,
  analysis6_results,
  analysis7_results,
  age_ge65_results,
  age_lt65_results
) %>%
  mutate(
    N_Cases = paste0(Cases, "/", N),
    HR_95CI = ifelse(
      is.na(p_value),
      "1.000 (Reference)",
      sprintf("%.2f (%.2f–%.2f)", HR, HR_lower, HR_upper)
    ),
    p_value_scientific = ifelse(
      is.na(p_value),
      "",
      format(p_value, scientific = TRUE, digits = 3)
    )
  ) %>%
  select(
    Analysis, Model, Exposure, Comparison,
    N_Cases, HR_95CI, p_value, p_value_scientific
  )
library(openxlsx)
write.xlsx(
  final_results,
  file = "./output/ELSA_vital8_results.xlsx",
  rowNames = FALSE
)


############################################
# Generate final model datasets
# (consistent with HR and PH tests)
############################################
generate_model_data <- function(data,
                                analysis_name,
                                covariates_formula,
                                exposure_var = "BrainVital8",
                                id_var = "idauniqc",
                                save_dir = "./output/model_datasets") {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  covariate_vars <- all.vars(as.formula(paste("~", covariates_formula)))
  required_cols <- unique(
    c(id_var, exposure_var, "time", "status", covariate_vars)
  )
  
  final_model_data <- data %>%
    select(all_of(required_cols)) %>%
    na.omit() %>%
    rename(
      id = all_of(id_var)
    )
  
  write.csv(
    final_model_data,
    paste0(save_dir, "/", analysis_name, "_final_dataset.csv"),
    row.names = FALSE
  )
  
  cat(
    "✔", analysis_name,
    "| N =", nrow(final_model_data),
    "| Cases =", sum(final_model_data$status), "\n"
  )
  
  return(final_model_data)
}

############################################
# Generate datasets for each analysis
############################################
main_model_data <- generate_model_data(
  final_data_clean,
  "Main_Analysis",
  "age + sex + bmi + drink + hypertension"
)

male_model_data <- generate_model_data(
  male_data,
  "Male_Subgroup",
  "age + bmi + drink + hypertension"
)

female_model_data <- generate_model_data(
  female_data,
  "Female_Subgroup",
  "age + bmi + drink + hypertension"
)

exclude_2year_model_data <- generate_model_data(
  exclude_2year_data,
  "Exclude_2year",
  "age + sex + bmi + drink + hypertension"
)

lancet_model_data <- generate_model_data(
  final_data_clean,
  "Analysis6_LANCET",
  "age + sex + bmi + drink + hypertension + LANCET"
)

libra2_model_data <- generate_model_data(
  final_data_clean,
  "Analysis7_LIBRA2",
  "age + sex + bmi + drink + hypertension + LIBRA2"
)

age_ge65_model_data <- generate_model_data(
  age_ge65_data,
  "Age_ge65",
  "sex + bmi + drink + hypertension"
)

age_lt65_model_data <- generate_model_data(
  age_lt65_data,
  "Age_lt65",
  "sex + bmi + drink + hypertension"
)

############################################
# Proportional hazards (PH) tests
############################################
extract_ph_brainvital8_global <- function(
    data, analysis_name, covariates_formula,
    exposure_var = "BrainVital8"
) {
  
  cox_model <- coxph(
    as.formula(
      paste("Surv(time, status) ~", exposure_var, "+", covariates_formula)
    ),
    data = data
  )
  
  ph_test <- cox.zph(cox_model)
  
  df_bv8 <- data.frame(
    Analysis = analysis_name,
    Variable = exposure_var,
    Chi_square = ph_test$table[exposure_var, "chisq"],
    df = ph_test$table[exposure_var, "df"],
    p_value = ph_test$table[exposure_var, "p"]
  )
  
  df_global <- data.frame(
    Analysis = analysis_name,
    Variable = "GLOBAL",
    Chi_square = ph_test$table["GLOBAL", "chisq"],
    df = ph_test$table["GLOBAL", "df"],
    p_value = ph_test$table["GLOBAL", "p"]
  )
  
  bind_rows(df_bv8, df_global)
}

############################################
# PH test for quartile model + GLOBAL
############################################
extract_ph_quartile_global <- function(
    data, analysis_name, covariates_formula
) {
  
  data <- data %>%
    mutate(
      total_score_group_q = cut(
        BrainVital8,
        breaks = quantile(
          BrainVital8,
          probs = c(0, 0.25, 0.5, 0.75, 1),
          na.rm = TRUE
        ),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE
      )
    )
  
  cox_model <- coxph(
    as.formula(
      paste("Surv(time, status) ~ total_score_group_q +", covariates_formula)
    ),
    data = data
  )
  
  ph_test <- cox.zph(cox_model)
  
  df_q <- data.frame(
    Analysis = analysis_name,
    Variable = "total_score_group_q",
    Chi_square = ph_test$table["total_score_group_q", "chisq"],
    df = ph_test$table["total_score_group_q", "df"],
    p_value = ph_test$table["total_score_group_q", "p"]
  )
  
  df_global <- data.frame(
    Analysis = analysis_name,
    Variable = "GLOBAL",
    Chi_square = ph_test$table["GLOBAL", "chisq"],
    df = ph_test$table["GLOBAL", "df"],
    p_value = ph_test$table["GLOBAL", "p"]
  )
  
  bind_rows(df_q, df_global)
}

############################################
# List of all analysis models
############################################
analysis_list <- list(
  list(main_model_data, "Main_Analysis",
       "age + sex + bmi + drink + hypertension"),
  list(male_model_data, "Male_Subgroup",
       "age + bmi + drink + hypertension"),
  list(female_model_data, "Female_Subgroup",
       "age + bmi + drink + hypertension"),
  list(exclude_2year_model_data, "Exclude_2year",
       "age + sex + bmi + drink + hypertension"),
  list(lancet_model_data, "Analysis6_LANCET",
       "age + sex + bmi + drink + hypertension + LANCET"),
  list(libra2_model_data, "Analysis7_LIBRA2",
       "age + sex + bmi + drink + hypertension + LIBRA2"),
  list(age_ge65_model_data, "Age_ge65",
       "sex + bmi + drink + hypertension"),
  list(age_lt65_model_data, "Age_lt65",
       "sex + bmi + drink + hypertension")
)

############################################
# PH tests: continuous BrainVital8 + GLOBAL
############################################
ph_continuous_all <- bind_rows(
  lapply(
    analysis_list,
    function(x)
      extract_ph_brainvital8_global(x[[1]], x[[2]], x[[3]])
  )
)

############################################
# PH tests: quartiles + GLOBAL
############################################
ph_quartile_all <- bind_rows(
  lapply(
    analysis_list,
    function(x)
      extract_ph_quartile_global(x[[1]], x[[2]], x[[3]])
  )
)

############################################
# Combine and export PH results
############################################
ph_combined <- bind_rows(ph_continuous_all, ph_quartile_all)

write.csv(
  ph_combined,
  "./output/ELSA_vital8_dementia_PH_AllModels_Combined.csv",
  row.names = FALSE
)

cat(
  "\n✔ Continuous BrainVital8 + quartile models + GLOBAL PH tests exported\n"
)
