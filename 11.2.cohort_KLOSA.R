# =========================================================
# Load required packages
# =========================================================
library(readxl)
library(dplyr)
library(rio)
library(survival)

# =========================================================
# Import datasets
# =========================================================
aligned <- read_excel("./input/BrainVital8_scores_aligned.xlsx")
behavior <- read_excel("./input/behavior_score0827.xlsx")
lancet <- read.csv("./input/lancet_scored.csv")
libra2 <- read.csv("./input/libra2_scored.csv")
klosa <- read_excel("./input/KLoSA_data.xlsx")
cov <- read_excel("./input/Covariates.xlsx")

# =========================================================
# Data preprocessing workflow
# =========================================================
klosa <- klosa %>% filter(pid %in% behavior$pid)              # Retain participants with behavioral data
klosa <- klosa %>% mutate(pid = as.numeric(pid))              # Harmonize identifier type
cov <- cov %>% mutate(pid = as.numeric(pid))
klosa <- klosa %>% left_join(cov %>% select(pid, hibpe, depression), by = "pid")
c <- klosa %>% left_join(aligned %>% select(pid, BrainVital8_aligned), by = "pid")

# Removed median-based dichotomization of itot
# c$itot <- ifelse(c$itot >= median(c$itot, na.rm=TRUE), 1, 0)   # removed

c_lancet <- c %>% left_join(lancet %>% select(pid, lancet), by = "pid")
c_libra2 <- c %>% left_join(libra2 %>% select(pid, libra2), by = "pid")

# =========================================================
# Categorical variable encoding
# =========================================================
c$gender1 <- factor(ifelse(c$gender1 == 5, 0, 1), levels=c(0,1), labels=c("Female","Male"))
c$smoken <- factor(c$smoken, levels=c(0,1), labels=c("No","Yes"))
c$drink  <- factor(c$drink,  levels=c(0,1), labels=c("No","Yes"))
c$hibpe  <- factor(c$hibpe,  levels=c(0,1), labels=c("No","Yes"))
c$diabe  <- factor(c$diabe,  levels=c(0,1), labels=c("No","Yes"))
c$depression <- factor(c$depression, levels=c(0,1), labels=c("No","Yes"))
c$edu <- factor(
  c$edu,
  levels=c(1,2,3,4),
  ordered=TRUE,
  labels=c("Primary or below","Middle school","High school","College or above")
)

# =========================================================
# Missing-data imputation strategy
# =========================================================
impute_missing <- function(df){
  cat_vars <- c("gender1","edu","smoken","drink","hibpe","diabe","depression")
  
  # Encode missing categories to preserve sample size
  for(v in cat_vars){
    if(v %in% names(df)){
      if(is.factor(df[[v]]) || is.ordered(df[[v]])){
        levels(df[[v]]) <- c(levels(df[[v]]), "Missing")
        df[[v]][is.na(df[[v]])] <- "Missing"
      } else { 
        df[[v]][is.na(df[[v]])] <- 9 
      }
    }
  }
  
  # Median-based imputation for continuous traits
  cont_vars <- c("a002_age","bmi","BrainVital8_aligned")
  for(v in cont_vars){
    if(v %in% names(df)){
      med <- median(df[[v]], na.rm=TRUE)
      df[[v]][is.na(df[[v]])] <- med
    }
  }
  return(df)
}

c <- impute_missing(c)
c_lancet <- impute_missing(c_lancet)
c_libra2 <- impute_missing(c_libra2)

# =========================================================
# Stratified datasets
# =========================================================
c_male <- c %>% filter(gender1=="Male")               # Male participants
c_female <- c %>% filter(gender1=="Female")           # Female participants
c_2year <- c %>% filter(!(dementia_diagnose==1 & survival_years <= 2))  
# Exclude early-onset dementia within 2 years to mitigate reverse causation

# =========================================================
# Re-factorization for derived datasets
# =========================================================
factorize_vars <- function(df){
  df$gender1 <- factor(ifelse(df$gender1==5,0,1), levels=c(0,1), labels=c("Female","Male"))
  df$smoken <- factor(df$smoken, levels=c(0,1), labels=c("No","Yes"))
  df$drink  <- factor(df$drink,  levels=c(0,1), labels=c("No","Yes"))
  df$hibpe  <- factor(df$hibpe,  levels=c(0,1), labels=c("No","Yes"))
  df$diabe  <- factor(df$diabe,  levels=c(0,1), labels=c("No","Yes"))
  df$depression <- factor(df$depression, levels=c(0,1), labels=c("No","Yes"))
  df$edu <- factor(
    df$edu,
    levels=c(1,2,3,4),
    ordered=TRUE,
    labels=c("Primary or below","Middle school","High school","College or above")
  )
  return(df)
}

c_lancet <- factorize_vars(c_lancet)
c_libra2 <- factorize_vars(c_libra2)

c_libra2 <- c_libra2[!is.na(c_libra2$libra2), ]      # Retain complete LIBRA-2 scores
c_lancet <- c_lancet[!is.na(c_lancet$lancet), ]      # Retain complete Lancet dementia scores

# =========================================================
# Covariate specification for modeling
# =========================================================
cov_c <- c("a002_age","bmi","gender1","drink","hibpe")
cov_c_mf <- cov_c[cov_c != "gender1"]               # Covariates for sex-stratified analyses
cov_c_lancet <- c("a002_age","bmi","gender1","drink","hibpe","lancet")
cov_c_libra2 <- c("a002_age","bmi","gender1","drink","hibpe","libra2")

# =========================================================
# Quartile construction and Cox analysis utilities
# =========================================================

# Derive quartiles based on the empirical distribution of BrainVital8 scores
generate_quartile <- function(df){
  df$BrainVital8_aligned_quartile <- cut(
    df$BrainVital8_aligned,
    breaks = quantile(df$BrainVital8_aligned, probs=c(0,0.25,0.5,0.75,1), na.rm=TRUE),
    labels = c("Q1","Q2","Q3","Q4"), 
    include.lowest = TRUE
  )
  return(df)
}

# Extract HR estimates, confidence intervals, and PH-test statistics from a Cox model
cox_result_table <- function(cox_model, var_name="BrainVital8_aligned"){
  model_summary <- summary(cox_model)
  coef_idx <- which(names(cox_model$coefficients)==var_name)
  
  hr <- exp(cox_model$coefficients[coef_idx])
  ci_lower <- model_summary$conf.int[coef_idx,"lower .95"]
  ci_upper <- model_summary$conf.int[coef_idx,"upper .95"]
  p_value <- model_summary$coefficients[coef_idx,"Pr(>|z|)"]
  
  # Schoenfeld residual test for proportional hazards assumption
  ph_test <- tryCatch(cox.zph(cox_model), error=function(e) NULL)
  if(!is.null(ph_test)){
    ph_p <- ph_test$table[var_name,"p"]
    ph_global_p <- ph_test$table["GLOBAL","p"]
  } else { 
    ph_p <- ""
    ph_global_p <- ""
  }
  
  data.frame(
    Variable = var_name,
    HR = sprintf("%.2f", hr),
    `95% CI` = sprintf("%.2f-%.2f", ci_lower, ci_upper),
    `P value` = format(p_value, scientific=TRUE, digits=2),
    `PH P value` = ph_p,
    `PH Global P value` = ph_global_p,
    Cases = model_summary$nevent,
    `Total N` = model_summary$n,
    stringsAsFactors = FALSE
  )
}

# Generate quartile-specific Cox estimates relative to the reference quartile
extract_quartile_results <- function(cox_model, data, quartile_var="BrainVital8_aligned_quartile"){
  model_summary <- summary(cox_model)
  coef_names <- names(cox_model$coefficients)
  quartile_coefs <- coef_names[grepl(paste0("^", quartile_var), coef_names)]
  
  results <- data.frame()
  
  # Reference quartile summary
  ref_level <- levels(data[[quartile_var]])[1]
  ref_cases <- sum(data$dementia_diagnose[data[[quartile_var]] == ref_level])
  ref_total <- sum(data[[quartile_var]] == ref_level)
  
  results <- rbind(
    results,
    data.frame(
      Comparison = paste0(ref_level, " (Ref)"),
      HR = "1.00",
      `95% CI` = "Reference",
      P_value = "",
      Cases = ref_cases,
      `Total N` = ref_total,
      stringsAsFactors = FALSE
    )
  )
  
  # Estimates for non-reference quartiles
  for(coef in quartile_coefs){
    comparison <- gsub(quartile_var, "", coef)
    hr <- exp(cox_model$coefficients[coef])
    ci_lower <- model_summary$conf.int[coef,"lower .95"]
    ci_upper <- model_summary$conf.int[coef,"upper .95"]
    p_value <- model_summary$coefficients[coef,"Pr(>|z|)"]
    
    cases <- sum(data$dementia_diagnose[data[[quartile_var]]==comparison])
    total <- sum(data[[quartile_var]]==comparison)
    
    results <- rbind(
      results,
      data.frame(
        Comparison = comparison,
        HR = sprintf("%.2f", hr),
        `95% CI` = sprintf("%.2f-%.2f", ci_lower, ci_upper),
        P_value = format(p_value, scientific=TRUE, digits=2),
        Cases = cases,
        `Total N` = total,
        stringsAsFactors = FALSE
      )
    )
  }
  return(results)
}

# Combine continuous-model and quartile-model output into a unified results table
combine_cont_quartile <- function(cont, quartile){
  cont$Model <- "Continuous"
  quartile$Model <- "Quartile"
  
  # Insert a blank separator row for readability
  empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=max(ncol(cont), ncol(quartile))))
  colnames(empty_row) <- colnames(cont)
  
  bind_rows(cont, empty_row, quartile)
}

# Standardize Cox regression tables (harmonize naming, merge columns)
format_cox_output <- function(df){
  colnames(df) <- trimws(colnames(df))
  
  # Consolidate case count columns
  cases_col <- if("Cases" %in% colnames(df)) "Cases" else if("Cases." %in% colnames(df)) "Cases." else NULL
  total_col <- if("Total.N" %in% colnames(df)) "Total.N" else if("Total N" %in% colnames(df)) "Total N" else NULL
  
  if(!is.null(cases_col) & !is.null(total_col)){
    df$N_Cases <- paste0(df[[cases_col]], "/", df[[total_col]])
    df[[cases_col]] <- NULL
    df[[total_col]] <- NULL
  }
  
  # Combine HR and CI into a single reporting column
  hr_col <- if("HR" %in% colnames(df)) "HR" else NULL
  ci_col <- if("X95..CI" %in% colnames(df)) "X95..CI" else if("95% CI" %in% colnames(df)) "95% CI" else NULL
  
  if(!is.null(hr_col) & !is.null(ci_col)){
    df$`HR(95% CI)` <- paste0(df[[hr_col]], " (", df[[ci_col]], ")")
    df[[hr_col]] <- NULL
    df[[ci_col]] <- NULL
  }
  
  return(df)
}
cox_ph_continuous <- function(cox_model, var_name){
  # Schoenfeld residual test for continuous exposure
  ph <- cox.zph(cox_model)
  tab <- ph$table
  data.frame(
    Variable = var_name,
    PH_P_value = if(var_name %in% rownames(tab)) tab[var_name,"p"] else NA,
    Chi_sq = if(var_name %in% rownames(tab)) tab[var_name,"chisq"] else NA,
    df = if(var_name %in% rownames(tab)) tab[var_name,"df"] else NA,
    PH_Global_P_value = tab["GLOBAL","p"],
    Chi_sq_Global = tab["GLOBAL","chisq"],
    df_Global = tab["GLOBAL","df"],
    stringsAsFactors = FALSE
  )
}

cox_ph_quartile <- function(cox_model, quartile_var){
  # Schoenfeld residual test for quartile-coded exposure
  ph <- cox.zph(cox_model)
  tab <- ph$table
  data.frame(
    Variable = quartile_var,
    PH_P_value = if(quartile_var %in% rownames(tab)) tab[quartile_var,"p"] else NA,
    Chi_sq = if(quartile_var %in% rownames(tab)) tab[quartile_var,"chisq"] else NA,
    df = if(quartile_var %in% rownames(tab)) tab[quartile_var,"df"] else NA,
    PH_Global_P_value = tab["GLOBAL","p"],
    Chi_sq_Global = tab["GLOBAL","chisq"],
    df_Global = tab["GLOBAL","df"],
    stringsAsFactors = FALSE
  )
}

# =========================================================
# Cox regression wrapper for continuous and quartile models
# =========================================================
run_cox <- function(df, covs){
  df <- generate_quartile(df)   # Construct quartiles prior to model fitting
  
  # Continuous exposure model
  cont_model <- coxph(
    as.formula(paste0(
      "Surv(survival_years, dementia_diagnose) ~ BrainVital8_aligned + ",
      paste(covs, collapse = "+")
    )), 
    data = df
  )
  
  # Quartile exposure model
  quart_model <- coxph(
    as.formula(paste0(
      "Surv(survival_years, dementia_diagnose) ~ BrainVital8_aligned_quartile + ",
      paste(covs, collapse = "+")
    )), 
    data = df
  )
  
  # Compile model outputs
  res_cont <- cox_result_table(cont_model)
  res_quart <- extract_quartile_results(quart_model, df)
  res <- combine_cont_quartile(res_cont, res_quart)
  res <- format_cox_output(res)
  
  # Proportional hazards diagnostics
  ph_cont <- cox_ph_continuous(cont_model, "BrainVital8_aligned")
  ph_quart <- cox_ph_quartile(quart_model, "BrainVital8_aligned_quartile")
  ph <- bind_rows(ph_cont, ph_quart)
  
  list(cox = res, ph = ph)
}

# =========================================================
# Model execution across all analytic subsets
# =========================================================
res_c <- run_cox(c, cov_c)                       # Full sample
res_c_male <- run_cox(c_male, cov_c_mf)          # Male subgroup
res_c_female <- run_cox(c_female, cov_c_mf)      # Female subgroup
res_c_2year <- run_cox(c_2year, cov_c)           # Excluding early dementia cases
res_c_lancet <- run_cox(c_lancet, cov_c_lancet)  # Model adjusted for Lancet risk score
res_c_libra2 <- run_cox(c_libra2, cov_c_libra2)  # Model adjusted for LIBRA-2 score

# Age-stratified analyses
c_ge65 <- c %>% filter(a002_age >= 65)           # Older participants
c_lt65 <- c %>% filter(a002_age < 65)            # Younger participants

cov_age <- cov_c[cov_c != "a002_age"]            # Remove age from covariates in age-stratified models

res_c_ge65 <- run_cox(c_ge65, cov_age)           # Age ≥ 65
res_c_lt65 <- run_cox(c_lt65, cov_age)           # Age < 65

# =========================================================
# Run all datasets and append Dataset label
# =========================================================
res_c <- run_cox(c, cov_c)
res_c$cox$Dataset <- "Main"
res_c$ph$Dataset <- "Main"

res_c_male <- run_cox(c_male, cov_c_mf)
res_c_male$cox$Dataset <- "Male"
res_c_male$ph$Dataset <- "Male"

res_c_female <- run_cox(c_female, cov_c_mf)
res_c_female$cox$Dataset <- "Female"
res_c_female$ph$Dataset <- "Female"

res_c_2year <- run_cox(c_2year, cov_c)
res_c_2year$cox$Dataset <- "Ex2years"
res_c_2year$ph$Dataset <- "Ex2years"

res_c_lancet <- run_cox(c_lancet, cov_c_lancet)
res_c_lancet$cox$Dataset <- "Lancet"
res_c_lancet$ph$Dataset <- "Lancet"

res_c_libra2 <- run_cox(c_libra2, cov_c_libra2)
res_c_libra2$cox$Dataset <- "Libra2"
res_c_libra2$ph$Dataset <- "Libra2"

res_c_ge65 <- run_cox(c_ge65, cov_c[cov_c != "a002_age"])
res_c_ge65$cox$Dataset <- "Ge_65"
res_c_ge65$ph$Dataset <- "Ge_65"

res_c_lt65 <- run_cox(c_lt65, cov_c[cov_c != "a002_age"])
res_c_lt65$cox$Dataset <- "Lt_65"
res_c_lt65$ph$Dataset <- "Lt_65"

# =========================================================
# Merge Cox results for the first worksheet
# =========================================================
all_cox_results <- bind_rows(
  res_c$cox,
  res_c_male$cox,
  res_c_female$cox,
  res_c_2year$cox,
  res_c_lancet$cox,
  res_c_libra2$cox,
  res_c_ge65$cox,
  res_c_lt65$cox
)

# Merge PH test results for the second worksheet
all_ph_results <- bind_rows(
  res_c$ph,
  res_c_male$ph,
  res_c_female$ph,
  res_c_2year$ph,
  res_c_lancet$ph,
  res_c_libra2$ph,
  res_c_ge65$ph,
  res_c_lt65$ph
)


# =========================================================
# Export to Excel
# =========================================================
export(
  list(
    Cox_Results = all_cox_results,
    PH_Results = all_ph_results
  ),
  file="./output/BrainVital8_all_cox.xlsx"
)
