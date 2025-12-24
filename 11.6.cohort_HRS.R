#### 1. MAIN ANALYSIS ####
### Setup: Load Libraries ###
library(survival)
library(dplyr)
library(tidyr)


### Data Preparation ##
baseline <- read.csv("./input/data.csv",header = T) 

# Change the variable name BrainVital8 to score
colnames(baseline)[colnames(baseline) == "BrainVital8"] <- "score"

### Cox Proportional-Hazards Models ###
### Model 1: Score as a continuous variable ###
cox_continuous <- coxph(
  Surv(time, status) ~ score + age + sex + bmi + drink + hypertension,
  data = baseline
)

# Extract the results of Model 1
continuous_summary <- summary(cox_continuous)
continuous_hr <- exp(cox_continuous$coefficients["score"])
continuous_ci <- exp(confint(cox_continuous)["score", ])
continuous_p <- continuous_summary$coefficients["score", "Pr(>|z|)"]
total_n <- nrow(baseline)
cases <- sum(baseline$status)

### Model 2: Score as a categorical (quartile) variable ###
baseline <- baseline %>%
  mutate(
    score_quartile = cut(score,
                         breaks = quantile(score, probs = c(0, 0.25, 0.5, 0.75, 1)),
                         include.lowest = TRUE,
                         labels = c("Q1", "Q2", "Q3", "Q4"))
  )

cox_quartile <- coxph(
  Surv(time, status) ~ score_quartile + age + sex + bmi + drink + hypertension,
  data = baseline
)

# Extract the results of Model 2
quartile_summary <- summary(cox_quartile)

# Construct the cese/n ratio
group_counts <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    cases = sum(status),
    total = n()
  )

# Extract HR, CI, and P-value of Model 2
quartile_hr <- exp(cox_quartile$coefficients[grep("score_quartile", names(cox_quartile$coefficients))])
quartile_ci <- exp(confint(cox_quartile)[grep("score_quartile", rownames(confint(cox_quartile))), ])
quartile_p <- quartile_summary$coefficients[grep("score_quartile", rownames(quartile_summary$coefficients)), "Pr(>|z|)"]

# Obtain the score range for each group
score_ranges <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    min_score = min(score),
    max_score = max(score)
  ) %>%
  mutate(score_range = paste0(round(min_score, 1), " - ", round(max_score, 1)))

### Process and Combine Results ###
# Create the result data frame for Model 1
continuous_result <- data.frame(
  Analysis_Typ = "Continuous",
  Comparison = "Full range",
  Score_Range = "Full range",
  HR_CI = paste0(
    round(continuous_hr, 2), 
    " (", 
    round(continuous_ci[1], 2), 
    "-", 
    round(continuous_ci[2], 2), 
    ")"
  ),
  P_Value = formatC(continuous_p, format = "e", digits = 2),
  Cases_TotalN = paste0(cases, "/", total_n)
)

# Create the result data frame for Model 2
quartile_results <- data.frame(
  Analysis_Typ = "Quartile",
  Comparison = c("Q1 (Ref)", paste0("Q", 2:4)),
  Score_Range = score_ranges$score_range,
  HR_CI = c(
    "Reference",
    paste0(
      round(quartile_hr, 2),
      " (",
      round(quartile_ci[, 1], 2),
      "-",
      round(quartile_ci[, 2], 2),
      ")"
    )
  ),
  P_Value = c("", formatC(quartile_p, format = "e", digits = 2)),
  Cases_TotalN = c(
    paste0(group_counts$cases[1], "/", group_counts$total[1]),
    paste0(group_counts$cases[2:4], "/", group_counts$total[2:4])
  )
)

# Integrate results from both models
final_results <- bind_rows(continuous_result, quartile_results)

### Final Output ###
# View the final integrated results table
print(final_results)

# Save the final results 
write.csv(final_results, "./output/hrs_mainanalysis_result.csv", row.names = FALSE)


### Proportional Hazards (PH) Assumption Tests ###
# Model 1 
ph_cont_score <- cox.zph(cox_continuous)

# Print PH test results
print(ph_cont_score)

# Save results
ph_cont_score_df <- as.data.frame(ph_cont_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)

write.csv(ph_cont_score_df, "./output/hrs_mainanalysis_ph_continuous.csv", row.names = FALSE)

# Model 2
ph_quart_score <- cox.zph(cox_quartile)

# Print PH test results
print(ph_quart_score)

# Save results
ph_quart_score_df <- as.data.frame(ph_quart_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)

write.csv(ph_quart_score_df, "./output/hrs_mainanalysis_ph_quartile.csv", row.names = FALSE)



#### 2. SENSITIVITY ANALYSES ####
#### 2.1 SENSITIVITY ANALYSIS A: Excluding events in the first 2 years of follow-up####
### Data Preparation ##
data <- read.csv("./input/data.csv",header = T)

# Exclude subjects who developed the condition within the past two years
baseline <- data[!(data$time <= 2 & data$status == 1), ]

# Change the variable name BrainVital8 to score
colnames(baseline)[colnames(baseline) == "BrainVital8"] <- "score"

### Cox Proportional-Hazards Models ###
### Model 1: Score as a continuous variable ###
cox_cont <- coxph(
  Surv(time, status) ~ score + age + sex + bmi + drink + hypertension,
  data = baseline
)

# Extract the results of Model 1
cont_sum <- summary(cox_cont)
cont_hr <- exp(cox_cont$coefficients["score"])
cont_ci <- exp(confint(cox_cont)["score", ])
cont_p <- cont_sum$coefficients["score", "Pr(>|z|)"]
total_n <- nrow(baseline)
cases <- sum(baseline$status)

continuous_result <- data.frame(
  Analysis_Typ = "Continuous",
  Comparison = "Full range",
  Score_Range = "Full range",
  HR_CI = paste0(round(cont_hr, 2), " (", round(cont_ci[1], 2), "-", round(cont_ci[2], 2), ")"),
  P_Value = formatC(cont_p, format = "e", digits = 2),
  Cases_TotalN = paste0(cases, "/", total_n)
)

### Model 2: Score as a categorical (quartile) variable ###
baseline<- baseline%>%
  mutate(score_quartile = cut(score,
                              breaks = quantile(score, probs = c(0, 0.25, 0.5, 0.75, 1)),
                              include.lowest = TRUE,
                              labels = c("Q1", "Q2", "Q3", "Q4")))

cox_quart <- coxph(
  Surv(time, status) ~ score_quartile + age + sex + bmi + drink + hypertension,
  data = baseline
)

# Extract the results of Model 2
quart_sum <- summary(cox_quart)

group_counts <- baseline%>%
  group_by(score_quartile) %>%
  summarise(cases = sum(status), total = n(), .groups = "drop")

quart_hr <- exp(cox_quart$coefficients[grep("score_quartile", names(cox_quart$coefficients))])
quart_ci <- exp(confint(cox_quart)[grep("score_quartile", rownames(confint(cox_quart))), ])
quart_p <- quart_sum$coefficients[grep("score_quartile", rownames(quart_sum$coefficients)), "Pr(>|z|)"]

score_ranges <- baseline%>%
  group_by(score_quartile) %>%
  summarise(min_score = min(score), max_score = max(score), .groups = "drop") %>%
  mutate(score_range = paste0(round(min_score, 1), " - ", round(max_score, 1)))

quartile_result <- data.frame(
  Analysis_Typ = "Quartile",
  Comparison = c("Q1 (Ref)", paste0("Q", 2:4)),
  Score_Range = score_ranges$score_range,
  HR_CI = c("Reference",
            paste0(round(quart_hr, 2), " (", round(quart_ci[, 1], 2), "-", round(quart_ci[, 2], 2), ")")),
  P_Value = c("", formatC(quart_p, format = "e", digits = 2)),
  Cases_TotalN = c(
    paste0(group_counts$cases[1], "/", group_counts$total[1]),
    paste0(group_counts$cases[2:4], "/", group_counts$total[2:4])
  )
)

### Process and Combine Results ###
final_results_filtered <- bind_rows(continuous_result, quartile_result)

# View the final integrated results table
print(final_results_filtered)

#  Final Output
write.csv(final_results_filtered, "./output/HRS_excl2yr_result.csv", row.names = FALSE)


### Proportional Hazards (PH) Assumption Tests ###
# Model 1
ph_cont_score <- cox.zph(cox_cont)

# Print PH test results
print(ph_cont_score)

# Save results
ph_cont_score_df <- as.data.frame(ph_cont_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)
write.csv(ph_cont_score_df, "./output/hrs_excl2yr_ph_continuous.csv", row.names = FALSE)

# Model 2
ph_quart_score <- cox.zph(cox_quart)

# Print PH test results
print(ph_quart_score)

# Save results
ph_quart_score_df <- as.data.frame(ph_quart_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)

write.csv(ph_quart_score_df, "./output/hrs_excl2yr_ph_quartile.csv", row.names = FALSE)



#### 2.2 SENSITIVITY ANALYSIS B: Stratified by Sex and Age Group ---------------------------------------------------------------------
### Data Preparation ##
baseline <- read.csv("./input/data.csv",header = TRUE)

# Change the variable name BrainVital8 to score
colnames(baseline)[colnames(baseline) == "BrainVital8"] <- "score"

### Cox Proportional-Hazards Models and Proportional Hazards (PH) Assumption Tests ###
run_cox_analysis <- function(data, strat_label, exclude_vars = NULL){
  
# File name sanitization: Replace non-alphanumeric characters in the stratification label with underscores.
  file_label <- gsub("[^A-Za-z0-9_]", "_", strat_label)
  
 # -- Model 1 --
  covariates <- c("score", "age", "sex", "bmi", "drink","hypertension")
  if(!is.null(exclude_vars)){
    covariates <- setdiff(covariates, exclude_vars)
  }
  
  formula_cont <- as.formula(paste0("Surv(time, status) ~ ", paste(covariates, collapse = " + ")))
  
  cox_cont <- coxph(formula_cont, data = data)
  cont_sum <- summary(cox_cont)
  
  # -- PH Assumption Tests of Model 1 --
  ph_cont <- cox.zph(cox_cont)
  
  # Convert to a data frame and ensure the column names match requirements (Variable, chisq, df, p)
  ph_cont_df <- as.data.frame(ph_cont$table) %>%
    mutate(Variable = rownames(.)) %>%
    select(Variable, chisq, df, p) 
  
  # Save results
  ph_cont_filename <- paste0("./output/PH_Cont_", file_label, ".csv")
  write.csv(ph_cont_df, ph_cont_filename, row.names = FALSE)
  cat(paste0("Continuous pH test results saved to: ", ph_cont_filename, "\n"))
  
  # -- Continuous HR Result Extraction --
  # Check whether ‘score’ is included in the model (it may not be present if “score” is listed in exclude_vars)
  if("score" %in% names(cox_cont$coefficients)){
    cont_hr <- exp(cox_cont$coefficients["score"])
    cont_ci <- exp(confint(cox_cont)["score", ])
    cont_p <- cont_sum$coefficients["score", "Pr(>|z|)"]
  } else {
    cont_hr <- NA; cont_ci <- c(NA, NA); cont_p <- NA
  }
  
  total_n <- nrow(data)
  cases <- sum(data$status)
  
  continuous_result <- data.frame(
    Stratification = strat_label,
    Analysis_Typ = "Continuous",
    Comparison = "Full range",
    Score_Range = "Full range",
    HR_CI = paste0(round(cont_hr, 2), " (", round(cont_ci[1], 2), "-", round(cont_ci[2], 2), ")"),
    P_Value = formatC(cont_p, format = "e", digits = 2),
    Cases_TotalN = paste0(cases, "/", total_n)
  )
  
  # -- Model 2 --
  data <- data %>%
    mutate(score_quartile = cut(score,
                                breaks = quantile(score, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                                include.lowest = TRUE,
                                labels = c("Q1", "Q2", "Q3", "Q4")))
  
  formula_quart <- as.formula(paste0("Surv(time, status) ~ score_quartile + ", 
                                     paste(setdiff(covariates, "score"), collapse = " + ")))
  
  cox_quart <- coxph(formula_quart, data = data)
  quart_sum <- summary(cox_quart)
  
  # -- PH Assumption Tests of Model 2 --
  ph_quart <- cox.zph(cox_quart)
  
  # Convert to a data frame and ensure the column names match requirements (Variable, chisq, df, p)
  ph_quart_df <- as.data.frame(ph_quart$table) %>%
    mutate(Variable = rownames(.)) %>%
    select(Variable, chisq, df, p)
  
  # Save results
  ph_quart_filename <- paste0("./output/PH_Quart_", file_label, ".csv")
  write.csv(ph_quart_df, ph_quart_filename, row.names = FALSE)
  cat(paste0("The results of the quartile pH test have been saved to: ", ph_quart_filename, "\n"))
  
  # -- Extraction of Quartile HR Results --
  group_counts <- data %>%
    group_by(score_quartile) %>%
    summarise(cases = sum(status), total = n(), .groups = "drop")
  
  quart_hr <- exp(cox_quart$coefficients[grep("score_quartile", names(cox_quart$coefficients))])
  quart_ci <- exp(confint(cox_quart)[grep("score_quartile", rownames(confint(cox_quart))), ])
  quart_p <- quart_sum$coefficients[grep("score_quartile", rownames(quart_sum$coefficients)), "Pr(>|z|)"]
  
  score_ranges <- data %>%
    group_by(score_quartile) %>%
    summarise(min_score = min(score), max_score = max(score), .groups = "drop") %>%
    mutate(score_range = paste0(round(min_score, 1), " - ", round(max_score, 1)))
  
  quartile_result <- data.frame(
    Stratification = strat_label,
    Analysis_Typ = "Quartile",
    Comparison = c("Q1 (Ref)", paste0("Q", 2:4)),
    Score_Range = score_ranges$score_range,
    HR_CI = c("Reference",
              paste0(round(quart_hr, 2), " (", round(quart_ci[, 1], 2), "-", round(quart_ci[, 2], 2), ")")),
    P_Value = c("", formatC(quart_p, format = "e", digits = 2)),
    Cases_TotalN = c(
      paste0(group_counts$cases[1], "/", group_counts$total[1]),
      paste0(group_counts$cases[2:4], "/", group_counts$total[2:4])
    )
  )
  
  # Return the merged HR results
  return(bind_rows(continuous_result, quartile_result))
}

# -- Main Program: Execute four hierarchical analyses and generate PH test files --
# Age < 65
data_age_lt65 <- subset(baseline, age < 65)
res_age_lt65 <- run_cox_analysis(data_age_lt65, "Age_LT65", exclude_vars = "age")

# Age ≥ 65
data_age_ge65 <- subset(baseline, age >= 65)
res_age_ge65 <- run_cox_analysis(data_age_ge65, "Age_GE65", exclude_vars = "age")

# Female
data_female <- subset(baseline, sex == 0)
res_female <- run_cox_analysis(data_female, "Female", exclude_vars = "sex")

# Male
data_male <- subset(baseline, sex == 1)
res_male <- run_cox_analysis(data_male, "Male", exclude_vars = "sex")

# Merge all hierarchical results
final_stratified_results <- bind_rows(res_age_lt65, res_age_ge65, res_female, res_male)

# Save the final merged HR/CI results
write.csv(final_stratified_results, "./output/hrs_age_gender_ph_result.csv", row.names = FALSE)


#### 2.3 SENSITIVITY ANALYSIS C: Additionally adjusted for LIBRA2 ---------------------------------------------------------------------
### Data Preparation ##
baseline <- read.csv("./input/data.csv",header = TRUE)

# Change the variable name BrainVital8 to score
colnames(baseline)[colnames(baseline) == "BrainVital8"] <- "score"

### Cox Proportional-Hazards Models ###
### Model 1: Score as a continuous variable ###
cox_continuous <- coxph(
  Surv(time, status) ~ score + age + sex + bmi + drink + hypertension + libra2,
  data = baseline
)

# Extract the results of Model 1
continuous_summary <- summary(cox_continuous)
continuous_hr <- exp(cox_continuous$coefficients["score"])
continuous_ci <- exp(confint(cox_continuous)["score", ])
continuous_p <- continuous_summary$coefficients["score", "Pr(>|z|)"]
total_n <- nrow( baseline)
cases <- sum( baseline$status)

### Model 2: Score as a categorical (quartile) variable ###
baseline <-  baseline %>%
  mutate(
    score_quartile = cut(score,
                         breaks = quantile(score, probs = c(0, 0.25, 0.5, 0.75, 1)),
                         include.lowest = TRUE,
                         labels = c("Q1", "Q2", "Q3", "Q4"))
  )

cox_quartile <- coxph(
  Surv(time, status) ~ score_quartile + age + sex + bmi + drink + hypertension + libra2,
  data = baseline
)

# Extract the results of Model 2
quartile_summary <- summary(cox_quartile)

# Construct the cese/n ratio
group_counts <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    cases = sum(status),
    total = n()
  )

# Extract HR, CI, and P-value
quartile_hr <- exp(cox_quartile$coefficients[grep("score_quartile", names(cox_quartile$coefficients))])
quartile_ci <- exp(confint(cox_quartile)[grep("score_quartile", rownames(confint(cox_quartile))), ])
quartile_p <- quartile_summary$coefficients[grep("score_quartile", rownames(quartile_summary$coefficients)), "Pr(>|z|)"]

# Obtain the score range for each group
score_ranges <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    min_score = min(score),
    max_score = max(score)
  ) %>%
  mutate(score_range = paste0(round(min_score, 1), " - ", round(max_score, 1)))

### Process and Combine Results ###
# Create the result data frame for Model 1
continuous_result <- data.frame(
  Analysis_Typ = "Continuous",
  Comparison = "Full range",
  Score_Range = "Full range",
  HR_CI = paste0(
    round(continuous_hr, 2), 
    " (", 
    round(continuous_ci[1], 2), 
    "-", 
    round(continuous_ci[2], 2), 
    ")"
  ),
  P_Value = formatC(continuous_p, format = "e", digits = 2),
  Cases_TotalN = paste0(cases, "/", total_n)
)

# Create the result data frame for Model 2
quartile_results <- data.frame(
  Analysis_Typ = "Quartile",
  Comparison = c("Q1 (Ref)", paste0("Q", 2:4)),
  Score_Range = score_ranges$score_range,
  HR_CI = c(
    "Reference",
    paste0(
      round(quartile_hr, 2),
      " (",
      round(quartile_ci[, 1], 2),
      "-",
      round(quartile_ci[, 2], 2),
      ")"
    )
  ),
  P_Value = c("", formatC(quartile_p, format = "e", digits = 2)),
  Cases_TotalN = c(
    paste0(group_counts$cases[1], "/", group_counts$total[1]),
    paste0(group_counts$cases[2:4], "/", group_counts$total[2:4])
  )
)

# Integrate results from both models
final_results <- bind_rows(continuous_result, quartile_results)

### Final Output ###
# View the final integrated results table
print(final_results)

# Save the final results table to a file (e.g., CSV)
write.csv(final_results, "./output/hrs_adj_libra2_result.csv", row.names = FALSE)


### Proportional Hazards (PH) Assumption Tests ###
# Model 1
ph_cont_score <- cox.zph(cox_continuous)

# Print PH test results
print(ph_cont_score)

# Save results
ph_cont_score_df <- as.data.frame(ph_cont_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)
write.csv(ph_cont_score_df, "./output/hrs_adj_libra2_ph_continuous.csv", row.names = FALSE)

# Model 2
ph_quart_score <- cox.zph(cox_quartile)

# Print PH test results
print(ph_quart_score)

# Save results
ph_quart_score_df <- as.data.frame(ph_quart_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)
write.csv(ph_quart_score_df, "./output/hrs_adj_libra2_ph_quartile.csv", row.names = FALSE)



#### 2.4 SENSITIVITY ANALYSIS D: Additionally adjusted for Lancet ---------------------------------------------------------------------
### Data Preparation ##
baseline <- read.csv("./input/data.csv",header = TRUE)

# Change the variable name BrainVital8 to score
colnames(baseline)[colnames(baseline) == "BrainVital8"] <- "score"

### Cox Proportional-Hazards Models ###
### Model 1: Score as a continuous variable ###
cox_continuous <- coxph(
  Surv(time, status) ~ score + age + sex + bmi + drink + hypertension + lancet,
  data = baseline
)

# Extract the results of Model 1
continuous_summary <- summary(cox_continuous)
continuous_hr <- exp(cox_continuous$coefficients["score"])
continuous_ci <- exp(confint(cox_continuous)["score", ])
continuous_p <- continuous_summary$coefficients["score", "Pr(>|z|)"]
total_n <- nrow( baseline)
cases <- sum( baseline$status)

### Model 2: Score as a categorical (quartile) variable ###
baseline <-  baseline %>%
  mutate(
    score_quartile = cut(score,
                         breaks = quantile(score, probs = c(0, 0.25, 0.5, 0.75, 1)),
                         include.lowest = TRUE,
                         labels = c("Q1", "Q2", "Q3", "Q4"))
  )

cox_quartile <- coxph(
  Surv(time, status) ~ score_quartile + age + sex + bmi + drink + hypertension + lancet,
  data = baseline
)

# Extract the results of Model 2
quartile_summary <- summary(cox_quartile)

# Construct the cese/n ratio
group_counts <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    cases = sum(status),
    total = n()
  )

# Extract HR, CI, and P-value
quartile_hr <- exp(cox_quartile$coefficients[grep("score_quartile", names(cox_quartile$coefficients))])
quartile_ci <- exp(confint(cox_quartile)[grep("score_quartile", rownames(confint(cox_quartile))), ])
quartile_p <- quartile_summary$coefficients[grep("score_quartile", rownames(quartile_summary$coefficients)), "Pr(>|z|)"]

# Obtain the score range for each group
score_ranges <- baseline %>%
  group_by(score_quartile) %>%
  summarise(
    min_score = min(score),
    max_score = max(score)
  ) %>%
  mutate(score_range = paste0(round(min_score, 1), " - ", round(max_score, 1)))

### Process and Combine Results ###
# Create the result data frame for Model 1
continuous_result <- data.frame(
  Analysis_Typ = "Continuous",
  Comparison = "Full range",
  Score_Range = "Full range",
  HR_CI = paste0(
    round(continuous_hr, 2), 
    " (", 
    round(continuous_ci[1], 2), 
    "-", 
    round(continuous_ci[2], 2), 
    ")"
  ),
  P_Value = formatC(continuous_p, format = "e", digits = 2),
  Cases_TotalN = paste0(cases, "/", total_n)
)

# Create the result data frame for Model 2
quartile_results <- data.frame(
  Analysis_Typ = "Quartile",
  Comparison = c("Q1 (Ref)", paste0("Q", 2:4)),
  Score_Range = score_ranges$score_range,
  HR_CI = c(
    "Reference",
    paste0(
      round(quartile_hr, 2),
      " (",
      round(quartile_ci[, 1], 2),
      "-",
      round(quartile_ci[, 2], 2),
      ")"
    )
  ),
  P_Value = c("", formatC(quartile_p, format = "e", digits = 2)),
  Cases_TotalN = c(
    paste0(group_counts$cases[1], "/", group_counts$total[1]),
    paste0(group_counts$cases[2:4], "/", group_counts$total[2:4])
  )
)

# Integrate results from both models
final_results <- bind_rows(continuous_result, quartile_results)

### Final Output ###
print(final_results)

# View the final integrated results table
write.csv(final_results, "./output/hrs_adj_lancet_result.csv", row.names = FALSE)


### Proportional Hazards (PH) Assumption Tests ###
# Model 1 
ph_cont_score <- cox.zph(cox_continuous)

# Print PH test results
print(ph_cont_score)

# Save results
ph_cont_score_df <- as.data.frame(ph_cont_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)
write.csv(ph_cont_score_df, "./output/hrs_adj_lancet_ph_continuous.csv", row.names = FALSE)

# Model 2
ph_quart_score <- cox.zph(cox_quartile)

# Print PH test results
print(ph_quart_score)

# Save results
ph_quart_score_df <- as.data.frame(ph_quart_score$table) %>%
  mutate(Variable = rownames(.)) %>%
  select(Variable, chisq, df, `p`)
write.csv(ph_quart_score_df, "./output/hrs_adj_lancet_ph.csv", row.names = FALSE)

