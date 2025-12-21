library(dplyr)
library(survival)

# Load data
BB <- read.csv("./input/brainvital8.csv")

# Check missing values in each variable
missing_count <- colSums(is.na(BB))
print(missing_count)

# View column names
names(BB)

# Create analysis function
run_analysis <- function(data, analysis_name, adjust_sex = TRUE) {
  # Calculate quartiles within current data
  data <- data %>%
    mutate(
      total_score_group_q = cut(
        BrainVital8,
        breaks = quantile(BrainVital8, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
        labels = c("Q1", "Q2", "Q3", "Q4"),
        include.lowest = TRUE
      )
    )
  
  # Calculate group statistics
  group_stats <- data %>%
    group_by(total_score_group_q) %>%
    summarise(
      N = n(),
      Cases = sum(status),
      .groups = "drop"
    )
  
  total_n <- nrow(data)
  total_cases <- sum(data$status)
  
  # Define covariates (does not include APOEe4_carrier)
  if(adjust_sex) {
    covariates_formula <- "age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression"
  } else {
    covariates_formula <- "age + bmi + smoke + drink + TDI + education + T2D + hypertension + depression"
  }
  
  # Continuous variable Cox regression
  formula_cont <- as.formula(paste("Surv(time, status) ~ BrainVital8 +", covariates_formula))
  cox_cont <- coxph(formula_cont, data = data)
  cont_summary <- summary(cox_cont)
  
  # Quartile group Cox regression
  formula_q <- as.formula(paste("Surv(time, status) ~ total_score_group_q +", covariates_formula))
  cox_q <- coxph(formula_q, data = data)
  q_summary <- summary(cox_q)
  
  # Extract results
  results_list <- list()
  
  # Continuous variable results
  if("BrainVital8" %in% rownames(cont_summary$coefficients)) {
    results_list[[1]] <- data.frame(
      Analysis = analysis_name,
      Model = "Continuous",
      Exposure = "BrainVital8",
      Comparison = "Per unit increase",
      N = total_n,
      Cases = total_cases,
      HR = exp(cont_summary$coefficients["BrainVital8", "coef"]),
      HR_lower = exp(confint(cox_cont)["BrainVital8", 1]),
      HR_upper = exp(confint(cox_cont)["BrainVital8", 2]),
      p_value = cont_summary$coefficients["BrainVital8", "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
  }
  
  # Quartile group results
  # Q1 reference group
  q1_stats <- group_stats %>% filter(total_score_group_q == "Q1")
  results_list[[2]] <- data.frame(
    Analysis = analysis_name,
    Model = "Quartile",
    Exposure = "BrainVital8",
    Comparison = "Q1 (Reference)",
    N = q1_stats$N,
    Cases = q1_stats$Cases,
    HR = 1.000,
    HR_lower = NA,
    HR_upper = NA,
    p_value = NA,
    stringsAsFactors = FALSE
  )
  
  # Q2-Q4 vs Q1
  for(i in 2:4) {
    term_name <- paste0("total_score_group_qQ", i)
    comp_name <- paste0("Q", i, " vs Q1")
    q_stats <- group_stats %>% filter(total_score_group_q == paste0("Q", i))
    
    if(term_name %in% rownames(q_summary$coefficients)) {
      results_list[[i+1]] <- data.frame(
        Analysis = analysis_name,
        Model = "Quartile",
        Exposure = "BrainVital8",
        Comparison = comp_name,
        N = q_stats$N,
        Cases = q_stats$Cases,
        HR = exp(q_summary$coefficients[term_name, "coef"]),
        HR_lower = exp(confint(cox_q)[term_name, 1]),
        HR_upper = exp(confint(cox_q)[term_name, 2]),
        p_value = q_summary$coefficients[term_name, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(bind_rows(results_list))
}

# Analysis 1: Main analysis in full population
cat("1. Main analysis in full population...\n")
main_results <- run_analysis(BB, "Main_Analysis", adjust_sex = TRUE)

# Analysis 2: Male subgroup analysis (do not adjust for sex)
cat("\n2. Male subgroup analysis...\n")
cat("Distribution of sex variable:\n")
print(table(BB$sex))

# Assuming sex = 1 is male
male_data <- BB %>% filter(sex == 1)
cat("Male sample size:", nrow(male_data), "Number of events:", sum(male_data$status), "\n")
male_results <- run_analysis(male_data, "Male_Subgroup", adjust_sex = FALSE)

# Analysis 3: Female subgroup analysis (do not adjust for sex)
cat("\n3. Female subgroup analysis...\n")
# Assuming sex = 0 is female
female_data <- BB %>% filter(sex == 0)
cat("Female sample size:", nrow(female_data), "Number of events:", sum(female_data$status), "\n")
female_results <- run_analysis(female_data, "Female_Subgroup", adjust_sex = FALSE)
                               
# Analysis 4: Exclude cases diagnosed within 2 years
cat("\n4. Excluding cases diagnosed within 2 years...\n")
exclude_2year_data <- BB %>%
  filter(!(status == 1 & time <= 730))
cat("Sample size after excluding 2-year cases:", nrow(exclude_2year_data), "\n")
cat("Number of excluded cases:", sum(BB$status == 1 & BB$time <= 730), "\n")
exclude_2year_results <- run_analysis(exclude_2year_data, "Exclude_2year", adjust_sex = TRUE)

# Analysis 5: APOEe4 = 0 population analysis
cat("\n5. APOEe4 = 0 population analysis...\n")
APOE_0_data <- BB %>% filter(APOEe4_carrier == 0)
cat("APOEe4 = 0 sample size:", nrow(APOE_0_data), "Number of events:", sum(APOE_0_data$status), "\n")
APOE_0_results <- run_analysis(APOE_0_data, "APOE_0_Only", adjust_sex = TRUE)

# Analysis 6: APOEe4 = 1 population analysis
cat("\n6. APOEe4 = 1 population analysis...\n")
APOE_1_data <- BB %>% filter(APOEe4_carrier == 1)
cat("APOEe4 = 1 sample size:", nrow(APOE_1_data), "Number of events:", sum(APOE_1_data$status), "\n")
APOE_1_results <- run_analysis(APOE_1_data, "APOE_1_Only", adjust_sex = TRUE)

# Combine all results
all_results <- bind_rows(
  main_results,
  male_results,
  female_results,
  exclude_2year_results,
  APOE_0_results,
  APOE_1_results
)

# Format results
final_results <- all_results %>%
  mutate(
    N_Cases = paste0(Cases, "/", N),
    HR_95CI = ifelse(is.na(p_value), "1.000 (Reference)",
                     sprintf("%.2f (%.2f-%.2f)", HR, HR_lower, HR_upper)),
    p_value_scientific = ifelse(is.na(p_value), "", 
                                format(p_value, scientific = TRUE, digits = 3)),
    p_value_raw = p_value
  ) %>%
  select(Analysis, Model, Exposure, Comparison, N_Cases, HR_95CI, p_value_raw, p_value_scientific)


# Export to CSV
write.csv(final_results, "./output/brainvital8_dementia_sensitivity.csv", row.names = FALSE)
cat("Results saved to: brainvital8_dementia_sensitivity.csv\n")





## Genetic Risk Score Analysis


# Load data
BB <- read.csv("./input/brainvital8.csv")
B <- read.csv("./input/FinnGen.csv")
B <- read.csv("./input/MVP.csv")
# Convert eid to character type for proper merging
B$eid <- as.character(B$eid)
BB$eid <- as.character(BB$eid)

# Load dplyr package
library(dplyr)

# Merge datasets by eid
A <- left_join(BB, B, by = "eid")

# Check missing values in each variable
missing_count <- colSums(is.na(A))
print(missing_count)

# Remove rows with any missing values
A <- na.omit(A)

# View column names
names(A)

# Categorize genetic risk score into tertiles (T1, T2, T3)
A$genetic_tertile <- cut(A$SCORESUM,
                         breaks = quantile(A$SCORESUM, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                         labels = c("T1", "T2", "T3"),
                         include.lowest = TRUE)

# Categorize BrainVital8 score into quartiles (Q1, Q2, Q3, Q4)
A$score_quartile <- cut(A$BrainVital8,
                        breaks = quantile(A$BrainVital8, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE),
                        labels = c("Q1", "Q2", "Q3", "Q4"),
                        include.lowest = TRUE)

# Create combined interaction variable with "Q1.T3" as reference
A$combined_group <- interaction(A$score_quartile, A$genetic_tertile)
A$combined_group <- relevel(A$combined_group, ref = "Q1.T3")

# Check distribution of combined groups
cat("Distribution of combined groups:\n")
table_data <- with(A, table(score_quartile, genetic_tertile))
print(table_data)

# Calculate detailed statistics for each combined group
cat("\nDetailed statistics for each combined group:\n")
combined_stats <- A %>%
  group_by(score_quartile, genetic_tertile) %>%
  summarise(
    Total = n(),
    Cases = sum(status),
    Incidence = paste0(round(mean(status) * 100, 2), "%")
  ) %>%
  arrange(genetic_tertile, score_quartile)

print(combined_stats)

# Run Cox regression analysis with combined groups
cox_combined <- coxph(Surv(time, status) ~ combined_group +
                        age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                      data = A)

cat("\n=== Cox regression results (reference: Q1.T3) ===\n")
summary(cox_combined)

# Extract detailed Cox regression results
cox_summary <- summary(cox_combined)

# Create results dataframe
results_df <- data.frame(
  Variable = names(cox_combined$coefficients),
  Coef = cox_combined$coefficients,
  HR = exp(cox_combined$coefficients),
  SE = cox_summary$coefficients[, "se(coef)"],
  Z = cox_summary$coefficients[, "z"],
  P_value = cox_summary$coefficients[, "Pr(>|z|)"],
  CI_lower = exp(cox_combined$coefficients - 1.96 * cox_summary$coefficients[, "se(coef)"]),
  CI_upper = exp(cox_combined$coefficients + 1.96 * cox_summary$coefficients[, "se(coef)"])
)

# Add formatted 95% CI column for readability
results_df$CI_95 <- paste0(round(results_df$CI_lower, 2), "-", round(results_df$CI_upper, 2))

# Select required columns
final_results <- results_df[, c("Variable", "HR", "CI_95", "P_value", "Coef", "SE", "Z")]

# Export to CSV
write.csv(final_results, "./output/combined_snp3_and_score.csv", row.names = FALSE)

# Method 1: Calculate trend using numerical scores
# Assign numerical scores representing "distance from reference group"
A$trend_score <- as.numeric(A$score_quartile) + as.numeric(A$genetic_tertile)

# Run Cox model for trend test
cox_trend <- coxph(Surv(time, status) ~ trend_score +
                     age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                   data = A)

# Extract and format trend p-value
trend_p <- summary(cox_trend)$coefficients["trend_score", "Pr(>|z|)"]
cat("\nP for trend:", format(trend_p, scientific = TRUE, digits = 3), "\n")

# Stratified analysis: Analyze lifestyle score effect within each genetic risk stratum
library(survival)

# Analysis within T1 genetic risk stratum
subset_T1 <- A %>% filter(genetic_tertile == "T1")
cox_T1 <- coxph(Surv(time, status) ~ score_quartile + age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                data = subset_T1)

# Analysis within T2 genetic risk stratum
subset_T2 <- A %>% filter(genetic_tertile == "T2")
cox_T2 <- coxph(Surv(time, status) ~ score_quartile + age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                data = subset_T2)

# Analysis within T3 genetic risk stratum
subset_T3 <- A %>% filter(genetic_tertile == "T3")
cox_T3 <- coxph(Surv(time, status) ~ score_quartile + age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                data = subset_T3)

# Output stratified analysis results
cat("\n=== Effect of lifestyle score within genetic risk T1 stratum ===\n")
cat("Sample size:", nrow(subset_T1), "Cases:", sum(subset_T1$status), "\n")
summary(cox_T1)

cat("\n=== Effect of lifestyle score within genetic risk T2 stratum ===\n")
cat("Sample size:", nrow(subset_T2), "Cases:", sum(subset_T2$status), "\n")
summary(cox_T2)

cat("\n=== Effect of lifestyle score within genetic risk T3 stratum ===\n")
cat("Sample size:", nrow(subset_T3), "Cases:", sum(subset_T3$status), "\n")
summary(cox_T3)

# Simplified version - only output basic statistics
cat("\n=== Basic statistics for each stratum ===\n")

for(layer in c("T1", "T2", "T3")) {
  subset_data <- get(paste0("subset_", layer))
  n <- nrow(subset_data)
  cases <- sum(subset_data$status)
  incidence <- round(cases / n * 100, 3)
  
  cat(layer, "stratum: N =", n, "Cases =", cases, "Incidence =", incidence, "%\n")
}

cat("\n=== Basic statistics for each group within strata ===\n")
basic_stats <- A %>%
  group_by(genetic_tertile, score_quartile) %>%
  summarise(
    N = n(),
    Cases = sum(status),
    Incidence = paste0(round(mean(status) * 100, 3), "%")
  ) %>%
  arrange(genetic_tertile, score_quartile)

print(basic_stats)

# Re-extract stratum results ensuring p-values match summary output exactly
strata_results <- data.frame()

# Precise extraction function - ensures p-values match summary exactly
extract_strata_results_exact <- function(cox_model, strata_name, data_subset) {
  model_summary <- summary(cox_model)
  coef_names <- names(cox_model$coefficients)
  score_terms <- coef_names[grepl("score_quartile", coef_names)]
  
  for(term in score_terms) {
    hr <- exp(cox_model$coefficients[term])
    ci <- exp(confint(cox_model)[term, ])
    p <- model_summary$coefficients[term, "Pr(>|z|)"]  # Extract directly from summary
    
    group <- gsub("score_quartile", "", term)
    strata_results <<- rbind(strata_results, data.frame(
      Strata = strata_name,
      Group = group,
      HR = round(hr, 2),
      CI_lower = round(ci[1], 4),
      CI_upper = round(ci[2], 4),
      CI_95 = paste0(round(ci[1], 4), "-", round(ci[2], 4)),
      P_value = p,  # Raw p-value, no rounding
      N = nrow(data_subset),
      Cases = sum(data_subset$status)
    ))
  }
}

# Re-extract results
extract_strata_results_exact(cox_T1, "T1", subset_T1)
extract_strata_results_exact(cox_T2, "T2", subset_T2)
extract_strata_results_exact(cox_T3, "T3", subset_T3)

# Add reference groups
strata_results <- rbind(
  data.frame(
    Strata = c("T1", "T2", "T3"),
    Group = "Q1 (Ref)",
    HR = 1.000,
    CI_lower = NA,
    CI_upper = NA,
    CI_95 = "Reference",
    P_value = NA,
    N = c(nrow(subset_T1), nrow(subset_T2), nrow(subset_T3)),
    Cases = c(sum(subset_T1$status), sum(subset_T2$status), sum(subset_T3$status))
  ),
  strata_results
)

# Display final table
cat("\n=== Final table for stratified analysis ===\n")
print(strata_results)

# Save results
write.csv(strata_results, "./output/snp3_stratified.csv", row.names = FALSE)

# Check column names
names(A)

# Method 1: Calculate interaction p-value using likelihood ratio test
# First, build additive model (without interaction term)
model_additive <- coxph(Surv(time, status) ~ score_quartile + genetic_tertile + 
                          age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                        data = A)

# Then build interaction model (with interaction term)
model_interaction <- coxph(Surv(time, status) ~ score_quartile * genetic_tertile + 
                             age + sex + bmi + smoke + drink + TDI + education + T2D + hypertension + depression, 
                           data = A)

# Compare models using likelihood ratio test
lr_test <- anova(model_additive, model_interaction)
cat("\n=== Interaction test (Likelihood Ratio Test) ===\n")
print(lr_test)


