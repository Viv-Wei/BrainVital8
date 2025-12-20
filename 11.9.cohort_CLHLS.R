library(tidyverse)
library(nhanesR)


# =========================
# CLHLS analysis
# =========================

# Load CLHLS dataset
d <- read_rds("./input/CLHLS_two_scores.rds")

# Create BrainVital8 quartiles for visualization
d$BrainVital8Q <- quant(d$BrainVital8, n = 4, Q = TRUE, round = 3)

# Preserve original dataset
data <- d

# Stratified analyses
d <- select_row(data, data$sex == "male")
d <- select_row(data, data$sex == "female")
d <- select_row(data, data$age >= 65)

# Recreate quartiles after subsetting
d$BrainVital8Q <- quant(d$BrainVital8, n = 4, Q = TRUE, round = 3)

# =========================
# Main models (CLHLS)
# =========================

# Model using BrainVital8 quartiles
model <- glm(
  adl_cog_dual_0826 ~ BrainVital8Q + age + sex + bmi + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

# Model using continuous BrainVital8 score
model <- glm(
  adl_cog_dual_0826 ~ BrainVital8 + age + sex + bmi + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n,
  data = d,
  family = binomial
) %>% reg_table(round = 3)

# =========================
# Sample size summaries
# =========================

table(d$BrainVital8Q, d$adl_cog_dual_0826)
table(d$BrainVital8Q)
table(d$adl_cog_dual_0826)
table(d$cesd10_n)

# =========================
# Sensitivity analyses (additional risk scores)
# =========================

# --- Adjustment for LIBRA2 ---
model <- glm(
  adl_cog_dual_0826 ~ BrainVital8Q + age + sex + bmiQ + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n + LIBRA2Q,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

model <- glm(
  adl_cog_dual_0826 ~ BrainVital8 + age + sex + bmi + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n + LIBRA2Q,
  data = d,
  family = binomial
) %>% reg_table(round = 3)

# --- Adjustment for Lancet risk score ---
model <- glm(
  adl_cog_dual_0826 ~ BrainVital8Q + age + sex + bmiQ + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n + lancet_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

model <- glm(
  adl_cog_dual_0826 ~ BrainVital8 + age + sex + bmi + smoking + drinking +
    years_schooling + total_income + suffer_diabetes +
    hypertension + cesd10_n + lancet_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 3)
