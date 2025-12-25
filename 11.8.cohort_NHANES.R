library(tidyverse)
library(nhanesR)

# =========================
# NHANES analysis
# =========================

# Load NHANES dataset with additional adjustment scores
d <- read_rds("./input/NHANES_BrainVital8_alignedQ.rds")

# Create quartiles of BrainVital8 for visualization and categorical analysis
d$BrainVital8Q <- quant(d$BrainVital8, n = 4, Q = TRUE, round = 3)

# Preserve original dataset
data <- d

# Check sex distribution
table(d$sex)

# Stratified analyses by sex and age group
d <- select_row(data, data$sex == "1")   # males
d <- select_row(data, data$sex == "2")   # females
d <- select_row(data, data$age >= 65)    # age ≥ 65
d <- select_row(data, data$age < 65)     # age < 65

# =========================
# Main logistic regression models
# =========================

# Model using BrainVital8 quartiles
model <- glm(
  high_risk ~ BrainVital8_alignedQ + age + sex + BMI +
    alcohol.user + Hypertension,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

# Model using continuous BrainVital8 score
model <- glm(
  high_risk ~ BrainVital8_aligned + age + sex + BMI +
    alcohol.user + Hypertension,
  data = d,
  family = binomial
) %>% reg_table(round = 3)

# =========================
# Sample size summaries
# =========================

table(d$BrainVital8_alignedQ, d$high_risk)
table(d$BrainVital8_alignedQ)
table(d$high_risk)

# =========================
# Sensitivity analyses
# =========================

# --- Additional adjustment for LIBRA2 score ---
model <- glm(
  high_risk ~ BrainVital8_alignedQ+ age + sex + BMI +
    alcohol.user + Hypertension +
    LIBRA2_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

model <- glm(
  high_risk ~ BrainVital8_aligned+ age + sex + BMI +
    alcohol.user + Hypertension+
    LIBRA2_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 3)

# --- Additional adjustment for Lancet risk score ---
model <- glm(
  high_risk ~ BrainVital8_alignedQ+ age + sex + BMI +
    alcohol.user + Hypertension+
    lancet_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 2)

model <- glm(
  high_risk ~ BrainVital8_aligned+ age + sex + BMI +
    alcohol.user + Hypertension +
    lancet_scores,
  data = d,
  family = binomial
) %>% reg_table(round = 3)
