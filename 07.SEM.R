library(lavaan)
library(dplyr)

# ============================================================
# 1) Load SHARE data
# ============================================================
file_path <- "./input/SHARE_original_data.csv"
dat_raw <- read.csv(file_path, check.names = FALSE, fileEncoding = "UTF-8")
dat <- dat_raw
names(dat) <- make.names(names(dat), unique = TRUE)

# ============================================================
# 2) Robust column getter (Compatible column names)
# ============================================================
norm_name <- function(x) gsub("[\\._\\s]+", "", tolower(x))

get_col2 <- function(df, target, allow_partial = TRUE) {
  nn_raw <- names(df)
  nn <- norm_name(nn_raw)
  tt <- norm_name(target)
  
  idx <- which(nn == tt)
  if (length(idx) == 1) return(df[[idx]])
  if (length(idx) > 1) {
    message("⚠️ Multiple exact matches, take the first one:", nn_raw[idx[1]])
    return(df[[idx[1]]])
  }
  
  if (allow_partial) {
    idx2 <- which(grepl(tt, nn, fixed = TRUE) | startsWith(nn, tt) | startsWith(tt, nn))
    if (length(idx2) == 1) {
      message("✅ using fuzzy matching:", target, " -> ", nn_raw[idx2])
      return(df[[idx2]])
    }
    if (length(idx2) > 1) {
      message("⚠️ Fuzzy matches to multiple columns:", target, " -> ",
              paste(nn_raw[idx2], collapse = ", "),
              "\n take the first one:", nn_raw[idx2[1]])
      return(df[[idx2[1]]])
    }
  }
  
  cand <- nn_raw[grepl(substr(tt, 1, min(5, nchar(tt))), nn, fixed = TRUE)]
  stop(paste0(
    "❌ Column not found:,",target,
    "\n Possible candidates:", ifelse(length(cand) == 0, "(No obvious candidate)", paste(cand, collapse = ", ")),
    "\n Check the true column name of names(dat)."
  ))
}

# ============================================================
# 3) Helpers
# ============================================================
z <- function(x) as.numeric(scale(x))

to_num <- function(x) {
  if (is.factor(x) || is.character(x)) return(as.numeric(as.character(x)))
  as.numeric(x)
}

ph010_sleep_fix <- function(x) {
  v <- to_num(x)
  v[v %in% c(-2, -1)] <- 0
  v
}

# 0/1 problem indicator -> health = 1 - v
sel_problem_health01 <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         ifelse(v %in% c(0,1), 1 - v, NA_real_))
}

# smoke br002_: 1 yes (worse), 5 no (better)
smoke_health01 <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         ifelse(v == 5, 1,
                ifelse(v == 1, 0, NA_real_)))
}

house_health01 <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         ifelse(v == 1, 1,
                ifelse(v %in% c(2,3,4,5), 0, NA_real_)))
}

marry_health01 <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         ifelse(v %in% c(1,2), 1,
                ifelse(v %in% c(3,4,5,6), 0, NA_real_)))
}

edu_health <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_, v)
}

job_health_proxy <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         dplyr::case_when(
           v == 2 ~ 1,
           v == 1 ~ 0.5,
           v %in% c(3,4,5,97) ~ 0,
           TRUE ~ NA_real_
         ))
}

# EURO: symptom item assumed 0/1 (1 symptom present worse) -> health=1-x
euro_symptom_health01 <- function(x) sel_problem_health01(x)
# euro2 hopefulness assumed 0/1 (1 hopeful better) -> health=x
euro_positive_health01 <- function(x) {
  v <- to_num(x)
  ifelse(is.na(v), NA_real_,
         ifelse(v %in% c(0,1), v, NA_real_))
}

# ============================================================
# 4) Build domains (higher = healthier)
# ============================================================
dat2 <- dat %>%
  mutate(
    mergeid = get_col2(dat, "mergeid"),
    
    # Sleep
    ph010d6_fix = ph010_sleep_fix(get_col2(dat, "ph010d6")),
    sleep_prob_h = z(sel_problem_health01(ph010d6_fix)),
    sleep_med_h  = z(sel_problem_health01(get_col2(dat, "ph011d9"))),
    
    # Diabetes
    diab_dx_h   = z(sel_problem_health01(get_col2(dat, "ph006d5"))),
    diab_drug_h = z(sel_problem_health01(get_col2(dat, "ph011d6"))),
    a1c_h       = -z(to_num(get_col2(dat, "a1c_x4"))),
    
    # Grip
    gs006_h = z(to_num(get_col2(dat, "gs006_"))),
    gs007_h = z(to_num(get_col2(dat, "gs007_"))),
    gs008_h = z(to_num(get_col2(dat, "gs008_"))),
    gs009_h = z(to_num(get_col2(dat, "gs009_"))),
    
    # Activity
    activity_h = z(to_num(get_col2(dat, "PA_score"))),
    
    # SDOH
    house_h   = z(house_health01(get_col2(dat, "ho002_"))),
    married_h = z(marry_health01(get_col2(dat, "dn014_"))),
    edu_h     = z(edu_health(get_col2(dat, "raeducl"))),
    
    income_raw = to_num(get_col2(dat, "Hhitothhinc")),
    income_q   = ifelse(is.na(income_raw), NA_real_, ntile(income_raw, 4)),
    income_h   = z(income_q),
    
    job_h = z(job_health_proxy(get_col2(dat, "ep005_"))),
    
    # ✅ Breath (Fixed to puff1/puff2)
    puff1_h = z(to_num(get_col2(dat, "puff1"))),
    puff2_h = z(to_num(get_col2(dat, "puff2"))),
    
    # EURO 1-12
    euro1_h  = z(euro_symptom_health01(get_col2(dat, "euro1"))),
    euro2_h  = z(euro_positive_health01(get_col2(dat, "euro2"))),
    euro3_h  = z(euro_symptom_health01(get_col2(dat, "euro3"))),
    euro4_h  = z(euro_symptom_health01(get_col2(dat, "euro4"))),
    euro5_h  = z(euro_symptom_health01(get_col2(dat, "euro5"))),
    euro6_h  = z(euro_symptom_health01(get_col2(dat, "euro6"))),
    euro7_h  = z(euro_symptom_health01(get_col2(dat, "euro7"))),
    euro8_h  = z(euro_symptom_health01(get_col2(dat, "euro8"))),
    euro9_h  = z(euro_symptom_health01(get_col2(dat, "euro9"))),
    euro10_h = z(euro_symptom_health01(get_col2(dat, "euro10"))),
    euro11_h = z(euro_symptom_health01(get_col2(dat, "euro11"))),
    euro12_h = z(euro_symptom_health01(get_col2(dat, "euro12"))),
    
    # Smoking
    smoke_h = z(smoke_health01(get_col2(dat, "br002_")))
  )

dat3 <- dat2 %>%
  mutate(
    sleep_comp   = rowMeans(cbind(sleep_prob_h, sleep_med_h), na.rm = TRUE),
    glucose_comp = rowMeans(cbind(diab_dx_h, diab_drug_h, a1c_h), na.rm = TRUE),
    grip_comp    = rowMeans(cbind(gs006_h, gs007_h, gs008_h, gs009_h), na.rm = TRUE),
    activity_comp = activity_h,
    SDOH_comp    = rowMeans(cbind(house_h, married_h, edu_h, income_h, job_h), na.rm = TRUE),
    
    # ✅ respiratory_comp  to puff1/puff2
    respiratory_comp = rowMeans(cbind(puff1_h, puff2_h), na.rm = TRUE),
    
    emotional_comp = rowMeans(
      cbind(euro1_h, euro2_h, euro3_h, euro4_h, euro5_h, euro6_h,
            euro7_h, euro8_h, euro9_h, euro10_h, euro11_h, euro12_h),
      na.rm = TRUE
    ),
    smoke_comp = smoke_h
  )

# ============================================================
# 5) SEM
# ============================================================
model_share <- '
  BrainVital8 =~ respiratory_comp + glucose_comp + sleep_comp +
                 grip_comp + activity_comp + SDOH_comp +
                 emotional_comp + smoke_comp

  grip_comp ~~ activity_comp
  respiratory_comp ~~ glucose_comp
  sleep_comp ~~ emotional_comp
'

fit_share <- cfa(
  model_share,
  data      = dat3,
  estimator = "MLR",
  missing   = "fiml",
  std.lv    = TRUE,
  control   = list(iter.max = 5000)
)

summary(fit_share, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
fitMeasures(fit_share, c("cfi","tli","rmsea","srmr","chisq","df","pvalue"))

# ============================================================
# 6) Export for alignment
# ============================================================
scores <- as.data.frame(lavPredict(fit_share, type = "lv"))

export_df <- dat3 %>%
  transmute(
    mergeid,
    respiratory_comp, glucose_comp, sleep_comp, grip_comp,
    activity_comp, SDOH_comp, emotional_comp, smoke_comp
  ) %>%
  bind_cols(scores)

out_path <- "./output/SHARE_BrainVital8_domains_for_alignment.csv"
write.csv(export_df, out_path, row.names = FALSE, fileEncoding = "UTF-8")
cat("\n✅ exported：", out_path, "\n")

# ============================================================
# 7) Export loadings & residual variances
# ============================================================
pe <- parameterEstimates(fit_share, standardized = TRUE)

loadings <- pe %>% filter(op == "=~") %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.all)
write.csv(loadings,
          "./output/SHARE_factor_loadings.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

man_vars <- lavNames(fit_share, type = "ov")
resid_var <- pe %>% filter(op == "~~", lhs == rhs, lhs %in% man_vars) %>%
  select(lhs, op, rhs, est, se, z, pvalue, std.all)
write.csv(resid_var,
          "./output/SHARE_residual_variances.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
