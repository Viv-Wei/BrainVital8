############################## Dementia Cox analysis ##############################
############################## Data input #########################################

# Load dementia outcomes and harmonize participant identifiers
allcause_dementia <- read.csv("./input/allcause_dementia.csv") |> mutate(eid = as.character(eid))
F00 <- read.csv("F00.csv") |> mutate(eid = as.character(eid))
F01 <- read.csv("F01.csv") |> mutate(eid = as.character(eid))
G30 <- read.csv("G30.csv") |> mutate(eid = as.character(eid))

# Load demographic, lifestyle and genetic covariates
cov1 <- read.csv("./input/M45_M51+cov.csv", header = T) |> 
  mutate(eid = as.character(eid)) |> 
  select(eid, age, sex, race, bmi, smoke, drink, APOEe4_carrier) 

# Recode covariates to harmonized categorical representations
{
  cov1 <- cov1 %>%
    mutate(race = case_when(
      race == "-3" ~ "9",
      race == "-1" ~ "9",
      race %in% c("1","2","3","4","5","6") ~ "2",
      race == "1001" ~ "1",
      race %in% c("1002","1003","2001","2002","2003","2004",
                  "3001","3002","3003","3004",
                  "4001","4002","4003") ~ "2",
      is.na(race) ~ "9",
      TRUE ~ as.character(race)
    ))
  
  cov1 <- cov1 %>%
    mutate(smoke = case_when(
      smoke == "-3" ~ "9",
      smoke %in% c("0","1","2") ~ smoke,
      is.na(smoke) ~ "9",
      TRUE ~ as.character(smoke)
    ))
  
  cov1 <- cov1 %>%
    mutate(drink = case_when(
      drink == "-3" ~ "9",
      drink %in% c("0","1","2") ~ drink,
      is.na(drink) ~ "9",
      TRUE ~ as.character(drink)
    ))
  
  cov1 <- cov1 %>%
    mutate(APOEe4_carrier = case_when(
      APOEe4_carrier %in% c("0","1","9") ~ APOEe4_carrier,
      is.na(APOEe4_carrier) ~ "9",
      TRUE ~ as.character(APOEe4_carrier)
    ))
  }

# Convert covariates to factor variables
cov1$race <- as.factor(cov1$race)
cov1$sex <- as.factor(cov1$sex)
cov1$smoke <- as.factor(cov1$smoke)
cov1$drink <- as.factor(cov1$drink)
cov1$APOEe4_carrier <- as.factor(cov1$APOEe4_carrier)

# Merge covariates with dementia outcomes
mydata   <- left_join(allcause_dementia, cov1, by = "eid")
mydata00 <- left_join(F00, cov1, by = "eid")
mydata01 <- left_join(F01, cov1, by = "eid")
mydata30 <- left_join(G30, cov1, by = "eid")

# Remove redundant index columns
mydata$X <- NULL
mydata30$X <- NULL

{
  # Load behavioral and psychosocial phenotype dataset
  shuzhen_df_out <- read.csv("./input/202emotion.csv")  
  shuzhen_df_out <- shuzhen_df_out[ ,-1]  
  colnames(shuzhen_df_out)[colnames(shuzhen_df_out) == "f.eid"] <- "eid"  
  
  {
    # Extract continuous behavioral and psychosocial measures
    selected_columns <- c(
      "eid",
      "f.2050.0.0","f.2060.0.0","f.2070.0.0","f.2080.0.0","f.4526.0.0",
      "f.4537.0.0","f.4548.0.0","f.4559.0.0","f.4570.0.0","f.4581.0.0",
      "f.4609.0.0","f.4620.0.0","f.20127.0.0","p28735","p28736","p28737",
      "p28738","p28739","p28740","p28741","p28742","p28743","p28744",
      "p28745","p28746","p28747","p28748","p28749","p28750","p28751",
      "p28752","f.20438.0.0","f.20436.0.0","f.20439.0.0","f.20440.0.0",
      "f.20442.0.0","f.20518.0.0","f.20510.0.0","f.20507.0.0","f.20519.0.0",
      "f.20514.0.0","f.20511.0.0","f.20513.0.0","f.20508.0.0","f.20517.0.0",
      "f.20505.0.0","f.20512.0.0","f.20506.0.0","f.20509.0.0","f.20516.0.0",
      "f.20515.0.0","f.20520.0.0","f.20414.0.0","f.20403.0.0","f.20416.0.0",
      "f.20413.0.0","f.20407.0.0","f.20412.0.0","f.20409.0.0","f.20408.0.0",
      "f.20453.0.0","f.20489.0.0","f.20488.0.0","f.20487.0.0","f.20490.0.0",
      "f.20491.0.0","f.20522.0.0","f.20523.0.0","f.20521.0.0","f.20524.0.0",
      "f.20525.0.0","f.20497.0.0","f.20498.0.0","f.20495.0.0","f.20496.0.0",
      "f.20494.0.0","f.20479.0.0","f.20485.0.0","f.20458.0.0","f.20459.0.0",
      "f.20460.0.0"
    )
    
    available_columns <- intersect(selected_columns, names(shuzhen_df_out))
    shuzhen_continues <- shuzhen_df_out[available_columns]
  }
  
  {
    # Extract categorical behavioral and lifestyle indicators
    selected_columns <- c(
      "eid",
      "f.20499.0.0","f.20500.0.0","f.20446.0.0","f.20441.0.0","f.20447.0.0",
      "f.20532.0.0","f.20435.0.0","f.20449.0.0","f.20450.0.0","f.20448.0.0",
      "f.20534.0.0","f.20437.0.0","f.20533.0.0","f.20535.0.0","f.20536.0.0",
      "f.20502.0.0","f.20501.0.0","f.20421.0.0","f.20425.0.0","f.20401.0.0",
      "f.20411.0.0","f.20405.0.0","f.20468.0.0","f.20474.0.0","f.20463.0.0",
      "f.20471.0.0","f.20531.0.0","f.20529.0.0","f.20526.0.0","f.20530.0.0",
      "f.20528.0.0","f.20527.0.0","f.20480.0.0","f.1920.0.0","f.1930.0.0",
      "f.1940.0.0","f.1950.0.0","f.1960.0.0","f.1970.0.0","f.1980.0.0",
      "f.1990.0.0","f.2000.0.0","f.2010.0.0","f.2020.0.0","f.2030.0.0",
      "f.2040.0.0","f.2090.0.0","f.2100.0.0","f.4598.0.0","f.4631.0.0",
      "f.4642.0.0","f.4653.0.0","f.20126.0.0","new_6145_1","new_6145_2",
      "new_6145_3","new_6145_4","new_6145_5","new_6145_6","p29016",
      "p29019","p29022","p29028","p29035","new_p28753_0","new_p28753_4",
      "new_p28753_6","new_p28753_8","new_p28753_1","new_p28753_7",
      "new_p28753_2","new_p28753_3","new_p28753_5","new_p28753_9",
      "new_p29017_0","new_p29017_1","new_p29017_2","new_p29020_2",
      "new_p29020_0","new_p29020_1","new_p29021_1","new_p29021_3",
      "new_p29021_0","new_p29021_2","new_p29032_0","new_p29032_2",
      "new_p29032_1","new_p29038_2","new_p29038_1","new_p29038_3",
      "new_p29038_0","new_p29039_5","new_p29039_4","new_p29039_3",
      "new_p29039_1","new_p29039_2","new_p29039_7","new_p29039_6",
      "new_p29051_3","new_p29051_6","new_p29051_1","new_p29051_8",
      "new_p29051_0","new_p29051_2","new_p29051_4","new_p29051_5",
      "new_p29051_7","new_p29065_1","new_p29065_4","new_p29065_8",
      "new_p29065_11","new_p29065_13","new_p29065_0","new_p29065_9",
      "new_p29065_10","new_p29065_2","new_p29065_6","new_p29065_3",
      "new_p29065_12","new_p29065_7","new_p29065_5"
    )
    
    available_columns <- intersect(selected_columns, names(shuzhen_df_out))
    shuzhen_category <- shuzhen_df_out[available_columns]
  }
  
  # Harmonize binary categorical indicators
  {
    shuzhen_category$f.20536.0.0[shuzhen_category$f.20536.0.0 %in% c(1,2,3)] <- 1
    shuzhen_category$f.20411.0.0[shuzhen_category$f.20411.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20405.0.0[shuzhen_category$f.20405.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20531.0.0[shuzhen_category$f.20531.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20529.0.0[shuzhen_category$f.20529.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20526.0.0[shuzhen_category$f.20526.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20530.0.0[shuzhen_category$f.20530.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20528.0.0[shuzhen_category$f.20528.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20527.0.0[shuzhen_category$f.20527.0.0 %in% c(1,2)] <- 1
    shuzhen_category$f.20126.0.0[shuzhen_category$f.20126.0.0 %in% c(1,2,3,4,5)] <- 1
  }
  
  shuzhen_category$eid <- as.character(shuzhen_category$eid)
  shuzhen_continues$eid <- as.character(shuzhen_continues$eid)
}

# Construct final analysis datasets for categorical and continuous exposures
full_data_cat   <- left_join(mydata,   shuzhen_category,  by = "eid")
full_data00_cat <- left_join(mydata00, shuzhen_category,  by = "eid")
full_data01_cat <- left_join(mydata01, shuzhen_category,  by = "eid")
full_data30_cat <- left_join(mydata30, shuzhen_category,  by = "eid")

full_data_con   <- left_join(mydata,   shuzhen_continues, by = "eid")
full_data00_con <- left_join(mydata00, shuzhen_continues, by = "eid")
full_data01_con <- left_join(mydata01, shuzhen_continues, by = "eid")
full_data30_con <- left_join(mydata30, shuzhen_continues, by = "eid")

############### Cox proportional hazards regression (loop-based analysis) #########

# Load libraries for survival analysis and data manipulation
library(survival)   
library(tidyverse) 

# Define covariates included in all Cox models
covariates <- c("age", "sex", "race", "bmi", "smoke", "drink", "APOEe4_carrier")

# Initialize container for regression results
all_results <- list()

# Define datasets with categorical exposures
cat_dataframes <- c("full_data_cat", "full_data00_cat", "full_data01_cat", "full_data30_cat")

# Define datasets with continuous exposures
con_dataframes <- c("full_data_con", "full_data00_con", "full_data01_con", 
                    "full_data30_con")

# Check availability of required datasets
missing_cat_dfs <- setdiff(cat_dataframes, ls())
missing_con_dfs <- setdiff(con_dataframes, ls())

# Stop analysis if required datasets are missing
if (length(missing_cat_dfs) > 0) {
  stop(paste("Missing categorical dataframes:", paste(missing_cat_dfs, collapse = ", ")))
}
if (length(missing_con_dfs) > 0) {
  stop(paste("Missing continuous dataframes:", paste(missing_con_dfs, collapse = ", ")))
}

# Loop over datasets containing continuous exposures
for (df_name in con_dataframes) {
  
  df <- get(df_name)
  
  # Specify columns corresponding to exposure variables
  var_cols <- 12:92   
  na_values <- c(NA)  
  
  var_names <- names(df)[var_cols]
  
  # Initialize results table for current dataset
  df_results <- data.frame(
    Dataset = character(),
    Variable = character(),
    N = integer(),
    Coef = numeric(),
    SE = numeric(),
    HR = numeric(),
    CI_lower = numeric(),
    CI_upper = numeric(),
    p_value = numeric(),
    pFDR = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Fit Cox models for each exposure variable
  for (var in var_names) {
    
    # Exclude observations with missing exposure values
    rows_to_keep <- !(df[[var]] %in% na_values)
    df_clean <- df[rows_to_keep, ]
    
    n <- sum(rows_to_keep)
    if (n == 0) next  
    
    # Construct Cox proportional hazards model
    formula <- as.formula(
      paste("Surv(time, status) ~", var, "+", paste(covariates, collapse = "+"))
    )
    
    fit <- tryCatch(
      coxph(formula, data = df_clean),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    summary_fit <- summary(fit)
    
    coef_row <- which(rownames(summary_fit$coefficients) == var)
    confint_row <- which(rownames(summary_fit$conf.int) == var)
    
    # Extract hazard ratios and statistical estimates
    df_results <- rbind(df_results, data.frame(
      Dataset = df_name,
      Variable = var,
      N = n,
      Coef = summary_fit$coefficients[coef_row, "coef"],
      SE = summary_fit$coefficients[coef_row, "se(coef)"],
      HR = summary_fit$conf.int[confint_row, "exp(coef)"],
      CI_lower = summary_fit$conf.int[confint_row, "lower .95"],
      CI_upper = summary_fit$conf.int[confint_row, "upper .95"],
      p_value = summary_fit$coefficients[coef_row, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    ))
  }
  
  # Adjust p values for multiple testing using FDR
  if (nrow(df_results) > 0) {
    df_results$pFDR <- p.adjust(df_results$p_value, method = "fdr")
  }
  
  all_results[[df_name]] <- df_results
}

# Combine results across all datasets
final_results <- bind_rows(all_results)

# Load packages for survival analysis and data handling
library(survival)
library(broom)
library(tidyverse)

# Define covariates included in all Cox models
covariates <- c("age", "sex", "race", "bmi", "smoke", "drink", "APOEe4_carrier")

# Specify datasets containing categorical exposures
df_names <- c("full_data_cat", "full_data00_cat", "full_data01_cat", "full_data30_cat")

# Initialize list to store regression results
all_results1 <- list()

# Perform Cox proportional hazards regression for categorical exposures
for (df_name in df_names) {
  df <- get(df_name)
  
  current_results <- data.frame()
  
  # Columns corresponding to categorical exposure variables
  var_cols <- 12:132
  var_names <- colnames(df)[var_cols]
  
  # Iterate over each categorical exposure
  for (var in var_names) {
    
    # Exclude missing or invalid category codes
    temp_df <- df %>%
      filter(!(!!sym(var) %in% c(999, 9999)), 
             !is.na(!!sym(var)))
    
    # Skip variables without sufficient variability
    if (length(unique(temp_df[[var]])) < 2) {
      next
    }
    
    # Construct Cox proportional hazards model
    formula <- as.formula(
      paste("Surv(time, status) ~", var, "+", 
            paste(covariates, collapse = " + "))
    )
    
    cox_model <- tryCatch(
      coxph(formula, data = temp_df),
      error = function(e) NULL
    )
    
    # Extract model estimates if fitting succeeded
    if (!is.null(cox_model)) {
      n <- nrow(temp_df)
      model_summary <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
      var_result <- model_summary[1, ]
      
      # Store hazard ratio estimates and confidence intervals
      result_row <- data.frame(
        dataset = df_name,
        variable = var,
        n = n,
        coef = var_result$estimate,
        se = var_result$std.error,
        HR = var_result$estimate,
        HR_95CI_lower = var_result$conf.low,
        HR_95CI_upper = var_result$conf.high,
        p_value = var_result$p.value,
        stringsAsFactors = FALSE
      )
      
      current_results <- bind_rows(current_results, result_row)
    }
  }
  
  # Adjust p values for multiple testing using FDR
  if (nrow(current_results) > 0) {
    current_results$fdr_p_value <- p.adjust(current_results$p_value, method = "fdr")
  }
  
  all_results1[[df_name]] <- current_results
}

# Combine categorical Cox regression results across datasets
final_results1 <- bind_rows(all_results1, .id = "dataset")

################# Integration with continuous-exposure results ###################

# Harmonize column names and merge categorical and continuous analyses
colnames(final_results1) <- colnames(final_results)
final_combined <- rbind(final_results, final_results1)

# Annotate exposures using external mapping file
matchdata <- read.csv("./input/202category.csv")

final_combined <- left_join(
  final_combined, 
  matchdata, 
  by = c("Variable" = "outcome")
)

######################################################
# Assemble dementia outcome dataset with covariates
mydata <- left_join(allcause_dementia, cov1, by = "eid")

{
  # Load behavioral and actigraphy-derived sleep variables
  file_cox_all <- read.csv("./input/file_cox_all.csv")
  names(file_cox_all)[names(file_cox_all) == "f.eid"] <- "eid"
  file_cox_all$eid <- as.character(file_cox_all$eid)
  
  file_cox_all <- file_cox_all %>%
    select(
      Sleep_duration_touchscreen,
      Getting_up_in_morning_touchscreen,
      chronotype_touchscreen,
      Nap_during_day_touchscreen,
      Sleeplessness_insomnia_touchscreen,
      Snoring_touchscreen,
      Daytime_dozing_sleeping_touchscreen,
      morning_person_derived,
      AD_L5hr_ENMO_mg_0.24hr,
      AD_L5_ENMO_mg_0.24hr,
      AD_M5hr_ENMO_mg_0.24hr,
      AD_M5_ENMO_mg_0.24hr,
      IS_interdailystability,
      IV_intradailyvariability,
      MVPA_min_morning,
      MVPA_min_noon,
      MVPA_min_evening,
      ENMO_morning,
      ENMO_noon,
      ENMO_evening,
      sleeponset_AD_T5A5_mn,
      sleeponset_AD_T5A5_sd,
      SptDuration_AD_T5A5_mn,
      SptDuration_AD_T5A5_sd,
      SleepDurationInSpt_AD_T5A5_mn,
      SleepDurationInSpt_AD_T5A5_sd,
      WASO_AD_T5A5_mn,
      WASO_AD_T5A5_sd,
      SleepRegularityIndex_AD_T5A5_mn,
      SleepRegularityIndex_AD_T5A5_sd,
      number_sib_sleepperiod_AD_T5A5_mn,
      number_sib_sleepperiod_AD_T5A5_sd,
      sleep_efficiency_mean,
      Nbouts_day_IN_bts_30_mean,
      Nbouts_day_IN_bts_20_30_mean,
      Nbouts_day_IN_bts_10_20_mean,
      mid_L5,
      mid_M5,
      mid_SPT,
      delta_wd_we_sptduration,
      delta_wd_we_sleepINspt,
      Job_involves_shift_work,
      Job_involves_night_shift_work,
      sleep_too_much_SptDuration,
      sleep_too_low_SptDuration,
      sleep_too_much_SleepDurationInSpt,
      sleep_too_low_SleepDurationInSpt,
      Nbouts_day_IN_bts_mean_all,
      activity_pattern,
      eid
    )
  
  # Load biochemical markers
  bio <- read.csv("./input/bio.csv") |> mutate(eid = as.character(eid))
  
  # Load lifestyle and environmental variables
  art <- read.csv("./input/article_data.csv", header = TRUE) |> 
    mutate(eid = as.character(eid))
  
  # Rename selected variables to descriptive labels
  name_map <- c(
    "p1160" = "Sleep_duration",
    "p1239" = "Current_tobacco_smoking",
    "p1249" = "Past_tobacco_smoking",
    "p1488" = "Tea_intake",
    "p1498" = "Coffee_intake",
    "p1528" = "Water_intake",
    "p24503" = "Greenspace_percentage_buffer_300m",
    "p24504" = "Domestic_garden_percentage_buffer_300m",
    "p24505" = "Water_percentage_buffer_300m",
    "p24506" = "Natural_environment_percentage_buffer_1000m",
    "p24507" = "Natural_environment_percentage_buffer_300m"
  )
  
  matched_names <- names(name_map)[names(name_map) %in% colnames(art)]
  colnames(art)[match(matched_names, colnames(art))] <- name_map[matched_names]

  ## --- Lifestyle factors ---
  # Extract lifestyle-related variables from baseline questionnaire data
  cols <- c("eid","Average_total_household_income_before_tax",
            "Current_employment_status",
            "Qualifications",
            "Breastfed_as_a_baby",
            "Comparative_body_size_at_age_10",
            "Comparative_height_size_at_age_10",
            "Maternal_smoking_around_birth",
            "Above_moderate_vigorous_walking_recommendation",
            "Age_first_had_sexual_intercourse",
            "Coffee_intake",
            "Current_tobacco_smoking", 
            "Daytime_dozing_sleeping_narcolepsy",
            "Ever_had_same_sex_intercourse",
            "Exposure_to_tobacco_smoke_at_home",
            "Exposure_to_tobacco_smoke_outside_home",
            "Hands_free_device_speakerphone_use_with_mobile_phone_in_last_3_month",
            "Leisure_social_activities",
            "Lifetime_number_of_sexual_partners",                     
            "MET_minutes_per_week_for_moderate_activity",
            "MET_minutes_per_week_for_vigorous_activity",
            "MET_minutes_per_week_for_walking",
            "Mineral_and_other_dietary_supplements",
            "Nap_during_day",
            "Past_tobacco_smoking",
            "Sleep_duration",
            "Sleeplessness_insomnia",
            "Snoring",
            "Summed_MET_minutes_per_week_for_all_activity",             
            "Tea_intake",                                                
            "Time_spend_outdoors_in_summer",                            
            "Time_spent_outdoors_in_winter",                            
            "Time_spent_using_computer",                                 
            "Time_spent_watching_television_TV",                         
            "Vitamin_and_mineral_supplements",
            "Water_intake",
            "Weekly_usage_of_mobile_phone_in_last_3_months")
  Lifestyles <- art[, intersect(cols, colnames(art))]
  colnames(Lifestyles)
  
  ## --- Psychosocial factors ---
  # Extract psychosocial and mental well-being variables
  cols <- c("eid",
            "Able_to_confide",
            "Fed_up_feelings",
            "Guilty_feelings",
            "Irritability",
            "Loneliness_isolation",
            "Miserableness",
            "Mood_swings",
            "Nervous_feelings",
            "Sensitivity_hurt_feelings",
            "Suffer_from_nerves",
            "Tense_highly_strung",
            "Worrier_anxious_feelings",
            "Worry_too_long_after_embarrassment")
  Psychosocial_factors <- art[, intersect(cols, colnames(art))]
  colnames(Psychosocial_factors)
  
  ## --- Local environmental exposures ---
  # Extract residential environmental and pollution-related variables
  cols <- c("eid",
            "Water_percentage_buffer_1000m",
            "Water_percentage_buffer_300m",
            "Nitrogen_dioxide_air_pollution_a", 
            "Nitrogen_dioxide_air_pollution_b",                         
            "Nitrogen_dioxide_air_pollution_c", 
            "Nitrogen_dioxide_air_pollution_d",
            "Nitrogen_oxides_air_pollution",                            
            "Particulate_matter_air_pollution_2_5_10um", 
            "Particulate_matter_air_pollution_pm10_a",                
            "Particulate_matter_air_pollution_pm10_b",                  
            "Particulate_matter_air_pollution_pm2_5",
            "Natural_environment_percentage_buffer_1000m",                 
            "Natural_environment_percentage_buffer_300m",
            "Greenspace_percentage_buffer_1000m",                       
            "Greenspace_percentage_buffer_300m",
            "Domestic_garden_percentage_buffer_1000m",             
            "Domestic_garden_percentage_buffer_300m",
            "Average_16_hour_sound_level_of_noise_pollution",        
            "Average_24_hour_sound_level_of_noise_pollution",           
            "Average_daytime_sound_level_of_noise_pollution",           
            "Average_evening_sound_level_of_noise_pollution",            
            "Average_night_time_sound_level_of_noise_pollution")
  Local_environment <- art[, intersect(cols, colnames(art))]
  colnames(Local_environment)
  
  # Derive composite air pollution metrics
  Local_environment <- Local_environment %>%
    mutate(Nitrogen_air_pollution = rowMeans(select(., Nitrogen_dioxide_air_pollution_a,
                                                    Nitrogen_dioxide_air_pollution_b,
                                                    Nitrogen_dioxide_air_pollution_c,
                                                    Nitrogen_dioxide_air_pollution_d,
                                                    Nitrogen_oxides_air_pollution),
                                             na.rm = TRUE)) %>%
    rename(Nitrogen_oxides = Nitrogen_air_pollution) %>%
    mutate(Particulate_matter_air_pollution_pm10 = rowMeans(select(.,
                                                                   Particulate_matter_air_pollution_pm10_a,
                                                                   Particulate_matter_air_pollution_pm10_b),
                                                            na.rm = TRUE)) %>%
    select(-all_of(c("Particulate_matter_air_pollution_pm10_a",
                     "Particulate_matter_air_pollution_pm10_b",
                     "Nitrogen_dioxide_air_pollution_a",
                     "Nitrogen_dioxide_air_pollution_b",
                     "Nitrogen_dioxide_air_pollution_c",
                     "Nitrogen_dioxide_air_pollution_d",
                     "Nitrogen_oxides_air_pollution")))
  
  # Load dietary intake variables
  Food_Intake <- read.csv("./input/food_take_1229.csv") |>
    mutate(eid = as.character(eid))
  sort(colnames(Food_Intake))
  
  Food_Intake <- Food_Intake %>%
    select(-c("p26000_i0", "p26001_i0", "p26004_i0", 
              "p26042_i4", "p26112_i4"))
  
  # Load anthropometric measurements
  Body_measure <- read.csv("./input/39_body_measure.csv", header = TRUE) |>
    mutate(eid = as.character(eid))
  
  # Harmonize variable names
  new_names <- colnames(Body_measure)
  new_names <- gsub("[^A-Za-z0-9]", "_", new_names)
  new_names <- gsub("_+", "_", new_names)
  new_names <- gsub("^_+|_+$", "", new_names)
  colnames(Body_measure) <- new_names
  colnames(Body_measure)
  
  # Load food preference data
  liking <- read.csv("./input/food_liking.tsv", sep = "\t") |>
    mutate(eid = as.character(eid))
  
  liking <- liking %>%
    select(-c("f.20614.0.0", "f.20656.0.0", "f.20657.0.0", 
              "f.20668.0.0", "f.20669.0.0", "f.20670.0.0", 
              "f.20741.0.0", "f.20749.0.0"))
  
  # Merge multi-domain exposure data by participant identifier
  data_list <- list(
    bio, Lifestyles, Local_environment, liking,
    Psychosocial_factors, Food_Intake, Body_measure, file_cox_all
  ) |>
    map(~ mutate(., eid = as.character(eid)))
  
  lapply(data_list, function(df) class(df$eid))
  outcome_data <- reduce(data_list, full_join, by = "eid")
}

outcome_data <- outcome_data %>% 
  mutate(eid = as.character(eid))
full_data <- left_join(mydata, outcome_data, by = "eid")

# Load libraries for survival analysis, data manipulation, and progress monitoring
library(progress)
library(survival)
library(dplyr)

# Define survival outcome variables and candidate exposure set
time_var <- "time"     
status_var <- "status"  
exposures <- setdiff(names(outcome_data), c("eid", "X"))  

# Initialize progress bar for high-throughput Cox regression
pb <- progress_bar$new(
  total = length(exposures),  
  format = "[:bar] :percent :elapsed",  
  clear = FALSE  
)

results <- list()

# Loop over each exposure variable
for (expo in exposures) {
  
  # Specify exposure, outcome, and adjustment covariates
  required_vars <- c(expo, time_var, status_var,
                     "age", "sex", "race", "bmi", "smoke", "drink", "APOEe4_carrier")
  
  # Skip exposure if required variables are unavailable
  if (!all(required_vars %in% names(full_data))) next
  
  # Subset data and exclude missing observations
  df <- full_data |> 
    select(all_of(required_vars)) |>  
    filter(
      !is.na(!!sym(expo)),          
      !is.na(!!sym(time_var)),      
      !is.na(!!sym(status_var))     
    )
  
  # Skip exposure if no valid samples remain
  if (nrow(df) == 0) next
  
  # Construct Cox proportional hazards model with covariate adjustment
  formula_str <- paste(
    "Surv(", time_var, ", ", status_var, ") ~",
    expo, "+ age + sex + race + bmi + smoke + drink + APOEe4_carrier"
  )
  
  # Fit Cox model with error handling
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = df),
    error = function(e) NULL
  )
  
  # Extract exposure-specific hazard ratio estimates
  if (!is.null(fit)) {
    tidy_fit <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == expo) %>%  
      mutate(
        exposure = expo,        
        n_samples = nrow(df),   
        hazard_ratio = exp(estimate),  
        conf.low = exp(conf.low),      
        conf.high = exp(conf.high)     
      )
    
    results[[expo]] <- tidy_fit
  }
  
  # Update progress indicator
  pb$tick()  
}

# Combine results across all exposures
final_result <- bind_rows(results)

# Adjust P values for multiple testing using FDR
final_result$p_fdr <- p.adjust(final_result$p.value, method = "fdr")

# Report any warnings generated during model fitting
warnings()
warnings()

# ========================== Result Processing & Annotation ==========================

# Define variable category lists (extract exposure variables by dataset, exclude ID)
fi_vars <- setdiff(names(Food_Intake), "eid")          
bm_vars <- setdiff(names(Body_measure), "eid")        
li_vars <- setdiff(names(liking), "eid")              
pf_vars <- setdiff(names(Psychosocial_factors), "eid")
ls_vars <- setdiff(names(Lifestyles), "eid")          
le_vars <- setdiff(names(Local_environment), "eid")   
bio_vars <- setdiff(names(bio), "eid")                
fix_vars <- setdiff(names(file_cox_all), "eid")       

final_result <- final_result %>%
  mutate(Category = case_when(
    exposure %in% fi_vars ~ "Food Intake",
    exposure %in% bm_vars ~ "Body measure",
    exposure %in% li_vars ~ "Food liking",
    exposure %in% pf_vars ~ "Psychosocial factors",
    exposure %in% ls_vars ~ "Lifestyles",
    exposure %in% le_vars ~ "Local environment",
    exposure %in% bio_vars ~ "Biochemical markers",
    exposure %in% fix_vars ~ "file_cox_all",
    TRUE ~ "Unknown"  # Label unclassified variables as "Unknown"
  ))

library(readxl)

data1 <- read_excel("./input/match_data.xlsx", sheet = "Sheet1")

names(data1)[names(data1) == "exprosure"] <- "exposure"

phewas_data <- left_join(final_result, data1, by = "exposure")

phewas_data <- phewas_data %>%
  mutate(exprosure1 = if_else(is.na(exprosure1), exposure, exprosure1))

# Check number of missing values in the exposure name column (quality control)
sum(is.na(phewas_data$exprosure1))
