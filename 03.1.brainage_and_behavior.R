##############################Behavioral group1##################################
# Load required libraries
library(dplyr)
library(purrr)
library(progress)
library(broom)
library(dplyr)

file_cox_all<-read.csv("./input/file_cox_all.csv")
names(file_cox_all)[names(file_cox_all) == "f.eid"] <- "eid"
file_cox_all$eid<-as.character(file_cox_all$eid)

file_cox_all <- file_cox_all %>% select("Sleep_duration_touchscreen" ,        
                                        "Getting_up_in_morning_touchscreen",  
                                        "chronotype_touchscreen"            , 
                                        "Nap_during_day_touchscreen"         ,
                                        "Sleeplessness_insomnia_touchscreen" ,
                                        "Snoring_touchscreen"                ,
                                        "Daytime_dozing_sleeping_touchscreen",
                                        "morning_person_derived"   ,
                                        "AD_L5hr_ENMO_mg_0.24hr"   ,          
                                        "AD_L5_ENMO_mg_0.24hr"      ,         
                                        "AD_M5hr_ENMO_mg_0.24hr"     ,        
                                        "AD_M5_ENMO_mg_0.24hr"        ,       
                                        "IS_interdailystability"      ,       
                                        "IV_intradailyvariability"     ,      
                                        "MVPA_min_morning"               ,    
                                        "MVPA_min_noon"                   ,   
                                        "MVPA_min_evening"                 ,  
                                        "ENMO_morning"                      , 
                                        "ENMO_noon"                          ,
                                        "ENMO_evening"                       ,
                                        "sleeponset_AD_T5A5_mn"              ,
                                        "sleeponset_AD_T5A5_sd"              ,
                                        "SptDuration_AD_T5A5_mn"             ,
                                        "SptDuration_AD_T5A5_sd"             ,
                                        "SleepDurationInSpt_AD_T5A5_mn"      ,
                                        "SleepDurationInSpt_AD_T5A5_sd"      ,
                                        "WASO_AD_T5A5_mn"                    ,
                                        "WASO_AD_T5A5_sd"                    ,
                                        "SleepRegularityIndex_AD_T5A5_mn"    ,
                                        "SleepRegularityIndex_AD_T5A5_sd"    ,
                                        "number_sib_sleepperiod_AD_T5A5_mn"  ,
                                        "number_sib_sleepperiod_AD_T5A5_sd"  ,
                                        "sleep_efficiency_mean"              ,
                                        "Nbouts_day_IN_bts_30_mean"          ,
                                        "Nbouts_day_IN_bts_20_30_mean"       ,
                                        "Nbouts_day_IN_bts_10_20_mean"       ,
                                        "mid_L5"                             ,
                                        "mid_M5"                             ,
                                        "mid_SPT"                            ,
                                        "delta_wd_we_sptduration"            ,
                                        "delta_wd_we_sleepINspt"             ,
                                        "Job_involves_shift_work"            ,
                                        "Job_involves_night_shift_work"      ,
                                        "sleep_too_much_SptDuration"         ,
                                        "sleep_too_low_SptDuration"          ,
                                        "sleep_too_much_SleepDurationInSpt"  ,
                                        "sleep_too_low_SleepDurationInSpt"   ,
                                        "Nbouts_day_IN_bts_mean_all"  ,
                                        "activity_pattern",
                                        "eid")

# Load biochemical marker data
bio <- read.csv("./input/bio.csv")

# Convert eid to character type
bio$eid<-as.character(bio$eid)

# Load base article data and ensure eid is character type
art <- read.csv("./input/article_data.csv", header = TRUE) |> mutate(eid = as.character(eid)) 

# Define mapping from field IDs to descriptive names
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

# Identify existing columns to rename
matched_names <- names(name_map)[names(name_map) %in% colnames(art)]

# Replace column names with descriptive names
colnames(art)[match(matched_names, colnames(art))] <- name_map[matched_names]

## --- Lifestyle Factors ---
# Select lifestyle-related columns from article data
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

## --- Psychosocial Factors ---
# Select psychosocial-related columns from article data
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

## --- Local Environment Variables ---
# Select environment-related columns from article data
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

# Combine pollutants into average scores
Local_environment <- Local_environment %>%
  mutate(Nitrogen_air_pollution = rowMeans(select(., Nitrogen_dioxide_air_pollution_a,
                                                  Nitrogen_dioxide_air_pollution_b,
                                                  Nitrogen_dioxide_air_pollution_c,
                                                  Nitrogen_dioxide_air_pollution_d,
                                                  Nitrogen_oxides_air_pollution), na.rm = TRUE)) %>%
  rename(Nitrogen_oxides = Nitrogen_air_pollution) %>%
  mutate(Particulate_matter_air_pollution_pm10 = rowMeans(select(.,
                                                                 Particulate_matter_air_pollution_pm10_a,
                                                                 Particulate_matter_air_pollution_pm10_b), na.rm = TRUE)) %>%
  select(-all_of(c("Particulate_matter_air_pollution_pm10_a",
                   "Particulate_matter_air_pollution_pm10_b",
                   "Nitrogen_dioxide_air_pollution_a",
                   "Nitrogen_dioxide_air_pollution_b",
                   "Nitrogen_dioxide_air_pollution_c",
                   "Nitrogen_dioxide_air_pollution_d",
                   "Nitrogen_oxides_air_pollution")))

# Load food intake data
Food_Intake <- read.csv("./input/food_take_1229.csv") |> mutate(eid = as.character(eid))
sort(colnames(Food_Intake))

# Remove specific columns from food intake data
Food_Intake <- Food_Intake %>%
  select(-c("p26000_i0", "p26001_i0", "p26004_i0", 
            "p26042_i4", "p26112_i4"))

# Load body measurement data
Body_measure <- read.csv("./input/39_body_measure.csv", header = TRUE) |> mutate(eid = as.character(eid)) 
new_names <- colnames(Body_measure)

# Step 1: Replace all non-alphanumeric characters with "_"
new_names <- gsub("[^A-Za-z0-9]", "_", new_names)

# Step 2: Replace consecutive "_" with single "_"
new_names <- gsub("_+", "_", new_names)

# Step 3: Remove leading/trailing "_"
new_names <- gsub("^_+|_+$", "", new_names)

# Apply new column names
colnames(Body_measure) <- new_names
colnames(Body_measure)

# Load food liking data
liking <- read.csv("./input/food_liking.tsv", sep="\t") |> 
  mutate(eid = as.character(eid))

# Remove non-food-liking columns from food liking data
liking <- liking %>%
  select(-c("f.20614.0.0", "f.20656.0.0", "f.20657.0.0", 
            "f.20668.0.0", "f.20669.0.0", "f.20670.0.0", 
            "f.20741.0.0", "f.20749.0.0"))

data_list <- list(
  bio, Lifestyles, Local_environment, liking,
  Psychosocial_factors, Food_Intake, Body_measure, file_cox_all
) |>
  map(~ mutate(., eid = as.character(eid))) 

lapply(data_list, function(df) class(df$eid))
outcome_data <- reduce(data_list, full_join, by = "eid")

# Load exposure protein data
mydata <- read.csv("./input/brain_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         brain_difference = Bias_Corrected_Age - age)
mydata1 <- read.csv("./input/baizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         baizhi_difference = Bias_Corrected_Age - age)
mydata2 <- read.csv("./input/huizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         huizhi_difference = Bias_Corrected_Age - age)
head(mydata1)
head(mydata2)
# Extract required columns
baizhi <- mydata1[, c("eid", "baizhi_difference")]
huizhi <- mydata2[, c("eid", "huizhi_difference")]
# Merge by eid
mydata_merged <- merge(baizhi, huizhi, by = "eid", all = TRUE)
mydata <- mydata %>%
  left_join(mydata_merged, by = "eid")
mydata$age<-NULL
mydata$Bias_Corrected_Age<-NULL

# Load covariates data
cov1 <- read.csv("./input/M45_M51+cov.csv", header = T) |> 
  mutate(eid = as.character(eid)) |> 
  select(eid, age, sex, race, bmi, smoke, drink, APOEe4_carrier) 
# Process education variable
# Process race variable
{
  cov1 <- cov1 %>%
    mutate(race = case_when(
      race == "-3" ~ "9",
      race == "-1" ~ "9",
      race == "1" ~ "2",
      race == "2" ~ "2",
      race == "3" ~ "2",
      race == "4" ~ "2",
      race == "5" ~ "2",
      race == "6" ~ "2",
      race == "1001" ~ "1",
      race == "1002" ~ "2",
      race == "1003" ~ "2",
      race == "2001" ~ "2",
      race == "2002" ~ "2",
      race == "2003" ~ "2",
      race == "2004" ~ "2",
      race == "3001" ~ "2",
      race == "3002" ~ "2",
      race == "3003" ~ "2",
      race == "3004" ~ "2",
      race == "4001" ~ "2",
      race == "4002" ~ "2",
      race == "4003" ~ "2",
      is.na(race) ~ "9",
      TRUE ~ as.character(race)
    ))
  
  cov1 <- cov1 %>%
    mutate(smoke = case_when(
      smoke == "-3" ~ "9",
      smoke == "0" ~ "0",
      smoke == "1" ~ "1",
      smoke == "2" ~ "2",
      is.na(smoke) ~ "9",
      TRUE ~ as.character(smoke)
    ))
  
  cov1 <- cov1 %>%
    mutate(drink = case_when(
      drink == "-3" ~ "9",
      drink == "0" ~ "0",
      drink == "1" ~ "1",
      drink == "2" ~ "2",
      is.na(drink) ~ "9",
      TRUE ~ as.character(drink)
    ))
  
  cov1 <- cov1 %>%
    mutate(APOEe4_carrier = case_when(
      APOEe4_carrier == "9" ~ "9",
      APOEe4_carrier == "0" ~ "0",
      APOEe4_carrier == "1" ~ "1",
      is.na(APOEe4_carrier) ~ "9",
      TRUE ~ as.character(APOEe4_carrier)
    ))
  }

# Convert variables to factors
cov1$race<-as.factor(cov1$race)
cov1$sex<-as.factor(cov1$sex)
cov1$smoke<-as.factor(cov1$smoke)
cov1$drink<-as.factor(cov1$drink)
cov1$APOEe4_carrier<-as.factor(cov1$APOEe4_carrier)

# Merge all into final dataset
mydata <- left_join(mydata, cov1, by = "eid")
full_data <- left_join(mydata, outcome_data, by = "eid")
head(full_data)


# Define outcome and exposure variables
outcomes <- c("brain_difference", "baizhi_difference","huizhi_difference")
exposures <- setdiff(names(outcome_data), "eid")

library(progress)
library(broom)
# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      df <- df |> filter(!is.na(.data[[expo]]) & !is.na(.data[[outcome]]) & !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    
    n <- nrow(df) 
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ ", expo, " + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      # Calculate standard deviations for standardization
      sd_expo <- sd(df[[expo]], na.rm = TRUE)
      sd_outcome <- sd(df[[outcome]], na.rm = TRUE)
      this_result <- tidy_fit |> 
        filter(term == expo) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          sd_exposure = sd_expo,
          sd_outcome = sd_outcome,
          n = n, 
          beta_std = estimate * (sd_expo / sd_outcome)
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value,
               sd_exposure, sd_outcome, n,beta_std)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}

# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

# Define variable categories for each dataset
fi_vars <- setdiff(names(Food_Intake), "eid")
bm_vars <- setdiff(names(Body_measure), "eid")
li_vars <- setdiff(names(liking), "eid")
pf_vars <- setdiff(names(Psychosocial_factors), "eid")
ls_vars <- setdiff(names(Lifestyles), "eid")
le_vars <- setdiff(names(Local_environment), "eid")
bio_vars <- setdiff(names(bio), "eid")
fix_vars <- setdiff(names(file_cox_all), "eid")
# Assign categories to each variable
library(dplyr)
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
    TRUE ~ "Unknown"
  ))

#exposure %in% met_vars ~ "Metabolite",

max(final_result$n, na.rm = TRUE)
min(final_result$n, na.rm = TRUE)

final_result %>%
  group_by(outcome) %>%
  summarise(
    max_n = max(n, na.rm = TRUE),
    min_n = min(n, na.rm = TRUE)
  )

library(readxl)
# Load additional mapping data
data1 <- read_excel("./input/match_data.xlsx", sheet = "Sheet1")

names(data1)[names(data1) == "exprosure"] <- "exposure"
# Join with mapping data
phewas_data<-left_join(final_result,data1, by="exposure")

# Fill in missing exposure names
phewas_data <- phewas_data %>%
  mutate(exprosure1 = if_else(is.na(exprosure1), exposure, exprosure1))

# Check for missing values
sum(is.na(phewas_data$exprosure1))

# Clean up food intake variable names
phewas_data$outcome[phewas_data$Category == "Food Intake"] <- gsub("_0$", "", phewas_data$outcome[phewas_data$Category == "Food Intake"])

# Save final results
#write.csv(phewas_data,"./output/three_imaging_brain_age.csv",row.names = FALSE)


mydata <- read.csv("./input/Predicted_Age.csv")
mydata <- mydata |> 
  rename(eid = SampleID) |>  
  mutate(eid = as.character(eid)) 
# Merge all into final dataset
mydata <- left_join(mydata, cov1, by = "eid")

mydata <- mydata %>%
  mutate(protein_difference = PredictedAge - age, .after = 2) 
full_data <- left_join(mydata, outcome_data, by = "eid")

# Define outcome and exposure variables
outcomes <- c("protein_difference") 
exposures <- setdiff(names(outcome_data), "eid")

library(progress)
library(broom)
# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      df <- df |> filter(!is.na(.data[[expo]]) & !is.na(.data[[outcome]]) & !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    n <- nrow(df)  
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ ", expo, " + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      # Calculate standard deviations for standardization
      sd_expo <- sd(df[[expo]], na.rm = TRUE)
      sd_outcome <- sd(df[[outcome]], na.rm = TRUE)
      this_result <- tidy_fit |> 
        filter(term == expo) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          sd_exposure = sd_expo,
          sd_outcome = sd_outcome,
          n = n, 
          beta_std = estimate * (sd_expo / sd_outcome)
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value,
               sd_exposure, sd_outcome, n, beta_std)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}

# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

# Define variable categories for each dataset
fi_vars <- setdiff(names(Food_Intake), "eid")
bm_vars <- setdiff(names(Body_measure), "eid")
li_vars <- setdiff(names(liking), "eid")
pf_vars <- setdiff(names(Psychosocial_factors), "eid")
ls_vars <- setdiff(names(Lifestyles), "eid")
le_vars <- setdiff(names(Local_environment), "eid")
bio_vars <- setdiff(names(bio), "eid")
fix_vars <- setdiff(names(file_cox_all), "eid")

# Assign categories to each variable
library(dplyr)
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
    TRUE ~ "Unknown"
  ))


# Display max and min sample sizes
max(final_result$n, na.rm = TRUE)
min(final_result$n, na.rm = TRUE)

# Summarize sample sizes per outcome
final_result %>%
  group_by(outcome) %>%
  summarise(
    max_n = max(n, na.rm = TRUE),
    min_n = min(n, na.rm = TRUE)
  )

# Load additional mapping data
data1 <- read_excel("./input/match_data.xlsx", sheet = "Sheet1")

names(data1)[names(data1) == "exprosure"] <- "exposure"
# Join with mapping data
phewas_data<-left_join(final_result,data1,by="exposure")

# Fill in missing exposure names
phewas_data <- phewas_data %>%
  mutate(exprosure1 = if_else(is.na(exprosure1), exposure, exprosure1))

# Check for missing values
sum(is.na(phewas_data$exprosure1))

# Clean up food intake variable names
phewas_data$outcome[phewas_data$Category == "Food Intake"] <- gsub("_0$", "", phewas_data$outcome[phewas_data$Category == "Food Intake"])

# Save final results
#write.csv(phewas_data,"蛋白脑龄.csv",row.names = FALSE)

##############################Behavioral group2##################################
# Read the phenotype dataset 
shuzhen_df_out <- read.csv("./input/202emotion.csv")  
# Remove the first column 
shuzhen_df_out <- shuzhen_df_out[ ,-1]  
# Rename the column "f.eid" to standard "eid" (unique identifier for participants)
colnames(shuzhen_df_out)[colnames(shuzhen_df_out) == "f.eid"] <- "eid"  

{
  # Define column names to extract
  selected_columns <- c(
    "eid",  
    "f.2050.0.0", "f.2060.0.0", "f.2070.0.0", "f.2080.0.0", "f.4526.0.0", 
    "f.4537.0.0", "f.4548.0.0", "f.4559.0.0", "f.4570.0.0", "f.4581.0.0", 
    "f.4609.0.0", "f.4620.0.0", "f.20127.0.0", "p28735", "p28736", 
    "p28737", "p28738", "p28739", "p28740", "p28741", 
    "p28742", "p28743", "p28744", "p28745", "p28746", 
    "p28747", "p28748", "p28749", "p28750", "p28751", 
    "p28752", "f.20438.0.0", "f.20436.0.0", "f.20439.0.0", "f.20440.0.0", 
    "f.20442.0.0", "f.20518.0.0", "f.20510.0.0", "f.20507.0.0", "f.20519.0.0", 
    "f.20514.0.0", "f.20511.0.0", "f.20513.0.0", "f.20508.0.0", "f.20517.0.0", 
    "f.20505.0.0", "f.20512.0.0", "f.20506.0.0", "f.20509.0.0", "f.20516.0.0", 
    "f.20515.0.0", "f.20520.0.0", "f.20414.0.0", "f.20403.0.0", "f.20416.0.0", 
    "f.20413.0.0", "f.20407.0.0", "f.20412.0.0", "f.20409.0.0", "f.20408.0.0", 
    "f.20453.0.0", "f.20489.0.0", "f.20488.0.0", "f.20487.0.0", "f.20490.0.0", 
    "f.20491.0.0", "f.20522.0.0", "f.20523.0.0", "f.20521.0.0", "f.20524.0.0", 
    "f.20525.0.0", "f.20497.0.0", "f.20498.0.0", "f.20495.0.0", "f.20496.0.0", 
    "f.20494.0.0", "f.20479.0.0", "f.20485.0.0", "f.20458.0.0", "f.20459.0.0",  
    "f.20460.0.0")
  
  # Check which of the defined columns actually exist in the dataframe
  available_columns <- intersect(selected_columns, names(shuzhen_df_out))
  
  shuzhen_continues <- shuzhen_df_out[available_columns]
  
}

{
  # Define column names to extract
  selected_columns <- c(
    "eid", 
    "f.20499.0.0", "f.20500.0.0", "f.20446.0.0", "f.20441.0.0", "f.20447.0.0",
    "f.20532.0.0", "f.20435.0.0", "f.20449.0.0", "f.20450.0.0", "f.20448.0.0", 
    "f.20534.0.0", "f.20437.0.0", "f.20533.0.0", "f.20535.0.0", "f.20536.0.0",
    "f.20502.0.0", "f.20501.0.0", "f.20421.0.0", "f.20425.0.0", "f.20401.0.0", 
    "f.20411.0.0", "f.20405.0.0", "f.20468.0.0", "f.20474.0.0", "f.20463.0.0", 
    "f.20471.0.0", "f.20531.0.0", "f.20529.0.0", "f.20526.0.0", "f.20530.0.0",
    "f.20528.0.0", "f.20527.0.0", "f.20480.0.0", "f.1920.0.0", "f.1930.0.0", 
    "f.1940.0.0", "f.1950.0.0", "f.1960.0.0", "f.1970.0.0", "f.1980.0.0",
    "f.1990.0.0", "f.2000.0.0", "f.2010.0.0", "f.2020.0.0", "f.2030.0.0",
    "f.2040.0.0", "f.2090.0.0", "f.2100.0.0", "f.4598.0.0", "f.4631.0.0",
    "f.4642.0.0", "f.4653.0.0", "f.20126.0.0", "new_6145_1", "new_6145_2", 
    "new_6145_3", "new_6145_4", "new_6145_5", "new_6145_6", "p29016", 
    "p29019", "p29022", "p29028", "p29035", "new_p28753_0", 
    "new_p28753_4", "new_p28753_6", "new_p28753_8", "new_p28753_1", "new_p28753_7", 
    "new_p28753_2", "new_p28753_3", "new_p28753_5", "new_p28753_9", "new_p29017_0",
    "new_p29017_1", "new_p29017_2", "new_p29020_2", "new_p29020_0", "new_p29020_1", 
    "new_p29021_1", "new_p29021_3", "new_p29021_0", "new_p29021_2", "new_p29032_0",
    "new_p29032_2", "new_p29032_1", "new_p29038_2", "new_p29038_1", "new_p29038_3",
    "new_p29038_0", "new_p29039_5", "new_p29039_4", "new_p29039_3", "new_p29039_1",
    "new_p29039_2", "new_p29039_7", "new_p29039_6", "new_p29051_3", "new_p29051_6",
    "new_p29051_1", "new_p29051_8", "new_p29051_0", "new_p29051_2", "new_p29051_4",
    "new_p29051_5", "new_p29051_7", "new_p29065_1", "new_p29065_4", "new_p29065_8",
    "new_p29065_11", "new_p29065_13", "new_p29065_0", "new_p29065_9", "new_p29065_10",
    "new_p29065_2", "new_p29065_6", "new_p29065_3", "new_p29065_12", "new_p29065_7",
    "new_p29065_5")
  
  available_columns <- intersect(selected_columns, names(shuzhen_df_out))
  shuzhen_category <- shuzhen_df_out[available_columns]
}

{
  shuzhen_category$f.20536.0.0[shuzhen_category$f.20536.0.0 %in% c(1, 2, 3)] <- 1
  shuzhen_category$f.20411.0.0[shuzhen_category$f.20411.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20405.0.0[shuzhen_category$f.20405.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20531.0.0[shuzhen_category$f.20531.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20529.0.0[shuzhen_category$f.20529.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20526.0.0[shuzhen_category$f.20526.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20530.0.0[shuzhen_category$f.20530.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20528.0.0[shuzhen_category$f.20528.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20527.0.0[shuzhen_category$f.20527.0.0 %in% c(1, 2)] <- 1
  shuzhen_category$f.20126.0.0[shuzhen_category$f.20126.0.0 %in% c(1, 2, 3,4,5)] <- 1
}
shuzhen_category$eid <- as.character(shuzhen_category$eid)
shuzhen_continues$eid <- as.character(shuzhen_continues$eid)

matchdata <- read.csv("./input/202category.csv")

#########################four brain age #########################

# Load exposure protein data
mydata <- read.csv("./input/brain_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         brain_difference = Bias_Corrected_Age - age)
mydata1 <- read.csv("./input/baizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         baizhi_difference = Bias_Corrected_Age - age)
mydata2 <- read.csv("./input/huizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         huizhi_difference = Bias_Corrected_Age - age)
head(mydata1)
head(mydata2)
# Extract required columns
baizhi <- mydata1[, c("eid", "baizhi_difference")]
huizhi <- mydata2[, c("eid", "huizhi_difference")]
# Merge by eid
mydata_merged <- merge(baizhi, huizhi, by = "eid", all = TRUE)
mydata <- mydata %>%
  left_join(mydata_merged, by = "eid")
mydata$age<-NULL
mydata$Bias_Corrected_Age<-NULL

# Merge all into final dataset
mydata <- left_join(mydata, cov1, by = "eid")
full_data <- left_join(mydata, shuzhen_continues, by = "eid")


# Define outcome and exposure variables
outcomes <- c("brain_difference", "baizhi_difference","huizhi_difference")
exposures <- setdiff(names(shuzhen_continues), "eid")

library(progress)
library(broom)
# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      df <- df |> filter(!is.na(.data[[expo]]) & !is.na(.data[[outcome]]) & !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    
    n <- nrow(df)  
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ ", expo, " + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      # Calculate standard deviations for standardization
      sd_expo <- sd(df[[expo]], na.rm = TRUE)
      sd_outcome <- sd(df[[outcome]], na.rm = TRUE)
      this_result <- tidy_fit |> 
        filter(term == expo) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          sd_exposure = sd_expo,
          sd_outcome = sd_outcome,
          n = n,
          beta_std = estimate * (sd_expo / sd_outcome)
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value,
               sd_exposure, sd_outcome, n,beta_std)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}

# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

final_result <- left_join(
  final_result, 
  matchdata, 
  by = c("exposure" = "outcome") 
)
#write.csv(final_result,"81+个表型与三种影像脑龄.csv",row.names = FALSE)



mydata <- read.csv("./input/Predicted_Age.csv")
mydata <- mydata |> 
  rename(eid = SampleID) |>  
  mutate(eid = as.character(eid))  
# Merge all into final dataset
mydata <- left_join(mydata, cov1, by = "eid")

mydata <- mydata %>%
  mutate(protein_difference = PredictedAge - age, .after = 2)  
full_data <- left_join(mydata, shuzhen_continues, by = "eid")

# Define outcome and exposure variables
outcomes <- c("protein_difference")  
exposures <- setdiff(names(shuzhen_continues), "eid")

library(progress)
library(broom)
# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      # Remove rows with missing values (only for current variable)
      df <- df |> filter(!is.na(.data[[expo]]) & !is.na(.data[[outcome]]) & !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    n <- nrow(df)  
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ ", expo, " + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      # Calculate standard deviations for standardization
      sd_expo <- sd(df[[expo]], na.rm = TRUE)
      sd_outcome <- sd(df[[outcome]], na.rm = TRUE)
      this_result <- tidy_fit |> 
        filter(term == expo) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          sd_exposure = sd_expo,
          sd_outcome = sd_outcome,
          n = n, 
          beta_std = estimate * (sd_expo / sd_outcome)
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value,
               sd_exposure, sd_outcome, n, beta_std)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}

# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

final_result <- left_join(
  final_result, 
  matchdata, 
  by = c("exposure" = "outcome")  
)
#write.csv(final_result,"./output/81+phenotype&brain_age.csv",row.names = FALSE)

###############################################
mydata <- mydata %>%
  mutate(protein_difference = PredictedAge - age, .after = 2)  
full_data <- left_join(mydata, shuzhen_category, by = "eid")
# Define outcome and exposure variables
outcomes <- c("protein_difference")  
exposures <- setdiff(names(shuzhen_category), "eid")  

library(progress)
library(broom)
library(dplyr)

# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      # Remove missing values: keep 0/1, drop 999/9999/NA
      df <- df |> 
        filter(.data[[expo]] %in% c(0, 1),  
               !is.na(.data[[outcome]]),
               !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    n <- nrow(df)  
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ factor(", expo, ") + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      this_result <- tidy_fit |> 
        filter(grepl(paste0("factor\\(", expo, "\\)"), term)) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          n = n, 
          reference_level = 0,  
          effect_level = 1      
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value, n,
               reference_level, effect_level)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}


# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

if(exists("matchdata")) {
  final_result <- left_join(
    final_result, 
    matchdata, 
    by = c("exposure" = "outcome")
  )
}
#write.csv(final_result,"121+个表型与蛋白影像脑龄.csv",row.names = FALSE)

####################################################
mydata <- read.csv("./input/brain_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         brain_difference = Bias_Corrected_Age - age)
mydata1 <- read.csv("./input/baizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         baizhi_difference = Bias_Corrected_Age - age)
mydata2 <- read.csv("./input/huizhi_age.csv") |> 
  select(eid, age, Bias_Corrected_Age) |> 
  mutate(eid = as.character(eid),
         huizhi_difference = Bias_Corrected_Age - age)
head(mydata1)
head(mydata2)
# Extract required columns
baizhi <- mydata1[, c("eid", "baizhi_difference")]
huizhi <- mydata2[, c("eid", "huizhi_difference")]
# Merge by eid
mydata_merged <- merge(baizhi, huizhi, by = "eid", all = TRUE)
mydata <- mydata %>%
  left_join(mydata_merged, by = "eid")
mydata$age<-NULL
mydata$Bias_Corrected_Age<-NULL

mydata <- left_join(mydata, cov1, by = "eid")
full_data <- left_join(mydata, shuzhen_category, by = "eid")
# Define outcome and exposure variables
outcomes <- c("brain_difference", "baizhi_difference","huizhi_difference")
exposures <- setdiff(names(shuzhen_category), "eid")  

library(progress)
library(broom)
library(dplyr)

# Initialize progress bar
pb <- progress_bar$new(total = length(exposures) * length(outcomes), 
                       format = "[:bar] :percent :elapsed", clear = FALSE)

results <- list()

# Loop through all exposure-outcome combinations
for (expo in exposures) {
  for (outcome in outcomes) {
    
    # Select relevant columns
    selected_columns <- c(expo, outcome, "age","sex","bmi","race","smoke","drink","APOEe4_carrier")
    selected_columns <- selected_columns[selected_columns %in% colnames(full_data)]
    df <- full_data |> select(eid, all_of(selected_columns))
    
    # Check if outcome exists in selected data
    if (outcome %in% colnames(df)) {
      df <- df |> 
        filter(.data[[expo]] %in% c(0, 1),  
               !is.na(.data[[outcome]]),
               !is.na(bmi))
    } else {
      pb$tick()
      next
    }
    
    n <- nrow(df) 
    
    # Create formula string for regression
    formula_str <- paste0(outcome, " ~ factor(", expo, ") + age + sex + bmi + race + smoke + drink + APOEe4_carrier")
    
    # Run linear regression with error handling
    fit <- tryCatch(
      lm(as.formula(formula_str), data = df),
      error = function(e) NULL
    )
    
    # Process results if model ran successfully
    if (!is.null(fit)) {
      tidy_fit <- tidy(fit)
      this_result <- tidy_fit |> 
        filter(grepl(paste0("factor\\(", expo, "\\)"), term)) |>  
        mutate(
          exposure = expo,
          outcome = outcome,
          n = n, 
          reference_level = 0,  
          effect_level = 1      
        ) |> 
        select(exposure, outcome, term, estimate, std.error, statistic, p.value, n,
               reference_level, effect_level)
      results[[length(results) + 1]] <- this_result
    }
    
    pb$tick()
  }
}


# Combine all results and calculate FDR-adjusted p-values
final_result <- bind_rows(results) |> 
  group_by(outcome) |> 
  mutate(p_fdr = p.adjust(p.value, method = "fdr")) |> 
  ungroup()

if(exists("matchdata")) {
  final_result <- left_join(
    final_result, 
    matchdata, 
    by = c("exposure" = "outcome")
  )
}
#write.csv(final_result,"./output/121+phenotype&brain_age.csv",row.names = FALSE)
