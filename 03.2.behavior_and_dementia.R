##############################dementia cox###########################
######################Load data########################
allcause_dementia <- read.csv("./input/allcause_dementia.csv") |> mutate(eid = as.character(eid))
F00 <- read.csv("F00.csv") |> mutate(eid = as.character(eid))
F01 <- read.csv("F01.csv") |> mutate(eid = as.character(eid))
G30 <- read.csv("G30.csv") |> mutate(eid = as.character(eid))

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
mydata <- left_join(allcause_dementia, cov1, by = "eid")
mydata00 <- left_join(F00, cov1, by = "eid")
mydata01 <- left_join(F01, cov1, by = "eid")
mydata30 <- left_join(G30, cov1, by = "eid")

mydata$X <- NULL
mydata30$X <- NULL

{
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
}

full_data_cat <- left_join(mydata, shuzhen_category, by = "eid")
full_data00_cat <- left_join(mydata00, shuzhen_category, by = "eid")
full_data01_cat <- left_join(mydata01, shuzhen_category, by = "eid")
full_data30_cat <- left_join(mydata30, shuzhen_category, by = "eid")
full_data_con <- left_join(mydata, shuzhen_continues, by = "eid")
full_data00_con <- left_join(mydata00, shuzhen_continues, by = "eid")
full_data01_con <- left_join(mydata01, shuzhen_continues, by = "eid")
full_data30_con <- left_join(mydata30, shuzhen_continues, by = "eid")

############### Loop Analysis for Cox Proportional Hazards Models ##############
# Load required libraries for survival analysis and data manipulation
library(survival)   # Core library for Cox regression and survival analysis
library(tidyverse)  # For data manipulation (dplyr) and visualization (ggplot2)

# Define covariates for adjustment in Cox models
# These are confounders to be controlled in the regression analysis
covariates <- c("age", "sex", "race", "bmi", "smoke", "drink", "APOEe4_carrier")

# Initialize empty list to store results across all datasets
all_results <- list()

# Define lists of dataframe names for categorical and continuous variables
# Categorical variable datasets (different dementia subtypes + categorical phenotypes)
cat_dataframes <- c("full_data_cat", "full_data00_cat", "full_data01_cat", "full_data30_cat")

# Continuous variable datasets (different dementia subtypes + continuous phenotypes)
con_dataframes <- c("full_data_con", "full_data00_con", "full_data01_con", 
                    "full_data30_con")

# Validate existence of all required dataframes (prevent runtime errors)
# Identify missing categorical dataframes
missing_cat_dfs <- setdiff(cat_dataframes, ls())
# Identify missing continuous dataframes
missing_con_dfs <- setdiff(con_dataframes, ls())

# Stop execution if categorical dataframes are missing
if (length(missing_cat_dfs) > 0) {
  stop(paste("Missing categorical dataframes:", paste(missing_cat_dfs, collapse = ", ")))
}
# Stop execution if continuous dataframes are missing
if (length(missing_con_dfs) > 0) {
  stop(paste("Missing continuous dataframes:", paste(missing_con_dfs, collapse = ", ")))
}

# Process continuous variable datasets (Cox regression for each variable)
for (df_name in con_dataframes) {
  
  # Retrieve the current dataframe by name from the global environment
  df <- get(df_name)
  
  # Define column range for independent variables (continuous phenotypes)
  var_cols <- 12:92   # Columns 12-92 contain continuous independent variables
  na_values <- c(NA)  # Only remove rows with NA values (no other missing codes)
  
  # Extract names of the continuous independent variables
  var_names <- names(df)[var_cols]
  
  # Initialize empty dataframe to store results for current dataset
  # Columns include key metrics from Cox regression output
  df_results <- data.frame(
    Dataset = character(),    # Name of the dataset being analyzed
    Variable = character(),   # Name of the independent variable
    N = integer(),            # Sample size after removing missing values
    Coef = numeric(),         # Regression coefficient (log hazard ratio)
    SE = numeric(),           # Standard error of the coefficient
    HR = numeric(),           # Hazard ratio (exponentiated coefficient)
    CI_lower = numeric(),     # Lower bound of 95% confidence interval for HR
    CI_upper = numeric(),     # Upper bound of 95% confidence interval for HR
    p_value = numeric(),      # Unadjusted p-value for the variable
    pFDR = numeric(),         # FDR-adjusted p-value (to control multiple testing)
    stringsAsFactors = FALSE  # Disable factor conversion for character columns
  )
  
  # Iterate over each continuous independent variable for Cox regression
  for (var in var_names) {
    # Filter rows: keep only non-missing values for the current variable
    # Preserve data for other variables to maximize sample size per analysis
    rows_to_keep <- !(df[[var]] %in% na_values)
    df_clean <- df[rows_to_keep, ]
    
    # Check remaining sample size (skip if no valid observations)
    n <- sum(rows_to_keep)
    if (n == 0) next  # Skip variable if all rows are removed
    
    # Build Cox model formula: Survival outcome ~ target variable + covariates
    # Surv(time, status) = survival object (time = follow-up time, status = event indicator)
    formula <- as.formula(paste("Surv(time, status) ~", var, "+", paste(covariates, collapse = "+")))
    
    # Fit Cox proportional hazards model with error handling
    # Return NULL if model fitting fails (e.g., separation, insufficient events)
    fit <- tryCatch(
      coxph(formula, data = df_clean),
      error = function(e) NULL
    )
    
    # Skip to next variable if model fitting failed
    if (is.null(fit)) next
    
    # Extract detailed model summary for result extraction
    summary_fit <- summary(fit)
    
    # Locate rows corresponding to the target variable in summary output
    coef_row <- which(rownames(summary_fit$coefficients) == var)
    confint_row <- which(rownames(summary_fit$conf.int) == var)
    
    # Append results to the current dataset's results dataframe
    df_results <- rbind(df_results, data.frame(
      Dataset = df_name,
      Variable = var,
      N = n,
      Coef = summary_fit$coefficients[coef_row, "coef"],          # Log HR
      SE = summary_fit$coefficients[coef_row, "se(coef)"],       # SE of log HR
      HR = summary_fit$conf.int[confint_row, "exp(coef)"],       # Hazard ratio
      CI_lower = summary_fit$conf.int[confint_row, "lower .95"], # 95% CI lower bound
      CI_upper = summary_fit$conf.int[confint_row, "upper .95"], # 95% CI upper bound
      p_value = summary_fit$coefficients[coef_row, "Pr(>|z|)"],  # Unadjusted p-value
      stringsAsFactors = FALSE
    ))
  }
  
  # Calculate FDR-adjusted p-values (Benjamini-Hochberg method)
  # Adjust p-values separately for each dataset to control false discovery rate
  if (nrow(df_results) > 0) {
    df_results$pFDR <- p.adjust(df_results$p_value, method = "fdr")
  }
  
  # Store results for current dataset in the global results list
  all_results[[df_name]] <- df_results
}

# Combine results from all continuous variable datasets into a single dataframe
# Bind rows to create a consolidated results table
final_results <- bind_rows(all_results)


# Load required libraries
# survival: Core package for Cox proportional hazards regression
# broom: Tidies model output into data frames for easy extraction
# tidyverse: Includes dplyr for data manipulation (filter/bind_rows/left_join)
library(survival)
library(broom)
library(tidyverse)

# Define list of covariates for adjustment in Cox regression models
# These are confounders to control for in the analysis (demographics/bio markers)
covariates <- c("age", "sex", "race", "bmi", "smoke", "drink", "APOEe4_carrier")

# Define list of dataframe names for categorical variable analysis
# Each dataframe corresponds to a dementia subtype + categorical phenotype data
df_names <- c("full_data_cat", "full_data00_cat", "full_data01_cat", "full_data30_cat")

# Initialize empty list to store results across all categorical datasets
all_results1 <- list()

# Iterate over each categorical dataset for Cox regression analysis
for (df_name in df_names) {
  # Retrieve the current dataframe by name from the global environment
  df <- get(df_name)
  
  # Initialize empty dataframe to store results for the current dataset
  current_results <- data.frame()
  
  # Define column range for independent variables (categorical phenotypes)
  # Columns 12-132 contain the categorical independent variables
  var_cols <- 12:132
  var_names <- colnames(df)[var_cols]
  
  # Iterate over each categorical independent variable
  for (var in var_names) {
    # Create temporary dataframe with missing values removed for current variable
    # Exclude coded missing values (999, 9999) and NA values
    # !!sym(var) converts string variable name to a symbol for dplyr filter
    temp_df <- df %>%
      filter(!(!!sym(var) %in% c(999, 9999)), 
             !is.na(!!sym(var)))
    
    # Check for sufficient data diversity (at least 2 unique groups required)
    # Skip variable if only one category remains (no variability to test)
    if (length(unique(temp_df[[var]])) < 2) {
      next
    }
    
    # Build formula for Cox proportional hazards model
    # Surv(time, status) = survival object (time = follow-up time, status = event indicator)
    # Formula structure: Survival outcome ~ target variable + adjustment covariates
    formula <- as.formula(paste("Surv(time, status) ~", var, "+", 
                                paste(covariates, collapse = " + ")))
    
    # Fit Cox model with error handling (return NULL if model fails to fit)
    # Common failures: perfect separation, insufficient events, collinearity
    cox_model <- tryCatch(
      {
        coxph(formula, data = temp_df)
      },
      error = function(e) {
        return(NULL)
      }
    )
    
    # Extract results only if model fitting was successful
    if (!is.null(cox_model)) {
      # Get sample size after removing missing values for current variable
      n <- nrow(temp_df)
      
      # Tidy model output using broom (returns standardized dataframe)
      # conf.int = TRUE: Include 95% confidence intervals
      # exponentiate = TRUE: Return hazard ratios (instead of log HR)
      model_summary <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
      
      # Extract results for the target variable (always first row in summary)
      # Covariates appear in subsequent rows and are not extracted here
      var_result <- model_summary[1, ]
      
      # Create a single row of results for the current variable
      result_row <- data.frame(
        dataset = df_name,          # Name of the dataset being analyzed
        variable = var,             # Name of the categorical independent variable
        n = n,                      # Sample size for this analysis
        coef = var_result$estimate, # Hazard ratio (exponentiated coefficient)
        se = var_result$std.error,  # Standard error of the hazard ratio
        HR = var_result$estimate,   # Hazard ratio (duplicated for consistency)
        HR_95CI_lower = var_result$conf.low,  # Lower bound of 95% CI for HR
        HR_95CI_upper = var_result$conf.high, # Upper bound of 95% CI for HR
        p_value = var_result$p.value,         # Unadjusted p-value for the variable
        stringsAsFactors = FALSE    # Disable automatic factor conversion
      )
      
      # Append current variable results to the dataset's results dataframe
      current_results <- bind_rows(current_results, result_row)
    }
  }
  
  # Calculate FDR-adjusted p-values (Benjamini-Hochberg method)
  # Adjust p-values per dataset to control false discovery rate
  if (nrow(current_results) > 0) {
    current_results$fdr_p_value <- p.adjust(current_results$p_value, method = "fdr")
  }
  
  # Store results for the current dataset in the global results list
  all_results1[[df_name]] <- current_results
}

# Combine results from all categorical datasets into a single dataframe
# .id = "dataset": Add a column indicating the source dataset name
final_results1 <- bind_rows(all_results1, .id = "dataset")

################# Merge Categorical and Continuous Results ###############
# Standardize column names of categorical results to match continuous results
# Ensures compatibility for row-binding (same column structure)
colnames(final_results1) <- colnames(final_results)

# Merge categorical and continuous results into one consolidated dataframe
final_combined <- rbind(final_results, final_results1)

# Read matching dataset (phenotype name mapping/annotation file)
# File contains 203 emotional indicators (6 categories, duplicates removed)
matchdata <- read.csv("./input/202category.csv")

# Merge results with phenotype annotation data
# Left join: Keep all rows from final_combined, match by Variable = outcome
# Preserves all analysis results even if no match is found in matchdata
final_combined <- left_join(
  final_combined, 
  matchdata, 
  by = c("Variable" = "outcome")  # Specify join columns (Variable in results = outcome in matchdata)
)


######################################################
# Merge all into final dataset
mydata <- left_join(allcause_dementia, cov1, by = "eid")
# Replace "allcause_dementia" with specific dementia datasets (F00, F01, G30) for each analysis

{
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
}

outcome_data <- outcome_data %>% 
  mutate(eid = as.character(eid))
full_data <- left_join(mydata, outcome_data, by = "eid")

# Load required libraries
library(progress)
library(survival)
library(dplyr)

# Define key variables for survival analysis
time_var <- "time"      # Column name for follow-up time (survival time)
status_var <- "status"  # Column name for event status (1=event occurred, 0=censored)
# Extract exposure variables: exclude ID column (eid) and irrelevant column (X)
exposures <- setdiff(names(outcome_data), c("eid", "X"))  

# Initialize progress bar (track progress of exposure variable loop)
pb <- progress_bar$new(
  total = length(exposures),  # Total number of exposure variables to process
  format = "[:bar] :percent :elapsed",  # Progress bar display format
  clear = FALSE  # Keep progress bar visible after completion
)

# Initialize empty list to store Cox model results for each exposure
results <- list()

# Iterate over each exposure variable for Cox regression analysis
for (expo in exposures) {
  # Define required variables for the model (exposure + survival + covariates)
  required_vars <- c(expo, time_var, status_var, "age", "sex", "race", "bmi", "smoke","drink","APOEe4_carrier")
  
  # Skip variable if any required column is missing from the dataset
  if (!all(required_vars %in% names(full_data))) next
  
  # Filter dataset to retain only required variables and remove missing values
  # Keep complete cases for exposure, time, and status variables
  df <- full_data |> 
    select(all_of(required_vars)) |>  # Select only required columns
    filter(
      !is.na(!!sym(expo)),          # Remove missing values for current exposure
      !is.na(!!sym(time_var)),      # Remove missing values for survival time
      !is.na(!!sym(status_var))     # Remove missing values for event status
    )
  
  # Skip variable if no valid observations remain after filtering
  if (nrow(df) == 0) next
  
  # Build formula for Cox proportional hazards model
  # Formula structure: Surv(time, status) ~ exposure + adjustment covariates
  formula_str <- paste(
    "Surv(", time_var, ", ", status_var, ") ~",
    expo, "+ age + sex + race + bmi + smoke + drink + APOEe4_carrier"
  )
  
  # Fit Cox model with error handling (return NULL if model fails to converge)
  # Common failures: perfect separation, insufficient events, collinearity
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = df),
    error = function(e) NULL
  )
  
  # Extract and store results only if model fitting was successful
  if (!is.null(fit)) {
    # Tidy model output (extract coefficients, p-values, confidence intervals)
    # Exponentiate results to get hazard ratios (HR) instead of log HR
    tidy_fit <- broom::tidy(fit, conf.int = TRUE) %>%
      filter(term == expo) %>%  # Keep only results for the target exposure
      mutate(
        exposure = expo,        # Add exposure variable name for tracking
        n_samples = nrow(df),   # Record sample size for this analysis
        hazard_ratio = exp(estimate),  # Calculate hazard ratio (HR)
        conf.low = exp(conf.low),      # Lower bound of 95% CI for HR
        conf.high = exp(conf.high)     # Upper bound of 95% CI for HR
      )
    
    # Store tidy results in the results list (named by exposure variable)
    results[[expo]] <- tidy_fit
  }
  
  pb$tick()  # Update progress bar (increment by 1 for each exposure processed)
}

# Combine all individual exposure results into a single dataframe
final_result <- bind_rows(results)

# Calculate FDR-adjusted p-values (Benjamini-Hochberg method)
# Control false discovery rate for multiple testing correction
final_result$p_fdr <- p.adjust(final_result$p.value, method = "fdr")

# Check for any warnings generated during model fitting/processing
warnings()
warnings()

# ========================== Result Processing & Annotation ==========================

# Define variable category lists (extract exposure variables by dataset, exclude ID)
fi_vars <- setdiff(names(Food_Intake), "eid")          # Food Intake variables
bm_vars <- setdiff(names(Body_measure), "eid")        # Body Measurement variables
li_vars <- setdiff(names(liking), "eid")              # Food Liking variables
pf_vars <- setdiff(names(Psychosocial_factors), "eid")# Psychosocial Factors variables
ls_vars <- setdiff(names(Lifestyles), "eid")          # Lifestyles variables
le_vars <- setdiff(names(Local_environment), "eid")   # Local Environment variables
bio_vars <- setdiff(names(bio), "eid")                # Biochemical Markers variables
fix_vars <- setdiff(names(file_cox_all), "eid")       # file_cox_all variables

# Assign categorical labels to each exposure variable
# Classify variables into predefined categories for downstream analysis/visualization
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

# Load readxl library to import Excel mapping data
library(readxl)
# Load additional variable mapping data (e.g., variable labels/descriptions)
data1 <- read_excel("./input/match_data.xlsx", sheet = "Sheet1")

# Standardize column name for merging (correct typo: "exprosure" → "exposure")
names(data1)[names(data1) == "exprosure"] <- "exposure"

# Merge Cox regression results with mapping data (left join to preserve all results)
# Add descriptive labels/annotations to exposure variables
phewas_data <- left_join(final_result, data1, by = "exposure")

# Fill missing exposure names: use original exposure name if mapped name is NA
# Ensure all variables have a readable name (no missing labels)
phewas_data <- phewas_data %>%
  mutate(exprosure1 = if_else(is.na(exprosure1), exposure, exprosure1))

# Check number of missing values in the exposure name column (quality control)
sum(is.na(phewas_data$exprosure1))
