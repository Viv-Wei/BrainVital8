#################
## 代码结束 — 逻辑和方法论已验证
################
图书馆（生存）
库（dplyr）
图书馆（阅读器）

##---------------
##-------------
##-------------
df <- read_csv("./输入/最终_cox_合并_处理2。csv”)

##-------------
## 2。数据预处理
##-------------

# 排除基线痴呆症
df <- df %>%
 过滤器（失智_状态！= "基线”)

# 构建生存时间和结果（将天数转换为年数）
df <- df %>%
 变异(
 time_days = dementia_to_jinzu_days, # 生存时间（天）
 时间_年 = 痴呆症_到_津祖_天 / 365。25,  # 生存时间（年）
 状态 = 组合_诊断 # 1=直呆症,0=直查
  )

# 确保收入是一个因素
df <- df %>%
 变异(
 收入 = 因素（收入,水平 = c(1, 2, 3, 4))
  )

# 一 BrainVital8 一千零八千零八（Q1 四千,Q4 四千）
df <- df %>%
 变异(
 BrainVital8_四分位数 = 切割(
 BrainVital8_脑活,
 breaks = brainVital8_aligned,probs = c(0, 0。25, 0。5, 0。75, 1),娜。rm = 真实的）,
 标签 = c("Q1”, "Q2”, "Q3”, "Q4”),
 包括。最低= 真实的
    )
  )

##-------------
## 3。添加组样本统计数据,包括 Q1 参考组
##-------------
run_cox_analysis_with_groups <- 功僽（数据、公式_str、分析_名称、 
 required_vars = 必填项, 必填_必填= "一者”) {
  
  # 如果指定了所需的变量,则过滤没有缺失值的样本
  如果 (！是。null（required_vars）){
 data_filtered <- 数据
    为了 （变量 在 所需_变量）{
 数据_已过滤 <- 数据_已过滤 %>% 
 筛选(！是。na(！！sym（var）))
    }
  } 得的 {
 data_filtered <- 数据
  }
  
  # 计算总样本量和病例数
  total_n <- nrow（数据_已过滤）
  total_cases <- sum（data_filtered$status == 1,na。rm = TRUE）
  
  # 计算每个 BrainVital8 四分位数的样本量和病例数
  group_stats <- 数据_已过滤 %>%
    group_by（BrainVital8_四分位数）%>%
    总结(
      Group_N = n(),
      Group_Case = sum（状态 == 1,na。rm = TRUE）,
      . 。。组 = 'drop'
    ) %>%
    脑活（BrainVital8_脑活）
  
  猫("\n不析:", 分析_名称）
  猫("\n 总计 N:", 总计_n）
  猫("\n 总案例数:", 总计_案例）
  猫("\n 按组划分的样本量:")
  对于（i 在 1:nrow（group_stats）中）{
    猫（粘贴0("\n ", group_stats$BrainVital8_quartile[i], ": N=", 
               group_stats$Group_N[i], ", 案例=", group_stats$Group_Case[i]))
  }
  猫("\n")
  
  # 检查样本量是否足够
  如果（total_n < 10 || total_cases < 5){
    警告（粘贴("分析", 分析_名称, "样本量不足"))
    return(NULL)
  }
  
  # Run Cox model
  model <- coxph(as.formula(formula_str), data = data_filtered)
  sm <- summary(model)
  
  # Extract results for all variables
  results_full <- data.frame(
    Analysis = analysis_name,
    Variable = rownames(sm$coefficients),
    HR       = exp(sm$coefficients[, "coef"]),
    CI_lower = exp(sm$coefficients[, "coef"] - 1.96 * sm$coefficients[, "se(coef)"]),
    CI_upper = exp(sm$coefficients[, "coef"] + 1.96 * sm$coefficients[, "se(coef)"]),
    P_value  = sm$coefficients[, "Pr(>|z|)"],
    N        = total_n,           # Total sample size
    Case     = total_cases,       # Total cases
    stringsAsFactors = FALSE
  )
  
  # Extract BrainVital8 related variables
  brainvital8_vars <- c("BrainVital8_aligned", "BrainVital8_quartileQ2", 
                        "BrainVital8_quartileQ3", "BrainVital8_quartileQ4")
  
  # Create a dataframe for all BrainVital8 variables (incl. Q1 ref)
  brainvital8_results <- data.frame()
  
  # 1. Add continuous variable results
  continuous_row <- results_full %>% 
    filter(grepl("BrainVital8_aligned", Variable))
  
  if (nrow(continuous_row) > 0) {
    continuous_row$Group_N <- total_n
    continuous_row$Group_Case <- total_cases
    brainvital8_results <- bind_rows(brainvital8_results, continuous_row)
  }
  
  # 2. Add categorical variable results (incl. Q1 ref)
  
  # First find statistics for Q1
  q1_stats <- group_stats %>% filter(BrainVital8_quartile == "Q1")
  q1_n <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_N[1], NA)
  q1_case <- ifelse(nrow(q1_stats) > 0, q1_stats$Group_Case[1], NA)
  
  # Add Q1 as the reference group
  q1_row <- data.frame(
    Analysis = analysis_name,
    Variable = "BrainVital8_quartileQ1 (ref)",
    HR = 1.00,
    CI_lower = NA,
    CI_upper = NA,
    P_value = NA,
    N = total_n,
    Case = total_cases,
    Group_N = q1_n,
    Group_Case = q1_case,
    stringsAsFactors = FALSE
  )
  
  brainvital8_results <- bind_rows(brainvital8_results, q1_row)
  
  # Add Q2, Q3, Q4
  for (quartile in c("Q2", "Q3", "Q4")) {
    var_name <- paste0("BrainVital8_quartile", quartile)
    
    # Extract results for this variable from Cox output
    quartile_row <- results_full %>% filter(Variable == var_name)
    
    # Extract statistics for this group from group_stats
    quartile_stats <- group_stats %>% filter(BrainVital8_quartile == quartile)
    quartile_n <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_N[1], NA)
    quartile_case <- ifelse(nrow(quartile_stats) > 0, quartile_stats$Group_Case[1], NA)
    
    if (nrow(quartile_row) > 0) {
      quartile_row$Group_N <- quartile_n
      quartile_row$Group_Case <- quartile_case
      brainvital8_results <- bind_rows(brainvital8_results, quartile_row)
    } else {
      # If variable is not in Cox model (e.g. small sample), still add group info
      missing_row <- data.frame(
        Analysis = analysis_name,
        Variable = var_name,
        HR = NA,
        CI_lower = NA,
        CI_upper = NA,
        P_value = NA,
        N = total_n,
        Case = total_cases,
        Group_N = quartile_n,
        Group_Case = quartile_case,
        stringsAsFactors = FALSE
      )
      brainvital8_results <- bind_rows(brainvital8_results, missing_row)
    }
  }
  
  # Reorder: Continuous variable, then Q1-Q4
  brainvital8_results <- brainvital8_results %>%
    mutate(order = case_when(
      Variable == "BrainVital8_aligned" ~ 1,
      Variable == "BrainVital8_quartileQ1 (ref)" ~ 2,
      Variable == "BrainVital8_quartileQ2" ~ 3,
      Variable == "BrainVital8_quartileQ3" ~ 4,
      Variable == "BrainVital8_quartileQ4" ~ 5,
      TRUE ~ 6
    )) %>%
    arrange(order) %>%
    select(-order)
  
  # PH test
  ph_test <- cox.zph(model)
  
  return(list(
    results_full = results_full,
    results = brainvital8_results,
    ph_test = ph_test,
    model   = model,
    data_used = data_filtered,
    group_stats = group_stats,
    model_type = model_type
  ))
}

##----------------------------------------------------------
## 4. Create function for both Cox models
##----------------------------------------------------------
run_both_models <- function(data, analysis_name, required_vars) {
  
  # Model 1: Using continuous variable BrainVital8_aligned
  formula_continuous <- "Surv(time_years, status) ~ BrainVital8_aligned + baseline_age + gender + bmi + education + smoking + drinking + income + hypertension + diabates + depression"
  
  res_continuous <- run_cox_analysis_with_groups(
    data = data,
    formula_str = formula_continuous,
    analysis_name = paste0(analysis_name, "_continuous"),
    required_vars = required_vars,
    model_type = "continuous"
  )
  
  # Model 2: Using categorical variable BrainVital8_quartile
  formula_categorical <- "Surv(time_years, status) ~ BrainVital8_quartile + baseline_age + gender + bmi + education + smoking + drinking + income + hypertension + diabates + depression"
  
  res_categorical <- run_cox_analysis_with_groups(
    data = data,
    formula_str = formula_categorical,
    analysis_name = paste0(analysis_name, "_categorical"),
    required_vars = required_vars,
    model_type = "categorical"
  )
  
  return(list(
    continuous = res_continuous,
    categorical = res_categorical
  ))
}

##----------------------------------------------------------
## 5. Main analysis (Total population) - Two models
##----------------------------------------------------------
cat("=== Starting sample screening for each analysis ===\n")

# Variables required for main analysis
main_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
               "baseline_age", "gender", "bmi", "education", "smoking", "drinking",
               "income", "hypertension", "diabates", "depression")

res_main_both <- run_both_models(
  data = df,
  analysis_name = "Main",
  required_vars = main_vars
)

res_main_continuous <- res_main_both$continuous
res_main_categorical <- res_main_both$categorical

##----------------------------------------------------------
## 6. Sex-stratified analysis - Two models
##----------------------------------------------------------

# Variables for male analysis (excluding gender)
male_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
               "baseline_age", "bmi", "education", "smoking", "drinking",
               "income", "hypertension", "diabates", "depression")

# Filter male samples first
df_male <- df %>% 
  filter(gender == 1) %>%
  filter(!is.na(gender))

res_male_both <- run_both_models(
  data = df_male,
  analysis_name = "Male",
  required_vars = male_vars
)

res_male_continuous <- res_male_both$continuous
res_male_categorical <- res_male_both$categorical

# Female analysis
df_female <- df %>% 
  filter(gender == 2) %>%
  filter(!is.na(gender))

res_female_both <- run_both_models(
  data = df_female,
  analysis_name = "Female",
  required_vars = male_vars
)

res_female_continuous <- res_female_both$continuous
res_female_categorical <- res_female_both$categorical

##----------------------------------------------------------
## 6.5 Age-stratified analysis - Two models
##----------------------------------------------------------

# Age stratification (<65 years)
df_age_lt65 <- df %>% 
  filter(baseline_age < 65) %>%
  filter(!is.na(baseline_age))

# Note: baseline_age is no longer included as a covariate here
age_strat_vars <- c("time_years", "status", "BrainVital8_aligned", "BrainVital8_quartile",
                    "gender", "bmi", "education", "smoking", "drinking",
                    "income", "hypertension", "diabates", "depression")

# Analysis for Age < 65
res_age_lt65_both <- run_both_models(
  data = df_age_lt65,
  analysis_name = "Age_lt65",
  required_vars = age_strat_vars
)

res_age_lt65_continuous <- res_age_lt65_both$continuous
res_age_lt65_categorical <- res_age_lt65_both$categorical

# Analysis for Age >= 65
df_age_ge65 <- df %>% 
  filter(baseline_age >= 65) %>%
  filter(!is.na(baseline_age))

res_age_ge65_both <- run_both_models(
  data = df_age_ge65,
  analysis_name = "Age_ge65",
  required_vars = age_strat_vars
)

res_age_ge65_continuous <- res_age_ge65_both$continuous
res_age_ge65_categorical <- res_age_ge65_both$categorical

##----------------------------------------------------------
## 7. Exclude cases within 2 years - Two models
##----------------------------------------------------------

if (!is.null(res_main_continuous)) {
  df_ex2yr <- res_main_continuous$data_used %>%
    filter(!(status == 1 & time_years < 2))
  
  res_ex2yr_both <- run_both_models(
    data = df_ex2yr,
    analysis_name = "Exclude_2yr",
    required_vars = main_vars
  )
  
  res_ex2yr_continuous <- res_ex2yr_both$continuous
  res_ex2yr_categorical <- res_ex2yr_both$categorical
}

##----------------------------------------------------------
## 8. Analysis including LIBRA2 - Two models
##----------------------------------------------------------
libra2_vars <- c(main_vars, "LIBRA2")

res_libra2_both <- run_both_models(
  data = df,
  analysis_name = "LIBRA2",
  required_vars = libra2_vars
)

res_libra2_continuous <- res_libra2_both$continuous
res_libra2_categorical <- res_libra2_both$categorical

##----------------------------------------------------------
## 9. Analysis including Lancet - Two models
##----------------------------------------------------------
lancet_vars <- c(main_vars, "lancet")

res_lancet_both <- run_both_models(
  data = df,
  analysis_name = "Lancet",
  required_vars = lancet_vars
)

res_lancet_continuous <- res_lancet_both$continuous
res_lancet_categorical <- res_lancet_both$categorical

##----------------------------------------------------------
## 10. Merge and save results
##----------------------------------------------------------
# Collect all valid results
all_results <- list()

# Add continuous model results
if (!is.null(res_main_continuous)) all_results[["Main_continuous"]] <- res_main_continuous
if (!is.null(res_male_continuous)) all_results[["Male_continuous"]] <- res_male_continuous
if (!is.null(res_female_continuous)) all_results[["Female_continuous"]] <- res_female_continuous
# Add age stratification results
if (!is.null(res_age_lt65_continuous)) all_results[["Age_lt65_continuous"]] <- res_age_lt65_continuous
if (!is.null(res_age_ge65_continuous)) all_results[["Age_ge65_continuous"]] <- res_age_ge65_continuous
if (!is.null(res_ex2yr_continuous)) all_results[["Exclude_2yr_continuous"]] <- res_ex2yr_continuous
if (!is.null(res_libra2_continuous)) all_results[["LIBRA2_continuous"]] <- res_libra2_continuous
if (!is.null(res_lancet_continuous)) all_results[["Lancet_continuous"]] <- res_lancet_continuous

# Add categorical model results
if (!is.null(res_main_categorical)) all_results[["Main_categorical"]] <- res_main_categorical
if (!is.null(res_male_categorical)) all_results[["Male_categorical"]] <- res_male_categorical
if (!is.null(res_female_categorical)) all_results[["Female_categorical"]] <- res_female_categorical
# Add age stratification results
if (!is.null(res_age_lt65_categorical)) all_results[["Age_lt65_categorical"]] <- res_age_lt65_categorical
if (!is.null(res_age_ge65_categorical)) all_results[["Age_ge65_categorical"]] <- res_age_ge65_categorical
if (!is.null(res_ex2yr_categorical)) all_results[["Exclude_2yr_categorical"]] <- res_ex2yr_categorical
if (!is.null(res_libra2_categorical)) all_results[["LIBRA2_categorical"]] <- res_libra2_categorical
if (!is.null(res_lancet_categorical)) all_results[["Lancet_categorical"]] <- res_lancet_categorical

# Merge BrainVital8 results (including Q1 ref)
brainvital8_results <- data.frame()

for (analysis_name in names(all_results)) {
  res <- all_results[[analysis_name]]
  brainvital8_results <- bind_rows(brainvital8_results, res$results)
}

# Merge full results (excluding group statistics)
full_results <- data.frame()

for (analysis_name in names(all_results)) {
  res <- all_results[[analysis_name]]
  full_results <- bind_rows(full_results, res$results_full)
}

##----------------------------------------------------------
## 11. Extract and save Global PH test results for both models
##----------------------------------------------------------
extract_global_ph_test_results <- function(results_list) {
  ph_results <- data.frame()
  
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    
    if (!is.null(res$ph_test)) {
      # Extract PH test results
      ph_table <- res$ph_test$table
      
      # Determine model type
      model_type <- ifelse(grepl("_continuous$", analysis_name), "continuous", "categorical")
      model_label <- ifelse(model_type == "continuous", 
                            "BrainVital8_aligned (continuous)", 
                            "BrainVital8_quartile (categorical)")
      
      # Extract global PH test results
      if ("GLOBAL" %in% rownames(ph_table)) {
        global_ph <- ph_table["GLOBAL", ]
        global_ph_df <- data.frame(
          Analysis = gsub("_(continuous|categorical)$", "", analysis_name),
          Model_Type = model_type,
          Model_Label = model_label,
          Global_chisq = global_ph["chisq"],
          Global_df = global_ph["df"],
          Global_p = global_ph["p"],
          stringsAsFactors = FALSE
        )
        ph_results <- bind_rows(ph_results, global_ph_df)
      }
    }
  }
  
  return(ph_results)
}

# Extract Global PH test results for all analyses
global_ph_results <- extract_global_ph_test_results(all_results)

##----------------------------------------------------------
## 12. Extract and save variable-specific PH test results
##----------------------------------------------------------
extract_variable_specific_ph_results <- function(results_list) {
  ph_results <- data.frame()
  
  for (analysis_name in names(results_list)) {
    res <- results_list[[analysis_name]]
    
    如果(！is。null（res$ph_test）) {
      # 提取 PH 测试结果
      ph_table <- res$ph_test$table
      
      # 确定模型类型
      model_type <- ifelse（grepl("_连续$", 分析_名称）, "连续的", “分类")
      
      # 脑Vital8 脑Vital8 脑Vital8 PH 脑Vital8 PH
 brainvital8_vars <- grep("BrainVital8”, rownames（ph_table）,子 = 真实的）
      
 智（智（brainvital8_vars）> 0){
 对于（brainvital8_vars 中元 var）{
 var_ph <- ph_table[var,]
 var_ph_df <- 数据框(
 分析=gsub("_（连续的|分类）$", “", 分析_名称）,
            模型_类型 = 模型_类型,
            变量 = var,
            chisq = var_ph["奇斯克"],
            df = var_ph["df"],
            p = var_ph["p"],
            stringsAsFactors = FALSE
          )
          ph_results <- bind_rows（ph_results,var_ph_df）
        }
      }
    }
  }
  
  返回（ph_结果）
}

# PH 测试结果
变量_ph_结果 <- 提取_变量_特定_ph_结果（所有_结果）

##--------------------------------
## 13。保存结果
##--------------------------------

# 保存全局 PH 测试结果
写入。csv(
 global_ph_results,
  ". 。/大脑/BrainVital8_Global_PH_test_results。csv”,
 行。得 = 得
)

# 保存特定于变量的 PH 测试结果
写入。csv(
 变量_ph_结果,
  ". 。/大脑/BrainVital8_变量_PH_测试_结果。csv”,
 行。得 = 得
)

# 保存 Cox 分析结果
写入。csv(
 brainvital8_脑果,
  ". 。。/BrainVital8_脑呆_Cox_脑果_with_group_stats_incl_Q1x。csv”,
 行。名称 = 错误的
)

##-------------------------------
## 14。创建摘要表
##-------------------------------
create_ph_summary_table <- 全局（global_results,variable_results）{
  # 结合全局结果和变量特定结果
 摘要_表 <- 数据。框架()
  
 对于（唯一（global_results$Analysis）一目一然）{
    # 获取此特定分析的所有结果
 analysis_global <- global_results %>% 分析结果（分析 == 分析）
 分析_变量 <- 变量_结果 %>% 过滤器（分析 == 分析）
    
    #Continuous 模型
 cont_global <- analysis_global %>% 模型类型（模型_类型 == "连续的")
 cont_variable <- 分析_变量 %>% 过滤器（模型_类型 == "连续的" & 变量 == "BrainVital8_大脑”)
    
 如果（nrow（cont_global）> 0){
 cont_row <- 数据行(
 分析=分析,
 模型="脑活（BrainVital8_aligned）",
 BrainVital8_Variable_PH = ifelse（nrow（cont_variable）> 0, 
 粘贴0("χ²=", round（cont_variable$chisq, 3), 
                                                ", df=", cont_variable$df, 
                                                ", p=", round（cont_variable$p, 4)),
                                         "不适用"),
 Global_PH = 粘贴0("χ²=", round（cont_global$Global_chisq, 3), 
                           ", df=", cont_global$Global_df, 
                           ", p=", round（cont_global$Global_p, 4)),
 字符串作为因素 = 错误的
      )
 summary_table <- bind_rows（summary_table,cont_row）
    }
    
    # 分类模型
 cat_global <- analysis_global %>% 分析了（分析_分析 =="分类”)
 cat_variables <- 分析_变量 %>% 过滤器（模型_类型 == "分类” & 格雷普尔("BrainVital8_脑力万能”, 变量）)
    
 如垜（nrow（cat_global）> 0){
      #Create PH PH PH PH PH
 var_ph_strings <- c()
 如果（nrow（cat_variables）> 0){
 对于（i 在 1:nrow（cat_variables）中）{
 var_name <- gsub("BrainVital8_脑力万能”, “", cat_variables$变量[i])
          var_ph_字符一 <- c（var_ph_字符一, 
                              粘名0（var_name,“: χ²=", round（cat_variables$chisq[i], 3),
                                     ", p=", round（cat_variables$p[i], 4)))
        }
        var_ph_combined <- 粘贴（var_ph_strings,折叠 = "; ")
      } 否则{
        var_ph_combined <- "不适用"
      }
      
      cat_row <- 数据框(
        分析=分析,
        模型="分类（BrainVital8_四分位数）",
        BrainVital8_Variable_PH = var_ph_combined,
        Global_PH = 粘贴0("χ²=", round（cat_global$Global_chisq, 3), 
                           ", df=", cat_global$Global_df, 
                           ", p=", round（cat_global$Global_p, 4)),
        stringsAsFactors = FALSE
      )
      summary_table <- bind_rows（summary_table,cat_row）
    }
  }
  
  返回（摘要_表）
}

# 创建汇总表
ph_summary_table <- create_ph_summary_table（global_ph_results,variable_ph_results）

# 保存汇总表
写入。csv(
  ph_摘要_表,
  ". ./输出/BrainVital8_PH_test_summary_table。csv",
  行。名称 = FALSE
)


