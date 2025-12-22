library(dplyr)
library(meta)
library(stringr)
library(readxl)
library(tidyr)
library(tibble)
library(tools)

# ===============================
# Read Excel file + automatically create file prefix
# ===============================
input_file <- "./input/meta.xlsx"
df <- read_excel(input_file)

# Automatically use input file name as output file prefix
file_base <- file_path_sans_ext(basename(input_file))

df <- df %>% rename(HR_CI = `HR (95% CI)`)

# ===============================
# Clean HR and calculate logHR / SE
# ===============================
df_clean <- df %>%
  filter(HR_CI != "Reference") %>%
  mutate(
    HR = as.numeric(str_extract(HR_CI, "^[0-9.]+")),
    lower_CI = as.numeric(str_extract(HR_CI, "(?<=\\()[0-9.]+")),
    upper_CI = as.numeric(str_extract(HR_CI, "(?<=, )[0-9.]+(?=\\))")),
    logHR = log(HR),
    SE = (log(upper_CI) - log(lower_CI)) / (2 * 1.96)
  )

# ===============================
# Meta-analysis
# ===============================
d_continuous <- df_clean %>% filter(Comparison == "Continuous")
meta_continuous <- metagen(
  TE = logHR,
  seTE = SE,
  studlab = Database,
  data = d_continuous,
  sm = "HR",
  method.tau = "DL",
  method.random.ci = "HK"
)

d_quartiles <- df_clean %>% filter(Comparison %in% c("Q2","Q3","Q4"))

meta_q2 <- metagen(
  TE = logHR, seTE = SE, studlab = Database,
  data = filter(d_quartiles, Comparison == "Q2"),
  sm = "HR", method.tau = "DL", method.random.ci = "HK"
)

meta_q3 <- metagen(
  TE = logHR, seTE = SE, studlab = Database,
  data = filter(d_quartiles, Comparison == "Q3"),
  sm = "HR", method.tau = "DL", method.random.ci = "HK"
)

meta_q4 <- metagen(
  TE = logHR, seTE = SE, studlab = Database,
  data = filter(d_quartiles, Comparison == "Q4"),
  sm = "HR", method.tau = "DL", method.random.ci = "HK"
)

# ===============================
# Extract pooled HR and 95% CI
# ===============================
get_meta_HR_CI <- function(meta_obj){
  data.frame(
    HR = exp(meta_obj$TE.random),
    CI_lower = exp(meta_obj$lower.random),
    CI_upper = exp(meta_obj$upper.random)
  )
}

df_meta_summary <- tibble(
  Comparison = c("Continuous", "Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1")
) %>% bind_cols(
  rbind(
    get_meta_HR_CI(meta_continuous),
    get_meta_HR_CI(meta_q2),
    get_meta_HR_CI(meta_q3),
    get_meta_HR_CI(meta_q4)
  )
) %>% mutate(
  HR_CI_formatted = sprintf("%.2f (%.2f, %.2f)", HR, CI_lower, CI_upper)
)

# ===============================
# ★ Correct total counts: include only Continuous, Q1, Q2, Q3, Q4
# ===============================
df_counts <- df %>%
  mutate(
    Cases = as.numeric(str_extract(`Cases/N`, "^[0-9]+")),
    N     = as.numeric(str_extract(`Cases/N`, "(?<=/)[0-9]+"))
  )

count_cont <- df_counts %>% filter(Comparison=="Continuous") %>% summarise(Cases=sum(Cases), N=sum(N))
count_q1   <- df_counts %>% filter(Comparison=="Q1") %>% summarise(Cases=sum(Cases), N=sum(N))
count_q2   <- df_counts %>% filter(Comparison=="Q2") %>% summarise(Cases=sum(Cases), N=sum(N))
count_q3   <- df_counts %>% filter(Comparison=="Q3") %>% summarise(Cases=sum(Cases), N=sum(N))
count_q4   <- df_counts %>% filter(Comparison=="Q4") %>% summarise(Cases=sum(Cases), N=sum(N))

# ===============================
# Build final results table
# ===============================
df_final <- tibble(
  Comparison = c("Continuous", "Q1", "Q2", "Q3", "Q4"),
  
  HR = c(
    df_meta_summary$HR[df_meta_summary$Comparison=="Continuous"],
    NA,
    df_meta_summary$HR[df_meta_summary$Comparison=="Q2 vs Q1"],
    df_meta_summary$HR[df_meta_summary$Comparison=="Q3 vs Q1"],
    df_meta_summary$HR[df_meta_summary$Comparison=="Q4 vs Q1"]
  ),
  
  HR_CI = c(
    df_meta_summary$HR_CI_formatted[df_meta_summary$Comparison=="Continuous"],
    NA,
    df_meta_summary$HR_CI_formatted[df_meta_summary$Comparison=="Q2 vs Q1"],
    df_meta_summary$HR_CI_formatted[df_meta_summary$Comparison=="Q3 vs Q1"],
    df_meta_summary$HR_CI_formatted[df_meta_summary$Comparison=="Q4 vs Q1"]
  ),
  
  CasesN = c(
    paste0(count_cont$Cases, "/", count_cont$N),
    paste0(count_q1$Cases, "/", count_q1$N),
    paste0(count_q2$Cases, "/", count_q2$N),
    paste0(count_q3$Cases, "/", count_q3$N),
    paste0(count_q4$Cases, "/", count_q4$N)
  )
)

print(df_final)

# ===============================
# Automatically write CSV (based on input file name)
# ===============================
output_csv <- paste0(file_base, "_results.csv")
write.csv(df_final, file = output_csv, row.names = FALSE)
cat("Saved result file:", output_csv, "\n")

# ===============================
# Save forest plots as PDFs (4 plots)
# ===============================
save_forest_pdf <- function(meta_obj, title, file_base, suffix="") {
  out_pdf <- paste0(file_base, "_", suffix, ".pdf")
  pdf(out_pdf, width = 10, height = 8)
  forest(meta_obj,
         main = title,
         print.pval.Q = TRUE,
         digits.pval.Q = 2,
         scientific.pval = TRUE)
  dev.off()
  cat("Saved PDF:", out_pdf, "\n")
}

save_forest_pdf(meta_continuous, "Continuous", file_base, "Continuous")
save_forest_pdf(meta_q2, "Q2 vs Q1", file_base, "Q2_vs_Q1")
save_forest_pdf(meta_q3, "Q3 vs Q1", file_base, "Q3_vs_Q1")
save_forest_pdf(meta_q4, "Q4 vs Q1", file_base, "Q4_vs_Q1")

cat("\nAll forest plots have been successfully saved as PDFs!\n")
