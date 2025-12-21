library(readr)
library(readxl)
library(dplyr)

# ===============================
# Input files
# ===============================
protein_path <- "./input/olink_int_new.txt"
coef_path <- "./input/Brain.xlsx"

# ===============================
# Data loading
# ===============================
protein_df <- read.delim(protein_path, check.names = FALSE)
coef_df <- read_excel(coef_path)

# ===============================
# Model parameters
# ===============================
intercept_val <- coef_df$intercept[1]
coef_proteins <- colnames(coef_df)[-(1:2)]

# ===============================
# Protein matching
# ===============================
protein_cols <- colnames(protein_df)

coef_proteins_lower <- tolower(trimws(coef_proteins))
protein_cols_lower <- tolower(trimws(protein_cols))

common_lower <- intersect(coef_proteins_lower, protein_cols_lower)

idx_in_protein <- match(common_lower, protein_cols_lower)
protein_cols_matched <- protein_cols[idx_in_protein]

idx_in_coef <- match(common_lower, coef_proteins_lower)
coef_proteins_matched <- coef_proteins[idx_in_coef]

# ===============================
# Normalization
# ===============================
protein_z <- protein_df
protein_z[protein_cols_matched] <- scale(protein_df[protein_cols_matched])
protein_z[protein_cols_matched][is.na(protein_z[protein_cols_matched])] <- 0

# ===============================
# Prediction
# ===============================
beta_values <- as.numeric(coef_df[1, coef_proteins_matched])
pred_ages <- as.matrix(protein_z[protein_cols_matched]) %*% beta_values + intercept_val

# ===============================
# Output
# ===============================
sample_id_col <- colnames(protein_df)[1]

result_df <- data.frame(
  SampleID = protein_df[[sample_id_col]],
  PredictedAge = as.numeric(pred_ages)
)

write.csv(result_df, "./output/Predicted_Age.csv", row.names = FALSE)

cat("Predicted age calculation completed. Results saved to Predicted_Age.csv\n")
