library(readr)
library(readxl)
library(dplyr)

# File paths
protein_path <- "./input/olink_int_new.txt"
coef_path <- "./input/Brain.xlsx"

# 1. Read protein expression data
protein_df <- read.delim(protein_path, check.names = FALSE)

# 2. Read coefficient data (wide format, columns: organ, intercept, multiple proteins)
coef_df <- read_excel(coef_path)

# 3. Extract intercept (assuming intercept is in the first row of the 'intercept' column)
intercept_val <- coef_df$intercept[1]

# 4. Get protein coefficient column names (columns from the 3rd onward in coef_df are protein names)
coef_proteins <- colnames(coef_df)[-(1:2)]

# 5. Protein expression data column names
protein_cols <- colnames(protein_df)

# 6. Standardize to lowercase and trim whitespace to resolve case sensitivity issues
coef_proteins_lower <- tolower(trimws(coef_proteins))
protein_cols_lower <- tolower(trimws(protein_cols))

# 7. Find intersecting proteins (case-insensitive matching)
common_lower <- intersect(coef_proteins_lower, protein_cols_lower)

# 8. Retrieve the original column names from the protein expression data
idx_in_protein <- match(common_lower, protein_cols_lower)
protein_cols_matched <- protein_cols[idx_in_protein]

# 9. Retrieve the original column names from the coefficient data
idx_in_coef <- match(common_lower, coef_proteins_lower)
coef_proteins_matched <- coef_proteins[idx_in_coef]

# 10. Apply z-score normalization to the matched columns in the protein expression data
protein_z <- protein_df
protein_z[protein_cols_matched] <- scale(protein_df[protein_cols_matched])

# ✅ 10.5 Replace NA values after z-score normalization with 0
protein_z[protein_cols_matched][is.na(protein_z[protein_cols_matched])] <- 0

# 11. Extract beta coefficients for the corresponding proteins (order matches the protein columns)
beta_values <- as.numeric(coef_df[1, coef_proteins_matched])

# 12. Calculate predicted age = intercept + (protein z-score matrix × beta coefficient vector)
pred_ages <- as.matrix(protein_z[protein_cols_matched]) %*% beta_values + intercept_val

# 13. Save results, assuming the first column of protein expression data is the sample ID
sample_id_col <- colnames(protein_df)[1]

result_df <- data.frame(
  SampleID = protein_df[[sample_id_col]],
  PredictedAge = as.numeric(pred_ages)
)

# 14. Write to CSV file
write.csv(result_df, "./output/Predicted_Age.csv", row.names = FALSE)

cat("Predicted age calculation completed. Results saved to Predicted_Age.csv\n")