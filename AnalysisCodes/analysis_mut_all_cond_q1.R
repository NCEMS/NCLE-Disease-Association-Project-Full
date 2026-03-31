# Q1
# Are natively entangled regions of proteins positively associated
# with disease causing missense mutations?

## parsing af_control_mutation.R ##
cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("type: af \n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Libraries
library(dplyr)

# Mutation data
source("DataExtractionCode/m1_data.R")

############################################################################
# Contingency table + chi-square
# Compare Mis_Pathogenic vs in_region at residue level
############################################################################
contingency_table <- table(res_merged_dt$in_region, res_merged_dt$Mis_Pathogenic)

rownames(contingency_table) <- c("Not in Entangled Region", "In Entangled Region")
colnames(contingency_table) <- c("Missense Benign", "Missense Pathogenic")

print(contingency_table)

print(table(res_merged_dt$in_region, res_merged_dt$Mis_Pathogenic, res_merged_dt$Essential))

chi_sq_result <- chisq.test(contingency_table)
print(chi_sq_result)

# Manual odds ratio from 2x2 table:
# (InRegion & Pathogenic * NotInRegion & Benign) / (InRegion & Benign * NotInRegion & Pathogenic)
odds_ratio <- (as.numeric(contingency_table["In Entangled Region", "Missense Pathogenic"]) *
                 as.numeric(contingency_table["Not in Entangled Region", "Missense Benign"])) /
  (as.numeric(contingency_table["In Entangled Region", "Missense Benign"]) *
     as.numeric(contingency_table["Not in Entangled Region", "Missense Pathogenic"]))
print(odds_ratio)

############################################################################
# Logistic regression: Mis_Pathogenic ~ in_region + Length + Essential
############################################################################
res_merged_dt$Mis_Pathogenic <- factor(res_merged_dt$Mis_Pathogenic, levels = c("No", "Yes"))
res_merged_dt$Essential <- factor(res_merged_dt$Essential, levels = c("No", "NT", "Yes"))

logit_model_len_res <- glm(
  Mis_Pathogenic ~ in_region + scale(Length) + Essential,
  data = res_merged_dt,
  family = binomial
)

summary(logit_model_len_res)

coef_tab_len_res <- coef(summary(logit_model_len_res))
coef_tab_len_res

odds_ratios <- exp(coef(logit_model_len_res))

glm_len_res_summary <- summary(logit_model_len_res)

coef_table_len_res <- as.data.frame(coef(glm_len_res_summary))
coef_table_len_res$Variable <- rownames(coef_table_len_res)
coef_table_len_res$Odds_Ratios <- exp(coef_table_len_res$Estimate)

# 95% CI for estimates ORs
z_crit <- qnorm(0.975)  # 1.96 for 95% CI`

coef_table_len_res$CI_Lower_Est <- coef_table_len_res$Estimate - z_crit * coef_table_len_res$`Std. Error`
coef_table_len_res$CI_Upper_Est <- coef_table_len_res$Estimate + z_crit * coef_table_len_res$`Std. Error`

coef_table_len_res$CI_Lower_OR <- exp(coef_table_len_res$CI_Lower_Est)
coef_table_len_res$CI_Upper_OR <- exp(coef_table_len_res$CI_Upper_Est)

rownames(coef_table_len_res) <- NULL
coef_table_len_res <- coef_table_len_res[, c("Variable", "Estimate", "Odds_Ratios",
                                             "CI_Lower_OR", "CI_Upper_OR",
                                             "Std. Error", "z value", "Pr(>|z|)")]

# Build per-type output file name (not written because write.csv is commented out)
fileName <- paste("mutation_q1_type_af_pthr_", plddt_thresh, sep = "")
final_output_file <- paste0(
  "results/dataframes/", fileName,
  "_results_len_res_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"
)

# Save
combined_file <- paste0(
  "Results/Dataframes/mutation_q1_alltypes_pthr_",
  plddt_thresh, "_results_len_res_",
  format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"
)
write.csv(coef_table_len_res, combined_file, row.names = FALSE)
