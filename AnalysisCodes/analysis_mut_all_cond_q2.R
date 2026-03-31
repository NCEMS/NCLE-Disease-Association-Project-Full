# Q2
# Are entangled proteins positively associated with
# disease causing missense mutations?

cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input:\n")
cat("types: af\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Mutation data
source("DataExtractionCode/m2_data.R")

# If disease associated only
if (disease_only == "yes") {
  protein_summary <- protein_summary %>%
    filter(`50th_percentile` == "Yes")
}
###############################################################################
# MODEL
###############################################################################
glm_model_len <- glm(
  cbind(Path_Mutations_m, Total_Mutations_n - Path_Mutations_m) ~
    Entanglement+ scale(Length) + Essential,
  data = protein_summary,
  family = binomial
)

glm_len_summary <- summary(glm_model_len)

# Extract coefficient table and add odds ratios
coef_table_len <- as.data.frame(coef(glm_len_summary))
coef_table_len$Variable <- rownames(coef_table_len)
coef_table_len$Odds_Ratios <- exp(coef_table_len$Estimate)
coef_table_len$CI_Lower <- exp(coef_table_len$Estimate - 1.96 * coef_table_len$`Std. Error`)
coef_table_len$CI_Upper <- exp(coef_table_len$Estimate + 1.96 * coef_table_len$`Std. Error`)
rownames(coef_table_len) <- NULL
coef_table_len <- coef_table_len[, c(
  "Variable","Estimate","Std. Error","Odds_Ratios",
  "CI_Lower", "CI_Upper", "z value","Pr(>|z|)"
)]

# File naming
fileName <- paste("mutation_q2_type_af_pthr_",plddt_thresh,sep="")
final_output_file <- paste0(
  "Results/Dataframes/", fileName,
  "_results_len_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"
)

# Save
write.csv(coef_table_len, final_output_file, row.names = FALSE)
