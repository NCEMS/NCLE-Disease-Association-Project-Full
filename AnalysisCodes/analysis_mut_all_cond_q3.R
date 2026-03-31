# Q3
# Do certain disease classes show greater association between
# entanglements and disease-causing mutations?

cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input: Controlm1_3\n")
cat("types: af\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Data Checks
# Libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(svglite)

# Mutation data
source("DataExtractionCode/m3_data.R")

# Get list of disease classes (excluding uniprotids and Entanglement columns)
disease_classes <- setdiff(names(UniScoDis50), c("uniprotids", "Entanglement"))

# Initialize list to store model results
results <- list()

for (class in disease_classes) {
  cat("\n--- Analyzing Disease Class:", class, "---\n")

  # Subset UniScoDis50 to only proteins associated with this disease class
  disease_subset_ids <- UniScoDis50 %>%
    filter(.data[[class]] == "Yes") %>%
    pull(uniprotids)

  # Subset final_df to those proteins
  df_subset <- final_df %>%
    filter(uniprotids %in% disease_subset_ids)

  # Count number of unique proteins in this subset
  protein_count <- n_distinct(df_subset$uniprotids)

  # Summarize mutation and pathogenicity info
  protein_summary <- df_subset %>%
    group_by(uniprotids) %>%
    summarise(
      Entanglement = first(Entanglement),
      Length = first(Length),
      Essential = first(Essential),
      Total_Mutations_n = sum(Mutation == "Yes"),
      Path_Mutations_m = sum(Pathogenic == "Yes", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(p_hat_m_div_n = if_else(Total_Mutations_n > 0,
                                   Path_Mutations_m / Total_Mutations_n,
                                   NA_real_))

  # Save the final dataframe used in the regression for this disease class
  dir.create("Data/DataRanInAnalysis", recursive = TRUE, showWarnings = FALSE)

  write.csv(
    protein_summary,
    file = paste0(
      "Data/DataRanInAnalysis/",
      "mutation_q3_type_af_pthr_", plddt_thresh,
      "_", make.names(class),
      "_protein_summary_used_in_glm",
      format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
      ".csv"
    ),
    row.names = FALSE
  )

  # Make sure Entanglement and Essential are a factor
  protein_summary$Entanglement <- factor(protein_summary$Entanglement, levels = c("No", "Yes"))
  protein_summary$Essential <- factor(protein_summary$Essential, levels = c("No", "NT", "Yes"))

  # Skip if only one level of Entanglement of Essential
  if (nlevels(droplevels(protein_summary$Entanglement)) < 2) {
    cat("Skipped: Only one level of Entanglement\n")
    next
  }

  if (nlevels(droplevels(protein_summary$Essential)) < 2) {
    cat("Skipped:", class, "- only one level of Essential\n")
    next
  }

  # prepare warning flags
  warn_len  <- FALSE

  # length + essential model, catch warnings
  glm_model_len <- withCallingHandlers(
    glm(
      cbind(Path_Mutations_m, Total_Mutations_n - Path_Mutations_m) ~
        Entanglement + scale(Length) + Essential,
      data = protein_summary,
      family = binomial
    ),
    warning = function(w) {
      warn_len <<- TRUE
      # print each warning as it happens:
      cat(conditionMessage(w),
          "in length + essential model for class", class, "\n")
      invokeRestart("muffleWarning")
    }
  )

  # Extract coefficient info
  coef_len  <- coef(summary(glm_model_len))["EntanglementYes", ]

  # Store results
  results[[class]] <- list(
    class = class,
    protein_count = protein_count,
    p_value_len = coef_len["Pr(>|z|)"],
    model_with_length = summary(glm_model_len),
    odds_ratios_len = exp(coef(glm_model_len)),
    conf_int_len = exp(confint(glm_model_len)),
    warn_len         = warn_len
  )

}

# Create summary of odds ratios including the length+SASA model
summary_table <- do.call(rbind, lapply(results, function(res) {
  # helper to safely extract a named value (or NA)
  safe_extract <- function(x, name) if (!is.null(x) && name %in% names(x)) x[name] else NA

  data.frame(
    Disease_Class                   = res$class,
    Protein_Count                   = res$protein_count,

    ## Length + Essentiality model
    OR_Entangled           = round(safe_extract(res$odds_ratios_len,     "EntanglementYes"), 3),
    CI_Lower               = round(safe_extract(res$conf_int_len[,1],    "EntanglementYes"), 3),
    CI_Upper               = round(safe_extract(res$conf_int_len[,2],    "EntanglementYes"), 3),
    OR_Length              = round(safe_extract(res$odds_ratios_len,     "scale(Length)"),   3),
    OR_ESS                 = round(safe_extract(res$odds_ratios_len,     "EssentialYes"),3),
    P_Value                = if (!is.na(res$p_value_len)) signif(res$p_value_len, 3) else NA,
    Coef_Estimate_Entangled = round(safe_extract(res$model_with_length$coefficients[, "Estimate"], "EntanglementYes"), 4),
    Std_Error_Entangled     = round(safe_extract(res$model_with_length$coefficients[, "Std. Error"], "EntanglementYes"), 4)
  )
}))


print(summary_table)

# Create summary
len_df <- summary_table %>%
  select(Disease_Class, Protein_Count, Coef_Estimate_Entangled, Std_Error_Entangled,
         OR_Entangled, CI_Lower, CI_Upper, P_Value) %>%
  rename(
    OR      = OR_Entangled
  )



# At least 20 proteins per class needed
len_df <- len_df[len_df$Protein_Count >= 20, ]

# Create a BH-adjusted p-value across ALL rows (all types & disease classes)
len_df$P_Value_BH_all <- NA_real_

non_na_idx <- which(!is.na(len_df$P_Value))
len_df$P_Value_BH_all[non_na_idx] <- p.adjust(
  len_df$P_Value[non_na_idx],
  method = "BH"
)

rownames(len_df) <- NULL


final_results_filtered <- len_df %>%
  select(
    Disease_Class,
    Coef_Estimate_Entangled,
    OR,
    CI_Lower,
    CI_Upper,
    P_Value_BH_all,
    Protein_Count
  ) %>%
  mutate(
    Significant = P_Value_BH_all <= 0.05
  )

# File naming
fileName <- paste("mutation_q3_type_af_pthr_",plddt_thresh,sep="")
final_output_file <- paste0(
  "Results/Dataframes/", fileName,
  "_results_len_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv"
)

# Save
write.csv(final_results_filtered, final_output_file, row.names = FALSE)

# Shorten name for plot
final_results_filtered$Disease_Class[
  final_results_filtered$Disease_Class ==
    "Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)"
] <- "Congenital and Hereditary Disorders (C16)"

# Plot of ORs by disease class
x_min <- 0.6
x_max <- 2.4

p <- ggplot(
  final_results_filtered,
  aes(
    y = reorder(Disease_Class, OR),
    x = OR,
    shape = Significant,
    color = Significant
  )
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 1.2) +
  geom_point(size = 4.2, stroke = 1.2) +
  geom_errorbarh(
    aes(xmin = CI_Lower, xmax = CI_Upper),
    height = 0.22,
    linewidth = 1.2
  ) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "red")) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = seq(x_min, x_max, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 28),
    expand = expansion(add = 0.35)
  ) +
  labs(
    y = "Disease Class",
    x = "Odds Ratio"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),

    axis.text.x  = element_text(size = 14, color = "black"),
    axis.text.y  = element_text(size = 14, color = "black", lineheight = 1),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),

    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(8, "pt"),

    legend.position = "none",

    plot.margin = ggplot2::margin(t = 12, r = 20, b = 12, l = 35, unit = "pt")
  )

ggsave(
  filename = paste0(
    "Results/Plots/mut_q3_af_noncov_50th",
    format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
    "_plot.svg"
  ),
  plot   = p,
  device = svglite::svglite,
  width  = 11,
  height = 10,
  units  = "in"
)

ggsave(
  filename = paste0(
    "Results/Plots/mut_q3_af_noncov_50th",
    format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
    "_plot.svg"
  ),
  plot   = p + theme(text = element_text(family = "sans")),
  device = pdf,
  width  = 11,
  height = 10,
  units  = "in"
)
