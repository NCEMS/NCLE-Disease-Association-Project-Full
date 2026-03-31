# Analysis 1.3
# logistic regression disease status ~ complexity-related metric values + Length + Essentiality

## parsing control_analysis1.3.R ##
cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input: Control1_3\n")
cat("type(s):", paste(types, collapse = "_and_"), "\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# output file names
fileName <- paste("a1_3_type_", paste(types, collapse = "_and_"), "_pthr_", plddt_thresh, sep = "")

# Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(rlang)
library(svglite)
library(grid)

result4slen_all <- list()   # Store results per type (af/crystal) before combining

# Loop over each dataset type requested (e.g., "af", "crystal")
for (cur_type in types) {

  type <- cur_type  # Current dataset type (kept as 'type' for downstream functions)

  # disease and entanglement association per disease MSH class data
  source("Functions/dataProcessing_df.R") # will load function that will load ent and disease data
  dfs <- dataProcessing_df(type, plddt_thresh) # function loaded to retrieve ent and disease data
  final_dfA <- dfs$final_dfA # disease data
  final_dfB <- dfs$final_dfB # entanglement data

  # remove values that were unable to be calculated
  final_dfB <- final_dfB %>% filter(!(Entanglement == "Yes" & (Travatos_G == 0 | is.na(Travatos_G))))

  final_dfA <- final_dfA %>%
    filter(uniprotids %in% final_dfB$gene)

  # If analyzing covalent, rebuild Entanglement status from covalent column
  if (identical(cov_vec, "covalent")) {

    keep_cols <- c("gene", "Length", "plddt", "Essential")

    rows_to_blank <- !(final_dfB$CCBond %in% "True") # Rows where CCBond is NOT true
    cols_to_blank <- setdiff(names(final_dfB), keep_cols) # All metric columns that should be blanked out

    final_dfB[rows_to_blank, cols_to_blank] <- NA # Blank metrics for rows failing the covalent bond criterion

    final_dfB <- final_dfB %>% select(-Entanglement) # Drop existing Entanglement column

    # Recreate Entanglement based on Gn/Gc thresholds (>= 0.6 in magnitude)
    final_dfB <- final_dfB %>%
      mutate(
        Entanglement = if_else(
          (abs(Gn) >= 0.6) | (abs(Gc) >= 0.6),
          "Yes", "No"
        )
      )

    # Aggregate to gene-level Entanglement: if any row is "Yes", gene is "Yes"; else "No"
    ent_map <- final_dfB %>%
      mutate(Entanglement = as.character(Entanglement)) %>%
      group_by(gene) %>%
      summarize(Entanglement = if_else(any(Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
                .groups = "drop")

    # Replace Entanglement in final_dfA by joining gene-level ent_map (default to "No" when missing)
    final_dfA <- final_dfA %>%
      select(-Entanglement) %>%                         # remove old column
      left_join(ent_map, by = c("uniprotids" = "gene")) %>%
      mutate(Entanglement = coalesce(Entanglement, "No"))

    final_dfB[["ENT.ID"]] <- final_dfB[["cc_count"]]

  }

  # Function where metric like Gmax and Gsum are created
  source("Functions/dataProcessing_df_analysis_3.R")
  percentiles <- rev(percentiles) # Reverse order so first iteration is the 50th threshold
  ######################################
  # Confounders per gene from final_dfB
  confounders <- final_dfB %>%
    group_by(gene) %>%
    summarise(
      Essential = first(Essential),
      .groups = "drop"
    ) %>%
    mutate(
      Essential = factor(Essential, levels = c("No", "NT", "Yes"))
    )

  # Create gene-level disease score: max(score) per UniProt ID (gene)
  gene_scores <- final_dfA %>%
    group_by(uniprotids) %>%
    summarise(
      score = max(score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(gene = uniprotids)

  result4slen_list <- list()

  # Define which metrics will be tested in the logistic regression
  all_metrics <- c("Gn", "Gc", "Gmax", "Gsum",
                   "N_term_thread",
                   "C_term_thread",
                   "ENT.ID",
                   "num_zipper_nc",
                   "num_loop_contacting_res",
                   "num_cross_nearest_neighbors",
                   "min_N_thread_depth_left",
                   "min_N_prot_depth_left",
                   "min_N_thread_slippage_left",
                   "min_C_thread_depth_right",
                   "min_C_prot_depth_right",
                   "min_C_thread_slippage_right",
                   "Travatos_G",
                   "Length")

  # Define which metrics will be tested in the logistic regression (no length)
  metrics_no_length <- setdiff(all_metrics, "Length")

  # Collapse final_df_updated to gene-level "max metric" data set
  # - For each metric column: use max per gene (or 0 if the whole group is NA)
  # - Length: max Length per gene
  final_dfB_max <- final_df_updated %>%
    group_by(gene) %>%
    summarise(
      across(
        all_of(metrics_no_length),
        \(x) if (all(is.na(x))) 0 else max(x, na.rm = TRUE) # If every value in that group/column is NA, return 0. Otherwise, return the maximum non-NA value.
      ),
      Length = max(Length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(all_of(metrics_no_length), ~replace_na(., 0)))

  # Add confounders and disease scores to the gene-level metrics table
  final_dfB_max <- final_dfB_max %>%
    left_join(confounders, by = "gene") %>%
    left_join(gene_scores, by = "gene")

  # Save the final gene-level dataframe used in the regression
  dir.create("Data/DataRanInAnalysis", recursive = TRUE, showWarnings = FALSE)

  write.csv(
    final_dfB_max,
    file = paste0(
      "Data/DataRanInAnalysis/",
      fileName,
      "_final_df_used_in_regression",
      format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
      ".csv"
    ),
    row.names = FALSE
  )

  # Create a named list of data frames, one per entanglement metric
  df_G_ent_slen <- setNames(
    lapply(metrics_no_length, function(m) { # Loop over each metric
      final_dfB_max %>%
        select(gene, !!sym(m), Length, Essential, score) # Select the gene ID, the current metric column, protein length, essentiality, and disease score
    }),
    metrics_no_length # Name each list element using the corresponding metric name
  )

  # - If remove_mode is TRUE, then for i>1 you set disease_association to NA for genes that were positive at i=1
  #   but are currently labeled negative, effectively removing those negatives from later models.
  remove_mode <- identical(tolower(get0("option", ifnotfound = "keep")), "remove")
  pos_uniprots_i1 <- NULL
  # --------------------------------------

  i=1
  for (i in 1:length(percentiles)) { # Loop over each percentile threshold (e.g., 95%, 75%, 50%)

    # Logistic Regression for each metric + Length
    logit_length <- function(df, metric, score_thresh) {
      df$disease_association <- ifelse(df$score >= score_thresh, 1, 0)

      # Apply "Remove" rule for i > 1
      if (remove_mode && i != 1L) {
        if ("gene" %in% names(df)) { # Check that df has a gene column
          drop_idx <- df$gene %in% pos_uniprots_i1 & df$disease_association == 0L # Flags rows where the gene was positive in iteration 1 but is currently labeled disease negative
          if (any(drop_idx)) df$disease_association[drop_idx] <- NA_integer_ # Set disease_association to NA for flagged rows
        }
      }

      mask0 <- (df$disease_association == 0L)

      # Counts
      cat("Unique genes not disease-associated at", names(percentiles)[i], ":",
          dplyr::n_distinct(df$gene[mask0]),
          "\n")
      print(table(df$disease_association))

      # Build model formula: outcome ~ (metric) + scaled Length + Essentiality
      formula <- as.formula(paste("disease_association ~", metric,"+ scale(Length) + Essential"))
      model <- glm(formula, family = "binomial", data = df) # Fit logistic regression
      conf_int <- suppressMessages(confint(model))  # Compute confidence intervals

      # Extract coefficient for the metric
      out <- c(Estimate = as.numeric(coef(model)[2]),
               OR = as.numeric(exp(as.numeric(coef(model)[2]))),
               Std_Error = summary(model)$coefficients[2, "Std. Error"],
               p_value = summary(model)$coefficients[2, "Pr(>|z|)"],
               CI_Lower = exp(conf_int[2, 1]),
               CI_Upper = exp(conf_int[2, 2]))
      return(out)
    }


    # At i == 1, cache which genes are disease-associated under the first (50th) threshold
    if (i == 1L) {
      thr <- as.numeric(percentiles[ i ][[1]])
      fdu <- final_df_updated
      pos_uniprots_i1 <- unique(fdu$gene[fdu$score >= thr])
      pos_uniprots_i1_not <- unique(fdu$gene[fdu$score < thr])
      message(length(pos_uniprots_i1), " disease associated genes at ", names(percentiles)[i], " (thr=", thr, ")")
      message(length(pos_uniprots_i1_not), " not disease associated genes at ", names(percentiles)[i], " (thr=", thr, ")")
    }
    # -------------------------------------------------------------------------------
    # For each metric-specific data frame in df_G_ent_slen,
    # run a logistic regression that includes scaled length
    out_logit_slength <- lapply(names(df_G_ent_slen), function(metric) {
      logit_length(df_G_ent_slen[[metric]], metric, score_thresh = percentiles[i])
    })

    # Combine the list of regression outputs into a single data frame
    out_logit_slength <- do.call(rbind,out_logit_slength) %>% data.frame()


    # Add a column identifying which metric each row of results corresponds to
    out_logit_slength$metric <- names(df_G_ent_slen)


    # Add a column indicating which percentile threshold was used
    out_logit_slength$percentile <- names(percentiles[i])

    # Store the combined results for this percentile iteration in the results list
    result4slen_list[[i]] <- out_logit_slength
  }

  # Stack results across percentiles into one data frame for this type
  result4slen <- do.call(rbind, result4slen_list)


  # Keep core columns
  result4slen <- result4slen[,c("metric", "percentile" , "p_value",
                                "Estimate", "OR", "Std_Error", "CI_Lower","CI_Upper")]


  # After creating result4slen at the end of your analysis add the type
  result4slen$type <- type
  result4slen_all[[type]] <- result4slen

}

# Combine both types
result4slen_combined <- bind_rows(result4slen_all)

# BH across types + metrics ONLY (i.e., within each percentile)
result4slen_combined <- result4slen_combined %>%
  group_by(percentile) %>%
  mutate(p_value_corrected_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Save to file
write.csv(result4slen_combined, file = paste0("Results/Dataframes/", fileName, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

# Create a "50%" filtered table for the OR plot
result4slen_filtered <- result4slen_combined %>%
  filter(percentile == "50%") %>%
  select(
    metric,
    OR,
    CI_Lower,
    CI_Upper,
    p_value_corrected_BH
  ) %>%
  mutate(
    Significant = p_value_corrected_BH <= 0.05
  )%>%
  arrange(desc(OR))

# Save the 50% filtered table
write.csv(result4slen_filtered, file = paste0("Results/Dataframes/", fileName,
                                     format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), "_plot_values", ".csv"), row.names = FALSE)

# Map raw metric names to plot-ready expressions (parsed into italics/subscripts)
label_map <- c(
  "min_C_prot_depth_right"      = "italic(d)[P](C)",
  "min_N_prot_depth_left"       = "italic(d)[P](N)",
  "min_N_thread_depth_left"     = "italic(d)[T](N)",
  "min_N_thread_slippage_left"  = "italic(d)[S](N)",
  "min_C_thread_depth_right"    = "italic(d)[T](C)",
  "min_C_thread_slippage_right" = "italic(d)[S](C)",
  "Travatos_G" = "\"max\"*\"|\"*italic(G)[Intra]*\"|\"",
  "Gsum"                        = "italic(G)[Sum]",
  "Gn"                          = "italic(G)[n]",
  "Gc"                          = "italic(G)[c]",
  "Gmax"                        = "italic(G)[Max]",
  "ENT.ID"                      = "italic(E)[NCLE]",
  "N_term_thread"               = "italic(N)[threads]",
  "C_term_thread"               = "italic(C)[threads]",
  "num_zipper_nc"               = "italic(N)[zipper]",
  "num_loop_contacting_res"     = "italic(N)[loop-cont]",
  "num_cross_nearest_neighbors" = "italic(N)[cross-cont]"
  )


# OR plot
x_min <- 0.8
x_max <- 1.85

p <- ggplot(
  result4slen_filtered,
  aes(
    y = reorder(metric, OR),
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
    labels = function(x) {
      lab <- label_map[x]
      lab[is.na(lab)] <- x[is.na(lab)]
      parse(text = lab)
    },
    expand = expansion(add = 0.35)
  ) +
  labs(
    y = "Entanglement Metric",
    x = "Odds Ratio"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),

    axis.text.x  = element_text(size = 18, color = "black"),
    axis.text.y  = element_text(size = 18, color = "black", lineheight = 1),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),

    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(8, "pt"),

    legend.position = "none",

    plot.margin = ggplot2::margin(t = 12, r = 20, b = 12, l = 35, unit = "pt")
  )

ggsave(
  filename = paste0(
    "Results/Plots/", fileName,
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
    "Results/Plots/", fileName,
    format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
    "_plot.svg"
  ),
  plot   = p + theme(text = element_text(family = "sans")),
  device = pdf,
  width  = 11,
  height = 10,
  units  = "in"
)
