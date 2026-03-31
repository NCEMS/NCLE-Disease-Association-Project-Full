# Analysis 1.2

## parsing Control1_2.R ##
cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input: Control1_2\n")
cat("type(s):", paste(types, collapse = "_and_"), "\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# output file names
(fileName = paste("a1_2_type_", paste(types, collapse = "_and_"), "_pthr_", plddt_thresh, "_", option, sep=""))

# libraries
library(dplyr)
library(broom)
library(stringr)
library(ggplot2)
library(svglite)

# where results will be saved
res_list  <- list()
data_list <- list()
fileNames <- character(0)

# function to run analysis
run_one <- function(type, cov_type, plddt_thresh) {
  # disease and entanglement association per disease MSH class data
  source("Functions/dataProcessing_df.R") # will load function that will load ent and disease data
  dfs <- dataProcessing_df(type, plddt_thresh) # function loaded to retrieve ent and disease data
  final_dfA <- dfs$final_dfA # disease data
  final_dfB <- dfs$final_dfB # entanglement data

  # Swap in alternative “entanglement-like” predictors depending on cov_type
  if (identical(cov_type, "covalent")) {
    final_dfA$Entanglement <- final_dfA$Cov_Entanglement
    final_dfB$Entanglement <- final_dfB$Cov_Entanglement
  } else if (identical(cov_type, "knot")) {
    final_dfA$Entanglement <- final_dfA$Knot
    final_dfB$Entanglement <- final_dfB$Knot
  }
  # Otherwise (e.g., "noncovalent"), keep Entanglement as-is

  #### Disgenet Score Percentiles ####
  # Used to calculate percentiles from "raw" disgenet results
  sub_finalA <- final_dfA[,c("Association_ID","score")]
  sub_finalA <- unique(sub_finalA) # Drop duplicate Association_ID-score pairs
  sub_finalA <- sub_finalA[sub_finalA$Association_ID>0,] # Keep DisGeNET rows
  # Calculate score percentiles
  percentiles <- quantile(sub_finalA$score, c(0.95, 0.75, 0.50))
  rm(sub_finalA) # Remove temporary object to save memory

  ### Prepare a dataframe for contingency tables (one per disease class) ####
  # Subset necessary columns
  UniScoDis <- final_dfA[, c("uniprotids", "score","diseaseClasses_MSH", "Entanglement")]

  # Remove duplicates
  UniScoDis <- distinct(UniScoDis)

  # pull out lengths and essentiality per uniprot ID
  length_df <- final_dfB %>%
    select(gene, Length, Essential) %>%
    distinct()

  # join onto final_df
  UniScoDis <- UniScoDis %>%
    left_join(length_df, by = c("uniprotids" = "gene")) # Add Length/Essential per protein

  # List of all disease classes
  all_disease_classes <- c(
    "Infections (C01)",
    "Neoplasms (C04)",
    "Musculoskeletal Diseases (C05)",
    "Digestive System Diseases (C06)",
    "Stomatognathic Diseases (C07)",
    "Respiratory Tract Diseases (C08)",
    "Otorhinolaryngologic Diseases (C09)",
    "Nervous System Diseases (C10)",
    "Eye Diseases (C11)",
    "Urogenital Diseases (C12)",
    "Cardiovascular Diseases (C14)",
    "Hemic and Lymphatic Diseases (C15)",
    "Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)",
    "Skin and Connective Tissue Diseases (C17)",
    "Nutritional and Metabolic Diseases (C18)",
    "Endocrine System Diseases (C19)",
    "Immune System Diseases (C20)",
    "Disorders of Environmental Origin (C21)",
    "Pathological Conditions, Signs and Symptoms (C23)",
    "Chemically-Induced Disorders (C25)",
    "Behavior and Behavior Mechanisms (F01)",
    "Psychological Phenomena (F02)",
    "Mental Disorders (F03)"
  )

  # Function to create new columns based on disease class codes
  create_disease_columns <- function(data, disease_classes) {
    # For each disease class, create a column:
    #   - value = score if that class code appears in diseaseClasses_MSH
    #   - value = 0 otherwise
    for (disease_class in disease_classes) {
      # extract code part of disease_class e.g., (C01)
      code <- str_extract(disease_class, "\\([CF]\\d+\\)")
      # create a new column
      # where each element = score if the code is in diseaseClasses_MSH, o/w 0
      data[[disease_class]] <- ifelse(grepl(code, data$diseaseClasses_MSH), data$score, 0)
    }
    return(data)
  }

  # Call the function to add columns
  UniScoDis <- create_disease_columns(UniScoDis, all_disease_classes)


  # Remove the score and diseaseClasses_MSH columns
  UniScoDis <- subset(UniScoDis, select = -c(score, diseaseClasses_MSH))

  # Group by uniprotids and calculate maximum scores for each disease class
  max_scores <- UniScoDis %>%
    group_by(uniprotids) %>%
    mutate(
      `Digestive System Diseases (C06)` = max(`Digestive System Diseases (C06)`),
      `Neoplasms (C04)` = max(`Neoplasms (C04)`),
      `Pathological Conditions, Signs and Symptoms (C23)` = max(`Pathological Conditions, Signs and Symptoms (C23)`),
      `Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)` = max(`Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)`),
      `Endocrine System Diseases (C19)` = max(`Endocrine System Diseases (C19)`),
      `Urogenital Diseases (C12)` = max(`Urogenital Diseases (C12)`),
      `Respiratory Tract Diseases (C08)` = max(`Respiratory Tract Diseases (C08)`),
      `Nervous System Diseases (C10)` = max(`Nervous System Diseases (C10)`),
      `Nutritional and Metabolic Diseases (C18)` = max(`Nutritional and Metabolic Diseases (C18)`),
      `Stomatognathic Diseases (C07)` = max(`Stomatognathic Diseases (C07)`),
      `Eye Diseases (C11)` = max(`Eye Diseases (C11)`),
      `Musculoskeletal Diseases (C05)` = max(`Musculoskeletal Diseases (C05)`),
      `Cardiovascular Diseases (C14)` = max(`Cardiovascular Diseases (C14)`),
      `Infections (C01)` = max(`Infections (C01)`),
      `Immune System Diseases (C20)` = max(`Immune System Diseases (C20)`),
      `Skin and Connective Tissue Diseases (C17)` = max(`Skin and Connective Tissue Diseases (C17)`),
      `Otorhinolaryngologic Diseases (C09)` = max(`Otorhinolaryngologic Diseases (C09)`),
      `Hemic and Lymphatic Diseases (C15)` = max(`Hemic and Lymphatic Diseases (C15)`),
      `Mental Disorders (F03)` = max(`Mental Disorders (F03)`),
      `Behavior and Behavior Mechanisms (F01)` = max(`Behavior and Behavior Mechanisms (F01)`),
      `Chemically-Induced Disorders (C25)` = max(`Chemically-Induced Disorders (C25)`),
      `Psychological Phenomena (F02)` = max(`Psychological Phenomena (F02)`),
      `Disorders of Environmental Origin (C21)` = max(`Disorders of Environmental Origin (C21)`)
    ) %>%
    distinct() # Keep one row per protein after max-ing

  #### Prepare data frames for Contingency tables ####

  # if score>=thresh, disease == "Yes"
  convert_to_yes_no <- function(value, thresh) {
    ifelse(value >= thresh, "Yes", "No")
  }

  # Apply the function to all disease class columns columns
  UniScoDis95 <- max_scores %>%
    mutate(across(all_of(all_disease_classes), ~convert_to_yes_no(., percentiles[1])))
  UniScoDis75 <- max_scores %>%
    mutate(across(all_of(all_disease_classes), ~convert_to_yes_no(., percentiles[2])))
  UniScoDis50 <- max_scores %>%
    mutate(across(all_of(all_disease_classes), ~convert_to_yes_no(., percentiles[3])))

  # vector of class-column names
  class_cols <- all_disease_classes

  # function to run one batch (pct_df + label)
  run_logistics <- function(df, pct_label) {
    # Fits one logistic regression per disease class column in df for a given percentile definition.
    # Returns a tibble with OR, CI, p-values, convergence flag, and raw cell counts.

    res_list <- vector("list", length(class_cols))

    for (i in seq_along(class_cols)) {
      cls <- class_cols[i] # Current disease class column name
      df2 <- df %>%
        mutate(
          Disease      = factor(.data[[cls]], levels = c("No","Yes")), # Outcome for this class (Yes/No)
          Entanglement = factor(Entanglement, levels = c("No","Yes")),
          Essential    = factor(Essential, levels = c("No", "NT", "Yes")),
          Length       = as.numeric(Length)
        )

      # count the 4 cells directly
      ent_yes_d_yes <- sum(df2$Entanglement == "Yes" & df2$Disease == "Yes", na.rm = TRUE)
      ent_no_d_yes  <- sum(df2$Entanglement == "No"  & df2$Disease == "Yes", na.rm = TRUE)
      ent_yes_d_no  <- sum(df2$Entanglement == "Yes" & df2$Disease == "No", na.rm = TRUE)
      ent_no_d_no   <- sum(df2$Entanglement == "No"  & df2$Disease == "No", na.rm = TRUE)

      ## ---- skip if Entanglement or Disease are all one value ----
      # If either variable has no variation, logistic regression is not identifiable.
      ent_levels <- unique(as.character(df2$Entanglement))
      dis_levels <- unique(as.character(df2$Disease))
      if (length(ent_levels) < 2 || length(dis_levels) < 2) {
        warn_msg <- sprintf(
          "Skipped: no variation in %s (Entanglement levels=%d, Disease levels=%d)",
          if (length(ent_levels) < 2 && length(dis_levels) < 2) "both Entanglement & Disease"
          else if (length(ent_levels) < 2) "Entanglement" else "Disease",
          length(ent_levels), length(dis_levels)
        )
        warning(sprintf("Skipping '%s' at %s — %s", cls, pct_label, warn_msg))

        # Return an NA row but still include the cell counts + warning text
        res_list[[i]] <- tibble(
          Percentile                   = pct_label,
          Disease_Class                = cls,
          Estimate                     = NA_real_,
          OR                           = NA_real_,
          CI_lower                     = NA_real_,
          CI_upper                     = NA_real_,
          p_value                      = NA_real_,
          DidNotConverge               = NA_real_,
          Entangled_and_Disease        = ent_yes_d_yes,
          NotEntangled_and_Disease     = ent_no_d_yes,
          Entangled_and_NoDisease      = ent_yes_d_no,
          NotEntangled_and_NoDisease   = ent_no_d_no,
          Warning = warn_msg
        )
        next  # Move to next disease class
      }

      # fit the model
      fit <- glm(
        Disease ~ Entanglement + scale(Length) + Essential,
        data    = df2,
        family  = binomial()
      )

      # Wald CIs
      ci_prof <- suppressWarnings(confint.default(fit)) # Compute coefficient CI (Wald)
      ci_term <- ci_prof["EntanglementYes", ] # Extract CI row for EntanglementYes

      # Extract and format model results for the entanglement term only
      td <- tidy(fit) %>%
        filter(term == "EntanglementYes") %>%
        transmute(
          Percentile                   = pct_label,
          Disease_Class                = cls,
          Estimate                     = estimate,
          OR                           = exp(estimate),
          CI_lower                     = exp(ci_term[1]),
          CI_upper                     = exp(ci_term[2]),
          p_value                      = p.value,
          DidNotConverge               = !fit$converged,  # Convergence flag
          Entangled_and_Disease        = ent_yes_d_yes,
          NotEntangled_and_Disease     = ent_no_d_yes,
          Entangled_and_NoDisease      = ent_yes_d_no,
          NotEntangled_and_NoDisease   = ent_no_d_no
        )

      res_list[[i]] <- td  # Store this class’s result
    }

    bind_rows(res_list) # Combine all classes into one tibble
  }

  ## === Run models; include the "Remove" behavior for 75/95 inside the if ======
  opt <- tolower(get0("option", ifnotfound = "keep"))

  if (opt == "remove") {
    # Copy frames to modify
    UniScoDis75_mod <- UniScoDis75
    UniScoDis95_mod <- UniScoDis95

    # For each disease class: if 50th=="Yes" but 75th/95th=="No", set to NA
    for (cls in class_cols) {
      UniScoDis75_mod[[cls]] <- ifelse(
        UniScoDis50[[cls]] == "Yes" & UniScoDis75_mod[[cls]] == "No",
        NA_character_,
        UniScoDis75_mod[[cls]]
      )
      UniScoDis95_mod[[cls]] <- ifelse(
        UniScoDis50[[cls]] == "Yes" & UniScoDis95_mod[[cls]] == "No",
        NA_character_,
        UniScoDis95_mod[[cls]]
      )
    }

    # Run models USING the modified 75/95 and the original 50
    res95_spec <- run_logistics(UniScoDis95_mod, "95th")
    res75_spec <- run_logistics(UniScoDis75_mod, "75th")
    res50_spec <- run_logistics(UniScoDis50,     "50th")

    # combine
    all_df <- bind_rows(
      UniScoDis95_mod %>% mutate(Percentile = "95th"),
      UniScoDis75_mod %>% mutate(Percentile = "75th"),
      UniScoDis50 %>% mutate(Percentile = "50th")
    )

  } else {
    # No removals: run models on the original tables
    res95_spec <- run_logistics(UniScoDis95, "95th")
    res75_spec <- run_logistics(UniScoDis75, "75th")
    res50_spec <- run_logistics(UniScoDis50, "50th")

    # combine
    all_df <- bind_rows(
      UniScoDis95 %>% mutate(Percentile = "95th"),
      UniScoDis75 %>% mutate(Percentile = "75th"),
      UniScoDis50 %>% mutate(Percentile = "50th")
    )

  }

  # combine all percentile-specific model results
  all_results_spec <- bind_rows(res95_spec, res75_spec, res50_spec)

  # add type & cov_type as leading columns
  all_results_spec <- all_results_spec %>%
    mutate(type = type, cov_type = cov_type) %>%
    relocate(type, cov_type, .before = everything())

  all_df <- all_df %>%
    mutate(type = type, cov_type = cov_type) %>%
    relocate(type, cov_type, .before = everything())

   # return outputs for this one run
  list(
    results  = all_results_spec,
    data     = all_df,
    fileName = fileName
  )
}

# ---- run all combos ----
for (tt in types) { # Loop over dataset types (e.g., af, crystal)
  for (cc in cov_vec) { # Loop over entanglement variants (e.g., noncovalent, covalent, knot)
    out <- run_one(tt, cc, plddt_thresh = plddt_thresh) # Run analysis for this combination
    res_list[[length(res_list) + 1L]] <- out$results  # Append model results
    data_list[[length(data_list) + 1L]] <- out$data  # Append run data tables
    fileNames <- unique(c(fileNames, out$fileName))

    # delete everything except the required objects
    to_keep <- c("run_one", "res_list", "data_list", "types", "cov_vec", "plddt_thresh", "option", "fileName", "fileNames", "tt", "cc")
    rm(list = setdiff(ls(), to_keep))
    invisible(gc()) # frees memory that is no longer being used
  }
}


# bind all into the final data frame and assign to `fileName`
final_results <- bind_rows(res_list)
all_data <- bind_rows(data_list)

# BH correction across types (af, crystal) for the same cov_type and threshold (Percentile)
final_results <- final_results %>%
  group_by(cov_type, Percentile) %>%  # Correct within each cov_type x percentile family
  mutate(
    n_tested_BH = sum(!is.na(p_value)), # Number of non-NA tests in this group
    p_value_BH  = { # Compute BH-adjusted p-values (keeping NAs)
      pv <- p_value
      idx <- which(!is.na(pv)) # Indices of valid p-values
      adj <- rep(NA_real_, length(pv))  # Initialize adjusted vector
      if (length(idx) > 0) adj[idx] <- p.adjust(pv[idx], method = "BH") # Apply BH only to non-NA
      adj # Return adjusted p-values aligned to rows
    }
  ) %>%
  ungroup() %>%
  # keep identifiers up front
  relocate(type, cov_type, Percentile, Disease_Class, .before = everything())

assign(fileName, final_results)  # Assign results into an object named by fileName

rm(list = setdiff(ls(), c("final_results", "all_data", "fileName")))

# Save:
write.csv(final_results, file = paste0("Results/Dataframes/", fileName, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
write.csv(all_data, file = paste0("Data/DataRanInAnalysis/", fileName, format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

# Plot
# Build a plotting subset: AF + noncovalent + 50th percentile
final_results_filtered <- final_results %>%
  filter(
    type == "af",
    cov_type == "noncovalent",
    Percentile == "50th"
  ) %>%
  filter(
    is.finite(CI_lower),
    is.finite(CI_upper)
  ) %>%
  select(
    Disease_Class,
    OR,
    CI_lower,
    CI_upper,
    p_value_BH
  ) %>%
  mutate(
    Significant = p_value_BH <= 0.05
  )

final_results_filtered$Disease_Class[
  final_results_filtered$Disease_Class ==
    "Congenital, Hereditary, and Neonatal Diseases and Abnormalities (C16)"
] <- "Congenital and Hereditary Disorders (C16)"


# Plot of ORs by disease class
x_min <- 0.2
x_max <- 2.2

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
    aes(xmin = CI_lower, xmax = CI_upper),
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
    "Results/Plots/a1_2_af_noncov_50th",
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
    "Results/Plots/a1_2_af_noncov_50th",
    format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"),
    "_plot.svg"
  ),
  plot   = p + theme(text = element_text(family = "sans")),
  device = pdf,
  width  = 11,
  height = 10,
  units  = "in"
)

write.csv(final_results_filtered, file = paste0("Results/Dataframes/", fileName,"_plot_values", format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)


