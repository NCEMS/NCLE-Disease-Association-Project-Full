# Score
type = "af"
plddt_thresh = 70

cat("type:", type, "\n")
cat("plddt_thresh:", plddt_thresh, "\n")

library(dplyr)
library(tidyr)
library(readr)

# process dfA, dfB using type, dfA, dfB
source("Functions/dataProcessing_df.R")
dfs <- dataProcessing_df(type, plddt_thresh)
final_dfA <- dfs$final_dfA
final_dfB <- dfs$final_dfB

# Will add aditional metrics needed
source("Functions/dataProcessing_df_analysis_3.R")

# Essentiality data
essentiality <- read_csv("Data/essentiality.csv")

# All metrics
all_metrics <- c("Gn", "Gc", "Gmax", "Gsum",
                 "N_term_thread",
                 "C_term_thread",
                 "ENT.ID",
                 "num_zipper_nc",
                 "num_loop_contacting_res",
                 "num_cross_nearest_neighbors",
                 "min_N_prot_depth_left",
                 "min_C_prot_depth_right",
                 "Travatos_G",
                 "Length")

# Entanglement metrics
ent_metrics <- c("Gn", "Gc", "Gmax", "Gsum",
                 "N_term_thread",
                 "C_term_thread",
                 "ENT.ID",
                 "num_zipper_nc",
                 "num_loop_contacting_res",
                 "num_cross_nearest_neighbors",
                 "min_N_prot_depth_left",
                 "min_C_prot_depth_right",
                 "Travatos_G")

# Max metric value
final_dfB_max <- final_df_updated %>%
  group_by(gene) %>%
  summarise(
    across(
      all_of(all_metrics),
      \(x) if (all(is.na(x))) 0 else max(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  mutate(across(all_of(all_metrics), ~replace_na(., 0)))

# Collapse to one row per gene
dfB_collapse <- final_dfB %>%
  group_by(gene) %>%
  summarise(
    Entanglement     = if_else(any(Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
    Cov_Entanglement = if_else(any(Cov_Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
    .groups = "drop"
  )

# Join to final_df_max
final_dfB_max <- final_dfB_max %>%
  left_join(dfB_collapse, by = "gene")


# Columns we will keep from B as-is; others become 0 for missing genes
keep_cols <- c("gene", "Length", "Entanglement", "Cov_Entanglement")

# Identify missing genes
missing_genes <- setdiff(unique(final_dfB$gene), unique(final_dfB_max$gene))

# Build rows for those missing genes
df_missing <- final_dfB %>%
  filter(gene %in% missing_genes) %>%
  group_by(gene) %>%
  summarise(
    Length            = suppressWarnings(max(Length, na.rm = TRUE)),
    Entanglement      = if_else(any(Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
    Cov_Entanglement  = if_else(any(Cov_Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
    .groups = "drop"
  )

# Add zero columns for all other fields that exist in final_dfB_max
cols_to_zero <- setdiff(colnames(final_dfB_max), keep_cols)

if (length(missing_genes) > 0) {
  df_missing <- df_missing %>%
    mutate(!!!setNames(rep(list(0), length(cols_to_zero)), cols_to_zero)) %>%
    select(all_of(colnames(final_dfB_max)))  # match column order
}

# Append to final_dfB_max (no-op if no missing genes)
final_dfB_max <- bind_rows(final_dfB_max, df_missing) %>%
  arrange(gene)

# Mutations
mutation_counts_af_70 <- read.csv("Data/mutation_counts_af_70_2.csv")

# Total Mutation Column
mutation_counts_af_70$Total_Mutations_n <-
  mutation_counts_af_70$pathogenic_count +
  mutation_counts_af_70$benign_count

# SASA data
af_sasa <- read.csv("Data/Human_AF_SASA.csv")

# Calculate mean SASA
unique_sasa_summary <- af_sasa %>%
  group_by(uniprotids) %>%
  summarise(
    total_SASA = sum(SASA, na.rm = TRUE),
    mean_SASA  = mean(SASA, na.rm = TRUE)
  )

# Bring in mutations
mut_join <- mutation_counts_af_70 %>%
  select(uniprotids, Total_Mutations_n) %>%
  rename(gene = uniprotids)

# Bring in mean_SASA from sasa summary
sasa_join <- unique_sasa_summary %>%
  select(uniprotids, mean_SASA) %>%
  rename(gene = uniprotids)

# Join to final_dfB_max
final_dfB_max <- final_dfB_max %>%
  left_join(mut_join, by = "gene") %>%
  left_join(sasa_join, by = "gene") %>%
  mutate(Total_Mutations_n = replace_na(Total_Mutations_n, 0))


#### Disgenet Score Percentiles ####
# Used to calculate percentiles from "raw" disgenet results
sub_finalA <- final_dfA[,c("Association_ID","score")]
sub_finalA <- unique(sub_finalA)
sub_finalA <- sub_finalA[sub_finalA$Association_ID>0,]

# Calculate score percentiles
percentiles <- quantile(sub_finalA$score, c(0.95, 0.75, 0.50))
rm(sub_finalA)

#### Contingency tables ####
# Subset necessary columns
UniScoDis <- final_dfA[, c("uniprotids", "score","Entanglement", "Cov_Entanglement")]

# Remove duplicates
UniScoDis <- distinct(UniScoDis)
# Keep Uniprots with the highest score
max_df <- UniScoDis %>%
  group_by(uniprotids) %>%
  filter(score == max(score)) %>%
  ungroup()

# Determine if the score falls in the percentile
final_df <- max_df %>%
  mutate(`50th_percentile` = ifelse(is.na(score), "No", ifelse(score >= percentiles[3], "Yes", "No")))

# Extract 50th_percentile from final_df and rename to disease
disease_join <- final_df %>%
  select(uniprotids, `50th_percentile`) %>%
  rename(gene = uniprotids, disease = `50th_percentile`)

# Join into final_dfB_max
final_dfB_max <- final_dfB_max %>%
  left_join(disease_join, by = "gene")

# Add Essential column to final_dfB_max
final_dfB_max <- final_dfB_max %>%
  left_join(
    essentiality,
    by = c("gene" = "Entry")
  ) %>%
  mutate(
    Essential = ifelse(is.na(Essential), "NT", Essential)
  )

final_dfB_max <- final_dfB_max %>% filter(!(Entanglement == "Yes" & (Travatos_G == 0 | is.na(Travatos_G))))

write.csv(final_dfB_max, "Data/final_dfB_max_af_70_raw_l_02_20.csv", row.names = FALSE)
