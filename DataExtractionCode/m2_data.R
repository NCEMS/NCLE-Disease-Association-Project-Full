# Data preprocessing

cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input:\n")
cat("types: af\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Data Checks
# Libraries
library(dplyr)
library(readr)

# Mutation data (pre-cleaned, missense only)
combined_mapped_clean <- read_csv("Data/combined_mapped_clean_missense.csv")

################################################################################
# Define mutation identity by gene, position, codon, and amino acid change
# to check for duplicates (same mutation recorded more than once)
################################################################################
mutation_counts <- combined_mapped_clean %>%
  group_by(uniprotids, User_input, Codon_change, Amino_acid_change, Amino_acid_position) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Keep only duplicated mutations (count > 1)
mutation_counts <- mutation_counts %>% filter(count > 1)

# Pull all rows in the main mutation table that correspond to those duplicates
# majority reported twice because they come from different sources but are the same
# removed in downstream analysis
repeated_mutations <- combined_mapped_clean %>%
  inner_join(mutation_counts, by = c("uniprotids", "User_input", "Codon_change", "Amino_acid_change", "Amino_acid_position"))

# Among duplicates, extract those from ClinVar pathogenic dataset
# (checking for contradictions where same User_input appears benign vs pathogenic)
repeated_clinvar_pathogenic <- unique(repeated_mutations %>%
                                        filter(source == "ClinVarPathogenic"))

# Save contradictory mutation identifiers (User_input) to remove later
user_inputs_to_remove <- unique(repeated_clinvar_pathogenic$User_input)


################################################################################
# analysis1_1.R -> final_dfA, final_dfB
################################################################################

final_df <- read_csv("Data/final_df_af_02_22_2026.csv") # Disease data
final_dfB <- read.csv("Data/final_dfB_af_02_22_2026.csv") # Entanglement data

# Filter final_dfB by using plddt
final_dfB <- final_dfB %>%
  filter(plddt >= plddt_thresh)

# Keep only matching proteins in final_df
final_df <- final_df %>%
  filter(uniprotids %in% final_dfB$gene)

# disease and entanglement data (protein-level)
final_dfA <- final_df
final_dfA <- final_dfA[, c("uniprotids", "score", "Entanglement", "95th_percentile", "75th_percentile", "50th_percentile")]

# length data (protein-level, from final_dfB)
length_data <- final_dfB %>%
  select(gene, Length) %>%
  distinct()

# essentiality data (protein-level)
ess_data <- final_dfB %>%
  select(gene, Essential) %>%
  distinct()

# Join length and essentiality into final_dfA (protein-level length and essentiality for each uniprotid)
final_dfA <- final_dfA %>%
  left_join(length_data, by = c("uniprotids" = "gene")) %>%
  left_join(ess_data,    by = c("uniprotids" = "gene"))

# Join mutation data onto final_dfA by uniprotids
final_dfA <- merge(final_dfA,
                   combined_mapped_clean[, c("uniprotids", "User_input", "Consequences","Amino_acid_change", "Amino_acid_position", "Mutation", "Pathogenic")],
                   by = "uniprotids",
                   all.x = TRUE)
# Keep unique rows
final_dfA <- final_dfA %>% distinct()

# Remove mutations if they are located outside of the "length" of the protein
final_dfA <- final_dfA %>%
  filter(Amino_acid_position <= Length)

# Remove mutations flagged as contradictory (same User_input with pathogenic vs benign)
final_dfA <- final_dfA %>% filter(!User_input %in% user_inputs_to_remove)

###############################################################################
# Build final_dfA for modeling
# - Restrict to missense mutations
# - Define Mis_Pathogenic at mutation level
# - Restrict to proteins with at least 10 missense mutations (pathogenic + benign)
###############################################################################
final_dfA <- final_dfA %>% filter(Mutation=="Yes" & Consequences=="missense")

final_dfA$Mis_Pathogenic <- ifelse(
  final_dfA$Pathogenic == "Yes",
  "Yes","No"
)

# Count pathogenic and benign missense mutations per protein
path_counts <- final_dfA %>%
  group_by(uniprotids) %>%
  summarise(
    pathogenic_count = sum(Pathogenic=="Yes"),
    benign_count = sum(Pathogenic=="No")
  )

# write.csv(path_counts, "Data/mutation_counts_af_70_2.csv", row.names = FALSE)

# Count pathogenic and benign missense mutations per protein
valid_uniprotids <- path_counts %>%
  filter(pathogenic_count + benign_count >= 10) %>%
  pull(uniprotids)

# Filter for proteins with at least 10 missense mutations
final_dfA <- final_dfA %>% filter(uniprotids %in% valid_uniprotids)

###############################################################################
# Aggregate protein summary
# - One row per uniprotid
# - Entanglement = "Yes" if any row for that protein is entangled
# - Total_Mutations_n = # missense mutations
# - Path_Mutations_m = # pathogenic missense mutations
# - p_hat_m_div_n = proportion of missense mutations that are pathogenic
###############################################################################
protein_summary <- final_dfA %>%
  group_by(uniprotids) %>%
  summarise(
    Entanglement = ifelse(any(Entanglement == "Yes", na.rm = TRUE), "Yes", "No"),
    Length = first(Length),
    Essential = first(Essential),
    `50th_percentile` = first(`50th_percentile`),
    Total_Mutations_n = sum(Mutation == "Yes"),
    Path_Mutations_m = sum(Mis_Pathogenic == "Yes", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    p_hat_m_div_n = if_else(
      Total_Mutations_n > 0,
      Path_Mutations_m / Total_Mutations_n,
      NA_real_
    )
  ) %>%
  mutate(
    Pathogenic = ifelse(Path_Mutations_m > 0, "Yes", "No")
  )



# Set factors for analysis
protein_summary$Entanglement <- factor(
  protein_summary$Entanglement, levels=c("No","Yes")
)

protein_summary$Essential <- factor(protein_summary$Essential, levels = c("No", "NT", "Yes"))

# Check counts
table(Essential = protein_summary$Essential, Entanglement = protein_summary$Entanglement)

protein_summary %>%
  group_by(Pathogenic, Essential, Entanglement) %>%
  summarise(N = n(), .groups = "drop")

# Save
write.csv(
  protein_summary,
  paste0("Data/DataRanInAnalysis/protein_summary_af_", plddt_thresh, "_", Sys.Date(), ".csv"),
  row.names = FALSE
)

rm(list = setdiff(ls(), c("protein_summary", "plddt_thresh", "disease_only")))


