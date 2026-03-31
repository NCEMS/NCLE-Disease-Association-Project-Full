# Data preprocessing

cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("type: af\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Libraries
library(dplyr)                        # data manipulation (pipes, group_by, summarise, joins)
library(stringr)                      # string utilities (not directly used in this snippet, but used in sourced functions)
library(readr)                        # fast CSV reading (read_csv)
library(data.table)                   # fast data.table operations (as.data.table, merges, by-group computations)


# Entanglement and Disease data
source("Functions/dataProcessing_df_mutation.R")

# Mutation data (pre-cleaned, missense only)
combined_mapped_clean <- read_csv("Data/combined_mapped_clean_missense.csv")

################################################################################
# Define mutation identity by gene, position, codon, and amino acid change
# to check for duplicates (same mutation recorded more than once)
################################################################################
mutation_counts <- combined_mapped_clean %>%                              # start from full mutation table
  group_by(uniprotids, User_input, Codon_change, Amino_acid_change, Amino_acid_position) %>% # define “same mutation” key
  summarise(count = n()) %>%                                              # count how many times each mutation key appears
  arrange(desc(count))                                                    # sort so most repeated are at the top

# Keep only duplicated mutations (count > 1)
mutation_counts <- mutation_counts %>% filter(count > 1)

# Pull all rows in the main mutation table that correspond to those duplicates
# majority reported twice because they come from different sources but are the same
# removed in downstream analysis
repeated_mutations <- combined_mapped_clean %>%
  inner_join(mutation_counts,                                             # keep only rows matching duplicated keys
             by = c("uniprotids", "User_input", "Codon_change", "Amino_acid_change", "Amino_acid_position"))

# Among duplicates, extract those from ClinVar pathogenic dataset
# (checking for contradictions where same User_input appears benign vs pathogenic)
repeated_clinvar_pathogenic <- unique(repeated_mutations %>%
                                        filter(source == "ClinVarPathogenic"))

# Save contradictory mutation identifiers (User_input) to remove later
user_inputs_to_remove <- unique(repeated_clinvar_pathogenic$User_input)


ent_dat_clean <- af_ent_dat_clean


################################################################################
# analysis1_1.R -> final_dfA, final_dfB
################################################################################

final_dfA <- read_csv("Data/final_df_af_02_22_2026.csv") # Disease data
final_dfB <- read.csv("Data/final_dfB_af_02_22_2026.csv") # Entanglement data

# Filter final_dfB by using plddt
final_dfB <- final_dfB %>%
  filter(plddt >= plddt_thresh)

# Keep only matching proteins in final_df
final_dfA <- final_dfA %>%
  filter(uniprotids %in% final_dfB$gene)

# disease and entanglement data (protein-level)
final_dfA <- final_dfA[, c("uniprotids", "score", "Entanglement", "95th_percentile", "75th_percentile", "50th_percentile", "Length", "Essential")]

# length and essentiality data (protein-level)
length_data <- final_dfA %>%
  select(uniprotids, Length, Essential) %>%
  distinct()

# Join mutation data onto final_dfA by uniprotids, matching rows from combined_mapped_clean are added when available
final_dfA <- merge(final_dfA,
                   combined_mapped_clean[, c("uniprotids", "User_input", "Consequences","Amino_acid_change", "Amino_acid_position", "Mutation", "Pathogenic")],
                   by = "uniprotids",
                   all.x = TRUE)                                        # keep all proteins even if they have no mutations

# Join entanglement region
final_dfA <- final_dfA %>%
  left_join(ent_dat_clean[, c("uniprotids", "numeric_EntRegion")],       # bring in list/vector of entangled residue positions per protein
            by = "uniprotids")

final_dfA <- final_dfA %>% distinct()

# Remove mutations flagged as contradictory (same User_input with pathogenic vs benign)
final_dfA <- final_dfA %>%
  filter(!User_input %in% user_inputs_to_remove)                        # exclude any rows whose User_input is in the removal list

# Quick check: show User_input values with >1 rows remaining (if any)
final_dfA %>%
  group_by(User_input) %>%
  summarise(count = n()) %>%
  filter(count > 1)

# Print all duplicated User_input rows, amino acid change is different so not duplicates
rows_dup <- final_dfA %>%
  filter(!is.na(User_input)) %>%                                        # ignore missing User_input rows
  group_by(User_input) %>%                                              # group by User_input
  filter(n() > 1) %>%                                                   # keep groups with >1 row
  ungroup() %>%                                                         # remove grouping
  print(n = Inf)                                                        # print all rows to console

# Create missense pathogenic flag at the mutation level
final_dfA$Mis_Pathogenic <- ifelse(final_dfA$Pathogenic == "Yes",
                                   "Yes",
                                   "No")
# Count pathogenic and benign missense mutations per uniprotid
path_counts <- final_dfA %>%
  group_by(uniprotids) %>%
  summarise(
    pathogenic_count = sum(Pathogenic == "Yes"),
    benign_count = sum(Pathogenic == "No" & Mutation == "Yes")
  )

# Keep proteins with at least 10 total missense mutations (pathogenic + benign)
valid_uniprotids <- path_counts %>%
  filter(pathogenic_count + benign_count >= 10) %>%
  pull(uniprotids)                                                      # extract vector of UniProt IDs

final_dfA <- final_dfA %>%
  filter(uniprotids %in% valid_uniprotids)

############################################################################
# Add residue data: build a full residue grid (1:Length) per valid protein
############################################################################
dt <- as.data.table(length_data)                                        # convert length_data to data.table for fast operations

# One row per residue per uniprotid: resid = 1,2,...,Length
res_dt <- dt[, .(resid = 1:Length), by = uniprotids]                    # expand each protein into a residue-level grid

# Keep only valid_uniprotids (≥10 missense mutations)
res_dt <- res_dt[uniprotids %in% valid_uniprotids]

# Convert final_dfA to data.table for fast merging
final_dfA_dt <- as.data.table(final_dfA)

# Collapse to one protein-level row per uniprotid
# (score, entanglement flags, percentiles, length, entangled region list)
gene_dt <- final_dfA_dt[
  ,
  .(
    score             = score[1],
    Entanglement      = Entanglement[1],
    `95th_percentile` = `95th_percentile`[1],
    `75th_percentile` = `75th_percentile`[1],
    `50th_percentile` = `50th_percentile`[1],
    Length            = Length[1],
    Essential         = Essential[1],
    numeric_EntRegion = list(unique(unlist(numeric_EntRegion)))         # store entangled region residue positions as a unique vector list
  ),
  by = uniprotids
]

# Attach protein-level info to every residue row
res_gene_dt <- gene_dt[                                                 # protein-level table with entanglement regions
  res_dt,                                                               # residue-level grid
  on = .(uniprotids)                                                    # join key: uniprotids
]

# Keep only mutation/residue-level columns needed later
mut_cols <- c("uniprotids", "User_input", "Consequences",
              "Amino_acid_change", "Amino_acid_position",
              "Mutation", "Pathogenic", "Mis_Pathogenic")

mut_dt <- final_dfA_dt[, ..mut_cols]                                    # subset mutation-level columns from final_dfA_dt

# Collapse to one mutation per (uniprotids, Amino_acid_position) with tie-breaking:
# 1) If any Mis_Pathogenic == "Yes", keep the first such row
# 2) Otherwise, keep the first row at that position
mut_dt_unique <- mut_dt[                                                # start from mutation rows
  ,
  {
    idx_yes <- which(Mis_Pathogenic == "Yes")                           # indices of pathogenic mutations at this protein+position
    if (length(idx_yes) > 0L) {                                         # if any pathogenic mutation exists at that position
      .SD[idx_yes[1L]]                                                  # keep the first pathogenic row
    } else {                                                            # otherwise (all benign/No)
      .SD[1L]                                                           # keep the first row
    }
  },
  by = .(uniprotids, Amino_acid_position)                               # group by protein and residue position
]

# Merge residue grid (all residues) with mutation info (per residue, if present)
res_merged_dt <- mut_dt_unique[                                         # left table: unique mutation info by residue position
  res_gene_dt,                                                          # right table: full residue grid with protein-level info
  on = .(uniprotids, Amino_acid_position = resid)                       # join on protein ID and map resid -> Amino_acid_position
]

# For residues with no mutation, fill mutation-related fields with "No"
res_merged_dt[, c("Mutation", "Pathogenic", "Mis_Pathogenic") :=        # update these three columns in place
                lapply(.SD, function(x) fifelse(is.na(x), "No", x)),     # replace NA with "No"
              .SDcols = c("Mutation", "Pathogenic", "Mis_Pathogenic")]  # only apply to these columns

# Filter to entangled proteins only
res_merged_dt <- res_merged_dt[Entanglement == "Yes"]

# Compute in_region
res_merged_dt[, in_region := {                                          # create in_region column by protein
  region <- numeric_EntRegion[[1L]]                                     # extract region vector/list (stored as list column) for this protein
  region <- unlist(region)                                              # flatten
  region <- region[!is.na(region)]                                      # remove NA residue indices
  fifelse(Amino_acid_position %in% region, "Yes", "No")                 # mark residue as in-region if its position is in entangled region
}, by = uniprotids]

saveRDS(
  res_merged_dt,
  paste0("Data/DataRanInAnalysis/res_merged_dt_af_", plddt_thresh, "_", Sys.Date(), ".rds")
)

res_merged_dt_csv <- res_merged_dt
res_merged_dt_csv$numeric_EntRegion <- sapply(
  res_merged_dt_csv$numeric_EntRegion,
  function(x) paste(x, collapse = "; ")
)

write.csv(
  res_merged_dt_csv,
  paste0(
    "Data/DataRanInAnalysis/res_merged_dt_af_", plddt_thresh, "_", Sys.Date(), ".csv"
  ),
  row.names = FALSE
)

rm(list = setdiff(ls(), c("res_merged_dt", "plddt_thresh")))
