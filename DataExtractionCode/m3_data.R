# Data preprocessing

cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("control file input: Controlm1_3\n")
cat("types: af\n")
cat("plddt_thresh:", plddt_thresh, "\n")

# Data Checks
# Libraries
library(dplyr)
library(stringr)
library(readr)

# Mutation data (pre-cleaned, missense only)
combined_mapped_clean <- read_csv("Data/combined_mapped_clean_missense.csv")

############################################################
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
final_dfA <- read.csv("Data/final_dfA02_22_2026.csv") # Disease data
final_dfB <- read.csv("Data/final_dfB_af_02_22_2026.csv") # Entanglement data

# Filter final_dfB by using plddt
final_dfB <- final_dfB %>%
  filter(plddt >= plddt_thresh)

# Keep only matching proteins in final_dfA
final_dfA <- final_dfA %>%
  filter(uniprotids %in% final_dfB$gene)

#### Disgenet Score Percentiles ####
# Used to calculate percentiles from "raw" disgenet results
sub_finalA <- final_dfA[,c("Association_ID","score")]
sub_finalA <- unique(sub_finalA)
sub_finalA <- sub_finalA[sub_finalA$Association_ID>0,]
# Calculate score percentiles
percentiles <- quantile(sub_finalA$score, c(0.95, 0.75, 0.50))
rm(sub_finalA)

### Prepare a dataframe for contingency tables (one per disease class) ####
# Subset necessary columns
UniScoDis <- final_dfA[, c("uniprotids", "score","diseaseClasses_MSH", "Entanglement")]

# Remove duplicates
UniScoDis <- distinct(UniScoDis)

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
  for (disease_class in disease_classes) {
    # extract code part of disease_class e.g., (C01)
    code <- str_extract(disease_class, "\\([CF]\\d+\\)")
    # extract non-code part of disease_class e.g., "Infections"
    name <- gsub("\\s*\\(.*\\)", "",disease_class)
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
  distinct()

#### Prepare dataframes for Contingency tables ####

# if score>=thresh, disease == "Yes"
convert_to_yes_no <- function(value, thresh) {
  ifelse(value >= thresh, "Yes", "No")
}

# Apply the function to all columns except the first and second one (uniprotids)
UniScoDis50 <- max_scores %>%
  mutate(across(-c(0,1), ~convert_to_yes_no(., percentiles[3])))

# Summarized Disease data
final_df <- read_csv("Data/final_df_af_0_2026-2-6.csv")

# Keep only matching proteins in final_df
final_df <- final_df %>%
  filter(uniprotids %in% final_dfB$gene)

# Add Length and Essentiality
final_df <- final_df %>%
  left_join(
    final_dfB %>%
      select(gene, Length, Essential) %>%
      distinct(),
    by = c("uniprotids" = "gene")
  )

# Merge with final_df using uniprotids as the matching key
final_df <- merge(final_df,
                  combined_mapped_clean[, c("uniprotids", "User_input", "Consequences","Amino_acid_change", "Amino_acid_position", "Mutation", "Pathogenic")],
                  by = "uniprotids",
                  all.x = TRUE)

# Keep unique rows
final_df <- final_df %>% distinct()

# Remove mutations if they are located outside of the "length" of the protein
final_df <- final_df %>%
  filter(Amino_acid_position <= Length)

# Remove mutations flagged as contradictory (same User_input with pathogenic vs benign)
final_df <- final_df %>%
  filter(!User_input %in% user_inputs_to_remove)

# Keep rows with Mutations only
final_df <- final_df %>% filter(Mutation == "Yes" & Consequences == "missense")

# Count pathogenic and benign entries per uniprotid
path_counts <- final_df %>%
  group_by(uniprotids) %>%
  summarise(
    pathogenic_count = sum(Pathogenic == "Yes"),
    benign_count = sum(Pathogenic == "No")
  )

# Identify uniprotids that meet the criteria (at least 10 missense mutations)
valid_uniprotids <- path_counts %>%
  filter(pathogenic_count + benign_count >= 10) %>%
  pull(uniprotids)

# Filter the original dataset
final_df <- final_df %>%
  filter(uniprotids %in% valid_uniprotids)
