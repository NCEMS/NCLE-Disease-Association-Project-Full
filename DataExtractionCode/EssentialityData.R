# Essentiality data extraction

# Libraries
library(dplyr)
library(tidyr)

# Load CRISPR Essential Gene Data
cris_ess <- read.csv("Data/CRISPRInferredCommonEssentials.csv")
tested_genes <- unique(colnames(read.csv("Data/CRISPRGeneEffect.csv",
                                         nrows = 0, check.names = FALSE)))

# Check overlap between essential genes and tested genes
length(intersect(cris_ess$Essentials,tested_genes))

# Clean and Parse Essential Gene Data
cris_ess <- cris_ess %>%
  separate(Essentials, into = c("GeneName", "GeneID"), sep = " \\(", fill = "right") %>%
  mutate(
    GeneID = gsub("\\)", "", GeneID),         # clean closing parenthesis
    # if GeneName is just digits, move it into GeneID
    GeneID = ifelse(grepl("^[0-9]+$", GeneName) & (is.na(GeneID) | GeneID == ""),
                    GeneName, GeneID),
    GeneName = ifelse(grepl("^[0-9]+$", GeneName), NA, GeneName),
    Essential = "Yes"
  )

# Clean and Parse Tested Gene Data
tested <- tibble(tested_genes = tested_genes) %>%
  separate(tested_genes, into = c("GeneName", "GeneID"), sep = " \\(", fill = "right") %>%
  mutate(
    GeneID = gsub("\\)", "", GeneID),         # remove closing parenthesis
    GeneID = ifelse(grepl("^[0-9]+$", GeneName) & (is.na(GeneID) | GeneID == ""),
                    GeneName, GeneID),
    GeneName = ifelse(grepl("^[0-9]+$", GeneName), NA, GeneName)
  )

# Extract GeneID and save (to map to Uniprot ID)
# write.table(
#   cris_ess$GeneID,
#   file = "results/dataframes/gene_ids_cris.txt",
#   row.names = FALSE,
#   col.names = FALSE,
#   quote = FALSE
# )

# write.table(
#   tested$GeneID,
#   file = "results/dataframes/gene_tested.txt",
#   row.names = FALSE,
#   col.names = FALSE,
#   quote = FALSE
# )

# Load UniProt ID Mapping Files
idmapping2 <- read.delim("Data/idmapping_2_reviewed_true_AND_model_organ_2025_10_03.tsv",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

idmapping3 <- read.delim("Data/idmapping_2025_10_14.tsv",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

idmapping3<- idmapping3 %>%
  filter(tolower(Reviewed) == "reviewed")

# Map essential genes to UniProt entries
essentiality_data_full <- merge(
  cris_ess,
  idmapping2[, c("From", "Entry", "Gene.Names")],
  by.x = "GeneID",
  by.y = "From",
  all.x = TRUE
)

# Map tested genes to UniProt entries
essentiality_data_full2 <- merge(
  tested,
  idmapping3[, c("From", "Entry", "Gene.Names")],
  by.x = "GeneID",
  by.y = "From",
  all.x = TRUE
)

# Check whether GeneName appears in UniProt Gene.Names
essentiality_data_full$GeneName_match <- ifelse(
  mapply(function(g, names) grepl(paste0("\\b", g, "\\b"), names), 
         essentiality_data_full$GeneName, 
         essentiality_data_full$Gene.Names),
  "Yes", 
  "No"
)

essentiality_data_full2$GeneName_match <- ifelse(
  mapply(function(g, names) grepl(paste0("\\b", g, "\\b"), names), 
         essentiality_data_full2$GeneName, 
         essentiality_data_full2$Gene.Names),
  "Yes", 
  "No"
)

# Mark tested genes as Essential if their UniProt Entry
# appears in the essential gene list
essentiality_data_full2 <- essentiality_data_full2 %>%
  mutate(
    Essential = ifelse(Entry %in% essentiality_data_full$Entry, "Yes", "No")
  )

# If multiple mappings exist for one Entry,
# mark Essential as "Yes" if ANY row is "Yes"
essentiality_data_full2_collapsed <- essentiality_data_full2 %>%
  group_by(Entry) %>%
  summarise(Essential = ifelse(any(Essential == "Yes"), "Yes", "No"), .groups = "drop")

# Save
write.csv(essentiality_data_full2_collapsed, "Data/essentiality.csv", row.names = FALSE)
