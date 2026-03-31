# Updated data set code

# Gather unique Uniprot IDs in the data frame and create a list. Save the Uniprots to a txt file and search
# their Entrez IDs in Uniprot. Include Uniprot IDs with Entrez IDs in the data frame (list: Uni_to_Ent).
# Query Entrez IDs through Disgenet and add all other columns. Start with the Disgenet output, subset
# relevant columns (Uniprot, Entrez ID, Disease Class Association, Score, disease_name, diseaseUMLSCUI) as
# B2 data frame, and add an association ID. Modify B2 to have one Uniprot ID per row and remove rows not
# in Uni_to_Ent. Add a blank Entanglement column to B2. Add Uniprot IDs in Uni_to_Ent but not in B2 with
# NA columns and Association ID = -1. If Uniprot ID has Gn or Gc value, set Entanglement to Yes, otherwise No.


# Gather unique Uniprot IDs in the data frame and create a list. Save the Uniprots to a txt file and
# search their Entrez IDs in Uniprot. Include Uniprot IDs with Entrez IDs in the data frame (list: Uni_to_Ent).
# Query Entrez IDs through Disgenet and add all other columns. Start with the new Entanglement file,
# subset gene, ENT.ID, Gn, Gc, N_term_thread, and C_term_thread columns, and remove genes without Entrez ID.

# Clearing environment
remove(list=ls())

# Loading libraries
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(gmodels)
library(readxl)
library(readr)
library(progress)
library(disgenet2r)

## Files
# Reading entanglement features data and length data
data <- read_delim(
  file = "Data/Human_AF_combined_20250521 (1).txt",
  delim = "|",
  trim_ws = TRUE,
  col_names = TRUE
)

# Make CCBond either True or False only
data <- data %>%
  mutate(CCBond = case_when(
    CCBond %in% c("True","1.0", "1") ~ "True",
    CCBond %in% c("False","0.0", "0") ~ "False",
    TRUE ~ NA
  ))

# CC count summary
cc_map <- data %>%
  select(gene, CCBond) %>%
  mutate(
    CCBond_flag = str_to_lower(str_trim(as.character(CCBond))) %in%
      "true"
  ) %>%
  group_by(gene) %>%
  summarise(cc_count = sum(CCBond_flag, na.rm = TRUE), .groups = "drop")

# Join back on gene and subtract from ENT.ID (because they will be used for counts in later analysis)
data <- data %>%
  left_join(cc_map, by = "gene") %>%
  mutate(
    cc_count = replace_na(cc_count, 0),
    `ENT-ID`   = ifelse(is.na(`ENT-ID`), NA_real_, as.numeric(`ENT-ID`) - cc_count + 1)
  )

crystal_data <- read_delim(
  "Data/Human_EXP_combined_20250614.txt",
  delim = "|",
  trim_ws = TRUE,
  col_names = TRUE,
  col_types = cols(
    Knot_type = col_character()
  )
)

# Make CCBond either True or False only
crystal_data <- crystal_data %>%
  mutate(CCBond = case_when(
    CCBond %in% c("True","1.0", "1") ~ "True",
    CCBond %in% c("False","0.0", "0") ~ "False",
    TRUE ~ NA
  ))

# CC count summary
cc_map_crys <- crystal_data %>%
  select(gene, CCBond) %>%
  mutate(
    CCBond_flag = str_to_lower(str_trim(as.character(CCBond))) %in%
      "true"
  ) %>%
  group_by(gene) %>%
  summarise(cc_count = sum(CCBond_flag, na.rm = TRUE), .groups = "drop")

# Join back on gene and subtract from ENT.ID
crystal_data <- crystal_data %>%
  left_join(cc_map_crys, by = "gene") %>%
  mutate(
    cc_count = replace_na(cc_count, 0),
    `ENT-ID`   = ifelse(is.na(`ENT-ID`), NA_real_, as.numeric(`ENT-ID`) - cc_count)
  )
# Load length and Trovato data
length <- read_excel("Data/Human_mappedCorrected_Lengths.xlsx")
travato_af <- read_excel("Data/Trovato.xlsx")
travato_crystal <- read.csv("Data/human_Xtal_maxGLN_Trovato.csv", sep = "|")

# Rename columns
data <- data %>%
  rename(
    ENT.ID = `ENT-ID`,
    Gn     = gn,
    Gc     = gc
  )

crystal_data <- crystal_data %>%
  rename(
    ENT.ID = `ENT-ID`,
    Gn     = gn,
    Gc     = gc
  )

# Add Entanglemet and Knot column
data$Entanglement <- ifelse((!is.na(data$Gn) & abs(data$Gn) >= 0.6 & (data$CCBond == "False" | is.na(data$CCBond))) |
                              (!is.na(data$Gc) & abs(data$Gc) >= 0.6 & (data$CCBond == "False" | is.na(data$CCBond))),
                            "Yes", "No")
data <- data %>%
  mutate(Knot = ifelse(!is.na(Knot_type) & Knot_type != "", "Yes", "No"))

crystal_data$Entanglement <- ifelse((!is.na(crystal_data$Gn) & abs(crystal_data$Gn) >= 0.6 & (crystal_data$CCBond == "False" | is.na(crystal_data$CCBond))) |
                                      (!is.na(crystal_data$Gc) & abs(crystal_data$Gc) >= 0.6 & (crystal_data$CCBond == "False" | is.na(crystal_data$CCBond))),
                                    "Yes", "No")

crystal_data <- crystal_data %>%
  mutate(Knot = ifelse(!is.na(Knot_type) & Knot_type != "", "Yes", "No"))

# Change name to Travtos_G and subset necessary columns
names(travato_af)[names(travato_af) == "maxGLN"] <- "Travatos_G"
names(travato_af)[names(travato_af) == "gene"] <- "uniprot"
travato_af <- travato_af[, c("uniprot", "Travatos_G")]

names(travato_crystal)[names(travato_crystal) == "maxGLN"] <- "Travatos_G"
names(travato_crystal)[names(travato_crystal) == "gene"] <- "uniprot"
travato_crystal <- travato_crystal[, c("uniprot", "Travatos_G")]


# Data Frame A

## Gather all unique Uniprot IDs that will be in the data frame.
# Create a list of all Uniprot IDs

uniprots <- unique(data$gene)
uniprots_crystal <- unique(crystal_data$gene)

# Save the Uniprots in a txt file and search their Entrez ID in Uniprot

#write.table(uniprots, file = "uniprots2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(uniprots_crystal, file = "uniprots_crystal.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

TransData <- read.table("Data/idmapping_2024_05_21.tsv", header = T)
TransData_Crystal <- TransData %>%
  filter(From %in% crystal_data$gene)


# Include Uniprot IDs that have an Entrez ID in the data frame (call them list Uni_to_Ent).

Uni_to_Ent <- unique(TransData$From)
Uni_to_Ent_Crystal <- unique(TransData_Crystal$From)

# Query Entrez IDs through Disgenet

api_key <- "6c0567e3-c666-4d3c-9a50-9d027e1e6517"
Sys.setenv(DISGENETPLUS_API_KEY= api_key)

# Function to query genes in batches and combine the results
query_genes <- function(gene_ids) {
  all_gene_results <- list()

  # Split gene_ids into batches
  gene_batches <- split(gene_ids, ceiling(seq_along(gene_ids) / 100))

  # Initialize progress bar
  pb <- progress_bar$new(total = length(gene_batches), format = "[:bar] :percent :elapsed - ETA: :eta")

  # Iterate over each batch of genes
  for (i in seq_along(gene_batches)) {
    batch <- gene_batches[[i]]
    # Convert gene IDs to numeric
    batch <- as.numeric(batch)

    # Query genes in the batch
    results <- gene2disease(gene = batch, vocabulary = "ENTREZ", database = "CURATED")

    # Check if all results are "no results for query"
    if (all(class(results) == "character")) {
      # If all genes have no results, create placeholder rows for all genes
      for (gene_id in batch) {
        placeholder_row <- data.frame(
          gene_symbol = NA,
          geneid = gene_id,
          ensemblid = NA,
          geneNcbiType = NA,
          geneDSI = NA,
          geneDPI = NA,
          uniprotids = NA,
          protein_classid = NA,
          protein_class_name = NA,
          disease_name = NA,
          diseaseType = NA,
          diseaseUMLSCUI = NA,
          diseaseClasses_MSH = NA,
          diseaseClasses_UMLS_ST = NA,
          diseaseClasses_DO = NA,
          diseaseClasses_HPO = NA,
          score = NA,
          yearInitial = NA,
          yearFinal = NA,
          diseaseid = NA,
          DisGeNET = "No"
        )
        all_gene_results <- c(all_gene_results, list(placeholder_row))
      }
    } else {
      results <- results@qresult
      # Extract only the columns we need from results
      required_columns <- c(
        "gene_symbol", "geneid", "ensemblid", "geneNcbiType",
        "geneDSI", "geneDPI", "uniprotids", "protein_classid",
        "protein_class_name", "disease_name", "diseaseType",
        "diseaseUMLSCUI", "diseaseClasses_MSH", "diseaseClasses_UMLS_ST",
        "diseaseClasses_DO", "diseaseClasses_HPO", "score",
        "yearInitial", "yearFinal", "diseaseid"
      )
      results <- results[, required_columns]
      results$DisGeNET <- "Yes"

      # Append results to all_gene_results
      all_gene_results <- c(all_gene_results, list(results))
    }

    # Update progress bar
    pb$tick()
  }

  # Combine all the results into a single data frame
  final_results <- do.call(rbind, all_gene_results)

  return(final_results)
}

# Track the time taken to execute the query_genes function
start_time <- Sys.time()

final_results <- as.data.frame(query_genes(as.numeric(unique(TransData$To))))
final_results_Crystal <- as.data.frame(query_genes(as.numeric(unique(TransData_Crystal$To))))

end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Output the elapsed time
cat("Elapsed time:", elapsed_time, "\n")

# Convert columns in list form to character vector
convert_list_columns <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) paste(x, collapse = ", "))
    }
  }
  return(df)
}

final_results <- convert_list_columns(final_results)

# Save file
# write.csv(final_results, file = "Data/final_resultsRAW05_21_2024.csv", row.names = FALSE)

# Open file
# final_results <- read.csv("Data/final_resultsRAW05_21_2024.csv")

# Remove pseudo and NA
final_results <- subset(final_results, geneNcbiType == "protein-coding")

# Replace NA and empty values in the diseaseClasses_MSH column
final_results <- final_results %>%
  mutate(diseaseClasses_MSH = ifelse(is.na(diseaseClasses_MSH) | diseaseClasses_MSH == "", "Uncategorized", diseaseClasses_MSH))

# Convert columns with list-like objects to character vectors
convert_list_columns <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) paste(x, collapse = ", "))
    }
  }
  return(df)
}

# Apply the function to convert list-like columns
final_results <- convert_list_columns(final_results)

## Add all other columns
# Start with output file from Disgenet

head(final_results)

# Subset the Uniprot, Entrez ID, Disease Class Association, Score, disease_name, and diseaseUMLSCUI from the Disgenet file (call it B2 data frame)
# Add an association id

B2 <- final_results %>%
  select(uniprotids, geneid, diseaseClasses_MSH, score,disease_name,diseaseUMLSCUI,protein_class_name,protein_classid)
maxB2 <- nrow(B2)
B2$Association_ID <- seq(1, maxB2)
colnames(B2)[colnames(B2) == "geneid"] <- "EntreID"


# Modify B2 in order for the data frame to have one Uniprot ID per row

B2 <- separate_rows(B2, uniprotids, sep = ", ")

# Identify duplicate rows based on all columns except the last one
dup_rows <- duplicated(B2[, -ncol(B2)]) | duplicated(B2[, -ncol(B2)], fromLast = TRUE)

# Subset duplicate rows
duplicate_B2 <- B2[dup_rows, ]

# Remove rows with Uniprots that are not in list Uni_to_Ent

B2 <- B2 %>%
  filter(uniprotids %in% Uni_to_Ent)
B2_Crystal <- B2 %>%
  filter(uniprotids %in% Uni_to_Ent_Crystal)


# Add the Uniprot IDs that are in list Uni_to_Ent  and but not in B2. Make the Entanglement column empty for these, add their corresponding Entrez ID/s from TransData, make Disease Class Association NA and Score from Disgenet 0 for the Uniprots we are adding.

names(TransData) <- c("uniprotids", "EntreID")
names(TransData_Crystal) <- c("uniprotids", "EntreID")

# Create a lookup table with concatenated EntreIDs for each uniprotid
lookup_table <- TransData %>%
  group_by(uniprotids) %>%
  summarise(EntreID = paste(unique(EntreID), collapse = ", "))

lookup_table_Crystal <- TransData_Crystal %>%
  group_by(uniprotids) %>%
  summarise(EntreID = paste(unique(EntreID), collapse = ", "))

# Merge the lookup table with B2 to replace the EntreID column
B2 <- B2 %>%
  left_join(lookup_table, by = "uniprotids") %>%
  select(-EntreID.x) %>%
  rename(EntreID = EntreID.y)

B2_Crystal <- B2_Crystal %>%
  left_join(lookup_table_Crystal, by = "uniprotids") %>%
  select(-EntreID.x) %>%
  rename(EntreID = EntreID.y)

# Find Uniprot IDs in `Uni_to_Ent` that are not in `B2`
uniprots_to_add <- setdiff(Uni_to_Ent, B2$uniprotids)
uniprots_to_add_Crystal <- setdiff(Uni_to_Ent_Crystal, B2_Crystal$uniprotids)

# Create a data frame for the additional Uniprots
additional_uniprots <- data.frame(
  uniprotids = uniprots_to_add,
  EntreID = sapply(uniprots_to_add, function(id) {
    lookup_table$EntreID[lookup_table$uniprotids == id]
  }),
  diseaseClasses_MSH = NA,
  score = 0,
  disease_name = NA,
  diseaseUMLSCUI = NA,
  Association_ID = -1 ,
  Entanglement = "",
  protein_class_name = "Other/None",
  protein_classid = "Other/None"
)

additional_uniprots_Crystal <- data.frame(
  uniprotids = uniprots_to_add_Crystal,
  EntreID = sapply(uniprots_to_add_Crystal, function(id) {
    lookup_table_Crystal$EntreID[lookup_table_Crystal$uniprotids == id]
  }),
  diseaseClasses_MSH = NA,
  score = 0,
  disease_name = NA,
  diseaseUMLSCUI = NA,
  Association_ID = -1 ,
  Entanglement = "",
  protein_class_name = "Other/None",
  protein_classid = "Other/None"
)

# Combine the additional Uniprots with the B2 data frame
B2$EntreID <- as.character(B2$EntreID)
B2 <- bind_rows(B2, additional_uniprots)

B2_Crystal$EntreID <- as.character(B2_Crystal$EntreID)
B2_Crystal <- bind_rows(B2_Crystal, additional_uniprots_Crystal)

# If the Uniprot ID has Gn or Gc value, then make the Entanglement column equal Yes, otherwise No.
B2 <- B2 %>%
  select(-Entanglement)
B2_Crystal <- B2_Crystal %>%
  select(-Entanglement)

# Summary of the Entanglement column in 'data' by uniprotids
entanglement_summary <- data %>%
  group_by(gene) %>%
  summarise(Entanglement = ifelse(any(Entanglement == "Yes"), "Yes", "No"))
entanglement_summary_crystal <- crystal_data %>%
  group_by(gene) %>%
  summarise(Entanglement = ifelse(any(Entanglement == "Yes"), "Yes", "No"))

# Merge summary with B2
B2 <- B2 %>%
  left_join(entanglement_summary, by = c("uniprotids" = "gene"))
B2_Crystal <- B2_Crystal %>%
  left_join(entanglement_summary_crystal, by = c("uniprotids" = "gene"))


final_dfA <- B2
final_dfA_Crystal <- B2_Crystal

# Data Frame B

# Remove genes not in List Uni_to_Ent

subset_df <- data %>%
  filter(gene %in% Uni_to_Ent)
subset_df_Crystal <- crystal_data %>%
  filter(gene %in% Uni_to_Ent_Crystal)

# Add length column
subset_df <- subset_df %>%
  left_join(
    length %>% select(gene, Length = `AF Length`),
    by = "gene"
  )
subset_df_Crystal <- subset_df_Crystal %>%
  left_join(length %>% select(gene, Length = `Crystal Mapped Length`),
            by = "gene")

# Add Travatos G column
subset_df <- subset_df %>%
  left_join(travato_af, by = c("gene" = "uniprot"))

subset_df_Crystal <- subset_df_Crystal %>%
  left_join(travato_crystal, by = c("gene" = "uniprot"))

# Update Travatos_G to 0 if NA and no entanglement
subset_df <- subset_df %>%
  mutate(Travatos_G = ifelse(Entanglement == "No" & is.na(Travatos_G), 0, Travatos_G))

subset_df_Crystal <- subset_df_Crystal %>%
  mutate(Travatos_G = ifelse(Entanglement == "No" & is.na(Travatos_G), 0, Travatos_G))


final_dfB <- subset_df
final_dfB_Crystal <- subset_df_Crystal

final_dfB <- final_dfB %>%
  mutate(
    Cov_Entanglement = if_else(
      (CCBond %in% c(TRUE, 1, "1", "TRUE", "True", "Yes", "yes")) &
        (abs(Gn) >= 0.6 | abs(Gc) >= 0.6),
      "Yes", "No", missing = "No"
    )
  )

final_dfB_Crystal <- final_dfB_Crystal %>%
  mutate(
    Cov_Entanglement = if_else(
      (CCBond %in% c(TRUE, 1, "1", "TRUE", "True", "Yes", "yes")) &
        (abs(Gn) >= 0.6 | abs(Gc) >= 0.6),
      "Yes", "No", missing = "No"
    )
  )

# Summary of the Cov Entanglement column in 'data' by uniprotids
cov_entanglement_summary <- final_dfB %>%
  group_by(gene) %>%
  summarise(Cov_Entanglement = ifelse(any(Cov_Entanglement == "Yes"), "Yes", "No"))
cov_entanglement_summary_crystal <- final_dfB_Crystal %>%
  group_by(gene) %>%
  summarise(Cov_Entanglement = ifelse(any(Cov_Entanglement == "Yes"), "Yes", "No"))

# Merge summary with B2
final_dfA <- final_dfA %>%
  left_join(cov_entanglement_summary, by = c("uniprotids" = "gene"))
final_dfA_Crystal <- final_dfA_Crystal %>%
  left_join(cov_entanglement_summary_crystal, by = c("uniprotids" = "gene"))

# Summary of the Knot column in 'data' by uniprotids
knot_summary <- final_dfB %>%
  group_by(gene) %>%
  summarise(Knot = ifelse(any(Knot == "Yes"), "Yes", "No"))
knot_summary_crystal <- final_dfB_Crystal %>%
  group_by(gene) %>%
  summarise(Knot = ifelse(any(Knot == "Yes"), "Yes", "No"))

# Merge summary with B2
final_dfA <- final_dfA %>%
  left_join(knot_summary, by = c("uniprotids" = "gene"))
final_dfA_Crystal <- final_dfA_Crystal %>%
  left_join(knot_summary_crystal, by = c("uniprotids" = "gene"))

final_dfB_Crystal$CCBond <- as.character(final_dfB_Crystal$CCBond)

# essential data
essential <- read.csv("Data/essentiality.csv")

# Join the essential summary by uniprotids
final_dfB <- final_dfB %>%
  # Join the Essential info by UniProt Entry
  left_join(essential[, c("Entry", "Essential")], by = c("gene" = "Entry")) %>%
  # Replace NA in Essential with "NT"
  mutate(Essential = ifelse(is.na(Essential), "NT", Essential))

final_dfB_Crystal <- final_dfB_Crystal %>%
  # Join the Essential info by UniProt Entry
  left_join(essential[, c("Entry", "Essential")], by = c("gene" = "Entry")) %>%
  # Replace NA in Essential with "NT"
  mutate(Essential = ifelse(is.na(Essential), "NT", Essential))

# Remove due to updated canonical sequence lengths in Uniprot
final_dfB <- final_dfB %>%
  group_by(gene) %>%
  filter(
    all(
      (min_N_prot_depth_left  > 0 | is.na(min_N_prot_depth_left)) &
        (min_C_prot_depth_right > 0 | is.na(min_C_prot_depth_right))
    )
  ) %>%
  ungroup()

final_dfB_Crystal <- final_dfB_Crystal %>%
  group_by(gene) %>%
  filter(
    all(
      (min_N_prot_depth_left  > 0 | is.na(min_N_prot_depth_left)) &
        (min_C_prot_depth_right > 0 | is.na(min_C_prot_depth_right))
    )
  ) %>%
  ungroup()

final_dfA<- final_dfA %>%
  filter(uniprotids %in% final_dfB$gene)

final_dfA_Crystal <- final_dfA_Crystal %>%
  filter(uniprotids %in% final_dfB_Crystal$gene)

# Save
# write.csv(final_dfA, file = "Data/final_dfA02_22_2026.csv", row.names = FALSE)
# write.csv(final_dfB, file = "Data/final_dfB02_22_2026.csv", row.names = FALSE)
# write.csv(final_dfA_Crystal, file = "Data/final_dfA_Crystal02_22_2026.csv", row.names = FALSE)
# write.csv(final_dfB_Crystal, file = "Data/final_dfB_Crystal02_22_2026.csv", row.names = FALSE)
