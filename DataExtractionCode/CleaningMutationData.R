# Libraries
library(dplyr)

# 1. ClinVar benign (TSV)
clinvar_benign <- read.table(
  file = "Data/MutationData/clinvar_benign.tsv/clinvar_benign.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 2. ClinVar benign mapped (CSV)
clinvar_benign_mapped <- read.csv(
  file = "Data/MutationData/clinvar_benign_mapped.csv/clinvar_benign_mapped.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# 3. ClinVar pathogenic (TSV)
clinvar_pathogenic <- read.table(
  file = "Data/MutationData/clinvar_pathogenic.tsv/clinvar_pathogenic.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 4. ClinVar pathogenic mapped (CSV)
clinvar_pathogenic_mapped <- read.csv(
  file = "Data/MutationData/clinvar_pathogenic_mapped.csv/clinvar_pathogenic_mapped.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# 5. gnomAD benign (TSV)
gnomad_benign <- read.table(
  file = "Data/MutationData/gnomad_benign/gnomad_benign.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 6. gnomAD benign mapped (CSV)
gnomad_benign_mapped <- read.csv(
  file = "Data/MutationData/gnomad_benign_mapped.csv/benign_mapped.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)


############################################################################
# Mapped variant data: combined_mapped_clean
############################################################################
# list of mutation data frames
mapped_files <- list(
  ClinVarBenign = clinvar_benign_mapped,
  ClinVarPathogenic = clinvar_pathogenic_mapped,
  GnomADBenign = gnomad_benign_mapped
)

# extract column names from each data frame
mapped_colnames <- lapply(mapped_files, colnames)

# check they have the same column names before joining
all_same <- function(x) {
  first <- x[[1]]
  all(sapply(x, function(y) identical(first, y)))
}

all_mapped_same <- all_same(mapped_colnames)

if (!all_mapped_same) {
  print("Differences in column names across mapped files:")
  print(mapped_colnames)
}

# join the mutation data frames
combined_mapped <- bind_rows(
  ClinVarBenign = clinvar_benign_mapped,
  ClinVarPathogenic = clinvar_pathogenic_mapped,
  GnomADBenign = gnomad_benign_mapped,
  .id = "source"
)

# remove empty columns
remove_single_value_columns <- function(df) {
  df[, sapply(df, function(col) length(unique(col)) > 1)]
}

combined_mapped_clean <- remove_single_value_columns(combined_mapped)

# reorder column to have uniprotid as the first column
combined_mapped_clean <- combined_mapped_clean[, c("Uniprot_canonical_isoform_.non_canonical.",
                                                   setdiff(names(combined_mapped_clean), "Uniprot_canonical_isoform_.non_canonical."))]
colnames(combined_mapped_clean)[1] <- "uniprotids"

# Mutation column
combined_mapped_clean$Mutation <- "Yes"

# Pathogenic mutation column
combined_mapped_clean$Pathogenic <- ifelse(combined_mapped_clean$source == "ClinVarPathogenic",
                                           "Yes",
                                           ifelse(combined_mapped_clean$source %in% c("GnomADBenign", "ClinVarBenign"),
                                                  "No",
                                                  NA))

# Filter to only have missense mutations
combined_mapped_clean <- combined_mapped_clean %>% filter(grepl("missense", Consequences, ignore.case = TRUE))

write.csv(combined_mapped_clean,
          file = "Data/combined_mapped_clean_missense.csv",
          row.names = FALSE)
