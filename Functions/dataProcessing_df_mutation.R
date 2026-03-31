# Entanglement Data
af_ent_raw <- read.csv(
  file = "Data/Human_AF_combined_20250521 (1).csv",
  header = TRUE,
  sep = "|",
  stringsAsFactors = FALSE
)

af_ent_raw <- af_ent_raw %>%
  rename(
    Gn = gn,
    Gc = gc
  )

af_ent_raw <- af_ent_raw %>%
  mutate(CCBond = case_when(
    CCBond %in% c("True","1.0", "1") ~ "True",
    CCBond %in% c("False","0.0", "0") ~ "False",
    TRUE ~ NA
  ))

af_ent_raw <- af_ent_raw %>%
  mutate(crossings = str_replace(crossings, "\\.0$", ""))

af_ent_raw <- distinct(af_ent_raw)

af_ent_dat <- af_ent_raw %>%
  rename(EntRegion = ent_region)

af_ent_dat <- af_ent_dat %>%
  filter(CCBond != "True",             # Keep rows where CCBond is not TRUE
         !is.na(EntRegion),          # Exclude rows with NA in EntRegion
         trimws(EntRegion) != "")    # Exclude rows with empty (or whitespace) EntRegion


af_ent_dat <- af_ent_dat[, c("gene", "EntRegion")]


af_ent_dat_clean <- af_ent_dat %>%
  # Rename the first column "gene" to "uniprotids"
  rename(uniprotids = gene) %>%
  # Create a new column "numeric_EntRegion" that converts the comma separated 
  # EntRegion string into a numeric vector.
  mutate(numeric_EntRegion = ifelse(
    str_trim(EntRegion) == "", 
    list(NA),  # If empty, store as NA (within a list)
    # Otherwise, split the text by comma and convert to numeric.
    lapply(EntRegion, function(x) as.numeric(str_trim(unlist(str_split(x, ",")))))
  )) %>%
  # Group by the uniprotids and merge the numeric lists:
  group_by(uniprotids) %>%
  summarise(
    numeric_EntRegion = list(na.omit(unlist(numeric_EntRegion))),
    n_EntRegion = length(numeric_EntRegion[[1]])
  )


rm(af_ent_raw, af_ent_dat)


