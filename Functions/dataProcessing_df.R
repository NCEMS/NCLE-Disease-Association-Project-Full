# Loading libraries
library(dplyr)
library(tidyr)

dataProcessing_df<- function(type, plddt_thresh){
  # load corresponding Data according to type
  # merge with plddt and keep observations with plddt >= plddt_thresh
  type = match.arg(type, c('af', 'crystal', 'crystal_yes_af_2nd', 'af_not_crystal', 'af_crystal'))

  # 1. All AlphaFold
  af_A <- read.csv("Data/final_dfA02_22_2026.csv")
  af_B <- read.csv("Data/final_dfB02_22_2026.csv")

  # 2. All Crystal
  crystal_A <- read.csv("Data/final_dfA_Crystal02_22_2026.csv")
  crystal_B <- read.csv("Data/final_dfB_Crystal02_22_2026.csv")

  #3. Data set when Crystal available, if not AF
  crystal_yes_af_2nd_A <- bind_rows(crystal_A, af_A %>% filter(!(uniprotids %in% crystal_A$uniprotids)))
  crystal_yes_af_2nd_B <- bind_rows(crystal_B %>% dplyr::select(-dplyr::any_of("CCBond")), af_B %>% filter(!(gene %in% crystal_B$gene)))

  # 4. Data set with Crystal uniprots but AF Data
  af_crystal_A <- af_A %>% filter(uniprotids %in% crystal_A$uniprotids)
  af_crystal_B <- af_B %>% filter(gene %in% crystal_B$gene)

  # 5. AF not in Crystal
  af_not_crystal_A <- anti_join(af_A, crystal_A, by = "uniprotids")
  af_not_crystal_B <- anti_join(af_B, crystal_B, by = "gene")

  # Assign `final_dfA` and `final_dfB` based on selected input type
  if (type == "af") {
    final_dfA <- af_A
    final_dfB <- af_B
  } else if (type == "crystal") {
    final_dfA <- crystal_A
    final_dfB <- crystal_B
  } else if (type == "crystal_yes_af_2nd") {
    final_dfA <- crystal_yes_af_2nd_A
    final_dfB <- crystal_yes_af_2nd_B
  } else if (type == "af_crystal") {
    final_dfA <- af_crystal_A
    final_dfB <- af_crystal_B
  } else if (type == "af_not_crystal") {
    final_dfA <- af_not_crystal_A
    final_dfB <- af_not_crystal_B
  } else {
    stop("Invalid input_type. Choose from 'af', 'crystal', 'crystal_yes_af_2nd', or 'af_crystal'.")
  }


  # Remove extra space
  final_dfA$uniprotids <- trimws(final_dfA$uniprotids)
  final_dfB$gene <- trimws(final_dfB$gene)


  # Load plddt Data
  plddt <- read.csv("Data/Avg_pLDDTs.csv", header = TRUE, sep = ",")
  plddt <- plddt %>% rename(uniprotids = gene)
  plddt <- plddt %>% rename(plddt = X.pLDDT.)

  # Join plddt Data
  final_dfA <- final_dfA %>%
    left_join(plddt, by = "uniprotids")
  final_dfB <- final_dfB %>%
    left_join(plddt, by = c("gene" = "uniprotids"))

  # Select Thresholds
  if (type %in% c("af", "af_crystal")) {
    final_dfA <- final_dfA %>%
      filter(plddt >= plddt_thresh)

    final_dfB <- final_dfB %>%
      filter(plddt >= plddt_thresh)
  }

  return(list(final_dfA = final_dfA, final_dfB= final_dfB))

}



