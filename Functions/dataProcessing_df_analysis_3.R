# Loading libraries
library(dplyr)
library(stringr)


#### Disgenet Score Percentiles ####
# Used to calculate percentiles from "raw" disgenet results
sub_finalA <- final_dfA[,c("Association_ID","score")]
sub_finalA <- unique(sub_finalA)
sub_finalA <- sub_finalA[sub_finalA$Association_ID>0,]
# Calculate score percentiles
percentiles <- quantile(sub_finalA$score, c(0.95, 0.75, 0.50))
rm(sub_finalA)

#### Create a data frame final_df_updated ####

# Absolute value of Gn and Gc
final_dfB$Gn <- abs(final_dfB$Gn)
final_dfB$Gc <- abs(final_dfB$Gc)

# Subset necessary columns
UniScoDis <- final_dfA[, c("uniprotids", "score")]

# Remove duplicates
UniScoDis <- distinct(UniScoDis)

# Keep Uniprots with the highest score
max_df <- UniScoDis %>%
  group_by(uniprotids) %>%
  filter(score == max(score)) %>%
  ungroup()

# Perform a left join to include the entanglement data
colnames(max_df)[1] <- "gene"
final_df_updated <- left_join(max_df, final_dfB, by = c("gene"))
rm(UniScoDis,max_df)

# keep only entangled proteins
final_df_updated <- final_df_updated %>%
  filter(Entanglement == "Yes")

# check whether any of Gn or Gc values are missing
stopifnot(sum(is.na(final_df_updated$Gn))==0)
stopifnot(sum(is.na(final_df_updated$Gc))==0)

# Create columns for Gsum: |Gn|+|Gc| and Gmax: Max(|Gn|, |Gc|)
final_df_updated$Gsum <- final_df_updated$Gc + final_df_updated$Gn
final_df_updated$Gmax <- pmax(final_df_updated$Gc, final_df_updated$Gn)
