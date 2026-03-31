cat('date:',format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
cat("type: af", "\n")
cat("plddt_thresh: 70", "\n")

# libraries
library(grf)
library(dplyr)
library(randomForest)
library(dplyr)
library(tidyr)
library(iml)
library(ggplot2)
library(gridExtra)

# loading data
# ent data
final_dfB_max<- read.csv("Data/final_dfB_max_af_70_raw_l_02_20.csv")

# transmemebrane proteins
membrane <- read.delim("Data/uniprotkb_organism_id_9606_AND_keyword_2025_11_16.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

# Role columns for U (post-treatment, entanglement-related) and Z (other covariates)
U_cols <- c("Gn", "Gc", "Gmax", "Gsum",
            "N_term_thread", "C_term_thread",
            "ENT.ID",
            "num_zipper_nc", "num_loop_contacting_res", "num_cross_nearest_neighbors",
            "min_N_prot_depth_left", "min_C_prot_depth_right",
            "Travatos_G")

Z_cols <- c("Length","Total_Mutations_n","mean_SASA","Essential")

# making categorical variables in to factors
df <- final_dfB_max
df$disease <- factor(df$disease)
df$Entanglement <- factor(df$Entanglement)
df$Essential <- factor(df$Essential)

# separating the data frames btw entangled and nonentangled proteins
df1 <- df %>% filter(Entanglement=="Yes")
df0 <- df %>% filter(Entanglement=="No")

# check
# all.equal(nrow(df0)+nrow(df1),nrow(df))

# predictors for entangled and nonentangled proteins
m1_predictors <- c(Z_cols, U_cols) # entangled predictors
m0_predictors <- Z_cols # nonentangled predictors

# datasets for entanglement == 1 and 0
df1_fit <- df1 %>% dplyr::select(disease, tidyselect::all_of(m1_predictors)) # entangled
df0_fit <- df0 %>% dplyr::select(disease, tidyselect::all_of(m0_predictors)) # nonentangled

# check missing values
all.equal(df1_fit,na.omit(df1_fit))
all.equal(df0_fit,na.omit(df0_fit))

# sapply(df1_fit,class)
# sapply(df0_fit,class)

# fit two random forests
set.seed(1)
m1 <- randomForest(disease ~ ., data = df1_fit,
                   ntree = 1000,
                   na.action = na.omit,
                   importance=TRUE) # Computes permutation (how much accuracy drops when
                                    # each variable is randomly permuted) and Gini feature importance
                                    # (measures how much each variable reduces node impurity across all trees).

m0 <- randomForest(disease ~ ., data = df0_fit,
                   ntree = 1000,
                   na.action = na.omit,
                   importance=TRUE)


# for entangled proteins, evaluate P(D=1|E=1,(U,Z)) - P(D=1|E=0,(Z))
ent_df <- df[df$Entanglement=="Yes",] # entangled proteins only

# estimate probabilities
phat_1<-predict(m1, newdata = ent_df, type="prob") # predicted disease probabilities from the model trained on entangled proteins
phat_0<-predict(m0, newdata = ent_df, type="prob") # predicted disease probabilities from the model trained on nonentangled proteins

ent_df$p1 <- phat_1[,2] # probability of disease if entangled
ent_df$p0 <- phat_0[,2] # probability of disease if nonentangled
ent_df$Delta_2RF <- phat_1[,2]-phat_0[,2] # estimated increase (or decrease) in disease probability that is attributable
                                          # to being entangled, after accounting for its Length, mutations,
                                          # mean SASA, and Essential status.

write.csv(ent_df, file = paste0("Results/Dataframes/delta_ent_", format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

####################################################################################################

# remove transmembrame proteins for simulations
ent_df <- ent_df %>%
  filter(!gene %in% membrane$Entry) # remove transmembrane proteins

# ordered ent+disease group
ent_disease_df <- filter(ent_df, disease =="Yes") %>% # only entangled and disease associated proteins
  arrange(desc(Delta_2RF)) # largest to smallest delta

n_treated <- nrow(ent_disease_df) # number of disease associated entangled proteins
treated_gene_ordered <- ent_disease_df$gene # disease associated entangled proteins from largest to smallest delta

ent_disease_df <- ent_disease_df[1:50,] # only top 50 delta
top50_gene_ordered <- ent_disease_df$gene

# scale "mean_SASA","Length","Total_Mutations_n" variables for matching
df <- df %>% mutate(s_SASA = scale(mean_SASA),
                    s_Length = scale(Length),
                    s_nmut = scale(Total_Mutations_n))

df <- df %>%
  filter(!gene %in% membrane$Entry) # remove transmembrane proteins

# top ent & disease group  + no ent & no disease group
match_df1 <- df[df$gene %in% treated_gene_ordered,c("gene","s_SASA","s_Length","s_nmut","Essential")] %>% mutate(group=1)
match_df1<- match_df1[match(treated_gene_ordered,match_df1$gene),] # ordered based on delta_2rf

match_df2 <- df[df$disease=="No"&df$Entanglement=="No",c("gene","s_SASA","s_Length","s_nmut","Essential")] %>% mutate(group=0)
match_df = rbind(match_df1,match_df2) # join

rownames(match_df) = NULL

# top50 ent & disease group  + no ent & no disease group
match_df1_50 <- df[df$gene %in% top50_gene_ordered,c("gene","s_SASA","s_Length","s_nmut","Essential")] %>% mutate(group=1)
match_df1_50<- match_df1_50[match(top50_gene_ordered,match_df1_50$gene),] # ordered based on delta_2rf

match_df2_50 <- df[df$disease=="No"&df$Entanglement=="No",c("gene","s_SASA","s_Length","s_nmut","Essential")] %>% mutate(group=0)
match_df_50 = rbind(match_df1_50,match_df2_50)

rownames(match_df_50) = NULL

## 1) Match with ratio = 1 for all and top 50
library(MatchIt)
set.seed(123)
m <- matchit(
  group ~ s_SASA + s_Length + s_nmut,
  data     = match_df,
  method   = "nearest",
  distance = "euclidean",
  ratio    = 1,
  replace  = FALSE,
  exact    = ~ Essential
)

idx_control = as.numeric(m$match.matrix[,1]) # selected index for "control" proteins (no disease nonetangled)
ctrl_gene <- match_df[idx_control,"gene"] # selected "control" proteins (no disease nonetangled)

m_50 <- matchit(
  group ~ s_SASA + s_Length + s_nmut,
  data     = match_df_50,
  method   = "nearest",
  distance = "euclidean",
  ratio    = 1,
  replace  = FALSE,
  exact    = ~ Essential
)

idx_control_50 = as.numeric(m_50$match.matrix[,1]) # selected "control" genes
c50_gene <- match_df_50[idx_control_50,"gene"]


#### outputs ####
# out1 and out1_50: gene, delta_2rf, dist, scaled variables (used for matching)
# out2 and out2_50: gene, delta_2rf, dist, variables (used for matching)
# out3 and out3_50: subset of df with 100 top genes + 100 matched genes

# check exact matching
all.equal(match_df[1:n_treated,"Essential"], match_df[idx_control,"Essential"]) # NA because not enough proteins to match all
all.equal(match_df_50[1:50,"Essential"], match_df_50[idx_control_50,"Essential"])

# distance
euclidean_dist <- apply(
  match_df[1:n_treated,2:4] - match_df[idx_control,2:4],
  1,
  function(x) sqrt(sum(x^2))
)

euclidean_dist_50 <- apply(
  match_df_50[1:50,2:4] - match_df_50[idx_control_50,2:4],
  1,
  function(x) sqrt(sum(x^2)))

# combine
out1 <- data.frame(cbind(match_df[1:n_treated,], match_df[idx_control,], euclidean_dist))
colnames(out1) = gsub(colnames(out1),pattern=".1",replacement="_ctrl")
out1 %>% head()

out1_50 <- data.frame(cbind(match_df_50[1:50,], match_df_50[idx_control_50,], euclidean_dist_50))
colnames(out1_50) = gsub(colnames(out1_50),pattern=".1",replacement="_ctrl")
out1_50 %>% head()

# summary table 2 based on raw covariates
out2.1 = ent_df[ent_df$gene %in% treated_gene_ordered,
                c("gene","Delta_2RF","mean_SASA","Length","Total_Mutations_n","Essential")]
out2.1 = out2.1[match(treated_gene_ordered,out2.1$gene),]

out2.2 = df[df$gene %in% ctrl_gene,
            c("gene","mean_SASA","Length","Total_Mutations_n","Essential")]
out2.2 = out2.2[match(ctrl_gene,out2.2$gene),]

all.equal(out2.1$gene, treated_gene_ordered)
all.equal(out2.2$gene, ctrl_gene)

out2 = data.frame(cbind(out2.1, out2.2, euclidean_dist))
colnames(out2) = gsub(colnames(out2),pattern=".1",replacement="_ctrl")
out2 %>% head()

all.equal(out1$gene,out2$gene)
all.equal(out1$gene_ctrl, out2$gene_ctrl)
out1$Delta_2RF<- out2$Delta_2RF

out2.1_50 = ent_df[ent_df$gene %in% top50_gene_ordered,
                   c("gene","Delta_2RF","mean_SASA","Length","Total_Mutations_n","Essential")]
out2.1_50 = out2.1_50[match(top50_gene_ordered,out2.1_50$gene),]

out2.2_50 = df[df$gene %in% c50_gene,
               c("gene","mean_SASA","Length","Total_Mutations_n","Essential")]
out2.2_50 = out2.2_50[match(c50_gene,out2.2_50$gene),]

all.equal(out2.1_50$gene, top50_gene_ordered)
all.equal(out2.2_50$gene, c50_gene)

out2_50 = data.frame(cbind(out2.1_50, out2.2_50, euclidean_dist_50))
colnames(out2_50) = gsub(colnames(out2_50),pattern=".1",replacement="_ctrl")
out2_50 %>% head()

all.equal(out1_50$gene,out2_50$gene)
all.equal(out1_50$gene_ctrl, out2_50$gene_ctrl)
out1_50$Delta_2RF<- out2_50$Delta_2RF

# summary table 3 (all covariates)
out3_tmp <- df[df$gene %in% c(treated_gene_ordered, ctrl_gene),] %>%
  select(-starts_with("s_")) # remove scaled variables
out3_tmp <- out3_tmp[match(c(treated_gene_ordered, ctrl_gene), out3_tmp$gene),]

out3 <- data.frame(
  cbind(
    out3_tmp[1:n_treated,],
    out3_tmp[(n_treated + 1):(2*n_treated),]
  )
)
colnames(out3) = gsub(colnames(out3),pattern=".1",replacement="_ctrl")

all.equal(out3$gene, treated_gene_ordered)
all.equal(out3$gene_ctrl, ctrl_gene)

out3$Delta_2RF<-out2$Delta_2RF
out3$euclidean_dist = out2$euclidean_dist

out3_50 <- df[df$gene%in%c(top50_gene_ordered,c50_gene),]%>%
  select(-starts_with("s_"))
out3_50 <- out3_50[match(c(top50_gene_ordered,c50_gene), out3_50$gene),]

out3_50 <- data.frame(
  cbind(
    out3_50[1:50,],
    out3_50[51:100,]
    )
  )
colnames(out3_50) = gsub(colnames(out3_50),pattern=".1",replacement="_ctrl")

all.equal(out3_50$gene,top50_gene_ordered)
all.equal(out3_50$gene_ctrl, c50_gene)

out3_50$Delta_2RF<-out2_50$Delta_2RF
out3_50$euclidean_dist_50 = out2_50$euclidean_dist_50

# filter to only rows with a matched control
out1 <- out1[!is.na(out1$gene_ctrl), ]
out2 <- out2[!is.na(out2$gene_ctrl), ]
out3 <- out3[!is.na(out3$gene_ctrl), ]

# reorder some columns
out1<- out1 %>% relocate(Delta_2RF, .after = "gene") %>% relocate(euclidean_dist, .after = "Delta_2RF")
out2<- out2 %>%  relocate(euclidean_dist, .after = "Delta_2RF")
out3<- out3 %>% relocate(Delta_2RF, .after = "gene") %>% relocate(euclidean_dist, .after = "Delta_2RF")

out1_50<- out1_50 %>% relocate(Delta_2RF, .after = "gene") %>% relocate(euclidean_dist_50, .after = "Delta_2RF")
out2_50<- out2_50 %>%  relocate(euclidean_dist_50, .after = "Delta_2RF")
out3_50<- out3_50 %>% relocate(Delta_2RF, .after = "gene") %>% relocate(euclidean_dist_50, .after = "Delta_2RF")

write.csv(out1,     paste0("Results/Dataframes/geneall_scaledcov_",    format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
write.csv(out2,     paste0("Results/Dataframes/geneall_unscaledcov_",  format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
write.csv(out3,     paste0("Results/Dataframes/geneall_allvar_",       format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

write.csv(out1_50,  paste0("Results/Dataframes/gene50_scaledcov_",     format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
write.csv(out2_50,  paste0("Results/Dataframes/gene50_unscaledcov_",   format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
write.csv(out3_50,  paste0("Results/Dataframes/gene50_allvar_",        format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)

# dataframes for entangled and nonentangled proteins with their predictors
X1 <- df1 %>% dplyr::select(gene, tidyselect::all_of(m1_predictors))
X0 <- df0 %>% dplyr::select(gene, tidyselect::all_of(m0_predictors))

mod1 <- Predictor$new(
  model = m1,
  data=X1,
  type = "prob",
  class="Yes"
)
mod0 <- Predictor$new(
  model = m0,
  data=X0,
  type = "prob",
  class="Yes"
)

# Compute Shapley values for observations
set.seed(123)
g50 <- out3_50[1:50,"gene"]
i=2
x1 <- X1[X1$gene==g50[i], , drop = FALSE]
shap <- Shapley$new(mod1, x.interest = x1) # for this protein, how much did each feature contribute to the predicted disease probability?

set.seed(123)
x0 <- X1[X1$gene==g50[i], c("gene", m0_predictors), drop = FALSE]
shap0 <- Shapley$new(mod0, x.interest = x0)


pdf(paste0("Results/Plots/shap_top50_pairs", format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".pdf"), width = 14, height = 8)

# Initialize storage dataframe
shap_results <- data.frame()

for (i in 1:50) {
  set.seed(123)

  # ---------- p1 ----------
  x1 <- X1[X1$gene==g50[i], , drop = FALSE]
  shap  <- Shapley$new(mod1, x.interest = x1)
  actual1  <- round(shap$y.hat.interest, 3)
  avgpred1 <- round(shap$y.hat.average, 3)

  # Extract p1 shapley values: wide format
  shap_p1 <- shap$results %>%
    select(feature, phi) %>%
    pivot_wider(names_from = feature, values_from = phi, names_prefix = "p1_")

  # Plot p1
  p1 <- plot(shap) +
    ggtitle("p1",
            subtitle = paste0("Actual prediction: ", actual1, "\nAverage prediction: ", avgpred1)) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 6))


  # ---------- p0 ----------
  set.seed(123)
  x0 <- X1[X1$gene==g50[i], c("gene", m0_predictors), drop = FALSE]
  shap0 <- Shapley$new(mod0, x.interest = x0)
  actual2  <- round(shap0$y.hat.interest, 3)
  avgpred2 <- round(shap0$y.hat.average, 3)

  # Extract p0 shapley values: wide format
  shap_p0 <- shap0$results %>%
    select(feature, phi) %>%
    pivot_wider(names_from = feature, values_from = phi, names_prefix = "p0_")

  # Plot p0
  p2 <- plot(shap0) +
    ggtitle("p0",
            subtitle = paste0("Actual prediction: ", actual2, "\nAverage prediction: ", avgpred2)) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10),
          axis.text.y = element_text(size = 6))


  # ---------- Combine into shap_results ----------
  shap_results <- bind_rows(
    shap_results,
    cbind(
      data.frame(
        gene = g50[i],
        actual_p1 = actual1,
        avgpred_p1 = avgpred1,
        actual_p0 = actual2,
        avgpred_p0 = avgpred2
      ),
      shap_p1,
      shap_p0
    )
  )


  # Display plots
  grid.arrange(p1, p2, ncol = 2)
}

dev.off()

# identify p1 and p0 shapley columns
p1_cols <- grep("^p1_", names(shap_results), value = TRUE)
p0_cols <- grep("^p0_", names(shap_results), value = TRUE)

shap_results <- shap_results %>%
  rowwise() %>%
  mutate(
    top_p1_contributor = sub("^p1_", "", p1_cols[which.max(abs(c_across(all_of(p1_cols))))]),
    top_p0_contributor = sub("^p0_", "", p0_cols[which.max(abs(c_across(all_of(p0_cols))))])
  ) %>%
  ungroup()

table(shap_results$top_p1_contributor)
table(shap_results$top_p0_contributor)

rownames(shap_results) = NULL

write.csv(shap_results, paste0("Results/Dataframes/shap_results", format(Sys.time(), "_%Y-%m-%d_%H-%M-%S"), ".csv"), row.names = FALSE)
