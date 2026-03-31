Data Extraction Code
================

# Overview

This folder contains scripts used to download, clean, harmonize, and
assemble the main datasets used across the project:

- Essentiality mapping to UniProt

- Missense mutation datasets (ClinVar / gnomAD), cleaned and harmonized

- Preprocessed mutation datasets used in residue-level and protein-level
  analyses

- Modeling-ready feature tables used for delta / RF workflows

These scripts produce analysis-ready output files in `Data/`, which are
then used by downstream analysis scripts.

# Data Extraction Codes

## 1) EssentialityData.R

### What it does

This script creates a UniProt-level essentiality label by combining
DepMap CRISPR essential genes with the set of genes tested in the CRISPR
screens, parsing gene identifiers, mapping those identifiers to reviewed
UniProt entries, and then assigning each UniProt entry an
essential/non-essential status based on whether it appears in the
inferred common essential list, collapsing across multiple mappings when
necessary so each UniProt entry has a single, consistent essentiality
classification for use as a modeling covariate.

### Primary outputs

`Data/essentiality.csv`

## 2) CleaningMutationData.R

### What it does

This script consolidates multiple mutation sources into a single
harmonized mutation dataset by loading ClinVar and gnomAD variant
tables, aligning the mapped outputs to a consistent schema, adding
standardized mutation and pathogenicity indicators, removing
uninformative columns, and filtering to missense variants only,
resulting in a clean, unified mutation table that can be used
consistently across residue-level and protein-level mutation analyses.

### Primary outputs

`Data/combined_mapped_clean_missense.csv`

## 3) m1_data.R (Mutation dataset for residue-level analysis)

### What it does

This script prepares the mutation data for residue-level analyses by
joining cleaned missense variants to the protein-level
disease/entanglement dataset, removing contradictory duplicate mutation
records, filtering proteins and variants according to analysis criteria
(including quality thresholds and minimum mutation burden), expanding
each protein into a full residue-level grid, merging mutations onto
residue positions, and labeling whether each residue falls inside an
entangled region, producing a residue-resolved dataset focused on
entangled proteins for downstream residue-level analyses.

### Primary outputs

- `Data/res_merged_dt_af_<plddt_thresh>_<Sys.Date()>.rds`

- `Data/res_merged_dt_af_<plddt_thresh>_<Sys.Date()>.csv`

## 4) m2_data.R (Mutation dataset for protein-level analysis)

### What it does

This script prepares the mutation data for protein-level analyses by
joining cleaned missense variants to the protein-level
disease/entanglement dataset, resolving contradictory duplicate mutation
records, applying quality and validity filters (including restricting
mutations to canonical protein length and requiring a minimum number of
missense variants per protein), and aggregating to one row per protein
with mutation burden summaries and pathogenic proportions, producing the
protein-level mutation summary dataset used in downstream
binomial/proportion-based modeling.

### Primary outputs

`Data/protein_summary_af_<plddt_thresh>_<Sys.Date()>.csv`

## 5) m3_data.R (Mutation dataset for protein-level within classes)

### What it does

This script prepares class-stratified mutation analysis inputs by
deriving disease-class-specific disease indicators from DisGeNET scores
(using percentile-based thresholds), and then joining cleaned missense
variants while applying the same duplicate-resolution and mutation
burden filters used previously, resulting in a dataset setup that
supports downstream protein-level mutation modeling within disease
classes.

### Primary outputs

`final_df`

## 6) RF_Preprocessing.R (Preprocessing for delta analysis)

### What it does

This script constructs the modeling-ready feature table used in the
delta/RF workflow by collapsing entanglement-related metrics to a single
row per protein, ensuring consistent handling of missing values, and
integrating additional covariates such as mutation burden, structural
exposure summaries, essentiality labels, and a disease outcome derived
from DisGeNET scores, producing a harmonized dataset suitable for
downstream predictive modeling and delta-score computation.

### Primary outputs

`Data/final_dfB_max_af_70_raw_l_02_20.csv`

## 7) DisgenetExtraction.R

### What it does

This script downloads and constructs the core disease–protein
association dataset by mapping UniProt IDs to Entrez IDs, querying
DisGeNET for gene–disease associations, and assembling the main
analysis-ready data frames used throughout the project (including
disease annotations, scores, and protein-level mappings). It also
integrates entanglement, structural, and annotation information to
produce the final datasets used in downstream analyses.

### Important note on DisGeNET access

This script uses the disgenet2r package to query DisGeNET. However, this
package is no longer actively maintained and may not function reliably.

To reproduce the data:

- You will need to create a DisGeNET account and obtain API access.
- If possible, request access to DisGeNET version 24.1, which was used
  for this project.
- Otherwise, you may need to use the latest available version, noting
  that results may differ due to updates in the database.

### Primary outputs

- final_dfA\*.csv
- final_dfB\*.csv
- Crystal equivalents of the above datasets
