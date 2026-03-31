Data
================

# Data Folder

This folder contains all raw, intermediate, and finalized datasets used
in the Disease Entanglement project. Each file is described below along
with the scripts where it is used.

# Data availability notice

Due to licensing restrictions associated with DisGeNET, this repository
does not include certain raw data required to fully reproduce all
analyses. Specifically, raw DisGeNET files with fields such as disease
names, scores, and identifiers have been removed to comply with the
DisGeNET user agreement.

As a result, some scripts in this repository will not run as-is because
of these missing data components.

This repository is intended to:

- Provide the original analysis code and workflow
- Document how the data were processed and analyzed
- Enable reproducibility for users who independently obtain access to
  the required DisGeNET data

Users who wish to fully reproduce the analyses must obtain the
appropriate data directly from DisGeNET and place it into the expected
directories as described in the repository.

## Data access (CyVerse)

The datasets required to run these analyses are not included in this
repository due to their large size and data-sharing restrictions. Only
the allowed subset of data is provided externally on CyVerse, as certain
components cannot be distributed in accordance with the DisGeNET user
agreement. As a result, full reproduction of the results is not possible
using the provided data alone, and access to DisGeNET is required to
obtain the complete dataset used in this study.

1.  Create an account

Create a free CyVerse account:
[CyVerse](https://user.cyverse.org/signup)

2.  Access and download the data

Once your account is created and verified, you will be able to access
the project data through the CyVerse Data Store at the following
location:

[Download Data_Cyverse folder]()

3.  Download the data

You can download the data in two ways:

**Option A: Direct download (via browser)**  
Download the full `Data_Cyverse` folder from the CyVerse interface.

**Option B: Using CyVerse GO commands (recommended for large files)**  
Instructions for installing and using GO commands:  
[CyVers Go
commands](https://learning.cyverse.org/ds/gocommands/installation/)

4.  Place data in the repository

After downloading:

- Locate the folder named `Data_Cyverse`
- Ensure the `Data/` folder in this repository is empty
- Move all contents from `Data_Cyverse` into the `Data/` folder

The directory structure inside `Data/` must match exactly what is
provided in `Data_Cyverse` before running any analyses.

### Human_AF_combined_20250521 (1).csv

**Used in:** `Functions/dataProcessing_df_mutation.R`  
Raw AlphaFold structural entanglement output (csv). Has entanglement
metrics for human proteins, including G metrics and entangled region
annotations.

| Column Name | Definition |
|----|----|
| gene | UniProt identifier for the protein |
| PDB | Structure identifier (AlphaFold model ID) |
| chain | Chain identifier within the structure |
| ENT.ID | Sequential identifier for each non-covalent lasso entanglement (NCLE); the maximum value per protein equals the number of entanglements (E<sub>NCLE</sub>) |
| Gn | Gaussian linking value for the N-terminal segment of the entanglement (G<sub>n</sub>) |
| N_term_thread | Number of N-terminal thread crossings of the loop plane (N<sub>threads</sub>) |
| Gc | Gaussian linking value for the C-terminal segment of the entanglement (G<sub>c</sub>) |
| C_term_thread | Number of C-terminal thread crossings of the loop plane (C<sub>threads</sub>) |
| i | Start residue index of the loop region |
| j | End residue index of the loop region |
| NC | Residue pair (i, j) forming the non-covalent contact that closes the loop |
| NC_wbuff | Non-covalent contact including ± buffer residues |
| NC_region | Region classification of the non-covalent contact |
| crossings | Residues where the threading segment crosses the loop plane |
| crossings_wbuff | Crossings including ± buffer residues |
| crossings_region | Region classification of crossings |
| ent_region | Region of the protein defined as entangled based on loop-closing and threading residues |
| loopsize | Number of residues in the representative NCLE loop |
| num_zipper_nc | Number of loop-closing native contacts forming the entanglement cluster (*N*<sub>zipper</sub>) |
| perc_bb_loop | Percentage of the primary sequence encompassed by the loop |
| num_loop_contacting_res | Total number of contacts made by loop residues with the rest of the protein (contacts defined within 8Å and ≥3 residues apart) (N<sub>loop-cont</sub>) |
| num_cross_nearest_neighbors | Total number of contacts made by crossing residues and their ±3 residue buffer with the rest of the protein (N<sub>cross-cont</sub>) |
| ent_coverage | Fraction of residues in the entangled region, defined as residues within 8Å of key loop-closing and threading residues |
| min_N_prot_depth_left | Minimum depth of N-terminal crossings from the N-terminus, normalized by protein length (*d*<sub>P</sub>(N)) |
| min_N_thread_depth_left | Minimum depth of N-terminal crossings from the N-terminus, normalized by thread length (*d*<sub>T</sub>(N)) |
| min_N_thread_slippage_left | Minimum depth of N-terminal crossings from the N-terminus in residues (*d*<sub>S</sub>(N)) |
| min_C_prot_depth_right | Minimum depth of C-terminal crossings from the C-terminus, normalized by protein length (*d*<sub>P</sub>(C)) |
| min_C_thread_depth_right | Minimum depth of C-terminal crossings from the C-terminus, normalized by thread length (*d*<sub>T</sub>(C)) |
| min_C_thread_slippage_right | Minimum depth of C-terminal crossings from the C-terminus in residues (*d*<sub>S</sub>(C)) |
| prot_size | Total number of residues in the protein |
| ACO | Absolute contact order |
| RCO | Relative contact order: ACO normalized by protein length |
| CCBond | Indicator of disulfide (cysteine–cysteine) bonds |
| Knot_type | Classification of knot topology in the protein |
| pLDDT_knotcore | Average AlphaFold pLDDT confidence score within the knot core region |

The entanglement variables correspond to entanglement metrics defined in
the Supplementary Material and described in detail in the associated
publication: Natively Entangled Proteins Are Linked to Human Disease and
Pathogenic Mutations Likely Due to Misfolding

### Human_AF_combined_20250521 (1).txt

**Used in:** `DataExtractionCode/DisgenetExtraction.R`  
Raw AlphaFold structural entanglement output (txt) used during disease
association integration and preprocessing.

This file contains the same variables and definitions as
`Human_AF_combined_20250521 (1).csv`, but is provided in text format for
compatibility with different preprocessing workflows.

### Human_EXP_combined_20250614.txt

**Used in:** `DataExtractionCode/DisgenetExtraction.R` Combined
experimental (crystal) structure entanglement dataset used to generate
crystal-based disease association datasets.

This file contains the same variables and definitions as
`Human_AF_combined_20250521 (1).csv`, but derived from experimentally
determined (crystal) protein structures instead of AlphaFold
predictions.

### Human_AF_SASA.csv

**Used in:** `DataExtractionCode/RF_Preprocessign.R`  
Contains solvent-accessible surface area (SASA) metrics for AlphaFold
structures. Used as input features for random forest modeling.

| Column Name | Definition |
|----|----|
| uniprotids | UniProt identifier for the protein |
| resid | Residue index (position in the protein sequence) |
| resname | Amino acid identity at the residue position (three-letter code) |
| SASA | Solvent-accessible surface area of the residue, representing the extent of exposure to solvent |

### Trovato.xlsx

**Used in:** `DataExtractionCode/DisgenetExtraction.R`  
Contains max \|G<sub>Intra</sub>\| entanglement metric values that are
merged with structural datasets.

| Column Name | Definition |
|----|----|
| gene | UniProt identifier for the protein |
| PDB | Structure identifier (AlphaFold model) |
| chain | Chain identifier within the structure |
| ENT-ID | Sequential identifier for each non-covalent lasso entanglement (NCLE); see above |
| maxGLN | Maximum absolute intrachain Gaussian linking value (max |
| loopsize | Number of residues in the representative entanglement loop |
| threadsize | Number of residues in the threading segment passing through the loop |

### human_Xtal_maxGLN_Trovato.csv

**Used in:** `DataExtractionCode/DisgenetExtraction.R`  
Crystal structure dataset containing
"max"*"\|"*italic(G)\[Intra\]\*"\| entanglement metric values.

This dataset contains the same entanglement metric (maxGLN) as
Trovato.xlsx, but computed from crystal structures instead of AlphaFold
predictions.

### final_resultsRAW05_21_2024.csv

**Used in:** `DataExtractionCode/DisgenetExtraction.R`  
Raw disease association results prior to filtering and merging with
structural data.

| Column Name | Definition |
|----|----|
| gene_symbol | Official gene symbol |
| geneid | NCBI Entrez Gene identifier |
| ensemblid | Ensembl gene identifier |
| geneNcbiType | Gene type classification from NCBI |
| geneDSI | Disease Specificity Index (DSI), measuring how specific a gene is to particular diseases |
| geneDPI | Disease Pleiotropy Index (DPI), measuring how broadly a gene is associated with multiple diseases |
| uniprotids | UniProt identifier(s) mapped to the gene |
| protein_classid | Identifier for the protein functional class |
| protein_class_name | Name of the protein functional class (e.g., enzyme, transcription factor) |
| disease_name | Name of the associated disease |
| diseaseType | Type/category of disease |
| diseaseUMLSCUI | Unique identifier for the disease in the UMLS system |
| diseaseClasses_MSH | Disease classification based on MeSH (Medical Subject Headings) |
| diseaseClasses_UMLS_ST | Disease classification based on UMLS Semantic Types |
| diseaseClasses_DO | Disease classification based on Disease Ontology (DO) |
| diseaseClasses_HPO | Disease classification based on Human Phenotype Ontology (HPO) |
| score | DisGeNET gene–disease association score representing the strength of evidence for each gene–disease pair |
| yearInitial | Year the gene–disease association was first reported |
| yearFinal | Most recent year the association was reported |
| diseaseid | Internal DisGeNET disease identifier |
| DisGeNET | Binary indicator (“Yes”/“No”) of whether the gene was present in DisGeNET at the time of query |

Note: This file is not included in this repository due to the data
restrictions described above related to DisGeNET’s user agreement.

### essentiality.csv

**Used in:**  
- `DataExtractionCode/EssentialityData.R` -
`DataExtractionCode/DisgenetExtraction.R` -
`DataExtractionCode/RF_Preprocessign.R`  
Processed gene essentiality annotations mapped to UniProt IDs.

| Column Name | Definition |
|----|----|
| Entry | UniProt identifier for the protein |
| Essential | Gene essentiality classification: “Yes” (essential), “No” (non-essential), or “NT” (not tested) |

### CRISPRGeneEffect.csv

**Used in:** `DataExtractionCode/EssentialityData.R` Raw CRISPR gene
effect dataset where rows represent cell lines/samples and columns
represent genes. Each value corresponds to a gene effect score, which
quantifies the impact of gene knockout on cell viability (more negative
values indicate higher essentiality).

| Column Name | Definition |
|----|----|
| ACH-XXXXX | DepMap cell line/sample identifier |
| \[Gene columns\] | Each column corresponds to a gene (formatted as “GeneSymbol (EntrezID)”); values represent gene effect scores |

### CRISPRInferredCommonEssentials.csv

**Used in:** `DataExtractionCode/EssentialityData.R`  
List of genes identified as commonly essential across multiple cell
lines, derived from CRISPR screening data.

| Column Name | Definition |
|----|----|
| Essentials | Gene identifier (formatted as “GeneSymbol (EntrezID)”) classified as commonly essential |

### final_dfA02_22_2026.csv

**Used in:**  
- `DataExtractionCode/DisgenetExtraction.R` -
`Functions/dataProcessing_df.R`  
Primary AlphaFold-based protein-level dataset used in disease
association analyses. Each row represents a protein–disease association.

This dataset contains columns previously defined in: -
`Human_AF_SASA.csv` (uniprotids) -
`final_resultsRAW05_21_2024.csv`(diseaseClasses_MSH, score,
disease_name, diseaseUMLSCUI, protein_class_name, protein_classid,
EntreID/geneid)

Additional columns:

| Column Name | Definition |
|----|----|
| Association_ID | Unique identifier for each gene–disease association record in DisGeNET (greater than 0) |
| Entanglement | Binary indicator of whether the protein contains at least one non-covalent lasso entanglement (NCLE) |
| Cov_Entanglement | Binary indicator of whether the protein contains covalent lasso entanglement |
| Knot | Binary indicator of whether the protein contains a topological knot (derived from structural topology analysis) |

Note: This file is not included in this repository due to the data
restrictions described above related to DisGeNET’s user agreement.

### final_dfA_Crystal02_22_2026.csv

**Used in:** - `DataExtractionCode/DisgenetExtraction.R` -
`Functions/dataProcessing_df.R` Crystal-based version of `final_dfA`
used for experimental structure comparisons.

This dataset contains the same columns and definitions as
`final_dfA02_22_2026.csv`, but is derived from experimentally determined
(crystal) structures instead of AlphaFold predictions.

Note: This file is not included in this repository due to the data
restrictions described above related to DisGeNET’s user agreement.

### final_dfB02_22_2026.csv

**Used in:** - `DataExtractionCode/DisgenetExtraction.R` -
`Functions/dataProcessing_df.R` Protein-level dataset integrating
structural features.

This dataset contains columns previously defined in:

- `Human_AF_combined_20250521 (1).csv` (entanglement and structural
  metrics)
- `Trovato.xlsx` (Travatos_G / maxGLN metric)
- `essentiality.csv` (Essential column)
- `final_dfA02_22_2026.csv` (Entanglement, Cov_Entanglement, Knot)

Additional columns:

| Column Name | Definition                                                   |
|-------------|--------------------------------------------------------------|
| cc_count    | Number of disulfide (cysteine–cysteine) bonds in the protein |
| Length      | Total protein length (number of residues)                    |

### final_dfB_Crystal02_22_2026.csv

**Used in:** - `DataExtractionCode/DisgenetExtraction.R` -
`Functions/dataProcessing_df.R` Crystal structure version of
`final_dfB`.

This dataset contains the same columns and definitions as
`final_dfB02_22_2026.csv`, but is derived from experimentally determined
(crystal) structures instead of AlphaFold predictions.

### final_dfB_af_02_22_2026.csv

**Used in:**  
- `DataExtractionCode/m1_data.R`  
- `DataExtractionCode/m2_data.R`  
- `DataExtractionCode/m3_data.R`  
AlphaFold protein-level dataset integrating structural features (no
pLDDT threshold).

This dataset contains columns previously defined in:

- `Human_AF_combined_20250521 (1).csv` (entanglement and structural
  metrics)
- `Trovato.xlsx` (Travatos_G / maxGLN metric)
- `essentiality.csv` (Essential column)
- `final_dfA02_22_2026.csv` (Entanglement, Cov_Entanglement, Knot)
- `final_dfB02_22_2026.csv` (cc_count, Length)

Additional columns:

| Column Name | Definition |
|----|----|
| organism | Organism name (e.g., Human) |
| plddt | Mean pLDDT score across the protein (AlphaFold confidence metric) |
| std | Standard deviation of pLDDT scores across the protein |

### final_dfB_max_af_70_raw_l_02_20.csv

**Used in:**  
- `AnalysisCodes/rf_and_match.R`  
- `DataExtractionCode/RF_Preprocessign.R`  
AlphaFold proteins filtered to pLDDT ≥ 70 with structural and length
metrics.

This dataset contains columns previously defined in:

- `Human_AF_combined_20250521 (1).csv` (entanglement metrics including
  Gn, Gc, Gmax, Gsum, threading, depth metrics)
- `Trovato.xlsx` (Travatos_G)
- `Human_AF_SASA.csv` (SASA-derived features)
- `final_dfB02_22_2026.csv` (Length, Entanglement, Cov_Entanglement)
- `essentiality.csv` (Essential)

| Column Name | Definition |
|----|----|
| Total_Mutations_n | Total number of missense mutations mapped to the protein |
| mean_SASA | Mean solvent-accessible surface area across all residues in the protein |
| disease | Binary indicator (“Yes”/“No”) of whether the protein is classified as disease-associated based on DisGeNET score (50th quantile) thresholds |

### final_df_af_02_22_2026.csv

**Used in:**  
- `DataExtractionCode/m1_data.R`  
- \`DataExtractionCode/m2_data.R Intermediate AlphaFold-only processed
dataset used during mutation data preparation.

This dataset contains columns previously defined in:

- `final_dfA02_22_2026.csv` (uniprotids, score, Entanglement,
  Cov_Entanglement, Knot)
- `final_dfB02_22_2026.csv` (Length)
- `essentiality.csv` (Essential)

Additional columns:

| Column Name | Definition |
|----|----|
| 95th_percentile | Binary indicator (“Yes”/“No”) of whether the protein’s DisGeNET score is above the 95th percentile threshold |
| 75th_percentile | Binary indicator (“Yes”/“No”) of whether the protein’s DisGeNET score is above the 75th percentile threshold |
| 50th_percentile | Binary indicator (“Yes”/“No”) of whether the protein’s DisGeNET score is above the 50th percentile threshold |

Note: This file is not included in this repository due to the data
restrictions described above related to DisGeNET’s user agreement.

### final_df_af_0_2026-2-6.csv

**Used in:** `DataExtractionCode/m3_data.R` Intermediate AlphaFold-only
processed dataset used during mutation data preparation.

This dataset contains the same columns and definitions as
`final_df_af_02_22_2026.csv`, except that it does not include the Length
and Essential columns.

Note: This file is not included in this repository due to the data
restrictions described above related to DisGeNET’s user agreement.

### combined_mapped_clean_missense.csv

**Used in:**  
- `DataExtractionCode/CleaningMutationData.R`  
- `DataExtractionCode/m1_data.R`  
- `DataExtractionCode/m2_data.R`  
- `DataExtractionCode/m3_data.R`  
Fully cleaned and mapped missense mutation dataset combining ClinVar and
gnomAD variants aligned to UniProt residues.

Variant annotation fields (e.g., genomic position, amino acid changes,
transcript mappings) are derived from the ProtVar annotation pipeline
(EMBL-EBI). See ProtVar documentation for full definitions of these
fields.

This dataset contains columns previously defined in:

- ProtVar annotation outputs (variant-level features)
- `final_dfA02_22_2026.csv` (uniprotids)

Additional columns:

| Column Name | Definition |
|----|----|
| source | Source database of the variant (e.g., ClinVar or gnomAD) |
| Mutation | Standardized mutation identifier (amino acid change format) |
| Pathogenic | Binary indicator (“Yes”/“No”) of whether the variant is classified as pathogenic |

### mutation_counts_af_70_2.csv

**Used in:**  
- `DataExtractionCode/RF_Preprocessign.R`  
- `DataExtractionCode/m2_data.R`  
Aggregated mutation counts per protein after filtering to pLDDT ≥ 70.

| Column Name | Definition |
|----|----|
| uniprotids | UniProt identifier for the protein |
| pathogenic_count | Total number of pathogenic missense variants mapped to the protein |
| benign_count | Total number of benign missense variants mapped to the protein |

### idmapping_2024_05_21.tsv

**Used in:** `DataExtractionCode/DisgenetExtraction.R`

### idmapping_2025_10_14.tsv

**Used in:** `DataExtractionCode/EssentialityData.R`

### idmapping_2_reviewed_true_AND_model_organ_2025_10_03.tsv

**Used in:** `DataExtractionCode/EssentialityData.R`

These files contain identifier mappings and protein metadata retrieved
from the UniProt ID mapping service.

They are used to map between:

- UniProt identifiers
- NCBI Entrez Gene identifiers
- Gene names and protein annotations

| Column Name | Definition |
|----|----|
| From | Input identifier (e.g., UniProt ID or Entrez Gene ID depending on query) |
| To | Mapped identifier (e.g., Entrez Gene ID or UniProt ID) |
| Entry | UniProt identifier for the protein |
| Reviewed | UniProt review status (“reviewed” for Swiss-Prot, “unreviewed” for TrEMBL) |
| Entry Name | UniProt entry name |
| Protein names | Full protein name(s) |
| Gene Names | Associated gene name(s) |
| Organism | Organism name (e.g., *Humans*) |
| Length | Protein length (number of residues); see definition in `final_dfB02_22_2026.csv` |

These mappings are used to harmonize identifiers across structural,
mutation, and disease datasets.

### uniprotkb_organism_id_9606_AND_keyword_2025_11_16.tsv

**Used in:** `AnalysisCodes/rf_and_match.R`  
UniProt metadata file containing basic protein annotations for human
proteins (organism ID: 9606). This dataset is used to retrieve protein
length.

| Column Name | Definition |
|----|----|
| Entry | UniProt identifier for the protein (see `uniprotids` above) |
| Length | Total protein length (number of residues); see definition in `final_dfB02_22_2026.csv` |
| Organism | Organism name (e.g., *Homo sapiens*) |

### cg_simulation_results_2026_03_09.xlsx

**Used in:** `AnalysisCodes/sim_mis_prop.R`  
Contains coarse-grained simulation results used to quantify misfolding
propensity.

This dataset contains columns previously defined in:

- `Human_AF_combined_20250521 (1).csv` (entanglement and structural
  metrics)
- `Trovato.xlsx` (Travatos_G)
- `final_dfB02_22_2026.csv` (Length, Entanglement, Cov_Entanglement)
- `final_dfB_max_af_70_raw_l_02_20.csv` (Total_Mutations_n, mean_SASA,
  disease)
- `essentiality.csv` (Essential)

Additional columns:

| Column Name | Definition |
|----|----|
| …1 | Row index or identifier |
| Delta_2RF | Difference in predicted disease probability between two random forest models (with vs. without entanglement features) |
| euclidean_dist_50 | Euclidean distance used for matching proteins in the top 50 selection |
| F_IDP | Fraction of intrinsically disordered residues in the protein |
| misfolding_prop | Estimated misfolding propensity from coarse-grained simulations |
| gene_ctrl | Matched control protein (UniProt identifier) |
| \*\_ctrl | Corresponding control protein features (same definitions as above, for matched non-entangled/non-disease proteins) |
| misfolding_prop_ctrl | Misfolding propensity for the matched control protein |

### Avg_pLDDTs.csv

**Used in:** `Functions/dataProcessing_df.R`  
Contains average pLDDT values per protein used for structural confidence
filtering.

| Column Name | Definition |
|----|----|
| organism | Organism name (e.g., *Humans*) |
| gene | UniProt identifier for the protein |
| <pLDDT> | Mean pLDDT score across the protein (structural confidence metric) |
| std | Standard deviation of pLDDT scores across the protein |

### Human_mappedCorrected_Lengths.xlsx

**Used in:** `DataExtractionCode/DisgenetExtraction.R` Contains
canonical protein lengths.

| Column Name | Definition |
|----|----|
| gene | UniProt identifier for the protein |
| AF Length | Protein length from AlphaFold model |
| Crystal Length | Protein length from experimental (crystal) structure |
| Crystal Mapped Length | Length of the region successfully mapped between crystal and AlphaFold structures |
| %mapped | Percentage of the crystal structure mapped to the AlphaFold sequence |
| AF - mappedCrystal | Difference between AlphaFold length and mapped crystal length |

## MutationData Subfolder

Contains raw and mapped mutation datasets used in
`DataExtractionCode/CleaningMutationData.R`:

### Raw datasets

- `clinvar_benign.tsv`
- `clinvar_pathogenic.tsv`
- `gnomad_benign`

These files contain raw missense variant data obtained from:

- ClinVar (clinical variant classifications)
- gnomAD (population variation data)

### Mapped datasets

- `clinvar_benign_mapped.csv`
- `clinvar_pathogenic_mapped.csv`
- `gnomad_benign_mapped.csv`

These files contain variants mapped to UniProt canonical protein
sequences.

Mapping and annotation are derived from the ProtVar pipeline, which
aligns genomic variants to protein residues and provides standardized
amino acid–level annotations.

## DataRanInAnalysis Folder

This folder contains the fully preprocessed datasets used immediately
before each analysis script is executed.

Each file represents the final cleaned, filtered, and harmonized
dataframe corresponding to a specific analysis (e.g., Q1, Q2, Q3, random
forest, etc.). These datasets ensure reproducibility by preserving the
exact inputs used to generate reported results.
