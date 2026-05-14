# New_Analysis

## Overview

This directory contains the manuscript-ready analysis pipeline for single-cell-resolution spatial transcriptomics analysis of TSC samples.

The maintained pipeline is the numbered script series:

- `01_Data_Preprocessing.py`
- `02_Manual_Annotation.py`
- `03_Manual_Annotation_Comb2.py`
- `04_Check_Gene_Expression.py`
- `05_SingleCell_DEGs.py`
- `06_Functional_Enrichment.py`
- `07_Spatial_Niche_NMF.py`
- `08_NMF_Downstream.py`
- `09_Venn_Factor7_DNExcit.py`
- `10_Core_Genes_Visualization.py`
- `11_GC_DN_Neighborhood_Communication.py`
- `12_CellChat_Analysis.R`
- `13_Spatial_Communication_Visualization.py`
- `14_NicheNet_Analysis.R`
- `15_Spatial_GRN_Squidpy.py`
- `16_Spatial_GRN_Comparison.py`
- `17_WGCNA_Analysis.R`
- `18_Trajectory_Pseudotime.py`
- `19_Trajectory_CNV.py`
- `20_Targeted_SCENIC_Upstream.py`
- `21_Targeted_NicheNet_Downstream.R`

## Repository Policy

- Only the scripts listed above are considered the maintained manuscript pipeline.
- Standard paths are centrally managed in `analysis_config.py` and `analysis_paths.R`.
- Downstream scripts are expected to consume standardized outputs from earlier steps instead of hard-coded file paths.
- Results are organized as `Results/<module>/<step_script>/...`.
- See `MANIFEST.md` for step-by-step commands, inputs, outputs, and environment entry points.
- See `ENVIRONMENT.md` for Python and R environment notes.
- See `renv.lock` for the current R-side package lock reference used for GitHub release preparation.
- Python scripts use the project-local virtual environment `ST_analysis/.venv`.
- R scripts use the user's global R installation and R package library.
- `17_WGCNA_Analysis.R` is part of the R workflow; it should be understood as an R script, not a Python script.

## Repository Layout

- `01-21` numbered scripts: maintained manuscript pipeline
- `tools/`: reusable helper utilities outside the core execution chain
- `supplementary/`: manuscript support scripts and notebooks not required for core reproduction
- `archive/`: historical or one-off materials kept only for reference
- `archive_dev/`: development-time checks, logs, temporary scripts, and test files

## Results Layout

The `Results` directory now uses continuous module numbering:

| Module | Meaning |
| --- | --- |
| `01_Basis` | preprocessing and base AnnData object |
| `02_Annotation` | annotation and expression validation |
| `03_Pathology_DEG_Enrichment` | pathology-centered DEG and enrichment analyses |
| `04_Niche_Neighborhood` | spatial niche, NMF, and core-gene analyses |
| `05_Communication` | neighborhood export, CellChat, NicheNet, communication visualization |
| `06_Regulation_Coexp` | spatial GRN and WGCNA |
| `07_Trajectory` | pseudotime and trajectory inference |
| `08_CNV` | CNV-related analyses |
| `09_Targeted_Networks` | targeted SCENIC and targeted NicheNet follow-up |

Within each module, each maintained script has its own subdirectory, for example:

```text
Results/
  01_Basis/
    01_Data_Preprocessing/
  02_Annotation/
    02_Manual_Annotation/
    03_Manual_Annotation_Comb2/
    04_Check_Gene_Expression/
  03_Pathology_DEG_Enrichment/
    05_SingleCell_DEGs/
    06_Functional_Enrichment/
  04_Niche_Neighborhood/
    07_Spatial_Niche_NMF/
    08_NMF_Downstream/
    09_Venn_Factor7_DNExcit/
    10_Core_Genes_Visualization/
```

## Standard Pipeline Order

Recommended execution order:

| Step | Script | Primary Input | Primary Output |
| --- | --- | --- | --- |
| 01 | `01_Data_Preprocessing.py` | raw spatial transcriptomics data | processed AnnData and cluster markers |
| 02 | `02_Manual_Annotation.py` | step 01 AnnData | basic annotation preview figures |
| 03 | `03_Manual_Annotation_Comb2.py` | step 01 AnnData | canonical annotated AnnData |
| 04 | `04_Check_Gene_Expression.py` | step 03 annotated AnnData | marker and spatial validation plots |
| 05 | `05_SingleCell_DEGs.py` | step 03 annotated AnnData | DE tables and comparison figures |
| 06 | `06_Functional_Enrichment.py` | step 03 annotated AnnData | enrichment and Moran results |
| 07 | `07_Spatial_Niche_NMF.py` | step 03 annotated AnnData | standardized NMF artifacts |
| 08 | `08_NMF_Downstream.py` | step 03 annotated AnnData and step 07 NMF artifacts | factor interpretation and enrichment |
| 09 | `09_Venn_Factor7_DNExcit.py` | step 05 DEGs and step 08 factor genes | overlap gene lists and Venn plots |
| 10 | `10_Core_Genes_Visualization.py` | step 03 annotated AnnData, step 05 DEGs, step 07 NMF artifacts, step 09 overlap genes | core-gene figures |
| 11 | `11_GC_DN_Neighborhood_Communication.py` | step 03 annotated AnnData | neighborhood composition and communication input export |
| 12 | `12_CellChat_Analysis.R` | step 11 communication export | CellChat results |
| 13 | `13_Spatial_Communication_Visualization.py` | step 03 annotated AnnData and step 12 CellChat output | spatial communication plots |
| 14 | `14_NicheNet_Analysis.R` | step 11 communication export | NicheNet ligand-receptor analyses |
| 15 | `15_Spatial_GRN_Squidpy.py` | step 03 annotated AnnData | pathological niche GRN output |
| 16 | `16_Spatial_GRN_Comparison.py` | step 03 annotated AnnData | comparative niche GRN output |
| 17 | `17_WGCNA_Analysis.R` | step 03 annotated AnnData | WGCNA modules and hub genes |
| 18 | `18_Trajectory_Pseudotime.py` | step 03 annotated AnnData | pseudotime and trajectory output |
| 19 | `19_Trajectory_CNV.py` | step 03 annotated AnnData and step 18 trajectory output | CNV-related output |
| 20 | `20_Targeted_SCENIC_Upstream.py` | step 16 GRN comparison output | targeted regulatory network summary |
| 21 | `21_Targeted_NicheNet_Downstream.R` | step 11 communication export | targeted downstream ligand-target results |

## Standardized Key Files

Important standardized outputs are centrally defined in:

- `analysis_config.py`
- `analysis_paths.R`

Examples:

- processed AnnData: `Results/01_Basis/01_Data_Preprocessing/adata_tsc_processed.h5ad`
- canonical annotated AnnData: `Results/02_Annotation/03_Manual_Annotation_Comb2/adata_tsc_annotated.h5ad`
- DEG directory: `Results/03_Pathology_DEG_Enrichment/05_SingleCell_DEGs/`
- NMF directory: `Results/04_Niche_Neighborhood/07_Spatial_Niche_NMF/`
- communication export directory: `Results/05_Communication/11_GC_DN_Neighborhood_Communication/CellPhoneDB_Input/`
- GRN comparison directory: `Results/06_Regulation_Coexp/16_Spatial_GRN_Comparison/`

## Legacy Content

The following content is retained for record-keeping only and is not part of the maintained manuscript pipeline:

- `原始代码/`
- `原始代码-260513/`
- `archive_dev/`
- `tools/`
- `supplementary/`
- `archive/`
- ad hoc checking scripts such as `check_*.py`
- temporary helper scripts such as `temp_*.py`
- exploratory notebooks and logs not referenced by the numbered maintained pipeline

Please do not mix outputs from legacy scripts with outputs from the maintained `01-21` pipeline.

## Reproducibility Notes

- Re-run the pipeline from step 01 when regenerating manuscript results after code updates.
- Use the standardized outputs from earlier steps rather than manually copying files across directories.
- Python and R steps share the same result layout through `analysis_config.py` and `analysis_paths.R`.
- Python environment management and R environment management are intentionally separated in this project.
