# Pipeline Manifest

## Purpose

This manifest lists the maintained manuscript pipeline scripts, their standard execution order, key inputs, key outputs, recommended commands, and environment entry points.

Project root is assumed to be:

```text
ST_analysis/
```

## Environments

- Python scripts: use `ST_analysis/.venv`
- R scripts: use the global/system R installation with R packages listed in `ENVIRONMENT.md`
- `17_WGCNA_Analysis.R`: remains an R script and should be run with `Rscript`; if `reticulate` is used for H5AD access, it only bridges to the project Python interpreter and does not change the fact that the script belongs to the R environment

## Maintained Scripts

| Step | Script | Language | Standard command | Primary input | Primary output |
| --- | --- | --- | --- | --- | --- |
| 01 | `01_Data_Preprocessing.py` | Python | `python New_Analysis/01_Data_Preprocessing.py` | raw ST data under sample directories | `Results/01_Basis/01_Data_Preprocessing/` |
| 02 | `02_Manual_Annotation.py` | Python | `python New_Analysis/02_Manual_Annotation.py` | step 01 processed AnnData | `Results/02_Annotation/02_Manual_Annotation/` |
| 03 | `03_Manual_Annotation_Comb2.py` | Python | `python New_Analysis/03_Manual_Annotation_Comb2.py` | step 01 processed AnnData | `Results/02_Annotation/03_Manual_Annotation_Comb2/` |
| 04 | `04_Check_Gene_Expression.py` | Python | `python New_Analysis/04_Check_Gene_Expression.py` | step 03 annotated AnnData | `Results/02_Annotation/04_Check_Gene_Expression/` |
| 05 | `05_SingleCell_DEGs.py` | Python | `python New_Analysis/05_SingleCell_DEGs.py` | step 03 annotated AnnData | `Results/03_Pathology_DEG_Enrichment/05_SingleCell_DEGs/` |
| 06 | `06_Functional_Enrichment.py` | Python | `python New_Analysis/06_Functional_Enrichment.py` | step 03 annotated AnnData | `Results/03_Pathology_DEG_Enrichment/06_Functional_Enrichment/` |
| 07 | `07_Spatial_Niche_NMF.py` | Python | `python New_Analysis/07_Spatial_Niche_NMF.py` | step 03 annotated AnnData | `Results/04_Niche_Neighborhood/07_Spatial_Niche_NMF/` |
| 08 | `08_NMF_Downstream.py` | Python | `python New_Analysis/08_NMF_Downstream.py` | step 03 annotated AnnData and step 07 NMF artifacts | `Results/04_Niche_Neighborhood/08_NMF_Downstream/` |
| 09 | `09_Venn_Factor7_DNExcit.py` | Python | `python New_Analysis/09_Venn_Factor7_DNExcit.py` | step 05 DEGs and step 08 factor summaries | `Results/04_Niche_Neighborhood/09_Venn_Factor7_DNExcit/` |
| 10 | `10_Core_Genes_Visualization.py` | Python | `python New_Analysis/10_Core_Genes_Visualization.py` | steps 03, 05, 07, 09 outputs | `Results/04_Niche_Neighborhood/10_Core_Genes_Visualization/` |
| 11 | `11_GC_DN_Neighborhood_Communication.py` | Python | `python New_Analysis/11_GC_DN_Neighborhood_Communication.py` | step 03 annotated AnnData | `Results/05_Communication/11_GC_DN_Neighborhood_Communication/` |
| 12 | `12_CellChat_Analysis.R` | R | `Rscript New_Analysis/12_CellChat_Analysis.R` | step 11 `CellPhoneDB_Input` export | `Results/05_Communication/12_CellChat_Analysis/` |
| 13 | `13_Spatial_Communication_Visualization.py` | Python | `python New_Analysis/13_Spatial_Communication_Visualization.py` | step 03 annotated AnnData and step 12 CellChat output | `Results/05_Communication/13_Spatial_Communication_Visualization/` |
| 14 | `14_NicheNet_Analysis.R` | R | `Rscript New_Analysis/14_NicheNet_Analysis.R` | step 11 `CellPhoneDB_Input` export | `Results/05_Communication/14_NicheNet_Analysis/` |
| 15 | `15_Spatial_GRN_Squidpy.py` | Python | `python New_Analysis/15_Spatial_GRN_Squidpy.py` | step 03 annotated AnnData | `Results/06_Regulation_Coexp/15_Spatial_GRN_Squidpy/` |
| 16 | `16_Spatial_GRN_Comparison.py` | Python | `python New_Analysis/16_Spatial_GRN_Comparison.py` | step 03 annotated AnnData | `Results/06_Regulation_Coexp/16_Spatial_GRN_Comparison/` |
| 17 | `17_WGCNA_Analysis.R` | R | `Rscript New_Analysis/17_WGCNA_Analysis.R` | step 03 annotated AnnData | `Results/06_Regulation_Coexp/17_WGCNA_Analysis/` |
| 18 | `18_Trajectory_Pseudotime.py` | Python | `python New_Analysis/18_Trajectory_Pseudotime.py` | step 03 annotated AnnData | `Results/07_Trajectory/18_Trajectory_Pseudotime/` |
| 19 | `19_Trajectory_CNV.py` | Python | `python New_Analysis/19_Trajectory_CNV.py` | step 03 annotated AnnData and step 18 trajectory output | `Results/08_CNV/19_Trajectory_CNV/` |
| 20 | `20_Targeted_SCENIC_Upstream.py` | Python | `python New_Analysis/20_Targeted_SCENIC_Upstream.py` | step 16 GRN comparison output | `Results/09_Targeted_Networks/20_Targeted_SCENIC_Upstream/` |
| 21 | `21_Targeted_NicheNet_Downstream.R` | R | `Rscript New_Analysis/21_Targeted_NicheNet_Downstream.R` | step 11 `CellPhoneDB_Input` export | `Results/09_Targeted_Networks/21_Targeted_NicheNet_Downstream/` |

## Standardized Key Files

| Name | Path |
| --- | --- |
| processed AnnData | `Results/01_Basis/01_Data_Preprocessing/adata_tsc_processed.h5ad` |
| canonical annotated AnnData | `Results/02_Annotation/03_Manual_Annotation_Comb2/adata_tsc_annotated.h5ad` |
| DEG output root | `Results/03_Pathology_DEG_Enrichment/05_SingleCell_DEGs/` |
| NMF output root | `Results/04_Niche_Neighborhood/07_Spatial_Niche_NMF/` |
| communication export root | `Results/05_Communication/11_GC_DN_Neighborhood_Communication/CellPhoneDB_Input/` |
| GRN comparison root | `Results/06_Regulation_Coexp/16_Spatial_GRN_Comparison/` |
| WGCNA root | `Results/06_Regulation_Coexp/17_WGCNA_Analysis/` |
| trajectory root | `Results/07_Trajectory/18_Trajectory_Pseudotime/` |
| CNV root | `Results/08_CNV/19_Trajectory_CNV/` |
| targeted network root | `Results/09_Targeted_Networks/` |

## Notes

- The maintained pipeline is the numbered `01-21` series only.
- Path resolution is centralized in `analysis_config.py` and `analysis_paths.R`.
- Legacy scripts, debug scripts, temporary scripts, test files, and historical logs are archived and should not be mixed into manuscript runs.
- Do not interpret `.venv` as the runtime environment for R scripts; it is only the Python environment used by Python scripts and any explicit `reticulate` bridge calls.
- Non-pipeline helper files are organized under `tools/`, `supplementary/`, `archive/`, and `archive_dev/` to keep the repository front page focused on the maintained workflow.
