# Environment Notes

## Scope

This document describes the manuscript environment used by the maintained `New_Analysis` pipeline.

This project uses two distinct runtime environments:

- Python scripts use the project virtual environment `ST_analysis/.venv`
- R scripts use the user's global/system R installation

Python version and Python package versions are derived from the current project virtual environment:

- `.venv`
- path: `ST_analysis/.venv`
- Python version: `3.13.1`

The exact Python package lock file exported from the current `.venv` is:

- `New_Analysis/requirements-python-venv-freeze.txt`

## Python Setup

From the project root:

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r New_Analysis\requirements-python-venv-freeze.txt
```

If `py -3.13` is unavailable, install Python `3.13.1` first and then create the virtual environment again.

## Python Core Packages

The following packages are especially important for this pipeline and were read from the current `.venv`:

| Package | Version |
| --- | --- |
| `python` | `3.13.1` |
| `anndata` | `0.12.10` |
| `scanpy` | `1.12` |
| `squidpy` | `1.8.1` |
| `numpy` | `2.4.3` |
| `scipy` | `1.17.1` |
| `pandas` | `2.3.3` |
| `matplotlib` | `3.10.8` |
| `seaborn` | `0.13.2` |
| `scikit-learn` | `1.8.0` |
| `harmonypy` | `0.0.10` |
| `gseapy` | `1.1.12` |
| `infercnvpy` | `0.6.1` |
| `commot` | `0.0.3` |
| `omnipath` | `1.0.12` |
| `python-igraph` | `1.0.0` |
| `leidenalg` | `0.11.0` |
| `matplotlib-venn` | `1.1.2` |
| `adjustText` | `1.3.0` |

## R Setup

R scripts are executed with the global/system `Rscript`. Their package installation and library paths are managed on the R side and are not taken from `.venv`.

Current detected global R version:

- `R version 4.5.1 (2025-06-13 ucrt)`
- R package lock file prepared for GitHub upload: `New_Analysis/renv.lock`

The maintained R scripts require at least the following packages:

- `CellChat`
- `patchwork`
- `ggplot2`
- `dplyr`
- `Matrix`
- `nichenetr`
- `Seurat`
- `tidyverse`
- `data.table`
- `WGCNA`
- `stringr`
- `reticulate`

Detected key package versions in the current global R environment:

| Package | Version |
| --- | --- |
| `CellChat` | `1.6.1` |
| `patchwork` | `1.3.1` |
| `ggplot2` | `4.0.1` |
| `dplyr` | `1.1.4` |
| `Matrix` | `1.7.3` |
| `nichenetr` | `2.2.1.1` |
| `Seurat` | `5.3.0` |
| `tidyverse` | `2.0.0` |
| `data.table` | `1.17.6` |
| `WGCNA` | `1.73` |
| `stringr` | `1.5.1` |
| `reticulate` | `1.42.0` |

Example installation:

```r
install.packages(c(
  "patchwork", "ggplot2", "dplyr", "Matrix", "Seurat",
  "tidyverse", "data.table", "WGCNA", "stringr", "reticulate"
))
```

Some packages such as `CellChat` and `nichenetr` may need GitHub or Bioconductor-style installation depending on the local R environment.

## WGCNA Note

`17_WGCNA_Analysis.R` is an R script and belongs to the R-side workflow.

Its primary runtime is:

- global/system `Rscript`
- R packages such as `WGCNA`, `Seurat`, `reticulate`, `ggplot2`, and `dplyr`

When needed, it may use:

- `reticulate`
- Python interpreter: `ST_analysis/.venv/Scripts/python.exe`

This Python bridge is only used to support R-side access to Python objects such as H5AD data. It does not mean that the R scripts are managed by the Python virtual environment.

## Reproducibility Recommendation

- Use the exact `requirements-python-venv-freeze.txt` file for Python reproduction.
- Avoid mixing a different conda environment with the maintained pipeline.
- Keep R package installation and Python virtual environment management separate.
- `renv.lock` is included as a reproducibility reference for the current global R package state.
