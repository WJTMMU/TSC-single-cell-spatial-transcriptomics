# TSC Single-Cell-Resolution Spatial Transcriptomics Analysis

This repository provides the manuscript-ready computational workflow for single-cell-resolution spatial transcriptomics analysis of tuberous sclerosis complex (TSC) samples. The pipeline integrates preprocessing, cell-type annotation, pathology-centered differential expression, spatial niche deconvolution, cell-cell communication analysis, gene regulatory network inference, trajectory analysis, CNV-associated interpretation, and targeted downstream prioritization into a reproducible `01-21` step framework.

The codebase is organized for manuscript submission and GitHub-based reproduction: maintained analysis scripts are separated from legacy materials, outputs are standardized under module-level result directories, and Python and R environments are documented independently for transparent reuse.

## Citation

- First author: `Ju Wu`
- Please cite the associated manuscript together with this software repository.

## Repository Contents

- `.gitignore`: excludes local environments, results, logs, caches, and large intermediate files
- `LICENSE`: `BSD-3-Clause`
- `CITATION.cff`: citation metadata for the repository
- `New_Analysis/`: maintained analysis pipeline and documentation

## Maintained Workflow

The maintained workflow is the numbered `01-21` script series under `New_Analysis/`.

Key documentation:

- `New_Analysis/README.md`
- `New_Analysis/MANIFEST.md`
- `New_Analysis/ENVIRONMENT.md`

## Environment

- Python scripts use the project virtual environment definition captured in:
  - `New_Analysis/requirements-python-venv-freeze.txt`
- R scripts use the global/system R installation and the package lock reference:
  - `New_Analysis/renv.lock`

## Not Included

The following are intentionally excluded from the upload-ready repository:

- local virtual environments such as `.venv/`
- generated outputs under `New_Analysis/Results/`
- development checks, scratch logs, temporary files, and local caches
- legacy and archived materials not required for core manuscript reproduction

## Reproducibility

Please start from:

- `New_Analysis/README.md` for overview
- `New_Analysis/MANIFEST.md` for execution order
- `New_Analysis/ENVIRONMENT.md` for environment setup
