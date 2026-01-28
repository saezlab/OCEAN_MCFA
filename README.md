# OCEAN_MCFA — Multicellular Factor Analysis (Snakemake)

This repository contains a Snakemake workflow implementing a Multicellular Factor Analysis (MCFA) of 10x single-nucleus RNA-seq data from the OCEAN Human Glomerular Disease atlas.

The analysis is organized as a reproducible Snakemake workflow and was used in the study: "A Human Glomerular Disease Atlas defines the APOL1-JAK-STAT feed forward loop in focal segmental glomerulosclerosis" (preprint): https://www.medrxiv.org/content/10.1101/2025.09.12.25335572v3

## Repository overview

- `Snakefile` — the top-level Snakemake workflow.
- `workflow/` — Snakemake rules and scripts.
  - `rules/` — rule definitions (e.g., `mofacell.smk`, `preprocessing.smk`, `robustness.smk`).
  - `scripts/` — helper Python and R scripts used by rules (grouped by purpose: `mofacell`, `preprocessing`, `figures`, `downstream`, ...).
- `config/` — workflow configuration files
  - `config.yaml` - workflow-specific parameter configuration.
  - `slurm/` - subfolder with Slurm helpers and job scripts.
- `data/` — input data (note: some large raw files are not included in the repo; see below).
  - `data/OCEAN_v3_Nu_102025a_CZI.h5ad` — the primary processed single-nucleus dataset used in the workflow (if present).
  - `data/metadata.csv` and additional metadata/spreadsheets.
- `results/` — generated outputs (scores, loadings, heatmaps, downstream tables, etc.).
- `figures/`, `plots/`, `supplementary_tables/` — outputs and manuscript materials.
- `LICENSE` — repository license.

## Goal / Scope

Reproduce the MCFA analysis and downstream figures used in the manuscript. The workflow handles preprocessing, running MOFA/MCFA factor analysis across 10x snRNA-seq samples, producing factor loadings and scores, downstream statistics, and final figures.

## Requirements

- Snakemake v7.30.2 (this workflow was implemented and tested with this version).
- Cluster access (Slurm). A `config/slurm/` folder contains example jobscript and utilities for submitting to Slurm, which may need adapting to fulfill job submission requirements of your Slurm system.
- Singularity images in SIF format (see below).
- Data (see below).

## Configuration

Primary configuration files live under `config/`:

- `config/config.yaml` — main workflow parameters (samples, settings, thresholds, etc.).
- `config/slurm/` — helper scripts and `slurm-jobscript.sh` for cluster execution.

Before running, inspect and edit `config/config.yaml` to point to the correct input paths for your system.

## Data

The repository references several data files under `data/`. This repository will be updated with access links once the files are readily available on Zenodo.

- `data/OCEAN_v3_Nu_102025a_CZI.h5ad` — processed AnnData used for downstream analysis.
- `data/metadata.csv` — additional sample metadata.

## Singularity images

The respository references several Singularity images to ensure full reproducibility of this workflow. This repository will be updated with access links once the files are readily available on Zenodo.

## Running the workflow (examples)

1) Dry run (see what would be executed):

```zsh
snakemake -n
```

2) Run on a Slurm cluster using the provided Slurm scripts (example):

```zsh
# Adapt config/slurm/config.yaml for your cluster, then:
snakemake --profile config/slurm
```

## Outputs

Results are written under `results/` with subfolders for preprocessing, `mofacell`, downstream tables, and figure-ready files. Examples:

- `results/mofacell/` — loadings, scores, R2 statistics, and saved model instances.
- `results/downstream/` — factor-metadata associations, enrichment analysis.
- `results/preprocessing/` - centered log ratio-transformed cell type proportions.

## Reproducing figures

Figure generation scripts are in `workflow/scripts/figures/`. Many R scripts assume the results are in `results/` and will load precomputed tables. Run the Snakemake workflow up to the figure-producing targets first, then use the R scripts to create the figure PDFs.

## How to cite us

Read and cite:

https://www.medrxiv.org/content/10.1101/2025.09.12.25335572v3

## Acknowledgements
The authors acknowledge support by the state of Baden-Württemberg through bwHPC and the German Research Foundation (DFG) through grant INST 35/1597-1 FUGG, as well as the data storage service SDS@hd supported by the Ministry of Science, Research and the Arts Baden-Württemberg (MWK) and the German Research Foundation (DFG) through grant INST 35/1503-1 FUGG. Charlotte Boys gratefully acknowledges DFG funding through the Clinical Research Unit 5011 InteraKD (Project ID: 445703531).

## License

This repository is provided under the GNU General Public License v3.0.

