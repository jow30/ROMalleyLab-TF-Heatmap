# TF Target Heatmap (DAP-seq)

Interactive Shiny app to visualize transcription factor (TF) target genes occurrence and conservation scores in heatmaps. It shows which genes of interest appear in the TF target sets of interest, optional TF/gene annotations, expression or logFC profiles, and conservation scores from the latest DAP-seq target database (by 2026 May).

## Quick start

### CRI Bio-core version (online)

1. Open [https://biocoreapps.bsd.uchicago.edu/apps/app/tfheatmap](https://biocoreapps.bsd.uchicago.edu/apps/app/tfheatmap) in your browser.

2. Login using your cnetID and password.

### Docker version (offline)

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).

2. Search for `jonalin/geneset-heatmap` in the Docker Hub and click `Run` to start the app.

3. Set host port to `3838` and click `Run`. The Docker image will be pulled and the app will be started in a container.

4. Open **http://localhost:3838** in your browser.

## Try the app with example data

Example inputs are in the [**`test/`**](https://github.com/jow30/ROMalleyLab-TF-Heatmap/tree/main/test) folder. Use them to learn the expected file formats before uploading your own data.

Example inputs can be quickly loaded by clicking `Load all example files` at the top of the sidebar.

## Example files

| App field | Example file | Description |
|-----------|--------------|-------------|
| TF target database | `TF_targets_At_cscore.txt` | Tab-separated TF–gene database with conservation scores (50 MB max) |
| Gene sets to plot | `GSAD_results_all.txt` | One TF/gene-set ID per line |
| Genes of interest | `selected_DEGs.txt` | One gene ID per line (e.g. DEGs) |
| TF annotation | `ath-258-tf-info_simple_pp.txt` | TF ID, display label, family/group (optional) |
| TF profiles | `logFC.txt` | TF logFC/expression across samples/conditions; first column = TF ID (optional) |
| Gene annotation | `gene_id_name_mapping.txt` | Gene ID, display name, optional group (optional) |
| Gene profiles | `genes_tpm.txt` | Gene logFC/expression across samples/conditions; first column = gene ID (optional) |

## Step-by-step usage

### Required inputs

1. **TF target database** — Tab-separated file with at least `gene`, `tf`, and score columns (e.g. `n_cons_species_minfrac0`). May include a `Category` column for subsetting.
2. **Gene sets to plot** — Plain text, one TF (gene-set) ID per line. IDs must exist in the database.
3. **Genes of interest** — Plain text, one gene ID per line.

### Options

- **Subset by Category** — Comma-separated values (e.g. `BP,TFT`) if the database has a `Category` column.
- **Occurrence cutoff** — `0` = all genes in selected sets; `1` = only genes of interest in those sets; `>1` = genes in at least N sets.
- **Diet mode** — Randomly subsample to 100 genes when the heatmap is very large.
- **Conservation score column** — Column from the database for the score to plot. 

## Input file formats (summary)

- **Lists** (`geneset`, `goi`): one ID per line, no header required.
- **Database**: tab-separated, header row, columns `gene`, `tf`, plus score and optional `Category`.
- **Annotation tables**: tab-separated; column 1 = ID matching the database; column 2 = label; column 3 = group (optional).
- **Profile tables**: tab-separated; column 1 = ID matching the database; remaining columns = numeric values (logFC, TPM, etc.).

Maximum upload size per file: **50 MB**.

## Run without Docker (local R)

Requires R ≥ 4.4 with packages: `shiny`, `ComplexHeatmap`, `circlize`, `dplyr`, `tibble`, `tidyr` (see `install.R`).

```bash
cd TFHeatmap
Rscript install.R   # once
R -e "shiny::runApp('.')"
```

Open the URL shown in the terminal (usually http://127.0.0.1:3838).

## macOS app (optional)

If you use Docker Desktop, you can build a double-clickable launcher:

```bash
docker build -t tf-heatmap-shiny .
./macos/build_app.sh
```

Open **`TFHeatmap.app`** in the project folder. It starts the container and opens the app in your browser.

## Project layout

```
app.R                 # Shiny UI and server
R/heatmap_core.R      # Heatmap analysis logic
R/ui_helpers.R        # Example file paths and UI helpers
test/                 # Example input files (start here)
Dockerfile            # Container build
GSAD_DAPseq_TFs.R     # Original batch script (reference)
```

## Troubleshooting

| Issue | Suggestion |
|-------|------------|
| Upload fails / file too large | Per-file limit is 100 MB; compress or subset very large tables |
| Gene sets not in database | Check TF IDs in `geneset` file match the `tf` column in the database |
| Empty score / split tabs | Select a score column, or provide TF annotation for split-by-family views |
| Docker build network error | Run `docker pull rocker/r-ver:4.4.1` then `docker build --pull=false -t tf-heatmap-shiny .` |

## Citation / reference

Developed for visualizing TF target and DAP-seq–related analyses. Adjust inputs for your organism and experiment; the `test/` files use *Arabidopsis* TF naming from the bundled example database.
