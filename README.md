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

## Citation / reference

Developed for visualizing TF target and DAP-seq–related analyses. Adjust inputs for your organism and experiment; the `test/` files use *Arabidopsis* TF naming from the bundled example database.
