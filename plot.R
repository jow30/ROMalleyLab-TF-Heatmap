library(clusterProfiler)
library(AnnotationHub)
library(ggplot2)
library(dplyr)
library(readxl)

working_dir <- "/Volumes/project/gzy8899/qiaoshan/bulkRNAseq/experiments/Hannah/"
setwd(working_dir)

hub <- AnnotationHub()
all_sp <- unique(hub$species) # can be a drop-down option for users to select the species of interest for enrichment analysis, e.g., "Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella", "Brassica oleracea"

### user input
# mandatory input
species <- c("Arabidopsis thaliana") # user input species of interest, e.g., "Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella", "Brassica oleracea"
goi <- "GSAD/selected_DEGs.txt" # user input gene query (DEGs or genes of interest)
out_prefix  <- "~/Documents/Projects/Ronan/GO_enrichment" # user input output file prefix for enrichment results
# optional input
gsmap <- "~/Documents/Projects/Ronan/RNAseq/Ara_ALL_t2g.txt" # user input gene set mapping file (GeneSetID, GeneID); an optional third column with GeneSet category can be provided for gsoi-specific enrichment analysis (e.g., only analyze TF targets in the "BP" category). If null, enrichment analysis will be performed based on the species database (OrgDb) and the specified gene set (ontology).
gsoi <- "CRG" # user input ontology/geneset for enrichment analysis, e.g., "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component, "ALL" for all three ontologies, or a specific gene set category in the gsmap file (e.g., "BP" for only analyze gene sets in the "BP" category in the gsmap file), cannot be null if gsmap is null
bg <- "protein" # user input background genes, e.g., "protein" for all protein-coding genes in the species database, "all" for all genes in the species database, or a custom background gene list (a text file with one gene ID per line, with the same key type as the query genes)
pvalue_cutoff = 0.05
fdr_cutoff = 0.1
min_gs_size = 5
max_gs_size = 5000
### end of user input

# load species database from AnnotationHub
qry <- query(hub, species)
qry_mat <- mcols(qry)[, c("title", "rdataclass", "description")] %>% as.matrix()
hub_id <- rownames(qry_mat)[qry_mat[,2]=="OrgDb"]
sp_db <- hub[[hub_id]]

# determine key type for the species database
goi <- unique(trimws(readLines(goi)))
lapply(columns(sp_db), function(c) sum(goi %in% keys(sp_db, keytype = c))) %>% 
  unlist() %>% 
  setNames(columns(sp_db)) %>% 
  sort(decreasing = TRUE) %>% 
  head(1) -> top_key_type
if (top_key_type == 0) {
  # print out the top 5 IDs for each key type in the species database to help users identify the correct key type for their query genes
  message("None of the query genes match any key type in the species database. Here are the top 5 IDs for each key type in the species database to help you identify the correct key type for your query genes:")
  for (c in columns(sp_db)) {
    message("Key type: ", c)
    message("Top 5 IDs: ", paste(head(keys(sp_db, keytype = c)), collapse = ", "))
  }
  stop("Please check your gene identifiers.")
} else {
  message("Use key type ", names(top_key_type), " (with ", top_key_type, " out of ", length(goi), " genes of interest matches).")
  goi <- goi[goi %in% keys(sp_db, keytype = names(top_key_type))] # filter query genes to keep only those that match the top key type in the species database
}
names(top_key_type) -> key_type

# select background genes (all genes in the species database with the same key type as the query genes)
if (bg == "protein") {
  if ("REFSEQ" %in% columns(sp_db)) {
    bg_genes <- grep("[N|X]P_", keys(sp_db, keytype = "REFSEQ"), value = T) # assuming protein-coding genes have REFSEQ IDs starting with "NP_" or "XP_")
    bg_genes <- AnnotationDbi::select(sp_db, keys = bg_genes, columns = c("REFSEQ", key_type), keytype = "REFSEQ") %>% as.data.frame() %>% dplyr::select(all_of(key_type)) %>% unlist() %>% unique() %>% na.omit
    message("Using ", length(bg_genes), " protein-coding genes with a REFSEQ prefix of NP_/XP_ as background.")
  } else {
    bg_genes <- keys(sp_db, keytype = key_type)
    message("REFSEQ key type not found in the species database. Protein-coding genes cannot be extracted. Using all ", length(bg_genes), " genes with key type ", key_type, " as background instead.")
  }
} else if (bg == "all") {
  bg_genes <- keys(sp_db, keytype = key_type)
  message("Using all ", length(bg_genes), " genes with key type ", key_type, " as background.")
} else if (file.exists(bg)) {
  bg_genes <- unique(trimws(readLines(bg)))
  bg_genes_raw_count <- length(bg_genes)
  bg_genes <- bg_genes[bg_genes %in% keys(sp_db, keytype = key_type)]
  message("Using ", length(bg_genes), " out of ", bg_genes_raw_count, " custom background genes that matches ", key_type, " IDs.")
} else {
  stop("Invalid background gene option. Please choose either 'protein' for all protein-coding genes,  or 'all' for all genes in the species database, or provide a custom background gene list (a text file with one gene ID per line, with the same key type as the query genes).")
}

goi <- goi[goi %in% bg_genes] # make sure genes of interest are in the background gene set

enrich <- function(query_genes, species_db = NULL, gsmap = NULL, key_type = NULL, background_genes = NULL, gene_set = NULL, pvalue_cutoff = 0.05, fdr_cutoff = 0.1, minGSSize = 5, maxGSSize = 5000, out_file = NULL) {
  if (!is.null(gsmap)) {
    if (!is.null(gene_set)) {
      gsmap <- gsmap[gsmap[[3]] == gene_set, ]
    }
    ego <- enricher(gene = query_genes,
                    TERM2GENE = gsmap, 
                    pAdjustMethod = "BH",
                    pvalueCutoff = pvalue_cutoff,
                    qvalueCutoff = fdr_cutoff,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    universe = background_genes)
  } else if (!is.null(species_db) & !is.null(key_type) & !is.null(gene_set)) {
    ego <- enrichGO(gene = query_genes,
                    OrgDb = species_db,
                    keyType = key_type,
                    ont = gene_set,
                    pAdjustMethod = "BH",
                    pvalueCutoff = pvalue_cutoff,
                    qvalueCutoff = fdr_cutoff,
                    minGSSize = minGSSize,
                    maxGSSize = maxGSSize,
                    universe = background_genes,
                    readable = TRUE)
  } else {
    stop("Please provide either a species database (OrgDb) or a gene set mapping file (gsmap). When using a species database, please also specify the key type and gene set (ontology) for enrichment analysis.")
  }
  
  # Save results to Excel
  if (!is.null(ego) && nrow(ego) > 0) {
    if (!is.null(out_file)) {
      write.xlsx(as.data.frame(ego), file = out_file, row.names = FALSE)
      message("Enrichment results saved to ", out_file)
    } else {
      message(nrow(ego), " significant enrichments found.")
      print(ego)
    }
  } else {
    message("No significant enrichment found.")
  }
  return(ego)
}

# Perform enrichment analysis
if (!is.null(gsmap)) {
  gsmap <- read.delim(gsmap, header = T, stringsAsFactors = F)
  gsmap <- gsmap[gsmap[[2]] %in% bg_genes, ] # make sure genes in gsmap are in the background gene set
  
  ora_res <- enrich(goi, gsmap = gsmap, background_genes = bg_genes, gene_set = gsoi, 
                    pvalue_cutoff = pvalue_cutoff, fdr_cutoff = fdr_cutoff, minGSSize = min_gs_size, maxGSSize = max_gs_size, 
                    out_file = paste0(paste0(c(out_prefix, gsoi), collapse = "_"), ".xlsx"))
} else {
  ora_res <- enrich(goi, species_db = sp_db, key_type = key_type, background_genes = bg_genes, gene_set = gsoi, 
                    pvalue_cutoff = pvalue_cutoff, fdr_cutoff = fdr_cutoff, minGSSize = min_gs_size, maxGSSize = max_gs_size,
                    out_file = paste0(paste0(c(out_prefix, gsoi), collapse = "_"), ".xlsx"))
}

