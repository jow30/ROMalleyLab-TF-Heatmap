library(clusterProfiler)
library(AnnotationHub)
library(ggplot2)
library(dplyr)
library(readxl)
library(writexl)

working_dir <- "/Volumes/project/gzy8899/qiaoshan/bulkRNAseq/experiments/Hannah/"
setwd(working_dir)

hub <- AnnotationHub()
all_sp <- unique(hub$species) # can be a drop-down option for users to select the species of interest for enrichment analysis, e.g., "Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella", "Brassica oleracea"

### user input
# mandatory input
species <- c("Arabidopsis thaliana") # user input species of interest, e.g., "Arabidopsis thaliana", "Arabidopsis lyrata", "Capsella rubella", "Brassica oleracea"
goi_tb <- "~/Documents/Projects/Ronan/gsea_input.txt" # user input gene query (DEGs or genes of interest)
out_prefix  <- "~/Documents/Projects/Ronan/GSEA" # user input output file prefix for enrichment results
# optional input
gsmap <- "~/Documents/Projects/Ronan/RNAseq/Ara_ALL_t2g.txt" # user input gene set mapping file (GeneSetID, GeneID); an optional third column with GeneSet category can be provided for gsoi-specific enrichment analysis (e.g., only analyze TF targets in the "BP" category). If null, enrichment analysis will be performed based on the species database (OrgDb) and the specified gene set (ontology).
gsoi <- "BP" # user input ontology/geneset for enrichment analysis, e.g., "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component, "ALL" for all three ontologies, or a specific gene set category in the gsmap file (e.g., "BP" for only analyze gene sets in the "BP" category in the gsmap file), cannot be null if gsmap is null
bg <- "protein" # user input background genes, e.g., "protein" for all protein-coding genes in the species database, "all" for all genes in the species database, or a custom background gene list (a text file with one gene ID per line, with the same key type as the query genes)
pvalue_cutoff = 0.05 # adjusted p-value (FDR) cutoff for significant enrichment, default is 0.05
qvalue_cutoff = 0.1
min_gs_size = 5
max_gs_size = 5000
### end of user input

# load species database from AnnotationHub
qry <- query(hub, species)
qry_mat <- mcols(qry)[, c("title", "rdataclass", "description")] %>% as.matrix()
hub_id <- tail(rownames(qry_mat)[qry_mat[,2]=="OrgDb"], 1)
sp_db <- hub[[hub_id]]

# determine key type for the species database
goi_tb <- read.table(goi_tb, header = T, stringsAsFactors = F)
goi_tb[[1]] <- trimws(goi_tb[[1]]) # remove leading and trailing whitespace from gene IDs
lapply(columns(sp_db), function(c) sum(goi_tb[[1]] %in% keys(sp_db, keytype = c))) %>% 
  unlist() %>% 
  setNames(columns(sp_db)) %>% 
  sort(decreasing = TRUE) %>% 
  head(1) -> top_key_type
if (top_key_type == 0 & is.null(gsmap) & file.exists(bg) == FALSE) {
  # print out the top 5 IDs for each key type in the species database to help users identify the correct key type for their query genes
  message("None of the query genes match any key type in the species database. Here are the top 5 IDs for each key type in the species database to help you identify the correct key type for your query genes:")
  for (c in columns(sp_db)) {
    message("Key type: ", c)
    message("Top 5 IDs: ", paste(head(keys(sp_db, keytype = c)), collapse = ", "))
  }
  stop("Please check your gene identifiers. If your gene identifiers do not match any key type in the species database, you can provide a custom gene set mapping file (gsmap) with the same key type as your query genes for enrichment analysis, and optionally provide a custom background gene list (a text file with one gene ID per line, with the same key type as your query genes) to filter for the genes to be included in the enrichment analysis.")
} else if (top_key_type == 0 & !is.null (gsmap)) {
  message("None of the query genes match any key type in the species database, but a custom gene set mapping file (gsmap) is provided. Enrichment analysis will be performed based on the gsmap file instead of the species database. ")
  key_type <- NULL
} else {
  message("Use key type ", names(top_key_type), " (with ", top_key_type, " out of ", length(goi_tb[[1]]), " genes of interest matches).")
  goi_tb <- goi_tb[goi_tb[[1]] %in% keys(sp_db, keytype = names(top_key_type)),] # filter query genes to keep only those that match the top key type in the species database
  key_type <- names(top_key_type)
}

# select background genes (all genes in the species database with the same key type as the query genes)
if (bg == "protein" & !is.null(key_type)) {
  if ("REFSEQ" %in% columns(sp_db)) {
    bg_genes <- grep("[N|X]P_", keys(sp_db, keytype = "REFSEQ"), value = T) # assuming protein-coding genes have REFSEQ IDs starting with "NP_" or "XP_")
    bg_genes <- AnnotationDbi::select(sp_db, keys = bg_genes, columns = c("REFSEQ", key_type), keytype = "REFSEQ") %>% as.data.frame() %>% dplyr::select(all_of(key_type)) %>% unlist() %>% unique() %>% na.omit()
    message("Using ", length(bg_genes), " protein-coding genes with a REFSEQ prefix of NP_/XP_ as background.")
  } else {
    bg_genes <- keys(sp_db, keytype = key_type)
    warning("REFSEQ key type not found in the species database. Protein-coding genes cannot be extracted. Using all ", length(bg_genes), " genes with key type ", key_type, " as background instead.")
  }
} else if (bg == "all" & !is.null(key_type)) {
  bg_genes <- keys(sp_db, keytype = key_type)
  message("Using all ", length(bg_genes), " genes with key type ", key_type, " as background.")
} else if (file.exists(bg)) {
  bg_genes <- unique(trimws(readLines(bg)))
  bg_genes_raw_count <- length(bg_genes)
  if (!is.null(key_type)) {
    bg_genes <- bg_genes[bg_genes %in% keys(sp_db, keytype = key_type)]
    message("Using ", length(bg_genes), " out of ", bg_genes_raw_count, " custom background genes that matches ", key_type, " IDs.")
  } else {
    message("Using ", length(bg_genes), " custom background genes from the provided file.")
  }
} else {
  bg_genes <- NULL
}

# make sure ranked genes of interest are in the background gene set
if (!is.null(bg_genes)) {
  if (sum(unique(goi_tb[[1]]) %in% bg_genes) == 0) {
    stop("None of the query genes match the background gene set. Please check your query gene list or background gene set.")
  } else {
    if (sum(unique(goi_tb[[1]]) %in% bg_genes)/length(unique(goi_tb[[1]])) < 0.5) {
      warning("Less than 50% of the query genes are in the background gene set. This may lead to unreliable enrichment results. Please check your query gene list and background gene set to make sure they are compatible for enrichment analysis.")
    } else {
      message(sum(unique(goi_tb[[1]]) %in% bg_genes), " out of ", length(unique(goi_tb[[1]])), " query genes are in the background gene set and will be included in the enrichment analysis.")
    }
    goi_tb <- goi_tb[goi_tb[[1]] %in% bg_genes,] 
  }
}

gsea <- function(ranked_query_gene_list, species_db = NULL, gsmap = NULL, key_type = NULL, gene_set = NULL, pvalue_cutoff = 0.05, qvalue_cutoff = 0.1, minGSSize = 5, maxGSSize = 10000, out_file = NULL) {
  if (!is.null(gsmap)) {
    if (!is.null(gene_set)) {
      gsmap <- gsmap[gsmap[[3]] == gene_set, c(1,2)] # filter gsmap to keep only gene sets in the specified gene set category (e.g., "BP") if gene_set is provided
    }
    ego <- GSEA(geneList = ranked_query_gene_list,
                TERM2GENE = gsmap, 
                exponent = 1, # weight of each step in the ranked gene list, default is 1 (no weighting), can be adjusted based on the expected distribution of true positives in the ranked gene list
                pAdjustMethod = "BH",
                pvalueCutoff = pvalue_cutoff,
                minGSSize = minGSSize,
                maxGSSize = maxGSSize,
                verbose = F,
                by = "fgsea") # use fgsea algorithm for GSEA, which is faster
  } else if (!is.null(species_db) & !is.null(key_type) & !is.null(gene_set)) {
    ego <- gseGO(geneList = ranked_query_gene_list,
                 OrgDb = species_db,
                 keyType = key_type,
                 ont = gene_set,
                 exponent = 1, # weight of each step in the ranked gene list, default is 1 (no weighting), can be adjusted based on the expected distribution of true positives in the ranked gene list
                 pAdjustMethod = "BH",
                 pvalueCutoff = pvalue_cutoff,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize, 
                 verbose = F, 
                 by = "fgsea") # use fgsea algorithm for GSEA, which is faster than the default GSEA algorithm in clusterProfiler, especially for large gene sets and ranked gene lists
  } else {
    stop("Please provide either a species database (OrgDb) or a gene set mapping file (gsmap). When using a species database, please also specify the key type and gene set (ontology) for enrichment analysis.")
  }
  
  # filter results based on FDR cutoff (adjusted p-value)
  ego <- dplyr::filter(ego, qvalue < qvalue_cutoff)
  
  # Save results to Excel
  if (!is.null(ego) && nrow(ego) > 0) {
    if (!is.null(out_file)) {
      writexl::write_xlsx(as.data.frame(ego), path = out_file)
      message("GSEA results saved to ", out_file)
    } else {
      message(nrow(ego), " significant enrichments found.")
      print(ego)
    }
  } else {
    message("No significant enrichment found.")
  }
  return(ego)
}

# convert the gene query table to a ranked gene list for GSEA, with gene IDs as names and the ranking metric (e.g., log2 fold change) as values, sorted in decreasing order
goi_list <- goi_tb[[2]] %>% setNames(goi_tb[[1]]) %>% sort(decreasing = TRUE)

# Perform enrichment analysis
if (!is.null(gsmap)) {
  gsmap <- read.delim(gsmap, header = T, stringsAsFactors = F)
  if (!is.null(bg_genes)) {
    if (sum(unique(gsmap[[2]]) %in% bg_genes) == 0) {
      stop("None of the genes in the gene set mapping file (gsmap) match the background gene set. Please check your gsmap file or background gene set.")
    } else {
      message(sum(unique(gsmap[[2]]) %in% bg_genes), " out of ", length(unique(gsmap[[2]])), " genes in the gene set mapping file (gsmap) are in the background gene set and will be included in the enrichment analysis.")
      gsmap <- gsmap[gsmap[[2]] %in% bg_genes, ] # make sure genes in gsmap are in the background gene set
    }
  } else if (sum(unique(gsmap[[2]]) %in% names(goi_list)) == 0) {
    stop("None of the genes in the gene set mapping file (gsmap) match the query genes. Please check your gsmap file or query gene list.")
  }
  
  gsea_res <- gsea(goi_list, gsmap = gsmap, gene_set = gsoi, 
                   pvalue_cutoff = pvalue_cutoff, qvalue_cutoff = qvalue_cutoff, minGSSize = min_gs_size, maxGSSize = max_gs_size,
                   out_file = paste0(paste0(c(out_prefix, gsoi, "GSEA"), collapse = "_"), ".xlsx"))
} else {
  gsea_res <- gsea(goi_list, species_db = sp_db, key_type = key_type, gene_set = gsoi, 
                   pvalue_cutoff = pvalue_cutoff, qvalue_cutoff = qvalue_cutoff, minGSSize = min_gs_size, maxGSSize = max_gs_size,
                   out_file = paste0(paste0(c(out_prefix, gsoi, "GSEA"), collapse = "_"), ".xlsx"))
}

