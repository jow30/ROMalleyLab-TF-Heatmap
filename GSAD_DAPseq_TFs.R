############################################################
# This script aims to visualize the potential master genes or gene modules involved in a bunch of pathways (like DEG-enriched biological processes). 
# In TF analysis, gene modules regulated by the a TF family can be visualized along with binding site conserve scores.
############################################################

# library(httr)
# url <- "https://systemsbiology.cau.edu.cn/PlantGSEAv2/download_fun.php"
# res <- GET(
#   url,
#   query = list(
#     download = "download",
#     fname = "Ara_ALL.txt"
#   ),
#   add_headers(
#     "User-Agent" = "Mozilla/5.0",
#     "Referer" = "https://systemsbiology.cau.edu.cn/PlantGSEAv2/download.php"
#   )
# )
# writeBin(content(res, "raw"), "Ara_ALL.txt")

### preprocess TF annotation file
# gs_anno <- read.table("~/Documents/Projects/Ronan/latest_TF_targets/ath-258-tf-info_simple.txt", header = T)
# gs_anno$tf_name_id <- gs_anno$gene_id
# gs_anno$tf_name_id[gs_anno$gene_id != gs_anno$tf_name] <- paste0(gs_anno$tf_name, "(", gs_anno$gene_id, ")")[gs_anno$gene_id != gs_anno$tf_name]
# gs_anno$anno_name <- apply(gs_anno[,c("tf_name_id","tf_family","tf_clade_curated")] %>% collapse, 1, paste, collapse = "-")
# write.table(gs_anno[, c("gene_id", "anno_name", "tf_family")], file = "ath-258-tf-info_simple_pp.txt", quote = F, sep = "\t", row.names = F, col.names = T)

### logFC info for each TF as additional annotation
# gs_anno <- read.table("ath-258-tf-info_simple_pp.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
# merge<-read.csv("../DEGs/sigDEGs_ABAvsControl_merge.csv")
# column_to_rownames(merge, var = "X")->merge
# merge[,grep("logFC",colnames(merge))]->merge_logFC
# write.table(merge_logFC, file = "logFC.txt", quote = F, sep = "\t", row.names = T, col.names = T)

### preprocess gene annotation file
# df <- read.table("../results/star_salmon/salmon.merged.tx2gene.tsv", header = T, sep = "\t", quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
# df <- data.frame("gene_id"=gsub(".Araport11.447", "", df$gene_id), "gene_name"=gsub(".Araport11.447", "", df$gene_name))
# write.table(df, file = "gene_id_name_mapping.txt", quote = F, sep = "\t", row.names = F, col.names = T)

### preprocess gene expression file
# tmp <- read.table("../results/star_salmon/salmon.merged.gene_tpm.tsv", header = T, sep = "\t", quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
# tmp$gene_id <- gsub(".Araport11.447", "", tmp$gene_id)
# write.table(tmp[,-2], file = "genes_tpm.txt", quote = F, sep = "\t", row.names = F, col.names = T)

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
library(tidyr)

working_dir <- "/Volumes/project/gzy8899/qiaoshan/bulkRNAseq/experiments/Hannah/GSAD"
setwd(working_dir)

### user inputs
subset_cat <- NULL # subset db by the Category column (optional), e.g. subset_cat <- c("BP","TFT")
geneset2plot <- "GSAD_results_all.txt" # user input selected gene set
goi <- "selected_DEGs.txt" # user input gene query (DEGs or genes of interest)
out_file <- "heatmap.pdf"
occur_cutoff <- 1 # min occurance of genes in the selected gene sets to be shown in the heatmap; if 0, show all genes in the selected gene sets (TRUE); if 1, show only the genes of interest that are in the selected gene sets; if >1, only genes that are present in at least cutoff number of gene sets will be shown
diet <- F # whether to reduce the size of the heatmap by randomly subseting genes in the selected gene sets; if TRUE, only 100 genes will be shown in the heatmap
genes_anno <- "gene_id_name_mapping.txt" # the first column should be gene name (matching the gene names in db), the second column should be informative names to display; the third column (optional) should be the group to display, if NULL, no annotation will be added; if not NULL, the annotation will be added as row names in the heatmap; if the third column is provided, the annotation will be colored by group
genes_profile <- "genes_tpm.txt" # profiles for genes, e.g. logFC of genes or expression levels of genes in each sample; the first column should be gene name (matching the gene names in db), the other columns should be the value to display; if NULL, no additional annotation will be added


### TF analysis related inputs
db_file <- "~/Documents/Projects/Ronan/latest_TF_targets/TF_targets_At_cscore.txt"
score2plot <- "n_cons_species_minfrac0" # column name in db to use for coloring; if NULL, will use binary presence/absence (0/1) for coloring
gs_anno <- "ath-258-tf-info_simple_pp.txt" # the first column should be gene set name (matching the gene set names in db), the second column should be the informative names to display; the third column (optional) should be the group to display, if NULL, no annotation will be added
gs_profile <- "logFC.txt" # profiles for gene sets, e.g. logFC of TFs or expression levels of TFs in each sample; the first column should be gene set name (matching the gene set names in db), the other columns should be the value to display; if NULL, no additional annotation will be added

### colors (not intended to be changed by users)
color_bg <- "grey90"
color_gs <- "orange"
color_hl <- "red2"
color_mid <- "yellow"

colorset <- c("darkturquoise", "#E41A1C","gold","#BEBADA","#A6D854","lightpink","mediumslateblue","darkseagreen1", 
              "deeppink4","#FC8D62","cornflowerblue","#A65628","darkkhaki","#A6CEE3","#33A02C","navajowhite", 
              "antiquewhite4","darkgreen", "#FF7F00","midnightblue","#6A3D9A", 
              "dodgerblue4", "violetred1", "darkorchid1","darkslategrey","darkolivegreen",
              "coral4","red3","#FFED6F", "burlywood4", "thistle","chartreuse3","#FB9A99", "darkslateblue", "deeppink", 
              "antiquewhite",  "chartreuse4","chocolate3", "violet","mediumblue","turquoise4",
              "springgreen3", "cyan2", "thistle2", "#8DD3C7", "mediumpurple1", "khaki1","#FDB462","hotpink4",
              "dodgerblue2", "deepskyblue1",  "darkslategray4", "lemonchiffon1","#E78AC3",
              "deeppink3","#FFFFB3","burlywood", "firebrick1",  "#CAB2D6","darkolivegreen4", "springgreen2",
              "forestgreen","skyblue","#999999","#F781BF", "#80B1D3","yellow1", "olivedrab1", 
              "magenta1",  "blue", "azure2", "#B3B3B3", "#BC80BD","darkorchid4", 
              "#8DA0CB","#377EB8","purple", "antiquewhite3", 
              "#1F78B4","firebrick", "#FCCDE5","gold3","palegreen", "azure4",
              "azure3", "#B2DF8A","cadetblue1", "darkviolet", 
              "darksalmon", "chartreuse1","#E5C494","#4DAF4A",      "lightslateblue", 
              "tomato3", "#D9D9D9", "#CCEBC5","slateblue4", 
              "green2", "royalblue1", "#66C2A5","lightblue1", "mediumorchid"
)

# load gene set database
db <- read.table(db_file, header = T, sep = "\t", quote = "", stringsAsFactors = F)

if (!is.null(subset_cat)) {
  db <- db[db$Category %in% subset_cat, ]
}

geneset2plot <- readLines(geneset2plot)
geneset2plot <- unique(geneset2plot)

goi <- readLines(goi)

# check whether all geneset2plot are in db
if (!all(geneset2plot %in% db$tf)) {
  stop("Some user-selected gene sets are not in the database (or the selected category).")
}

gsg <- db$gene[db$tf %in% geneset2plot]
gsg <- unique(gsg)
goi_in_gsg <- goi[goi %in% gsg]
warning(paste0(length(goi_in_gsg), " out of ", length(goi), " DEGs are in the selected gene sets."))

if (occur_cutoff == 0) {
  genes2plot <- gsg
} else {
  genes2plot <- goi_in_gsg
}

# annotate genes if genes_anno is not NULL
if (!is.null(genes_anno)) {
  genes_anno <- read.table(genes_anno, header = T, sep = "\t", quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
  if (ncol(genes_anno) < 2) {
    stop("The annotation file should have at least two columns: the first column is geneID, the second column is annotation (gene name) to display.")
  } else {
    genes_anno[[2]] <- trimws(as.character(genes_anno[[2]]))
    genes_anno[[2]][is.na(genes_anno[[2]]) | genes_anno[[2]] == "" | tolower(genes_anno[[2]]) == "null"] <- genes_anno[[1]][is.na(genes_anno[[2]]) | genes_anno[[2]] == "" | tolower(genes_anno[[2]]) == "null"]
  }
  if (ncol(genes_anno) < 3) {
    gene_grp <- NULL
  } else {
    genes_anno[[3]] <- trimws(as.character(genes_anno[[3]]))
    genes_anno[[3]][is.na(genes_anno[[3]]) | genes_anno[[3]] == "" | tolower(genes_anno[[3]]) == "null"] <- "grp_unavail"
  }
  genes_anno <- genes_anno[genes_anno[[1]] %in% genes2plot, ]
  genes_anno <- genes_anno[match(genes2plot, genes_anno[[1]]), , drop = F]
  
  gene_labels <- genes_anno[[2]] %>% setNames(genes2plot)
  if (ncol(genes_anno) >= 3) {
    gene_grp <- factor(genes_anno[[3]], levels = unique(genes_anno[[3]]))
    names(gene_grp) <- genes_anno[[1]]
  }
} else {
  gene_labels <- genes2plot %>% setNames(genes2plot)
  gene_grp <- NULL
}

# annotate genesets if gs_anno is not NULL
if (!is.null(gs_anno)) {
  gs_anno <- read.table(gs_anno, header = T, sep = "\t", quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
  if (ncol(gs_anno) < 2) {
    stop("The annotation file should have at least two columns: the first column is gene set name, the second column is annotation to display.")
  } else {
    gs_anno[[2]] <- trimws(as.character(gs_anno[[2]]))
    gs_anno[[2]][is.na(gs_anno[[2]]) | gs_anno[[2]] == "" | tolower(gs_anno[[2]]) == "null"] <- "anno_unavail"
  }
  if (ncol(gs_anno) < 3) {
    gs_anno[[3]] <- "grp_unavail"
  } else {
    gs_anno[[3]] <- trimws(as.character(gs_anno[[3]]))
    gs_anno[[3]][is.na(gs_anno[[3]]) | gs_anno[[3]] == "" | tolower(gs_anno[[3]]) == "null"] <- "grp_unavail"
  }
  gs_anno <- gs_anno[gs_anno[[1]] %in% geneset2plot, ]
  gs_anno <- gs_anno[match(geneset2plot, gs_anno[[1]]), , drop = F]
  
  gs_labels <- gs_anno[[2]] %>% setNames(geneset2plot)

  gs_grp <- factor(gs_anno[[3]], levels = unique(gs_anno[[3]]))
  names(gs_grp) <- gs_anno[[1]]
} else {
  gs_labels <- geneset2plot %>% setNames(geneset2plot)
  gs_grp <- NULL
}

if (!is.null(gs_profile)) {
  gs_profile <- read.table(gs_profile, header = T, sep = "\t", row.names = 1, quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
  if (ncol(gs_profile) < 2) {
    stop("The additional annotation file should have at least two columns: the first column is gene set name, the other columns are annotation to display.")
  }
  gs_profile <- gs_profile[rownames(gs_profile) %in% geneset2plot, ]
  gs_profile <- gs_profile[match(geneset2plot, rownames(gs_profile)), , drop = F]
  rownames(gs_profile) <- geneset2plot
}

if (!is.null(genes_profile)) {
  genes_profile <- read.table(genes_profile, header = T, sep = "\t", row.names = 1, quote = "", stringsAsFactors = F, fill = T, na.strings = c("", "NA", "null", "NULL", "Null"))
  if (ncol(genes_profile) < 2) {
    stop("The additional annotation file should have at least two columns: the first column is gene name, the other columns are annotation to display.")
  }
  genes_profile <- genes_profile[rownames(genes_profile) %in% genes2plot, ]
  genes_profile <- genes_profile[match(genes2plot, rownames(genes_profile)), , drop = F]
  rownames(genes_profile) <- genes2plot
}

# create a matrix for heatmap
mat <- matrix(0, nrow = length(genes2plot), ncol = length(geneset2plot))
rownames(mat) <- genes2plot
colnames(mat) <- geneset2plot

for (gs in geneset2plot) {
  gsg_sel <- db$gene[db$tf == gs]
  if (occur_cutoff == 0) {
    mat[gsg_sel, gs] <- 1
    gsg_sel <- gsg_sel[gsg_sel %in% goi_in_gsg]
    mat[gsg_sel, gs] <- 2
  } else {
    gsg_sel <- gsg_sel[gsg_sel %in% genes2plot]
    mat[gsg_sel, gs] <- 2
  }
}

# add row annotation to show the number of gene sets each gene belongs to
gs_count <- rowSums(mat > 0)

if (occur_cutoff > 1) {
  genes2plot <- genes2plot[gs_count >= occur_cutoff]
  mat <- mat[genes2plot, , drop = FALSE]
  gs_count <- gs_count[genes2plot]
}

if (diet && nrow(mat) > 100) {
  genes2plot <- sample(genes2plot, 100)
  mat <- mat[genes2plot, , drop = FALSE]
  gs_count <- gs_count[genes2plot]
}

row_anno <- rowAnnotation(GeneSetCount = gs_count[rownames(mat)],
                          col = list(GeneSetCount = colorRamp2(c(min(gs_count), max(gs_count)), c("white", "blue"))),
                          show_legend = T, annotation_legend_param = list(title = "Gene Set Count"))
if (!is.null(gs_grp)) {
  top_anno <- HeatmapAnnotation(`Family/Group` = gs_grp[colnames(mat)],
                                col = list(`Family/Group` = structure(colorset[1:length(unique(gs_grp))], names = levels(gs_grp))),
                                show_legend = T, annotation_legend_param = list(title = "Family/Group"))
} else {
  top_anno <- NULL
}

if (!is.null(gs_profile)) {
  anno_df <- gs_profile[colnames(mat), , drop = FALSE]
  col_fun <- colorRamp2(c(min(anno_df, na.rm = TRUE), 0, max(anno_df, na.rm = TRUE)), c("darkgreen", "white", "red3"))
  col_list <- setNames(
    rep(list(col_fun), ncol(anno_df)),
    colnames(anno_df)
  )
  show_legend <- c(TRUE, rep(FALSE, ncol(anno_df) - 1))
  names(show_legend) <- colnames(anno_df)
  bottom_anno <- HeatmapAnnotation(df = anno_df,
                                   col = col_list,
                                   show_legend = show_legend, 
                                   annotation_legend_param = list(title = "Gene Set Profile"))
  profile_ncol <- ncol(gs_profile)
} else {
  bottom_anno <- NULL
  profile_ncol <- 0
}

if (!is.null(genes_profile)) {
  anno_df <- genes_profile[rownames(mat), , drop = FALSE]
  col_fun <- colorRamp2(c(min(anno_df, na.rm = TRUE), (min(anno_df, na.rm = TRUE)+max(anno_df, na.rm = TRUE))/2, max(anno_df, na.rm = TRUE)), c("darkgreen", "white", "red3"))
  col_list <- setNames(
    rep(list(col_fun), ncol(anno_df)),
    colnames(anno_df)
  )
  show_legend <- c(TRUE, rep(FALSE, ncol(anno_df) - 1))
  names(show_legend) <- colnames(anno_df)
  
  if (!is.null(gene_grp)) {
    right_anno <- rowAnnotation(`Gene Group` = gene_grp[rownames(mat)],
                                df = anno_df,
                                col = c(list(`Gene Group` = structure(colorset[1:length(unique(gene_grp))], names = levels(gene_grp))), col_list),
                                show_legend = c(`Gene Group` = TRUE, show_legend),
                                annotation_legend_param = list(title = "Gene Group"))
    profile_ncol_g <- ncol(genes_profile)+1
  } else {
    right_anno <- rowAnnotation(df = anno_df,
                                col = col_list,
                                show_legend = show_legend, 
                                annotation_legend_param = list(title = "Gene Profile"))
    profile_ncol_g <- ncol(genes_profile)
  }
} else if (is.null(genes_profile) & !is.null(gene_grp)) {
  right_anno <- rowAnnotation(`Gene Group` = gene_grp[rownames(mat)],
                              col = list(`Gene Group` = structure(colorset[1:length(unique(gene_grp))], names = levels(gene_grp))),
                              show_legend = T, annotation_legend_param = list(title = "Gene Group"))
  profile_ncol_g <- 1
} else {
  right_anno <- NULL
  profile_ncol_g <- 0
}

if(prod(dim(mat))>100000) {
  warning("The heatmap contains more than 100,000 cells, which may take a long time to render. Consider reducing the number of gene sets or genes to plot.")
}

cell_size <- unit(3, "mm")
pdf(out_file, width = 8 + ncol(mat) * 0.08 + profile_ncol_g * 0.2, height = 10 + nrow(mat) * 0.08 + profile_ncol * 0.2)
ht <- Heatmap(mat, name = "occurance",
              col = c("0" = color_bg, "1" = color_gs, "2" = color_hl),
              width  = ncol(mat) * cell_size,
              height = nrow(mat) * cell_size,
              show_row_names = T, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              row_labels = gene_labels[rownames(mat)],
              column_labels = gs_labels[colnames(mat)],
              column_title_rot = 45,
              column_title_gp = gpar(fontsize = 9, fontface = "bold"),
              rect_gp = gpar(col = "white", lwd = 1),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              show_heatmap_legend = F,
              left_annotation = row_anno,
              top_annotation = top_anno,
              bottom_annotation = bottom_anno,
              right_annotation = right_anno)
draw(ht)
dev.off()
if (!is.null(gs_grp)) {
  pdf(gsub(".pdf", "_separate_grp.pdf", out_file), width = 8 + ncol(mat) * 0.08 + profile_ncol_g * 0.2, height = 10 + nrow(mat) * 0.08 + profile_ncol * 0.2)
  ht <- Heatmap(mat, name = "occurance",
                col = c("0" = color_bg, "1" = color_gs, "2" = color_hl),
                width  = ncol(mat) * cell_size,
                height = nrow(mat) * cell_size,
                show_row_names = T, show_column_names = T,
                cluster_rows = T, cluster_columns = T,
                column_split = gs_grp[colnames(mat)],
                cluster_column_slices = F,
                row_labels = gene_labels[rownames(mat)],
                column_labels = gs_labels[colnames(mat)],
                column_title_rot = 45,
                column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                rect_gp = gpar(col = "white", lwd = 1),
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                show_heatmap_legend = F,
                left_annotation = row_anno,
                top_annotation = top_anno,
                bottom_annotation = bottom_anno,
                right_annotation = right_anno)
  draw(ht)
  dev.off()
}

if (!is.null(score2plot)) {
  # mat for scores
  mat <- db %>%
    dplyr::select(gene, tf, all_of(score2plot)) %>%
    pivot_wider(
      names_from = tf, #TF to columns
      values_from = all_of(score2plot), # genes to rows
      values_fill = 0
    ) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  mat <- mat[genes2plot, geneset2plot]
  
  if(prod(dim(mat))>100000) {
    warning("The heatmap contains more than 100,000 cells, which may take a long time to render. Consider reducing the number of gene sets or genes to plot.")
  }
  
  cell_size <- unit(3, "mm")
  pdf(gsub(".pdf", paste0("_", score2plot, ".pdf"), out_file), width = 8 + ncol(mat) * 0.08 + profile_ncol_g * 0.2, height = 10 + nrow(mat) * 0.08 + profile_ncol * 0.2)
  ht <- Heatmap(mat, name = score2plot,
                col = colorRamp2(
                  c(0, max(mat)/2, max(mat)),
                  c(color_bg, color_mid, color_hl)
                ),
                width  = ncol(mat) * cell_size,
                height = nrow(mat) * cell_size,
                show_row_names = T, show_column_names = T,
                cluster_rows = T, cluster_columns = T,
                row_labels = gene_labels[rownames(mat)],
                column_labels = gs_labels[colnames(mat)],
                column_title_rot = 45,
                column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                rect_gp = gpar(col = "white", lwd = 1),
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                show_heatmap_legend = T,
                left_annotation = row_anno,
                top_annotation = top_anno,
                bottom_annotation = bottom_anno,
                right_annotation = right_anno)
  draw(ht)
  dev.off()
  if (!is.null(gs_grp)) {
    pdf(gsub(".pdf", paste0("_", score2plot, "_separate_grp.pdf"), out_file), width = 8 + ncol(mat) * 0.08 + profile_ncol_g * 0.2, height = 10 + nrow(mat) * 0.08 + profile_ncol * 0.2)
    ht <- Heatmap(mat, name = score2plot,
                  col = colorRamp2(
                    c(0, max(mat)/2, max(mat)),
                    c(color_bg, color_mid, color_hl)
                  ),
                  width  = ncol(mat) * cell_size,
                  height = nrow(mat) * cell_size,
                  show_row_names = T, show_column_names = T,
                  cluster_rows = T, cluster_columns = T,
                  column_split = gs_grp[colnames(mat)],
                  cluster_column_slices = F,
                  row_labels = gene_labels[rownames(mat)],
                  column_labels = gs_labels[colnames(mat)],
                  column_title_rot = 45,
                  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
                  rect_gp = gpar(col = "white", lwd = 1),
                  row_names_gp = gpar(fontsize = 8),
                  column_names_gp = gpar(fontsize = 8),
                  show_heatmap_legend = T,
                  left_annotation = row_anno,
                  top_annotation = top_anno,
                  bottom_annotation = bottom_anno,
                  right_annotation = right_anno)
    draw(ht)
    dev.off()
  }
} 
