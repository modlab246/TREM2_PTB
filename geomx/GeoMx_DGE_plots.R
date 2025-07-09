
# Generate gene sets for contrasts and intersections as well as the

setwd("/Users/jonathanperrie/Documents/UCLA/projects/geomx")

library(DESeq2)
library(dplyr)
library(ggrepel)
library(openxlsx)
library(sva)
library(circlize)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(ggsignif)
library(umap)
library(eulerr)
library(readxl)

# raw geomx counts 
raw_count_input <- "data/geomx_ptb_raw_data.csv"
select_genes_input <- "data/geomx_select_genes_volc.xlsx"
mat <- read.table(raw_count_input,sep=",",row.names=1,header=T,check.names=F)

SampleAnnotationFile <- "dcc/geomx_ptb_annotations.xlsx"
sample_annotation <- read_excel(SampleAnnotationFile)

# Comprehensive row names: region, sample, ROI
batch_names1 <- sample_annotation[,c("Sub-classification","Scan name","ROI")]
batch_names1$ROI <- sprintf("%02d",batch_names1$ROI)
batch_names1$`Scan name` <- gsub("PTB","PTB ",batch_names1$`Scan name`)
batch_names1 <- apply(batch_names1,1,function(x) paste(x, collapse =" | "))

batch_names2 <- sample_annotation[,c("Segment tags","Scan name","ROI")]
batch_names2$ROI <- sprintf("%02d",batch_names2$ROI)
batch_names2$`Scan name` <- gsub("PTB","PTB ",batch_names2$`Scan name`)
batch_names2 <- apply(batch_names2,1,function(x) paste(x, collapse =" | "))

# may remove this in future but at least consistent 
roi_19_bool <- grepl("\\| PTB 22.2 \\| 19",colnames(mat))
ptb_bool <- grepl("\\| PTB ",colnames(mat)) & !roi_19_bool

mat <- mat[,ptb_bool]
# Ensure that batch_names1 and batch_names2 are properly ordered
batch_names1 <- batch_names1[match(colnames(mat), batch_names2)]

treatment <- trimws(sapply(strsplit(batch_names1,"\\|"),"[",1))
batch <- trimws(sapply(strsplit(batch_names1,"\\|"),function(x) paste(x[2:(length(x)-1)], collapse = " ")))

meta <- data.frame(treatment=trimws(treatment),batch=trimws(batch))
rownames(meta) <- colnames(mat)
meta$treatment[meta$treatment == "Giant"] = "Core"
meta$treatment <- factor(meta$treatment, levels=c("Alveoli","Core","Mantle", "Infiltrate","Fibrosis"))
meta$batch <- as.factor(meta$batch)

# treatment (region) and batch (patient) as additive effects modelling expression 
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = meta,
                              design = ~ treatment+batch)
dds <- DESeq(dds) 

# we only care about alveoli vs core vs mantle (ignore infiltrate and fibrosis)
res_mc <- results(dds,contrast=c("treatment","Mantle","Core"))
res_ma <- results(dds,contrast=c("treatment","Mantle","Alveoli"))
res_ca <- results(dds,contrast=c("treatment","Core","Alveoli"))

# function for filtering by pval and threshold while sorting by fc (+/-)
process_res <- function(res,thr,fc){
  res<-res[!is.na(res$padj),]
  res<-res[res$padj<=thr,]
  
  neg_res = res[res$log2FoldChange<(-fc),]
  plus_res = res[res$log2FoldChange>fc,]
  
  a<-rownames(neg_res[order(abs(neg_res$log2FoldChange),decreasing=TRUE),])
  b<-rownames(plus_res[order(abs(plus_res$log2FoldChange),decreasing=FALSE),])
  
  res<-rbind(res[a,],res[b,])
  return(res)
}

# volcano plot function 
plot_volc <- function(res, title, select_genes = NULL, top_thr=1E-3, top_fc=2, fname=""){
  # filter out res with nan for adjusted pval 
  res <- filter(data.frame(res), abs(padj) >= 0)
  
  # use top threshold/fc or select genes 
  if (is.null(select_genes)) {
    top_g <- rownames(filter(res, (padj <= top_thr & abs(log2FoldChange) >= top_fc) | abs(log2FoldChange) >= (top_fc + 1)))
  } else {
    top_g <- intersect(rownames(res), select_genes)
  }
  
  # gene markers that are blue colored 
  padj_thr <- 5E-2
  log2fc_thr <- 1
  
  gene_hit <- rownames(filter(res, padj <= padj_thr & abs(log2FoldChange) >= log2fc_thr))
  # default marker is grey
  res$marker = "grey"
  res$marker[unlist(lapply(rownames(res), function(x){ x %in% gene_hit}))] = "blue"
  res$Gene <- rownames(res)
  
  g <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj), label=Gene, size=marker)) +
    geom_point(aes(color=marker), alpha = 0.35, stroke = 0.5) + scale_color_identity() +
    # top genes have a red center 
    geom_point(data = subset(res, Gene %in% top_g), 
               aes(color = "red"), alpha = 1, size = 0.75) +
    scale_size_manual(values=c("grey" = 0.5, "blue" = 3)) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed") +
    geom_vline(xintercept = c(-log2fc_thr, log2fc_thr), linetype = "dashed") +
    # right side text annotations 
    geom_text_repel(data = subset(res, Gene %in% top_g & log2FoldChange > 0),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2,
                    max.overlaps = Inf, nudge_x = 0.5, nudge_y = 0.5,
                    segment.alpha = 0.25, force = 1.5) +
    # left side text annotations 
    geom_text_repel(data = subset(res, Gene %in% top_g & log2FoldChange < 0),
                    size = 4, point.padding = 0.15, color = "black",
                    min.segment.length = .1, box.padding = .2,
                    max.overlaps = Inf, nudge_x = -0.5, nudge_y = 0.5,
                    segment.alpha = 0.25, force = 1.5) +
    ggtitle(title) +
    ylab(expression(paste("-log"[10], "(", "adj P-value", ")"))) +
    xlab("") + 
    annotate(
      "text",
      x = 0, 
      y = -Inf,  
      label = expression(paste("log"[2], "FC")),
      vjust = 2,
      size= 20/2.845276
    ) + 
    coord_cartesian(xlim = c(-6, 6), clip = 'off') +
    scale_x_continuous(breaks = seq(-6, 6, by = 2)) + 
    theme_classic() +
    theme(legend.position="none",
          panel.background=element_blank(),
          text = element_text(size=20),
          plot.margin = margin(2, 3,0, 1, "lines"))
  
  if (endsWith(fname, ".png")){
    ggsave(file = fname, plot = g, device = "png", width=9, height=6, dpi = 400)
  } else if (endsWith(fname, ".pdf")) {
    ggsave(file = fname, plot = g, device = "pdf", width=9, height=6, dpi = 400)
  } else {
    g
  }
}

# read in select genes to highlight in volcano 
select_genes <- read.xlsx(select_genes_input)

# plot volcanoes with offset for titles 
plot_volc(res_mc, paste(sprintf("%*s", 10 + nchar("Core"), "Core"),"Mantle",sep=strrep(" ",43)),
          select_genes = c(select_genes$CvM, select_genes$MvC), top_thr=1E-20, top_fc=1, fname="plots/manuscript/mantle_core_dge.png")
plot_volc(res_ca, paste(sprintf("%*s", 10 + nchar("Alveoli"), "Alveoli"),"Core",sep=strrep(" ",40)),
          select_genes = c(select_genes$CvA, select_genes$AvC), top_thr=1E-30, top_fc=1, fname="plots/manuscript/alveoli_core_dge.png")
plot_volc(res_ma, paste(sprintf("%*s", 10 + nchar("Alveoli"), "Alveoli"),"Mantle",sep=strrep(" ",40)),
          select_genes = c(select_genes$AvM, select_genes$MvA), top_thr=1E-30, top_fc=1, fname="plots/manuscript/mantle_alveoli_dge.png")

# find union set between two contrast lists from deseq2 
# flip_res is used to flip log2fc of some list 
prepare_joint_df <- function(res1, res2, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = FALSE, flip_res2 = FALSE) {
  res1_copy <- data.frame(res1)
  res2_copy <- data.frame(res2)
  
  # flip log2fc 
  if (flip_res1) {
    res1_copy$log2FoldChange <- -res1_copy$log2FoldChange
  }
  if (flip_res2) {
    res2_copy$log2FoldChange <- -res2_copy$log2FoldChange
  }
  
  # filter by fc and pval
  df1 <- filter(res1_copy, log2FoldChange > log2FC_threshold, padj < padj_threshold)
  df2 <- filter(res2_copy, log2FoldChange > log2FC_threshold, padj < padj_threshold)
  
  geneNames <- union(rownames(df1), rownames(df2))
  combined_df <- cbind(res1_copy[geneNames, ], res2_copy[geneNames, ])
  colnames(combined_df) <- c(paste0(colnames(df1), ".1"), paste0(colnames(df2), ".2"))
  
  # joint pval (*) and fc (+)
  combined_df$jointP <- combined_df$padj.1 * combined_df$padj.2
  combined_df$jointFC <- abs(combined_df$log2FoldChange.2) + abs(combined_df$log2FoldChange.1)
  
  combined_df <- na.omit(combined_df)
  combined_df$gene <- rownames(combined_df)
  
  return(combined_df)
}

# find genes that are at least implicated in both sets and significant 
# res_mc means the default order is (+) fc for mantle and (-) fc for core 
core_df <- prepare_joint_df(res_ca, res_mc, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = FALSE, flip_res2 = TRUE)
alv_df <- prepare_joint_df(res_ca, res_ma, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = TRUE, flip_res2 = TRUE)
man_df <- prepare_joint_df(res_mc, res_ma, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = FALSE, flip_res2 = FALSE)

# granuloma genes are in core and mantle 
gran_df <- prepare_joint_df(res_ca, res_ma, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = FALSE, flip_res2 = FALSE)
# macrophage genes are in core and alveoli 
mph_df <- prepare_joint_df(res_mc, res_ma, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = TRUE, flip_res2 = TRUE)
# outzone genes are not in core 
outzone_df <- prepare_joint_df(res_mc, res_ca, log2FC_threshold = 0, padj_threshold = 0.05, flip_res1 = FALSE, flip_res2 = TRUE)

# gene lists from DGE contrasts ordered by fc 
ca_geneset <- rev(rownames(filter(arrange(as.data.frame(process_res(res_ca,5E-2,0)),log2FoldChange),log2FoldChange>0)))
ac_geneset <- rownames(filter(arrange(as.data.frame(process_res(res_ca,5E-2,0)),log2FoldChange),log2FoldChange<0))
mc_geneset <- rev(rownames(filter(arrange(as.data.frame(process_res(res_mc,5E-2,0)),log2FoldChange),log2FoldChange>0)))
cm_geneset <- rownames(filter(arrange(as.data.frame(process_res(res_mc,5E-2,0)),log2FoldChange),log2FoldChange<0))
ma_geneset <- rev(rownames(filter(arrange(as.data.frame(process_res(res_ma,5E-2,0)),log2FoldChange),log2FoldChange>0)))
am_geneset <- rownames(filter(arrange(as.data.frame(process_res(res_ma,5E-2,0)),log2FoldChange),log2FoldChange<0))

# gene lists from joint sets (now taking intersection by forcing positive fc and significance) ordered by fc
core_geneset <- rev(rownames(filter(arrange(core_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))
alv_geneset <- rev(rownames(filter(arrange(alv_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))
man_geneset <- rev(rownames(filter(arrange(man_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))
gran_geneset <- rev(rownames(filter(arrange(gran_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))
mph_geneset <- rev(rownames(filter(arrange(mph_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))
outzone_geneset <- rev(rownames(filter(arrange(outzone_df,jointFC),log2FoldChange.1>0 & log2FoldChange.2>0 & padj.1 < 5E-2 & padj.2 < 5E-2)))

maxlen <- max(
  length(ca_geneset), length(ac_geneset), length(mc_geneset), length(cm_geneset),
  length(ma_geneset), length(am_geneset), length(core_geneset), length(alv_geneset),
  length(man_geneset), length(gran_geneset), length(mph_geneset), length(outzone_geneset)
)

# pad for equal size 
pad_list <- function(x,maxlen=maxlen,rep_val="") {
  c(x, rep(rep_val, maxlen - length(x)))  
}

# DGE results 
gene_df_rbt <- data.frame(
  CvA = pad_list(ca_geneset,maxlen=maxlen),
  AvC = pad_list(ac_geneset,maxlen=maxlen),
  CvM = pad_list(cm_geneset,maxlen=maxlen),
  MvC = pad_list(mc_geneset,maxlen=maxlen),
  MvA = pad_list(ma_geneset,maxlen=maxlen),
  AvM = pad_list(am_geneset,maxlen=maxlen),
  Core = pad_list(core_geneset,maxlen=maxlen),
  Alveoli = pad_list(alv_geneset,maxlen=maxlen),
  Mantle = pad_list(man_geneset,maxlen=maxlen),
  Granuloma = pad_list(gran_geneset,maxlen=maxlen),
  Macrophage = pad_list(mph_geneset,maxlen=maxlen),
  Outzone = pad_list(outzone_geneset,maxlen=maxlen)
)

# function to stitch DGE result with indices ordered by each gene set 
order_gene_pos <- function(res_index, gene_df) {
  res_index <- data.frame(res_index)
  
  # mva, avm, ... core, mantle, ... outzone 
  for (col_name in names(gene_df)) {
    res_index[[col_name]] <- -1
    genes_to_match <- unname(unlist(gene_df[[col_name]]))
    matching_indices <- which(rownames(res_index) %in% genes_to_match)
    # the genes with the largest fc will have the biggest index (sort by largest, decreasing)
    res_index[matching_indices, col_name] <- nrow(res_index) - match(rownames(res_index)[matching_indices], genes_to_match)
  }
  
  return(res_index)
}

# Create a workbook
wb <- createWorkbook()

# Add the first sheet with gene_df_rbt
addWorksheet(wb, "Gene")
writeData(wb, "Gene", gene_df_rbt)

# Add sheets for res_mc, res_ca, and res_ma
addWorksheet(wb, "MvC")
res_mc_index <- order_gene_pos(res_mc,gene_df_rbt)
writeData(wb, "MvC", res_mc_index, rowNames=TRUE)

addWorksheet(wb, "CvA")
res_ca_index <- order_gene_pos(res_ca,gene_df_rbt)
writeData(wb, "CvA", res_ca_index, rowNames=TRUE)

addWorksheet(wb, "MvA")
res_ma_index <- order_gene_pos(res_ma,gene_df_rbt)
writeData(wb, "MvA", res_ma_index, rowNames=TRUE)

# Save the workbook to a file
saveWorkbook(wb, "data/DGE_genelist_RBT.xlsx", overwrite = TRUE)

# heatmap plot 
# we only want ACM ROIs
acm_bool <- meta$treatment %in% c("Core","Mantle","Alveoli")

# scale counts (like log)
rld <- vst(dds[,acm_bool])
rlog_mat <- assay(rld)
zscale_mat <- t(scale(t(rlog_mat)))

# color by region and patient 
ha <- HeatmapAnnotation(
  df = data.frame(
    Region = meta$treatment[acm_bool],
    Batch  = factor(meta$batch[acm_bool],
                    levels = c("PTB 21.1","PTB 22.2","PTB 22.3","PTB 22.1"))
  ),
  col = list(
    Region = c("Alveoli" = "#f8766d",
               "Core"    = "#00ba38",
               "Mantle"  = "#619cff"),
    Batch  = c("PTB 21.1" = "#635474",
               "PTB 22.2" = "#0974D8",
               "PTB 22.3" = "#FDF490",
               "PTB 22.1" = "#790100")
  ),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16)
)

micro_genes <- c("CAMP","CCL1","CCL11","CCL13","CCL17","CCL18","CCL19","CCL20",
                 "CCL21","CCL22","CCL24","CCL25","CCL26","CCL27","CCL28","CCL8",
                 "CD40LG","CXCL1","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14",
                 "CXCL2","CXCL3","CXCL6","CXCL9","DEFB1","GBP1","GNLY","GZMB","HAMP",
                 "IDO1","IFNG","IL1B","IL26","IL32","KLK5","KRT6A","LYZ","NPC2","P2RX7",
                 "PLA2G2A","PRF1","S100A12","S100A7","S100A7A","S100A8","S100A9","TXN",
                 "VDR","XCL1")
flat_select_genes <- micro_genes[micro_genes %in% unname(unlist(gene_df_rbt))]



hm_df <- zscale_mat[flat_select_genes, ]
colnames(hm_df) <- sub(
  ".*PTB\\s*([0-9\\.]+)\\s*\\|\\s*([0-9]+).*",
  "\\1-\\2",
  colnames(hm_df)
)

ht_exp <- Heatmap(
  hm_df, 
  name = "Expression",
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  cluster_rows = TRUE, 
  cluster_columns = TRUE, 
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  top_annotation = ha,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = FALSE
)

# we grab the order and clusters from a first call to recreate a row annotation 
ht_obj <- draw(ht_exp)
row_dend <- row_dend(ht_obj)
row_clusters <- cutree(as.hclust(row_dend), k = 4)
row_cluster_factor <- factor(row_clusters)
cluster_colors <- c("green","red","cyan","blue") 
names(cluster_colors) <- levels(row_cluster_factor)

ra <- rowAnnotation(
  Cluster = row_cluster_factor,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16),
  show_legend = FALSE
)

ht_exp <- Heatmap(
  hm_df, 
  name = "Expression",
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  cluster_rows = TRUE, 
  cluster_columns = TRUE, 
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  top_annotation = ha,
  right_annotation = ra,
  column_names_gp = gpar(fontsize = 9),
  row_names_gp = gpar(fontsize = 14),
  show_heatmap_legend = FALSE
)

# saving heatmap without legend 
output_file <- "plots/manuscript/hm_antimicrobial_nolegend_boxes.png"
png(output_file, width = 1800, height = 1800, res = 200)
draw(ht_exp)

highlight_block <- function(ht_name, col_start, col_end, row_start, row_end, border_col = "grey", border_lwd = 2, fill = NA) {
  decorate_heatmap_body(ht_name, {
    n_col <- ncol(hm_df)         
    n_row <- nrow(hm_df)         
    
    # find length of box to be drawn 
    width  <- (col_end - col_start) / n_col
    height <- (row_end - row_start) / n_row
    
    # coords are at center of box
    x_center <- (col_start + col_end) / (2 * n_col)
    y_center <- ((n_row - row_start) + (n_row - row_end)) / (2 * n_row)
    
    # draw rectangle 
    grid.rect(
      x = x_center, y = y_center,
      width  = width,  height = height,
      gp = gpar(col = border_col, lwd = border_lwd, fill = fill)
    )
  })
}

# antimicrobial, 4 clusters, robert does not like color coordinating 
# highlight_block("Expression", 0.05, 24.95, 0.05, 9.95, border_col = scales::alpha(cluster_colors[1], 0.8), border_lwd = 3.5)
# highlight_block("Expression", 0.05, 38.95, 10.05, 16.95, border_col = scales::alpha(cluster_colors[3], 0.8), border_lwd = 3.5)
# highlight_block("Expression", 25.05, 38.95, 17.05, 19.95, border_col = scales::alpha(cluster_colors[4], 0.8), border_lwd = 3.5)
# highlight_block("Expression", 39.05, 58.95, 20.05, 25.95, border_col = scales::alpha(cluster_colors[2], 0.8), border_lwd = 3.5)

highlight_block("Expression", 0.05, 24.95, 0.05, 9.95, border_col = scales::alpha("red", 0.8), border_lwd = 3.5)
highlight_block("Expression", 0.05, 38.95, 10.05, 16.95, border_col = scales::alpha("red", 0.8), border_lwd = 3.5)
highlight_block("Expression", 25.05, 38.95, 17.05, 19.95, border_col = scales::alpha("red", 0.8), border_lwd = 3.5)
highlight_block("Expression", 39.05, 58.95, 20.05, 25.95, border_col = scales::alpha("red", 0.8), border_lwd = 3.5)

dev.off()

# re-do clustering with legend, so we can strip just that element and place properly 
ha <- HeatmapAnnotation(
  df = data.frame(
    Region = meta$treatment[acm_bool],
    Batch  = factor(meta$batch[acm_bool],
                    levels = c("PTB 21.1","PTB 22.2","PTB 22.3","PTB 22.1"))
  ),
  col = list(
    Region = c("Alveoli" = "#f8766d",
               "Core"    = "#00ba38",
               "Mantle"  = "#619cff"),
    Batch  = c("PTB 21.1" = "#635474",
               "PTB 22.2" = "#0974D8",
               "PTB 22.3" = "#FDF490",
               "PTB 22.1" = "#790100")
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_position = "topcenter",    # put the title above and centered
    title_gp       = gpar(fontsize = 16),
    labels_gp      = gpar(fontsize = 12)
  )
)

ra <- rowAnnotation(
  Cluster = row_cluster_factor,
  col = list(Cluster = cluster_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 16),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_position = "topcenter",   
    title_gp       = gpar(fontsize = 16),
    labels_gp      = gpar(fontsize = 12),
    padding = unit(c(5, 80, 5, 10), "mm")
  )
)


ht_exp <- Heatmap(
  hm_df, 
  name = "Expression",
  col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
  cluster_rows = TRUE, 
  cluster_columns = TRUE, 
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  top_annotation = ha,
  right_annotation = ra,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Expression",               # optional override
    title_position = "topcenter",       # <-- center the title
    title_gp = gpar(fontsize = 16),
    labels_gp = gpar(fontsize = 10)
  )
)

# Define the output file path
output_file <- "plots/manuscript/hm_antimicrobial_legend_noboxes.png"

# Save the heatmap as a PNG file
png(output_file, width = 2700, height = 1800, res = 200)
draw(ht_exp,heatmap_legend_side="bottom",legend_grouping="original",annotation_legend_side="left",  padding = unit(c(10, 30, 10, 10), "mm"))
dev.off()

# PCA and UMAP 
counts <- assay(vst(dds[,dds$treatment %in% c("Core","Mantle","Alveoli")]))
geneVars <- rowVars(as.matrix(counts))

# select just most variable genes for ACM 
select <- rownames(counts)[order(geneVars, decreasing=T)[1:500]]
pca <- prcomp(t(counts[select,]))
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
group<-meta[meta$treatment %in% c("Core","Mantle","Alveoli"),]$treatment
intgroup.df <- data.frame("label"=group,"shape"=sub("^\\D+", "",meta[meta$treatment %in% c("Core","Mantle","Alveoli"),]$batch))

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))
g<-ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke=1.33) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  coord_fixed() +
  # square plot 
  xlim(-pca_range, pca_range) +  
  ylim(-pca_range, pca_range) +  
  theme_classic() +
  labs(color = "Region", shape = "Patient")+
  scale_shape_manual(values = c(16, 15, 17, 18))+
  scale_color_manual(values = c("#F8766D","#00BA38","#619CFF")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=20,family = 'sans'),
        legend.position.inside = c(0.9, 0.85),legend.box="horizontal",
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.spacing.y = unit(2, "mm")) + 
  guides(color = guide_legend(override.aes = list(size = 5)),  
         shape = guide_legend(override.aes = list(size = 3)))

g_main <- g + theme(legend.position = "none")
ggsave(file = "plots/manuscript/pca.png", plot = g_main, device = "png", width=9, height=9, dpi = 400)

# Build the ggplot object to extract the legend
legend <- ggplot_gtable(ggplot_build(g))
legend_grob <- legend$grobs[[which(sapply(legend$grobs, function(x) x$name) == "guide-box")]]

# Turn off the default device to avoid empty file saving
dev.off()

# Draw the legend on a new page using grid system
grid::grid.newpage()
grid::grid.draw(legend_grob)

# Save the extracted legend to a PNG file
ggsave("plots/manuscript/pca_legend.png", grid::grid.grabExpr(grid::grid.draw(legend_grob)), device = "png", width = 4, height = 4, dpi = 400)

# run using PCA as init 
set.seed(0)
umap_result <- umap(pca$x)

d <- data.frame(PC1=umap_result$layout[,1], PC2=umap_result$layout[,2], group=group, intgroup.df)
pca_range <- max(abs(range(d$PC1, d$PC2)))

g2<-ggplot(data = d, aes(x = PC1, y = PC2, color = group, shape = shape)) +
  geom_point(size = 4, alpha = 0.6, stroke=1.33) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  coord_fixed() +
  # square plot 
  xlim(-pca_range, pca_range) +  
  ylim(-pca_range, pca_range) + 
  theme_classic() +
  labs(color = "Region", shape = "Patient")+
  scale_shape_manual(values = c(16, 15, 17, 18))+
  scale_color_manual(values = c("#F8766D","#00BA38","#619CFF")) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=20,family = 'sans'),
        legend.position = c(0.9, 0.85),legend.box="horizontal",
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.spacing.y = unit(2, "mm")) + 
  guides(color = guide_legend(override.aes = list(size = 5)),  
         shape = guide_legend(override.aes = list(size = 3)))
g2_main <- g2 + theme(legend.position = "none")
ggsave(file = "plots/manuscript/UMAP.png", plot = g2_main, device = "png", width=9,height=9, dpi = 400)

