# img 
This directory contains all image related work on the raw tiff files getting the bounding boxes and labels properly assigned 

# Manuscript plots
1. UMAP/PCA
2. Overlaps
3. Table #2
4. Volcano plots
5. Pathway viz
6. GeoMx heatmap

# GeoMx_DCC_QC.R
* Reads in DCC files and applies threshold to ROIs and genes 
* Saves just PTB files that pass initial QC

# GeoMx_DGE_plots.R
* Pairwise contrasts with DESeq2
* UMAP/PCA (1)
* Volcano plots (4)
* GeoMx heatmap (6)
* Generates a spreadsheet of all significant genes across every region

# GeoMx_generate_summary.R
* Tables from GeoMx data
* Table #2 (3)

# pathway_plots_dend.ipynb
* Reads in pathways highlighted by Robert
* Clusters them based on their gene ratio scores
  * hits / total gene set
* Generate a heatmap ordered by clustering (5)

# marker_overlap.ipynb
* Find Jaccard overlaps between reference sets and GeoMx genes within a region (ACM)
* Saves spreadsheets of all overlapping genes
* Generates heatmaps of Jaccard scores across ACM (2)
