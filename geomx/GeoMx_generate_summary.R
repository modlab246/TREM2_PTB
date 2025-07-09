setwd("/Users/jonathanperrie/Documents/UCLA/projects/geomx")

# Load libraries
library(dplyr)
library(knitr)

# geomx_dcc_qc inputs 
pdata_file <- "data/pdata_table.csv"
protdata_file <- "data/protdata_table.csv"
sdata_file <- "data/sdata_table.csv"
qc_file <- "data/qc_pass_table.csv" 
gdr_file <- "data/spot_detection_table.csv"

# output file locations 
qc_output <- "data/qc_count_table.csv"
geomx_output <- "data/geomx_summary_table.csv"
geomx_dge_output <- "data/dge_summary_table.csv"

# basic metadata, nuclei, area 
table_df <- read.table(pdata_file, sep=",", header=TRUE, row.names=1, check.names=FALSE)
roi_df <- read.table(protdata_file, sep=",", header=TRUE, row.names=1)
table_df$roi <- sprintf("%02d", roi_df$x)
table_df$`Scan name` <- gsub("PTB", "PTB ", table_df$`Scan name`)

# read alignment 
align_df <- read.table(sdata_file, sep=",", header=TRUE, row.names=1)

table_df <- table_df[rownames(align_df), ] %>% cbind(align_df)

# Generate batch names
table_df_names <- apply(table_df[c("Segment tags", "Scan name", "roi")], 1, function(x) paste(x, collapse = " | "))
table_df$batch_names <- table_df_names

QCResults <- read.table(qc_file, sep=",", header=TRUE, row.names=1)

# Add QC status and Gene Detection Rate 
table_df$QC <- QCResults$QCStatus
table_df$gdr <- 1

gdr_df <- read.table(gdr_file, header = TRUE, sep = ",", row.names=1)
GeneDetectionRate <- setNames(gdr_df$gdr, gdr_df$name)
table_df$gdr[table_df$batch_names %in% names(GeneDetectionRate)] <- GeneDetectionRate
table_df$gdr <- table_df$gdr * 100


# exclude 22.2 - 19 (core-mantle) and other annotations 
other_bool <- grepl("Other", table_df$batch_names)
roi_19_bool <- grepl("\\| PTB 22.2 \\| 019", table_df$batch_names)
ptb_bool <- grepl("^PTB", table_df$`Scan name`)
table_df$post_qc <- ptb_bool & (other_bool | roi_19_bool)

# table 1: counts and QC/GDR summary (79 spots)
res1 <- table_df %>%
  filter(Tissue == "Lung") %>%
  group_by(`Scan name`, Tissue, `Segment tags`) %>%
  summarise(
    Region_count = n(),
    QC_count = sum(QC != "PASS", na.rm = TRUE),
    GDR_count = sum(gdr < 5, na.rm = TRUE),
    post_QC = sum(post_qc, na.rm = TRUE)
  ) %>%
  mutate(
    Max_QC = pmax(QC_count, GDR_count, post_QC, na.rm = TRUE) # Max of the three columns
  )

summary_row <- res1 %>%
  ungroup() %>%  
  summarise(
    `Scan name` = "Total",
    Tissue = "Lung",
    `Segment tags` = "All",
    Region_count = sum(Region_count, na.rm = TRUE),
    QC_count = sum(QC_count, na.rm = TRUE),
    GDR_count = sum(GDR_count, na.rm = TRUE),
    post_QC = sum(post_QC, na.rm = TRUE),
    Max_QC = sum(Max_QC, na.rm = TRUE)
  )

res1 <- bind_rows(res1, summary_row)

# table 2: averages and summaries for GDR > 0.05 and QC == "PASS" (70 spots)
res2 <- table_df %>%
  filter(Tissue == "Lung", gdr >= 5, gdr <= 100, QC == "PASS") %>%
  group_by(`Scan name`, Tissue, `Segment tags`) %>%
  summarise(
    Avg_NucleiCount = mean(AOINucleiCount, na.rm = TRUE),
    Avg_SurfaceArea = mean(AOISurfaceArea, na.rm = TRUE),
    Sum_AlignedReads = sum(Aligned_Count, na.rm = TRUE),
    Avg_AlignedP = mean(Aligned_Percentage, na.rm = TRUE),
    Region_Count = n(),
    GDR_mean = mean(gdr, na.rm = TRUE),
    .groups = "drop"
  )

summary_row <- table_df %>%
  filter(
    Tissue == "Lung",
    gdr    >=  5,
    gdr    <=100,
    QC     == "PASS"
  ) %>%
  summarise(
    `Scan name`      = "Total",
    Tissue           = "Lung",
    `Segment tags`   = "All",
    Avg_NucleiCount  = mean(AOINucleiCount,   na.rm = TRUE),
    Avg_SurfaceArea  = mean(AOISurfaceArea,   na.rm = TRUE),
    Sum_AlignedReads = sum(Aligned_Count,      na.rm = TRUE),
    Avg_AlignedP     = mean(Aligned_Percentage, na.rm = TRUE),
    Region_Count     = n(),
    GDR_mean         = mean(gdr,               na.rm = TRUE)
  )
res2 <- bind_rows(res2, summary_row)

res2 <- res2 %>%
  mutate(
    Avg_NucleiCount = as.integer(round(Avg_NucleiCount)),          
    Avg_SurfaceArea = gsub("e", "E", formatC(Avg_SurfaceArea, format = "e", digits = 2)), 
    Sum_AlignedReads = gsub("e", "E", formatC(Sum_AlignedReads, format = "e", digits = 2)), 
    Avg_AlignedP = as.character(sprintf("%.1f", round(Avg_AlignedP, 1))),       
    GDR_mean = as.character(sprintf("%.1f", round(GDR_mean, 1)))
  )

# table 3: averages and summaries excluding Post-QC Regions
res3 <- table_df %>%
  filter(Tissue == "Lung", gdr >= 5, gdr <= 100, QC == "PASS", !post_qc) %>%
  group_by(`Scan name`, Tissue, `Segment tags`) %>%
  summarise(
    Avg_NucleiCount = mean(AOINucleiCount, na.rm = TRUE),
    Avg_SurfaceArea = mean(AOISurfaceArea, na.rm = TRUE),
    Sum_AlignedReads = sum(Aligned_Count, na.rm = TRUE),
    Avg_AlignedP = mean(Aligned_Percentage, na.rm = TRUE),
    Region_Count = n(),
    GDR_mean = mean(gdr, na.rm = TRUE),
    .groups = "drop"
  ) 

summary_row <- table_df %>%
  filter(Tissue == "Lung", gdr >= 5, gdr <= 100, QC == "PASS", !post_qc) %>%
  summarise(
    `Scan name`      = "Total",
    Tissue           = "Lung",
    `Segment tags`   = "All",
    Avg_NucleiCount  = mean(AOINucleiCount,   na.rm = TRUE),
    Avg_SurfaceArea  = mean(AOISurfaceArea,   na.rm = TRUE),
    Sum_AlignedReads = sum(Aligned_Count,      na.rm = TRUE),
    Avg_AlignedP     = mean(Aligned_Percentage, na.rm = TRUE),
    Region_Count     = n(),
    GDR_mean         = mean(gdr,               na.rm = TRUE)
  )
res3 <- bind_rows(res3, summary_row)

res3 <- res3 %>%
  mutate(
    Avg_NucleiCount = as.integer(round(Avg_NucleiCount)),         
    Avg_SurfaceArea = gsub("e", "E", formatC(Avg_SurfaceArea, format = "e", digits = 2)), 
    Sum_AlignedReads = gsub("e", "E", formatC(Sum_AlignedReads, format = "e", digits = 2)), 
    Avg_AlignedP = as.character(sprintf("%.1f", round(Avg_AlignedP, 1))),     
    GDR_mean = as.character(sprintf("%.1f", round(GDR_mean, 1)))
  )

# Display the results
print("Res1: Counts and QC/GDR Summary")
colnames(res1) <- c("Patient", "Tissue", "Annotation", "Count", "Read QC", "Detection QC", "Post QC", "Total QC")
print(res1)

print("Res2: Averages and Summaries for gdr > 0.05 and QC == 'PASS'")
colnames(res2) <- c("Patient", "Tissue", "Annotation", "Nuclei per ROI", "Area per ROI", "Total aligned", "Average aligned %", "Count", "Gene detection %")
print(res2)

print("Res3: Excluding Post-QC Regions")
colnames(res3) <- c("Patient", "Tissue", "Annotation", "Nuclei per ROI", "Area per ROI", "Total aligned", "Average aligned %", "Count", "Gene detection %")
print(res3)

# .csv table format 
write.table(res1[,c("Patient","Annotation","Count","Read QC","Detection QC","Post QC", "Total QC")], qc_output, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(res2[,c("Patient","Annotation","Count","Average aligned %","Total aligned","Nuclei per ROI","Gene detection %")], geomx_output, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(res3[,c("Patient","Annotation","Count","Average aligned %","Total aligned","Nuclei per ROI","Gene detection %")], geomx_dge_output, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)