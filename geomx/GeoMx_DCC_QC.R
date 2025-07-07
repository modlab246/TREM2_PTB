
# This script processes NanoString GeoMx data for spatial genomic 
# analysis. It includes steps for loading, preprocessing, quality 
# control, and normalization of DCC files across multiple batches. 
# Outputs include raw and normalized count data, nuclei metadata
# and negative probe counts. 


setwd("/Users/jonathanperrie/Documents/UCLA/projects/geomx")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(readxl)
library(openxlsx)
library(knitr)
library(ggplot2)

# Input and output parameters 

# digital count conversion filepaths 
dcc_b1 <- "dcc/DCC_batch_1"
dcc_b2 <- "dcc/DCC_batch_2"
dcc_b3 <- "dcc/DCC_batch_3"
dcc_b4 <- "dcc/DCC_batch_4"
dcc_b5 <- "dcc/DCC_batch_5"

# b1 = PTB 21.1
# b2 = TBL 7, 16, 42
# b3 = PTB 22.2, TBL 60
# b4 = PTB 22.1, 22.3
# b5 = EPTB 22.2, TBL 62
DCCFiles <- c(
  dir(dcc_b1, pattern = ".dcc$", full.names = TRUE, recursive = TRUE),
  dir(dcc_b2, pattern = ".dcc$", full.names = TRUE, recursive = TRUE),
  dir(dcc_b3, pattern = ".dcc$", full.names = TRUE, recursive = TRUE),
  dir(dcc_b4, pattern = ".dcc$", full.names = TRUE, recursive = TRUE),
  dir(dcc_b5, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
)

# Panel kit configuration 
PKCFiles <- "dcc/Hs_R_NGS_WTA_v1.0.pkc"

# GeoMx annotation metadata 
SampleAnnotationFile <- "dcc/geomx_ptb_annotations.xlsx"

# Output file locations
nuclei_output <- "data/geomx_ptb_metadata.csv"
negprobes_output <- "data/geomx_ptb_negprobes.csv"
raw_data_output <- "data/geomx_ptb_raw_data.csv"
q3_data_output <- "data/geomx_ptb_q3_data.csv"
qc_output <- "data/qc_pass_table.csv"
pdata_output <- "data/pdata_table.csv"
sdata_output <- "data/sdata_table.csv"
protdata_output <- "data/protdata_table.csv"
gdr_output <- "data/spot_detection_table.csv"

# data processing script 

# Load data 
batched_data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                       pkcFiles = PKCFiles,
                                       phenoDataFile = SampleAnnotationFile,
                                       phenoDataSheet = "Sheet1",
                                       phenoDataDccColName = "Sample_ID",
                                       protocolDataColNames = c("ROI","Segment label"))

# Offset and flag 
batched_data <- shiftCountsOne(batched_data, useDALogic = TRUE)

QC_params <- list(minSegmentReads = 1000, 
       percentTrimmed = 80,    
       percentStitched = 80,  
       percentAligned = 75,    
       percentSaturation = 50, # Minimum sequencing saturation 
       minNegativeCount = 1,   # Minimum negative control counts 
       maxNTCCount = 9000,     # Maximum counts observed in NTC well 
       minNuclei = 20,        
       minArea = 1000)        

batched_data <- setSegmentQCFlags(batched_data, 
                                  qcCutoffs = QC_params)  

write.table(pData(batched_data)[c("AOINucleiCount", "AOISurfaceArea", "Scan name", "Segment tags", "Tissue")], pdata_output, sep=",", col.names=NA)
write.table(protocolData(batched_data)[["ROI"]], protdata_output, sep=",", col.names=NA)

sdata_clean <- sData(batched_data)[c("Aligned (%)", "Aligned")]
colnames(sdata_clean) <- c("Aligned_Percentage","Aligned_Count")
sdata_clean$Aligned_Percentage <- as.numeric(unlist(sdata_clean$Aligned_Percentage))
sdata_clean$Aligned_Count <- as.numeric(unlist(sdata_clean$Aligned_Count))
write.table(sdata_clean, sdata_output, sep=",", col.names=NA)


# Summarize QC results
QCResults <- protocolData(batched_data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1, function(x) ifelse(sum(x) == 0, "PASS", "WARNING"))

write.table(QCResults, qc_output, sep=",", col.names=NA)

QC_Summary["TOTAL FLAGS", ] <- c(sum(QCResults[, "QCStatus"] == "PASS"),
                                 sum(QCResults[, "QCStatus"] == "WARNING"))

# Calculate negative geometric means
negativeGeoMeans <- esBy(negativeControlSubset(batched_data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(batched_data)[["NegGeoMean"]] <- negativeGeoMeans

pkcs <- annotation(batched_data)
modules <- gsub(".pkc", "", pkcs)
negCols <- paste0("NegGeoMean_", modules)
pData(batched_data)[, negCols] <- sData(batched_data)[["NegGeoMean"]]

# Detatch negative geometric mean column ahead of aggregateCounts call
pData(batched_data) <- pData(batched_data)[, !colnames(pData(batched_data)) %in% negCols]
batched_data <- batched_data[, QCResults$QCStatus == "PASS"]


# Set BioProbe QC flags based on probe ratio and Grubb’s test outlier detection
batched_data <- setBioProbeQCFlags(batched_data, 
                                   qcCutoffs = list(minProbeRatio = 0.1,
                                                    percentFailGrubbs = 20), 
                                   removeLocalOutliers = TRUE)

# Define QC metrics for each probe based on outcomes of QC tests
ProbeQCResults <- fData(batched_data)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# Subset data to include only probes that passed QC checks
ProbeQCPassed <- subset(batched_data,
                        fData(batched_data)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
                          fData(batched_data)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE)


batched_data <- ProbeQCPassed 

# Aggregate counts for the final data set
target_batched_data <- aggregateCounts(batched_data)
pData(target_batched_data)$roi <- protocolData(target_batched_data)[["ROI"]]

# filter out data sets we have ignored and DCC control files 
index <- !grepl("TBL 7", pData(target_batched_data)$`Scan name`) & !grepl("TBL 60", pData(target_batched_data)$`Scan name`) & 
  !grepl("No Template Control", pData(target_batched_data)$`Segment tags`)


# Separate out metadata and counts 
batch_meta<-cbind(pData(target_batched_data)[index,],sData(target_batched_data)[index,])
batch_count<-exprs(target_batched_data)[,index]

# Separate out negative probes 
negativeProbefData <- subset(fData(batched_data), CodeClass == "Negative")
batch_ngpb <- exprs(batched_data[fData(batched_data)$TargetName %in% negativeProbefData$TargetName,])
batch_ngpb <- batch_ngpb[,index]

# Comprehensive row names: region, sample, ROI
batch_names <- pData(target_batched_data)[colnames(batch_count)[grepl("dcc",colnames(batch_count))],c("Segment tags","Scan name","roi")]
batch_names$roi <- sprintf("%02d",batch_names$roi)
batch_names$`Scan name` <- gsub("PTB","PTB ",batch_names$`Scan name`)
batch_names <- apply(batch_names,1,function(x) paste(x, collapse =" | "))

# update names 
colnames(batch_count) <- batch_names[colnames(batch_count)]
colnames(batch_ngpb) <- batch_names[colnames(batch_ngpb)]

# Data structure for nuclei counts 
nuclei_names <- c("Sub-classification","AOINucleiCount","AOISurfaceArea","ROICoordinateX","ROICoordinateY",
                  "umiQ30","rtsQ30","Segment label","SoftwareVersion",
                  "Raw","Trimmed","Stitched","Aligned","DeduplicatedReads")
nuclei <- rbind(batch_meta[,nuclei_names])
nuclei["TrimmedPercent"] <- unname(unlist(batch_meta["Trimmed (%)"]))
nuclei["StitchedPercent"] <- unname(unlist(batch_meta["Stitched (%)"]))
nuclei["AlignedPercent"] <- unname(unlist(batch_meta["Aligned (%)"]))
nuclei["SaturatedPercent"] <- unname(unlist(batch_meta["Saturated (%)"]))
nuclei["DCC"] <- rownames(nuclei)
rownames(nuclei) <- c(batch_names[rownames(batch_meta["AOINucleiCount"])])
colnames(nuclei)[colnames(nuclei) == "SoftwareVersion"] <- "GeoMxPipelineNGS"
colnames(nuclei)[colnames(nuclei) == "DeduplicatedReads"] <- "Deduplicated"
colnames(nuclei)[colnames(nuclei) == "SegmentLabel"] <- "Segment"

# Calculate LOQ 
pkcs <- annotation(batched_data)
modules <- gsub(".pkc", "", pkcs)
module <- modules[1]

# LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

LOQ <- data.frame(row.names = colnames(batch_count))

NegGeoMean <- apply(batch_ngpb,2,ngeoMean)
NegGeoSD <- apply(batch_ngpb,2,ngeoSD)

LOQ <- data.frame(row.names = colnames(batch_count))
LOQ[, module] <-pmax(minLOQ, NegGeoMean * NegGeoSD ^ cutoff)

# Threshold counts against LOQ 
LOQ_Mat <- c()
Mat_i <- apply(batch_count,1,FUN = function(x) { x > LOQ})
LOQ_Mat <- t(rbind(LOQ_Mat, Mat_i))
colnames(LOQ_Mat) <- colnames(batch_count)

# Filter out spots by number of genes detected 
GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
GeneDetectionRate <- GenesDetected / nrow(batch_count)

# Convert GeneDetectionRate to a data frame
gdr_df <- data.frame(
  name = names(GeneDetectionRate),
  gdr = as.numeric(GeneDetectionRate)
)
write.table(gdr_df, gdr_output, sep=",", col.names=NA)

DetectionThreshold <- cut(GeneDetectionRate,
                          breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
                          labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

batch_count <- batch_count[,GeneDetectionRate >= .05]
batch_ngpb <- batch_ngpb[,GeneDetectionRate >= .05]
nuclei <- nuclei[GeneDetectionRate >= .05,]

# Filter out genes by spot detection 
LOQ_Mat <- LOQ_Mat[, colnames(batch_count)]
DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
DetectionRate <- DetectedSegments / ncol(batch_count)

batch_count <- batch_count[DetectionRate >= 0.05,]

# Q3 normalization
qs <- apply(batch_count, 2, function(x) stats::quantile(x, 0.75))
q3 <- sweep(batch_count, 2L, qs / ngeoMean(qs), FUN = "/")

# Save to CSV
write.table(nuclei[grepl(" PTB",rownames(nuclei)),], nuclei_output, sep=",", col.names=NA)
write.table(batch_ngpb[,grepl(" PTB",colnames(batch_ngpb))], negprobes_output, sep=",", col.names=NA)
write.table(batch_count[,grepl(" PTB",colnames(batch_count))], raw_data_output, sep=",", col.names=NA)
write.table(q3[,grepl(" PTB",colnames(q3))], q3_data_output, sep=",", col.names=NA)