setwd("C:\\AI_Omics_Internship_2025")
getwd()

if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
install.packages("dplyr")
BiocManager::install("limma")
BiocManager::install("curl")


library(GEOquery)           
library(affy)                 
library(arrayQualityMetrics)  
library(dplyr)
library(curl)

gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)

View(gse_data)

expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)

feature_data <-  fData(gse_data$GSE79973_series_matrix.txt.gz)
phenotype_data <-  pData(gse_data$GSE79973_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1)) 

getGEOSuppFiles("GSE79973", baseDir = "Raw_data_h", makeDirectory = TRUE)
dir.create("raw_data_h")

untar("raw_data_h\\GSE79973_RAW.tar", exdir = "raw_data_h\\CEL_Files")

raw_data <- ReadAffy(celfile.path = "Normaliztion-Class_3B.R/raw_data_h/CEL_Files")

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)

normalized_data <- rma(raw_data)
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)

processed_data <- as.data.frame(exprs(normalized_data))

dim(processed_data) 
row_median <- rowMedians(as.matrix(processed_data))
row_median
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2) 
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ] 
colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data 
class(phenotype_data$source_name_ch1) 
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))
class(groups)
levels(groups)
