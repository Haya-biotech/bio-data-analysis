# Group Project
# Haya Nahhas: Differential Gene Expression

#### Install and Load Required Packages ####
# Check if BiocManager is installed; install if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages required for microarray analysis
BiocManager::install(c("limma", "AnnotationDbi", "hgu133a.db"))
# Install CRAN packages for data manipulation and visualization
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))
# Load Bioconductor packages
library(AnnotationDbi)   
library(hgu133a.db)  
library(limma)           
library(dplyr)           
library(tibble)          
library(ggplot2)         
library(pheatmap)    
getwd()
list.files("C:/AI_Omics_Internship_2025/Differential_Expression_Gene_2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# Display objects available in the annotation package
ls("package:hgu133a.db")
columns(hgu133a.db)
keytypes(hgu133a.db)
getwd()
ls()
# Load preprocessed expression and phenotype data
data <- read.csv("C:/AI_Omics_Internship_2025/GSE15852_normalized_expression.csv", row.names = 1)

# Extract probe IDs from processed microarray data
probe_ids <- rownames(data)

# Map probe IDs to gene symbols using the platform annotation database
gene_symbols <- mapIds(
  hgu133a.db,          
  keys = probe_ids,        
  keytype = "PROBEID",    
  column = "SYMBOL",      
  multiVals = "first"      
)
# Convert mapping to a data frame and rename columns
gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)
# Summarize number of probes per gene symbol
duplicate_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))
# Identify genes associated with multiple probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)
# Verify if probe IDs in mapping correspond to expression data
all(gene_map_df$PROBEID == row.names(data))

# Merge annotation (SYMBOL) with expression matrix
processed_data_df <- data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = gene_symbols[PROBEID]) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)
# Remove probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

# Select only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)
# limma::avereps() computes the average for probes representing the same gene
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)
# Example to demonstrate how avereps works
x <- matrix(rnorm(8*3), 8, 3)
colnames(x) <- c("S1", "S2", "S3")
rownames(x) <- c("b", "a", "a", "c", "c", "b", "b", "b")
head(x)
avereps(x)  

dim(averaged_data)
# Convert averaged expression data to matrix format
data_2 <- as.data.frame(averaged_data)
data_2 <- data.matrix(data)
str(data_2)        
is.numeric(data_2) 
#### Differential Gene Expression Analysis ####
# -------------------------------------------------------------
# Define sample groups based on phenotype data
# Adjust group labels according to dataset annotation
library(GEOquery)
gse <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)
groups <- factor(pheno$source_name_ch1,
                 levels = c("normal breast tissue", "breast tumor tissue"),
                 labels = c("normal", "tumor"))
table(groups)
levels(groups)
unique(pheno$source_name_ch1)
# Create design matrix for linear modeling
# -------------------------------------------------------------
# Using no intercept (~0 + groups) allows each group to have its own coefficient
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
# Fit linear model to expression data
ncol(data_2)
nrow(design)
fit_1 <- lmFit(data_2, design)
# Define contrast to compare cancer vs normal samples
contrast_matrix <- makeContrasts(tumor_vs_normal = tumor - normal,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# Extract list of differentially expressed genes (DEGs)
# -------------------------------------------------------------
deg_results <- topTable(fit_2,
                        coef = "tumor_vs_normal",  
                        number = Inf,              
                        adjust.method = "BH")       
# Classify DEGs into Upregulated, Downregulated, or Not Significant
# -------------------------------------------------------------
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated",
         "No")
))
# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
# Combine both sets of DEGs
dir.create("Result_4")
deg_updown <- rbind(upregulated, downregulated)

write.csv(deg_results, file = "Result_4/DEGs_Results.csv")
write.csv(upregulated, file = "Result_4/Upregulated_DEGs.csv")
write.csv(downregulated, file = "Result_4/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "Result_4/Updown_DEGs.csv")
#### Data Visualization ####
# -------------------------------------------------------------

# -------------------------------------------------------------
# Volcano Plot: visualizes DEGs by logFC and adjusted p-values
# -------------------------------------------------------------
dir.create("Result_Plots")
png("Result_Plots/volcano_plot.png", width = 2000, height = 1500, res = 300)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(P-value)",
       color = "Regulation")
dev.off()

# Heatmap of Top Differentially Expressed Genes
# -------------------------------------------------------------

# Select top genes with smallest adjusted p-values
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 10)

# Subset averaged expression matrix for selected genes
heatmap_data <- data_2[top_genes, ]
# Generate unique column names per sample group for display
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))

# Assign formatted names to heatmap columns
colnames(heatmap_data) <- heatmap_names
png("heatmap_top10_DEGs.png", width = 2000, height = 1500, res = 300)
# Generate heatmap without additional scaling
pheatmap(
  heatmap_data,
  scale = "none", # for already normalized data
  cluster_rows = FALSE,              
  cluster_cols = TRUE,              
  show_rownames = TRUE,              
  show_colnames = TRUE,              
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 10 Differentially Expressed Genes"
)

dev.off()

save.image(file = "Differential_Expression_Gene.RData")