if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Seurat", "data.table"))

library(GEOquery)
library(Seurat)

# 1. Download the GEO accession
gse <- getGEO("GSE115978", GSEMatrix = TRUE)
a
library(Seurat)
library(data.table)
library(dplyr)

# --- 1. Load the Count Matrix ---
# Read the large text file (fread is much faster than read.csv)
print("Loading count matrix...")
counts <- fread("data/GSE115978_counts.csv.gz")

# Convert to a standard matrix format (Genes as rows, Cells as columns)
# Note: Check if genes are rows or columns in the file. 
# Usually, rows=genes. If not, use t() to transpose.
counts_matrix <- as.matrix(counts[,-1]) 
rownames(counts_matrix) <- counts$V1 # Assuming first column is gene names

# --- 2. Create Seurat Object ---
# This is the container for your entire analysis
seu_obj <- CreateSeuratObject(
  counts = counts_matrix,
  project = "Melanoma_ICB",
  min.cells = 3,       # Filter genes expressed in <3 cells
  min.features = 200   # Filter cells with <200 genes (empty droplets)
)

# --- 3. Add Mitochondrial Percentage (QC) ---
# This is critical for the "Quality Control" bullet point
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")

# --- 4. Visualize QC Metrics ---
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

print("Data loaded successfully!")
print(seu_obj)
counts <- fread(file.choose())
library(Seurat)
library(data.table)

# 1. Clear old data
rm(list = ls()) 
gc()

# 2. POINT R TO YOUR NEW E: DRIVE FOLDER
# (Change "E:" to "G:" if you used the G drive)
setwd("E:/T_Cell_Project")

# 3. Load the data
print("Loading data from E: drive...")
counts <- fread("GSE115978_counts.csv.gz")

# 4. Verify dimensions
print(dim(counts))

# 5. Initialize Seurat Object
counts_matrix <- as.matrix(counts[,-1])
rownames(counts_matrix) <- counts$V1

seu_obj <- CreateSeuratObject(counts = counts_matrix, project = "Melanoma_ICB")
print(seu_obj)
seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 1. Check if there are any genes starting with "MT-" (Human standard)
print(grep("^MT-", rownames(seu_obj), value = TRUE))

# 2. Check if there are any genes starting with "mt-" (Mouse/Alternative standard)
print(grep("^mt-", rownames(seu_obj), value = TRUE))
# --- STEP 1: Normalization ---
# Log-normalize the data to make it comparable across cells
print("Normalizing data...")
seu_obj <- NormalizeData(seu_obj)

# --- STEP 2: Find Variable Features ---
# Identify the top 2000 genes that vary the most (these drive the biology)
print("Finding variable features...")
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)

# --- STEP 3: Scaling ---
# Scale data so all genes have the same weight (mean=0, variance=1)
# Note: This might take 1-2 minutes
print("Scaling data...")
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)

# --- STEP 4: PCA (Principal Component Analysis) ---
# Reduce dimensions to the top 50 "directions" of variance
print("Running PCA...")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))

# --- STEP 5: Clustering & UMAP ---
# This creates the "groups" of cells based on their gene expression
print("Running UMAP and Clustering...")
seu_obj <- FindNeighbors(seu_obj, dims = 1:20) # Use top 20 PCs
seu_obj <- FindClusters(seu_obj, resolution = 0.5) # Resolution 0.5 is standard
seu_obj <- RunUMAP(seu_obj, dims = 1:20)

# --- STEP 6: Visualizing the Result ---
# This is the "money shot" - your first map of the tumor!
# We look for specific "flag" genes to name the clusters.

# CD3D = All T-Cells (The ones we want!)
# CD8A = Killer T-Cells (The ones that get exhausted)
# CD4  = Helper T-Cells
# MS4A1 = B-Cells (Neighbors, not our focus)
# PMEL = Melanoma Tumor Cells (The enemy)

FeaturePlot(seu_obj, features = c("CD3D", "CD8A", "CD4", "MS4A1", "PMEL"))
DimPlot(seu_obj, reduction = "umap", label = TRUE)
t_cells <- subset(seu_obj, idents = c(0, 1, 11))
print("Re-processing T-cells...")
t_cells <- NormalizeData(t_cells)
t_cells <- FindVariableFeatures(t_cells)
t_cells <- ScaleData(t_cells)
t_cells <- RunPCA(t_cells)
t_cells <- RunUMAP(t_cells, dims = 1:20)
t_cells <- FindNeighbors(t_cells, dims = 1:20)
t_cells <- FindClusters(t_cells, resolution = 0.5)

# --- Show the T-Cell Map ---
DimPlot(t_cells, reduction = "umap", label = TRUE)
# --- Identify Exhaustion ---
# We look for 4 key markers:
# 1. CD8A   = Killer T-cells (The base type)
# 2. CD4    = Helper T-cells (Different lineage)
# 3. PDCD1  = PD-1 (The main exhaustion marker)
# 4. HAVCR2 = TIM-3 (Another strong exhaustion marker)

FeaturePlot(t_cells, features = c("CD8A", "CD4", "PDCD1", "HAVCR2", "LAG3"), ncol = 3)
new.cluster.ids <- c("Exhausted CD8+", "CD4+ Helper", "CD4+ Helper", 
                     "CD4+ Treg", "Functional CD8+", "Proliferating", "Other")
names(new.cluster.ids) <- levels(t_cells)
t_cells <- RenameIdents(t_cells, new.cluster.ids)
DimPlot(t_cells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
print("Finding exhaustion markers...")
exhaustion_markers <- FindMarkers(t_cells, ident.1 = "Exhausted CD8+", ident.2 = "Functional CD8+")

# Show the top 10 genes driving exhaustion
print(head(exhaustion_markers, n = 10))
print("Finding UPREGULATED exhaustion markers...")
up_markers <- FindMarkers(t_cells, 
                          ident.1 = "Exhausted CD8+", 
                          ident.2 = "Functional CD8+", 
                          only.pos = TRUE,
                          min.pct = 0.25)

# Show the top 10 genes that define exhaustion
print(head(up_markers, n = 10))
ggsave("Final_Exhaustion_Map.png", width = 10, height = 7)
library(ggplot2)
ggsave("Final_Exhaustion_Map.png", width = 10, height = 7)
ggsave("Final_Exhaustion_Map.png", width = 10, height = 7)
ggplot2::ggsave("Final_Exhaustion_Map.png", width = 10, height = 7)
print(head(t_cells@meta.data))
# --- Clinical Correlation (Patient vs. Exhaustion) ---

# 1. Calculate the percentage of each cell type per patient
library(ggplot2) # Loading it just in case

# Create a frequency table (Patient x Cell Type)
cell_counts <- table(t_cells$orig.ident, Idents(t_cells))

# Convert to percentages (so patients with more cells don't dominate)
cell_props <- prop.table(cell_counts, margin = 1) * 100
cell_props_df <- as.data.frame(cell_props)
colnames(cell_props_df) <- c("Patient_ID", "Cell_Type", "Percentage")

# 2. Plot the Stacked Bar Chart
ggplot(cell_props_df, aes(x = Patient_ID, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "T-Cell Composition per Patient",
       y = "Percentage of T-Cells",
       x = "Patient ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate labels so they fit

# 3. Save this final clinical chart
ggplot2::ggsave("Clinical_Correlation_Chart.png", width = 12, height = 6)
# Create a frequency table (Patient x Cell Type)
cell_counts <- table(t_cells$orig.ident, Idents(t_cells))

# Convert to percentages (so patients with more cells don't dominate)
cell_props <- prop.table(cell_counts, margin = 1) * 100
cell_props_df <- as.data.frame(cell_props)
colnames(cell_props_df) <- c("Patient_ID", "Cell_Type", "Percentage")

# 2. Plot the Stacked Bar Chart
ggplot(cell_props_df, aes(x = Patient_ID, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "T-Cell Composition per Patient",
       y = "Percentage of T-Cells",
       x = "Patient ID") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate labels so they fit

# 3. Save this final clinical chart
ggplot2::ggsave("Clinical_Correlation_Chart.png", width = 12, height = 6)
# --- Clinical Correlation (The "Explicit" Fix) ---

# 1. Prepare the data (This part is standard)
cell_counts <- table(t_cells$orig.ident, Idents(t_cells))
cell_props <- prop.table(cell_counts, margin = 1) * 100
cell_props_df <- as.data.frame(cell_props)
colnames(cell_props_df) <- c("Patient_ID", "Cell_Type", "Percentage")

# 2. Plot using EXPLICIT commands
# We put 'ggplot2::' before EVERYTHING to bypass the error
plot <- ggplot2::ggplot(cell_props_df, ggplot2::aes(x = Patient_ID, y = Percentage, fill = Cell_Type)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Clinical Correlation: T-Cell Exhaustion by Patient",
                y = "Percentage of T-Cells",
                x = "Patient ID") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

# 3. Display the chart
print(plot)

# 4. Save the chart
ggplot2::ggsave("Clinical_Correlation_Chart.png", plot = plot, width = 12, height = 6)
FeaturePlot(t_cells, features = c("TCF7", "PDCD1"), blend = TRUE)
FeaturePlot(t_cells, features = c("CD69", "ITGAE", "PDCD1", "TCF7"), ncol = 2)
FeaturePlot(t_cells, features = c("GAPDH", "FABP5"), blend = TRUE)
# --- CAR-T PIVOT: Signaling Domains ---
# We check for the "Survival Signals" used in CAR-T engineering.

# TNFRSF9 = 4-1BB (Promotes Fatty Acid Metabolism & Longevity)
# ICOS = Inducible Co-stimulator (Promotes persistence)
# CD28 = The "Fast & Furious" signal (Promotes Glycolysis)

FeaturePlot(t_cells, features = c("TNFRSF9", "ICOS", "CD28", "GAPDH"), ncol = 2)
saveRDS(t_cells, file = "t_cells_processed.rds")
