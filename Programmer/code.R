setwd("/projectnb/bf528/users/group1/xdhan/project4/programmer/")
install.packages("Seurat")
install.packages("dplyr")
library(Seurat)
library(tximport)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(Matrix.utils)
#Get the UMI matrix paht
files <- file.path("/projectnb/bf528/users/group1/project4_data/umi_matrix.csv")
file.exists(files)

#Load matrix
#txi <- tximport(files, type="alevin")
umi_matrix <- read.csv(file=files,header=T)
ds <-umi_matrix
rownames(ds) <- umi_matrix[,1]
ds <- ds[,-1]
#Convert Ensemble Id to gene symbol
ENSGID <- rownames(ds)
# Remove decimal
for (i in c(1:length(ENSGID))){
  ENSGID[i] = strsplit(as.character(ENSGID[i]),"\\.")[[1]][1]
}
symbolID <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ENSGID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# Remove the rows without symbol name
ds <- ds[!is.na(match(ENSGID,symbolID$GENEID)),]
ds$symbol <- symbolID$SYMBOL
ds <- as.data.frame(ds)
# Merge the rows with same name (gene symbol)
ds <- ds %>% group_by(symbol) %>% summarise_all(funs(sum))
dim(ds)
ds <- as.data.frame(ds)
rownames(ds)<-ds$symbol
ds <- ds[,colnames(ds)!="symbol"]
# Count the reads of matrix
pbmc <- CreateSeuratObject(counts = ds , min.cells = 3, min.features = 200, project = "10X_PBMC")
dim(pbmc@meta.data)
#mito.genes <- grep("^MT-", rownames(pbmc@meta), value = T)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Subset the qualified cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt<5)
dim(pbmc@meta.data)
head(pbmc@meta.data)

#Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
head(pbmc@meta.data)

# Filter out low variance genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = T)
plot2
# Scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = TRUE)

# Determine the dimensionality
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

# Cluster Cells
pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
# Non-Linear Dimensional Reduction
pbmc <- RunUMAP(pbmc, dims = 1:7)
DimPlot(pbmc, reduction = "umap")
pbmc$seurat_clusters
sum(pbmc$seurat_clusters==0)
sum(pbmc$seurat_clusters==1)
sum(pbmc$seurat_clusters==2)
sum(pbmc$seurat_clusters==3)
sum(pbmc$seurat_clusters==4)
saveRDS(pbmc, file = "pbmc3k_final.rds")
