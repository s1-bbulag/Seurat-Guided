library(dplyr)
library(Seurat)
library(patchwork)
library(httpgd)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
options(future.globals.maxSize = 1e9)


update.packages(oldPkgs = c("withr", "rlang"), type = "source")

if (!requireNamespace('remotes', quietly = TRUE)) {
  install.packages('remotes')
}
remotes::install_github('satijalab/azimuth', ref = 'master')

remotes::install_github('satijalab/seurat-wrappers')

gse756dt <- read.table("C:/Users/bbula/OneDrive/Desktop/VS Code DS/DSCode/Seurat-Guided/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", header = TRUE, sep = "\t")
dim(gse756dt)
# Read in entire file
gse756d <- read.delim("C:/Users/bbula/OneDrive/Desktop/VS Code DS/DSCode/Seurat-Guided/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", header = TRUE, stringsAsFactors = FALSE)
View(gse756d) 
head(rownames(gse756d))

gse756genes <- gse756d[,2]
View(gse756genes)

dim(gse756d)


#selecting to remove cols 1-18
cellgse756_columns = gse756d[, 18:ncol(gse756d)]
head(colnames(cellgse756_columns))
head(rownames(cellgse756_columns))
# Set the row names of cellgse756_columns to be gse756genes
# Make duplicate gene names unique
gse756genes <- make.unique(gse756genes)

# Now try setting the row names again
rownames(cellgse756_columns) <- gse756genes

# Create a Seurat gse756d_seuratect
gse756d_seurat <- CreateSeuratObject(counts = cellgse756_columns)
#gse756d_seurat <- CreateSeuratgse756d_seuratect(counts = gse756d)
#the above was a former attempt 

dim(gse756d_seurat)
summary(gse756d_seurat@meta.data)
summary(gse756d_seurat@meta.data$nCount_RNA)
summary(gse756d_seurat@meta.data$nFeature_RNA)
dim(gse756d_seurat@meta.data)

gse756d_seurat@meta.data$nCount=gse756d_seurat@meta.data$nFeature_RNA
gse756d_seurat[["percent.mt"]] <- PercentageFeatureSet(gse756d_seurat, pattern = "^MT-")
names(gse756d_seurat@meta.data)
View(gse756d_seurat)

gse756d_seurat


# The [[ operator can add columns to gse756d_seuratect metadata. This is a great place to stash QC stats
gse756d_seurat[["percent.mt"]] <- PercentageFeatureSet(gse756d_seurat, pattern = "^MT-")


head(gse756d_seurat@meta.data, 20)


gse756d_seurat@meta.data$nCount=NULL

gse756d_seurat@meta.data$nCount_RNA=NULL

gse756d_seurat@meta.data$percent.mt=NULL
VlnPlot(gse756d_seurat, features = "nFeature_RNA")


plot1 <- FeatureScatter(gse756d_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gse756d_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FindVariableFeatures(gse756d_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



gse756d_seurat <- subset(gse756d_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

gse756d_seurat <- NormalizeData(gse756d_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


gse756d_seurat <- NormalizeData(gse756d_seurat)


gse756d_seurat <- FindVariableFeatures(gse756d_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gse756d_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gse756d_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(gse756d_seurat)
gse756d_seurat <- ScaleData(gse756d_seurat, features = all.genes)

dim(gse756d_seurat)



gse756d_seurat <- RunPCA(gse756d_seurat, features = VariableFeatures(gse756d_seuratect = gse756d_seurat))

# Examine and visualize PCA results a few different ways
print(gse756d_seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(gse756d_seurat, dims = 1:2, reduction = "pca")

DimHeatmap(gse756d_seurat, dims = 1, cells = 500, balanced = TRUE)

ElbowPlot(gse756d_seurat)

gse756d_seurat <- FindNeighbors(gse756d_seurat, dims = 1:10)
gse756d_seurat <- FindClusters(gse756d_seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(gse756d_seurat), 5)

gse756d_seurat <- RunUMAP(gse756d_seurat, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(gse756d_seurat, reduction = "umap")

saveRDS(gse756d_seurat, file = "../output/gse756d_seurat_tutorial.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(gse756d_seurat, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(gse756d_seurat, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)



gse756d_seurat <- NormalizeData(gse756d_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
gse756d_seurat <- NormalizeData(gse756d_seurat)


gse756d_seurat <- FindVariableFeatures(gse756d_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gse756d_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(gse756d_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(gse756d_seurat)
gse756d_seurat <- ScaleData(gse756d_seurat, features = all.genes)
gse756d_seurat <- RunPCA(gse756d_seurat, features = VariableFeatures(gse756d_seuratect = gse756d_seurat))
VizDimLoadings(gse756d_seurat, dims = 1:2, reduction = "pca")

gse756d_seurat <- FindNeighbors(gse756d_seurat, dims = 1:10)
gse756d_seurat <- FindClusters(gse756d_seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(gse756d_seurat), 5)



gse756d_seurat <- RunUMAP(gse756d_seurat, dims = 1:10)




# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(gse756d_seurat, reduction = "umap")


# find all markers of cluster 2
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
gse756d_seurat.markers <- FindAllMarkers(gse756d_seurat, only.pos = TRUE)
gse756d_seurat.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)


cluster0.markers <- FindMarkers(gse756d_seurat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(gse756d_seurat, features = c("", ""))

# you can plot raw counts as well
VlnPlot(gse756d_seurat, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(gse756d_seurat, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))



gse756d_seurat.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(gse756d_seurat, features = top10$gene) + NoLegend()


new.cluster.ids <- c("")
names(new.cluster.ids) <- levels(gse756d_seurat)
gse756d_seurat <- RenameIdents(gse756d_seurat, new.cluster.ids)
DimPlot(gse756d_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



library(ggplot2)
plot <- DimPlot(gse756d_seurat, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)



table(gse756d_seurat@meta.data$orig.ident)


#now 3-16 columns are we want as the samples 
#this will be the dosage file 

#select_LN_cells=WhichCells(gse756d_seurat, idents = c('BC03LN','BC07LN'))
#gs756d_seurat_LN <- subset(gse756d_seurat, cells%in%select_LN_cells)_LN=subset(gse756d_seurat, idents = c('BC03LN','BC07LN')),cells=select_LN_cells)
#table(gse756d_seurat_LN@meta.data$orig.ident)

#control_list=unique(gse756d_seurat@meta.data$orig.ident)
#control_list

#control_list=c('BC01','BC02','BC03','BC03','BC03','BC03','BC03','BC03','BC03','BC03')



#control_list%in%control_list!%in%c('BC03LN','BC07LN')
#select_control_cells=WhichCells(gse756d_seurat,expression=orig.ident%in%)




saveRDS(gse756d_seurat, file = "C:/Users/bbula/OneDrive/Desktop/VS Code DS/DSCode/Seurat-Guided/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE756seurat.rds")
saveRDS(gse756d_seurat.markers, file = "C:/Users/bbula/OneDrive/Desktop/VS Code DS/DSCode/Seurat-Guided/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE756seuratmarkers.rds")


#Kevin's script 
names(gse756d_seurat@meta.data)
gse756d_seurat@meta.data$nCount=NULL
gse756d_seurat@meta.data$nCount_RNA=NULL
gse756d_seurat@meta.data$percent.mt=NULL


VlnPlot(gse756d_seurat, features = "nFeature_RNA")

# Normalize Data
gse756d_seurat <- NormalizeData(gse756d_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Variable Features
gse756d_seurat <- FindVariableFeatures(gse756d_seurat, selection.method = "vst", nfeatures = 500)
## Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(gse756d_seurat), 10)

top10 
## plot variable features with and without labels
plot1 <- VariableFeaturePlot(gse756d_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scale Data
all.genes <- rownames(gse756d_seurat)
gse756d_seurat <- ScaleData(gse756d_seurat, features = all.genes)

# PCA
gse756d_seurat <- RunPCA(gse756d_seurat, features = VariableFeatures(object = gse756d_seurat))

# Find Clusters
gse756d_seurat <- FindNeighbors(gse756d_seurat, dims = 1:10)
gse756d_seurat <- FindClusters(gse756d_seurat, resolution = 0.5)

# UMAP
gse756d_seurat <- RunUMAP(gse756d_seurat, dims = 1:10)
DimPlot(gse756d_seurat, reduction = "umap")

# Find marker genes
gse756d_seurat.markers <- FindAllMarkers(gse756d_seurat, only.pos = TRUE)
gse756d_seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


FeaturePlot(gse756d_seurat,features='ENSG00000127324.4')
FeaturePlot(gse756d_seurat,features="TSPAN8")



#Here, I am doing the work necessary in order to 
#categorize the cell markers into cell type labels
cluster0.markers <- FindMarkers(gse756d_seurat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers
head(cluster0.markers, n = 20)

#For 0, TSPAN8 is more toward intestinal tissues, 



gse756d_seurat.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
DoHeatmap(gse756d_seurat, features = top10$gene) + NoLegend()



levels(gse756d_seurat)


select_LN_cells=WhichCells(gse756d_seurat,expression=orig.ident%in%c('BC03LN','BC07LN'))
gse756d_seurat_LN=subset(gse756d_seurat,cells=select_LN_cells)
table(gse756d_seurat_LN@meta.data$orig.ident)

control_list=c('BC01','BC02','BC03','BC04','BC05','BC06','BC07','BC08','BC09','BC10','BC11')
select_control_cells=WhichCells(gse756d_seurat,expression=orig.ident%in%control_list)
gse756d_seurat_control=subset(gse756d_seurat,cells=select_control_cells)
table(gse756d_seurat_control@meta.data$orig.ident)


# Save your workspace to a file
save.image("my_workspace.RData3")
load("my_workspace.RData3")



save.image("my_workspace.RDataC")
load("my_workspace.RDataC")

# at this point, I need to select the genes from the original dataset,
# for the lymph node data, and then rerun the entire pipeline to 
# get a dataset for the lymph node stuff. At that point, we can
# run a gene set enrichment analysis between the two sets, to see
# if there are any significant differences between the two datasets.
# regarding the difference between the control and the experimental group 


#one thing to consider is; do I need to filter B3 and B7 out of the 
#control group (the non lymph node set)_? 

#Lets switch gears to the SeuratIntegration

#BiocManager::install("TFBSTools", type = "source", force = TRUE)
library(Matrix)
gse756d_seuratint <- gse756d_seurat

gse756d_seurat

gse756d_seuratint <- RunAzimuth(gse756d_seuratint, reference = "pbmcref")
gse756d_seuratint <- RunUMAP(gse756d_seuratint, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")

DimPlot(gse756d_seuratint, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))


#Caveat, the dimensionality ratio creates different plots based on the ratio
#the original ratio was 1:10, but the integrated analysis was 1:30, which yielded 
#different results. However, once we accounted for the dimension, the 1:10 yielded the same results. 



library(clusterProfiler)
search_kegg_organism('ece', by='kegg_code')

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)

mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)
head(mkk2)


browseKEGG(kk, 'hsa04110')








### Kevin Script P2


library(clusterProfiler)
library(tidyverse)

##### Trait Analysis: gene ~ LN_binary #####
# Create pseudobulk gene expresion objects
## Pseudobulk means 1 column as 1 subject; single cell data is 1 column as 1 cell
a <- read.table('C:/Users/bbula/OneDrive/Desktop/VS Code DS/DSCode/Seurat-Guided/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt', header=T)
pseudo_bulk <- a[, 4:17]


view(a)
tibble(a)
head(pseudo_bulk)
tibble(pseudo_bulk)
view(pseudo_bulk)


pseudo_bulk <- cbind(gene_names,pseudo_bulk)
dim(pseudo_bulk)
names(pseudo_bulk)

gene_names <- a[,1]


view(gene_names)


# Remove decimals from ENSEMBL ID
gene_names <- gsub("\\..*", "", gene_names)
gex <- a[,2:dim(pseudo_bulk)[2]]

gex2 <- pseudo_bulk[,2:15]
gex3 <- pseudo_bulk[,4-17]
view(gex2)
view(gex3)
# Create LN meta data (0 for no LN, 1, for LN)
LN_meta.data <- as.data.frame(c(0,0,0,0,1,0,0,0,0,1,0,0,0,0))


view(beta)
view(gex)
view(lm_df)

n <- length(gene_names)
result.i <- NULL
result <- NULL
for (i in 1:n){
  # For each gene in the data set, find if that gene is significantly different between 2 groups (control vs LN)
  gene <- gene_names[i]
  # Grab the expression values for the gene expression (gex) of gene "i"
  beta <- gex2[i,]
  # Create a data frame to input into the linear regression model (lm=linear model, df = dataframe)
  ## Rows are individuals, column 1 = gene expression, column 2 = binary LN
  lm_df <- cbind(t(beta),LN_meta.data)
  colnames(lm_df) <- c(gene,'LN_binary')
  fix.eff <- paste0(gene,' ~ LN_binary')
  fix.eff <- formula(fix.eff)
  #fit=glm(fix.eff,data=lm_df,family=binomial(link=logit))
  fit <- glm(fix.eff,data=lm_df)
  exposure <- 'LN_binary'
  result.anno <- gene # Specify column name for statistic summary table
  if (!inherits(fit, "try-error")) {
    BETA <- formatC(summary(fit)$coefficients[exposure,1],format="e",digits=2)
    STDERR <- formatC(summary(fit)$coefficients[exposure,2],format="e",digits=2)
    t <- formatC(summary(fit)$coefficients[exposure,3],format="e",digits=2)
    P <- formatC(summary(fit)$coefficients[exposure,4],format="e",digits=2)
    result.i <- cbind.data.frame(result.anno,exposure,BETA,STDERR,t,P)
    result <- as.data.frame(rbind(result,result.i))
  } else{
    BETA <- "NA"
    STDERR <- "NA"
    t <- "NA"
    P <- "NA"
    result.i <- cbind.data.frame(result.anno,exposure,BETA,STDERR,t,P)
    result <- as.data.frame(rbind(result,result.i))
  }
}

result$BETA <- as.numeric(result$BETA)
result$STDERR <- as.numeric(result$STDERR)
result$t <- as.numeric(result$t)
result$P <- as.numeric(result$P)


tibble(result$P)
# If you get an ERCC error, you can ignore it, the sumamry statistics data frame still comes out just fine

# Set P-value threshold to account for multiple testing
bonferroni_threshold <- 0.05/n
top_genes <- filter(result, P<bonferroni_threshold)
write.table(top_genes,'top_genes_summary_statistics2.txt')

##### Gene Pathway Analysis #####
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)

# Load in human gene set
hs <- org.Hs.eg.db
# Get a list of top genes
top_genes_ls <- top_genes$result.anno
tibble(top_genes)
tibble(top_genes_ls)
view(top_genes_ls)

view(top_genes)


ensembleP <- top_genes[1:85,]
view(ensembleP)

selected_rows <- top_genes[86:170, 1]
view(selected_rows)
## Here I had trouble replicating Kevin's CNET
## In the ensembl_entrez_mapping file, there is a second ENSG00000156273 between 29-30, with symbols BACH1 and GRIK1-AS2
#I see, I think the problem is that my genes are not being annotated properly, as many ENS IDS are resulting in errors



gnames <- a[,1:2]
view(gnames)
# Get the indices of the selected genes in the ENSG column
indices <- match(selected_rows, gnames$gene_id)

# Get the corresponding gene names
selected_gene_names <- gnames$gene_name[indices]
view(indices)





view(selected_gene_names)


###Reverse engineer to get ENTREZID
symindices <- match(selected_gene_names, gene.df$SYMBOL)
view(symindices)

selected_ENTREZ <- gene.df$ENTREZID[symindices]
view(selected_ENTREZ)



pannotate <- selected_gene_names
view(pannotate)
pannotate$ENTREZID <- selected_ENTREZ




view(ensembl_entrez_mapping)



head(top_genes_ls)
view(gene.df)

gene.df <- bitr(selected_gene_names, fromType = "SYMBOL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)


gene.df2 <- bitr(top_genes_ls, fromType = "ENSEMBL",
        toType = c("ENTREZID"),
        OrgDb = org.Hs.eg.db)


ego2 <- enrichGO(gene         = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENTREZID',
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(ego2)                





#top_genes_ls=gsub("\\..*", "", top_genes_ls)
# Annotate genes from symbol into entrezID
ensembl_entrez_mapping <- AnnotationDbi::select(hs,
                                                keys = top_genes_ls,
                                                keytype = "ENSEMBL",
                                                columns = c("ENTREZID", "SYMBOL"))




view(ensembl_entrez_mapping)
head(ensembl_entrez_mapping)


# Perform GO gene enrichment

view(ensembl_entrez_mapping$ENSEMBL)
view(ensembl_entrez_mapping$ENTREZID)

ggo <- enrichGO(gene=ensembl_entrez_mapping$ENTREZID,
             OrgDb=org.Hs.eg.db,
             ont='ALL',
             pAdjustMethod='fdr',
             pool=T,
             readable=T)


             
head(ggo)


head(bgo)
bgo <- enrichGO(gene         = ensembleP$result.anno,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                readable   = TRUE)  







ggoA <- enrichGO(gene=ensembl_entrez_mapping$ENSEMBL,
             OrgDb=org.Hs.eg.db,
             keyType = "ENSEMBL",
             ont='ALL',
             pAdjustMethod='fdr',
             pool=T,
             readable=T)


             


view(ggo)
view(ggoA)
head(ggoA)
# Create cnet plot (gene concept network plot)
edox <- setReadable(ggo,'org.Hs.eg.db','ENTREZ')
p1 <- cnetplot(edox)
p2 <- cnetplot(edox,circular=T,colorEdge=T)
#cowplot::plot_grid(p1,p2,ncol2=2,labels=LETTERS[1:3],rel_width=c(.8,.8,1.2))
p3 <- p1+p2
ggsave('cnet2.png',p3,height=8,width=14)

# Perform KEGG gene enrichment
kk=enrichKEGG(gene = ensembl_entrez_mapping$ENTREZID,
              organism = 'hsa',
              pvalueCutoff = 0.05)
# Extract KEGG pathway IDs
pathway_ids <- kk$ID
browseKEGG(kk,pathway_ids[1]) # Since KEGG only identified 1 cellular pathway, you can just set the index to 1







view(memes)
memes<- mapIds(org.Hs.eg.db, keys = selected_gene_names,
       column = "ENTREZID", keytype = "SYMBOL")

view(top_genes_ls)

library(AnnotationDbi)
library("org.Hs.eg.db")
#columns(org.Hs.eg.db) # returns list of available keytypes
ensembl_entrez_mapping$ENTREZID = mapIds(org.Hs.eg.db,
                    keys=ensembl_entrez_mapping$ENSEMBL, #Column containing Ensembl gene ids
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")


