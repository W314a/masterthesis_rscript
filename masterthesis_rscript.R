# Installation of all needed packages and libraries ----
# Packages
install.packages("Seurat") # used for single cell RNA sequencing analysis
install.packages("metap") # used for meta analysis of p-values
install.packages("devtools") # helps installing, developing and managing R packages
install.packages("remotes") # installing packages from GitHub or other remote sources
install.packages("cowplot") # enhances ggplot2 with better multi-panel plots
install.packages("gridExtra") # provides functions for arranging multiple plots on a grid
install.packages("BiocManager") # used to install and manage Bioconductor packages
BiocManager::install("EnhancedVolcano") # creates volcano plots for visualizing differential gene expression
BiocManager::install('multtest') # provides multiple hypothesis testing correction methods
BiocManager::install("clusterProfiler") # Used for functional enrichment analysis
# Libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(reticulate)
library(cowplot)
library(metap)
library(gridExtra)




# Colour palettes
sample_palette <- c("#d79001", "#02bf7d", "#9590ff")
heatmap_palette <- c("#423e9e", "white", "#f0a818")
featureplot_palette <- c("#423e9e", "grey", "#f0a818")
df_PC_palette <- c("IgG1" = "#f564e3", "IgG2b/c" = "#f8766d", "IgG3" = "#b79f00", "IgA 1" = "#02bfc4", "IgA 2" = "#609dff", "IgM" = "#00ba39")




# Commands applied to data beforehand ----
df@commands




# Preparing the data set ----
# Load the data set
df <- readRDS("/Users/pia/Documents/Ausbildung und Arbeit/Studium/Master/SS24/Analyses/Bioinformatics/scRNAseq Data/scManz.integrated_od_markers_annotated_24.10.2023.rds")
# Setting the resolution to 0.5 (can be set to 0.1, 0.3, 0.5 or 0.8)
Idents(df) <- "integrated_snn_res.0.5"



# Analysis of the original data set ----


## Quality control and UMAP ----
VlnPlot(df, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, 
        group.by = "Sample", cols = sample_palette)
# Visualiation of the dimentionality of the data
ElbowPlot(df, ndims = 100, reduction = "pca") + labs(x = "Principal Component", y = "Standard Deviation")
# PCA plot split by sample to compare the samples
DimPlot(df, reduction = "pca", group.by = "Sample", cols = sample_palette) + labs(x = "Principal Component 1", y = "Principal Component 2")
# UMAP Plot with resolution of 0.5
DimPlot(df, reduction = "umap", split.by = "Sample", label = TRUE) # UMAP for each sample to compare the clustering
DimPlot(df, reduction = "umap", label = TRUE) # UMAP with all samples combined


## Marker genes of the original data ----
# Finding markers
Markers <- FindAllMarkers(df)
# Top 10 DEGs
top10 <- Markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()
# Heatmap of the top10 DEGs
DoHeatmap(df, 
          features = top10$gene,
          size = 3, angle = 20) + theme(text = element_text(size = 9)) + guides(color = "none") + scale_fill_gradientn(colours = heatmap_palette)


## Expression of IgH genes ----
# Expression of all IgH genes visualised as violin plots for all 10 clusters
VlnPlot(df, 
        features = c("Ighm", "Ighd", "Igha", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3", "Ighe"), 
        same.y.lims = 7.5, pt.size = 0)
# Co-expression of Ighm and Ighd
FeaturePlot(df, 
            features = c("Ighm", "Ighd"), 
            blend = TRUE, blend.threshold = 0.5, # to show co-expression in the Feature Plot
            label = TRUE)
# Expression of Cd4 und Cd8
VlnPlot(df, 
        features = c("Cd4", "Cd8a", "Cd8b1"), 
        same.y.lims = 5, pt.size = 0)


## Removing immunoglobulin heavy chain genes ----
# Remove Ig heavy chain genes
df_no_Igs <- subset(df, 
                    features = setdiff(rownames(df), 
                                       c("Igha", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3", "Ighm", "Ighd", "Ighe")))
# renaming the Idents for clearer visualisation
df_no_Igs <- RenameIdents(df_no_Igs,
                          c("7" = "IgG1",
                            "1" = "IgG2b/c",
                            "3" = "IgG3",
                            "2" = "IgA",
                            "5" = "IgM"))
# Re-normalise and scale the data
df_no_Igs <- NormalizeData(df_no_Igs)
df_no_Igs <- ScaleData(df_no_Igs)
# Balloon plot for selected B-cell marker genes
B_cell_Markers <- DotPlot(df_no_Igs, 
                          features = c("Igkc", "Iglc2", "Iglc3",
                                       "Fcer2a","Cr2", "Ms4a1", "Cd19", "Ebf1",
                                       "Pou2af1","Xbp1", "Irf8", "Blnk", "Bach2", "Pax5", "Bcl6","Prdm1", "Fosb", "Irf4", "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17", 
                                       "Jchain")) & coord_flip()
B_cell_Markers <- B_cell_Markers + scale_color_gradient2(high = "#f0a818", low = "#423e9e", mid = "grey")
B_cell_Markers




# Analysis of the PC clusters ----


## Subset the original clusters for the identified PC clusters 1, 2, 3, 5 and 7 ----
df_PC <- subset(df, idents = c("1", "2", "3", "5", "7"))


## Reanalyse the subset (PCA, UMAP and clustering analysis) ----
df_PC <- FindVariableFeatures(df_PC)
df_PC <- ScaleData(df_PC)
df_PC <- RunPCA(df_PC)
df_PC <- RunUMAP(df_PC, dims = 1:50)
df_PC <- FindClusters(df_PC, 
                      resolution = .7, # sets the resolution to 0.7
                      algorithm = 4, # uses the Leiden algorithm
                      method = "igraph", 
                      graph.name = "integrated_snn")


## UMAP visualisation ----
DimPlot(df_PC, reduction = "umap")


## Marker gene expression ----
# Compare each cluster to all other clusters
Markers_PC <- FindAllMarkers(df_PC)
# find the top 1ß DEGs
top10_PC <- Markers_PC %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj) %>%
  slice_head(n = 10) %>%
  ungroup()
# Visualisation of the top 10 DEGs as a heatmap
DoHeatmap(df_PC, 
          features = top10_PC$gene,
          size = 4, angle = 20) + theme(text = element_text(size = 15)) + guides(color = "none") + scale_fill_gradientn(colours = heatmap_palette)
# Visualisation of the expression of IgH genes throughout the clusters
FeaturePlot(df_PC,
            features = c("Ighm", "Ighd", "Igha", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3"),
            label = TRUE, cols = featureplot_palette, ncol = 4)

# Renaming the clusters with their respective PC isotype
df_PC <- RenameIdents(df_PC, c("6" = "IgG1",
                               "1" = "IgG2b/c",
                               "2" = "IgG3",
                               "4" = "IgA 1",
                               "5" = "IgA 2",
                               "3" = "IgM"))



## Homing and survival molecules ----

### Visualisation as violin plots ----
VlnPlot(df_PC, features = c("Cxcr3", "Cxcr4", "Cxcr5", "Ccr7", "Ccr9", "Ccr10"), 
        pt.size = 0, y.max = 5, ncol = 3, cols = df_PC_palette)

### Visualisation as Volcano Plots ----

# Finding DEGs between the clusters
# Comparing each cluster with all other clusters
IgG1vsAll <- FindMarkers(df_PC, ident.1 = "IgG1", ident.2 = c("IgG2b/c", "IgG3", "IgA 1", "IgA 2", "IgM"))
IgG2bcvsAll <- FindMarkers(df_PC, ident.1 = "IgG2b/c", ident.2 = c("IgG1", "IgG3", "IgA 1", "IgA 2", "IgM"))
IgG3vsAll <- FindMarkers(df_PC, ident.1 = "IgG3", ident.2 = c("IgG1", "IgG2b/c", "IgA 1", "IgA 2", "IgM"))
IgA1vsAll <- FindMarkers(df_PC, ident.1 = "IgA 1", ident.2 = c("IgG2b/c", "IgG3", "IgM", "IgG1", "IgA 2"))
IgA2vsAll <- FindMarkers(df_PC, ident.1 = "IgA 2", ident.2 = c("IgG2b/c", "IgG3", "IgM", "IgG1", "IgA 1"))
IgMvsAll <- FindMarkers(df_PC, ident.1 = "IgM", ident.2 = c("IgG2b/c", "IgG3", "IgA 1", "IgA 2", "IgG1"))
# Comparing the two IgA PC clusters with one another
IgA1vsIgA2 <- FindMarkers(df_PC, ident.1 = "IgA 2", ident.2 = "IgA 1")
IgG1vsIgG <- FindMarkers(df_PC, ident.1 = "IgG1", ident.2 = c("IgG2b/c", "IgG3"))
IgG2bcvsIgG <- FindMarkers(df_PC, ident.1 = "IgG2b/c", ident.2 = c("IgG1", "IgG3"))
IgG3vsIgG <- FindMarkers(df_PC, ident.1 = "IgG3", ident.2 = c("IgG1", "IgG2b/c"))

#### Chemokines Receptor genes ----
ChemokineR_Genes <- c("Cxcr3", "Cxcr4", "Cxcr5", "Ccr7", "Ccr9", "Ccr10")

# IgA 1 vs PCs
# Subset to include only the genes of interest
IgA1vsAll_Chemokine_R <- IgA1vsAll[rownames(IgA1vsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA1_Chemo <- EnhancedVolcano(IgA1vsAll_Chemokine_R,
                              lab = rownames(IgA1vsAll_Chemokine_R), selectLab = rownames(IgA1vsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                              title = 'PC clusters and IgA 1 PCs', subtitle = '', caption = "") + theme_classic()

# IgA 2 vs PCs
# Subset to include only the genes of interest
IgA2vsAll_Chemokine_R <- IgA2vsAll[rownames(IgA2vsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA2_Chemo <- EnhancedVolcano(IgA2vsAll_Chemokine_R,
                              lab = rownames(IgA2vsAll_Chemokine_R), selectLab = rownames(IgA2vsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                              title = 'PC clusters and IgA 2 PCs', subtitle = '', caption = "") + theme_classic()


# IgM vs PCs
# Subset to include only the genes of interest
IgMvsAll_Chemokine_R <- IgMvsAll[rownames(IgMvsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgM_Chemo <- EnhancedVolcano(IgMvsAll_Chemokine_R,
                             lab = rownames(IgMvsAll_Chemokine_R), selectLab = rownames(IgMvsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                             x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                             title = 'PC clusters and IgM PCs', subtitle = '', caption = "") + theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(""))

# IgG1 vs PCs
# Subset to include only the genes of interest
IgG1vsAll_Chemokine_R <- IgG1vsAll[rownames(IgG1vsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG1_Chemo <- EnhancedVolcano(IgG1vsAll_Chemokine_R,
                              lab = rownames(IgG1vsAll_Chemokine_R), selectLab = rownames(IgG1vsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                              title = 'PC clusters and IgG1 PCs', subtitle = '', caption = "") + theme_classic()

# IgG2b/c vs PCs
# Subset to include only the genes of interest
IgG2bcvsAll_Chemokine_R <- IgG2bcvsAll[rownames(IgG2bcvsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG2bc_Chemo <- EnhancedVolcano(IgG2bcvsAll_Chemokine_R,
                                lab = rownames(IgG2bcvsAll_Chemokine_R), selectLab = rownames(IgG2bcvsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1, FCcutoff = 1,
                                title = 'PC clusters and IgG2b/c PCs', subtitle = '', caption = "") + theme_classic() + 
  theme(legend.title = element_blank()) 

# IgG3vs PCs
# Subset to include only the genes of interest
IgG3vsAll_Chemokine_R <- IgG3vsAll[rownames(IgG3vsAll) %in% ChemokineR_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG3_Chemo <- EnhancedVolcano(IgG3vsAll_Chemokine_R,
                              lab = rownames(IgG3vsAll_Chemokine_R), selectLab = rownames(IgG3vsAll_Chemokine_R), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1, 
                              title = 'PC clusters and IgG3 PCs', subtitle = '', caption = "") + theme_classic() 

# All Chemokine plots together
Chemo_Plots <- cowplot::plot_grid(IgG1_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim (-2, 115),
                                  IgG2bc_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim (-2, 115),
                                  IgG3_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim (-2, 115),
                                  IgA1_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 115), 
                                  IgA2_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 115),
                                  IgM_Chemo + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim (-2, 115),
                                  align = "vh")
# legend
legend_plot <- cowplot::get_legend(IgG2bc_Chemo)
# Combine legend plot and Volcano Plots
cowplot::plot_grid(Chemo_Plots, legend_plot, rel_widths = c(20, 4))



#### Survival factor genes ----
Survival_Genes <- c("Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17", "Il6ra")

# IgA 1 vs PCs
# Subset to include only the genes of interest
IgA1vsAll_Survival <- IgA1vsAll[rownames(IgA1vsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA1_Survival <- EnhancedVolcano(IgA1vsAll_Survival,
                                 lab = rownames(IgA1vsAll_Survival), selectLab = rownames(IgA1vsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                 x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                                 title = 'PC clusters and IgA 1 PCs', subtitle = '', caption = "") + theme_classic()

# IgA 2 vs PCs
# Subset to include only the genes of interest
IgA2vsAll_Survival<- IgA2vsAll[rownames(IgA2vsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA2_Survival <- EnhancedVolcano(IgA2vsAll_Survival,
                                 lab = rownames(IgA2vsAll_Survival), selectLab = rownames(IgA2vsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                 x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                                 title = 'PC clusters and IgA 2 PCs', subtitle = '', caption = "") + theme_classic()


# IgM vs PCs
# Subset to include only the genes of interest
IgMvsAll_Survival <- IgMvsAll[rownames(IgMvsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgM_Survival <- EnhancedVolcano(IgMvsAll_Survival,
                                lab = rownames(IgMvsAll_Survival), selectLab = rownames(IgMvsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                                title = 'PC clusters and IgM PCs', subtitle = '', caption = "") + theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(""))

# IgG1 vs PCs
# Subset to include only the genes of interest
IgG1vsAll_Survival <- IgG1vsAll[rownames(IgG1vsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG1_Survival<- EnhancedVolcano(IgG1vsAll_Survival,
                                lab = rownames(IgG1vsAll_Survival), selectLab = rownames(IgG1vsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                                title = 'PC clusters and IgG1 PCs', subtitle = '', caption = "") + theme_classic()

# IgG2b/c vs PCs
# Subset to include only the genes of interest
IgG2bcvsAll_Survival<- IgG2bcvsAll[rownames(IgG2bcvsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG2bc_Survival<- EnhancedVolcano(IgG2bcvsAll_Survival,
                                  lab = rownames(IgG2bcvsAll_Survival), selectLab = rownames(IgG2bcvsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                  x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1, 
                                  title = 'PC clusters and IgG2b/c PCs', subtitle = '', caption = "") + theme_classic()

# IgG3vs PCs
# Subset to include only the genes of interest
IgG3vsAll_Survival <- IgG3vsAll[rownames(IgG3vsAll) %in% Survival_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG3_Survival <- EnhancedVolcano(IgG3vsAll_Survival,
                                 lab = rownames(IgG3vsAll_Survival), selectLab = rownames(IgG3vsAll_Survival), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                                 x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1,
                                 title = 'PC clusters and IgG3 PCs', subtitle = '', caption = "") + theme_classic() 

# All Chemokine plots together
Survival_Plots <- cowplot::plot_grid(IgG1_Survival + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 200),
                                     IgG2bc_Survival + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 200),
                                     IgG3_Survival+ theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 200),
                                     IgA1_Survival + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 200),
                                     IgA2_Survival + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 200),
                                     IgM_Survival + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 210),
                                     align = "vh")
# legend
legend_plot <- cowplot::get_legend(IgG2bc_Chemo)
# Combine legend plot and Volcano Plots
cowplot::plot_grid(Survival_Plots, legend_plot, rel_widths = c(20, 4))


#### Integrin and additional genes ----
Int_Add_Genes <- c("Itga4", "Itgb2", "Cd44", "Cd28", "Cd37")

# IgA 1 vs PCs
# Subset to include only the genes of interest
IgA1vsAll_Int <- IgA1vsAll[rownames(IgA1vsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA1_Int <- EnhancedVolcano(IgA1vsAll_Int,
                            lab = rownames(IgA1vsAll_Int), selectLab = rownames(IgA1vsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                            x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                            title = 'PC clusters and IgA 1 PCs', subtitle = '', caption = "") + theme_classic()

# IgA 2 vs PCs
# Subset to include only the genes of interest
IgA2vsAll_Int<- IgA2vsAll[rownames(IgA2vsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA2_Int <- EnhancedVolcano(IgA2vsAll_Int,
                            lab = rownames(IgA2vsAll_Int), selectLab = rownames(IgA2vsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                            x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                            title = 'PC clusters and IgA 2 PCs', subtitle = '', caption = "") + theme_classic()


# IgM vs PCs
# Subset to include only the genes of interest
IgMvsAll_Int <- IgMvsAll[rownames(IgMvsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgM_Int <- EnhancedVolcano(IgMvsAll_Int,
                           lab = rownames(IgMvsAll_Int), selectLab = rownames(IgMvsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                           x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, 
                           title = 'PC clusters and IgM PCs', subtitle = '', caption = "") + theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(""))

# IgG1 vs PCs
# Subset to include only the genes of interest
IgG1vsAll_Int <- IgG1vsAll[rownames(IgG1vsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG1_Int<- EnhancedVolcano(IgG1vsAll_Int,
                           lab = rownames(IgG1vsAll_Int), selectLab = rownames(IgG1vsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                           x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                           title = 'PC clusters and IgG1 PCs', subtitle = '', caption = "") + theme_classic()

# IgG2b/c vs PCs
# Subset to include only the genes of interest
IgG2bcvsAll_Int<- IgG2bcvsAll[rownames(IgG2bcvsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG2bc_Int<- EnhancedVolcano(IgG2bcvsAll_Int,
                             lab = rownames(IgG2bcvsAll_Int), selectLab = rownames(IgG2bcvsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                             x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1, 
                             title = 'PC clusters and IgG2b/c PCs', subtitle = '', caption = "") + theme_classic()

# IgG3vs PCs
# Subset to include only the genes of interest
IgG3vsAll_Int <- IgG3vsAll[rownames(IgG3vsAll) %in% Int_Add_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG3_Int <- EnhancedVolcano(IgG3vsAll_Int,
                            lab = rownames(IgG3vsAll_Int), selectLab = rownames(IgG3vsAll_Int), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                            x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 1,
                            title = 'PC clusters and IgG3 PCs', subtitle = '', caption = "") + theme_classic() 

# All Chemokine plots together
Int_Add_Plots <- cowplot::plot_grid(IgG1_Int + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    IgG2bc_Int + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    IgG3_Int+ theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    IgA1_Int + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    IgA2_Int + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    IgM_Int + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 80),
                                    align = "vh")
# legend
legend_plot <- cowplot::get_legend(IgG2bc_Chemo)
# Combine legend plot and Volcano Plots
cowplot::plot_grid(Int_Add_Plots, legend_plot, rel_widths = c(20, 4))


#### IgA 1 vs IgA 2 ----
Homing_Genes <- c("Cxcr3", "Cxcr4", "Cxcr5", "Ccr7", "Ccr9", "Ccr10", "Tnfrsf13c", "Tnfrsf13b", "Tnfrsf17", "Il6ra", "Itga4", "Itgb2", "Cd44", "Cd28", "Cd37")

# Subset to include only the genes of interest
IgA1vsIgA2_Homing <- IgA1vsIgA2[rownames(IgA1vsIgA2) %in% Homing_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgA_Homing <- EnhancedVolcano(IgA1vsIgA2_Homing,
                              lab = rownames(IgA1vsIgA2_Homing), selectLab = rownames(IgA1vsIgA2_Homing), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                              title = 'IgA 1 and IgA 2 PCs', subtitle = '', caption = "") + theme_classic()
IgA_Homing + theme(legend.position = "right", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4)



#### IgG-isotypes comparison ----

# IgG1 vs IgG
# Subset to include only the genes of interest
IgG1vsIgG_Homing <- IgG1vsIgG[rownames(IgG1vsIgG) %in% Homing_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG1_Homing <- EnhancedVolcano(IgG1vsIgG_Homing,
                               lab = rownames(IgG1vsIgG_Homing), selectLab = rownames(IgG1vsIgG_Homing), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                               x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                               title = 'IgG1 and IgG PCs', subtitle = '', caption = "") + theme_classic()
IgG1_Homing + theme(legend.position = "right", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4)

# IgG2bc vs IgG
# Subset to include only the genes of interest
IgG2bcvsIgG_Homing <- IgG2bcvsIgG[rownames(IgG2bcvsIgG) %in% Homing_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG2bc_Homing <- EnhancedVolcano(IgG2bcvsIgG_Homing,
                              lab = rownames(IgG2bcvsIgG_Homing), selectLab = rownames(IgG2bcvsIgG_Homing), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                              x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                              title = 'IgG2b/c and IgG PCs', subtitle = '', caption = "") + theme_classic()
IgG2bc_Homing + theme(legend.position = "right", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4)

# IgG3 vs IgG
# Subset to include only the genes of interest
IgG3vsIgG_Homing <- IgG3vsIgG[rownames(IgG3vsIgG) %in% Homing_Genes, ]
# Plot the volcano plot using EnhancedVolcano
IgG3_Homing <- EnhancedVolcano(IgG3vsIgG_Homing,
                               lab = rownames(IgG3vsIgG_Homing), selectLab = rownames(IgG3vsIgG_Homing), drawConnectors = TRUE, arrowheads = FALSE, labSize = 4, pointSize = 3,
                               x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05,
                               title = 'IgG3 and IgG PCs', subtitle = '', caption = "") + theme_classic()
IgG3_Homing + theme(legend.position = "right", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4)

# IgG Plots
IgG_Plots <- cowplot::plot_grid(IgG2bc_Homing + theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 10),
                                IgG3_Homing+ theme(legend.position = "none", plot.title = element_text(size = 12), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14)) + xlim(-4, 4) + ylim(-2, 10), 
                                align = "vh")
# legend
legend_plot <- cowplot::get_legend(IgG2bc_Chemo)
# Combine legend plot and Volcano Plots
cowplot::plot_grid(IgG_Plots, legend_plot, rel_widths = c(20, 4))
