# Basic single cell processing using Seurat 


#### Table of Contents
1. [Quality control](#quality-control)<br/>
2. [Clustering and dimensionality reduction](#clustering-and-dimensionality-reduction)<br/>
3. [Integration](#integration)<br/>
4. [Marker Identification](#marker-identification)<br/>
5. [Cell type annotation](#cell-type-annotation)<br/>
6. [Differential gene expression analysis](#differential-gene-expression-analysis)<br/>
7. [Pathway enrichment](#pathway-enrichment)<br/>

#### Notes

This workflow in its current form is built on Seurat functions with some of the below functions serving as wrappers for Seurat functions. Seurat vignettes should be explored before running this to learn about the different steps.<br/>
https://satijalab.org/seurat/articles/get_started.html<br/>
https://satijalab.org/seurat/articles/integration_introduction.html

This is a living document that will be modified and added to.

## Quality control

### Create Seurat object from 10X output

In the current working directory, each sample has its own directory containing barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz from 10X.
```R
library(Seurat)
source("Source/scFunctions.R")
setwd("~/Documents/Single_Cell_Pipeline_Example/data/")
samples <- list.dirs("./", full.names = FALSE)
samples = samples[samples!=""]
```

This function assumes sample directory names are formatted genotype-diet-tissue-replicate in the below example. (AD-F-hypo-1, must be separated by dash)
```R
seuratObject <- CreateSeuratObjectFrom10X(sample_dirs = samples,
                                          meta_data = c("genotype","diet",
                                                        "tissue","replicate"))
```

Set working directory to project folder (preceding directory)
```R
setwd("../")
```
### Check data quality and filter cells
```R
# if working with human data, change the appropriate prefixes
seuratObject = getCellQuality(seuratObject = seuratObject,
                              feature_patterns=list("percent.mito"="^mt-",
                                                    "percent.ribo"=c("^Rps","^Rpl"),
                                                    "percent.pred"=c("^Gm1","^Gm2","^Gm3","^Gm4","^Gm5","^Gm6","^Gm7","^Gm8","^Gm9"),
                                                    "percent.Hb"=c("^Hba-","^Hbb-","^Hbq")))
cellQualityPlot(seuratObject=seuratObject,
                fileName="./Plots/QCPlots/ViolinPlots/PreFilter/HYP_2020.pdf",
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito",
                               "percent.ribo","percent.pred","percent.Hb"),
                identPlot = "orig.ident",
                H=9,W=20,
                pointSize=0)
```
![Prefilter QC plot](https://github.com/jessicading/single-cell-utils/blob/master/Figures/HYP_2020_Pre_QC.png?raw=true)
```R
# thresholds can be changed as needed
filtered_sample = subset(x = seuratObject, 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < 7500 & 
                           nCount_RNA > 500 &
                           nCount_RNA < 35000 &
                           percent.mito < .15 & 
                           percent.ribo < .6 & # not really needed
                           percent.pred < .10 & 
                           percent.Hb < 0.05)
# redo cell quality plot to see QC after filtering
cellQualityPlot(seuratObject=filtered_sample,
                fileName=paste0("./Plots/QCPlots/ViolinPlots/PostFilter/HYP_2020.pdf"),
                H=9,W=20,
                featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito",
                               "percent.ribo","percent.pred","percent.Hb"),
                identPlot = "orig.ident",
                pointSize=0)
```
![Prefilter QC plot](https://github.com/jessicading/single-cell-utils/blob/master/Figures/HYP_2020_Post_QC.png?raw=true)


## Clustering and dimensionality reduction
First cluster cells without integration (batch correction). In order words, cluster using the RNA assay. This function will create plots in the Plots directory of the current working directory. Please review parameters of this function and adjust accordingly. Defaults should
work fine for most cases.

```R
filtered_sample <- processSeurat(seuratObject=filtered_sample, 
                                 name="AD_HYP_2020", # should be tissue or descriptive name
                                 plot_meta=c("orig.ident", "diet"))
# each cell now has cluster assignments at different resolutions
grep("^RNA",colnames(filtered_sample@meta.data), value = TRUE)
```
```
[1] "RNA_snn_res.0.5" "RNA_snn_res.1"   "RNA_snn_res.1.5" "RNA_snn_res.2"   "RNA_snn_res.2.5"
[6] "RNA_snn_res.3"   "RNA_snn_res.3.5" "RNA_snn_res.4" 
```
```R
dir.create("Saves")
saveRDS(filtered_sample, "Saves/RNA_Clust_AD_ADF_2020.rds")
```

This plot is generated by `ClusterAnalysis()` used below. It is shown here for explanation purposes. The individual plots from `processSeurat()` can be seen in Plots/rnaUMAP/resolution/PCA_UMAPres0.5.pdf, Plots/rnaUMAP/AD_HYP_2020_diet.pdf, and AD_HYP_2020_orig.ident.pdf. 
![RNA res 0.5 results](https://github.com/jessicading/single-cell-utils/blob/master/Figures/AD_HYP_2020_RNA_res.0.5_umap.png?raw=true)
Here we see separation of fructose treated and untreated cells for some clusters. For example, cluster 3 and 5 look to be the same cell type but are designated as different clusters. For this reason, we choose to run integration analysis to get clusters representing conserved cell types across conditions. If a condition/treatment is hypothesized to have an effect, ideally in the RNA clustering, there should be large "shifts" in the cells (in the reduction plots) and replicates of the same condition are overlapping (consistent).

## Integration
(Optional) Integration is highly recommended such as in the case described above.
```R
# save umap and tsne reductions from RNA clustering to add to integrated
# seuratObject (if applicable) since it will be overwritten by Integration
rna_pca_red <- filtered_sample@reductions$pca
rna_umap_red <- filtered_sample@reductions$umap
rna_tsne_red <- filtered_sample@reductions$tsne

filtered_sample = runSampleCCA(seuratObject=filtered_sample,
                               combineLevel="orig.ident", # which "batch" you want to correct for
                               # can leave these parameters as is
                               numFeatures=2000,numDims=30,kParam=25,
                               # change for your specific analysis
                               tissue="AD_HYP_2020",
                               resUse=seq(0.5,4,by=0.5))
# each cell now has integrated cluster assignments at different resolutions
grep("^integrated",colnames(filtered_sample@meta.data), value = TRUE)
```
```
[1] "integrated_snn_res.0.5" "integrated_snn_res.1"   "integrated_snn_res.1.5"
[4] "integrated_snn_res.2"   "integrated_snn_res.2.5" "integrated_snn_res.3"  
[7] "integrated_snn_res.3.5" "integrated_snn_res.4" 
```
```R
# Write integrated clustering to new slots
filtered_sample@reductions$IntegrateSamplePCA <- filtered_sample@reductions$pca
filtered_sample@reductions$IntegrateSampleUMAP <- filtered_sample@reductions$umap
filtered_sample@reductions$IntegrateSampleTSNE <- filtered_sample@reductions$tsne
# Add RNA clustering back to seuratObject (gets overwritten from integration)
filtered_sample@reductions$pca <- rna_pca_red
filtered_sample@reductions$umap <- rna_umap_red
filtered_sample@reductions$tsne <- rna_tsne_red
saveRDS(filtered_sample, "Saves/Integrated_AD_ADF_2020.rds")
```
Optionally, plot integration cluster distribution in RNA clustering and vice versa to observe where the integrated clusters were projected on in the RNA clustering
```R
ClusterDistribitionAcrossMethods(seuratObject = filtered_sample,
                                 # can change to higher resolutions
                                 clustering1="RNA_snn_res.0.5",
                                 clustering2="integrated_snn_res.0.5", 
                                 reduction1 = "IntegrateSampleUMAP",
                                 reduction2 = "umap",
                                 output_name="Plots/HYP_RNA_Integrated_Cluster_Distribution.png",
                                 plot_titles_size = 20)
```
![RNA and Integrated cluster distribution](https://github.com/jessicading/single-cell-utils/blob/master/Figures/HYP_RNA_Integrated_Cluster_Distribution.png?raw=true)

### Run Harmony
Harmony is a much faster batch correction method that performs similarly to Seurat Integration. 
```R
# run harmony on non integrated, did not work on integrated for some reason
harmonized <- readRDS("Saves/RNA_Clust_AD_ADF_2020.rds")
harmonized <- RunHarmony(harmonized, 
                         "orig.ident", # can be other metadata and multiple metadata
                         plot_convergence = TRUE)
harmonized <- RunUMAP(harmonized, reduction = "harmony", dims = 1:30)
# add harmony reduction to filtered_sample (main seurat)
filtered_sample@reductions$harmony_umap <- harmonized@reductions$umap
harmonized <- FindNeighbors(harmonized, reduction = "harmony", dims = 1:30)
# save RNA clustering clusters
harmonized <- FindClusters(harmonized, resolution = 0.5) 
filtered_sample$harmony_snn_res.0.5 <- harmonized$RNA_snn_res.0.5
# compare against Seurat integration
ggarrange(DimPlot(filtered_sample, group.by = "Sample", reduction = "IntegrateSampleUMAP") +
            ggtitle("Seurat Integration") + theme(title = element_text(size = 25)),
          DimPlot(filtered_sample, group.by = "integrated_snn_res.0.5", reduction = "IntegrateSampleUMAP", 
                  label = TRUE, label.size = 8)+ggtitle("Seurat Integration") + 
            NoLegend()+theme(title = element_text(size = 25)),
          DimPlot(filtered_sample, group.by = "Sample", reduction = "harmony_umap") +
            ggtitle("Harmony") + theme(title = element_text(size = 25)),
          DimPlot(filtered_sample, group.by = "harmony_snn_res.0.5", reduction = "harmony_umap",
                  label = TRUE, label.size = 8)+ggtitle("Harmony") + NoLegend() + 
            theme(title = element_text(size = 25)), 
          ncol = 4)
```
From a quick qualitative assessment, Seurat Integration may have performed slightly better. 
![Harmony and Seurat Integration](https://github.com/jessicading/single-cell-utils/blob/master/Figures/Seurat_Integration_Harmony_Comparison.png?raw=true)

Plot to visualize the consistency between the two methods. A minor inconsistency includes Seurat Integration cluster 0 to be split between Harmony cluster 0 and 15 (and cluster 15 is comprised of a single sample).
```
ClusterDistribitionAcrossMethods(seuratObject = filtered_sample,
                                 # can change to higher resolutions
                                 clustering1="integrated_snn_res.0.5",
                                 clustering2="harmony_snn_res.0.5", 
                                 reduction1="harmony_umap", 
                                 reduction2="IntegrateSampleUMAP", 
                                 output_name="Plots/HYP_Integrated_Harmony_Cluster_Distribution.png",
                                 plot_titles_size = 20)
```
![Harmony and Integrated cluster distribution](https://github.com/jessicading/single-cell-utils/blob/master/Figures/HYP_Integrated_Harmony_Cluster_Distribution.png?raw=true)

## Marker Identification
```R
# Just adding Sample meta column
filtered_sample$Sample <- paste0(filtered_sample$genotype, 
                                 filtered_sample$diet, 
                                 filtered_sample$replicate)
filtered_sample$Sample <- gsub("W","",filtered_sample$Sample)
```
This function plots visualizations of cluster identity and calculates marker genes and visualizes them in a combined plot (output in Figures directory). Number of markers to show is 10 and number of clusters to show on one image is 5 and these can be changed by setting the `nMarkers` and `nClustersOneImage` parameters. This function returns a Seurat object with markers in `filtered_sample@misc$markers`.
```R
filtered_sample <- ClusterAnalysis(seuratObject = filtered_sample, 
                                   condition = c("diet","Sample"), # Sample can be orig.ident
                                   sample_colors = NULL, # can pass in named vector where names are samples and values are colors
                                                         # see function definition for example
                                   show_clusters = c("integrated_snn_res.0.5", # can pass in just one
                                                     "RNA_snn_res.0.5"),
                                   name="AD_HYP_2020",
                                   sample_proportion_plot_width=10)
```
integrated_snn_res.0.5 clustering results shown below.
![UMAP integrated res 0.5](https://github.com/jessicading/single-cell-utils/blob/master/Figures/AD_HYP_2020_integrated_res.0.5_UMAP.png?raw=true)

![integrated res 0.5 sample proportions](https://github.com/jessicading/single-cell-utils/blob/master/Figures/AD_HYP_2020_integrated_res.0.5_SampleProportions.png?raw=true)

![markers plot example](https://github.com/jessicading/single-cell-utils/blob/master/Figures/AD_HYP_2020_integrated_res.0.5_UMAP_Markers_Set1.png?raw=true)

```R
# Additional markers can be calculated and plotted using this function.
# There may be biologically relevant clusters at higher resolutions so it's good to check.
filtered_sample <- makeMarkerFeaturePlots(seuratObject=filtered_sample,
                                          calculate_clusters=c("integrated_snn_res.1",
                                                               "integrated_snn_res.1.5"),
                                          output_name="Figures/Markers/AD_HYP_2020_markers.png",
                                          red=c("IntegrateSampleUMAP"))
# save seuratObject with markers
saveRDS(filtered_sample, "Saves/Integrated_AD_ADF_2020.rds")   

# this function can also be called using precalculated markers by setting the "markers"
# parameter to just make marker plots
integrated_res.0.5_markers <- filtered_sample@misc$markers[filtered_sample@misc$markers$res=="integrated_snn_res.0.5",]
makeMarkerFeaturePlots(seuratObject = filtered_sample, 
                       markers = integrated_res.0.5_markers, 
                       clusters_show = list("integrated_snn_res.0.5"=c(0,2,5,7,11)), # Optional, can choose to visualize
                                                                                     # only select clusters
                       output_name="Figures/Markers/AD_HYP_2020_Select_Cluster_Markers.png",
                       red = "IntegrateSampleUMAP", 
                       png_resolution_factor = 2) # Optional, twice the resolution
```
(Optional but recommended) Find conserved markers. Conserved markers are not so much a concern with single cell data sets with very distinct cell types and "successful" integration analysis, but it is still good to run to check. For subclustering, it is more highly recommended.
```R
Idents(seuratObject) <- "integrated_snn_res.0.5"
conserved_markers <- data.frame()
for(clust in unique(seuratObject$integrated_snn_res.0.5)){
  temp <- FindConservedMarkers(seuratObject,
                               ident.1 = clust, 
                               grouping.var = "orig.ident")
  temp$gene <- rownames(temp)
  temp$cluster <- clust
  conserved_markers <- rbind(conserved_markers, temp)
}
```

## Cell type annotation
Shown is a prediction of cell types by simple marker overlap from reference markers which can be either from Seurat `FindMarkers()` output or a two column module gene file with cell type in module column and gene in gene column.
```R
cluster_predictions <- markerGeneOverlap(query_markers=integrated_res.0.5_markers,                                        
                                         reference_markers="Markers/Hypothalamus_WT_markers_2019.txt",
                                         overlap_plot_wid=8,
                                         overlap_plot_hgt=5)
head(cluster_predictions$predictions)
```
```
  Cluster ByOverlapPrediction ByEnrichmentPrediction ByFDRPrediction
1       0           Astrocyte              Astrocyte       Astrocyte
2       1        Endothelium1           Endothelium2    Endothelium1
3       2     Oligodendrocyte        Oligodendrocyte Oligodendrocyte
4       3            Ependyma               Ependyma        Ependyma
5       4           Microglia             Macrophage       Microglia
6       5              Neuron                 Neuron          Neuron
```
Overlap results are contained in `cluster_predictions$overlap_results`.

Generated in InitialCellTypeNaming/ <br/>
![cell type prediction](https://github.com/jessicading/single-cell-utils/blob/master/Figures/heatmap.png?raw=true)

Other resources:<br/>
   https://panglaodb.se/search.html<br/>
   https://singlecell.broadinstitute.org/single_cell<br/>
   http://biocc.hrbmu.edu.cn/CellMarker/search.jsp

Plotting a cluster tree may also help with allocating clusters to certain cell types. Subtypes can be explored in subclustering analysis. You may also want to explore trees using higher resolution clusters.
```R
Idents(filtered_sample) <- "integrated_snn_res.0.5"
filtered_sample <- BuildClusterTree(object = filtered_sample)
pdf("Figures/Cluster_DimPlot_Tree.pdf")
DimPlot(filtered_sample, 
        group.by = "integrated_snn_res.0.5", 
        label = TRUE, label.size = 5) + NoLegend()
PlotClusterTree(object = filtered_sample)
dev.off()
```

The reference markers used were from a very similar dataset (and same tissue) so I feel comfortable relying on the overlap results (and prior experience) to label cell types. In other cases, more research is needed (or just calling the cluster by its top marker). This can be a very time consuming step.
```R
View(cluster_predictions$predictions)
```
The "ByOverlapPrediction" is a pretty good estimation of the cell type except for cluster 11 where the "ByEnrichmentPrediction" seems to be a more accurate cell type label (Tanycyte instead of Ependyma). Generally all three prediction methods should agree on the same cell type but review those that do not.
```R
predictions <- cluster_predictions$predictions
new.cluster.ids <- c(predictions$ByOverlapPrediction[1:11],
                     predictions$ByEnrichmentPrediction[12],
                     predictions$ByOverlapPrediction[13:19])
names(new.cluster.ids) <- levels(filtered_sample)
filtered_sample <- RenameIdents(filtered_sample, new.cluster.ids)
filtered_sample$Cell_type <- filtered_sample@active.ident
DimPlot(filtered_sample, reduction = "IntegrateSampleUMAP", label = TRUE)
saveRDS(filtered_sample, "Saves/Integrated_AD_ADF_2020.rds")
```
![UMAP Cell types](https://github.com/jessicading/single-cell-utils/blob/master/Figures/Cell_type_annotation.png?raw=true)

## Differential gene expression analysis
```R
filtered_sample$condition <- gsub("[0-9]","",filtered_sample$Sample)
DEGs <- runDEGs(seuratObject = filtered_sample, 
                comparisons = list("F effect"=c("ADF","AD")),
                Cell_type_column="Cell_type",
                Condition_column="condition", 
                logfc_threshold = 0.1 # I usually set this to 0 for plotting aesthetics but it will be slow
                )
# more comparisons can be added to the comparisons parameter like so:
comparisons = list("F effect"=c("ADF","AD"),
                   "DHA effect"=c("ADDHA","AD"))
```                
### Run sample-wise DEGs (Highly recommended)
DEG analysis implemented in Seurat does not account for multiple samples/replicates in condition analysis. DEGs that are consistent across samples are ideal. This can help with choosing which genes to discuss/validate. `runSampleWiseDEGs()` is similar to `FindConservedMarkers()`.

```R
Idents(filtered_sample) <- "Cell_type"
SampleWise_DEGs <- runSampleWiseDEGs(seuratObject = filtered_sample,
                                     comparisons = list("F effect"=c("ADF","AD")),
                                     Cell_type_column = "Cell_type",
                                     Condition_column = "condition",
                                     sampleNameMetaEntry = "Sample", 
                                     logfc_threshold=0.05)
Global_Sample_Wise_DEGs_Compared <- compareGlobalAndSampleWiseDEGs(global_DEGs = DEGs, 
                                                                   sample_wise_DEGs = SampleWise_DEGs,
                                                                   output_name = "Global_SampleWise_DEG_Comparison.txt")
filtered_DEGs <- sampleWiseFilteredDEGsConsistent(global_DEGs = DEGs, 
                                                  sample_wise_DEGs = SampleWise_DEGs, 
                                                  nTotalSamples = 6,
                                                  nSampleSignificance = 5, # one less than total
                                                                           # can do 6 to be even more stringent
                                                  fdr_threshold = 0.1, 
                                                  lfc_threshold = .1)
DEGs <- addSampleWiseConsistency(DEG_df = DEGs, 
                                 stringentDEG_df = filtered_DEGs, 
                                 match_columns = c("Cell_type") # if there are multiple comparisons
                                                                 # and tissues, add the columns here
                                  ) 
sum(DEGs$SampleWiseConsistent=="Yes")==nrow(filtered_DEGs)

DEGs_Sample_Wise_Info <- merge(global_DEGs, 
                               SampleWise_DEGs,
                               by = c("Cell_type","GENE","Comparison"))
```

### Summarize and visualize DEGs
See scUtils/DEGSummaryWorkflow for more examples.

```R
DEG_df <- addInfoAndTrim(DEG_df = DEGs, 
                         FDR_cutoff = 0.05, 
                         # if analyzing multiple comparisons, put Comparison last (for the summarizeTopDegs function if you
                         # are adding the mean normalized gene expression)
                         meta_data = c("Cell_type","Comparison"), # can put just one, i.e. "Cell_type"
                         lfc_cutoff = .25)
meanNormExprDf <- makeMeanNormExprDf(seuratObject = filtered_sample, 
                                     meta_data=c("Cell_type","condition"), 
                                     genes = unique(DEG_df$GENE))
summarizedDEGs <- summarizeTopDegs(DEG_df = DEG_df, MODULE = "MODULE_direction", 
                                   comparisons = c("ADF.v.AD"="F effect"), # names of comparisons must match names of conditions list
                                   conditions = list("ADF.v.AD"=c("ADF","AD")),
                                   meanNormExprDf = meanNormExprDf,
                                   convertToHuman = TRUE,
                                   celltypes = unique(DEG_df$Cell_type), # don't need to set this is already have Cell_type column
                                   gene_info = "Resources/10090.protein.info.v11.0.txt",
                                   GWAS_info = "Resources/gwas-association-downloaded_2021-01-27-EFO_1001870-withChildTraits.tsv", 
                                   output_name="DEGs/AD_2020_Fructose_Effect_DEG_summary.xls") # must be .xls
                                   
# Plot top consistent DEGs
summary <- summarizedDEGs[["F effect"]]
DEGs$Tissue <- "Hypothalamus"
# Visualize top consistent DEGs
DEGs_logFCMinMax <- DEGs
# so that large DEGs do not dominate
DEGs_logFCMinMax$avg_log2FC <- MinMax(DEGs_logFCMinMax$avg_log2FC, 
                                     min = -1, max = 1) 
gene_dot_plot(DEG_df = DEGs_logFCMinMax, 
              vars = c("GENE","Cell_type"),
              var_list = list("GENE"= summary$GENE[!grepl("^Rpl|^Rps",summary$GENE)][1:30]), 
              cluster = c("Cell_type"), # for hierarchical clustering. if multiple tissues, can
                                        # choose "GENE" to make less hairy
              facet_by = "Comparison") # helpful for multiple comparisons                                   
```
![Cell type consistent](https://github.com/XiaYangLabOrg/scUtils/blob/master/Basic_pipeline/DEGs/Ct_Consistent_DEGs.png?raw=true)

```R
# Plot celltype specific DEGs
specificDEGs <- summary[summary$nCellTypes==1,]
# This table generation to identify "truly" specific DEGs
SpecificDEGsLogFCs <- MakeSpecificDEGsLogFCsDf(specificDEGs = specificDEGs, 
                                               DEG_df = DEGs, # use unfiltered DEGs for comprehensive data  
                                               column = "Cell_type", 
                                               output_name = "DEGs/AD_2020_Fructose_Ct_Specific_DEGs.xls")
# pick ct specific DEGs
ct_specific_DEGs = c()
cts <- names(SpecificDEGsLogFCs)[!(grepl("lfc|fdr",names(SpecificDEGsLogFCs)))]
for(ct in cts){
  specific_DEGs <- SpecificDEGsLogFCs[[paste0(ct,"_fdr_lfc_filtered")]]$GENE
  if(length(specific_DEGs)<3){
    nGenes = length(specific_DEGs)
  } else if(length(specific_DEGs)==0){
    next
  } else {
    nGenes = 3
  }
  ct_specific_DEGs <- c(ct_specific_DEGs, specific_DEGs[1:nGenes])
}
ct_specific_DEGs <- ct_specific_DEGs[!is.na(ct_specific_DEGs)]

# not setting cluster as I want to retain order as inputted
gene_dot_plot(DEG_df = DEGs_logFCMinMax, 
              vars = c("GENE","Cell_type"),
              var_list = list("Cell_type"=cts,
                              "GENE"= ct_specific_DEGs),  
              facet_by = "Comparison")
```
![Cell type specific](https://github.com/XiaYangLabOrg/scUtils/blob/master/Basic_pipeline/DEGs/Ct_Specific_DEGs.png?raw=true)


## Pathway enrichment
Functions for pathway enrichment of either a two-column "MODULE" "GENE" file or output from FindMarkers from Seurat using two-column "module" "gene" pathway resource data.

### Example Uses

```R

# Examples ------------------------------------------------------------------------------------------
# two column data, with group name in MODULE column and genes in GENE column
load("./Final_Microglia_clustering/Slingshot/HP_pseudotime_groups.RData")
colnames(pseudotime_clusters) <- c("MODULE","GENE") 
# returns comprehensive pathway enrichment results, can be used in ggplot2
pathway_df <- makePathwayEnrichmentDf(DEG_df = pseudotime_clusters, 
                        resources_path = "../Resources/", # containing two-column "module" "gene" .txt pathway data
                        output_Dir = "./PseudotimeGroupPathwayEnrichment", # must be new for each different enrichment analysis!
                        convertToHuman = TRUE) # set FALSE if query genes already human

# Seurat FindMarkers output 
# Must have genes in "GENE" column and if not setting the 'MODULE_column' in makePathwayEnrichmentDf, cell type in "Cell_type" column
logFC <- read.delim("./DEGs/HP_HYP_DEGs_new.txt")
logFC <- logFC[logFC$Tissue=="Hippocampus",]
logFC <- logFC[logFC$Cell_type %in% c("Astrocyte","Microglia","Pericyte","Neuron","Oligodendrocyte","Endothelium"),]

# returns comprehensive pathway enrichment results with summed logFC values for each pathway 
pathway_df <- makePathwayEnrichmentDf(DEG_df = logFC, 
                                      resources_path = "../Resources/", 
                                      output_Dir = "./DEGPathways", 
                                      convertToHuman = TRUE, logFC_threshold = .1,
                                      addlogFC = TRUE)

# setting the 'MODULE_column' because have another level of data (tissue)
logFC <- read.delim("./DEGs/HP_HYP_DEGs_new.txt")
logFC <- logFC[logFC$Cell_type %in% c("Astrocyte","Microglia","Pericyte","Neuron","Oligodendrocyte","Endothelium"),]
logFC$Ct_tissue <- paste(logFC$Cell_type, logFC$Tissue, sep = "_")
pathway_df <- makePathwayEnrichmentDf(DEG_df = logFC, 
                                      resources_path = "../Resources/", 
                                      output_Dir = "./DEGPathwaysTwoTissues", 
                                      convertToHuman = TRUE, logFC_threshold = .1,
                                      addlogFC = TRUE, 
                                      MODULE_column = "Ct_tissue")
```

