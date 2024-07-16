# Decoding the immune checkpoint signatures in human atherosclerosis: implications of type 2 diabetes and lipid-lowering intervention
DOI: add paper DOI

***

## Human carotid plaque scRNA-seq analysis
Raw data and processed count matrices for this study have been deposited in GEO under the accession numbers [GSE246315](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246315), [GSE235437](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235437), and [GSE224273](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224273)

[Code and associated files](Human_carotid_plaque_scRNAseq/)

### Description of analysis pipeline and associated files:
1.	[Create seurat object](Human_carotid_plaque_scRNAseq/R_scripts/create_seurat.R)
      Read in cellranger generated count matrices and generates a list of seurat objects 
2.  [CITE-seq sample processing](Human_carotid_plaque_scRNAseq/R_scripts/CITEseq_GSM7502475_Sample38_filter_process.R)
	    Process and QC the single CITE-seq sample (GSM7502475_Sample38)
3.  [Filter and QC](Human_carotid_plaque_scRNAseq/R_scripts/filter_process.Rmd)
			QC and filtering of low quality barcodes for each sample and exports a list of filtered seurat objects
4.  [Integration and reference annotation](Human_carotid_plaque_scRNAseq/R_scripts/annotate_integrate.R)
	    Map each sample to [CITE-seq PBMC reference](https://doi.org/10.1016/j.cell.2021.04.048) and integrate
5.  [Clustering](Human_carotid_plaque_scRNAseq/R_scripts/cluster_all_cells.R)
      Dimensionality reduction, clustering, marker gene identification for each cluster
6.  [Subset major cell types](Human_carotid_plaque_scRNAseq/R_scripts/Subset_major_celltypes.R)
      Subset major cell types (myeloid, NK, B, T, CD4 T, CD8 T, DN T, DP T)
7.  [Subclustering](Human_carotid_plaque_scRNAseq/R_scripts/subclustering/)
      Dimensionality reduction, clustering, and marker gene identification of cells in each major population
8.  [Reclassification of DNT cells](Human_carotid_plaque_scRNAseq/R_scripts/DNT_clustering_diagnostics.Rmd)
      T-cell clustering diagnostics and reclassification of double negative T-cells (DNT) as CD4 or CD8
9.  [Celltypist](Human_carotid_plaque_scRNAseq/celltypist/celltypist_run.py)
      Annotate cells with [celltypist](https://www.celltypist.org/)
10. [Cellchat](Human_carotid_plaque_scRNAseq/cellchat/)
      Run [cellchat](https://github.com/jinworks/CellChat)
      
      
***  


## Human PBMC CITE-seq analysis
Raw data and processed count matrices for this study have been deposited in GEO under the accession number [GSE246317](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246317)

[Code and associated files](Human_PBMC_CITE-seq/)

### Description of analysis pipeline and associated files:
1.	[Pre-processing and filtering of 4h post-treatment samples](Human_PBMC_CITE-seq/PBMC_CITEseq_analysis_preprocessing_4h.Rmd)
      Read in cellranger generated count matrices, inspect quality, and filter low quality barcodes
2.	[Pre-processing and filtering of 24h post-treatment samples](Human_PBMC_CITE-seq/PBMC_CITEseq_analysis_preprocessing_24h.Rmd)
      Read in cellranger generated count matrices, inspect quality, and filter low quality barcodes
3.	[Integration and clustering](Human_PBMC_CITE-seq/anti_CTLA4_PD1_integrate.Rmd)
      Integrate, and cluster anti-CTLA4 and anti-PD1 samples. Identify marker genes for each cluster
4.  [Protein data clustering](Human_PBMC_CITE-seq/Cluster_and_annotate_protein_data.Rmd)
      Use protein expression and clustering to annotate WNN clusters with major cell type identities
5.  [T-cell gating](Human_PBMC_CITE-seq/subclustering/T_cells_wnn_subcluster.Rmd)
      Divide T-cells into CD4, CD8, DN T, and DP T, based on expression of CD4 and CD8 protein
6.  [Subclustering](Human_PBMC_CITE-seq/subclustering/)
      Dimensionality reduction, clustering, and marker gene identification of cells in each major population
      
      
***



## Human PBMC scRNA-seq analysis
Raw data and processed count matrices for this study have been deposited in GEO under the accession number [GSE272294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272294)

[Code and associated files](Human_PBMC_scRNAseq/)

### Description of analysis pipeline and associated files:
1.	[Create seurat object and filter](Human_PBMC_scRNAseq/R_scripts/preprocess_filter.R)
      Read in cellranger generated count matrices, filter low quality barcodes, and save a list of seurat objects
2.	[Mark multiplets](Human_PBMC_scRNAseq/R_scripts/detect_doublets.R)
      Mark multiplets using [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder/)
3.	[Remove proliferation associated genes](Human_PBMC_scRNAseq/R_scripts/Remove_cc_genes.R)
      Remove cell-cycle/proliferative genes for integration and clustering
4.	[Integration](Human_PBMC_scRNAseq/R_scripts/Integrate.R)
      Integrate data
5.  [Dimensionality reduction and clustering](Human_PBMC_scRNAseq/R_scripts/Cluster.R)
      Remove marked multiplets and cluster
6.	[Celltypist](Human_PBMC_scRNAseq/R_scripts/celltypist/celltypist_run.py)
      Annotate cells with [celltypist](https://www.celltypist.org/)
7.  [Differential expression analysis](Human_PBMC_scRNAseq/R_scripts/Differential_expression.R)
      Marker gene identification for each cluster. Differential expression analysis of immune checkpoint genes between type 2 and no diabetes NT samples for each cell type. The latter was run after subclustering and annotation
8.  [Subset major cell types](Human_PBMC_scRNAseq/R_scripts/Add_celltypist_labels_&_subset_cells.R)
      Add celltypist labels and subset major cell types (myeloid, NK, B, T, CD4 T, CD8 T, DN T, DP T)
9.  [Subclustering](Human_PBMC_scRNAseq/R_scripts/subclustering/)
      Dimensionality reduction, clustering, and marker gene identification of cells in each major population
10. [Cellchat](Human_PBMC_scRNAseq/cellchat/)
      Run [cellchat](https://github.com/jinworks/CellChat)
      
      
***   

## Mouse scRNA-seq analysis
Raw data and processed count matrices for this study have been deposited in GEO under the accession numbers [GSE272294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272294), [GSE141038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141038), [GSE161494](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161494), [GSE168389](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168389), and [GSE253555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253555)

[Code and associated files](Mouse_plaque_scRNAseq/)

### Description of analysis pipeline and associated files:
1.	[Process and integrate](Mouse_plaque_scRNAseq/Process_integrate.R)
      Read in cellranger generated count matrices, filter low quality barcodes, integrate, cluster, identify marker genes per cluster
2.  [Identify proliferation associated genes](Mouse_plaque_scRNAseq/Cellcycle_genes.R)
      Identify cell-cycle/proliferative genes for removal prior to integration and clustering
3.  [Differential expression analysis](Mouse_plaque_scRNAseq/FindMarkersCondition.Rmd)
      Differential expression analysis of immune checkpoint genes between NT and LL. The latter was run after subclustering and annotation
4.  [Cellchat](Mouse_plaque_scRNAseq/cellchat/)
      Run [cellchat](https://github.com/jinworks/CellChat)

***

## Human coronary scRNA-seq analysis
Raw data and processed count matrices for this study have been deposited in GEO under the accession number [GSE264666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264666)

[Code and associated files](Human_coronary_scRNAseq/)

### Description of analysis pipeline and associated files:
1.	[Create seurat object and filter](Human_coronary_scRNAseq/R_scripts/QC)
      Read in cellranger generated count matrices, filter low quality barcodes, detect and mark multiplets via [scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
2.	[Integrate](Human_coronary_scRNAseq/R_scripts/Pre_process_integrate.Rmd)
      Remove multiplets, integrate, cluster, and identify marker genes per cluster. Annotate CD45- clusters and subcluster CD45+ cells
3.	[Subset major cell types](Human_coronary_scRNAseq/R_scripts/CD45+_analysis&subsetting.Rmd)
      Subset major cell types (myeloid, NK, B, T, CD4 T, CD8 T, DN T, DP T) from CD45+ cells
4.  [Subcluster](Human_coronary_scRNAseq/R_scripts/subclustering)
      Dimensionality reduction, clustering, and marker gene identification of cells in each major population
5.  [Plaque vs control](Human_coronary_scRNAseq/R_scripts/Plaque_vs_control.Rmd)
      Differential expression analysis of plaque vs control