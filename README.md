# CD8_asthma_2025

## Requirements
For scRNA-seq data the project has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.1.3)  
• [Cell Ranger](https://cran.r-project.org/) (v3.1.0)  
• [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)  
• [Seurat](https://satijalab.org/seurat/) (v4.3.0)  
• [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) (v1.2.1)  
• [Cellchat](https://github.com/jinworks/CellChat) (v2.2.0)  

Regarding Spatial transcriptomics analysis, data has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.3.3)  
• [Cell Ranger](https://cran.r-project.org/) (v3.1.0)  
• [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)  
• [Seurat](https://satijalab.org/seurat/) (v5.1.0)  
• [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) (v1.2.1)

## Raw data access
The single-cell RNA-seq raw and processed files from scRNA-seq and spatial transcriptomics data can be downloaded through the following GEO accession number:  

## Raw data pre-processing
# Single-cell RNA-seq
• To do the 10x demultiplexing and mapping pull our in-house pipeline using Cell Ranger.  
• To do the donor demultiplexing pull our in-house pipeline.  
• To do the single-cell quality control pull our in-house pipeline.  
• To do the doublet detection use our in-house pipeline using Scrublet.  
• To generate the clustering of single-cell data just pull our in-house pipeline using Seurat.  
• To do the aggregation of VDJ libraries pull our in-house pipeline.  
> Relevant scripts all located in: ./downstream_analysis

# Spatial Transcriptomics
• To generate the clustering of Spatial Transcriptomics data just pull our in-house pipeline using Seurat.  
> Relevant scripts all located in: ./downstream_analysis  
For more specific information about the data generation and processing, please check the "methods" section within the manuscript.  

