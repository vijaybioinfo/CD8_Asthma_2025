# CD8_asthma_2025

## Requirements
For scRNA-seq data the project has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.1.3)  
• [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) (v7.1.0)  
• [Scrublet](https://github.com/swolock/scrublet/blob/master/README.md) (v0.2.3)  
• [Seurat](https://satijalab.org/seurat/) (v4.3.0)  
• [vdjtools](https://vdjtools-doc.readthedocs.io/en/master/) (v1.2.1)  
• [Cellchat](https://github.com/jinworks/CellChat) (v2.2.0)

Regarding Spatial transcriptomics analysis, data has been processed following modules/programs:  
• [R](https://cran.r-project.org/) (v4.3.3)  
• [Cell Ranger](https://www.10xgenomics.com/support/software/xenium-ranger/latest)) (v4.0.1.1)  
• [Seurat](https://satijalab.org/seurat/) (v5.1.0)  
• [phenoptr](https://akoyabio.github.io/phenoptr/) (v0.3.2)

## Raw data access
The single-cell RNA-seq raw and processed files from scRNA-seq and spatial transcriptomics data can be downloaded through the following GEO accession number:  

## Raw data pre-processing 
For more specific information about the data generation and processing, please check the "methods" section within the manuscript.  

> Relevant scripts for both 5'v2 scRNA-seq and Spatial Transcriptomics data are all located in: ./pre-processing 

### Single-cell RNA-seq
• To do the 10x demultiplexing and mapping pull [our in-house pipeline](https://github.com/vijaybioinfo/cellranger_wrappeR) using Cell Ranger.  
• To do the donor demultiplexing pull [our in-house pipeline](https://github.com/vijaybioinfo/ab_capture).  
• To do the single-cell quality control pull [our in-house pipeline](https://github.com/vijaybioinfo/quality_control).  
• To do the doublet detection use [our in-house pipeline](https://github.com/vijaybioinfo/doublet_detection) using Scrublet.  
• To generate the clustering of single-cell data just pull [our in-house pipeline](https://github.com/vijaybioinfo/clustering) using Seurat.  
• To do the aggregation of VDJ libraries pull [our in-house pipeline](https://github.com/vijaybioinfo/VDJ_aggr).  

### Spatial Transcriptomics
• To generate the clustering of Spatial Transcriptomics data just pull [our in-house pipeline](https://github.com/vijaybioinfo/Xenium-clustering/) using Seurat.  

## Figures 
> Relevant scripts are located in: ./figures

DGEA - You can follow [our in-house pipeline](https://github.com/vijaybioinfo/dgea)  
CellChat - You can use our in-house cellchat implementation.  

## Usage & Citation
If you want to clone this repository run:  
> git clone https://github.com/vijaybioinfo/Asthma_2025.git  


## Maintainers
Current maintainers:  

Francisco Emmanuel Castañeda-Castro (fcastaneda@lji.org)  
Donaldo Sosa-Garcia (donaldo@lji.org)  
Vijayanand Lab.  
Division of Vaccine Discovery La Jolla Institute for Immunology La Jolla, CA 92037, USA  

## Contact
Please email Francisco Emmanuel Castaneda-Castro (fcastaneda@lji.org), Donaldo Sosa-Garcia (donaldo@lji.org) and/or Vijayanand Pandurangan (vijay@lji.org).
