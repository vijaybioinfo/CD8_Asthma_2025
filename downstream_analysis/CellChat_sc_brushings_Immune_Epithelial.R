# brushings_immune
    library(CellChat)
    library(patchwork)
    library(dplyr)
    options(stringsAsFactors = FALSE)
    # Step 1: load seurat objects 

    #remove bad quality cells from brushings immune
    brushings_immune2 <- subset(brushings_immune, orig.celltype_tcells !=  "Tcells_Badquality")

    brushings_immune2$orig.celltype_cellchat <- paste0("I_", as.character(brushings_immune2$orig.celltype_tcells))
    table(brushings_immune2$orig.celltype_cellchat, brushings_immune2$orig.celltype_tcells)

    brushings_epi2 <- brushings_epi

    brushings_epi2$orig.celltype_cellchat <- paste0("S_", as.character(brushings_epi2$orig.celltype))
    table(brushings_epi2$orig.celltype, brushings_epi2$orig.celltype_cellchat)

    brushings_immune2$all <- "all_immune"
    brushings_epi2$all <- "all_epi"
    
    setwd("/home/donaldo/bioadhoc/SaHe/SaHe_bru/results/cellchat/brushings_cellchat/")

    meta_brushings_immune2 <- brushings_immune2@meta.data[,c("orig.donor_paper", "orig.disease", "orig.celltype_cellchat", "UMAP_1", "UMAP_2", "origlib")]
    # brushings_trms$orig.celltype <- NULL
    
    meta_brushings_epi2 <- brushings_epi2@meta.data[, c("orig.donor_paper", "orig.disease", "orig.celltype_cellchat", "UMAP_1", "UMAP_2", "origlib")]

    brushings_immune2@meta.data <- meta_brushings_immune2
    
    brushings_immune2 <- RenameCells(brushings_immune2, new.names = paste0("immune-", Cells(brushings_immune2)))

    brushings_epi2@meta.data <- meta_brushings_epi2

    brushings_epi2 <- RenameCells(brushings_epi2, new.names = paste0("epi-", Cells(brushings_epi2)))

    # The overlap should be minimal
    table(colnames(brushings_immune2) %in% colnames(brushings_epi2), meta_brushings_immune2$origlib, useNA="always")

    table(colnames(brushings_epi2) %in%  colnames(brushings_immune2), meta_brushings_epi2$origlib, useNA="always")

    brushings_cellchat <- merge(brushings_immune2, brushings_epi2, add.cells.ids = c("immune", "epi"), project="brushings_cellchat_immune_all_epi")
    nrow(brushings_cellchat@meta.data)
    length(unique(rownames(brushings_cellchat@meta.data)))
    
    brushings_cellchat$samples <- factor(brushings_cellchat$orig.donor_paper, levels=unique(brushings_cellchat$orig.donor_paper))

    table(brushings_cellchat$orig.donor_paper, useNA="always")
    table(brushings_cellchat$orig.disease)
    table(brushings_cellchat$orig.celltype_cellchat, brushings_cellchat$orig.disease)

    outdir <- "/home/donaldo/bioadhoc/SaHe/SaHe_bru/results/cellchat/brushings_cellchat"
 
    trm_colors_cellchat <- c(
    'I_CD8' = '#DF594C',
    'I_CD4' = '#F5C04A',
    'I_Int_macro' = '#DEC6F5',
    'I_Class_macro' = '#A06BFF',
    'I_MAST' = '#D6A57C',
    'I_Bcells' = '#00A9FF',
    'I_NK' = '#313C5C',
    'I_pDC' = '#4CAF50',
    'S_Ciliated' = '#F4A07A',
    'S_Suprabasal' = '#7AC8D9',  
    'S_Basal' = '#F7D678',        
    'S_Secretory' = '#F194C4',
    'S_Ionocytes' = '#9BE3D7'  
    )

  trm_colors_cellchat2 <- c(
        'I_Bcells' = '#bfbfbf',
        'I_Class_macro' = '#bfbfbf',
        'I_Int_macro' = '#bfbfbf',
        'I_MAST' = '#bfbfbf',
        'I_NK' = '#bfbfbf',
        'I_pDC' = '#bfbfbf',
        'I_CD8' = '#DF594C',
        'I_CD4' = '#bfbfbf',
        'S_Basal' = '#F7D678',    
        'S_Ciliated' = '#F4A07A',
        'S_Ionocytes' = '#9BE3D7',
        'S_Secretory' = '#F194C4',
        'S_Suprabasal' = '#7AC8D9')
    
    
    source("/home/donaldo/bioadhoc/SaHe/SaHe_bru/scripts/cellchat/cellchat_v1.R") #cellchat function for prcessing
    
    running_cellchat(
      object_cellchat = brushings_cellchat,
      pop_size = TRUE, 
      pval_cells=0.05,
      fc_cells=0,
      mean_method = list("truncatedMean", 0.1),
      dgea_method = list("cellchat"),
      comparing_groups = c("Severe.asthma", "Mild.asthma", "Healthy"),
      donor_info = "orig.donor_paper",
      outdir=outdir,
      senders_celltypes=c('I_CD8', 'I_CD4', 'I_Int_macro', 'I_Class_macro', 'I_MAST', 'I_Bcells','I_NK','I_pDC'),
      receivers_celltypes=c("S_Secretory", "S_Ciliated", "S_Basal", "S_Suprabasal", "S_Ionocytes"),
      celltype_selected = 'I_CD8',
      colors_celltypes=trm_colors_cellchat,
      color_edges=trm_colors_cellchat2,
      include_non_protein_signalling=FALSE
      )
