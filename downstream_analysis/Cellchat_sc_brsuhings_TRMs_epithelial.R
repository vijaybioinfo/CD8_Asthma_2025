library(CellChat)
    library(patchwork)
    library(dplyr)
    options(stringsAsFactors = FALSE)

    bru_trms_wo4$orig.celltype_cellchat <- paste0("I_", as.character(bru_trms_wo4$RNA_snn_res.0.4))
    table(bru_trms_wo4$orig.celltype_cellchat, bru_trms_wo4$RNA_snn_res.0.4)

    brushings_epi2 <- brushings_epi

    brushings_epi2$orig.celltype_cellchat <- paste0("S_", as.character(brushings_epi2$orig.celltype))
    table(brushings_epi2$orig.celltype, brushings_epi2$orig.celltype_cellchat)

    bru_trms_wo4$all <- "all_TRMs"
    brushings_epi2$all <- "all_epi"
    
    setwd("/home/donaldo/bioadhoc/SaHe/SaHe_bru/results/cellchat/brushings_cellchat_TRM/")

    meta_bru_trms_wo4 <- bru_trms_wo4@meta.data[,c("orig.donor_paper", "orig.disease", "orig.celltype_cellchat", "UMAP_1", "UMAP_2", "origlib")]
    # brushings_trms$orig.celltype <- NULL
    
    meta_brushings_epi2 <- brushings_epi2@meta.data[, c("orig.donor_paper", "orig.disease", "orig.celltype_cellchat", "UMAP_1", "UMAP_2", "origlib")]

    bru_trms_wo4_cc <- bru_trms_wo4

    bru_trms_wo4_cc@meta.data <- meta_bru_trms_wo4
    
    bru_trms_wo4_cc <- RenameCells(bru_trms_wo4_cc, new.names = paste0("TRM-", Cells(bru_trms_wo4_cc)))

    brushings_epi2@meta.data <- meta_brushings_epi2

    brushings_epi2 <- RenameCells(brushings_epi2, new.names = paste0("epi-", Cells(brushings_epi2)))

    # The overlap should be minimal
    table(colnames(bru_trms_wo4_cc) %in% colnames(brushings_epi2), meta_bru_trms_wo4$origlib, useNA="always")

    table(colnames(brushings_epi2) %in%  colnames(bru_trms_wo4_cc), meta_brushings_epi2$origlib, useNA="always")

    brushings_cellchat <- merge(bru_trms_wo4_cc, brushings_epi2, add.cells.ids = c("TRM", "epi"), project="brushings_cellchat_TRM_all_epi")
    nrow(brushings_cellchat@meta.data)
    length(unique(rownames(brushings_cellchat@meta.data)))
    
    brushings_cellchat$samples <- factor(brushings_cellchat$orig.donor_paper, levels=unique(brushings_cellchat$orig.donor_paper))

    table(brushings_cellchat$orig.donor_paper, useNA="always")
    table(brushings_cellchat$orig.disease)
    table(brushings_cellchat$orig.celltype_cellchat, brushings_cellchat$orig.disease)

    outdir <- "/home/donaldo/bioadhoc/SaHe/SaHe_bru/results/cellchat/brushings_cellchat_TRM"
 
      trm_colors_cellchat <- c(
        'I_0' = '#91D0A3',
        'I_1' = '#9B59B6',
        'I_2' = '#E27D60',
        'I_3' = '#F6C85F',
        'S_Basal' = '#F7D678',                  
        'S_Ciliated' = '#F4A07A',
        'S_Ionocytes' = '#9BE3D7',
        'S_Secretory' = '#F194C4',              
        'S_Suprabasal' = '#7AC8D9')
      
      trm_colors_cellchat2 <- c(
        'I_0' = '#bfbfbf',
        'I_1' = '#9B59B6',
        'I_2' = '#bfbfbf',
        'I_3' = '#bfbfbf',
        'S_Basal' = '#F7D678',                  
        'S_Ciliated' = '#F4A07A',
        'S_Ionocytes' = '#9BE3D7',
        'S_Secretory' = '#F194C4',              
        'S_Suprabasal' = '#7AC8D9')
    
    source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/scripts/cellchat_v1.R") #cellchat function for prcessing
    
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
      senders_celltypes=c('I_0', 'I_1', 'I_2', 'I_3'),
      receivers_celltypes=c("S_Secretory", "S_Ciliated", "S_Basal", "S_Suprabasal", "S_Ionocytes"),
      celltype_selected = 'I_1',
      colors_celltypes=trm_colors_cellchat,
      color_edges=trm_colors_cellchat2,
      include_non_protein_signalling=FALSE
      )
