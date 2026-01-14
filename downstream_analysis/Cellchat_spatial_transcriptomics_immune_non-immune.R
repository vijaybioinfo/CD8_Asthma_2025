# Cellchat analysis in Spatial Transcriptomic data

# Set annotations for cellchat
  SaHe_all_slide2@meta.data$orig.celltype_cellchat <- SaHe_all_slide2$orig.celltype_tc_dc
  table(SaHe_all_slide2$Xenium_snn_res.0.6, SaHe_all_slide2@meta.data$orig.celltype_cellchat, useNA="always")
  table(SaHe_all_slide2@meta.data$orig.celltype_cellchat, useNA="always")

  source_cellchat <- c('CD8', 'CD4', 'monocytes', 'Mast', 'bcells', 'DCs', 'NK') 

  #colors for visualizations
  pop_cl_colors <- c('basal' = '#F7D678',
                      'endothelial' = '#3792c1',
                      'fibroblasts' = '#24cc87',
                      'smooth_muscle' = '#00c3af',
                      'ciliated' = '#F4A07A',
                      'Tcells' = '#C0332C',
                      'secretory1' = '#F194C4',
                      'monocytes' = '#A06BFF',
                      'DCs' = '#DEC6F5',
                      'bcells' = '#00A9FF',
                      'secretory2' = '#c92e8c',
                      'pericytes' = '#f4c17aff',
                      'Mast' = '#D6A57C',
                      'schwann' = '#333eb2ff',
                      'Chondrocyte' = '#98a164ff',
                      'CD8' = '#DF594C',
                      'CD4' = '#F5C04A',
                      'NK' = '#313C5C')

  outdir <- "/home/donaldo/bioadhoc/spatial/SaHe/slide_5k_SaHe2/re_segmentation/results/figures/cellchat_spatial/slide2_new_seg_no_imputation_corrected/"
  dir.create(outdir)

  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/scripts/cellchat_v1.R")

  running_cellchat(
    object_cellchat = SaHe_all_slide2,
    pop_size = TRUE ,
    pval_cells = 0.05,
    fc_cells = 0,
    mean_method = list("truncatedMean", 0.05),
    dgea_method = list("cellchat"),
    comparing_groups = c("Severe.asthma", "Mild.asthma", "Healthy"),
    outdir = outdir,
    senders_celltypes = source_cellchat,
    celltype_selected = "CD8",
    receivers_celltypes = names(pop_cl_colors)[!names(pop_cl_colors) %in% source_cellchat & !names(pop_cl_colors) == 'Tcells'],
    colors_celltypes=pop_cl_colors,
    include_non_protein_signalling=FALSE,
    column_group= "orig.disease",
    identities_clusters = "orig.celltype_cellchat",
    exp_type = list(type="Spatial", params="distance.use = TRUE, raw.use = TRUE, contact.range = 10, contact.dependent = TRUE, interaction.range = 250, scale.distance = 0.5"),
    donor_info = "orig.donor"
  )
