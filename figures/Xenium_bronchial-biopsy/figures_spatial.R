
  #  module load R/4.3.3
  #  R

  .libPaths(c("/mnt/biohome/donaldo/R/for_el8_kernel/4.3", "/home/fcastaneda/R/x86_64-pc-linux-gnu-library/4.3"))


resources = c(
  "/home/ciro/scripts/handy_functions/devel/file_reading.R", # readfile
  "/home/ciro/scripts/handy_functions/devel/utilities.R",
  # filters_columns is.file.finished show_commas (cluster_reports)
  "/home/ciro/scripts/clustering/R/plotting.R", # cluster_reports
  "/home/ciro/scripts/clustering/R/utilities.R", # get_top_n
  "/home/ciro/scripts/handy_functions/devel/filters.R", # sample_even
  "/home/ciro/scripts/handy_functions/devel/plots.R", # plot_pct getlegend mytheme
  "/home/ciro/scripts/handy_functions/R/stats_summary_table.R",
  "/home/fcastaneda/bin/clustering/R/plotting.R" #plots of clustering pipeline
)
for(i in resources){ source(i) }

  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R")
  # source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  source("/home/fcastaneda/bin/clustering/R/stats_summary_table.R")
  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/gpr25/gsea_signature.R")
  source("/home/ciro/scripts/handy_functions/R/gsea_tests.R") # gsea_matrix, gsea_plot_summary, gsea_process_list")
  source("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/gpr25/gsea_signature.R") # clean_feature_list signature_scoring
  # source("/home/ciro/scripts/handy_functions/R/stats_summary_table.R")
  source("/home/fcastaneda/bin/clustering/R/stats_summary_table.R")
  source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
  source("/home/ciro/scripts/figease/figease.R")
  library(Seurat)
  # library(vcd)
  library(ggridges)
  library(ggplot2)
  library(dplyr)
  library(phenoptr)
  library(matrixStats)
  library(ggplot2)
  library(ggforce)
  library(stringr)
  library(rstatix)
  library(gridExtra)
  library(ggpubr)
  library(reshape2)
  library(pheatmap)
  library(dplyr)
  library(devtools)
  library(scGSVA)
  library(scales)
  # library(scales)

  Image_plot_e<- function(
         object = SaHe_all_slide2_2_temp,
          fov = "fov",
          group.by = "",
          cols = pop_cl_colors,
          size = pop_cl_size_temp[SaHe_all_slide2_2_temp$orig.celltype_cellchat],
          border.size=NA,
          shuffle.cols = TRUE,
          coord.fixed=TRUE,
          crop=TRUE,
          alpha = 1,
          overlap = FALSE,
          split.by = NULL,
          outdir_image = paste0(out_dir1, "colon_clusters_v3_size1.pdf")
  ){        
              object$column_to_color <- object@meta.data[,group.by]
              cells <- Cells(x = object)
              fov <- fov %||% DefaultFOV(object = object)
              fov <- Filter(f = function(x) {
                  return(x %in% Images(object = object) && inherits(x = object[[x]],
                      what = "FOV"))
              }, x = fov)

              boundaries <- sapply(X = fov, FUN = function(x) {
                  return(DefaultBoundary(object = object[[x]]))
              }, simplify = FALSE, USE.NAMES = TRUE)

              boundaries <- Seurat:::.BoundariesByImage(object = object, fov = fov, boundaries = boundaries)
              fov <- names(x = boundaries)

              overlap <- rep_len(x = overlap, length.out = length(x = fov))
              crop <- rep_len(x = crop, length.out = length(x = fov))
              names(x = crop) <- fov
              group.by <- boundaries %!NA% group.by %||% "ident"
              vars <- c(group.by, split.by)
              md <- FetchData(object = object, vars = vars[!is.na(x = vars)], cells = cells)
            
              pnames <- unlist(x = lapply(X = seq_along(along.with = fov),
                  FUN = function(i) {
                      return(if (isTRUE(x = overlap[i])) {
                          fov[i]
                      } else {
                          paste(fov[i], boundaries[[i]], sep = "_")
                      })
                  }))

              pdata <- vector(mode = "list", length = length(x = pnames))
              names(x = pdata) <- pnames

              for (i in names(x = pdata)) {
                i<-names(x = pdata)
                  ul <- unlist(x = strsplit(x = i, split = "_"))
                  img <- paste(ul[1:length(ul) - 1], collapse = "_")
                  lyr <- ul[length(ul)]
                  if (is.na(x = lyr)) {
                      lyr <- boundaries[[img]]
                  }

                  pdata[[i]] <- lapply(X = lyr, FUN = function(l) {
                      if (l == "NA") {
                          return(NA)
                      }
                      df <- fortify(model = object[[img]][[l]])
                      df <- df[df$cell %in% cells, , drop = FALSE]
                      if (!is.null(x = md)) {
                          df <- merge(x = df, y = md, by.x = "cell", by.y = 0,
                            all.x = TRUE)
                      }
                      df$cell <- paste(l, df$cell, sep = "_")
                      df$boundary <- l
                      return(df)
                  })
                  
                  pdata[[i]] <- do.call(what = "rbind", args = pdata[[i]])
              
              }

                  mdata <- NULL
              
              plots <- vector(mode = "list", length = length(x = pdata) *
                  ifelse(test = length(x = group.by), yes = length(x = group.by),
                      no = 1L))
              idx <- 1L
              # for (group in group.by) {
                group <- group.by
                  # for (i in seq_along(along.with = pdata)) {
                    i<-1
                      img <- unlist(x = strsplit(x = names(x = pdata)[i],
                          split = "_"))[1L]

              order_colon<-names(pop_cl_size)
              order_cells<-gsub("centroids_", "", pdata[[i]] %>% arrange(factor(object@meta.data[,group.by], levels = rev(order_colon))) %>% pull(cell))
              sizes2<-pop_cl_size_temp[object@meta.data[,group.by]]
              names(sizes2) <-  rownames(object@meta.data)

          # sizes2[order_cells]

            # pdf("aver.pdf", width=10, height=15)
             image_tissue <- SingleImagePlot(data = pdata[[i]] %>% arrange(factor(object@meta.data[,group.by], levels = rev(order_colon))), col.by = pdata[[i]] %!NA%
                group, molecules = mdata[[img]], cols = cols,
                shuffle.cols = FALSE, size =sizes2[order_cells], alpha = alpha,
                mols.size = 0.1, mols.cols = NULL,
                mols.alpha = 1, border.color = "white",
                border.size = NA, na.value = "grey50",
                dark.background = TRUE) + coord_fixed() + NoAxes(panel.background = element_blank())
              return(image_tissue)

  }

  setwd("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset")

  outdir_spa <- "/mnt/bioadhoc/Groups/vd-vijay/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/spatial_redo"
  dir.create("spatial_redo")
  dir.create("spatial_redo/distance")

  slide1 <- readRDS("/mnt/BioAdHoc/Groups/vd-vijay/donaldo/spatial/SaHe/slide_5k_SaHe1/re_segmentation/info/seurat_objects/slide1_all_cells_50pct_30pc_with_populations_and_tcells_DCs_defined_Dec_15_2025.rds")
  # slide1$orig.donor <- slide1$orig.donor_paper
  table(slide1$orig.donor)
  slide1$cellnames <- rownames(slide1@meta.data)
  slide1 ### delate macrophages
  slide2 <- readRDS("/mnt/BioAdHoc/Groups/vd-vijay/donaldo/spatial/SaHe/slide_5k_SaHe2/re_segmentation/info/seurat_objects/slide2_all_cells_50pct_30pc_with_populations_and_tcells_fib_DC_defined_dec_15_2025.rds")
  slide2$cellnames <- rownames(slide2@meta.data)
  table(slide2$orig.donor)

  keep_cells2 <- slide2$cellnames[! slide2$orig.donor %in% c("NIHHC_326", "NIHMA_322")]
  keep_cells1 <- slide1$cellnames[! slide1$orig.donor %in% c("NIHW_319")]

  slide1 <- subset(slide1, cells = as.character(keep_cells1))
  slide2 <- subset(slide2, cells = as.character(keep_cells2))
  # disease_slide1 <- c("HC"= "Healthy", "MA"= "Mild.asthma",  "SA"= "Severe.asthma")
  # slide1$orig.disease <- as.character(disease_slide1[slide1$orig.disease])

  seuobj <- list(
    "slide1"= slide1, 
    "slide2"= slide2)

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

  ## 1st step calculate distances from T cells to all other celltypes
  source("/home/fcastaneda/bin/clustering/R/stats_summary_table.R")

  ## Calculating the distance matrix per each seurat object
  all_data <- lapply(names(seuobj), function(slide_test){
    # slide_test <- names(seuobj)[1]
    slide_xe <- seuobj[[slide_test]]
    slide_xe@meta.data$orig.celltype_cellchat <- slide_xe@meta.data$orig.celltype_tc_dc

    coords <- Seurat::GetTissueCoordinates(slide_xe)
    colnames(coords)[3] <- "cellname"

    slide_xe@meta.data$cellname <- rownames(slide_xe@meta.data)

    all(coords$cellname %in% slide_xe@meta.data$cellname)

    coords2 <- left_join(coords, slide_xe@meta.data[, c("cellname", "orig.celltype_cellchat", "orig.donor")], by = 'cellname')
    
    colnames(coords2) <- c("Cell X Position", "Cell Y Position", "Cell ID", "Phenotype", "sample")
    coords2 <- coords2 %>% select("Cell ID", "Cell X Position", "Cell Y Position", Phenotype, sample) 
    rownames(coords2) <- coords2$'Cell ID'

    head(coords2)

     list_dist_mtx <- lapply(unique(coords2$sample), function(donor) {
        # donor <- unique(coords2$sample)[1]
        donor_distances <- coords2 %>% filter(sample == donor) 
        dist_mtx <- distance_matrix(donor_distances)
        return(dist_mtx)
    })

    names(list_dist_mtx) <- unique(coords2$sample)

    celltypes_comp_orig <- unique(slide_xe$orig.celltype_cellchat)
    # celltypes_comp_orig <- c("fibroblasts", "CD8", "CD4", "MAIT")

    dir_spatial <- paste0("spatial_redo/distance/", slide_test, "/")
    dir.create(dir_spatial)

    all_list <- lapply(celltypes_comp_orig, function(from_celltype){ 

        # from_celltype <-  "fibroblasts"
        cat(from_celltype, "\n")
          
        matrix_donor <- lapply(names(list_dist_mtx), function(donor){   
          # donor <- names(list_dist_mtx)[1]
          # donor <- "NIHW_577"

          # 
          # cat(donor, "\n")
          sample_matrix <- list_dist_mtx[[donor]]
          celltypes_comp <- as.character(celltypes_comp_orig)
          celltypes_comp <- grep(from_celltype, celltypes_comp, value=TRUE, invert=TRUE)
          # celltypes_comp <- gsub("^c[0-9]{1,2}_", "", celltypes_comp)
          # celltypes_comp <- grep(celltypes_comp, )

          distance_to_closest_cell <- lapply(celltypes_comp, function(to_celltype){ 
            # to_celltype <- "DCs"
            # to_celltype <- "c12_schwann_cells"

            # cat(to_celltype, "\n")
          
            from_celltype_cells <- slide_xe@meta.data %>% filter(orig.celltype_cellchat == from_celltype  & orig.donor ==  donor) %>% pull(cellname)
            to_celltype_cells <- slide_xe@meta.data %>% filter(orig.celltype_cellchat == to_celltype & orig.donor ==  donor) %>% pull(cellname)

            if(length(from_celltype_cells) > 1 & length(to_celltype_cells) > 1) {

              sample_matrix2 <- sample_matrix[to_celltype_cells, from_celltype_cells]
              to_ret <- data.frame(
                distance = rowMins(sample_matrix2),
                to = to_celltype, 
                cell= names(rowMins(sample_matrix2)),
                donor= donor)

            } else {
              to_ret <- data.frame(
                distance = 0,
                to = to_celltype, 
                cell= "no_cells",
                donor= donor)
            }
            return(to_ret) 
          })

          distance_from_tcells <- Reduce(rbind, distance_to_closest_cell)
          return(distance_from_tcells)
        })

        donor_matrix <-  Reduce(rbind,matrix_donor)
        
        dist_out <- paste0(dir_spatial, from_celltype)
        # dist_out <- paste0("distance_closest/", from_celltype  ,"/")

        # dir.create(dist_out)
        # cat("Plotting .... ")
        
        # pdf(paste0(dist_out, "/Per_donor.pdf"))

        #   distance_plot <- ggplot(donor_matrix, aes(x = donor, y = distance)) + 
        #     geom_jitter(position = position_jitter(0.2), alpha=0.7) + 
        #     geom_violin(width = 0.5) +
        #     geom_boxplot(width = 0.1) + facet_wrap_paginate(~ to, ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance to ", from_celltype)) 

        #   for(i in 1:n_pages(distance_plot)){
        #       p_save <-  distance_plot + 
        #         facet_wrap_paginate(~ to, ncol = 1, nrow = 1, scales = "free", page = i)
        #       print(p_save)
        #     }

        # dev.off()

        #   pdf(paste0(dist_out, "/Per_celltype.pdf"))

        #   distance_plot <- ggplot(donor_matrix, aes(x = to, y = distance)) + 
        #     geom_jitter(position = position_jitter(0.2), alpha=0.7) + 
        #     geom_violin(width = 0.5) +
        #     geom_boxplot(width = 0.1) + facet_wrap_paginate(~donor, ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance to ", from_celltype))

        #   for(i in 1:n_pages(distance_plot)){
        #       p_save <-  distance_plot + 
        #         facet_wrap_paginate(~ donor, ncol = 1, nrow = 1, scales = "free", page = i)
        #       print(p_save)
        #     }

        # dev.off()

        # pdf(paste0(dist_out, "/Per_disease.pdf"))
        #   donor_matrix$disease <- str_extract(donor_matrix$donor, "HC|MA|W")
        #   # donor_matrix
        #   distance_plot <- ggplot(donor_matrix, aes(x = disease, y = distance)) + 
        #     geom_jitter(position = position_jitter(0.2), alpha=0.7) + 
        #     geom_violin(width = 0.5) +
        #     geom_boxplot(width = 0.1) + facet_wrap_paginate(~ to, ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance to ", from_celltype))

        #   for(i in 1:n_pages(distance_plot)){
        #       p_save <-  distance_plot + 
        #         facet_wrap_paginate(~ to, ncol = 1, nrow = 1, scales = "free", page = i)
        #       print(p_save)
        #     }

        # dev.off()

        # pdf(paste0(dist_out, "/Per_disease_all.pdf"), width=25)
        #   donor_matrix$disease <- str_extract(donor_matrix$donor, "HC|MA|W")
        #   # donor_matrix
        #   distance_plot <- ggplot(donor_matrix, aes(x = disease, y = distance)) + 
        #     geom_jitter(position = position_jitter(0.2), alpha=0.7) + 
        #     geom_violin(width = 0.5) +
        #     geom_boxplot(width = 0.1) + facet_wrap_paginate(~ to, ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance to ", from_celltype))+ facet_wrap(~ to, nrow=1)

        #   print(distance_plot)
        # dev.off()

        # pdf(paste0(dist_out, "/Per_disease_celltype.pdf"))
        #   donor_matrix$disease <- str_extract(donor_matrix$donor, "HC|MA|W")
        #   # donor_matrix
        #   distance_plot <- ggplot(donor_matrix, aes(x = to, y = distance)) + 
        #     geom_jitter(position = position_jitter(0.2), alpha=0.7) + 
        #     geom_violin(width = 0.5) +
        #     geom_boxplot(width = 0.1) + facet_wrap_paginate(~ disease , ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance to ", from_celltype))

        #   for(i in 1:n_pages(distance_plot)){
        #       p_save <-  distance_plot + 
        #         facet_wrap_paginate(~ disease, ncol = 1, nrow = 1, scales = "free", page = i)
        #       print(p_save)
        #     }

        # dev.off()
      donor_matrix$disease <- str_extract(donor_matrix$donor, "HC|MA|W")
      write.csv(donor_matrix, paste0(dist_out, "_distance_matrix.csv"))
      return(donor_matrix)

    })
  
  })

  distance_files <- list.files("/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/spatial_redo/distance/", pattern="distance_matrix.csv", recursive=TRUE, full.names=TRUE)

  celltypes_comp_orig <-  c("CD8", "CD4", "MAST")


############################################################  
########### Distance   ################
############################################################  

## Code to calculate the distance to every population from a given celltype of interest 

  lapply(celltypes_comp_orig, function(per_celltype){ 
      per_celltype <- celltypes_comp_orig[1]
      dir.create(paste0(outdir_spa, "/analysis/"))
      dir.create(paste0(outdir_spa, "/analysis/", per_celltype))

      files_celltype <- grep(per_celltype, distance_files, value=TRUE)
      
      distance_to_celltypes <-Reduce(rbind, lapply(files_celltype, function(distances_files){
        distances_slides <- read.csv(distances_files)
        distances_slides$slide <- str_extract(distances_files, "slide1|slide2")
        return(distances_slides)
        }))

          # distance_to_celltypes
            distance_plot <- ggplot(distance_to_celltypes, aes(x = to, y = distance, fill=to)) + 
            geom_jitter(position = position_jitter(0.2), alpha=0.1, size=0.1) + 
            geom_violin(width = 1) +
            scale_fill_manual(values=pop_cl_colors) +
            geom_boxplot(width = 0.1, outlier.shape = NA) + facet_wrap_paginate(~ disease , ncol = 1, nrow = 1, scales = "free") + theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle(paste0("Distance from ", per_celltype)) + scale_y_continuous(trans = log2_trans()) 
            
      distance_to_celltypes$to_merge <- ifelse(distance_to_celltypes$to %in% c("secretory1", "secretory2"), "secretory", distance_to_celltypes$to)     
      write.csv(distance_to_celltypes, paste0(dist_sum, "/Per_disease_", per_celltype,".csv"))

      dist_sum <- "spatial_redo/distance/summary"
      dir.create(dist_sum)

      ## Distance per celltype summary
      pdf(paste0(dist_sum, "/Per_disease_", per_celltype,".pdf"), width=20, height=5)
    
        for(i in 1:n_pages(distance_plot)){
            p_save <-  distance_plot +  geom_hline(yintercept=50) +
               facet_wrap_paginate(~ disease, ncol = 1, nrow = 1, scales = "free", page = i)
            print(p_save)
          }

      dev.off()

  })

#
#
######################################################################################  
####### Concentric bins distance vs gene expression analysis.  #######################
######################################################################################

## Code to calculate the correlation between positive cells and distance from a reference celltype

  # setwd("donut_corr") 
  dir.create("spatial_redo/donut_corr/")

  seuobj <- list(
    "slide1"= slide1, 
    "slide2"= slide2)

  ## Analysis
  all_data <- lapply(names(seuobj), function(slide_test){
      # slide_test <- "slide1"
      slide_xe <- seuobj[[slide_test]]
      slide_xe@meta.data$orig.celltype_cellchat <- slide_xe@meta.data$orig.celltype_tc_dc

      coords <- Seurat::GetTissueCoordinates(slide_xe)
      colnames(coords)[3] <- "cellname"

      slide_xe@meta.data$cellname <- rownames(slide_xe@meta.data)

      all(coords$cellname %in% slide_xe@meta.data$cellname)

      coords2 <- left_join(coords, slide_xe@meta.data[, c("cellname", "orig.celltype_cellchat", "orig.donor")], by = 'cellname')
      
      colnames(coords2) <- c("Cell X Position", "Cell Y Position", "Cell ID", "Phenotype", "sample")
      coords2 <- coords2 %>% select("Cell ID", "Cell X Position", "Cell Y Position", Phenotype, sample) 
      rownames(coords2) <- coords2$'Cell ID'

      head(coords2)

      list_dist_mtx <- lapply(unique(coords2$sample), function(donor) {
          # donor <- unique(coords2$sample)[1]
          donor_distances <- coords2 %>% filter(sample == donor) 
          dist_mtx <- distance_matrix(donor_distances)
          return(dist_mtx)
      })

      names(list_dist_mtx) <- unique(coords2$sample)
  
      celltypes_comp <- as.character(unique(slide_xe$orig.celltype_cellchat))
      # celltypes_comp <- c("CD8", "MAST") ## Define them if your are interested just in some of them
      # celltypes_comp <- c("CD8") #, "MAST")
      # celltypes_comp <- c("fibroblasts", "smooth_muscle", "endothelial")

      lapply(celltypes_comp, function(from_celltype){ 
        # from_celltype <- celltypes_comp[1]
        cat(from_celltype, "\n")
        out_dir1<- paste0("spatial_redo/donut_corr/", slide_test, "/", from_celltype, "/")
        dir.create(paste0("spatial_redo/donut_corr/", slide_test, "/"))
        dir.create(out_dir1)

           ## Selecting gene lists
        files_name <- paste0("/home/fcastaneda/vd-vijay/donaldo/spatial/SaHe/slide_5k_SaHe2/re_segmentation/results/figures/crater/disease/merged_nocov_legend_3_dots/")
        to_dgea <- list.files(files_name, pattern="cpm_results", full.names = TRUE, recursive=TRUE)

        celltypes_to <- grep(from_celltype, unique(slide_xe$orig.celltype_cellchat), value=TRUE, invert=TRUE)
        celltype_dgea <- str_extract(to_dgea, paste0(celltypes_to, collapse="|") )
        celltypes_to<- celltypes_to[celltypes_to %in% celltype_dgea]

        # celltypes_to <- c("fibroblasts", "smooth_muscle", "endothelial")
        celltypes_to <- c("CD8")

        to_celltype<-celltypes_to

        matrix_donor <- lapply(names(list_dist_mtx), function(donor){   
            # donor <- names(list_dist_mtx)[1]
            
            cat(donor, "\n")
            out_dir1<- paste0(out_dir1, "/", donor, "/")
            dir.create(out_dir1)
            sample_matrix <- list_dist_mtx[[donor]] 

            slide_xe_temp <- subset(slide_xe, orig.donor == donor)

            from_celltype_cells <- slide_xe_temp@meta.data %>% filter(orig.celltype_cellchat == from_celltype  & orig.donor ==  donor) %>% pull(cellname)
            sample_matrix2 <- sample_matrix[, from_celltype_cells]

          ## donut approach
          
            out_dir1_donut <-paste0(out_dir1, "donut_approach/")
            dir.create(out_dir1_donut)

            sample_matrix2_temp <- sample_matrix2
            # sample_matrix2_temp[sample_matrix2_temp>200] <- NA
            try(rm(keep_cells_dist))
            keep_cells_dist<-c("no_cell")

            ### Correlation
            sample_matrix2_temp <- sample_matrix2
            # sample_matrix2_temp[sample_matrix2_temp>200] <- NA
            try(rm(keep_cells_dist))
            keep_cells_dist<-c("no_cell")

            bines_distances_corr <- lapply(c(15,35,55,75,95,115,135,155,175,195), function(dist){
                # dist<-15
                # cat(dist, "\n")
                sample_matrix2_r<-rowSums(sample_matrix2_temp<dist,na.rm = TRUE)
                # sample_matrix2_r <- rowSums(sample_matrix2<distances_ranges)
                keep_cells <- names(sample_matrix2_r[sample_matrix2_r>0])
                # keep_cells<- grep(paste(keep_cells_dist, collapse="|"), keep_cells, invert=TRUE, value=TRUE )
                keep_cells<- keep_cells[!keep_cells %in% keep_cells_dist]
                slide_xe_temp$keep_cells <- ifelse(slide_xe_temp$cellname %in% keep_cells, slide_xe_temp$orig.celltype_cellchat, "no_color")

                slide_xe_temp$keep_cells <-ifelse(slide_xe_temp$orig.celltype_cellchat == from_celltype, from_celltype, slide_xe_temp$keep_cells)

                distance_to_closest_cell_bin <- lapply(celltypes_to, function(to_celltype2){ 
                      # to_celltype2 <- celltypes_to[1]
                      cat(to_celltype2, "\n")

                      to_dgea_sub <- grep(to_celltype2, to_dgea, value=TRUE)

                      dgea_c <- lapply(to_dgea_sub, function(x){
                        # x <- to_dgea[1]
                        genes_select <- read.csv(x) %>% select(excel_name, log2FoldChange, padj, group)
                        genes_select$gene_name <- gsub("'", "", genes_select$excel_name)
                        genes_select$comp <- gsub(files_name, "", x)
                        genes_select$comp <- gsub("_.*|/", "", genes_select$comp)
                        return(genes_select)
                      })

                      dgea_c <- Reduce(rbind, dgea_c)

                      table(dgea_c$group, dgea_c$comp)
                      # table(dgea_c$log2FoldChange, dgea_c$group)
                      table(dgea_c$comp, useNA="always")

                      genes_severe_specific <- dgea_c %>% filter(group == "Severe.asthma") %>% group_by(gene_name) %>% summarize(rep=n()) %>% filter(rep >1 ) %>% pull(gene_name)

                      if(to_celltype2 == "CD8"){ 
                        genes_cd8 <- c("AREG", "ENTPD1", "HLA-DRA", "HLA-DQA1", "FABP5", "GLUL", "BATF", "RBPJ", "RUNX1", "RUNX2", "GZMB", "GZMH")
                        
                        genes_severe_specific <- genes_cd8[genes_cd8 %in% rownames(slide_xe_temp)]
                        } else { 
                          genes_severe_specific <- genes_severe_specific
                        }


                      genes_list <- list("1st"=unique(as.character(genes_severe_specific)))
                      names(genes_list) <- to_celltype2
  
                      to_celltype_cells <- slide_xe_temp@meta.data %>% filter(orig.celltype_cellchat == to_celltype2  & orig.donor ==  donor & keep_cells != "no_color") %>% pull(cellname)

                      genes_ins <- genes_list[[to_celltype2]]
                      tmp2<-add_gene_tag(genes_ins, slide_xe_temp@meta.data, as.matrix(slide_xe_temp@assays$Xenium$data[genes_ins, ]), thresh=0.01, tag = c('tag', 'p', 'n'))
                      slide_xe_temp@meta.data = joindf(slide_xe_temp@meta.data, as.data.frame(tmp2))

                      slide_xe_temp$cellname <- rownames(slide_xe_temp@meta.data)

                      genes_to_focus <- c("cellname", paste0("tag_", genes_ins))

                      # to_celltype_cells <- slide_xe_temp@meta.data %>% filter( == to_celltype  & orig.donor ==  donor) %>% pull(cellname)
                    
                      if(length(to_celltype_cells) > 1) { # Greater than 1 or 0 ???
                        
                        sample_matrix2_1 <- sample_matrix2_temp[to_celltype_cells, ]
                        sample_matrix2_1_temp <- sample_matrix2_1<dist & sample_matrix2_1>0
                        ncells_surrounding <- rowSums(sample_matrix2_1_temp,na.rm = TRUE)

                        sample_matrix2_1[sample_matrix2_1>dist]<- NA

                        closest_cell <- rowMins(sample_matrix2_1, na.rm = TRUE)
                          # sample_matrix2_1[sample_matrix2_1<1]<- NA Line documented because in this vector the same cells shouldn't be included (same cells means the diagonal is not 0)
                        mean_distance <- rowMeans(sample_matrix2_1, na.rm = TRUE)
                        if(!all(names(closest_cell) == names(ncells_surrounding))) stop("Wrong assigantion of cell IDs")
                          
                            df_distances <- data.frame(
                              row.names= names(ncells_surrounding),
                              surrounding=ncells_surrounding,
                              min_dist = closest_cell[names(ncells_surrounding)],
                              mean_dist = mean_distance[names(ncells_surrounding)]
                            )

                        keep_meta <- slide_xe_temp@meta.data[to_celltype_cells,genes_to_focus]
                        class_cells <- apply(keep_meta[,genes_to_focus], 2, function(patt) grepl("p$", patt))
                        rownames(class_cells) <- rownames(keep_meta)
                        pos_cells2 <- rowSums(class_cells)[rowSums(class_cells) > 0]
                          # neg_cells <- rowSums(class_cells)[rowSums(class_cells) < 0]
                        keep_meta$signature <- ifelse(rownames(keep_meta) %in% names(pos_cells2), "sig_pos", "sig_neg")

                        keep_meta2 <- keep_meta %>% select(signature, all_of(genes_to_focus), cellname)

                        # keep_meta2
                        # keep_meta$chosen_gene <- keep_meta[,genes_to_focus]

                        radius_plot_gene <- merge(df_distances, keep_meta2, by="row.names", all=TRUE)
                        if(!all(nrow(radius_plot_gene) == nrow(df_distances))) stop("Wrong assigantion of cell IDs")
                        radius_plot_gene$donor<-donor
                        radius_plot_gene$to_celltype<-to_celltype2
                        radius_plot_gene$distance_bin <- dist
                        radius_plot_gene$Row.names <- NULL

                      } else{

                        df_distances <- data.frame(
                            row.names= "no_cells",
                            surrounding=NA,
                            min_dist = NA,
                            mean_dist = NA
                          )

                        for(i in c("signature", genes_to_focus)){
                        df_distances[,i] <- NA   
                        }
                        
                        df_distances$donor<-donor
                        df_distances$to_celltype<-to_celltype2
                        df_distances$distance_bin <- dist

                        radius_plot_gene<-df_distances

                      }
                      return(radius_plot_gene)
                })

                names(distance_to_closest_cell_bin) <- celltypes_to

                colnames(sample_matrix2_temp)
                keep_cells_dist<<-c(keep_cells_dist,keep_cells)

                return(distance_to_closest_cell_bin)
            })

            names(bines_distances_corr) <- paste0("d", c(15,35,55,75,95,115,135,155,175,195))

            donut_corr_list <- lapply(celltypes_to, function(cellty){ 
                # cellty<-celltypes_to[1]
                bines_distances_celty <- lapply(bines_distances_corr, function(x) x[[cellty]]) 
                return(Reduce(rbind, bines_distances_celty)) 
            })
            names(donut_corr_list) <- celltypes_to
            
            saveRDS(donut_corr_list, paste0(out_dir1_donut, "distances_corr.rds"))
            # saveRDS(radius_list, paste0(out_dir1_radius, "distances_complete.rds"))

        })

      })

  })

  ### Summary

    source("/home/vfajardo/bioadhoc/HIPC/paper_developments/HIPC-Lung/paper_items/jobs_scripts/functions_collection.0.2.R")

        blank.complement.2e <- theme(
                  line=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none',
                  axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.6), axis.ticks.length=unit(0.2, "cm")
                )
        # files_disease_donut <- list.files(out_dir1, pattern="bines", recursive=TRUE, full.names = TRUE)
        out_dir0 <- "spatial_redo/donut_corr/"
        dir.create(paste0(out_dir0, "/summary_donut/"))

      # celltypes_from_sum <- c("fibroblasts", "smooth_muscle", "endothelial")
      # celltypes_from_sum <- c("CD8")

      ### Summary donuts correlations
      lapply(celltypes_from_sum, function(from_sum){ 
        # from_sum <- celltypes_from_sum[1]
        # out_dir1<-  "spatial_slide2_new_seg/new_analysis4_check/CD8//check/"
        files_disease_donut2 <- list.files(out_dir0, pattern="corr.rds", recursive=TRUE, full.names = TRUE)
        files_disease_donut2 <- grep(from_sum, files_disease_donut2, value=TRUE)
        dir.create(paste0(out_dir0, "/summary_donut/"))
        out_dir1<-paste0(out_dir0, "/summary_donut/", from_sum)
        dir.create(out_dir1)

        # celltypes_to <- c("fib1")
        # celltypes_to <- celltypes_to[-c(2, 7)]

        # celltypes_to <- c("fibroblasts", "smooth_muscle")

        # celltypes_to <- c("CD8")
        # celltypes_to <- c("fibroblasts", "smooth_muscle", "endothelial")
        
        ## Modify accdordingly
        corre_plot <- lapply(c("NIHW", "NIHHC", "NIHMA"), function(disease){
          
          # disease <- "NIHW"
          dir.create(paste0(out_dir1, "/", disease))
          cat("Summary for", disease, "\n")

          files_disease_donut_keep <- grep(disease, files_disease_donut2, value=TRUE)

          distance_data <- lapply(files_disease_donut_keep, readRDS)

          names(distance_data) <- gsub(paste0(out_dir1, "/|/donut_approach/distances_bines.rds"), "" ,files_disease_donut_keep)

          aver_summary <- lapply(celltypes_to, function(cellty){ 
              # cellty<-celltypes_to[2]
              # cellty <- "fibroblasts"
              cat("### Summary for", cellty, "\n")
              dir.create(paste0(out_dir1, "/", disease, "/", cellty))

              distance_data_celty <- lapply(distance_data, function(x) x[[cellty]]) 

              data_celltype<- Reduce(rbind, distance_data_celty)
              data_celltype$surrounding_logical <- ifelse(data_celltype$surrounding == 0, "no_cells", "cells")
              genes_plot<-grep("tag", colnames(data_celltype), value=TRUE)

              keep_gene <- unlist( lapply(data_celltype[,genes_plot], function(x){
                # x<-data_celltype[,genes_plot[1]]
                x<-paste0(x, "_", data_celltype$distance_bin)
                genes_values<-table(x)[grepl("p_", names(table(x)))]
                return(length(names(genes_values))>3)
              }) )

              genes_plot <- genes_plot[keep_gene]
              # genes_plot<-c("tag_TCF21", "tag_MAOA", "tag_CCN2")

              plot_all_corr <- lapply(c("signature", genes_plot), function(gene1){ 
                gene1 <- genes_plot[2]
                # gene1 <- "signature"
                cat(gene1, "\n")
                data_celltype <-  data_celltype %>% filter(!is.na(surrounding))

                # data_celltype %>% filter(is.na(surrounding))
                
                data_celltype$surrounding_logical <- ifelse(data_celltype$surrounding_logical == "cells", paste0(data_celltype$surrounding_logical, "_", data_celltype[,gene1]), data_celltype$surrounding_logical) #this is redundant

                data_celltype$surrounding_logical <- gsub("^cells_", "", data_celltype$surrounding_logical) #this is redundant

                table( data_celltype$surrounding_logical, data_celltype[,gene1])

                pos_gene <- grep("p$|pos$",unique(data_celltype[,gene1]), value=TRUE)
                neg_gene <- grep("n$|neg$",unique(data_celltype[,gene1]), value=TRUE)
                
                percentage_celltype_with_cells <- data_celltype %>% group_by(distance_bin, surrounding_logical) %>% summarize(cells=n()) 

                # bines_plot <- c("15","35","55","75","95","115","135","155","175","195", "215")

                bines_plot <- c("15","35","55","75","95","115","135","155")

                bines_plot <- bines_plot[bines_plot %in% unique(percentage_celltype_with_cells$distance_bin)]

                percentage_celltype_with_cells <- percentage_celltype_with_cells %>% filter(distance_bin %in% bines_plot)

                percentage_celltype_with_cells$distance_bin <- factor(as.character(percentage_celltype_with_cells$distance_bin), bines_plot)

                percentage_celltype_with_cells$surrounding_logical <- factor(as.character(percentage_celltype_with_cells$surrounding_logical), levels= c( pos_gene, neg_gene))

                number_cells <- percentage_celltype_with_cells %>% pull(cells) %>% sum()
                gene_name <- gsub("p$|_pos$", "", pos_gene)

                percentage_celltype_with_cells_corr <- percentage_celltype_with_cells %>% group_by(distance_bin) %>% mutate(prop=cells/sum(cells)*100)

                percentage_celltype_with_cells_corr$distance_bin2 <- as.numeric(as.character(percentage_celltype_with_cells_corr$distance_bin))

                table(percentage_celltype_with_cells_corr$distance_bin2, percentage_celltype_with_cells_corr$distance_bin)

                # table(percentage_celltype_with_cells_corr$distance_bin2, percentage_celltype_with_cells_corr$distance_bin)
                # write.csv(percentage_celltype_with_cells_corr, "")
                
                  percentage_celltype_with_cells_corr_pos <- percentage_celltype_with_cells_corr %>% filter(surrounding_logical == pos_gene)

                    corr.res <- cor.test(
                            x=percentage_celltype_with_cells_corr_pos$distance_bin2,
                            y=percentage_celltype_with_cells_corr_pos$prop, method='pearson')

                    correla_gene <- data.frame(
                      cor = as.character(corr.res$estimate),
                      pvalue = as.character(corr.res$p.value),
                      tot_cells = number_cells,
                      gene = gene1,
                      celltype = cellty
                    )

                    # plot_corr <-  
                choseen_color <- case_when(
                  gene1 == "tag_PDGFRA" ~ "#E432F7",
                  gene1 == "tag_ATP1A1" ~ "#E432F7",
                  gene1 == "tag_PECAM1" ~ "#E432F7",
                  gene1 == "tag_MAOA" ~ "#E432F7",
                  gene1 == "tag_SERPINA3" ~ "#F0A04F",
                  gene1 == "tag_SNAI2" ~ "#F0A04F",
                  gene1 == "tag_FZD4" ~ "#F0A04F",
                  gene1 == "tag_TCF21" ~ "#67DCD0",
                  gene1 == "tag_MUC5AC" ~ "#67DCD0",
                  gene1 == "tag_CCL14" ~ "#67DCD0",
                  gene1 == "tag_CCN2" ~ "#9DF15C",
                  gene1 == "tag_PEACAM1" ~ "#E432F7",
                  gene1 == "tag_TGFBR2" ~ "#F0A04F",
                  gene1 == "tag_EPAS1" ~ "#9DF15C",
                  TRUE ~ "black"
                )

                breaks_plot <- unique(c(floor(min(percentage_celltype_with_cells_corr_pos$prop))-3, ceiling(max(percentage_celltype_with_cells_corr_pos$prop))+3))

                breaks_add <- if ((breaks_plot[2]-breaks_plot[1])<10) seq(floor(min(percentage_celltype_with_cells_corr_pos$prop)), ceiling(max(percentage_celltype_with_cells_corr_pos$prop)), 3) else seq(floor(min(percentage_celltype_with_cells_corr_pos$prop)), ceiling(max(percentage_celltype_with_cells_corr_pos$prop)), 5)

                breaks_plot <- unique(c(breaks_plot, breaks_add))

                breaks_plot <- unique(c(breaks_add))

                breaks_plot <- breaks_plot[order(breaks_plot)]

                breaks_plot<- c(5, 15, 25)
                # lbreaks <- length(breaks_plot)
                # if(breaks_plot[lbreaks]-breaks_plot[lbreaks-1]==1) breaks_plot[lbreaks-1] <- NA

                # breaks_plot<-breaks_plot[!is.na(breaks_plot)]

                plot_corr <- ggplot(data=percentage_celltype_with_cells_corr_pos, aes(x=distance_bin2, y=prop)) +
                      # scale_colour_manual
                      ggpubr::stat_cor(method='pearson', label.x = 50, label.y = ceiling(max(percentage_celltype_with_cells_corr_pos$prop))+0.5, size=2) +
                      geom_smooth(method="lm", se = TRUE, level = 0.95, color="red", linetype = "dashed", fill="#bfbdbd") +
                      geom_point(size=3, shape=16, colour=choseen_color) +
                      # scale_y_continuous(breaks = breaks_plot,limits = c(floor(min(percentage_celltype_with_cells_corr_pos$prop)-3), ceiling(max(percentage_celltype_with_cells_corr_pos$prop)+3))) +
                      scale_y_continuous(breaks = breaks_plot,limits = c(5,25)) +
                      scale_x_continuous(breaks = unique(percentage_celltype_with_cells_corr_pos$distance_bin2), limits=c(0,ceiling(max(percentage_celltype_with_cells_corr_pos$distance_bin2)))) +
                      labs(
                        title = cellty,
                        x="Distance: Donuts bins",
                        y=pos_gene) + theme(legend.position = "none", axis.text.x = element_text(color = "black", size = 4), axis.text.y = element_text(color = "black", size = 6)) + blank.complement.2e 
                    pdf(paste0(out_dir1, "/", disease, "/", cellty, "/", gene1, "v3.pdf"), width=2, height=2)
                    print(plot_corr)
                    dev.off()
                    # # print(p)
                    # dev.off()

                    percentage_celltype_with_cells_corr_pos$gene <- pos_gene


                    cat("Enter \n")
                return(list(plots=plot_corr, corr_data=percentage_celltype_with_cells_corr_pos))
                
              })

              names(plot_all_corr) <- c("signature", genes_plot)
              
              plot_all_corr_plots <- lapply(plot_all_corr, function(x) x[["plots"]])
              corr_data_plot <- lapply(plot_all_corr, function(x) x[["corr_data"]])

              corr_data_plot_2 <- Reduce(rbind, corr_data_plot)
              corr_data_plot_2$disease <- disease
              # corr_data_plot_2$gene <- gsub("tag_", "", corr_data_plot_2$gene)
              # corr_data_plot_2$significant <- ifelse(corr_data_plot_2$pvalue < 0.05,corr_data_plot_2$cor, NA)

                pdf(paste0(out_dir1, "/", disease, "/", cellty, "/genes_correlations_155_color_SE.pdf"), width=6, height=length(plot_all_corr_plots)/1.5)
                print(ggarrange(plotlist=plot_all_corr_plots, ncol=3, nrow=ceiling((length(genes_plot)+1)/ 3)))
                dev.off()

                write.csv(corr_data_plot_2, paste0(out_dir1, "/", disease, "/", cellty, "/genes_correlations_155.csv"))
              return(corr_data_plot_2)

          })
          
          return(Reduce(rbind,aver_summary))

        })

        correp <- Reduce(rbind,corre_plot)
      })
#
#
############################################################
################# Cells surrounding  #######################
############################################################

## Code to count the cells surrounding (50 nm) a given populationof interest 

      seuobj <- list(
        "slide1"= slide1, 
        "slide2"= slide2)
      
      dir.create("spatial_redo/donut_surrounding/")
      
      ## Analysis
      all_data <- lapply(names(seuobj), function(slide_test){
          slide_xe <- seuobj[[slide_test]]
          slide_xe@meta.data$orig.celltype_cellchat <- slide_xe@meta.data$orig.celltype_tc_dc

          slide_xe@meta.data$orig.celltype_cellchat <- ifelse(slide_xe@meta.data$orig.celltype_cellchat %in% c("secretory1", "secretory2"), "secretory", slide_xe@meta.data$orig.celltype_cellchat)

          coords <- Seurat::GetTissueCoordinates(slide_xe)
          colnames(coords)[3] <- "cellname"

          slide_xe@meta.data$cellname <- rownames(slide_xe@meta.data)

          all(coords$cellname %in% slide_xe@meta.data$cellname)

          coords2 <- left_join(coords, slide_xe@meta.data[, c("cellname", "orig.celltype_cellchat", "orig.donor")], by = 'cellname')
          
          colnames(coords2) <- c("Cell X Position", "Cell Y Position", "Cell ID", "Phenotype", "sample")
          coords2 <- coords2 %>% select("Cell ID", "Cell X Position", "Cell Y Position", Phenotype, sample) 
          rownames(coords2) <- coords2$'Cell ID'

          head(coords2)

          list_dist_mtx <- lapply(unique(coords2$sample), function(donor) {
              # donor <- unique(coords2$sample)[1]
              donor_distances <- coords2 %>% filter(sample == donor) 
              dist_mtx <- distance_matrix(donor_distances)
              return(dist_mtx)
          })

          names(list_dist_mtx) <- unique(coords2$sample)
      
          celltypes_comp <- as.character(unique(slide_xe$orig.celltype_cellchat))
          celltypes_comp <- c("CD8", "Mast")
          # celltypes_comp <- c("CD8") #, "Mast")

          lapply(celltypes_comp, function(from_celltype){ 
            # from_celltype <- celltypes_comp[1]
            cat(from_celltype, "\n")
            out_dir1<- paste0("spatial_redo/donut_surrounding/", slide_test, "/", from_celltype, "/")
            dir.create(paste0("spatial_redo/donut_surrounding/", slide_test, "/"))
            dir.create(out_dir1)

              ## Selecting gene lists
            files_name <- paste0("/home/fcastaneda/vd-vijay/donaldo/spatial/SaHe/slide_5k_SaHe2/re_segmentation/results/figures/crater/disease/merged_nocov_legend_3_dots/")
            to_dgea <- list.files(files_name, pattern="cpm_results", full.names = TRUE, recursive=TRUE)

            celltypes_to <- grep(from_celltype, unique(slide_xe$orig.celltype_cellchat), value=TRUE, invert=TRUE)
            celltype_dgea <- str_extract(to_dgea, paste0(celltypes_to, collapse="|") )
            celltypes_to<- celltypes_to[celltypes_to %in% celltype_dgea]

            celltypes_to <- c("basal", "secretory", "ciliated", "fibroblasts", "endothelial", "smooth_muscle", "pericytes")

            to_celltype<-celltypes_to

            matrix_donor <- lapply(names(list_dist_mtx), function(donor){   
                # donor <- names(list_dist_mtx)[1]
                
                cat(donor, "\n")
                out_dir1<- paste0(out_dir1, "/", donor, "/")
                dir.create(out_dir1)
                sample_matrix <- list_dist_mtx[[donor]] 

                slide_xe_temp <- subset(slide_xe, orig.donor == donor)

                from_celltype_cells <- slide_xe_temp@meta.data %>% filter(orig.celltype_cellchat == from_celltype  & orig.donor ==  donor) %>% pull(cellname)
                sample_matrix2 <- sample_matrix[, from_celltype_cells]
      
              ## donut approach
              
                out_dir1_donut <-paste0(out_dir1, "donut_approach/")
                dir.create(out_dir1_donut)

                sample_matrix2_temp <- sample_matrix2
                # sample_matrix2_temp[sample_matrix2_temp>200] <- NA
                try(rm(keep_cells_dist))
                keep_cells_dist<-c("no_cell")

                ### Correlation
                sample_matrix2_temp <- sample_matrix2
                # sample_matrix2_temp[sample_matrix2_temp>200] <- NA
                try(rm(keep_cells_dist))
                keep_cells_dist<-c("no_cell")

                bines_distances_corr <- lapply(c(0,50,100,150), function(dist){
                    # dist<-10
                    # cat(dist, "\n")
                    sample_matrix2_r<-rowSums(sample_matrix2_temp<dist,na.rm = TRUE)
                    # sample_matrix2_r <- rowSums(sample_matrix2<distances_ranges)
                    keep_cells <- names(sample_matrix2_r[sample_matrix2_r>0])
                    # keep_cells<- grep(paste(keep_cells_dist, collapse="|"), keep_cells, invert=TRUE, value=TRUE )
                    keep_cells<- keep_cells[!keep_cells %in% keep_cells_dist]
                    slide_xe_temp$keep_cells <- ifelse(slide_xe_temp$cellname %in% keep_cells, slide_xe_temp$orig.celltype_cellchat, "no_color")

                    slide_xe_temp$keep_cells <-ifelse(slide_xe_temp$orig.celltype_cellchat == from_celltype, from_celltype, slide_xe_temp$keep_cells)

                    distance_to_closest_cell_bin <- lapply(celltypes_to, function(to_celltype2){ 
                          # to_celltype2 <- celltypes_to[11]
                          cat(to_celltype2, "\n")

                          to_dgea_sub <- grep(to_celltype2, to_dgea, value=TRUE)

                          dgea_c <- lapply(to_dgea_sub, function(x){
                            # x <- to_dgea[1]
                            genes_select <- read.csv(x) %>% select(excel_name, log2FoldChange, padj, group)
                            genes_select$gene_name <- gsub("'", "", genes_select$excel_name)
                            genes_select$comp <- gsub(files_name, "", x)
                            genes_select$comp <- gsub("_.*|/", "", genes_select$comp)
                            return(genes_select)
                          })

                          dgea_c <- Reduce(rbind, dgea_c)

                          table(dgea_c$group, dgea_c$comp)
                          # table(dgea_c$log2FoldChange, dgea_c$group)
                          table(dgea_c$comp, useNA="always")

                          genes_severe_specific <- dgea_c %>% filter(group == "Severe.asthma") %>% group_by(gene_name) %>% summarize(rep=n()) %>% filter(rep >1 ) %>% pull(gene_name)

                          genes_list <- list("1st"=unique(as.character(genes_severe_specific)))
                          names(genes_list) <- to_celltype2
      
                          to_celltype_cells <- slide_xe_temp@meta.data %>% filter(orig.celltype_cellchat == to_celltype2  & orig.donor ==  donor & keep_cells != "no_color") %>% pull(cellname)

                          genes_ins <- genes_list[[to_celltype2]]
                          tmp2<-add_gene_tag(genes_ins, slide_xe_temp@meta.data, as.matrix(slide_xe_temp@assays$Xenium$data[genes_ins, ]), thresh=0.01, tag = c('tag', 'p', 'n'))
                          slide_xe_temp@meta.data = joindf(slide_xe_temp@meta.data, as.data.frame(tmp2))

                          slide_xe_temp$cellname <- rownames(slide_xe_temp@meta.data)

                          genes_to_focus <- c("cellname", paste0("tag_", genes_ins))

                          # to_celltype_cells <- slide_xe_temp@meta.data %>% filter( == to_celltype  & orig.donor ==  donor) %>% pull(cellname)
                        
                          if(length(to_celltype_cells) > 1) { # Greater than 1 or 0 ???
                            
                            sample_matrix2_1 <- sample_matrix2_temp[to_celltype_cells, ]
                            sample_matrix2_1_temp <- sample_matrix2_1<dist & sample_matrix2_1>0
                            ncells_surrounding <- rowSums(sample_matrix2_1_temp,na.rm = TRUE)

                            sample_matrix2_1[sample_matrix2_1>dist]<- NA

                            closest_cell <- rowMins(sample_matrix2_1, na.rm = TRUE)
                              # sample_matrix2_1[sample_matrix2_1<1]<- NA Line documented because in this vector the same cells shouldn't be included (same cells means the diagonal is not 0)
                            mean_distance <- rowMeans(sample_matrix2_1, na.rm = TRUE)
                            if(!all(names(closest_cell) == names(ncells_surrounding))) stop("Wrong assigantion of cell IDs")
                              
                                df_distances <- data.frame(
                                  row.names= names(ncells_surrounding),
                                  surrounding=ncells_surrounding,
                                  min_dist = closest_cell[names(ncells_surrounding)],
                                  mean_dist = mean_distance[names(ncells_surrounding)]
                                )

                            keep_meta <- slide_xe_temp@meta.data[to_celltype_cells,genes_to_focus]
                            class_cells <- apply(keep_meta[,genes_to_focus], 2, function(patt) grepl("p$", patt))
                            rownames(class_cells) <- rownames(keep_meta)
                            pos_cells2 <- rowSums(class_cells)[rowSums(class_cells) > 0]
                              # neg_cells <- rowSums(class_cells)[rowSums(class_cells) < 0]
                            keep_meta$signature <- ifelse(rownames(keep_meta) %in% names(pos_cells2), "sig_pos", "sig_neg")

                            keep_meta2 <- keep_meta %>% select(signature, all_of(genes_to_focus), cellname)

                            # keep_meta2
                            # keep_meta$chosen_gene <- keep_meta[,genes_to_focus]

                            radius_plot_gene <- merge(df_distances, keep_meta2, by="row.names", all=TRUE)
                            if(!all(nrow(radius_plot_gene) == nrow(df_distances))) stop("Wrong assigantion of cell IDs")
                            radius_plot_gene$donor<-donor
                            radius_plot_gene$to_celltype<-to_celltype2
                            radius_plot_gene$distance_bin <- dist
                            radius_plot_gene$Row.names <- NULL

                          } else{

                            df_distances <- data.frame(
                                row.names= "no_cells",
                                surrounding=NA,
                                min_dist = NA,
                                mean_dist = NA
                              )

                            for(i in c("signature", genes_to_focus)){
                            df_distances[,i] <- NA   
                            }
                            
                            df_distances$donor<-donor
                            df_distances$to_celltype<-to_celltype2
                            df_distances$distance_bin <- dist

                            radius_plot_gene<-df_distances

                          }
                          return(radius_plot_gene)
                    })

                    names(distance_to_closest_cell_bin) <- celltypes_to

                    colnames(sample_matrix2_temp)
                    keep_cells_dist<<-c(keep_cells_dist,keep_cells)

                    return(distance_to_closest_cell_bin)
                })

                names(bines_distances_corr) <- paste0("d", c(0,50,100,150))

                donut_corr_list <- lapply(celltypes_to, function(cellty){ 
                    # cellty<-celltypes_to[1]
                    bines_distances_celty <- lapply(bines_distances_corr, function(x) x[[cellty]]) 
                    return(Reduce(rbind, bines_distances_celty)) 
                })
                names(donut_corr_list) <- celltypes_to
                
                saveRDS(donut_corr_list, paste0(out_dir1_donut, "distances_corr.rds"))
                # saveRDS(radius_list, paste0(out_dir1_radius, "distances_complete.rds"))

            })

          })

      })

      ## Summary

      #  out_dir1<-  "spatial_redo/donut_surrounding/"
        files_disease_donut2 <- list.files(out_dir1, pattern="corr.rds", recursive=TRUE, full.names = TRUE)
        dir.create(paste0(out_dir1, "/summary_donut/"))

        # /mnt/bioadhoc/Groups/vd-vijay/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/spatial_redo/donut_surrounding//slide1/CD8/NIHHC_310/donut_approach/distances_corr.rds

        all_files <- lapply(files_disease_donut2, function(files_nei){ 
          # files_nei <- files_disease_donut2[65]
          cat(files_nei, "\n")
          nei_files <-readRDS(files_nei)
          nei_files_select <- lapply(nei_files, function(cellty) { 
            cellty <- cellty %>% select(cellname, surrounding, min_dist, mean_dist, donor, to_celltype, distance_bin) %>% filter(!is.na(cellname)) %>% filter(distance_bin == 50)
            return(cellty)
          })

          nei_files_select<- Reduce(rbind, nei_files_select)

          if(nrow(nei_files_select) == 0){
            summ_donor <- data.frame(
              donor = rep(gsub("|/donut_approach/distances_corr.rds", "", gsub("spatial_redo/donut_surrounding/.*/NIH", "NIH", files_nei)), 7),
              to_celltype= c("basal", "ciliated", "endothelial", "fibroblasts", "pericytes", "secretory", "smooth_muscle"),
              distance_bin = rep(50,7),
              mean = 0
            )
          } else { 
            summ_donor <- nei_files_select %>% group_by(donor, to_celltype, distance_bin) %>% summarize("mean"= mean(surrounding))
          }
          summ_donor$from <- str_extract(files_nei, c("CD8|Mast"))
          summ_donor$slide <- str_extract(files_nei, "slide1|slide2")
          return(summ_donor)
        })

        all_files_summ <- Reduce(rbind, all_files)
        all_files_summ$mean2 <- all_files_summ$mean
        all_files_summ$donor_slide <- paste0(all_files_summ$donor, "_", all_files_summ$slide) 
        table(all_files_summ$donor_slide)

        write.csv(dcast(all_files_summ %>% filter(from == "CD8"), donor_slide~to_celltype, value.var = "mean2"), "spatial_redo/donut_surrounding/summary/CD8_neighborhods.csv")
        write.csv(dcast(all_files_summ %>% filter(from == "Mast"), donor_slide~to_celltype, value.var = "mean2"), "spatial_redo/donut_surrounding/summary/MAST_neighborhods.csv")

      all_files_summ$disease <- str_extract(all_files_summ$donor, "HC|MA|W")

    ## Summary plots. Modify accordingly
      tmp.ggplot <- ggplot(data=all_files_summ %>% filter(disease == "HC"), aes(x=factor(from, levels=c("CD8", "Mast")), y=mean)) +
          geom_line(aes(group=donor_slide), color="black", linewidth=0.7) +
          geom_jitter(shape=1, stroke=2, size=3, width=0, height=0) +
          geom_boxplot(alpha=0.4, width=0.7, linewidth=3, fatten=2.5, outlier.shape=NA) +
          scale_y_continuous(breaks=scales::pretty_breaks(n=3)) + facet_wrap(~ to_celltype, ncol=7)
          # scale_color_manual(values=order_tissue_color) +
          # scale_fill_manual(values=order_tissue_color) +

      dir.create("spatial_redo/donut_surrounding/summary")
      pdf("spatial_redo/donut_surrounding/summary/cells_surr_HC.pdf", width=10, height=8)
      print(tmp.ggplot)
      dev.off()

