#!/usr/bin/R
# ---
# Author: Francisco Emmanuel Castaneda-Castro
# Other contributors: 
    # Donaldo Sosa-Garcia
# This code is an implementation of Cellchat's functions to perform cell-cell communication analysis
# Date: 2025-06-24
# Version 1.0.0

#### Set of functions to run cellchat
library(circlize)
library(reshape2)
library(ggforce)
library(ggthemes)
source("/home/donaldo/bioadhoc/SaHe/SaHe_bru/scripts/cellchat/ranknet_mod_function.R")
# source("/home/ciro/scripts/handy_functions/devel/plots.R")
plot_blank <- function(gp, lays = "ext", ...){ # Ciro's blank plot function
  if(!is.null(gp$facet)) gp <- gp + theme(strip.text.x = element_blank())
  if("patchwork" %in% class(gp) && isTRUE(length(gp$patches$plots) > 0)){
    for(i in 1:length(gp$patches$plots)){
      y <- gp$patches$plots[[i]]; if(length(y$layers) == 0) next
      gp$patches$plots[[i]] = plot_rm_layer(y, lays = lays, ...) + theme_axes()
    }
    gp
  }else{
    plot_rm_layer(gp, lays = lays, ...) + theme_axes()
  }
}

plot_rm_layer <- function(gp, lays = "ext", verbose = FALSE){
  if(verbose) cat("Layers:", sapply(gp$layers, function(x) class(x$geom)[1] ), "\n")
  tvar <- sapply(gp$layers, function(x) !grepl(lays, class(x$geom)[1], ignore.case = TRUE) )
  gp$layers <- gp$layers[tvar]
  gp
}

theme_axes <- function (
  base_size = 11,
  base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22
){
  half_line <- base_size / 2
  t <- theme(
    line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white",colour = "black", size = base_rect_size, linetype = 1),
    text = element_text(
      family = base_family, face = "plain",
      colour = "black", size = base_size, lineheight = 0.9,
      hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
      debug = FALSE
    ),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    legend.position = "none",
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    strip.text = element_blank(), strip.background = element_blank()
  )
  # ggplot_global$theme_all_null %+replace% t
}

### Cellchat pipeline

netVisual_chord_gene_mod <- function(object, slot.name = "net", color.use = NULL,
                                 signaling = NULL, pairLR.use = NULL, net = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                 thresh = 0.05, edge_colors=NULL,
                                 ...){

  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use) | sum(c("interaction_name","pathway_name") %in% colnames(pairLR.use) == 0)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }

  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }

  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")
    cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
    cols.common <- intersect(cols.default,colnames(object@LR$LRsig))
    pairLR = dplyr::select(object@LR$LRsig, cols.common)
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx,]
    net <- cbind(net, temp)
  }

  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }

  if (!is.null(pairLR.use)){
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }

  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  } else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  } else {
    targets.use <- levels(object@idents)
  }
  # remove the interactions with zero values
  net$receptor <- paste0(net$receptor, ".") #### Avoids CD99-CD99 problem
  cat("CD99 problem", "\n")

  df <- subset(net, prob > 0)

  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }

  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }

  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }

  cell.order.sources <- sources.use[sources.use %in% as.character(levels(object@idents))]
  cat("ORDER AYUDA SOURCES : ", cell.order.sources, "\n")
  cell.order.targets <- targets.use[targets.use %in% as.character(levels(object@idents))]
  cat("ORDER AYUDA targets: ", cell.order.targets, "\n")

  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]

  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])

  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)

  # define cell type color
  if (is.null(color.use)){
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  } else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  }

  # define edge color
  if(!is.null(edge_colors)){
      edge.color <- edge_colors[as.character(df.ordered.source$source)]
      names(edge.color) <- as.character(df.ordered.source$source)
  } else {
    edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  }

  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)

  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)

  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector

  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]

  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df.plot,
               order = order.sector, 
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
              link.arr.type = link.arr.type, big.gap = 20, small.gap = 0.3,annotationTrack = "grid",reduce = reduce,annotationTrackHeight = annotationTrackHeight,  preAllocateTracks = list(track.height = min(strwidth(order.sector))))
  
    circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = 0.5)
  }, bg.border = NA)

  
  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    lgd <- ComplexHeatmap::Legend(at = names(color.use), type = "grid", legend_gp = grid::gpar(fill = color.use), title = "Cell State")
    ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
  }

  circos.clear()


  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()

  return(gg)
}

netVisual_chord_gene_mod_no_legend_blank <- function(object, slot.name = "net", color.use = NULL,
                                 signaling = NULL, pairLR.use = NULL, net = NULL,
                                 sources.use = NULL, targets.use = NULL,
                                 lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                 link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                 transparency = 0.4, link.border = NA,
                                 title.name = NULL, legend.pos.x = 20, legend.pos.y = 20, show.legend = TRUE,
                                 thresh = 0.05, edge_colors=NULL,
                                 ...){


  if (!is.null(pairLR.use)) {
    if (!is.data.frame(pairLR.use) | sum(c("interaction_name","pathway_name") %in% colnames(pairLR.use) == 0)) {
      stop("pairLR.use should be a data frame with a signle column named either 'interaction_name' or 'pathway_name' ")
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      message("slot.name is set to be 'netP' when pairLR.use contains signaling pathways")
      slot.name = "netP"
    }
  }

  if (!is.null(pairLR.use) & !is.null(signaling)) {
    stop("Please do not assign values to 'signaling' when using 'pairLR.use'")
  }

  if (is.null(net)) {
    prob <- slot(object, "net")$prob
    pval <- slot(object, "net")$pval
    prob[pval > thresh] <- 0
    net <- reshape2::melt(prob, value.name = "prob")
    colnames(net)[1:3] <- c("source","target","interaction_name")
    cols.default <- c("interaction_name_2", "pathway_name",  "ligand",  "receptor" ,"annotation","evidence")
    cols.common <- intersect(cols.default,colnames(object@LR$LRsig))
    pairLR = dplyr::select(object@LR$LRsig, cols.common)
    idx <- match(net$interaction_name, rownames(pairLR))
    temp <- pairLR[idx,]
    net <- cbind(net, temp)
  }

  if (!is.null(signaling)) {
    pairLR.use <- data.frame()
    for (i in 1:length(signaling)) {
      pairLR.use.i <- searchPair(signaling = signaling[i], pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
      pairLR.use <- rbind(pairLR.use, pairLR.use.i)
    }
  }

  if (!is.null(pairLR.use)){
    if ("interaction_name" %in% colnames(pairLR.use)) {
      net <- subset(net,interaction_name %in% pairLR.use$interaction_name)
    } else if ("pathway_name" %in% colnames(pairLR.use)) {
      net <- subset(net, pathway_name %in% as.character(pairLR.use$pathway_name))
    }
  }

  if (slot.name == "netP") {
    net <- dplyr::select(net, c("source","target","pathway_name","prob"))
    net$source_target <- paste(net$source, net$target, sep = "sourceTotarget")
    net <- net %>% dplyr::group_by(source_target, pathway_name) %>% dplyr::summarize(prob = sum(prob))
    a <- stringr::str_split(net$source_target, "sourceTotarget", simplify = T)
    net$source <- as.character(a[, 1])
    net$target <- as.character(a[, 2])
    net$ligand <- net$pathway_name
    net$receptor <- " "
  }

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- levels(object@idents)[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  } else {
    sources.use <- levels(object@idents)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- levels(object@idents)[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  } else {
    targets.use <- levels(object@idents)
  }
  # remove the interactions with zero values
  net$receptor <- paste0(net$receptor, ".") #### Avoids CD99-CD99 problem
  df <- subset(net, prob > 0)
  
  if (nrow(df) == 0) {
    stop("No signaling links are inferred! ")
  }

  if (length(unique(net$ligand)) == 1) {
    message("You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway")
  }

  df$id <- 1:nrow(df)
  # deal with duplicated sector names
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(' ',j-1),collapse = ''))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(' ',j-1),collapse = ''))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }

  cell.order.sources <- sources.use[sources.use %in% as.character(levels(object@idents))]
  cat("ORDER AYUDA SOURCES : ", cell.order.sources, "\n")
  cell.order.targets <- targets.use[targets.use %in% as.character(levels(object@idents))]
  cat("ORDER AYUDA targets: ", cell.order.targets, "\n")

  df$source <- factor(df$source, levels = cell.order.sources)
  df$target <- factor(df$target, levels = cell.order.targets)
  # df.ordered.source <- df[with(df, order(source, target, -prob)), ]
  # df.ordered.target <- df[with(df, order(target, source, -prob)), ]
  df.ordered.source <- df[with(df, order(source, -prob)), ]
  df.ordered.target <- df[with(df, order(target, -prob)), ]

  order.source <- unique(df.ordered.source[ ,c('ligand','source')])
  order.target <- unique(df.ordered.target[ ,c('receptor','target')])

  # define sector order
  order.sector <- c(order.source$ligand, order.target$receptor)

  # define cell type color
  if (is.null(color.use)){
    color.use = scPalette(nlevels(object@idents))
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  } else if (is.null(names(color.use))) {
    names(color.use) <- levels(object@idents)
    color.use <- color.use[levels(object@idents) %in% as.character(union(df$source,df$target))]
  }

  # define edge color
  # define edge color
  if(!is.null(edge_colors)){
      edge.color <- edge_colors[as.character(df.ordered.source$source)]
      names(edge.color) <- as.character(df.ordered.source$source)
  } else {
    edge.color <- color.use[as.character(df.ordered.source$source)]
  names(edge.color) <- as.character(df.ordered.source$source)
  }

  # define grid colors
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- as.character(order.source$source)
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- as.character(order.target$target)
  grid.col <- c(as.character(grid.col.ligand), as.character(grid.col.receptor))
  names(grid.col) <- order.sector

  df.plot <- df.ordered.source[ ,c('ligand','receptor','prob')]

  if (directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df.plot,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
              link.arr.type = link.arr.type, big.gap = 20, small.gap = 0.3,annotationTrack = "grid",reduce = reduce,annotationTrackHeight = annotationTrackHeight,  preAllocateTracks = list(track.height = min(strwidth(order.sector))))
  
  circos.clear()

  if(!is.null(title.name)){
    text(-0, 1.02, title.name, cex=1)
  }
  gg <- recordPlot()

  return(gg)
}

cellchat_comparisions<-function(
      cellchat_objects=cellchat_objects,
      groups_to_compare= order_compare,
      padj=0.05,
      fc=0.25,
      method = "wilcoxon",
      output_dir = outdir,
      senders_celltypes=senders_cellchat,
      receivers_celltypes=receivers_cellchat,
      colors_celltypes=trm_colors_cellchat,
      color_edges=trm_colors_cellchat2, 
      per_celltypes = TRUE
      ){ 

      ## Atention padj is really p-value. I need to change it

      dir.create(output_dir)

      library(ComplexHeatmap)
      library("tidyr")

      outdir2 <- paste0(output_dir, "/fc", fc, "_p" ,padj)
      dir.create(outdir2)
      
      cellchat_all <- mergeCellChat(cellchat_objects, add.names = names(cellchat_objects))

      ### Performing pair-wise comparisions
        combs_num <- combn(names(cellchat_objects),2,simplify=FALSE)
        names(combs_num)<-unlist(lapply(combs_num, function(x){
                    # x<-combs_num[1]
                    paste(unlist(x), collapse=".vs.")
                  }))

      lapply(names(combs_num), function(com){ 

        # com <- names(combs_num)[1]
        cat("Performing analysis for ", com, "comparision \n")
        outdir_com <- paste0(outdir2, "/", com, "/")
        dir.create(outdir_com)

        celltype_dgea_imm <- unique(cellchat_all@meta$orig.celltype_cellchat)
        senders_celltypes <- senders_celltypes[senders_celltypes %in% celltype_dgea_imm]
        receivers_celltypes <- receivers_celltypes[receivers_celltypes %in% celltype_dgea_imm]

        # Whether the cell-cell communication is enhanced or not  
        # The interaction between which cell types is significantly changed
        # How the major sources and targets change from one condition to another

        pdf(paste0(outdir_com, "differential_circus.pdf"), width=10, height=10)
        #where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, arrow.size = 0.2, vertex.label.color = "black")

        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, measure = "weight", arrow.size = 0.2, vertex.label.color = "black")
        
        dev.off()

        pdf(paste0(outdir_com, "differential_circus_onlySenders.pdf"), width=10, height=10)
        #where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, vertex.label.color = "#00000000", arrow.size = 0.2, sources.use = senders_celltypes, targets.use=receivers_celltypes)
        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, measure = "weight", vertex.label.color = "#00000000", arrow.size = 0.2, sources.use = senders_celltypes, targets.use=receivers_celltypes)
        dev.off()

        pdf(paste0(outdir_com, "differential_circus_onlySenders-labels.pdf"), width=10, height=10)
        #where red (or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
        par(mfrow = c(1,1), xpd=TRUE)
        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, vertex.label.color = "black", arrow.size = 0.2, sources.use = senders_celltypes, targets.use=receivers_celltypes)
        netVisual_diffInteraction(cellchat_all, comparison=rev(combs_num[[com]]), vertex.weight = 20,  weight.scale = T, color.use=colors_celltypes[levels(cellchat_all@idents$joint)], top=1, label.edge = FALSE, measure = "weight", vertex.label.color = "black", arrow.size = 0.2, sources.use = senders_celltypes, targets.use=receivers_celltypes)
        dev.off()

        pdf(paste0(outdir_com, "differential_heatmap.pdf"), width=15, height=10)
        # Therefore, the bar height indicates the degree of change in terms of the number of interactions or interaction strength between the two conditions. In the colorbar, red (or blue) represents increased (or decreased) signaling in the second dataset compared to the first one.
        gg1 <- netVisual_heatmap(cellchat_all, comparison=rev(combs_num[[com]]))
        #> Do heatmap based on a merged object
        gg2 <- netVisual_heatmap(cellchat_all, measure = "weight", comparison=rev(combs_num[[com]]))
        #> Do heatmap based on a merged object
        print(gg1 + gg2)
        dev.off()

        pdf(paste0(outdir_com, "barplot_signals.pdf"), width=15, height=10)
        gg1 <- rankNet(cellchat_all, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

        gg2 <- rankNet(cellchat_all, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
        print(gg1 + gg2)
        dev.off()

        pgrouop<-combs_num[[com]][1] #posgroup v 
        ngrouop<-combs_num[[com]][2] #i ngrouop
    
        pdf(paste0(outdir_com, "incomming_outgoing_pathways.pdf"), width = 15, height = 16)

          pathway.union <- union(cellchat_objects[[pgrouop]]@netP$pathways, cellchat_objects[[ngrouop]]@netP$pathways)

          ht1 = netAnalysis_signalingRole_heatmap(cellchat_objects[[pgrouop]], pattern = "outgoing", signaling = pathway.union, title = pgrouop, width = 15, height = 16)
          
          ht2 = netAnalysis_signalingRole_heatmap(cellchat_objects[[ngrouop]], pattern = "outgoing", signaling = pathway.union, title = ngrouop, width = 15, height = 16)

          draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

          ht1 = netAnalysis_signalingRole_heatmap(cellchat_objects[[pgrouop]], pattern = "incoming", signaling = pathway.union, title = pgrouop, width = 15, height = 16, color.heatmap = "GnBu")
          ht2 = netAnalysis_signalingRole_heatmap(cellchat_objects[[ngrouop]], pattern = "incoming", signaling = pathway.union, title = ngrouop, width = 15, height = 16, color.heatmap = "GnBu")
          draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
          
        dev.off()

        # Identify dysfunctional signaling by using differential expression analysis

        # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
        pos.dataset = pgrouop
        # define a char name used for storing the results of differential expression analysis
        features.name = paste0(pos.dataset, ".merged")
        rm(cellchat_dgea)
        cellchat_dgea <- mergeCellChat(cellchat_objects[c(pgrouop,ngrouop)], add.names = c(pgrouop,ngrouop))
      
        # perform differential expression analysis 
        # Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

        cellchat_dgea <- identifyOverExpressedGenes(cellchat_dgea, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = fc,thresh.p = padj, group.DE.combined = FALSE) 

        dir.create(paste0(outdir_com, "/dgea/"))
        # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
        net <- netMappingDEG(cellchat_dgea, features.name = features.name, variable.all = TRUE, thresh=padj)
        write.csv(net, paste0(outdir_com, "/dgea/net.csv"))
        # extract the ligand-receptor pairs with upregulated ligands in LS

        lapply(c("Senders", "Receivers"), function(direction_senders){ 
          # direction_senders <- "Senders"
            cat("########", direction_senders, "\n")
            senders_celltypes_dir <- if (direction_senders == "Senders") senders_celltypes else receivers_celltypes
            receivers_celltypes_dir <- if (direction_senders == "Senders") receivers_celltypes else senders_celltypes

            if(per_celltypes){
              lapply(receivers_celltypes_dir, function(celtype_target){
                  # celtype_target <- receivers_celltypes_dir[1]
                  # cellchat_dgea
                  cat("######", celtype_target, "\n")

                  net.up <- tryCatch({(subsetCommunication(cellchat_dgea, net = net, datasets =  pgrouop,ligand.logFC = fc, receptor.logFC = NULL, sources.use = senders_celltypes_dir, targets.use = celtype_target))}, error =  function(err) { cat("No significant pathway found"); return(NA)} )
                  # extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
                  net.down <- tryCatch({(subsetCommunication(cellchat_dgea, net = net, datasets =  ngrouop,ligand.logFC = -fc, receptor.logFC = NULL, sources.use = senders_celltypes_dir, targets.use = celtype_target))}, error =  function(err) { cat("No significant pathway found"); return(NA)} )

                  # gene.up <- extractGeneSubsetFromPair(net.up, cellchat_dgea)
                  # gene.down <- extractGeneSubsetFromPair(net.down, cellchat_dgea)
                try(lapply(c("up", "down"), function(dgea_res){
                  # dgea_res<- "up"
                  cat("###", dgea_res, "\n")

                  title_plot <- paste0(dgea_res, " in ", pos.dataset)
                  result_table <- if(dgea_res == "up") net.up else net.down
                  pairLR.use.up = result_table[, "interaction_name", drop = F]
                  # cellchat_objects
                  to_plot <- tryCatch({(netVisual_bubble(cellchat_all, pairLR.use = pairLR.use.up, sources.use = senders_celltypes_dir, targets.use = celtype_target, comparison = 1:length(cellchat_objects),  angle.x = 90, remove.isolate = F,title.name = paste0("Up-regulated signaling in ", pgrouop),  return.data = TRUE))}, error =  function(err) { cat("No significant pathway found \n"); return(NA)} )
                  # if(is.na(to_plot)) return(0)

                  # pdf(paste0(outdir_com, "/dgea/", celtype_target, "_checking.pdf"), width=12, height=12)
                  # p<-netVisual_bubble(cellchat_all, pairLR.use = pairLR.use.up, sources.use = senders_celltypes, targets.use = celtype_target, comparison = 1:length(cellchat_objects),  angle.x = 90, remove.isolate = F, grid.on= T, title.name = paste0("Up-regulated signaling in ",names(cellchat_objects)[pgrouop]))

                  # p1<-netVisual_bubble(cellchat_dgea, pairLR.use = pairLR.use.up, sources.use = senders_celltypes, targets.use = celtype_target, comparison = c(1,2),  angle.x = 90, remove.isolate = F, grid.on= T, title.name = paste0("Up-regulated signaling in ",names(cellchat_objects)[pgrouop]))
                  # print(p+p1)
                  # dev.off()

                  hjust.x = 1
                  hjust.y = 0.5
                  vjust.x = 0.3
                  angle.x = 90
                  to_plot <- to_plot$communication # Do not eliminate. It gives the error when to_plot do not have pathways
                  result_table$to_filter <- paste0(result_table$source, "_", result_table$target, "_", result_table$ligand, "_", result_table$receptor)
                  to_plot$to_filter <- paste0(to_plot$source
                  , "_", to_plot$target, "_", to_plot$ligand, "_", to_plot$receptor)

                  # to_plot %>% filter(!to_filter %in% result_table$to_filter)
                  to_plot$gray <- ifelse(is.na(to_plot$prob), 0,1)
                  to_plot$prob2 <- ifelse(to_plot$to_filter %in% result_table$to_filter, to_plot$prob, NA)

                  to_plot_order_y <- to_plot %>% arrange(desc(prob), pathway_name) %>% pull(interaction_name_2) %>% unique()
                  dataset_comp <- unique(to_plot$dataset)
                  source_cellchat_l <- length(unique(senders_celltypes_dir))
                  order_source.target <- to_plot %>% arrange(dataset, factor(source, levels=senders_celltypes_dir)) %>% pull(source.target) %>% unique()
                  y_axis_size <- ifelse(length(to_plot_order_y) < 3, 5, to_plot_order_y/1.5)

                  xintercept = seq(0.5+source_cellchat_l, source_cellchat_l*length(dataset_comp), by = source_cellchat_l)
                  library(RColorBrewer)
                  dir.create(paste0(outdir_com, "dgea/"))
                  write.csv(to_plot, paste0(outdir_com, "/dgea/", celtype_target, "-", dgea_res, "_", direction_senders, ".csv"))
                  write.csv(result_table, paste0(outdir_com, "/dgea/DGEA_table_", celtype_target, "-", dgea_res, "_", direction_senders, ".csv"))

                  
                  pdf(paste0(outdir_com, "/dgea/", celtype_target, "-", dgea_res, "_", direction_senders, "_emma.pdf"), width=source_cellchat_l, height=y_axis_size)
                  
                      if(sum(!is.na(to_plot$prob2)) > 1) { 
                      print(ggplot(to_plot, aes(x = factor(source.target, levels=order_source.target), y = factor(interaction_name_2, levels=rev(to_plot_order_y)), color = prob2, size = pval, alpha=gray)) +
                      geom_point(pch = 16) +
                      # geom_point(pch = 1, color="gray90") +
                      theme_linedraw() + theme(panel.grid.major = element_blank()) +
                      theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank()) +
                      scale_x_discrete(position = "bottom") + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = 0.5) + 
                      scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_plot$prob2, 0,na.rm= T), quantile(to_plot$prob2, 1,na.rm= T)), breaks = c(quantile(to_plot$prob2, 0,na.rm= T), quantile(to_plot$prob2, 1,na.rm= T)), labels = c("min","max")) +
                      scale_alpha_identity() +
                      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + labs(title=title_plot)
                      )
                  } else {
                      print(ggplot(to_plot, aes(x = factor(source.target, levels=order_source.target), y = factor(interaction_name_2, levels=rev(to_plot_order_y)), color = prob2, size = pval, alpha=gray)) +
                      geom_point(pch = 16) +
                      # geom_point(pch = 1, color="gray90") +
                      theme_linedraw() + theme(panel.grid.major = element_blank()) +
                      theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                          axis.title.x = element_blank(),
                          axis.title.y = element_blank()) +
                      scale_x_discrete(position = "bottom") + geom_vline(xintercept = xintercept, linetype="dashed", color = "grey60", size = 0.5)  +
                      scale_alpha_identity() +
                      guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) + labs(title=title_plot)
                      )
                  }
                  dev.off()
                  }))
              })
            } 

            net.up <- tryCatch({(subsetCommunication(cellchat_dgea, net = net, datasets =  pgrouop,ligand.logFC = fc, receptor.logFC = NULL, sources.use = senders_celltypes_dir, targets.use = receivers_celltypes_dir))}, error =  function(err) { cat("No significant pathway found"); return(NA)} )
                # extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
            net.down <- tryCatch({(subsetCommunication(cellchat_dgea, net = net, datasets =  ngrouop,ligand.logFC = -fc, receptor.logFC = NULL, sources.use = senders_celltypes_dir, targets.use = receivers_celltypes_dir))}, error =  function(err) { cat("No significant pathway found"); return(NA)} )

                # par(mfrow = c(1,2), xpd=TRUE)
                dir.create(paste0(outdir_com, "/dgea/chords/"))
                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_mod.pdf"), width=10, height=10)
                try(netVisual_chord_gene(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene(cellchat_objects[[ngrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_mod_ALL_gap.pdf"), width=10, height=10)
                try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_mod_ALL_gap_GRAY.pdf"), width=10, height=10)
                try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1, color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1, color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Down-regulated signaling in ", pgrouop)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_mod_ALL_gap_blank.pdf"), width=10, height=10)
                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = celltype_dgea_imm, targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_mod_ALL_gap_GRAY_blank.pdf"), width=10, height=10)
                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes,edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = celltype_dgea_imm, targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges,title.name = paste0("Down-regulated signaling in ", pgrouop)))

                dev.off()

                ##Donaldo's modification Dec 13 2025
                if(direction_senders == "Senders"){

                lapply(seq_along(senders_celltypes_dir), function(per_sender) {
                    # per_sender <- 1
                    per_sender_ct <- senders_celltypes_dir[[per_sender]]
                    per_receiver_ct_specific <- receivers_celltypes[[per_sender]]

                    #per celltype to group
                    pdf(paste0(outdir_com, "/dgea/chords/", "Senders-chord", "to", per_receiver_ct_specific, "_ALL_gap.pdf"), width=10, height=10)
                    print(try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = senders_celltypes_dir, targets.use = per_receiver_ct_specific, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop))))
                
                    print(try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = senders_celltypes_dir, targets.use = per_receiver_ct_specific, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop))))

                    dev.off()

                    pdf(paste0(outdir_com, "/dgea/chords/", "Senders-chord", "to", per_receiver_ct_specific, "_ALL_gap_blank.pdf"), width=10, height=10)
                    print(try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = senders_celltypes_dir, targets.use = per_receiver_ct_specific, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop))))
                
                    print(try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = senders_celltypes_dir, targets.use = per_receiver_ct_specific, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop))))

                    dev.off()

                  }
                )

                lapply(seq_along(senders_celltypes_dir), function(per_sender) {
                    # per_sender <- 1
                    per_sender_ct <- senders_celltypes_dir[[per_sender]]

                  lapply(seq_along(receivers_celltypes_dir), function(per_receiver) {
                    # per_receiver <- 3
                    per_receiver_ct <- receivers_celltypes_dir[[per_receiver]]

                    #per celltype to celltype
                    pdf(paste0(outdir_com, "/dgea/chords/", per_sender_ct, "to", per_receiver_ct, "_ALL_gap.pdf"), width=10, height=10)
                    print(try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = per_sender_ct, targets.use = per_receiver_ct, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop))))
                
                    print(try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = per_sender_ct, targets.use = per_receiver_ct, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop))))

                    dev.off()

                    pdf(paste0(outdir_com, "/dgea/chords/", per_sender_ct, "to", per_receiver_ct, "_ALL_gap_blank.pdf"), width=10, height=10)
                    print(try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = per_sender_ct, targets.use = per_receiver_ct, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Up-regulated signaling in ", pgrouop))))
                
                    print(try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = per_sender_ct, targets.use = per_receiver_ct, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes ,title.name = paste0("Down-regulated signaling in ", pgrouop))))

                    dev.off()

                    pdf(paste0(outdir_com, "/dgea/chords/", per_sender_ct, "to",per_receiver_ct, "-chord_mod_ALL_gap_GRAY.pdf"), width=10, height=10)
                    print(try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = rev(per_sender_ct), targets.use = per_receiver_ct, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1, color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop))))
                    
                    try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = rev(per_sender_ct), targets.use = per_receiver_ct, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1, color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Down-regulated signaling in ", pgrouop)))

                    dev.off()

                    pdf(paste0(outdir_com, "/dgea/chords/", per_sender_ct, "to", per_receiver_ct, "-chord_mod_ALL_gap_GRAY_blank.pdf"), width=10, height=10)
                    try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = rev(per_sender_ct), targets.use = per_receiver_ct, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes,edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                    
                    try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = per_sender_ct, targets.use = per_receiver_ct, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges,title.name = paste0("Down-regulated signaling in ", pgrouop)))

                    dev.off()
                    
                    return(NULL)
                  })

                  return(NULL)

                })

              }

              ## End od Donaldo's modification

                net.up <- net.up %>% arrange(desc(prob))
                first_quartile <- quantile(net.up$prob, probs = 0.8)
                net.up2 <- net.up %>% filter(prob > first_quartile)

                net.down <- net.down %>% arrange(desc(prob))
                first_quartile <- quantile(net.down$prob, probs = 0.8)
                net.down2 <- net.down %>% filter(prob > first_quartile)

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_filtered_20q.pdf"), width=10, height=10)
                try(netVisual_chord_gene(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up2, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene(cellchat_objects[[ngrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down2, lab.cex = 0.8, small.gap = 3.5, transparency=0.1,color.use=colors_celltypes, title.name = paste0("Down-regulated signaling in ", pgrouop,scale = FALSE)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_filtered_gap_20q.pdf"), width=10, height=10)
                try(netVisual_chord_gene_mod(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up2, lab.cex = 0.8, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod(cellchat_objects[[ngrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down2, lab.cex = 0.8, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Down-regulated signaling in ", pgrouop,scale = FALSE)))

                dev.off()

                pdf(paste0(outdir_com, "/dgea/chords/", direction_senders, "-chord_filtered_gap_blank_20q.pdf"), width=10, height=10)

                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[pgrouop]], sources.use = rev(senders_celltypes_dir), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.up2, lab.cex = 0.8, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Up-regulated signaling in ", pgrouop)))
                
                try(netVisual_chord_gene_mod_no_legend_blank(cellchat_objects[[ngrouop]], sources.use = rev(celltype_dgea_imm), targets.use = receivers_celltypes_dir, slot.name = 'net', net = net.down2, lab.cex = 0.8, transparency=0.1,color.use=colors_celltypes, edge_colors=color_edges, title.name = paste0("Down-regulated signaling in ", pgrouop,scale = FALSE)))

                dev.off()
                
                
          #> You may try the function `netVisual_chord_cell` for visualizing individual signaling pathway
      
        })
       
      })
}

running_cellchat<-function(
      object_cellchat = brushings_cellchat,
      pop_size = T, 
      pval_cells=pval,
      fc_cells=fc_cells,
      mean_method = mean_method,
      dgea_method = dgea_method,
      comparing_groups = comparing_groups,
      column_group= "orig.disease",
      identities_clusters = "orig.celltype_cellchat",
      donor_info = "orig.donor_paper",
      outdir=outdir,
      senders_celltypes=source_cellchat,
      celltype_selected = "Tcells",
      receivers_celltypes=receivers_celltypes,
      colors_celltypes=trm_colors_cellchat,
      color_edges=NULL, 
      include_non_protein_signalling=FALSE,
      exp_type = list(type="Standar", params=NULL) # if spatial use exp_type = list(type="Spatial", params="distance.use = TRUE, raw.use = FALSE, contact.range = 10, contact.dependent = TRUE, interaction.range = 150, scale.distance = 1") or add any parameter for the computeCommunProb function

 ){ 

    #  object_cellchat = SaHe_all_slide2_2
    #   pop_size = TRUE 
    #   pval_cells=0.05
    #   fc_cells=0
    #   mean_method = list("truncatedMean", 0.05)
    #   dgea_method = list("cellchat")
    #   comparing_groups = c("Severe.asthma", "Mild.asthma", "Healthy")
    #   outdir="/mnt/bioadhoc/Groups/vd-vijay/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/spatial_cellchat/slide2_new_seg/"
    #   senders_celltypes=source_cellchat  
    #   celltype_selected = "CD8"
    #   receivers_celltypes=names(pop_cl_colors)[!names(pop_cl_colors) %in% source_cellchat]
    #   colors_celltypes=pop_cl_colors
    #   include_non_protein_signalling=FALSE
    #   column_group= "orig.disease"
    #   identities_clusters = "orig.celltype_cellchat"
    #   exp_type = list(type="Spatial", params="distance.use = TRUE, raw.use = FALSE, contact.range = 10, contact.dependent = TRUE, interaction.range = 250, scale.distance = 0.5")
    #   donor_info = "orig.donor"
      
      # library(cellchat)
      library(CellChat)
      library(patchwork)
      library(dplyr)
      library(crayon)
      options(stringsAsFactors = FALSE)
      ## CHECKING THAT CONDITIONS ARE PRESENT IN THE CELLCHAT_OBJECT
      dir.create(outdir)
      cat(red("Attention. This function is only compatible with HUMAN data. If you want to use mouse, please modify it. It shouldn't be difficult, sorry. \n"))

      brushings_cellchat = object_cellchat ## Sorry for the variable name. This code was first developed for the "CD8 BRUSHINGS, asthma Herrera de la Mata paper"
      brushings_cellchat$orig.celltype_cellchat <- brushings_cellchat@meta.data[, identities_clusters]
      brushings_cellchat$samples <- brushings_cellchat@meta.data[, donor_info]

      cat("Using donors : ", names(table(brushings_cellchat$samples)) , "\n")

      folder_suffix <- paste0("DGEA_",unlist(dgea_method)[1], "_fc", fc_cells, "p", pval_cells, "MEAN-", paste0(unlist(mean_method), collapse="_"), "-pop_size_", as.character(pop_size), "_ProteinSignaling", include_non_protein_signalling)

      # folder_suffix <- paste0("DGEA_",unlist(dgea_method)[1], "_fc", fc_cells, "p", pval_cells, "MEAN-", paste0(unlist(mean_method), collapse="_"), "-pop_size_", as.character(pop_size)) Just run this line if you want to run previous versions
      folder_suffix <- paste0(folder_suffix, "_TYPE-", gsub(",|=", "_", gsub("\\s", "", paste0(unlist(exp_type), collapse = "-"))))
      # outdir <- paste0(outdir, "/", folder_suffix, "/SA")
      outdir <- paste0(outdir, "/", folder_suffix)
      dir.create(outdir)

      cat("Working in: ", red(outdir), "\n")

      comparing_groups <- comparing_groups[comparing_groups %in% unique(brushings_cellchat@meta.data[,column_group])]
    
      nCores <- 2

      if(!file.exists(paste0(outdir, "/hola.csv"))){

        cellchat_objects <-  parallel::mclapply(comparing_groups, function(disease_cellchat){
         
          # disease_cellchat <- comparing_groups[1]
          outdir <- paste0(outdir, "/", disease_cellchat, "/")
          dir.create(outdir)
          brushings_cellchat_SA<- subset(brushings_cellchat, !!sym(column_group) == disease_cellchat)

          cat("Performing cellchat for", unique(brushings_cellchat_SA@meta.data[,column_group]), "\n")
          table(brushings_cellchat_SA$orig.celltype_cellchat, useNA="always")

          if(exp_type[[1]] == "Spatial") { 
            cat("Creating spatial object ... \n")
            spatial.locs = GetTissueCoordinates(brushings_cellchat_SA[["fov"]], scale = NULL)
            data.input = GetAssayData(brushings_cellchat_SA, slot = "data", assay = "Xenium")

            table(colnames(data.input) %in% spatial.locs$cell)
            table(spatial.locs$cell %in% colnames(data.input))

            rownames(spatial.locs) <- spatial.locs$cell
            
            Idents(brushings_cellchat_SA) <- identities_clusters

            meta <- brushings_cellchat_SA@meta.data # manually create a dataframe consisting of the cell labels

            all(spatial.locs$cell %in% rownames(meta))
            all(rownames(meta) %in% spatial.locs$cell)

            spatial.locs$cell<-NULL

            conversion.factor = 1
            d = computeCellDistance(spatial.locs)
            spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
            spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
            saveRDS(spatial.factors, file = paste0(outdir,"/spatial_factors_colon.rds"))

            cat(red("Creating SPATIAL cellchat object \n"))
            cellchat <- createCellChat(object = data.input, meta = meta, group.by = identities_clusters, datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)
   
          } else { 
          cat(red("Creating STANDAR single cell cellchat object \n"))
          cellchat <- createCellChat(object = brushings_cellchat_SA, group.by = "orig.celltype_cellchat", assay = "RNA")
          }
          ## Data base

          CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
          
          pdf(paste0(outdir, "DB_categories.pdf"))
          print(showDatabaseCategory(CellChatDB))
          dev.off()

          # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
          if(isTRUE(include_non_protein_signalling)) CellChatDB.use <- CellChatDB else CellChatDB.use <- subsetDB(CellChatDB)
          
          cellchat@DB <- CellChatDB.use

          #By default CellChat uses the raw data (i.e., object@data.signaling) instead of the projected data. To use the projected data, users should run the function projectData before running computeCommunProb, and then set raw.use = FALSE when running computeCommunProb.

          cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

          cellchat <- updateCellChat(cellchat) # Just in case. 

          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")

          # future::plan("multisession", workers = 1) # do parallel

           if(unlist(dgea_method)[1] == "MAST") { 
            cat("Bypassing cellchat DGEA methods, instead, the genes provided on the second argument will be used. \n ")

            genes_to_use.info <- dgea_method[[2]] %>% filter(logFC >= fc_cells, pvalues.adj < pval_cells, features %in% rownames(cellchat@data.signaling), disease == disease_cellchat)
            genes_to_use.info[,c("disease","celltypes")] <- NULL
            genes_to_use <- unique(genes_to_use.info$features)

            cellchat@var.features$features.info <- genes_to_use.info
            cellchat@var.features$features <- genes_to_use
            
            cat("Using only: ", length(genes_to_use), " genes to perform the cellchat phase1? analysis (group independent) \n")

           } else {
              cellchat <- identifyOverExpressedGenes(cellchat,  thresh.fc = fc_cells, thresh.p =pval_cells)
           }

          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")
          # cat("Here222 \n")

          cat(red("Differentially expressed genes = ", length(unique(cellchat@var.features$features)), "out of", length(unique(c(CellChatDB.use$interaction$ligand, CellChatDB.use$interaction$receptor)))), "\n")
    

          cellchat <- identifyOverExpressedInteractions(cellchat)

          summary_genes_used<-list(
          "ngenes"= length(unique(cellchat@var.features$features)),
          "DEG_genes"= unique(cellchat@var.features$features),
          "ligands_used" = unique(cellchat@LR$LRsig$ligand),
          "receptor_used" = unique(cellchat@LR$LRsig$receptor),
          "total_ligands_in_DB"= unique(CellChatDB.use$interaction$ligand), 
          "total_receptors_in_DB"= unique(CellChatDB.use$interaction$receptor))

          paste0("DEG_", summary_genes_used$ngenes, "_used_", length(unique(c(summary_genes_used$ligands_used, summary_genes_used$receptor_used))), "_outof_", length(unique(c(summary_genes_used$total_ligands_in_DB, summary_genes_used$total_receptors_in_DB))))
            
          
          saveRDS(summary_genes_used, paste0(outdir, "summary_DEG_", summary_genes_used$ngenes, "_used_", length(unique(c(summary_genes_used$ligands_used, summary_genes_used$receptor_used))), "_outof_", length(unique(c(summary_genes_used$total_ligands_in_DB, summary_genes_used$total_receptors_in_DB))), ".rds"))

          if(exp_type[[1]] == "Spatial"){ 
            cat(red("Running Spatial: smoothData computeCommunProb ... \n"))
            # cellchat <- smoothData(cellchat, adj = PPI.human)

              # cat("Here222 \n")
              # cat("Here222 \n")
              # cat("Here222 \n")
              # cat("Here222 \n")
              # cat("Here222 \n")

            if(unlist(mean_method)[1] == "truncatedMean"){               
              compute_prob_args <- paste0("cellchat <- computeCommunProb(cellchat, type = unlist(mean_method)[1], trim= as.numeric(unlist(mean_method)[2]) ,population.size = pop_size, ", exp_type[[2]], ")" )
              cat(red("Using: \n ", compute_prob_args), "\n")
              eval(parse(text = compute_prob_args))
            } else {
              compute_prob_args <- paste0("cellchat <- computeCommunProb(cellchat, type = unlist(mean_method)[1], population.size = pop_size, ", exp_type[[2]], ")")
              cat(red("Using: \n ", compute_prob_args), "\n")
              eval(parse(text = compute_prob_args))
            }
          } else {   
            cat(red("Running standar single cell cellchat: computeCommunProb \n"), 
              if(!is.null(exp_type[["params"]])){ 
               paste0(red("Following parameters: will be ignored in computeCommunProb "), exp_type[[2]], "\n") 
              })
          
            if(unlist(mean_method)[1] == "truncatedMean"){               
                cellchat <- computeCommunProb(cellchat, type = unlist(mean_method)[1], trim= as.numeric(unlist(mean_method)[2]) ,population.size = pop_size)
              } else {
                cellchat <- computeCommunProb(cellchat, type = unlist(mean_method)[1], population.size = pop_size)
              }
          }

          #> The number of highly variable ligand-receptor pairs used for signaling inference is 692
          # CAUTION: The number of inferred ligand-receptor pairs clearly depends on the method for calculating the average gene expression per cell group. By default, CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. However, we find that CellChat performs well at predicting stronger interactions, which is very helpful for narrowing down on interactions for further experimental validations. In computeCommunProb, we provide an option for using other methods, such as 5% and 10% truncated mean, to calculating the average gene expression. Of note, ‘trimean’ approximates 25% truncated mean, implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%. To use 10% truncated mean, USER can set type = "truncatedMean" and trim = 0.1. To determine a proper value of trim, CellChat provides a function computeAveExpr, which can help to check the average expression of signaling genes of interest, e.g, computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1). Therefore, if well-known signaling pathways in the studied biological process are not predicted, users can try truncatedMean with lower values of trim to change the method for calculating the average gene expression per cell group.

          #When analyzing unsorted single-cell transcriptomes, under the assumption that abundant cell populations tend to send collectively stronger signals than the rare cell populations, CellChat can also consider the effect of cell proportion in each cell group in the probability calculation. USER can set population.size = TRUE. 
          ## This is the cas for brushings and biopsy samples

          ptm = Sys.time()
          # cellchat@data.signaling
          # mean_method <- list("truncatedMean", 0.1)
          # cellchat@data.signaling
      
          # cellchat@data.signaling


          # Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.

          cellchat <- filterCommunication(cellchat, min.cells = 10)

          # df.net <- subsetCommunication(cellchat) 

          cellchat <- computeCommunProbPathway(cellchat)
          slot.name = "netP"
          # df.net <- subsetCommunication(cellchat)
          # df.net <- subsetCommunication(cellchat, sources.use = c("c0","c1", "c2", "c3", "c4"), targets.use = c("Suprabasal", "Ciliated", "Basal", "Secretory", "Ionocytes"))

          cellchat <- aggregateNet(cellchat)
          execution.time = Sys.time() - ptm
          print(as.numeric(execution.time, units = "secs"))
          #> [1] 38.73308

          ptm = Sys.time()
          groupSize <- as.numeric(table(cellchat@idents))

          pdf(paste0(outdir, "pathways_inferred.pdf"))
          par(mfrow = c(1,2), xpd=TRUE)
          netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
          netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
          dev.off()

          pdf(paste0(outdir, "pathways_inferred_weight.pdf"), width=10, height=10)
          mat <- cellchat@net$weight
          par(mfrow = c(5,2), xpd=TRUE)
          for (i in 1:nrow(mat)) {
            mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
            mat2[i, ] <- mat[i, ]
            netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
          }

          dev.off()

          cat("Plotting... \n")
          ptm = Sys.time()
          # Compute the network centrality scores
          cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

          pdf(paste0(outdir, "Outgoing_and_incoming_signals.pdf"), width=length(unique(cellchat@netP$pathways))/6, height=length(unique(cellchat@netP$pathways))/4)
          # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
          ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = length(unique(cellchat@netP$pathways))/6, height = length(unique(cellchat@netP$pathways))/4)
          ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = length(unique(cellchat@netP$pathways))/6, height = length(unique(cellchat@netP$pathways))/4)
          print(ht1 + ht2)
          dev.off()

          # vertex.receiver = c("Suprabasal", "Ciliated", "Basal", "Secretory", "Ionocytes") # a numeric vector. 

          # Access all the signaling pathways showing significant communications
          pathways.show.all <- cellchat@netP$pathways
          
          saveRDS(pathways.show.all, paste0(outdir, "/", "Sig_pathways_", length(pathways.show.all), ".rds"))
      
          # check the order of cell identity to set suitable vertex.receiver
          dir.create(paste0(outdir, "/all_pathways/"))
          levels(cellchat@idents)
          # vertex.receiver = seq(1,4)
          try(for (i in 1:length(pathways.show.all)) {
            # i<-1
            setwd(paste0(outdir, "/all_pathways/"))
            # Visualize communication network associated with both signaling pathway and individual L-R pairs
            order_celltypes <- names(table(cellchat@meta$orig.celltype_cellchat)[table(cellchat@meta$orig.celltype_cellchat) > 10])

            order_celltypes <- grep(paste0(source_cellchat, collapse="|"), order_celltypes[order(order_celltypes)])

            netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = order_celltypes, layout = "hierarchy", out.format = c("svg"))
            # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
            gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
            ggsave(filename=paste0(outdir, "/all_pathways/", pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 7, height = 7, units = 'in', dpi = 300)
          })
          rm(brushings_cellchat_SA); cat("Finishing \n")
          saveRDS(cellchat, paste0(outdir, "/", disease_cellchat, "_cellchat_object.rds"))
          return(cellchat)
        }, mc.cores = 2, mc.preschedule = FALSE, mc.silent= FALSE)
        
        # Cellchat comparision of multiple datasets or conditions.
        #   https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

          #  cellchat_objects_orig <- cellchat_objects
        # outdir <- "/home/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/cellchat_pipeline2/w_immune_all_epi/fc0p0.05MEAN-truncatedMean_0.1-pop_size_TRUE/"
        # load(paste0(outdir, "/cellchat_object.list_brushings_Asthma_Immune.RData"))            
        names(cellchat_objects) <- comparing_groups
        # cellchat <- mergeCellChat(cellchat_objects, add.names = comparing_groups)
        
        cat("Saving Objects \n")
        save(cellchat_objects, file = paste0(outdir, "/cellchat_object_list.RData"))
        # outdir <- "/mnt/bioadhoc/Groups/vd-vijay/fcastaneda/fcastaneda-temp/rnaseq-sc-standar/biopBal_Sara_2024/brushings_all/results/figures/figures_mindset/cellchat_pipeline/w_immune_all_epi/fc0.25p0.05MEAN-triMean-pop_size_TRUE"
        write.csv("hola", paste0(outdir, "/hola.csv")) # Checking of I'm in the rigth path. Ignore

      } else { 
        cat("Loading previous run... \n")
        load(paste0(outdir, "/cellchat_object_list.RData"))
        cat("Load done \n ")
      }

        # cellchat_objects[["Severe.asthma"]]@pairLRsig
        cat("Names:", names(cellchat_objects), "\n")
        net_groups <- lapply(comparing_groups, function(x){
          # x <- comparing_groups[1]
            cat(x, "\n")
          get_pathways <- subsetCommunication(cellchat_objects[[x]], slot.name = "netP")
          net <- netMappingDEG(cellchat_objects[[x]], features.name = "features", variable.all = TRUE)
          write.csv(get_pathways, paste0(outdir, "/", x, "_significant_pathways.csv"))
          write.csv(net, paste0(outdir, "/", x, "_significant_genes.csv"))
          return(net)
        })

        # net_groups$Severe.asthma %>% filter(pathway_name == "IL4" & source ==  "CD8" & target == "endothelial")

        names(net_groups) <- comparing_groups
        
        cat("Comparision_analysis - nonDGEA \n")

        outdir2<-paste0(outdir, "/comparision_analysis/")     
        dir.create(outdir2)

        library(ComplexHeatmap)
        library("tidyr")

        cellchat_all <- mergeCellChat(cellchat_objects, add.names = names(cellchat_objects))
        save(cellchat_all, file = paste0(outdir, "/cellchat_merged.RData"))

        celltype_dgea_imm <- unique(cellchat_all@meta$orig.celltype_cellchat)
        senders_celltypes <- senders_celltypes[senders_celltypes %in% celltype_dgea_imm]
        receivers_celltypes <- receivers_celltypes[receivers_celltypes %in% celltype_dgea_imm]

        ### Ordering the cellchat object list in the given order
          # combs_num <- Reduce(c,lapply(2:length(cellchat_objects),
          #                 function(x) combn(1:length(cellchat_objects),x,simplify=FALSE)))
        groups_to_compare <- comparing_groups[comparing_groups %in% names(cellchat_objects)]
        if(length(groups_to_compare)<2) stop("Selected groups are not present or are less than 2")
        cellchat_objects <- cellchat_objects[groups_to_compare]

          ## Chords plots
          # groupSize <- as.numeric(table(cellchat_objects[[1]]@idents))
          # netVisual_circle(cellchat_objects[[1]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, color.use = colors_celltypes[levels(cellchat_objects[[1]]@idents)],  title.name = "Number of interactions")
          dir.create(paste0(outdir2, "/circle/"))

          lapply(c("count", "weight"), function(stats_size){
            # stats_size <- "count"
            weight.max <- getMaxWeight(cellchat_objects, attribute = c("idents",stats_size))
            
            group.size.max <- max(table(cellchat_all@meta$orig.celltype_cellchat, cellchat_all@meta[,column_group]))

            pdf(paste0(outdir2, "/circle/circle-ALL_", stats_size, ".pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)      
            for (i in 1:length(groups_to_compare)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], vertex.weight = groupSize, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0(stats_size, " of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max, color.use=colors_celltypes[order_colors]))
            }
            dev.off()

            # cellchat_all
            pdf(paste0(outdir2, "/circle/circle-ALL_", stats_size, "_blank.pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)
            for (i in 1:length(groups_to_compare)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], weight.scale = T, vertex.weight = groupSize, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0(stats_size," of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max, color.use=colors_celltypes[order_colors],  edge.label.color = "#00000000", vertex.label.color = "#00000000"))
            }
            dev.off()
          
            ### For senders

            weight.max2 <- max(unlist(lapply(names(cellchat_objects), function(x){
              # x<-names(cellchat_objects)[1]
              tomax <- reshape2::melt(cellchat_objects[[x]]@net[[stats_size]])
              colnames(tomax)[1:2] <- c("source","target")
              df.net <- filter(tomax, (source %in% senders_celltypes) & (target %in% receivers_celltypes))
              df.net$value[is.na(df.net$value)] <- 0
              return(max(df.net$value))
            })))
      
            pdf(paste0(outdir2, "/circle/circle-Senders_to_Receivers_", stats_size, ".pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)
            for (i in 1:length(groups_to_compare)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], weight.scale = T, vertex.weight = groupSize, label.edge= F, edge.weight.max = weight.max2, edge.width.max = 8, title.name = paste0("Number of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max, color.use=colors_celltypes[order_colors], sources.use = senders_celltypes, targets.use = receivers_celltypes))
              rm(groupSize)
            }
            dev.off()

            pdf(paste0(outdir2, "/circle/circle-Senders_to_Receivers_", stats_size, "_blank.pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)
            for (i in 1:length(groups_to_compare)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], weight.scale = T, vertex.weight = groupSize, label.edge= F, edge.weight.max = weight.max2, edge.width.max = 8, title.name = paste0("Number of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max, color.use=colors_celltypes[order_colors], sources.use = senders_celltypes, targets.use = receivers_celltypes, edge.label.color = "#00000000"))
            }
            dev.off()
          
            ### For others
        
            weight.max2 <- max(unlist(lapply(names(cellchat_objects), function(x){
              # x<-names(cellchat_objects)[1]
              tomax <- reshape2::melt(cellchat_objects[[x]]@net[[stats_size]])
              colnames(tomax)[1:2] <- c("source","target")
              df.net <- filter(tomax, (source %in% receivers_celltypes) & (target %in% senders_celltypes))
              df.net$value[is.na(df.net$value)] <- 0
              return(max(df.net$value))
            })))
      
            pdf(paste0(outdir2, "/circle/circle-Receiver_to_Sender_", stats_size, ".pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)
            for (i in 1:length(groups_to_compare)) {
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              order_colors <-levels(cellchat_objects[[i]]@idents)
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], weight.scale = T, vertex.weight = groupSize, label.edge= F, edge.weight.max = weight.max2, edge.width.max = 8, title.name = paste0("Number of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max, color.use=colors_celltypes[order_colors], sources.use = receivers_celltypes, targets.use = senders_celltypes))
              rm(groupSize)
            }
            dev.off()

            pdf(paste0(outdir2, "/circle/circle-Receiver_to_Sender_", stats_size, "_blank.pdf"), width=10, height=7)
            par(mfrow = c(1,1), xpd=TRUE)
            for (i in 1:length(groups_to_compare)) {
              groupSize <- as.numeric(table(cellchat_objects[[i]]@idents))
              order_colors <-levels(cellchat_objects[[i]]@idents)
              print(netVisual_circle(cellchat_objects[[i]]@net[[stats_size]], weight.scale = T, vertex.weight = groupSize, label.edge= F, edge.weight.max = weight.max2, edge.width.max = 8, title.name = paste0("Number of interactions - ", names(cellchat_objects)[i]), vertex.weight.max=group.size.max,  color.use=colors_celltypes[order_colors], sources.use = receivers_celltypes , targets.use = senders_celltypes, edge.label.color = "#00000000"))
            }
            dev.off()

          })

          ### Rank net plots ## This part takes almost 400 lines.. Please try to reduce it, at least at half.. I think is possible. 16-Sep-2025
          dir.create(paste0(outdir2, "/info_flow/"))
         
          message_plots<-paste0("Selected_means-", paste(celltype_selected, sep="&"))

          write.table(message_plots, paste0(outdir2,"/info_flow/",message_plots, ".txt"))

      ## Information flow
        ## All
          sig_pathways_all <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1))
              aver$dataset <- group2
              return(aver)
          })

          rank_all <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE)

          paths_select_all <- Reduce(rbind, sig_pathways_all) %>% pull(pathway_name) %>% unique()
          paths_select_all <- paths_select_all[paths_select_all %in% unique(rank_all$signaling.contribution$name)]

          rank_all <- rankNet(cellchat_all, signaling=paths_select_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE)

          pdf(paste0(outdir2, "/info_flow/information_flow_ALL.pdf"), width=10, height=ifelse(length(paths_select_all) < 100, 10, length(paths_select_all)/10))
          print(rank_all$gg.obj)
          dev.off()
          write.csv(rank_all$signaling.contribution, paste0(outdir2, "/info_flow/information_flow_ALL.csv"))

        ## Senders to Receivers
        rank_senders <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)

          sig_pathways_senders <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% senders_celltypes & target %in% receivers_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_senders <- Reduce(rbind, sig_pathways_senders) %>% pull(pathway_name) %>% unique()
          paths_select_senders <- paths_select_senders[paths_select_senders %in% unique(rank_senders$signaling.contribution$name)]

          rank_senders <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_senders, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)

          pdf(paste0(outdir2, "/info_flow/information_flow-Senders_to_Receivers.pdf"), width=10, height=ifelse(length(paths_select_senders) < 100, 10, length(paths_select_senders)/10))
          print(rank_senders$gg.obj)
          dev.off()
          write.csv(rank_senders$signaling.contribution, paste0(outdir2, "/info_flow/information_flow-Senders_to_Receivers.csv"))
       
        ## Receivers to Senders
          rank_others <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)

          sig_pathways_others <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% receivers_celltypes & target %in% senders_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_others <- Reduce(rbind, sig_pathways_others) %>% pull(pathway_name) %>% unique()
          paths_select_others <- paths_select_others[paths_select_others %in% unique(rank_others$signaling.contribution$name)]

          rank_others <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_others, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)

          pdf(paste0(outdir2, "/info_flow/information_flow-Receivers_to_Senders.pdf"), width=10, height=ifelse(length(paths_select_others) < 100, 10, length(paths_select_senders)/10))
          print(rank_others$gg.obj)
          dev.off()
          write.csv(rank_others$signaling.contribution, paste0(outdir2, "/info_flow/information_flow-Receivers_to_Senders.csv"))

        ## Selected to Receivers

          rank_cd8 <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)

          sig_pathways_cd8 <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% celltype_selected & target %in% receivers_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_cd8 <- Reduce(rbind, sig_pathways_cd8) %>% pull(pathway_name) %>% unique()
          paths_select_cd8 <- paths_select_cd8[paths_select_cd8 %in% unique(rank_cd8$signaling.contribution$name)]

          rank_cd8 <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_cd8, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)

          pdf(paste0(outdir2, "/info_flow/information_flow-Selected_to_Receivers.pdf"), width=10, height=ifelse(length(paths_select_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_cd8$gg.obj)
          dev.off()
          write.csv(rank_cd8$signaling.contribution, paste0(outdir2, "/info_flow/information_flow-Selected_to_Receivers.csv"))

          ## Selected to Receivers for each one
          dir.create(paste0(outdir2, "/info_flow/Selected_to_Receivers_individuals/"))
          lapply(receivers_celltypes, function(to_per_receiver){
              # to_per_receiver <- receivers_celltypes[2]
              cat("Info flow for Selected to ", to_per_receiver, "\n")
              rank_to_per_receiver <- "No_interactions"
              try(rank_to_per_receiver <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver))
              if(rank_to_per_receiver[1] == "No_interactions") return("No_interactions")

              sig_pathways_to_per_receiver <- lapply(comparing_groups, function(group2){
                  # group2 <- comparing_groups[3]
                  aver <- net_groups[[group2]]
                  aver$dataset <- group2
                  aver <- aver %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% celltype_selected & target %in% to_per_receiver)
                  return(aver)
              })

              paths_select_to_per_receiver <- Reduce(rbind, sig_pathways_to_per_receiver) %>% pull(pathway_name) %>% unique()
              # paths_select_to_per_receiver_d <- Reduce(rbind, sig_pathways_to_per_receiver)
              

              paths_select_to_per_receiver <- paths_select_to_per_receiver[paths_select_to_per_receiver %in% as.character(unique(rank_to_per_receiver$signaling.contribution$name))]

              # paths_select_to_per_receiver_d<-paths_select_to_per_receiver_d %>% filter(pathway_name %in% as.character(unique(rank_to_per_receiver$signaling.contribution$name)))

              if(length(paths_select_to_per_receiver) == 0 ) return("No_interactions")

              rank_to_per_receiver <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_to_per_receiver, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver)

              # rank_to_per_receiver_2 <- rankNet(cellchat_all, mode = "comparison", pairLR=unique(as.character(paths_select_to_per_receiver_d$interaction_name)), comparison=1:length(groups_to_compare), measure = "weight", stacked = T, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver, slot.name = 'net')
              

              pdf(paste0(outdir2, "/info_flow/Selected_to_Receivers_individuals/information_flow-Selected_to_", to_per_receiver, ".pdf"), width=10, height=ifelse(length(paths_select_to_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_to_per_receiver$gg.obj)
              # print(rank_to_per_receiver_2$gg.obj)
              dev.off()
              
              write.csv(rank_to_per_receiver$signaling.contribution, paste0(outdir2,"/info_flow/Selected_to_Receivers_individuals/information_flow-Selected_to_", to_per_receiver, ".csv"))
          })

        ## Receivers to Selected     
          rank_to_cd8 <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use =celltype_selected)

          sig_pathways_to_cd8 <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% receivers_celltypes & target %in% celltype_selected)
              return(aver)
          })

          paths_select_to_cd8 <- Reduce(rbind, sig_pathways_to_cd8) %>% pull(pathway_name) %>% unique()
          paths_select_to_cd8 <- paths_select_to_cd8[paths_select_to_cd8 %in% unique(rank_to_cd8$signaling.contribution$name)]

          rank_to_cd8 <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_to_cd8, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes , targets.use = celltype_selected)

          pdf(paste0(outdir2, "/info_flow/information_flow-Receivers_to_Selected.pdf"), width=10, height=ifelse(length(paths_select_to_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_to_cd8$gg.obj)
          dev.off()
          write.csv(rank_to_cd8$signaling.contribution, paste0(outdir2, "/info_flow/information_flow-Receivers_to_Selected.csv"))

          ## Receivers to Selected for each one
          dir.create(paste0(outdir2, "/info_flow/Receivers_to_Selected_individuals/"))
          lapply(receivers_celltypes, function(from_per_receiver){
              # from_per_receiver <- receivers_celltypes[2]
              cat("Info flow for ", from_per_receiver, "to Selected\n")
              rank_from_per_receiver <- "No_interactions"
              try(rank_from_per_receiver <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver, targets.use = celltype_selected))

              if(rank_from_per_receiver[1] == "No_interactions") return("No_interactions")

              sig_pathways_from_per_receiver <- lapply(comparing_groups, function(group2){
                  # group2 <- comparing_groups[3]
                  aver <- net_groups[[group2]]
                  aver$dataset <- group2
                  aver <- aver %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% from_per_receiver & target %in% celltype_selected )
                  return(aver)
              })

              paths_select_from_per_receiver <- Reduce(rbind, sig_pathways_from_per_receiver) %>% pull(pathway_name) %>% unique()

              paths_select_from_per_receiver <- paths_select_from_per_receiver[paths_select_from_per_receiver %in% as.character(unique(rank_from_per_receiver$signaling.contribution$name))]

              if(length(paths_select_from_per_receiver) == 0 ) return("No_interactions")

              rank_from_per_receiver <- rankNet(cellchat_all, mode = "comparison", signaling=paths_select_from_per_receiver, comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver , targets.use = celltype_selected)

              pdf(paste0(outdir2, "/info_flow/Receivers_to_Selected_individuals/information_flow-", from_per_receiver, "_to_Selected", ".pdf"), width=10, height=ifelse(length(paths_select_from_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_from_per_receiver$gg.obj)
              dev.off()
              
              write.csv(rank_from_per_receiver$signaling.contribution, paste0(outdir2,"/info_flow/Receivers_to_Selected_individuals/information_flow-", from_per_receiver, "_to_Selected.csv"))
          })

          #### Dot plots for immune to epi?? taking only cd8 significant pathways

          dir.create(paste0(outdir2, "/info_flow/dotplot_information_flow/"))
          ## Selecting only significant pathways
          to_dotplot <- Reduce(rbind, sig_pathways_all)
          to_dotplot <- to_dotplot %>% filter(pathway_name %in% paths_select_cd8 & source %in% senders_celltypes & target %in% receivers_celltypes)

          ## Selecting only significant pair of ligands - receptors
          to_dotplot_cd8 <- Reduce(rbind, sig_pathways_cd8)
          sig_pairs <- as.character(unique(to_dotplot_cd8$interaction_name))
          to_dotplot <- to_dotplot %>% filter(interaction_name %in% sig_pairs)
          table(to_dotplot$source)
          table(to_dotplot$target)

          to_dotplot$source_to_group <- ifelse(to_dotplot$source == celltype_selected, celltype_selected, "O_immune")
          to_dotplot$target_to_group <- "Epi"
          table(to_dotplot$source, to_dotplot$target_to_group)
          table(to_dotplot$target, to_dotplot$target_to_group)

          write.csv(to_dotplot, paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Selected_to_CombinedReceivers.csv"))

          to_dotplot_min<- to_dotplot %>% group_by(interaction_name, source_to_group, target_to_group, dataset) %>% summarize(sum_prob=sum(prob), mean_prob=mean(prob), p_mean=-log10(mean(pval)+1), pathway=unique(pathway_name), pathway2=unique(interaction_name_2))

          to_dotplot_min$interaction <- paste0(to_dotplot_min$source_to_group, " -> ", to_dotplot_min$target_to_group, ".(", to_dotplot_min$dataset, ")")

          to_dotplot_min <- to_dotplot_min %>% arrange(factor(pathway, levels=paths_select_cd8))

            hjust.x = 1
            hjust.y = 0.5
            vjust.x = 0.3
            angle.x = 90
            order_xaxis <- unique(to_dotplot_min$interaction)[unique(to_dotplot_min$interaction) %>% order()]

          pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Selected_to_CombinedReceivers.pdf"), width=10, height=20)
          print(ggplot(to_dotplot_min, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_min$pathway2))), color = mean_prob)) +
                geom_point(pch = 16, size=4) +
                # geom_point(pch = 1, color="gray90") +
                theme_linedraw() + theme(panel.grid.major = element_blank()) +
                theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) +
                scale_x_discrete(position = "bottom") +  
                scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_min$mean_prob, 0,na.rm= T), quantile(to_dotplot_min$mean_prob, 1,na.rm= T)), breaks = c(quantile(to_dotplot_min$mean_prob, 0,na.rm= T), quantile(to_dotplot_min$mean_prob, 1,na.rm= T)), labels = c("min","max")) +
                guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) ) # + labs(title=title_plot)
          dev.off()

       to_dotplot_maj<- to_dotplot %>% group_by(interaction_name, source, dataset) %>% summarize(
        source=unique(source), 
        target=unique(target_to_group),
         sum_prob=sum(prob), 
         mean_prob=mean(prob), 
         p_mean=-log10(mean(pval)+1), 
         pathway=unique(pathway_name), 
         pathway2=unique(interaction_name_2), 
         ligand.pct.1=unique(ligand.pct.1), 
         ligand.pct.2=unique(ligand.pct.2), 
         ligand.pvalues=mean(ligand.pvalues), 
         ligand.logFC=mean(ligand.logFC), 
         receptor.pct.1=mean(receptor.pct.1), 
         receptor.pct.2=mean(receptor.pct.2), 
         receptor.pvalues=mean(receptor.pvalues), 
         receptor.logFC=mean(receptor.logFC))

       to_dotplot_maj$interaction <- paste0(to_dotplot_maj$source, " -> ", to_dotplot_maj$target, ".(", to_dotplot_maj$dataset, ")")

      to_dotplot_maj <- to_dotplot_maj %>% arrange(factor(pathway, levels=paths_select_cd8))
    
      to_dotplot_maj <- to_dotplot_maj %>% group_by(pathway2) %>%
        mutate(x_centered = scale(mean_prob))
      
      order_xaxis <- unique(to_dotplot_maj$interaction)[unique(to_dotplot_maj$interaction) %>% order()]

        pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Senders_to_CombinedReceivers_percExpr.pdf"), width=15, height=15)

          print(ggplot(to_dotplot_maj, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj$pathway2))), color = x_centered, fill=ligand.pct.1, size = receptor.pct.2)) +
                geom_point(shape = 21, stroke = 2) +
                scale_size_area(max_size = 10) +
                # geom_point(pch = 1, color="gray90") +
                theme_linedraw() + theme(panel.grid.major = element_blank()) +
                theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) +
                scale_x_discrete(position = "bottom") +  
                scale_colour_gradientn(colors = c("#0a0a0a", "white", "#f24444", "#f50c0c", "#f50202"), na.value = "gray90", limits=c(quantile(to_dotplot_maj$x_centered, 0,na.rm= T), quantile(to_dotplot_maj$x_centered, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj$x_centered, 0,na.rm= T), quantile(to_dotplot_maj$x_centered, 1,na.rm= T)), labels = c("min","max")) + scale_fill_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#f29102"), na.value = "gray90") +
                guides(color = guide_colourbar(barwidth = 0.5, title = "z-score Commun. Prob.")))
          dev.off()
          
          pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Senders_to_CombinedReceivers.pdf"), width=15, height=10)

          print(ggplot(to_dotplot_maj, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj$pathway2))), color = mean_prob)) +
                geom_point(pch = 16, size=4) +
                # geom_point(pch = 1, color="gray90") +
                theme_linedraw() + theme(panel.grid.major = element_blank()) +
                theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) +
                scale_x_discrete(position = "bottom") +  
                scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), labels = c("min","max")) +
                guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) ) # + labs(title=title_plot)
          dev.off()

          write.csv(to_dotplot_maj, paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Senders_to_CombinedReceivers_percExpr.csv"))
          
          #### Facet per target

       to_dotplot_maj_target<- to_dotplot %>% group_by(interaction_name, source, dataset) %>% summarize(
        source=unique(source), 
        target=unique(target),
         sum_prob=sum(prob), 
         mean_prob=mean(prob), 
         p_mean=-log10(mean(pval)+1), 
         pathway=unique(pathway_name), 
         pathway2=unique(interaction_name_2), 
         ligand.pct.1=unique(ligand.pct.1), 
         ligand.pct.2=unique(ligand.pct.2), 
         ligand.pvalues=mean(ligand.pvalues), 
         ligand.logFC=mean(ligand.logFC), 
         receptor.pct.1=mean(receptor.pct.1), 
         receptor.pct.2=mean(receptor.pct.2), 
         receptor.pvalues=mean(receptor.pvalues), 
         receptor.logFC=mean(receptor.logFC))

      to_dotplot_maj_target$interaction <- paste0(to_dotplot_maj_target$source, " -> ",".(", to_dotplot_maj_target$dataset, ")")

      to_dotplot_maj_target <- to_dotplot_maj_target %>% arrange(factor(pathway, levels=paths_select_cd8))
    
      to_dotplot_maj_target <- to_dotplot_maj_target %>% group_by(pathway2, target) %>%
        mutate(x_centered = scale(mean_prob))
      
      order_xaxis <- unique(to_dotplot_maj_target$interaction)[unique(to_dotplot_maj_target$interaction) %>% order()]


       dotplot_maj <- ggplot(to_dotplot_maj_target, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj_target$pathway2))), color = x_centered, fill=ligand.logFC, size = receptor.logFC)) +
                geom_point(shape = 21, stroke = 2) +
                scale_size_area(max_size = 10) +
                # geom_point(pch = 1, color="gray90") +
                theme_linedraw() + theme(panel.grid.major = element_blank()) +
                theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) +
                scale_x_discrete(position = "bottom") +  
                scale_colour_gradientn(colors = c("#0a0a0a", "white", "#f24444", "#f50c0c", "#f50202"), na.value = "gray90", limits=c(quantile(to_dotplot_maj_target$x_centered, 0,na.rm= T), quantile(to_dotplot_maj_target$x_centered, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj_target$x_centered, 0,na.rm= T), quantile(to_dotplot_maj_target$x_centered, 1,na.rm= T)), labels = c("min","max")) + scale_fill_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#f29102"), na.value = "gray90") + 
                facet_wrap_paginate(~ target, ncol = 1, nrow = 1) + 
                guides(color = guide_colourbar(barwidth = 0.5, title = "z-score  Commun. Prob."))

          pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/pathways_from_Selected-displayed_in_Senders_to_Receivers_percExpr.pdf"), width=12, height=12)

          for(i in 1:n_pages(dotplot_maj)){
            p_save <-  dotplot_maj + 
              facet_wrap_paginate(~ target, ncol = 1, nrow = 1, page = i)
            print(p_save)
          }

          dev.off()
          
          ## What is this, I don't know.. I think is the same as dotplot_information_flow/pathways_from_Selected-displayed_in_Senders_to_CombinedReceivers.pdf
          # pdf(paste0(outdir2, "/info_flow/information_flow--selected_immune_comb_dotplot_major.pdf"), width=15, height=10)

          # print(ggplot(to_dotplot_maj, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj$pathway2))), color = mean_prob)) +
          #       geom_point(pch = 16, size=4) +
          #       # geom_point(pch = 1, color="gray90") +
          #       theme_linedraw() + theme(panel.grid.major = element_blank()) +
          #       theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
          #             axis.title.x = element_blank(),
          #             axis.title.y = element_blank()) +
          #       scale_x_discrete(position = "bottom") +  
          #       scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), labels = c("min","max")) +
          #       guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) ) # + labs(title=title_plot)
          # dev.off()

          # write.csv(to_dotplot_maj, paste0(outdir2, "/info_flow/dotplot_information_flow--selected_EPIcombined.csv"))


         #### Dot plots for epi to immune?? taking only cd8 significant pathways

          ## Selecting only significant pathways
          to_dotplot <- Reduce(rbind, sig_pathways_all)
          to_dotplot <- to_dotplot %>% filter(pathway_name %in% paths_select_to_cd8 & source %in% receivers_celltypes & target %in% senders_celltypes)

          ## Selecting only significant pair of ligands - receptors
          to_dotplot_cd8 <- Reduce(rbind, sig_pathways_to_cd8)
          sig_pairs <- as.character(unique(to_dotplot_cd8$interaction_name))
          to_dotplot <- to_dotplot %>% filter(interaction_name %in% sig_pairs)
          table(to_dotplot$source)
          table(to_dotplot$target)

          to_dotplot$source_to_group <- "Epi"
          to_dotplot$target_to_group <- ifelse(to_dotplot$target == celltype_selected, celltype_selected, "O_immune")
          table(to_dotplot$source, to_dotplot$target_to_group)
          table(to_dotplot$target, to_dotplot$target_to_group)

          write.csv(to_dotplot, paste0(outdir2, "/info_flow/dotplot_information_flow/Receivers_to_Selected_pathways-displayed_in_CombinedReceivers_to_Selected.csv"))

          to_dotplot_min<- to_dotplot %>% group_by(interaction_name, source_to_group, target_to_group, dataset) %>% summarize(sum_prob=sum(prob), mean_prob=mean(prob), p_mean=-log10(mean(pval)+1), pathway=unique(pathway_name), pathway2=unique(interaction_name_2))

          to_dotplot_min$interaction <- paste0(to_dotplot_min$source_to_group, " -> ", to_dotplot_min$target_to_group, ".(", to_dotplot_min$dataset, ")")
          table(to_dotplot_min$interaction)

          to_dotplot_min <- to_dotplot_min %>% arrange(factor(pathway, levels=sig_pathways_to_cd8))

            hjust.x = 1
            hjust.y = 0.5
            vjust.x = 0.3
            angle.x = 90
            order_xaxis <- unique(to_dotplot_min$interaction)[unique(to_dotplot_min$interaction) %>% order()]

          pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/Receivers_to_Selected_pathways-displayed_in_CombinedReceivers_to_Selected.pdf"), width=10, height=20)
          print(ggplot(to_dotplot_min, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_min$pathway2))), color = mean_prob)) +
                geom_point(pch = 16, size=4) +
                # geom_point(pch = 1, color="gray90") +
                theme_linedraw() + theme(panel.grid.major = element_blank()) +
                theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank()) +
                scale_x_discrete(position = "bottom") +  
                scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_min$mean_prob, 0,na.rm= T), quantile(to_dotplot_min$mean_prob, 1,na.rm= T)), breaks = c(quantile(to_dotplot_min$mean_prob, 0,na.rm= T), quantile(to_dotplot_min$mean_prob, 1,na.rm= T)), labels = c("min","max")) +
                guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) ) # + labs(title=title_plot)
          dev.off()

        to_dotplot_maj<- to_dotplot %>% group_by(interaction_name, source, dataset) %>% summarize(
            source=unique(source_to_group), 
            target=unique(target),
            sum_prob=sum(prob), 
            mean_prob=mean(prob), 
            p_mean=-log10(mean(pval)+1), 
            pathway=unique(pathway_name), 
            pathway2=unique(interaction_name_2), 
            ligand.pct.1=unique(ligand.pct.1), 
            ligand.pct.2=unique(ligand.pct.2), 
            ligand.pvalues=mean(ligand.pvalues), 
            ligand.logFC=mean(ligand.logFC), 
            receptor.pct.1=mean(receptor.pct.1), 
            receptor.pct.2=mean(receptor.pct.2), 
            receptor.pvalues=mean(receptor.pvalues), 
            receptor.logFC=mean(receptor.logFC))

        to_dotplot_maj$interaction <- paste0(to_dotplot_maj$source, " -> ", to_dotplot_maj$target, ".(", to_dotplot_maj$dataset, ")")

        to_dotplot_maj <- to_dotplot_maj %>% arrange(factor(pathway, levels=sig_pathways_to_cd8))
        
        to_dotplot_maj <- to_dotplot_maj %>% group_by(pathway2) %>%
            mutate(x_centered = scale(mean_prob))
        
        order_xaxis <- unique(to_dotplot_maj$interaction)[unique(to_dotplot_maj$interaction) %>% order()]

            pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/Receivers_to_Selected_pathways-displayed_in_CombinedReceivers_to_Senders.pdf"), width=15, height=20)

            print(ggplot(to_dotplot_maj, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj$pathway2))), color = x_centered)) +
                    geom_point(pch = 16, size=4) +
                    # geom_point(pch = 1, color="gray90") +
                    theme_linedraw() + theme(panel.grid.major = element_blank()) +
                    theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
                        axis.title.x = element_blank(),
                        axis.title.y = element_blank()) +
                    scale_x_discrete(position = "bottom") +  
                    scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_maj$x_centered, 0,na.rm= T), quantile(to_dotplot_maj$x_centered, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj$x_centered, 0,na.rm= T), quantile(to_dotplot_maj$x_centered, 1,na.rm= T)), labels = c("min","max")) +
                    guides(color = guide_colourbar(barwidth = 0.5, title = "z-score Commun. Prob."))) # + labs(title=title_plot)
            dev.off()
            
            # pdf(paste0(outdir2, "/info_flow/dotplot_information_flow/Receivers_to_Selected_pathways-displayed_in_CombinedReceivers_to_Senders2.pdf"), width=15, height=10)

            # print(ggplot(to_dotplot_maj, aes(x = factor(interaction, levels=order_xaxis), y = factor(pathway2, levels=rev(unique(to_dotplot_maj$pathway2))), color = mean_prob)) +
            #         geom_point(pch = 16, size=4) +
            #         # geom_point(pch = 1, color="gray90") +
            #         theme_linedraw() + theme(panel.grid.major = element_blank()) +
            #         theme(axis.text.x = element_text(angle = angle.x, hjust= hjust.x, vjust = vjust.x),
            #             axis.title.x = element_blank(),
            #             axis.title.y = element_blank()) +
            #         scale_x_discrete(position = "bottom") +  
            #         scale_colour_gradientn(colors = c("#502188", "white", "#FCD8B0", "#f29102", "#b56d02"), na.value = "gray90", limits=c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), breaks = c(quantile(to_dotplot_maj$mean_prob, 0,na.rm= T), quantile(to_dotplot_maj$mean_prob, 1,na.rm= T)), labels = c("min","max")) +
            #         guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob.")) ) # + labs(title=title_plot)
            # dev.off()

            write.csv(to_dotplot_maj, paste0(outdir2, "/info_flow/Receivers_to_Selected_pathways-displayed_in_CombinedReceivers_to_Senders.csv"))


      #### ------- revision Donaldo ------- #### ## CHEEECK!!
      dir.create(paste0(outdir2, "/info_flow/info_flow_LR/"))

            ## Information flow LR
        ## All LR
            # sig_pathways_all <- lapply(comparing_groups, function(group2){
            #       aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1))
            #       aver$dataset <- group2
            #       return(aver)
            #   })

        rank_all <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE)

        paths_select_all <- Reduce(rbind, sig_pathways_all) %>% pull(interaction_name, pathway_name)
        paths_select_all <- paths_select_all[names(paths_select_all) %in% unique(rank_all$signaling.contribution$name)]
        names(paths_select_all) <- NULL
        paths_select_all <- as.character(unique(paths_select_all))
        # rank_all_aver <- rank_all
        
        #stacked plot
        rank_all_aver <- rankNet(cellchat_all, pairLR=paths_select_all, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE)
        rank_all_aver_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_all, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE)

        pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow_ALL_LR.pdf"), width=10, height=ifelse(length(paths_select_all) < 100, 10, length(paths_select_all)/10))
        print(rank_all_aver$gg.obj)
        print(rank_all_aver_mod$gg.obj)
        dev.off()

        write.csv(rank_all_aver$signaling.contribution, paste0(outdir2, "/info_flow/info_flow_LR/information_flow_ALL_LR.csv"))
        
        #not stacked plot
        rank_all_aver_scaled <- rankNet(cellchat_all, pairLR=paths_select_all, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE)
        

        pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow_ALL_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_all) < 100, 10, length(paths_select_all)/10))
        print(rank_all_aver_scaled$gg.obj)
        dev.off()        
          
    ## Senders to Receivers LR
        rank_senders <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)

          sig_pathways_senders <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% senders_celltypes & target %in% receivers_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_senders <- Reduce(rbind, sig_pathways_senders) %>% pull(interaction_name, pathway_name)
          paths_select_senders <- paths_select_senders[names(paths_select_senders) %in% unique(rank_senders$signaling.contribution$name)]
          names(paths_select_senders) <- NULL
          paths_select_senders <- as.character(unique(paths_select_senders))

          #stacked plot
          rank_senders <- rankNet(cellchat_all, pairLR=paths_select_senders, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)
          rank_senders_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_senders, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Senders_to_Receivers_LR.pdf"), width=10, height=ifelse(length(paths_select_senders) < 100, 10, length(paths_select_senders)/10))
          print(rank_senders$gg.obj)
          print(rank_senders_mod$gg.obj)
          dev.off()

          write.csv(rank_senders$signaling.contribution, paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Senders_to_Receivers_LR.csv"))

          #not stacked plot
          rank_senders_scaled <- rankNet(cellchat_all, pairLR=paths_select_senders, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = senders_celltypes, targets.use = receivers_celltypes)
          
          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Senders_to_Receivers_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_senders) < 100, 10, length(paths_select_senders)/10))
          print(rank_senders_scaled$gg.obj)
          dev.off()
        
    ## Receivers to Senders LR
          rank_others <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)

          sig_pathways_others <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% receivers_celltypes & target %in% senders_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_others <- Reduce(rbind, sig_pathways_others) %>% pull(interaction_name, pathway_name)
          paths_select_others <- paths_select_others[names(paths_select_others) %in% unique(rank_others$signaling.contribution$name)]
          names(paths_select_others) <- NULL
          paths_select_others <- as.character(unique(paths_select_others))

          #stacked plot
          rank_others <- rankNet(cellchat_all, pairLR=paths_select_others, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)
          rank_others_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_others, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Senders_LR.pdf"), width=10, height=ifelse(length(paths_select_others) < 100, 10, length(paths_select_others)/10))
          print(rank_others$gg.obj)
          print(rank_others_mod$gg.obj)
          dev.off()

          write.csv(rank_others$signaling.contribution, paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Senders_LR.csv"))

          #not stacked plot
          rank_others <- rankNet(cellchat_all, pairLR=paths_select_others, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use = senders_celltypes)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Senders_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_others) < 100, 10, length(paths_select_others)/10))
          print(rank_others$gg.obj)
          dev.off()
        
        ## Selected to receivers LR
        dir.create(paste0(outdir2, "/info_flow/info_flow_LR/Selected_to_Receivers_individuals/"))

          rank_cd8 <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)

          sig_pathways_cd8 <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% celltype_selected & target %in% receivers_celltypes)
              aver$dataset <- group2
              return(aver)
          })

          paths_select_cd8 <- Reduce(rbind, sig_pathways_cd8) %>% pull(interaction_name, pathway_name)
          paths_select_cd8 <- paths_select_cd8[names(paths_select_cd8) %in% unique(rank_cd8$signaling.contribution$name)]
          names(paths_select_cd8) <- NULL
          paths_select_cd8 <- as.character(unique(paths_select_cd8))

          #stacked plot
          rank_cd8 <- rankNet(cellchat_all, pairLR=paths_select_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)
          rank_cd8_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)
          
          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Selected_to_Receivers_LR.pdf"), width=10, height=ifelse(length(paths_select_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_cd8$gg.obj)
          print(rank_cd8_mod$gg.obj)
          dev.off()

          write.csv(rank_cd8$signaling.contribution, paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Selected_to_Receivers_LR.csv"))

          #not stacked plot
          rank_cd8 <- rankNet(cellchat_all, pairLR=paths_select_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = receivers_celltypes)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Selected_to_Receivers_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_cd8$gg.obj)
          dev.off()
          
       ## Selected to Receivers for each one L-R
          dir.create(paste0(outdir2, "/info_flow/info_flow_LR/Selected_to_Receivers_individuals"))
          lapply(receivers_celltypes, function(to_per_receiver){
              # to_per_receiver <- receivers_celltypes[2]
              cat("Info flow for Selected to ", to_per_receiver, "\n")
              rank_to_per_receiver <- "No_interactions"
              try(rank_to_per_receiver <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver))
              if(rank_to_per_receiver[1] == "No_interactions") return("No_interactions")

              sig_pathways_to_per_receiver <- lapply(comparing_groups, function(group2){
                  # group2 <- comparing_groups[3]
                  aver <- net_groups[[group2]]
                  aver$dataset <- group2
                  aver <- aver %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% celltype_selected & target %in% to_per_receiver)
                  return(aver)
              })

              # paths_select_to_per_receiver <- Reduce(rbind, sig_pathways_to_per_receiver) %>% pull(pathway_name) %>% unique()
              # paths_select_to_per_receiver <- paths_select_to_per_receiver[paths_select_to_per_receiver %in% as.character(unique(rank_to_per_receiver$signaling.contribution$name))]

              paths_select_to_per_receiver <- Reduce(rbind, sig_pathways_to_per_receiver) %>% pull(interaction_name, pathway_name)
              paths_select_to_per_receiver <- paths_select_to_per_receiver[names(paths_select_to_per_receiver) %in% unique(rank_to_per_receiver$signaling.contribution$name)]
              names(paths_select_to_per_receiver) <- NULL
              paths_select_to_per_receiver <- as.character(unique(paths_select_to_per_receiver))

              # paths_select_to_per_receiver_d<-paths_select_to_per_receiver_d %>% filter(pathway_name %in% as.character(unique(rank_to_per_receiver$signaling.contribution$name)))

              if(length(paths_select_to_per_receiver) == 0 ) return("No_interactions")

              #stacked plot
              rank_to_per_receiver <- rankNet(cellchat_all, pairLR=paths_select_to_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver)
              rank_to_per_receiver_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_to_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver)

              pdf(paste0(outdir2, "/info_flow/info_flow_LR/Selected_to_Receivers_individuals/information_flow-Selected_to_", to_per_receiver, "_LR.pdf"), width=10, height=ifelse(length(paths_select_to_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_to_per_receiver$gg.obj)
              print(rank_to_per_receiver_mod$gg.obj)
              # print(rank_to_per_receiver_2$gg.obj)
              dev.off()
              
              write.csv(rank_to_per_receiver$signaling.contribution, paste0(outdir2,"/info_flow/info_flow_LR/Selected_to_Receivers_individuals/information_flow-Selected_to_", to_per_receiver, "_LR.csv"))

              #not stacked plot
              rank_to_per_receiver <- rankNet(cellchat_all, pairLR=paths_select_to_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = celltype_selected, targets.use = to_per_receiver)

              pdf(paste0(outdir2, "/info_flow/info_flow_LR/Selected_to_Receivers_individuals/information_flow-Selected_to_", to_per_receiver, "_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_to_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_to_per_receiver$gg.obj)
              # print(rank_to_per_receiver_2$gg.obj)
              dev.off()

          })
          
      ## Receivers to Selected LR
          rank_to_cd8 <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes, targets.use =celltype_selected)

          sig_pathways_to_cd8 <- lapply(comparing_groups, function(group2){
              aver <- net_groups[[group2]] %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% receivers_celltypes & target %in% celltype_selected)
              return(aver)
          })

          paths_select_to_cd8 <- Reduce(rbind, sig_pathways_to_cd8) %>% pull(interaction_name, pathway_name)
          paths_select_to_cd8 <- paths_select_to_cd8[names(paths_select_to_cd8) %in% unique(rank_to_cd8$signaling.contribution$name)]
          names(paths_select_to_cd8) <- NULL
          paths_select_to_cd8 <- as.character(unique(paths_select_to_cd8))
          
          #stacked plot
          rank_to_cd8 <- rankNet(cellchat_all, pairLR=paths_select_to_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes , targets.use = celltype_selected)
          rank_to_cd8_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_to_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes , targets.use = celltype_selected)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Selected_LR.pdf"), width=10, height=ifelse(length(paths_select_to_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_to_cd8$gg.obj)
          print(rank_to_cd8_mod$gg.obj)
          dev.off()

          write.csv(rank_to_cd8$signaling.contribution, paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Selected_LR.csv"))
  
          #not stacked plot
          rank_to_cd8 <- rankNet(cellchat_all, pairLR=paths_select_to_cd8, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = receivers_celltypes , targets.use = celltype_selected)

          pdf(paste0(outdir2, "/info_flow/info_flow_LR/information_flow-Receivers_to_Selected_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_to_cd8) < 100, 10, length(paths_select_senders)/10))
          print(rank_to_cd8$gg.obj)
          dev.off()
  
      ## Receivers to Selected for each one LR
          dir.create(paste0(outdir2, "/info_flow/info_flow_LR/Receivers_to_Selected_individuals"))
          lapply(receivers_celltypes, function(from_per_receiver){
              # from_per_receiver <- receivers_celltypes[2]
              cat("Info flow for ", from_per_receiver, "to Selected\n")
              rank_from_per_receiver <- "No_interactions"
              try(rank_from_per_receiver <- rankNet(cellchat_all, mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver, targets.use = celltype_selected))

              if(rank_from_per_receiver[1] == "No_interactions") return("No_interactions")

              sig_pathways_from_per_receiver <- lapply(comparing_groups, function(group2){
                  # group2 <- comparing_groups[3]
                  aver <- net_groups[[group2]]
                  aver$dataset <- group2
                  aver <- aver %>% filter(!is.na(ligand.pct.1) & !is.na(receptor.pct.1) & source %in% from_per_receiver & target %in% celltype_selected )
                  return(aver)
              })

              paths_select_from_per_receiver <- Reduce(rbind, sig_pathways_from_per_receiver) %>% pull(interaction_name, pathway_name)
              paths_select_from_per_receiver <- paths_select_from_per_receiver[names(paths_select_from_per_receiver) %in% unique(rank_from_per_receiver$signaling.contribution$name)]
              names(paths_select_from_per_receiver) <- NULL
              paths_select_from_per_receiver <- as.character(unique(paths_select_from_per_receiver))

              if(length(paths_select_from_per_receiver) == 0 ) return("No_interactions")
              
              # Stacked plot
              rank_from_per_receiver <- rankNet(cellchat_all, pairLR=paths_select_from_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver , targets.use = celltype_selected)
              rank_from_per_receiver_mod <- rankNet_mod(cellchat_all, pairLR=paths_select_from_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = T, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver , targets.use = celltype_selected)

              pdf(paste0(outdir2, "/info_flow/info_flow_LR/Receivers_to_Selected_individuals/information_flow-", from_per_receiver, "_to_Selected", "_LR.pdf"), width=10, height=ifelse(length(paths_select_from_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_from_per_receiver$gg.obj)
              print(rank_from_per_receiver_mod$gg.obj)
              dev.off()

              write.csv(rank_from_per_receiver$signaling.contribution, paste0(outdir2,"/info_flow/info_flow_LR/Receivers_to_Selected_individuals/information_flow-", from_per_receiver, "_to_Selected_LR.csv"))
              
              # Not stacked plot
              rank_from_per_receiver <- rankNet(cellchat_all, pairLR=paths_select_from_per_receiver, slot.name = "net", mode = "comparison", comparison=1:length(groups_to_compare), measure = "weight", stacked = F, do.stat = TRUE, return.data = TRUE, sources.use = from_per_receiver , targets.use = celltype_selected)

              pdf(paste0(outdir2, "/info_flow/info_flow_LR/Receivers_to_Selected_individuals/information_flow-", from_per_receiver, "_to_Selected", "_scaled_LR.pdf"), width=10, height=ifelse(length(paths_select_from_per_receiver) < 100, 10, length(paths_select_senders)/10))
              print(rank_from_per_receiver$gg.obj)
              dev.off()

          })

    ### End Donaldo's revision

          #### Barplots counting # of interactions
          pdf(paste0(outdir2, "/barplots_count_interactions.pdf"), width=10, height=10)
          gg1 <- compareInteractions(cellchat_all, show.legend = F, group = 1:length(groups_to_compare))
          gg2 <- compareInteractions(cellchat_all, show.legend = F, group = 1:length(groups_to_compare), measure = "weight")
          print(gg1 + gg2)
          dev.off()

         ## Scatter for all pathways

          dir.create(paste0(outdir2, "/scatter/"))

          num.link <- sapply(cellchat_objects, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})

          weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
          gg <- list()
          x_axis_scatter <- c()
          y_axis_scatter <- c()
          
          scatter_plots_d<-data.frame()
          for (i in 1:length(cellchat_objects)) {
            order_colors <-levels(cellchat_objects[[i]]@idents)
            scatter_plots <- netAnalysis_signalingRole_scatter(cellchat_objects[[i]], title = names(cellchat_objects)[i], weight.MinMax = weight.MinMax, color.use=colors_celltypes[order_colors])
            x_axis_scatter <- c(x_axis_scatter, scatter_plots$data$x)
            y_axis_scatter <- c(y_axis_scatter, scatter_plots$data$y)
            gg[[i]] <- scatter_plots
            scatter_plots_disease <- as.data.frame(scatter_plots$data)
            scatter_plots_disease$set <- unique(cellchat_objects[[i]]@meta[,column_group])
            scatter_plots_d <- rbind(scatter_plots_d,scatter_plots_disease)
          }

          for(i in 1:length(cellchat_objects)){
            gg[[i]] <- gg[[i]] + xlim(0, max(x_axis_scatter)) +  ylim(0, max(y_axis_scatter))
          } 

          pdf(paste0(outdir2, "/", "/scatter/Scatter_all_pathways.pdf"), width=15, height=3.5)
          print(patchwork::wrap_plots(plots = gg ))
          dev.off()

          write.csv(scatter_plots_d, paste0(outdir2, "/", "/scatter/Scatter_all_pathways.csv"))

          scatter_plots_d$condition <- ifelse(scatter_plots_d$labels %in% senders_celltypes, "Senders", "Receivers")

           scatter_plots_d$condition2 <- paste0(scatter_plots_d$condition, "_", scatter_plots_d$set)
            
          #scatter approach
          pdf(paste0(outdir2,"/scatter/Scatter_Senders_Receivers.pdf"), width=8, height=3)
            scatter_plots_d_r <- scatter_plots_d %>% filter(scatter_plots_d$condition == "Receivers")

            sca_r <- ggplot(scatter_plots_d_r, aes(x = x, y = y)) + 
            geom_point(aes(size = Count, fill = labels, color = labels), alpha = .8, shape = 21, stroke = .6) +
            scale_fill_manual(values = colors_celltypes) +
            scale_color_manual(values = colors_celltypes) +
            scale_size_continuous(
              range = c(1, 12),
              breaks = pretty(scatter_plots_d_r$Count, n = 4)  # show all sizes in legend
            ) +
            facet_wrap(~condition2, nrow = 2, ncol=3) +
            theme_tufte() + theme(axis.line=element_line()) + 
            scale_x_continuous(limits=c(0,(max(scatter_plots_d_r$x)+max(scatter_plots_d_r$x)/3))) + scale_y_continuous(limits=c(0,(max(scatter_plots_d_r$y)+max(scatter_plots_d_r$y)/3)))

            scatter_plots_d_s <- scatter_plots_d %>% filter(scatter_plots_d$condition == "Senders")
            sca_s <- ggplot(scatter_plots_d_s, aes(x = x, y = y)) + 
            geom_point(aes(size = Count, fill = labels, color = labels), alpha = .8, shape = 21, stroke = .6) +
            scale_fill_manual(values = colors_celltypes) +
            scale_color_manual(values = colors_celltypes) +
            scale_size_continuous(
              range = c(1, 12),
              breaks = pretty(scatter_plots_d_s$Count, n = 4)  # show all sizes in legend
            ) +
            facet_wrap(~condition2, nrow = 2, ncol=3) +
            theme_tufte() + theme(axis.line=element_line()) + 
            scale_x_continuous(limits=c(0,(max(scatter_plots_d_s$x)+max(scatter_plots_d_s$x)/3))) + scale_y_continuous(limits=c(0,(max(scatter_plots_d_s$y)+max(scatter_plots_d_s$y)/3)))
          print(sca_s)
          print(sca_r)

          dev.off()

          pdf(paste0(outdir2,"/scatter/Scatter_Senders_Receivers_blank.pdf"), width=9, height=6)
          print(plot_blank(sca_s))
          print(plot_blank(sca_r))
          dev.off()

        # Scatter for EGF pathway
        
          # gg <- list()
          # x_axis_scatter <- c()
          # y_axis_scatter <- c()
          
          # for (i in 1:length(cellchat_objects)) {
          #   scatter_plots <- netAnalysis_signalingRole_scatter(cellchat_objects[[i]], title = names(cellchat_objects)[i], color.use=colors_celltypes[order_colors], signaling = "EGF")
          #   x_axis_scatter <- c(x_axis_scatter, scatter_plots$data$x)
          #   y_axis_scatter <- c(y_axis_scatter, scatter_plots$data$y)
          #   gg[[i]] <- scatter_plots
          # }

          # for(i in 1:length(cellchat_objects)){
          #   gg[[i]] <- gg[[i]] + xlim(0, max(x_axis_scatter)) +  ylim(0, max(y_axis_scatter))
          # } 

          # pdf(paste0(outdir2, "/", "/scatter/Scatter_EGF_pathway.pdf"), width=15, height=4)
          # print(patchwork::wrap_plots(plots = gg ))
          # dev.off()


        ## Scatter for all pathways - no labels
          # num.link <- sapply(cellchat_objects, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})

          # weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
          # gg <- list()
          # x_axis_scatter <- c()
          # y_axis_scatter <- c()
          
          # for (i in 1:length(cellchat_objects)) {
          #   order_colors <-levels(cellchat_objects[[i]]@idents)
          #   scatter_plots <- netAnalysis_signalingRole_scatter(cellchat_objects[[i]], title = names(cellchat_objects)[i], weight.MinMax = weight.MinMax, color.use=colors_celltypes[order_colors],do.label=F)
          #   x_axis_scatter <- c(x_axis_scatter, scatter_plots$data$x)
          #   y_axis_scatter <- c(y_axis_scatter, scatter_plots$data$y)
          #   gg[[i]] <- scatter_plots
          # }

          # for(i in 1:length(cellchat_objects)){
          #   gg[[i]] <- gg[[i]] + xlim(0, max(x_axis_scatter)) +  ylim(0, max(y_axis_scatter))
          # } 

          # pdf(paste0(outdir2, "/", "/scatter/Scatter_all_pathways_nolabels.pdf"), width=15, height=3.5)
          # patchwork::wrap_plots(plots = gg )
          # dev.off()

        # Scatter for EGF pathway
        
          # gg <- list()
          # x_axis_scatter <- c()
          # y_axis_scatter <- c()
          
          # for (i in 1:length(cellchat_objects)) {
          #   scatter_plots <- netAnalysis_signalingRole_scatter(cellchat_objects[[i]], title = names(cellchat_objects)[i], color.use=colors_celltypes[order_colors], signaling = "EGF", do.label=F)
          #   x_axis_scatter <- c(x_axis_scatter, scatter_plots$data$x)
          #   y_axis_scatter <- c(y_axis_scatter, scatter_plots$data$y)
          #   gg[[i]] <- scatter_plots
          # }

          # for(i in 1:length(cellchat_objects)){
          #   gg[[i]] <- gg[[i]] + xlim(0, max(x_axis_scatter)) +  ylim(0, max(y_axis_scatter))
          # } 

          # pdf(paste0(outdir2, "/scatter/Scatter_EGF_pathway_nolabels.pdf"), width=15, height=4)
          # print(patchwork::wrap_plots(plots = gg ))
          # dev.off()

        #### Fibroblast plot
        # senders_celltypes <- senders_celltypes[senders_celltypes %in% celltype_dgea_imm]
        # receivers_celltypes <- receivers_celltypes[receivers_celltypes %in% celltype_dgea_imm]
        dir.create(paste0(outdir2, "/chord/"))

            pdf(paste0(outdir2, "/chord/Chord_Senders.pdf"), width=10, height=10)
            # par(mfrow = c(1, 3), xpd=TRUE)
            for (i in 1:length(cellchat_objects)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = senders_celltypes, targets.use = receivers_celltypes, lab.cex = 0.5, title.name = paste0("Signaling from Senders - ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
            }

            dev.off()

            pdf(paste0(outdir2, "/chord/Chord_Receiverz.pdf"), width=10, height=10)
            # par(mfrow = c(1, 3), xpd=TRUE)
            for (i in 1:length(cellchat_objects)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = receivers_celltypes, targets.use = senders_celltypes, lab.cex = 0.5, title.name = paste0("Signaling from Others - ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
            }
            dev.off()

            pdf(paste0(outdir2, "/chord/Chord_fromSelected.pdf"), width=10, height=10)
            # par(mfrow = c(1, 3), xpd=TRUE)
            for (i in 1:length(cellchat_objects)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = celltype_selected, targets.use = receivers_celltypes, lab.cex = 0.5, title.name = paste0("Signaling from ", celltype_selected ,"- ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
            }
            dev.off()

          pdf(paste0(outdir2, "/chord/Chord_toSelected.pdf"), width=10, height=10)
            # par(mfrow = c(1, 3), xpd=TRUE)
            for (i in 1:length(cellchat_objects)) {
              order_colors <-levels(cellchat_objects[[i]]@idents)
              try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = receivers_celltypes , targets.use = celltype_selected, lab.cex = 0.5, title.name = paste0("Signaling Others to ", celltype_selected ,"- ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
            }
            dev.off()
        
          lapply(c("Senders", "Receivers"), function(sender_chors) { 
            #  sender_chors <- "Senders"
             sender_chords <- if(sender_chors == "Senders") senders_celltypes else receivers_celltypes
             receivers_chords <- if(sender_chors == "Senders") receivers_celltypes else senders_celltypes
             to_receiver <- grep(sender_chors, c("Senders", "Receivers"), invert=TRUE, value=TRUE)

            lapply(sender_chords, function(per_sender) {
              # per_sender <- sender_chords[1]
              pdf(paste0(outdir2, "/chord/Chord_", per_sender, "_to_", to_receiver, ".pdf"), width=10, height=10)
              # par(mfrow = c(1, 3), xpd=TRUE)
              for (i in 1:length(cellchat_objects)) {
                order_colors <-levels(cellchat_objects[[i]]@idents)
                try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = per_sender, targets.use = receivers_chords, lab.cex = 0.5, title.name = paste0("Signaling from Senders - ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
              } # HERE

              dev.off()

              pdf(paste0(outdir2, "/chord/Chord_", to_receiver, "_to_", per_sender, ".pdf"), width=10, height=10)
              # par(mfrow = c(1, 3), xpd=TRUE)
              for (i in 1:length(cellchat_objects)) {
                order_colors <-levels(cellchat_objects[[i]]@idents)
                try(netVisual_chord_gene_mod(cellchat_objects[[i]], sources.use = receivers_chords, targets.use = per_sender, lab.cex = 0.5, title.name = paste0("Signaling from Others - ", names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors]))
              }
              dev.off()

            })

          })

        ## End Chords
      
        # ind_outputs <- paste0(outdir2, "/EGF_pathway_specific") 
        # dir.create(ind_outputs) 

        # pathways.show <- c("EGF") 
        # weight.max <- getMaxWeight(cellchat_objects, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

        # pdf(paste0(ind_outputs, "/circle.pdf"), width=10, height=10)
        # par(mfrow = c(1,1), xpd=TRUE)
        # for (i in 1:length(cellchat_objects)) {
        #   order_colors <-levels(cellchat_objects[[i]]@idents)
        #   netVisual_aggregate(cellchat_objects[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors])
        # }
        # dev.off()

        # pdf(paste0(ind_outputs, "/heatmap.pdf"), width=15, height=10)
        # par(mfrow = c(1,1), xpd=TRUE)
        # ht <- list()
        # for (i in 1:length(cellchat_objects)) {
        #   order_colors <-levels(cellchat_objects[[i]]@idents)
        #   ht[[i]] <- netVisual_heatmap(cellchat_objects[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors])
        # }

        # #> Do heatmap based on a single object 
        # #> 
        # #> Do heatmap based on a single object
        # print(ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm")))
        # dev.off()

        # pdf(paste0(ind_outputs, "/chord.pdf"), width=10, height=10)
        # # Chord diagram
        # par(mfrow = c(1,1), xpd=TRUE)
        # for (i in 1:length(cellchat_objects)) {
        #   order_colors <-levels(cellchat_objects[[i]]@idents)
        #   netVisual_aggregate(cellchat_objects[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(cellchat_objects)[i]), color.use=colors_celltypes[order_colors])
        # }
        # dev.off()

      ######## DGEA
      if(is.null(color_edges)){
        color_edges <- colors_celltypes
        color_edges[-which(names(color_edges) %in% celltype_selected)] <- "gray90"
      } 
      outdir2 <- paste0(outdir2, "/dgea")
      dir.create(outdir2)

      cellchat_comparisions(
        cellchat_objects=cellchat_objects,
        groups_to_compare= comparing_groups,
        padj=0.05,
        fc=0.25,
        method = "wilcoxon",
        output_dir = outdir2,
        senders_celltypes=senders_celltypes,
        receivers_celltypes=receivers_celltypes,
        colors_celltypes=colors_celltypes, 
        color_edges=color_edges)

      #https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
      #  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

      cellchat_comparisions(
        cellchat_objects=cellchat_objects,
        groups_to_compare= comparing_groups,
        padj=0.05,
        fc=0.05,
        method = "wilcoxon",
        output_dir = outdir2,
        senders_celltypes=senders_celltypes,
        receivers_celltypes=receivers_celltypes,
        colors_celltypes=colors_celltypes, 
        color_edges=color_edges)

               

}
