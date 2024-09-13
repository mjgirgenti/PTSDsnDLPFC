my_netVisual_chord <- function(net.diff, edge.transparency = FALSE, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                measure = c("weight","count"),
                                layout = c("circle", "chord"),
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,title.cex=1.1,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL, cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20, title.name = NULL,
                                ...) {
    gg <- my_netVisual_chord_cell_internal(net.diff, edge.transparency = edge.transparency, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce, title.cex = title.cex,
                                        title.name = title.name, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
    return(gg)
                                }



my_netVisual_chord_cell_internal <- function(net, edge.transparency = FALSE, color.use = NULL, group = NULL, cell.order = NULL,
                                          sources.use = NULL, targets.use = NULL,
                                          lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
                                          remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
                                          transparency = 0.4, link.border = NA,title.cex=1,
                                          title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
  if (inherits(x = net, what = c("matrix", "Matrix"))) {
    cell.levels <- union(rownames(net), colnames(net))
    net <- reshape2::melt(net, value.name = "prob")
    colnames(net)[1:2] <- c("source","target")
  } else if (is.data.frame(net)) {
    if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
      stop("The input data frame must contain three columns named as source, target, prob")
    }
    cell.levels <- as.character(union(net$source,net$target))
  }
  if (!is.null(cell.order)) {
    cell.levels <- cell.order
  }
  net$source <- as.character(net$source)
  net$target <- as.character(net$target)

  # keep the interactions associated with sources and targets of interest
  if (!is.null(sources.use)){
    if (is.numeric(sources.use)) {
      sources.use <- cell.levels[sources.use]
    }
    net <- subset(net, source %in% sources.use)
  }
  if (!is.null(targets.use)){
    if (is.numeric(targets.use)) {
      targets.use <- cell.levels[targets.use]
    }
    net <- subset(net, target %in% targets.use)
  }
  # remove the interactions with zero values
  net <- subset(net, prob != 0)
  if(dim(net)[1]<=0){message("No interaction between those cells")}
  # create a fake data if keeping the cell types (i.e., sectors) without any interactions
  if (!remove.isolate) {
    cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
    if (length(cells.removed) > 0) {
      net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
      colnames(net.fake) <- colnames(net)
      net <- rbind(net, net.fake)
      link.visible <- net[, 1:2]
      link.visible$plot <- FALSE
      if(nrow(net) > nrow(net.fake)){
        link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
      }
      # directional <- net[, 1:2]
      # directional$plot <- 0
      # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
      # link.arr.type = "big.arrow"
      # message("Set scale = TRUE when remove.isolate = FALSE")
      scale = TRUE
    }
  }

  df <- net
  cells.use <- union(df$source,df$target)

  # define grid order
  order.sector <- cell.levels[cell.levels %in% cells.use]

  # define grid color
  if (is.null(color.use)){
    color.use = scPalette(length(cell.levels))
    names(color.use) <- cell.levels
  } else if (is.null(names(color.use))) {
    names(color.use) <- cell.levels
  }

  grid.col <- color.use[order.sector]
  names(grid.col) <- order.sector
  
  # set grouping information
  if (!is.null(group)) {
    group <- group[names(group) %in% order.sector]
  }

  # # define edge color
  # edge.color <- ifelse(df$prob > 0, '#b2182b', '#2166ac')

    ###############################################
    # Create a custom color palette from red to blue
    # range_value <- max(abs(net$prob))
    range_value <- 0.0655241524091399
    ### Create a custom color palette from red to blue
    # color_ramp <- colorRamp(c("blue", "white", "red"))
    color_ramp <- colorRamp(c("#317EC2", "white", "#C03830")) # milder
    color_ramp2 <- colorRamp2(c(0,range_value),c( "white", "#C03830"))
    # color_ramp <- colorRamp(c("#333333", "white", "#333333")) # all gray
    
    edge.color <- rgb(color_ramp((net$prob + range_value) / (2 * range_value)), maxColorValue = 255)  
    edge.color2 <- color_ramp2
    ###############################################
  
  rgb_matrix <- col2rgb(edge.color) / 255
  if (edge.transparency){
  max_prob<-max(abs(df$prob))
    for (i in 1:dim(df)[1]){
        edge.color[i] <- rgb(rgb_matrix[1,i],rgb_matrix[2,i],rgb_matrix[3,i], abs(df$prob[i]/max_prob))
    }  
  }
  #edge.color <- color.use[as.character(df$source)]
  #return (rgb_matrix)
  if (directional == 0 | directional == 2) {
    link.arr.type = "triangle"
  } else {
    link.arr.type = "big.arrow"
  }

  circos.clear()
  chordDiagram(df,
               order = order.sector,
               col = edge.color,
               grid.col = grid.col,
               transparency = transparency,
               link.border = link.border,
               directional = directional,
               direction.type = c("diffHeight","arrows"),
               link.arr.type = link.arr.type, # link.border = "white",
               annotationTrack = "grid",
               annotationTrackHeight = annotationTrackHeight,
               preAllocateTracks = list(track.height = max(strwidth(order.sector))),
               small.gap = small.gap,
               big.gap = big.gap,
               link.visible = link.visible,
               scale = TRUE,
               group = group,
               link.target.prop = link.target.prop,
               diffHeight = mm_h(4), target.prop.height = mm_h(2),
               reduce = reduce,
               ...)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    xplot = get.cell.meta.data("xplot")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
  }, bg.border = NA)

  # https://jokergoo.github.io/circlize_book/book/legends.html
  if (show.legend) {
    # lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
    # ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))

    # lgd <- ComplexHeatmap::Legend(at = c(min(abs(net$prob)),max(abs(net$prob))), type = "lines", legend_gp = grid::gpar(fill = edge.color), title = "Cell State")
    # ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))

      lgd <- ComplexHeatmap::Legend(title = "prob", at = c(0,0.0655241524091399), col_fun = edge.color2)
      ComplexHeatmap::draw(lgd, x = unit(1, "npc"), y = unit(legend.pos.y, "mm"), just = c("right"))
      # draw(lgd,  just = "right")
  }

  if(!is.null(title.name)){
    # title(title.name, cex = 1)
    text(0, 0.9, title.name, cex=title.cex)
  }
  circos.clear()
  gg <- recordPlot()
  return(gg)
}




# ########################
# netVisual_chord_cell_internal <- function(net, color.use = NULL, group = NULL, cell.order = NULL,
#                                           sources.use = NULL, targets.use = NULL,
#                                           lab.cex = 0.8,small.gap = 1, big.gap = 10, annotationTrackHeight = c(0.03),
#                                           remove.isolate = FALSE, link.visible = TRUE, scale = FALSE, directional = 1, link.target.prop = TRUE, reduce = -1,
#                                           transparency = 0.4, link.border = NA,
#                                           title.name = NULL, show.legend = FALSE, legend.pos.x = 20, legend.pos.y = 20,...){
#   if (inherits(x = net, what = c("matrix", "Matrix"))) {
#     cell.levels <- union(rownames(net), colnames(net))
#     net <- reshape2::melt(net, value.name = "prob")
#     colnames(net)[1:2] <- c("source","target")
#   } else if (is.data.frame(net)) {
#     if (all(c("source","target", "prob") %in% colnames(net)) == FALSE) {
#       stop("The input data frame must contain three columns named as source, target, prob")
#     }
#     cell.levels <- as.character(union(net$source,net$target))
#   }
#   if (!is.null(cell.order)) {
#     cell.levels <- cell.order
#   }
#   net$source <- as.character(net$source)
#   net$target <- as.character(net$target)

#   # keep the interactions associated with sources and targets of interest
#   if (!is.null(sources.use)){
#     if (is.numeric(sources.use)) {
#       sources.use <- cell.levels[sources.use]
#     }
#     net <- subset(net, source %in% sources.use)
#   }
#   if (!is.null(targets.use)){
#     if (is.numeric(targets.use)) {
#       targets.use <- cell.levels[targets.use]
#     }
#     net <- subset(net, target %in% targets.use)
#   }
#   # remove the interactions with zero values
#   net <- subset(net, prob > 0)
#   if(dim(net)[1]<=0){message("No interaction between those cells")}
#   # create a fake data if keeping the cell types (i.e., sectors) without any interactions
#   if (!remove.isolate) {
#     cells.removed <- setdiff(cell.levels, as.character(union(net$source,net$target)))
#     if (length(cells.removed) > 0) {
#       net.fake <- data.frame(cells.removed, cells.removed, 1e-10*sample(length(cells.removed), length(cells.removed)))
#       colnames(net.fake) <- colnames(net)
#       net <- rbind(net, net.fake)
#       link.visible <- net[, 1:2]
#       link.visible$plot <- FALSE
#       if(nrow(net) > nrow(net.fake)){
#         link.visible$plot[1:(nrow(net) - nrow(net.fake))] <- TRUE
#       }
#       # directional <- net[, 1:2]
#       # directional$plot <- 0
#       # directional$plot[1:(nrow(net) - nrow(net.fake))] <- 1
#       # link.arr.type = "big.arrow"
#       # message("Set scale = TRUE when remove.isolate = FALSE")
#       scale = TRUE
#     }
#   }

#   df <- net
#   cells.use <- union(df$source,df$target)

#   # define grid order
#   order.sector <- cell.levels[cell.levels %in% cells.use]

#   # define grid color
#   if (is.null(color.use)){
#     color.use = scPalette(length(cell.levels))
#     names(color.use) <- cell.levels
#   } else if (is.null(names(color.use))) {
#     names(color.use) <- cell.levels
#   }
#   grid.col <- color.use[order.sector]
#   names(grid.col) <- order.sector

#   # set grouping information
#   if (!is.null(group)) {
#     group <- group[names(group) %in% order.sector]
#   }

#   # define edge color
#   ################################################################
#   # edge.color <- color.use[as.character(df$source)]
#    range_value <- max(abs(net$prob))
#     ### Create a custom color palette from red to blue
#     # color_ramp <- colorRamp(c("blue", "white", "red"))
#     color_ramp <- colorRamp(c("#317EC2", "white", "#C03830")) # milder
#     # color_ramp <- colorRamp(c("#333333", "white", "#333333")) # all gray
    
#     edge.color <- rgb(color_ramp((net$prob + range_value) / (2 * range_value)), maxColorValue = 255)  


#   ################################################################

#   if (directional == 0 | directional == 2) {
#     link.arr.type = "triangle"
#   } else {
#     link.arr.type = "big.arrow"
#   }

#   circos.clear()
#   chordDiagram(df,
#                order = order.sector,
#                col = edge.color,
#                grid.col = grid.col,
#                transparency = transparency,
#                link.border = link.border,
#                directional = directional,
#                direction.type = c("diffHeight","arrows"),
#                link.arr.type = link.arr.type, # link.border = "white",
#                annotationTrack = "grid",
#                annotationTrackHeight = annotationTrackHeight,
#                preAllocateTracks = list(track.height = max(strwidth(order.sector))),
#                small.gap = small.gap,
#                big.gap = big.gap,
#                link.visible = link.visible,
#                scale = scale,
#                group = group,
#                link.target.prop = link.target.prop,
#                reduce = reduce,
#                ...)
#   circos.track(track.index = 1, panel.fun = function(x, y) {
#     xlim = get.cell.meta.data("xlim")
#     xplot = get.cell.meta.data("xplot")
#     ylim = get.cell.meta.data("ylim")
#     sector.name = get.cell.meta.data("sector.index")
#     circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
#   }, bg.border = NA)

#   # https://jokergoo.github.io/circlize_book/book/legends.html
#   if (show.legend) {
#     lgd <- ComplexHeatmap::Legend(at = names(grid.col), type = "grid", legend_gp = grid::gpar(fill = grid.col), title = "Cell State")
#     ComplexHeatmap::draw(lgd, x = unit(1, "npc")-unit(legend.pos.x, "mm"), y = unit(legend.pos.y, "mm"), just = c("right", "bottom"))
#   }

#   if(!is.null(title.name)){
#     # title(title.name, cex = 1)
#     text(-0, 1.02, title.name, cex=1)
#   }
#   circos.clear()
#   gg <- recordPlot()
#   return(gg)
# }