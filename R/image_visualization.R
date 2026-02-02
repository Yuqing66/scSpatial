#' Visualize categorical metadata and transcripts in spatial context (centroids).
#'
#' @description Similar to Seurat's ImageDimPlot but allows highlighting specific
#' categories and optionally overlaying molecule locations. Supports faceting by
#' metadata column, cropping, and adding a scale bar.
#' @param object Seurat object with spatial images and centroids.
#' @param fov Name of the field of view to plot (must exist in `object@images`).
#' @param group.by Metadata column used to color cells. Defaults to `ident`.
#' @param split.by Metadata column used to facet the plot.
#' @param size Base point size for cells.
#' @param cols Named vector with group categories as names and colors as values.
#' @param alpha Transparency of the cells.
#' @param highlight.by Metadata column containing highlight category labels.
#'   Defaults to `group.by`.
#' @param highlight.groups Vector specifying categories to highlight.
#' @param highlight.size Point size for highlighted cells.
#' @param highlight.cols Named vector with highlight category labels as names and colors as values.
#' @param highlight.alpha Transparency of highlighted cell points.
#' @param molecules Character vector of molecules to plot.
#' @param molecules.size Point size of molecules.
#' @param molecules.cols Named vector with molecules as names and colors as values.
#' @param molecules.alpha Transparency of molecules to plot.
#' @param dark.background Logical; use a black panel background. By default TRUE.
#' @param crop If TRUE, crop to the extent of cells/molecules. Or a numeric vector
#'   with four numbers specifying `c(min.x, max.x, min.y, max.y)`.
#' @param flip Logical; kept for compatibility with Seurat output (currently unused).
#' @param scalebar.length Optional positive numeric length of the scale bar in image units.
#' @param scalebar.numConv Numeric conversion factor applied to the scale bar label.
#' @param scalebar.unit Unit label appended to the scale bar text.
#' @param scalebar.position One of `bottomright`, `bottomleft`, `topright`, `topleft`,
#'   or a numeric length-2 vector specifying the `(x, y)` start position.
#' @param scalebar.color Color for the scale bar and label. Defaults to white on
#'   dark background and black otherwise.
#' @param scalebar.text.size Text size for the scale bar label.
#' @param scalebar.margin Margin as a fraction of the plot span when positioning the scale bar.
#' @return A ggplot object.
#' @examples
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", split.by = "Disease")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", highlight.groups = c("Lymphocyte","KC_glandular_eccrine"))
#' srt$lowCount <- ifelse(srt$nCount_seqFISH < 100, "low", "high")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", highlight.by = "lowCount", highlight.groups = "low", highlight.size = 0.5)
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", molecules = c("IFNB1","IL17A"), molecules.size = 0.7)
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", crop = T)
#' @seealso Seurat::ImageDimPlot
#' @import ggplot2 dplyr Seurat
#' @export


ImageDimPlot.ssc <- function (object, fov, group.by = NULL, split.by = NULL, size = 0.1, 
                              cols = NULL, alpha = 1, highlight.by = NULL, highlight.groups = NULL, 
                              highlight.size = 0.2, highlight.cols = NULL, highlight.alpha = 1, 
                              molecules = NULL, molecules.size = 0.1, molecules.cols = NULL, 
                              molecules.alpha = 1, dark.background = T, crop = NULL, flip = F, 
                              scalebar.length = NULL, scalebar.numConv = 1, scalebar.unit = NULL, scalebar.position = "bottomright", 
                              scalebar.color = NULL, scalebar.text.size = 3, scalebar.margin = 0.03) 
{
  if (is.null(group.by)) {
    group.by <- "ident"
  }
  
  fov_image <- object@images[[fov]]
  coords <- fov_image$centroids@coords
  df <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    stringsAsFactors = FALSE
  )
  
  ind.fov <- match(fov_image$centroids@cells, colnames(object))
  df.meta <- FetchData(object = object, vars = c(group.by, split.by), cells = ind.fov)
  
  if (is.null(highlight.by)) {
    highlight.by <- group.by
  } else if (!identical(highlight.by, group.by)) {
    df.meta[[highlight.by]] <- object@meta.data[ind.fov, highlight.by]
  }
  
  df.meta$highlight <- ifelse(df.meta[[highlight.by]] %in% highlight.groups, "y", "n")
  df <- cbind(df, df.meta)
  
  if (!is.null(molecules)) {
    valid_molecules <- molecules[molecules %in% rownames(object)]
    missing <- setdiff(molecules, valid_molecules)
    if (length(missing) > 0) {
      message("Molecules not in the object: ", paste(missing, collapse = ","))
    }
    if (length(valid_molecules) > 0) {
      molecule_coords <- fov_image@molecules$molecules[valid_molecules]
      df.mol <- do.call(rbind, lapply(names(molecule_coords), function(name) {
        mol_coords <- molecule_coords[[name]]@coords
        data.frame(
          x = mol_coords[, 1],
          y = mol_coords[, 2],
          mol = name,
          stringsAsFactors = FALSE
        )
      }))
      colnames(df.mol)[3] <- group.by
      
      if (!is.null(split.by)) {
        split_values <- unique(df[[split.by]])
        df_mol_base <- df.mol
        df.mol <- do.call(rbind, lapply(split_values, function(val) {
          tmp <- df_mol_base
          tmp[[split.by]] <- val
          tmp
        }))
      }
      
      df.mol$highlight <- "m"
      df <- rbind(df, df.mol)
    }
  }
  
  plot_xlim <- range(df$x, na.rm = TRUE)
  plot_ylim <- range(df$y, na.rm = TRUE)
  
  g <- ggplot(df) + 
    geom_point(aes_string(x = "x", y = "y", alpha = "highlight", size = "highlight", colour = group.by), 
               shape = 16) + 
    theme_classic() + 
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
          axis.line = element_blank(), panel.grid = element_blank()) + 
    guides(colour = guide_legend(title = group.by), size = "none", alpha = "none") + 
    scale_size_manual(values = c(y = highlight.size, n = size, m = molecules.size)) + 
    scale_alpha_manual(values = c(y = highlight.alpha, n = alpha, m = molecules.alpha))
  
  g2 <- ggplot_build(g)
  color_scale <- g2$plot$scales$get_scales("colour")
  cols.default <- color_scale$palette.cache
  names(cols.default) <- sort(unique(df[[group.by]]))
  
  effective.cols <- cols
  if (!is.null(highlight.cols)) {
    if (is.null(effective.cols)) {
      effective.cols <- cols.default
    }
    effective.cols[names(highlight.cols)] <- highlight.cols
  }
  if (!is.null(molecules.cols)) {
    if (is.null(effective.cols)) {
      effective.cols <- cols.default
    }
    effective.cols[names(molecules.cols)] <- molecules.cols
  }
  if (!is.null(effective.cols)) {
    g <- g + scale_color_manual(values = effective.cols)
  }
  
  if (!is.null(split.by)) {
    g <- g + facet_wrap(reformulate(split.by))
  }
  
  coord_args <- list(ratio = 1)
  if (!is.null(crop)) {
    if (identical(crop, TRUE)) {
      crop <- c(min(df$x) - 1, max(df$x) + 1, min(df$y) - 1, max(df$y) + 1)
    } else if (length(crop) != 4) {
      stop("crop must be TRUE or a numeric vector c(min.x, max.x, min.y, max.y).")
    }
    plot_xlim <- crop[1:2]
    plot_ylim <- crop[3:4]
  }
  
  g <- g + coord_fixed(ratio = 1, xlim = plot_xlim, ylim = plot_ylim)
  
  if (dark.background) {
    g <- g + theme(panel.background = element_rect(fill = "black", colour = "black"))
  }
  
  if (!is.null(scalebar.length)) {
    if (!is.numeric(scalebar.length) || length(scalebar.length) != 1 || scalebar.length <= 0) {
      stop("scalebar.length must be a positive numeric value.")
    }
    
    span_x <- diff(plot_xlim)
    span_y <- diff(plot_ylim)
    margin <- max(scalebar.margin, 0)
    margin_x <- span_x * margin
    margin_y <- span_y * margin
    offset_y <- if (span_y > 0) span_y * 0.02 else span_x * 0.02
    bar_colour <- if (is.null(scalebar.color)) {
      if (isTRUE(dark.background)) "white" else "black"
    } else {
      scalebar.color
    }
    
    if (is.character(scalebar.position)) {
      pos <- match.arg(scalebar.position, c("bottomright", "bottomleft", "topright", "topleft"))
      if (scalebar.length > span_x) {
        warning("scalebar.length exceeds the x-range of the plot and may be clipped.", call. = FALSE)
      }
      
      from_left <- grepl("left", pos)
      from_bottom <- grepl("bottom", pos)
      
      x_start <- if (from_left) plot_xlim[1] + margin_x else plot_xlim[2] - margin_x - scalebar.length
      y_start <- if (from_bottom) plot_ylim[1] + margin_y else plot_ylim[2] - margin_y
      label_y <- if (from_bottom) y_start + offset_y else y_start - offset_y
      text_vjust <- if (from_bottom) 0 else 1
    } else if (is.numeric(scalebar.position) && length(scalebar.position) == 2) {
      x_start <- scalebar.position[1]
      y_start <- scalebar.position[2]
      label_y <- y_start + offset_y
      text_vjust <- 0
    } else {
      stop("scalebar.position must be 'bottomright', 'bottomleft', 'topright', 'topleft', or a numeric length-2 vector.")
    }
    
    x_end <- x_start + scalebar.length
    scalebar.label <- round(scalebar.length * scalebar.numConv, digits = 2)
    label <- if (is.null(scalebar.unit) || scalebar.unit == "") {
      scalebar.label
    } else {
      paste(scalebar.label, scalebar.unit)
    }
    
    g <- g +
      annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_start, colour = bar_colour, size = 0.5) +
      annotate("text", x = (x_start + x_end) / 2, y = label_y, label = label, colour = bar_colour, size = scalebar.text.size, vjust = text_vjust)
  }
  
  if (flip) {
    g <- g + coord_flip()
  }
  
  g
}



#' Transition animation between plots of cells
#'
#' @description Create a transition animation between different cell coordinate plots, e.g. fov, umap, pca
#' @param object Seurat object.
#' @param initial Name of cell coordinates to start with. Embedding name in srt@reductions, or "images" for spatial slide.
#' @param final Name of cell coordinates to transition into.
#' @param group.by The name of metadata column to color by.
#' @param fov Slides to include. NULL for all fovs or a vector of slide names.
#' @param ncol Number of fovs in a row.
#' @param fov.size c(width, height) to specify the size of fov panel grid.
#' @param point.size Point size for cells.
#' @param point.alpha Point transparency.
#' @param image.coord.flip To match the Seurat output, the x and y axises are by default flipped.
#' @param match.scale.by Scale the initial and final coordinates. NULL, x, y, or both.
#' @param plot.legend Logical. If TRUE, plot the legend. By default FALSE.
#' @param plot.boundaryFrames Logical. If TRUE, plot the boundary frames, and use list as the output. By default TRUE.
#' @examples
#' g <- plotTransition(srt, initial="umap", final="images", group.by="subcelltype", fov=NULL, point.size = 0.3, match.scale.by = "y")
#' animate(g, duration = 7, fps = 10, start_pause = 10, end_pause = 10,
#'         width = 16, height = 8, units = "in", res = 150,
#'         renderer = gifski_renderer())
#' anim_save("transition_animation_umap_to_images.gif")
#' @import ggplot2 dplyr gganimate
#' @export

plotTransition <- function(srt, initial="umap", final="images", group.by=NULL,
                           fov=NULL, ncol=NULL, fov.size=NULL,
                           point.size=0.5, point.alpha=1, image.coord.flip=T, 
                           match.scale.by="y", plot.legend=FALSE, plot.boundaryFrames=T,
                           cols=NULL){
  coords_initial <- getCoords(srt, name = initial, fov=fov, metadata=group.by)
  coords_final <- getCoords(srt, name = final, fov=fov, metadata=group.by)
  #
  if (nrow(coords_initial) != nrow(coords_final)){
    cells.intersect <- intersect(rownames(coords_initial), rownames(coords_final))
    coords_initial <- coords_initial[cells.intersect,]
    coords_final <- coords_final[cells.intersect,]
  }
  # special process for coordinates from images
  if ("images" %in% c(initial, final)){
    dfname <- paste0("coords_",c("initial","final")[c(initial,final) == "images"])
    df <- get(dfname)
    if (image.coord.flip){
      df[,c(1,2)] <- df[,c(2,1)]
    }
    # change image coordination to avoid overlap
    fovs <- unique(df$fov)
    if (length(fovs) > 1){
      if (is.null(ncol)){
        ncol <- ceiling(sqrt(length(fovs)))
      }
      nrow <- ceiling(length(fovs)/ncol)
      if (is.null(fov.size)){
        maxx <- max(tapply(df[,1], df$fov, max)) + 1
        maxy <- max(tapply(df[,2], df$fov, max)) + 1
        fov.size <- c(maxx, maxy)
      }
      for (i in 1:length(fovs)){
        fovname <- fovs[i]
        ind <- df$fov == fovname
        df[ind,1] <- df[ind,1] + fov.size[1]*((i-1) %% ncol)
        df[ind,2] <- df[ind,2] + fov.size[2]*(nrow - floor((i-1)/ncol))
      }
    }
    assign(dfname, df)
  }

  if (match.scale.by=="x"){
    i.min <- min(coords_initial[,1])
    i.max <- max(coords_initial[,1])
    f.min <- min(coords_final[,1])
    f.max <- max(coords_final[,1])
    i.d <- i.max - i.min
    f.d <- f.max - f.min
    i.min.y <- min(coords_initial[,2])
    f.min.y <- min(coords_final[,2])
    if (i.d > f.d){
      coords_final[,1] <- (coords_final[,1] - f.min)/f.d*i.d + i.min
      coords_final[,2] <- (coords_final[,2] - f.min.y)/f.d*i.d + i.min.y
    }else{
      coords_initial[,1] <- (coords_initial[,1] - i.min)/i.d*f.d + f.min
      coords_initial[,2] <- (coords_initial[,2] - i.min.y)/i.d*f.d + f.min.y
    }
  }else if (match.scale.by=="y"){
    i.min <- min(coords_initial[,2])
    i.max <- max(coords_initial[,2])
    f.min <- min(coords_final[,2])
    f.max <- max(coords_final[,2])
    i.d <- i.max - i.min
    f.d <- f.max - f.min
    i.min.x <- min(coords_initial[,1])
    f.min.x <- min(coords_final[,1])
    if (i.d > f.d){
      coords_final[,2] <- (coords_final[,2] - f.min)/f.d*i.d + i.min
      scale.final.x <- (coords_final[,1] - f.min.x)/f.d*i.d
      coords_final[,1] <- scale.final.x + (i.d-(max(scale.final.x) - min(scale.final.x)))/2
    }else{
      coords_initial[,2] <- (coords_initial[,2] - i.min)/i.d*f.d + f.min
      scale.initial.x <- (coords_initial[,1] - i.min.x)/i.d*f.d
      coords_initial[,1] <- scale.initial.x + ((max(coords_final[,1]) - min(coords_final[,1]))-(max(scale.initial.x) - min(scale.initial.x)))/2
    }
  }else if (match.scale.by=="both"){
    for (i in 1:2){
      i.min <- min(coords_initial[,i])
      i.max <- max(coords_initial[,i])
      f.min <- min(coords_final[,i])
      f.max <- max(coords_final[,i])
      i.d <- i.max - i.min
      f.d <- f.max - f.min
      if (i.d > f.d){
        coords_final[,i] <- (coords_final[,i] - f.min)/f.d*i.d + i.min
      }else{
        coords_initial[,i] <- (coords_initial[,i] - i.min)/i.d*f.d + f.min
      }
    }
  }



  # animation plot
  colnames(coords_initial)[c(1,2)] <- c("x","y")
  coords_initial$frame <- 1
  colnames(coords_final)[c(1,2)] <- c("x","y")
  coords_final$frame <- 2
  animation.data <- rbind(coords_initial, coords_final)

  aes_mapping <- aes(x = x, y = y)
  if (!is.null(group.by)) {
    aes_mapping <- modifyList(aes_mapping, aes(color = !!sym(group.by)))
  }
  
  g <- ggplot(animation.data, aes_mapping) +
    geom_point(size=point.size, alpha=point.alpha) +
    coord_fixed(ratio = 1) +
    transition_time(frame) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  g.list <- list()
  g.list[[1]] <- g
  
  if (plot.boundaryFrames){
    for (i in 2:3){
      ind.frame <- i-1
      g.list[[i]] <- ggplot(animation.data[animation.data$frame==ind.frame,], 
                            aes_mapping) +
        geom_point(size=point.size, alpha=point.alpha) +
        coord_fixed(ratio = 1) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())
    }
    names(g.list) <- c("animation","initial","final")
  }
  
  if (!plot.legend){
    g.list <- lapply(g.list, function(x) x + theme(legend.position = "none"))
  }
  if (!is.null(cols)) {
    g.list <- lapply(g.list, function(x) x + scale_color_manual(values = cols))
  }

  if (length(g.list) == 1){
    return(g.list[[1]])
  }
  
  return(g.list)
}



#' Plot the initial frame and the final frame of the transition animation created by plotTransition.
#'
#' @description
#' Use the same parameters as in the plotTransition that generated the animation.
#' @param object Seurat object.
#' @param initial Name of cell coordinates to start with. Embedding name in srt@reductions, or "images" for spatial slide.
#' @param final Name of cell coordinates to transition into.
#' @param group.by The name of metadata column to color by.
#' @param fov Slides to include. "all" or a vector of slide names.
#' @param ncol Number of fovs in a row.
#' @param fov.size c(width, height) to specify the size of fov panel grid.
#' @param point.size Point size for cells.
#' @param point.alpha Point transparency.
#' @param image.coord.flip To match the Seurat output, the x and y axises are by default flipped.
#' @param match.scale.by Scale the initial and final coordinates. NULL, x, y, or both.
#' @examples
#' g.list <- plotTransition.terminal(srt, initial="umap", final="images", group.by="subcelltype", fov=NULL, point.size = 0.3, match.scale.by = "y")
#' g.list[[1]]
#' ggsave(paste0("transition_animation_umap_to_images_start.png"), width = 16, height = 8)
#' g.list[[2]]
#' ggsave(paste0("transition_animation_umap_to_images_end.png"), width = 16, height = 8)
#' @import ggplot2 dplyr gganimate
#' @export

plotTransition.terminal <- function(srt, initial="umap", final="images", group.by=NULL,
                                    fov=NULL, ncol=NULL, fov.size=NULL,
                                    point.size=0.5, point.alpha=1, image.coord.flip=T, match.scale.by="y"){
  coords_initial <- getCoords(srt, name = initial, fov=fov, metadata=group.by)
  coords_final <- getCoords(srt, name = final, fov=fov, metadata=group.by)
  #
  if (nrow(coords_initial) != nrow(coords_final)){
    cells.intersect <- intersect(rownames(coords_initial), rownames(coords_final))
    coords_initial <- coords_initial[cells.intersect,]
    coords_final <- coords_final[cells.intersect,]
  }
  # special process for coordinates from images
  if ("images" %in% c(initial, final)){
    dfname <- paste0("coords_",c("initial","final")[c(initial,final) == "images"])
    df <- get(dfname)
    if (image.coord.flip){
      df[,c(1,2)] <- df[,c(2,1)]
    }
    # change image coordination to avoid overlap
    fovs <- unique(df$fov)
    if (length(fovs) > 1){
      if (is.null(ncol)){
        ncol <- ceiling(sqrt(length(fovs)))
      }
      nrow <- ceiling(length(fovs)/ncol)
      if (is.null(fov.size)){
        maxx <- max(tapply(df[,1], df$fov, max)) + 1
        maxy <- max(tapply(df[,2], df$fov, max)) + 1
        fov.size <- c(maxx, maxy)
      }
      for (i in 1:length(fovs)){
        fovname <- fovs[i]
        ind <- df$fov == fovname
        df[ind,1] <- df[ind,1] + fov.size[1]*((i-1) %% ncol)
        df[ind,2] <- df[ind,2] + fov.size[2]*(nrow - floor((i-1)/ncol))
      }
    }
    assign(dfname, df)
  }

  if (match.scale.by=="x"){
    i.min <- min(coords_initial[,1])
    i.max <- max(coords_initial[,1])
    f.min <- min(coords_final[,1])
    f.max <- max(coords_final[,1])
    i.d <- i.max - i.min
    f.d <- f.max - f.min
    i.min.y <- min(coords_initial[,2])
    f.min.y <- min(coords_final[,2])
    if (i.d > f.d){
      coords_final[,1] <- (coords_final[,1] - f.min)/f.d*i.d + i.min
      coords_final[,2] <- (coords_final[,2] - f.min.y)/f.d*i.d + i.min.y
    }else{
      coords_initial[,1] <- (coords_initial[,1] - i.min)/i.d*f.d + f.min
      coords_initial[,2] <- (coords_initial[,2] - i.min.y)/i.d*f.d + f.min.y
    }
  }else if (match.scale.by=="y"){
    i.min <- min(coords_initial[,2])
    i.max <- max(coords_initial[,2])
    f.min <- min(coords_final[,2])
    f.max <- max(coords_final[,2])
    i.d <- i.max - i.min
    f.d <- f.max - f.min
    i.min.x <- min(coords_initial[,1])
    f.min.x <- min(coords_final[,1])
    if (i.d > f.d){
      coords_final[,2] <- (coords_final[,2] - f.min)/f.d*i.d + i.min
      scale.final.x <- (coords_final[,1] - f.min.x)/f.d*i.d
      coords_final[,1] <- scale.final.x + (i.d-(max(scale.final.x) - min(scale.final.x)))/2
    }else{
      coords_initial[,2] <- (coords_initial[,2] - i.min)/i.d*f.d + f.min
      scale.initial.x <- (coords_initial[,1] - i.min.x)/i.d*f.d
      coords_initial[,1] <- scale.initial.x + ((max(coords_final[,1]) - min(coords_final[,1]))-(max(scale.initial.x) - min(scale.initial.x)))/2
    }
  }else if (match.scale.by=="both"){
    for (i in 1:2){
      i.min <- min(coords_initial[,i])
      i.max <- max(coords_initial[,i])
      f.min <- min(coords_final[,i])
      f.max <- max(coords_final[,i])
      i.d <- i.max - i.min
      f.d <- f.max - f.min
      if (i.d > f.d){
        coords_final[,i] <- (coords_final[,i] - f.min)/f.d*i.d + i.min
      }else{
        coords_initial[,i] <- (coords_initial[,i] - i.min)/i.d*f.d + f.min
      }
    }
  }


  # animation plot
  colnames(coords_initial)[c(1,2)] <- c("x","y")
  coords_initial$frame <- 1
  colnames(coords_final)[c(1,2)] <- c("x","y")
  coords_final$frame <- 2
  animation.data <- rbind(coords_initial, coords_final)

  g <- list()
  g[[1]] <- ggplot(animation.data[animation.data$frame==1,], aes(x=x,y=y,col=!!sym(group.by))) +
    geom_point(size=point.size, alpha=point.alpha) +
    coord_fixed(ratio = 1) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  g[[2]] <- ggplot(animation.data[animation.data$frame==2,], aes(x=x,y=y,col=!!sym(group.by))) +
    geom_point(size=point.size, alpha=point.alpha) +
    coord_fixed(ratio = 1) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())

  return(g)
}



#' Plot the path on top of a image plot.
#' @param ggplot.object A ggplot object.
#' @param path.coords A data frame with x and y coordinates of the path.
#' @param ends.close If TRUE, close the path by connecting the last point to the first point.
#' @param path.col Color of the path.
#' @param path.size Size of the path.
#' @param path.alpha Transparency of the path.
#' @return A ggplot object with the path added.
#' @examples
#' g <- ImageDimPlot.ssc(srt, fov = fov, group.by = NULL)
#' coords <- getClickCoordinates(g)
#' ImageDimPlot.path(g, path.coords = coords)
#' @import ggplot2
#' @export


ImageDimPlot.path <- function(ggplot.object, path.coords, ends.close = F,
                              path.col = "red", path.size = 0.5, path.alpha = 1){
  if (ends.close){
    path.coords <- rbind(path.coords, path.coords[1,])
  }
  g <- ggplot.object + geom_path(data = path.coords, aes(x=x,y=y), color = path.col, size = path.size, alpha = path.alpha)
  return(g)
}












library(spatstat.geom)
library(spatstat.explore)
library(dplyr)
library(ggplot2)
library(viridisLite)  # for palettes if needed
#' @title Contour plot of spatial feature expression
#'
#' @description
#' Plot filled contours of a feature's spatial expression by applying kernel
#' density estimation to point-marked expression values in a field of view.
#'
#' @param object Seurat object containing spatial coordinates and expression data.
#' @param feature Feature (e.g., gene) name to visualize.
#' @param fov Field of view (image) name to plot.
#' @param plot.all If TRUE, plot all cell centroids beneath the contours.
#' @param group.by Metadata column name for coloring cell points when `plot.all = TRUE`.
#' @param size Point size for cell centroids.
#' @param cols Named vector of colors for `group.by`.
#' @param alpha Transparency for cell centroids.
#' @param sigma Bandwidth for kernel density; defaults to `bw.diggle(pp)` when NULL.
#' @param n_levels Number of contour levels to draw.
#' @param color_palette Vector of colors to interpolate for filled contours.
#' @param threshold Quantile threshold (0-1) for contouring the highest values.
#' @param contour.alpha Transparency for contour fills.
#' @param dark.background If TRUE, use a dark plot background and palette.
#' @param flip If TRUE, flip axes to match Seurat plotting orientation.
#' @param scalebar.length Length of the scalebar in coordinate units; if NULL, no scalebar is drawn.
#' @param scalebar.numConv Multiplicative conversion from coordinate units to display units.
#' @param scalebar.unit Unit label to append to the scalebar.
#' @param scalebar.position Position of the scalebar; one of "bottomright", "bottomleft",
#'   "topright", "topleft", or a numeric length-2 vector.
#' @param scalebar.color Color of the scalebar; defaults to white on dark background and black otherwise.
#' @param scalebar.text.size Text size for the scalebar label.
#' @param scalebar.margin Fractional margin from plot edges for automatic positions.
#' @return A ggplot object.
#' @examples
#' g <- ImageFeaturePlot.contour(seqfish, feature = "IFNB1", fov = "HS009fov1")
#' 
ImageFeaturePlot.contour <- function(object, feature, fov, 
                                     plot.all = T, group.by = NULL, size = 0.1, cols = NULL, alpha = 0.5,
                                     sigma = NULL, n_levels = 6, color_palette = NULL, threshold = 0.9, contour.alpha = 0.7,
                                     dark.background = T, flip = F,
                                     scalebar.length = NULL, scalebar.numConv = 1, scalebar.unit = NULL, scalebar.position = "bottomright",
                                     scalebar.color = NULL, scalebar.text.size = 3, scalebar.margin = 0.03){
  
  df <- getCoords.cell(object, fov = fov)
  df$expr <- object@assays$seqFISH@layers$data[rownames(object) == feature, match(rownames(df), colnames(object)), drop = TRUE]
  
  if (!is.null(group.by)){
    df$group.by <- object@meta.data[,group.by]
  }
  
  # 1. Define observation window (rectangle; replace with polygon mask if you have tissue boundary)
  win <- owin(range(df$x), range(df$y))
  
  # 2. Point pattern with marks = expression
  pp <- ppp(df$x, df$y, window = win, marks = df$expr)
  
  # 3. Choose bandwidth (sigma). Automatic selectors often oversmooth / undersmooth for sparse genes.
    # or set manually, e.g. sigma <- 40  (in your spatial units)
  sigma <- sigma %||% bw.diggle(pp)
  
  # 4a. Local transcript *sum* surface (kernel-weighted)
  sum_surface <- density(pp, weights = pp$marks, sigma = sigma, at = "pixels", edge=TRUE)
  # units: (sum of expression) / area
  
  df_sum  <- as.data.frame(sum_surface)  # columns: x, y, value
  
  # only plot the top n% of the values
  thr <- quantile(df_sum$value, threshold, na.rm=TRUE)
  maxv <- max(df_sum$value, na.rm=TRUE)
  breaks <- seq(thr, maxv, length.out = n_levels + 1)
  
  plot_xlim <- range(df$x, na.rm = TRUE)
  plot_ylim <- range(df$y, na.rm = TRUE)

  g <- ggplot() + theme_void()
  
  if (plot.all){
    if (is.null(group.by)){
      g <- g + geom_point(data = df, aes(x, y), color="grey", size=size, alpha = alpha)
    }else{
      g <- g + geom_point(data = df, aes(x, y, col = group.by), size=size, alpha = alpha)
      if (!is.null(cols)){
        g <- g + scale_color_manual(values = cols, drop = FALSE)
      }
    }
  }
  
  g <- g + geom_contour_filled(data = df_sum, aes(x=x, y=y, z = value), alpha = contour.alpha, breaks = breaks)

  if (dark.background){
    g <- g + theme(panel.background = element_rect(fill = "black", color = "black"))
    color_palette = color_palette %||% c("#3B0047","#D700FF")
  }else{
    color_palette = color_palette %||% c("#F7E7E6","red")
  }
  g <- g + scale_fill_manual(
    values = colorRampPalette(color_palette)(length(breaks) - 1),
    name = paste0("value ≥ ", signif(thr,3)),
    drop = FALSE)
  
  g <- g + coord_fixed(ratio = 1, xlim = plot_xlim, ylim = plot_ylim)

  if (!is.null(scalebar.length)) {
    if (!is.numeric(scalebar.length) || length(scalebar.length) != 1 || scalebar.length <= 0) {
      stop("scalebar.length must be a positive numeric value.")
    }

    span_x <- diff(plot_xlim)
    span_y <- diff(plot_ylim)
    margin <- max(scalebar.margin, 0)
    margin_x <- span_x * margin
    margin_y <- span_y * margin
    offset_y <- if (span_y > 0) span_y * 0.02 else span_x * 0.02
    bar_colour <- if (is.null(scalebar.color)) {
      if (isTRUE(dark.background)) "white" else "black"
    } else {
      scalebar.color
    }

    if (is.character(scalebar.position)) {
      pos <- match.arg(scalebar.position, c("bottomright", "bottomleft", "topright", "topleft"))
      if (scalebar.length > span_x) {
        warning("scalebar.length exceeds the x-range of the plot and may be clipped.", call. = FALSE)
      }

      from_left <- grepl("left", pos)
      from_bottom <- grepl("bottom", pos)

      x_start <- if (from_left) plot_xlim[1] + margin_x else plot_xlim[2] - margin_x - scalebar.length
      y_start <- if (from_bottom) plot_ylim[1] + margin_y else plot_ylim[2] - margin_y
      label_y <- if (from_bottom) y_start + offset_y else y_start - offset_y
      text_vjust <- if (from_bottom) 0 else 1
    } else if (is.numeric(scalebar.position) && length(scalebar.position) == 2) {
      x_start <- scalebar.position[1]
      y_start <- scalebar.position[2]
      label_y <- y_start + offset_y
      text_vjust <- 0
    } else {
      stop("scalebar.position must be 'bottomright', 'bottomleft', 'topright', 'topleft', or a numeric length-2 vector.")
    }

    x_end <- x_start + scalebar.length
    scalebar.label <- round(scalebar.length * scalebar.numConv, digits = 2)
    label <- if (is.null(scalebar.unit) || scalebar.unit == "") {
      scalebar.label
    } else {
      paste(scalebar.label, scalebar.unit)
    }

    g <- g +
      annotate("segment", x = x_start, xend = x_end, y = y_start, yend = y_start, colour = bar_colour, size = 0.5) +
      annotate("text", x = (x_start + x_end) / 2, y = label_y, label = label, colour = bar_colour, size = scalebar.text.size, vjust = text_vjust)
  }

  if (flip) {
    g <- g + coord_flip()
  }

  return(g)
  
}












getImageSize <- function(object, fov){
  mat <- srt@images[[fov]]@boundaries$centroids@bbox
  w <- mat[1,2] - mat[1,1]
  h <- mat[2,2] - mat[2,1]
  
  res <- c(w,h)
  names(res) <- c("width","height")
  return(res)
}



