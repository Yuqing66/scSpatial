#' Visualize categorical metadata and transcripts in spatial context (centroids).
#'
#' @description Similar to ImageDimPlot in Seurat but allowing highlight of specific categories in group.by.
#' @param object Seurat object.
#' @param fov The fov to plot.
#' @param group.by Column name in the metadata to group color the cells
#' @param split.by Column name in the metadata to split the cells
#' @param size Size of the cells.
#' @param cols A named vector with group categories as names and colors as values.
#' @param alpha Transparency of the cells.
#' @param highlight.by The name of metadata column containing the highlight category labels. Same as group.by by default.
#' @param highlight.groups A vector specifying categories to highlight.
#' @param highlight.size Point size for highlighted cells.
#' @param highlight.cols A named vector with highlight category labels as names and colors as values.
#' @param highlight.alpha Transparency of highlighed cell points.
#' @param molecules Molecules to plot the location.
#' @param molecules.size Point size of molecules.
#' @param molecules.cols A named vector with molecules as names and colors as values.
#' @param molecules.alpha Transparency of molecules to plot.
#' @param dark.backgrount Black backgrount. By default TRUE.
#' @param crop If TRUE, crop off the region with no cells or molecules to plot. Or a vector with four numbers specifying c(min.x, max.x, min.y, max.y)
#' @param flip To match the Seurat output, the x and y axises are by default flipped.
#' @examples
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", split.by = "Disease")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", highlight.groups = c("Lymphocyte","KC_glandular_eccrine"))
#' srt$lowCount <- ifelse(srt$nCount_seqFISH < 100, "low", "high")
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", highlight.by = "lowCount", highlight.groups = "low", highlight.size = 0.5)
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", molecules = c("IFNB1","IL17A"), molecules.size = 0.7)
#' ImageDimPlot.ssc(srt, fov = "UV109fov1", group.by = "celltype_res0.1", crop = T)
#' @import ggplot2 dplyr Seurat
#' @export


ImageDimPlot.ssc <- function(object, fov, group.by=NULL, split.by=NULL, size=0.1, cols=NULL, alpha=1,
                             highlight.by=NULL, highlight.groups=NULL, highlight.size=0.2, highlight.cols=NULL, highlight.alpha=1,
                             molecules=NULL, molecules.size=0.1, molecules.cols=NULL, molecules.alpha=1,
                             dark.background=T, crop=NULL, flip=F){
  # get coordinates and metadata
  df <- data.frame(x=object@images[[fov]]$centroids@coords[,1],
                   y=object@images[[fov]]$centroids@coords[,2])

  if (is.null(group.by)){
    group.by <- "ident"
  }
  ind.fov <- match(object@images[[fov]]$centroids@cells,colnames(object))
  df.meta <- FetchData(object = object, vars = c(group.by, split.by), cells = ind.fov)
  if (is.null(highlight.by)){
    highlight.by <- group.by
  }else{
    df.meta[, highlight.by] <- object@meta.data[ind.fov, highlight.by]
  }
  df.meta$highlight <- ifelse(df.meta[, highlight.by] %in% highlight.groups, "y","n")
  df <- cbind(df, df.meta)

  # get molecule location
  if (!is.null(molecules)){
    ind <- !(molecules %in% rownames(object))
    if (any(ind)){
      message(paste0("Molecules not in the object: ", paste(molecules[ind], collapse = ",")))
      molecules <- molecules[!ind]
    }
    molecules.list <- object@images[[fov]]@molecules$molecules[molecules]
    df.mol <- data.frame()
    for (i in 1:length(molecules.list)){
      df.mol <- rbind(df.mol, data.frame(x=molecules.list[[i]]@coords[,1],
                                         y=molecules.list[[i]]@coords[,2],
                                         mol=names(molecules.list)[i]))
    }
    colnames(df.mol) <- c("x","y",group.by)
    if (!is.null(split.by)){
      splitbys <- unique(df.meta[,split.by])
      df.mol2 <- data.frame()
      for (i in 1:length(splitbys)){
        tmp <- df.mol
        tmp[,split.by] <- splitbys[i]
        df.mol2 <- rbind(df.mol2, tmp)
      }
      df.mol <- df.mol2
    }
    df.mol$highlight <- "m"
    df <- rbind(df, df.mol)
  }

  g <- ggplot() +
    geom_point(data = df, mapping = aes(x=x,y=y,alpha=highlight,size=highlight,col=get(group.by)), shape=16)+
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())  +
    theme(panel.grid = element_blank()) +
    guides(color = guide_legend(title = group.by),
           size = "none",
           alpha = "none") +
    scale_size_manual(values = c('y'=highlight.size,'n'=size,'m'=molecules.size)) +
    scale_alpha_manual(values = c('y'=highlight.alpha,'n'=alpha,'m'=molecules.alpha)) +
    coord_fixed(ratio = 1)

  # colors
  g2 <- ggplot_build(g)
  color_scale <- g2$plot$scales$get_scales("colour")
  cols.default <- color_scale$palette.cache
  names(cols.default) <- sort(unique(df[,group.by]))

  if (!is.null(highlight.cols)){
    if (is.null(cols)){
      cols <- cols.default
    }
    for (i in 1:length(highlight.cols)){
      cols[names(highlight.cols)[i]] <- highlight.cols[i]
    }
  }

  if (!is.null(molecules.cols)){
    if (is.null(cols)){
      cols <- cols.default
    }
    for (i in 1:length(molecules.cols)){
      cols[names(molecules.cols)[i]] <- molecules.cols[i]
    }
  }

  if (!is.null(cols)){
    g <- g + scale_color_manual(values = cols)
  }


  if (!is.null(split.by)){
    g <- g + facet_wrap(reformulate(split.by))
  }

  if (!is.null(crop)){
    if (crop == T){
      crop <- c(min(df$x)-1,max(df$x)+1, min(df$y)-1,max(df$y)+1)
    }else if (length(crop) != 4){
      message("Please provide crop in format: c(min.x, max.x, min.y, max.y)")
    }
    g <- g + coord_fixed(ratio = 1, xlim = crop[1:2], ylim = crop[3:4])
  }
  if (dark.background){
    g <- g + theme(panel.background = element_rect(fill = 'black', color = 'black'))
  }
  if (flip){
    g <- g + coord_flip()
  }
  return(g)
}


#' Transition animation between plots of cells
#'
#' @description Create a transition animation between different cell coordinate plots, e.g. fov, umap, pca
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
#' g <- plotTransition(srt, initial="umap", final="images", group.by="subcelltype", fov=NULL, point.size = 0.3, match.scale.by = "y")
#' animate(g, duration = 7, fps = 10, start_pause = 10, end_pause = 10,
#'         width = 16, height = 8, units = "in", res = 150,
#'         renderer = gifski_renderer())
#' anim_save("transition_animation_umap_to_images.gif")
#' @import ggplot2 dplyr gganimate
#' @export

plotTransition <- function(srt, initial="umap", final="images", group.by=NULL,
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

  g <- ggplot(animation.data, aes(x=x,y=y,col=!!sym(group.by))) +
    geom_point(size=point.size, alpha=point.alpha) +
    coord_fixed(ratio = 1) +
    transition_time(frame) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())

  return(g)
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




