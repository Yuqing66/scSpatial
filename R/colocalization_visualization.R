

#' @title Get df to plot the field
#'
#' @description
#' Calculate the value of the field over the 2d space. The result can be in long or wide format for downstream plotting.
#'
#' @param field A field object created by createField.
#' @param nbins Number of bins on each axis. If only one number, the same number is used for both x and y.
#' @param binlim x, y range for bin calculation, e.g. c(min(x),max(x),min(y),max(y)).
#' @param shape "long" by default, but can change to "wide".
#' @examples
#' binlim.srt <- getSize(srt, fov = "UV109fov1")
#' df.wide <- plotField(field, nbins = 100, binlim = binlim.srt, shape = "wide")
#' plotMt(z=df.wide)
#' df.long <- plotField(field, nbins = 100, binlim = binlim.srt, shape = "long")
#' colnames(df.long) <- c("x","y","z")
#' plotMt.long(df.long)
#' @import magrittr tibble dplyr tidyr
#' @export
#'


plotField <- function(field, nbins, binlim, shape="long"){
  ndim <- ncol(field$coords)
  if (length(nbins) == 1){
    nbins <- rep(nbins, ndim)
  }
  x.cuts <- seq(from = binlim[1], to = binlim[2], length = nbins[1] + 1)
  y.cuts <- seq(from = binlim[3], to = binlim[4], length = nbins[2] + 1)
  coords.plot <- expand.grid(x.cuts, y.cuts)

  values <- getValueInField(field, coords=coords.plot)
  coords.plot$value <- values
  if (shape=="wide"){
    res <- coords.plot %>% pivot_wider(names_from = "Var1", values_from = "value") %>% column_to_rownames("Var2")
  }else{
    res <- coords.plot
  }
  return(res)
}




#' @title Plot 2d raster plot
#'
#' @description
#' Plot 2d raster plot with data.frame in wide format.
#' macOS Preview will do some interpolation of the raster to make the plot in pdf smoother/blurry. Adobe is fine.
#' Raster might have some white line/gap in between.
#'
#' @param x,y Optional. Vectors for coordinates of z;
#' @param z is a count matrix.
#' @param flip Flip the plot to match seurat plot.
#' @param method "tile" by default, but can change to "raster".
#' @examples
#' binlim.srt <- getSize(srt, fov = "UV109fov1")
#' df.wide <- plotField(field, nbins = 100, binlim = binlim.srt, shape = "wide")
#' plotMt(z=df.wide)
#' @import ggplot2 reshape2
#' @export
#'

plotMt <- function(x=NULL,y=NULL,z,flip=F,method="tile"){
  wide <- data.frame(z)
  if (is.null(x)) x <- 1:ncol(z)
  if (is.null(y)) y <- 1:nrow(z)
  colnames(wide) <- x
  wide$y <- y
  long <- melt(wide, id.vars = c("y"),
               variable.name = "x",
               value.name = "z")
  long$x <- as.numeric(as.vector(long$x)) # factor to numeric

  if (method == "tile"){
    g <- ggplot(long,aes(x=x,y=y))+
      geom_tile(aes(fill=z,color=z),size=0,linewidth=0)
    # color and size are used to remove the white lines between tiles
  }else if (method == "raster"){
    g <- ggplot(long,aes(x=x,y=y))+
      geom_raster(aes(fill=z,color=z))
  }
  g <- g + theme(axis.title = element_blank(),axis.text = element_blank(),
                 axis.ticks = element_blank(),legend.title = element_blank())

  if (flip){
    g <- g + coord_flip()
  }
  return(g)
}



#' @title Plot 2d raster plot
#'
#' @description
#' Plot 2d raster plot with data.frame in long format
#'
#' @param df A data.frame with columns x, y, z for coordinates x,y and value z.
#' @param flip Flip the plot to match seurat plot.
#' @examples
#' colnames(df.long) <- c("x","y","z")
#' plotMt.long(df.long)
#' @import ggplot2
#' @export
#'

plotMt.long <- function(df, flip=F){
  g <- ggplot(df,aes(x=x,y=y))+
    geom_tile(aes(fill=z,color=z),linewidth=0.5)+
    # color and size are used to remove the white lines between tiles
    theme(axis.title = element_blank(),axis.text = element_blank(),
          axis.ticks = element_blank(),legend.title = element_blank())
  if (flip){
    g <- g + coord_flip()
  }
  return(g)
}





#' @title Plot the 3d contour plot
#'
#' @description
#' Plot the 3d contour plot.
#'
#' @param df A data.frame with columns x, y for coordinates z for values. Or a wide matrix with values.
#' @param flip Flip the plot to match seurat plot.
#' @param shape "long" by default, but can change to "wide".
#' @examples
#' g <- plotContour(df.long, shape = "long")
#' htmlwidgets::saveWidget(rglwidget(), "field_IFNB1_nbins100_contour.html")
#' browseURL("field_IFNB1_nbins100_contour.html")
#' @import rgl
#' @export
#'


plotContour <- function(df, flip=T, shape="long"){
  if (shape == "long"){
    df <- df %>% pivot_wider(names_from = "x", values_from = "z") %>% column_to_rownames("y")
  }
  if (flip){
    x=as.numeric(rownames(df))
    y=as.numeric(colnames(df))
  }else{
    x=as.numeric(colnames(df))
    y=as.numeric(rownames(df))
  }
  g <- persp3d(x=x,y=y,z=as.matrix(df), contour=T, col=terrain.colors(5))
  return(g)
}






#' @title ImageDimPlot color by provided values.
#' @description
#' Plot points on a 2d space color by values. If Seurat object is provided, plot the location of the cells beneath.
#' @param coords A data.frame with columns x, y.
#' @param values A vector of values for coloring.
#' @param alpha Transparency of the points.
#' @param size Size of the points.
#' @param flip Flip the plot to match seurat plot.
#' @param color Color of the points.
#' @param object Seurat object.
#' @examples
#' ind.moDC <- srt$subcelltype == "monocyte"
#' coords.moDC <- getCoords.cell(srt, ind = ind.moDC, fov = "roi1")
#' values <- getValueInField(field, coords.moDC)
#' ImageDimPlot.values(coords=coords.moDC, values, fov="roi1", srt=srt)
#' @import Seurat
#' @export
#'


ImageDimPlot.values <- function(coords, values, alpha = 1, size = 0.8, flip=F,
                                color=c("grey","purple","red","orange","yellow"),
                                srt=NULL, fov=NULL, alpha.srt=0.1, size.srt=0.1, color.srt="lightblue",dark.background=T){
  g <- ggplot()
  if(!is.null(srt)){
    df <- data.frame(x=srt@images[[fov]]$centroids@coords[,1],
                     y=srt@images[[fov]]$centroids@coords[,2])
    g <- g + geom_point(data = df, mapping = aes(x=x,y=y), alpha=alpha.srt, size=size.srt, color=color.srt)
  }
  coords$value <- values
  g <- g + geom_point(data = coords, mapping = aes(x=x,y=y,col=value), alpha = alpha, size = size) +
    theme(axis.title = element_blank())  +
    theme(panel.grid = element_blank()) +
    theme_classic() +
    coord_fixed(ratio = 1)

  if (!is.null(color)){
    g <- g + scale_color_gradientn(colors = color)
  }
  if (dark.background){
    g <- g + theme(panel.background = element_rect(fill = 'black', color = 'black'))
  }

  if (flip){
    g <- g + coord_flip()
  }
  return(g)
}




#' @title ImageDimPlot.sizeGuide
#'
#' @description
#' Plot circles with different radius on ImageDimPlot.
#' Serve as a size guide to choose sd for plotField colocalization analysis.
#' @param object Seurat object.
#' @param fov Name of the image.
#' @param group.by Column name in the metadata to group color the points.
#' @param size Size of the cells.
#' @param cols A named vector with group categories as names and colors as values.
#' @param alpha Transparency of the cells.
#' @param radii A vector of radii to plot.
#' @param gap Gap between circles.
#' @param position Position of the circles. "bottomleft" by default.
#' @param direction Aligning direction of the circles. "v" by default.
#' @param shape.col Color of the circles.
#' @param shape.alpha Transparency of the circles.
#' @param label Add radius label to the circles.
#' @param label.col Color of the label.
#' @param label.size Size of the label.
#' @param label.gap Gap between the circle and the label.
#' @param dark.background Black background. By default TRUE.
#' @param flip Flip coordinates to match the Seurat output. By default FALSE.
#' @examples
#' ImageDimPlot.sizeGuide(srt, fov = "UV238fov1", group.by = "subCellType.1", gap = 1000, position = "topleft", flip=T)
#' @import Seurat ggforce
#' @export


ImageDimPlot.sizeGuide <- function(object, fov, group.by=NULL, size = 0.2, cols = NULL, alpha = 1,
                       radii = c(100,200,300,400,500,700,1000), gap = NULL,
                       position = "bottomleft", direction = "v",
                       shape.col = "lightblue", shape.alpha = 0.8,
                       label = T, label.col = "red", label.size = 3, label.gap = NULL,
                       dark.background = T, flip = F){
  df <- getCoords.cell(object, fov = fov, meta.cols = group.by)
  if (flip){
    colnames(df) <- c("y","x",group.by)
  }
  binlim <- getSize(object, fov = fov)

  if (is.null(gap)){
    gap <- min(radii)
  }
  if (is.null(label.gap)){
    label.gap <- min(radii)
  }
  # circle center increments (r1 + gap + r2)
  cci <- c(0,radii[1:c(length(radii)-1)]) + gap + radii
  # circle center: cumulated increments
  cc <- data.frame(l=cumsum(cci), w=gap+max(radii), r=radii)
  cc$rl <- cc$l
  cc$rw <- cc$w + cc$r + label.gap

  if (direction %in% c("v", "vertical")){
    colnames(cc) <- c("y","x","r","ry","rx")
    hjust=0
    vjust=0.5
  }else if(direction %in% c("h", "horizontal")){
    colnames(cc) <- c("x","y","r","rx","ry")
    hjust=0.5
    vjust=0
  }else{
    stop("direction must be v or h")
  }


  if(grepl("top", position)){
    cc$y <- binlim[4] - cc$y
    cc$ry <- binlim[4] - cc$ry
    if (direction %in% c("h", "horizontal")){
      vjust <- 1 - vjust
    }
  }
  if(grepl("right", position)){
    cc$x <- binlim[2] - cc$x
    cc$rx <- binlim[2] - cc$rx
    if (direction %in% c("v", "vertical")){
      hjust <- 1 - hjust
    }
  }


  g <- ggplot() +
    geom_point(df, mapping = aes(x=x,y=y,col=.data[[group.by]]), size = size, alpha = alpha) +
    theme(panel.grid = element_blank()) +
    coord_fixed(ratio = 1)
  if (!is.null(cols)) {
    g <- g + scale_color_manual(values = cols)
  }

  g <- g + geom_circle(data = cc, mapping = aes(x0=x,y0=y,r=r), fill=shape.col, alpha=shape.alpha, linewidth=0)

  if (label){
    g <- g + geom_text(data = cc, mapping = aes(x=rx,y=ry,label=r), size = label.size, color = label.col, hjust=hjust, vjust=vjust)
  }
  if (dark.background){
    g <- g + theme(panel.background = element_rect(fill = 'black', color = 'black'))
  }

  return(g)
}




#' @title FeaturePlot.cont.loess
#' @description
#' Plot the gene expression along a continuous variable in meta.data. Smoothed by loess.
#' To check the gene expression pattern according to the colocalization values.
#' @param object Seurat object.
#' @param gene Gene name.
#' @param continuous Column name in the metadata containing the continuous values for x-axis.
#' @param group.by Column name in the metadata to group color the points.
#' @param split.by Column name in the metadata to facet the plot.
#' @param assay Assay to use.
#' @param cols A named vector with group categories as names and colors as values.
#' @param smooth.span Only used when there are fewer than 1,000 observations. Fraction of points used to fit each local regression. Larger numbers make a smoother curve.
#' @param smooth.se Display standard error of the regression. By default FALSE.
#' @param trans Transformation of the x-axis. By default "log1p".
#' @examples
#' srt.inf <- srt[,srt$subCellType.1 %in% c("MC_Inf","MC_Inf2","MC")]
#' FeaturePlot.cont.loess(srt.inf, gene = "LYVE1", continuous = "values_CD4_subCellType2",
#'                        group.by = "Sample", split.by = "Disease", smooth.span = 0.9, trans="log1p")
#' @import ggplot2 Seurat
#' @export

FeaturePlot.cont.loess <- function(object, gene, continuous, group.by=NULL, split.by=NULL, assay = NULL,
                                   cols=NULL, smooth.span=NULL, smooth.se=F, trans="log1p"){
  if (is.null(assay)){
    assay <- object@active.assay
  }
  df <- data.frame(x=object@meta.data[,continuous],
                   gene=object@assays[[assay]]$data[rownames(object) == gene,])
  if (!is.null(group.by) | !is.null(split.by)){
    df.meta <- object@meta.data[,c(group.by,split.by), drop=F]
    df <- cbind(df, df.meta)
  }

  g <- ggplot(df, aes(x=x, y=gene, col=.data[[group.by]])) +
    stat_smooth(se = smooth.se, span = smooth.span)
  if (!is.null(trans)){
    g <- g + scale_x_continuous(transform = trans)
  }
  if (!is.null(cols)) {
    g <- g + scale_color_manual(values = cols)
  }
  if (length(split.by) == 1) {
    g <- g + facet_grid(facets = reformulate(split.by),
                        scales = "free_x")
  }else if (length(split.by) == 2) {
    g <- g + facet_grid(facets = reformulate(split.by[1],
                                             split.by[2]), scales = "free_x")
  }else if (length(split.by) > 2) {
    stop("Parameter split.by needs to be a string with equal or less than two variables.")
  }

  return(g)
}






#' @title FeaturePlot.cont.rollmean
#' @description
#' Plot the gene expression along a continuous variable in meta.data. Smoothed by rolling average sliding across the x-axis.
#' To check the gene expression pattern according to the colocalization values.
#' @param object Seurat object.
#' @param gene Gene name.
#' @param continuous Column name in the metadata containing the continuous values for x-axis.
#' @param group.by Column name in the metadata to group color the points.
#' @param split.by Column name in the metadata to facet the plot.
#' @param assay Assay to use.
#' @param cols A named vector with group categories as names and colors as values.
#' @param window.prop Proportion of data used in the sliding window for smoothing.
#' @param window.n Number of points in the sliding window for smoothing.
#' @param trans Transformation of the x-axis. By default "log1p".
#' @examples
#' srt.inf <- srt[,srt$subCellType.1 %in% c("MC_Inf","MC_Inf2","MC")]
#' FeaturePlot.cont.rollmean(srt.inf, gene = "MMP9", continuous = "values_CD4_subCellType2",
#'                           group.by = "Sample", split.by = "Disease", window.prop = 0.05, trans="log1p")
#' @import ggplot2 Seurat zoo
#' @export

FeaturePlot.cont.rollmean <- function(object, gene, continuous, group.by=NULL, split.by=NULL, assay = NULL,
                                   cols=NULL, window.prop=0.05, window.n=NULL, trans="log1p"){
  if (is.null(assay)){
    assay <- object@active.assay
  }
  df <- data.frame(x=object@meta.data[,continuous],
                   gene=object@assays[[assay]]$data[rownames(object) == gene,])
  if (!is.null(group.by)){
    group.vec <- object@meta.data[,group.by]
  }else{
    group.vec <- NULL
  }
  if (!is.null(split.by)){
    split.vec <- object@meta.data[,split.by]
  }else{
    split.vec <- NULL
  }
  g <- rollMeanPlot(x=df$x, y=df$gene, group.vec=group.vec,split.vec=split.vec,
                    window.prop=window.prop, window.n=window.n, trans=trans)
  if (!is.null(cols)) {
    g <- g + scale_color_manual(values = cols)
  }
  return(g)
}





#' @title plot rolling mean
#' @description
#' Visualize the change of y along a continuous variable x, with moving window for smoothing.
#' The underlying function for FeaturePlot.cont.rollmean.
#' @param x,y A vector of x and y values.
#' @param group Optional. A vector of group for coloring.
#' @param split Optional. A vector of split for facetting.
#' @param window.prop Proportion of data used in the sliding window for smoothing.
#' @param window.n Number of points in the sliding window for smoothing. window.prop will be ignored if window.n is provided.
#' @param trans Transformation of the x-axis.
#' @examples
#' rollMeanPlot(x=srt.inf$values_CD4_subCellType2,
#'              y=srt.inf@assays$seqFISH@layers$data[rownames(srt) == "MMP9", ],
#'              group=srt.inf$Sample,
#'              trans="log1p")
#' @import ggplot2 zoo magrittr dplyr
#' @export

rollMeanPlot <- function(x, y, group.vec=NULL, split.vec=NULL, window.prop=0.05, window.n=NULL, trans=NULL){
  df <- data.frame(x=x, y=y)

  if (!is.null(group.vec) || !is.null(split.vec)){
    df <- addOptionalCols(df, group.vec, split.vec)
    df$tmp <- paste(df$group.vec, df$split.vec, sep="_")

    slices <- unique(df$tmp)
    df$rollingavg <- NA
    for (i in seq_along(slices)){
      ind <- df$tmp == slices[i]
      df$rollingavg[ind] <- rollMean(df$x[ind], df$y[ind], window.prop=window.prop, window.n=window.n)
    }
  }else{
    df$rollingavg <- rollMean(df$x[ind], df$y[ind], window.prop=window.prop, window.n=window.n)
  }

  g <- ggplot(df) +
    geom_line(aes(x = x, y = rollingavg, col = group.vec), size = 1)

  if (!is.null(trans)){
    g <- g + scale_x_continuous(transform = trans)
  }

  if (!is.null(split.vec)){
    g <- g + facet_wrap(~split.vec, scales = "free_y")
  }

  return(g)
}



#' @title Calculate the rolling mean vector
#' @description
#' Calculate the rolling mean of y along x, with moving window for smoothing.
#' Return the rolling mean vector in the same order as input x.
#' @param x,y A vector of x and y values.
#' @param window.prop Proportion of data used in the sliding window for smoothing.
#' @param window.n Number of points in the sliding window for smoothing.
#' @examples
#' x <- 1:1000
#' y <- x^0.5 + rnorm(1000, sd = 100)
#' rollMean(x,y,window.prop=0.1)
#' @import zoo
#' @export

rollMean <- function(x, y, window.prop=0.05, window.n=NULL){
  if(is.null(window.n)){
    window.n <- floor(length(x)*window.prop)
    if (window.n < 2){
      rollingavg <- rep(NA, length(x))
      return(rollingavg)
    }
  }
  padding <- rep(NA, window.n - 1)
  df <- data.frame(x=x, y=y)
  ind <- order(df$x)
  rollingavg <- c(padding, rollmean(df$y[ind], k = window.n))
  rollingavg <- rollingavg[order(ind)]
  return(rollingavg)
}


