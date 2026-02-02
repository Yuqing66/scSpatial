#' Extract cell coordinates in Seurat object.
#'
#' @description
#' Extract cell coordinates by name of reduction, or fov name of the slide image.
#'
#' @param object Seurat object.
#' @param name Reduction name as in names(srt@reductions), or "image" to extract image coordinates.
#' @param fov A vector of fovs to extract.
#' @param metadate Columns in meta.data to export together with cell coordinates.
#' @examples
#' coords <- getCoords(srt, name = "umap", fov=NULL, metadata=c("CellType"))
#' coords <- getCoords(srt, name = "images", fov=c("fov1","fov2"), metadata=c("CellType"))
#' @export

#
getCoords <- function(srt, name, fov=NULL, metadata=NULL){
  if (name %in% names(srt@reductions)){
    coords <- as.data.frame(srt@reductions[[name]]@cell.embeddings)
    coords$fov <- name
    coords <- cbind(coords, srt@meta.data[,metadata, drop=F])
  }else if (name == "images"){
    coords <- data.frame()
    if (is.null(fov)){
      fovs <- names(srt@images)
    }else{
      fovs <- intersect(fov, names(srt@images))
      if (length(fovs) < length(fov)) message("not all fov present in the seurat object")
    }
    for (i in 1:length(fovs)){
      fovname <- fovs[i]
      tmp <- as.data.frame(srt@images[[fovname]]@boundaries$centroids@coords)
      rownames(tmp) <- srt@images[[fovname]]@boundaries$centroids@cells
      tmp$fov <- fovname
      tmp <- cbind(tmp, srt@meta.data[match(rownames(tmp), colnames(srt)),metadata, drop=F])
      coords <- rbind(coords, tmp)
    }
  }
  return(coords)
}



#' @title Add optional columns to a data frame.
#' @description
#' Add optional columns to a data frame if not NULL.
#' @param df A data frame.
#' @param ... Vectors to add to the data frame. Format column_name=vector. Use vector name if column_name is not provided.
#' @return A data frame with optional columns added.
#' @examples
#' m1 <- data.frame(a=1:5, b=letters[1:5])
#' c <- 2:6
#' m2 <- addOptionalCols(df=m1, c, d=3:7)
#' @export


addOptionalCols <- function(df, ...){
  dots <- list(...)
  mc <- match.call(expand.dots = FALSE)[["..."]]

  for (i in 1:length(dots)){
    x <- dots[[i]]
    if (!is.null(x)){
      colname <- names(dots)[i]
      if (is.null(colname) || colname=="") colname <- deparse(mc[[i]])

      if (colname %in% colnames(df)){
        message(paste0(colname," already exists in the data frame."))
        next
      }
      if (length(x) != nrow(df)){
        stop(paste0("The length of ",colname," is not equal to the number of rows in the data frame."))
      }
      df[colname] <- x
    }
  }
  return(df)
}








#### get the coordinates (and values) ####

#' @title Get the coordinates of cells in a Seurat object
#'
#' @description
#' Extract the coordinates of cells of an fov.
#' Subset the Seurat object for cell types etc before applying this function.
#'
#' @param object Seurat object. Need to subset srt before applying this function.
#' @param fov The name of the image fov in the Seurat object.
#' @examples
#' coords <- getCoords.cell(srt[,srt$CellType == "Keratinocyte"], fov = "UV109fov1")
#' @import Seurat
#' @export
#'

getCoords.cell <- function(object, fov, meta.cols=NULL){
  fov <- as.character(fov)
  coords <- as.data.frame(object@images[[fov]]@boundaries$centroids@coords)
  rownames(coords) <- object@images[[fov]]@boundaries$centroids@cells
  if (!is.null(meta.cols)){
    meta.df <- object@meta.data[match(rownames(coords), colnames(object)), meta.cols, drop=FALSE]
    coords <- cbind(coords, meta.df)
  }
  return(coords)
}




#' @title Get the coordinates of transcripts in a Seurat object
#'
#' @description
#' Extract the coordinates of transcripts in a fov.
#'
#' @param object Seurat object.
#' @param transcript The name of the transcript.
#' @param fov The name of the image fov in the Seurat object.
#' @examples
#' coords <- getCoords.transcript(srt, transcript="IFNB1", fov="UV109fov1")
#' @import Seurat
#' @export
#'

getCoords.transcript <- function(object, transcript, fov){
  img <- object@images[[fov]]
  coords <- as.data.frame(img@molecules$molecules[[transcript]]@coords)
  return(coords)
}


#' @title Get the coordinates and values of features in a Seurat object
#'
#' @description
#' Get the coordinates of cells in an image fov together with the expression value of a gene.
#'
#' @param object Seurat object. Need to subset srt before applying this function.
#' @param feature The name of the gene.
#' @param fov The name of the image fov in the Seurat object.
#' @param assay The name of the assay.
#' @param slot The slot of the assay.
#' @examples
#' coords <- getCoords.feature(srt, feature="IFNB1", fov="UV109fov1", assay="seqFISH", slot="data")
#' @import Seurat
#' @export
#'

getCoords.feature <- function(object, feature, fov, assay=NULL, slot="data"){
  if (is.null(assay)){
    assay <- srt@active.assay
  }
  img <- object@images[[fov]]
  coords <- as.data.frame(img@boundaries$centroids@coords)
  ind.img <- colnames(object) %in% img@boundaries$centroids@cells
  coords$value <- object@assays[[assay]][slot][feature, ind.img]
  return(coords)
}
# coords <- getCoords.feature(srt, feature="IFNB1", fov="UV109fov1", assay="seqFISH", slot="data")










# plot 2d raster density plot
# x, y are vectors for coordinates of z; z is a count matrix
library(reshape2)
plot.mt <- function(x=NULL,y=NULL,z,flip=F){
  wide <- data.frame(z)
  if (is.null(x)) x <- 1:ncol(z)
  if (is.null(y)) y <- 1:nrow(z)
  colnames(wide) <- x
  wide$y <- y
  long <- melt(wide, id.vars = c("y"),
               variable.name = "x",
               value.name = "z")
  long$x <- as.numeric(as.vector(long$x)) # factor to numeric
  g <- ggplot(long,aes(x=x,y=y))+
    geom_tile(aes(fill=z,color=z),linewidth=0)+
    # color and size are used to remove the white lines between tiles
    theme(axis.title = element_blank(),axis.text = element_blank(),
          axis.ticks = element_blank(),legend.title = element_blank())
  if (flip){
    g <- g + coord_flip()
  }
  return(g)
}
# plot.mt(z=df.wide)

# colnames of df should be c("x","y","z")
plot.mt.long <- function(df, flip=F){
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














#' Rotate coordinate and bounding box data for a spatial component
#'
#' Internal helper used by `rotateCoords()` to update centroid or molecule
#' coordinates and their bounding boxes for 90-degree rotations.
#'
#' @param coords Matrix-like object with x/y columns to rotate.
#' @param another.bbox Two-row bounding box matrix with rows `x` and `y` and columns
#'   `min`/`max`. If provided, will use bbox for the original min max values.
#' @param angle Rotation angle in degrees. Must be 90, 180, or 270.
#'
#' @return A list containing `coords` (rotated coordinates) and `bbox`
#'   (updated bounding box).
#' @noRd
rotate_component_coords <- function(coords, another.bbox = NULL, angle) {
  x <- coords[, 1]
  y <- coords[, 2]
  
  if (!is.null(another.bbox)) {
    x0 <- another.bbox["x", "min"]
    y0 <- another.bbox["y", "min"]
    width <- another.bbox["x", "max"] - x0
    height <- another.bbox["y", "max"] - y0
  } else {
    x0 <- min(x)
    y0 <- min(y)
    width <- max(x) - x0
    height <- max(y) - y0
  }
  
  rel_x <- coords[, 1] - x0
  rel_y <- coords[, 2] - y0
  
  rotated <- switch(
    as.character(angle),
    "90" = cbind(rel_y + y0, width - rel_x + x0),
    "180" = cbind(width - rel_x + x0, height - rel_y + y0),
    "270" = cbind(height - rel_y + y0, rel_x + x0),
    stop("angle must be 90, 180, or 270.")
  )
  colnames(rotated) <- colnames(coords)
  rownames(rotated) <- rownames(coords)
  
  # Update bounding box
  bbox <- rbind(
    x = c(min = min(rotated[, 1], na.rm = TRUE), max = max(rotated[, 1], na.rm = TRUE)),
    y = c(min = min(rotated[, 2], na.rm = TRUE), max = max(rotated[, 2], na.rm = TRUE))
  )
  colnames(bbox) <- c("min", "max")
  rownames(bbox) <- c("x", "y")
  
  list(coords = rotated, bbox = bbox)
}

rotate_spatial_component <- function(component, another.bbox = NULL, angle) {
  if (is.null(component)) {
    return(component)
  }
  rotated <- rotate_component_coords(component@coords, another.bbox, angle)
  component@coords <- rotated$coords
  component@bbox <- rotated$bbox
  component
}

#' Rotate spatial coordinates stored in a Seurat field of view
#'
#' Rotates cell centroids and molecule coordinates for a selected Seurat
#' imaging field of view (`fov`) by an exact 90-degree increment. Both the
#' individual coordinate matrices and their bounding boxes are updated so
#' downstream plotting functions see the rotated geometry.
#'
#' @param srt A Seurat object containing the imaging field of view to rotate.
#' @param fov A single character value naming the entry in `srt@images` to
#'   rotate (for example `"UV109fov1"`).
#' @param angle Desired clockwise rotation in degrees. Must be one of
#'   `0`, `90`, `180`, or `270`.
#'
#' @return The input Seurat object with the requested field of view modified
#'   in-place.
#'
#' @examples
#' seqfish <- rotateCoords(seqfish, "UV109fov1", 90)
#' 
#'
#' @export
rotateCoords <- function(srt, fov, angle) {
  angle <- suppressWarnings(as.numeric(angle))
  if (is.na(angle)) {
    stop("angle must be numeric.")
  }
  angle <- angle %% 360
  valid_angles <- c(0, 90, 180, 270)
  if (!angle %in% valid_angles) {
    stop("angle must be one of 0, 90, 180, or 270 (clockwise).")
  }
  if (angle == 0) {
    return(srt)
  }
  
  if (!fov %in% names(srt@images)) {
    stop("Requested fov not found in the Seurat object.")
  }
  
  fov_image <- srt@images[[fov]]
  
  if (!is.null(fov_image$centroids)) {
    fov_image[["centroids"]] <- rotate_spatial_component(fov_image$centroids, another.bbox = NULL, angle = angle)
  }
  
  if (!is.null(fov_image@molecules)) {
    container <- fov_image@molecules
    for (i in seq_along(container$molecules)) {
      mol <- container$molecules[[i]]
      container[["molecules"]][[i]] <- rotate_spatial_component(mol, another.bbox = fov_image$centroids@bbox, angle = angle)
    }
    fov_image@molecules <- container
  }
  
  srt@images[[fov]] <- fov_image
  srt
}
