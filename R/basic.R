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






