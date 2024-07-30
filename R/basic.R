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
