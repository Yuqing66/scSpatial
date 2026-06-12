#' @title Calculate minimum distance to reference populations
#'
#' @description
#' For each slide (image) and each reference population, computes the
#' minimum Euclidean distance from every cell to the nearest reference cell.
#' Unlike colocalization which sums Gaussian field values from all reference
#' cells, minimum distance captures proximity to a tissue structure regardless
#' of how many reference cells define it.
#'
#' @param srt Seurat object with spatial images.
#' @param celltype.col Name of the metadata column containing cell type labels.
#' @param ref.celltypes Reference cell type(s). Can be:
#'   \itemize{
#'     \item A character vector: cells of any of these types are treated as one
#'       reference population (e.g., \code{c("CD4", "CD8")} for immune aggregates).
#'     \item A named list of character vectors: each element defines a separate
#'       reference population (e.g., \code{list(epidermis = c("KC_basal"),
#'       immune = c("CD4", "CD8"))}).
#'   }
#' @param images Character vector of image names to process. If \code{NULL},
#'   uses all images from \code{names(srt@@images)}.
#' @param save.plots If \code{TRUE}, saves per-reference spatial value plots
#'   via \code{\link{ImageDimPlot.values}}.
#' @param plot.dir Output directory for plots (used only when
#'   \code{save.plots = TRUE}).
#' @param verbose If \code{TRUE}, prints progress messages.
#'
#' @return A \code{data.frame} with rows matching \code{colnames(srt)} and
#'   one column per reference population. Each value is the Euclidean distance
#'   from that cell to the nearest reference cell. Cells in images where a
#'   reference population has fewer than 2 cells receive \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Single reference population
#' dist.vals <- calcMinDistance(srt, celltype.col = "CellType",
#'                              ref.celltypes = c("KC_basal"))
#'
#' # Multiple reference populations
#' dist.vals <- calcMinDistance(srt, celltype.col = "CellType",
#'                              ref.celltypes = list(
#'                                epidermis = c("KC_basal"),
#'                                immune = c("CD4", "CD8")
#'                              ))
#' }
#' @import Seurat
#' @import ggplot2
#' @export
calcMinDistance <- function(srt, celltype.col, ref.celltypes,
                            images = NULL, save.plots = FALSE, plot.dir = NULL,
                            verbose = TRUE) {
  # Default to all images
  if (is.null(images)) {
    images <- names(srt@images)
  }

  # Normalize ref.celltypes: plain character vector -> named list
  if (is.character(ref.celltypes)) {
    ref.name <- paste(ref.celltypes, collapse = "_")
    ref.celltypes <- setNames(list(ref.celltypes), ref.name)
  }

  # Initialize output data.frame
  df.values <- data.frame(row.names = colnames(srt))

  # Create plot directory if needed
  if (save.plots && !is.null(plot.dir)) {
    dir.create(plot.dir, showWarnings = FALSE, recursive = TRUE)
  }

  for (i.image in seq_along(images)) {
    image.name <- images[i.image]
    if (verbose) message("Processing image: ", image.name,
                         " (", i.image, "/", length(images), ")")

    # Image size for optional plots
    if (save.plots) {
      tmp <- srt@images[[image.name]]@boundaries$centroids@bbox
      w <- (tmp[1, 2] - tmp[1, 1]) / 5000
      h <- (tmp[2, 2] - tmp[2, 1]) / 5000
    }

    # Pre-extract coordinates for all cells in this image ONCE
    coords.all <- getCoords.cell(srt, fov = image.name)
    cells.in.image <- rownames(coords.all)

    # Get all cell types for cells in this image ONCE
    celltypes.all <- srt@meta.data[cells.in.image, celltype.col]

    for (i.ref in seq_along(ref.celltypes)) {
      ref.name <- names(ref.celltypes)[i.ref]
      ref.types <- ref.celltypes[[i.ref]]
      if (verbose) message("  Reference: ", ref.name,
                           " (", i.ref, "/", length(ref.celltypes), ")")

      # Identify reference cells in this image
      cells.ref <- cells.in.image[!is.na(celltypes.all) &
                                    celltypes.all %in% ref.types]

      if (length(cells.ref) < 2) {
        if (verbose) message("    Skipping: fewer than 2 reference cells")
        next
      }

      # Compute minimum Euclidean distance from every cell to nearest ref cell
      coords.ref <- as.matrix(coords.all[cells.ref, 1:2, drop = FALSE])
      coords.query <- as.matrix(coords.all[, 1:2])

      min_dists <- rep(Inf, nrow(coords.query))
      for (j in seq_len(nrow(coords.ref))) {
        d <- sqrt((coords.query[, 1] - coords.ref[j, 1])^2 +
                    (coords.query[, 2] - coords.ref[j, 2])^2)
        min_dists <- pmin(min_dists, d)
      }
      names(min_dists) <- rownames(coords.all)

      # Store in output data.frame
      col.name <- ref.name
      if (!col.name %in% colnames(df.values)) {
        df.values[, col.name] <- NA
      }
      df.values[names(min_dists), col.name] <- min_dists

      # Optional spatial plots
      if (save.plots && !is.null(plot.dir)) {
        p <- ImageDimPlot.values(
          coords = coords.all, values = min_dists, fov = image.name,
          bg.coords = coords.all, dark.background = TRUE, size = 0.1
        ) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank())
        ggsave(file.path(plot.dir,
                         paste0("mindist_", image.name, "_",
                                ref.name, ".png")),
               plot = p, width = w + 0.8, height = h)
      }
    }
  }

  return(df.values)
}
