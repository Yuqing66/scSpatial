#' Calculate colocalization values for all cell types across all slides
#'
#' @description
#' For each slide (image) and each reference cell type, creates a Gaussian
#' density field from the reference cell positions and evaluates it at every
#' cell position. Returns a data.frame of per-cell colocalization values that
#' can be fed into \code{\link{calcColocalizationHeatmap}}.
#'
#' @param srt Seurat object with spatial images.
#' @param celltype.col Name of the metadata column containing cell type labels.
#' @param sd Standard deviation of the Gaussian field (controls the spatial
#'   smoothing radius). Larger values produce smoother fields.
#' @param scale.factor Scale factor for field values (default 10000).
#' @param images Character vector of image names to process. If \code{NULL},
#'   uses all images from \code{names(srt@@images)}.
#' @param save.plots If \code{TRUE}, saves per-reference-celltype spatial value
#'   plots via \code{\link{ImageDimPlot.values}}.
#' @param plot.dir Output directory for plots (used only when
#'   \code{save.plots = TRUE}).
#' @param verbose If \code{TRUE}, prints progress messages.
#'
#' @return A \code{data.frame} with rows matching \code{colnames(srt)} and
#'   columns matching reference cell type names. Each value represents how
#'   strongly a cell colocalizes with the spatial distribution of the reference
#'   cell type. Cells in images where a reference type has fewer than 2 cells
#'   receive \code{NA}.
#'
#' @examples
#' \dontrun{
#' coloc <- calcColocalization(srt, celltype.col = "CellType", sd = 500)
#' }
#' @import Seurat
#' @import ggplot2
#' @export
calcColocalization <- function(srt, celltype.col, sd = 500, scale.factor = 10000,
                               images = NULL, save.plots = FALSE, plot.dir = NULL,
                               verbose = TRUE) {
  if (is.null(images)) {
    images <- names(srt@images)
  }

  celltypes.ref <- levels(factor(srt@meta.data[, celltype.col]))
  cutoff <- findDistance.euclidean(0.0001, sd = sd)

  df.values <- data.frame(row.names = colnames(srt))

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

    for (i.ref in seq_along(celltypes.ref)) {
      celltype.ref <- celltypes.ref[i.ref]
      if (verbose) message("  Reference: ", celltype.ref,
                           " (", i.ref, "/", length(celltypes.ref), ")")

      # Identify reference cells directly from pre-extracted celltypes
      cells.ref.in.image <- cells.in.image[!is.na(celltypes.all) & celltypes.all == celltype.ref]

      if (length(cells.ref.in.image) < 2) {
        if (verbose) message("    Skipping: fewer than 2 reference cells")
        next
      }

      # Create Gaussian field from reference cells using subset of coords
      coords.ref <- coords.all[cells.ref.in.image, , drop = FALSE]
      field <- createField(coords.ref, sd = sd, h = NULL, d.cutoff = cutoff,
                           scale.factor = scale.factor)

      # Evaluate field at all cell positions
      values <- getValueInField(field, coords.all)
      names(values) <- rownames(coords.all)

      # Store
      col.name <- celltype.ref
      if (!col.name %in% colnames(df.values)) {
        df.values[, col.name] <- NA
      }
      df.values[names(values), col.name] <- values

      # Optional spatial plots
      if (save.plots && !is.null(plot.dir)) {
        p <- ImageDimPlot.values(
          coords = coords.all, values = values, fov = image.name,
          bg.coords = coords.all, dark.background = TRUE, size = 0.1
        ) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank())
        ggsave(file.path(plot.dir,
                         paste0("coloc_", image.name, "_",
                                celltype.ref, "_sd", sd, ".png")),
               plot = p, width = w + 0.8, height = h)
      }
    }
  }

  return(df.values)
}




#' Compute per-slide colocalization heatmap matrices
#'
#' @description
#' Takes per-cell colocalization values from \code{\link{calcColocalization}},
#' optionally merges and removes cell types, aggregates values by cell type
#' mean, and scales so that self-colocalization = 100 and the overall mean = 0.
#'
#' @param srt Seurat object.
#' @param coloc.values A \code{data.frame} returned by
#'   \code{\link{calcColocalization}}.
#' @param celltype.col Name of the metadata column containing cell type labels.
#' @param merge.list Optional named list for merging cell types.
#'   Names are the new (merged) labels; values are character vectors of original
#'   labels to merge. For example,
#'   \code{list(CD4 = c("CD4", "CD4_CCL17", "CD4_CXCL13"))}.
#' @param rm.celltypes Character vector of cell types to remove from
#'   \strong{all} slides.
#' @param rm.celltypes.per.image Named list of cell types to remove per image.
#'   Names must match \code{names(srt@@images)}.
#'   For example,
#'   \code{list(UV109fov1 = c("B_cell", "KC_bulge"))}.
#' @param images Character vector of image names to process. If \code{NULL},
#'   uses \code{names(srt@@images)}.
#'
#' @return A named list (one element per image). Each element is itself a list
#'   with:
#'   \describe{
#'     \item{\code{heatmap}}{Aggregated mean colocalization matrix
#'       (cell type \eqn{\times} cell type), including an \code{"all"} row.}
#'     \item{\code{heatmap.scaled}}{Scaled version where self-colocalization
#'       = 100 and overall mean = 0. Rows reordered to match columns.}
#'   }
#'
#' @examples
#' \dontrun{
#' heatmap.list <- calcColocalizationHeatmap(
#'   srt, coloc.values,
#'   celltype.col = "CellType",
#'   merge.list   = list(CD4 = c("CD4", "CD4_CCL17")),
#'   rm.celltypes = c("Neutrophil"),
#'   rm.celltypes.per.image = list(UV109fov1 = c("B_cell"))
#' )
#' }
#' @import dplyr
#' @importFrom tibble column_to_rownames
#' @export
calcColocalizationHeatmap <- function(srt, coloc.values, celltype.col,
                                      merge.list = NULL,
                                      rm.celltypes = NULL,
                                      rm.celltypes.per.image = NULL,
                                      images = NULL) {
  if (is.null(images)) {
    images <- names(srt@images)
  }

  result <- list()

  for (img.name in images) {
    # Identify cells in this image
    cells.in.image <- srt@images[[img.name]]@boundaries$centroids@cells
    cells.in.image <- intersect(cells.in.image, rownames(coloc.values))

    df <- coloc.values[cells.in.image, , drop = FALSE]
    df.celltype <- df.celltype.org <-
      as.vector(srt@meta.data[cells.in.image, celltype.col])

    # ---- Merge cell types ----
    if (!is.null(merge.list)) {
      for (j in seq_along(merge.list)) {
        celltype.to     <- names(merge.list)[j]
        celltypes.from  <- merge.list[[j]]

        # Update labels
        df.celltype[df.celltype %in% celltypes.from] <- celltype.to

        # Merge columns (sum values from merged types)
        celltypes.from.ind <- which(colnames(df) %in% celltypes.from)
        if (length(celltypes.from.ind) > 1) {
          merged.values <- rowSums(df[, celltypes.from.ind], na.rm = TRUE)
        } else if (length(celltypes.from.ind) == 1) {
          merged.values <- df[, celltypes.from.ind]
        } else {
          next
        }

        df[, celltypes.from.ind] <- NULL
        df[, celltype.to] <- merged.values
      }
    }

    # ---- Remove cell types ----
    rm.list.combined <- c(rm.celltypes, rm.celltypes.per.image[[img.name]])
    if (length(rm.list.combined) > 0) {
      ind.keep <- !((df.celltype %in% rm.list.combined) |
                       (df.celltype.org %in% rm.list.combined))
      df <- df[ind.keep, !(colnames(df) %in% rm.list.combined), drop = FALSE]
      df.celltype     <- df.celltype[ind.keep]
      df.celltype.org <- df.celltype.org[ind.keep]
    }

    # ---- Remove cell types with NA values ----
    if (nrow(df) > 0 && ncol(df) > 0) {
      rm.list.na <- c(
        colnames(df)[is.na(df[1, ])],              # ref (columns) with NA
        unique(df.celltype.org[is.na(df[, 1])])     # query (rows) with NA
      )
      if (length(rm.list.na) > 0) {
        ind.keep <- !((df.celltype %in% rm.list.na) |
                         (df.celltype.org %in% rm.list.na))
        df <- df[ind.keep, !(colnames(df) %in% rm.list.na), drop = FALSE]
        df.celltype <- df.celltype[ind.keep]
      }
    }

    if (nrow(df) == 0 || ncol(df) == 0) {
      warning("No data remaining for image: ", img.name, ". Skipping.")
      next
    }

    # ---- Aggregate: mean per cell type ----
    df.heatmap <- df %>%
      mutate(celltype = df.celltype) %>%
      group_by(.data$celltype) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("celltype")

    # Add overall mean row
    query.all <- apply(df, 2, mean)
    df.heatmap <- rbind(df.heatmap, all = query.all)

    # ---- Scale: self = 100, overall mean = 0 ----
    df.heatmap.scaled <- df.heatmap
    ind.all <- grep("^all$", rownames(df.heatmap))
    for (j in seq_len(ncol(df.heatmap))) {
      celltype.ref.name <- colnames(df.heatmap)[j]
      ind.ref <- which(rownames(df.heatmap) == celltype.ref.name)
      if (length(ind.ref) == 0) next

      denom <- df.heatmap[ind.ref, j] - df.heatmap[ind.all, j]
      if (abs(denom) < .Machine$double.eps) {
        df.heatmap.scaled[, j] <- 0
      } else {
        sf <- 100 / denom
        df.heatmap.scaled[, j] <- (df.heatmap[, j] - df.heatmap[ind.all, j]) * sf
      }
    }

    # Reorder rows to match column names
    df.heatmap.scaled <- df.heatmap.scaled[colnames(df.heatmap.scaled), ]

    result[[img.name]] <- list(
      heatmap        = df.heatmap,
      heatmap.scaled = df.heatmap.scaled
    )
  }

  return(result)
}




#' Plot colocalization heatmaps
#'
#' @description
#' Plot colocalization heatmaps from the output of
#' \code{\link{calcColocalizationHeatmap}}. Supports three modes:
#' \itemize{
#'   \item \code{"average"} — average the scaled values across all slides.
#'   \item \code{"per.slide"} — plot each slide's scaled heatmap individually.
#'   \item \code{"difference"} — for each slide, subtract the cross-slide
#'     average to highlight slide-specific deviations.
#' }
#'
#' @param heatmap.list Output from \code{\link{calcColocalizationHeatmap}}.
#' @param type One of \code{"average"} (default), \code{"per.slide"}, or
#'   \code{"difference"}.
#' @param cluster.by Clustering mode controlling how cell types are ordered on
#'   each axis. One of:
#'   \itemize{
#'     \item \code{"reference"} (default) — cluster reference cell types
#'       (columns) by hierarchical clustering, then reorder query rows to match
#'       the same order. Recommended for comparing colocalization profiles across
#'       reference types. Applies consistently to all \code{type} outputs.
#'     \item \code{"query"} — cluster query cell types (rows) by hierarchical
#'       clustering, then reorder reference columns to the same order. Useful
#'       when you want to group query types with similar colocalization patterns.
#'       Applies consistently to all \code{type} outputs.
#'     \item \code{"both"} — cluster rows and columns independently.
#'     \item \code{"none"} — no clustering; preserve the input order.
#'     \item \code{"average.order"} — only for \code{type = "difference"}.
#'       Computes the reference-based clustering order from the cross-slide
#'       average heatmap, then applies that fixed order to both columns and rows
#'       of each difference heatmap. Reproduces the original workflow of ordering per-slide
#'       deviations by the average's structure.
#'   }
#' @param file PDF file path. If \code{NULL}, returns the Heatmap object(s)
#'   without saving. For \code{"per.slide"} and \code{"difference"}, the sample
#'   name is appended before the \code{.pdf} extension.
#' @param width PDF width in inches (default 10).
#' @param height PDF height in inches (default 9).
#'
#' @return For \code{type = "average"}: a \code{ComplexHeatmap::Heatmap} object
#'   (returned invisibly when saved to file).\cr
#'   For \code{type = "per.slide"} or \code{"difference"}: a named list of
#'   Heatmap objects.
#'
#' @examples
#' \dontrun{
#' # Average heatmap
#' plotColocalizationHeatmap(heatmap.list, type = "average")
#'
#' # Per-slide heatmaps saved to PDF
#' plotColocalizationHeatmap(heatmap.list, type = "per.slide",
#'                           file = "heatmap_per_slide.pdf")
#'
#' # Difference (each slide minus the average)
#' plotColocalizationHeatmap(heatmap.list, type = "difference")
#' }
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer separate pivot_wider
#' @importFrom rlang .data
#' @importFrom grDevices pdf dev.off
#' @export
plotColocalizationHeatmap <- function(heatmap.list,
                                     type = "average",
                                     cluster.by = "reference",
                                     file = NULL,
                                     width = 10, height = 9) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Install with:\n",
         "  BiocManager::install('ComplexHeatmap')")
  }

  type       <- match.arg(type,       c("average", "per.slide", "difference"))
  cluster.by <- match.arg(cluster.by, c("reference", "query", "both", "none",
                                         "average.order"))

  if (cluster.by == "average.order" && type != "difference") {
    stop("cluster.by = 'average.order' is only supported for type = 'difference'.")
  }

  # rows = query cell types; columns = reference cell types
  cluster_rows    <- cluster.by %in% c("both", "query")
  cluster_columns <- cluster.by %in% c("both", "reference")

  # ================================================================
  # Per-slide heatmaps
  # ================================================================
  if (type == "per.slide") {
    hmaps <- list()
    for (sample.name in names(heatmap.list)) {
      mat <- as.matrix(heatmap.list[[sample.name]]$heatmap.scaled)
      
      if (cluster.by %in% c("reference", "query")) {
        mat.sq <- mat[rownames(mat) != "all", , drop = FALSE]
        if (cluster.by == "reference") {
          hc            <- hclust(dist(t(mat.sq)))       # cluster columns (reference)
          ordered.names <- colnames(mat.sq)[hc$order]
          if ("all" %in% rownames(mat)) {
            mat.ordered <- rbind(mat.sq[ordered.names, ordered.names, drop = FALSE], all = mat["all", ordered.names])
          } else {
            mat.ordered <- mat.sq[ordered.names, ordered.names, drop = FALSE]
          }
        } else {                                         # "query"
          hc            <- hclust(dist(mat.sq))          # cluster rows (query)
          ordered.names <- rownames(mat.sq)[hc$order]
          mat.ordered <- mat.sq[ordered.names, ordered.names, drop = FALSE]
        }
        
        hmap <- ComplexHeatmap::Heatmap(
          mat.ordered,
          cluster_rows    = FALSE,
          cluster_columns = FALSE,
          show_row_dend   = FALSE,
          show_column_dend = FALSE,
          column_title    = sample.name
        )
      } else {
        hmap <- ComplexHeatmap::Heatmap(
          mat,
          cluster_rows    = cluster_rows,
          cluster_columns = cluster_columns,
          show_row_dend   = FALSE,
          show_column_dend = FALSE,
          column_title    = sample.name
        )
      }
      
      hmaps[[sample.name]] <- hmap

      if (!is.null(file)) {
        f <- sub("\\.pdf$", paste0("_", sample.name, ".pdf"), file)
        pdf(f, width = width, height = height)
        ComplexHeatmap::draw(hmap)
        dev.off()
      }
    }
    
    if (!is.null(file)) {
      return(invisible(hmaps))
    } else {
      return(hmaps)
    }
  }

  # ================================================================
  # Average / Difference — both need the cross-slide average
  # ================================================================
  avg.result <- .averageHeatmaps(heatmap.list)

  if (type == "average") {
    # avg.wide is square (cell types x cell types, no "all" row).
    # Keep mat.sq separate so "all" is never included in clustering.
    mat.sq <- as.matrix(avg.result$avg.wide)

    if (cluster.by %in% c("reference", "query")) {
      # Compute clustering explicitly with hclust/dist.
      # This avoids two problems with the draw()-then-extract approach:
      #   1. draw() renders an intermediate plot as a side effect.
      #   2. Including the "all" row (all zeros) distorts the dendrogram.
      if (cluster.by == "reference") {
        hc            <- hclust(dist(t(mat.sq)))       # cluster columns (reference)
        ordered.names <- colnames(mat.sq)[hc$order]
        # Mirror column order onto rows; append "all" row at the bottom.
        mat.ordered <- rbind(mat.sq[ordered.names, ordered.names, drop = FALSE],
                             all = 0)
      } else {                                         # "query"
        hc            <- hclust(dist(mat.sq))          # cluster rows (query)
        ordered.names <- rownames(mat.sq)[hc$order]
        # Mirror row order onto columns; omit "all" (matches original behavior).
        mat.ordered <- mat.sq[ordered.names, ordered.names, drop = FALSE]
      }

      hmap <- ComplexHeatmap::Heatmap(
        mat.ordered,
        cluster_rows     = FALSE,
        cluster_columns  = FALSE,
        show_row_dend    = FALSE,
        show_column_dend = FALSE,
        column_title     = "Average across slides"
      )

    } else {
      # "both" or "none": let ComplexHeatmap handle clustering directly.
      mat <- rbind(mat.sq, all = 0)
      hmap <- ComplexHeatmap::Heatmap(
        mat,
        cluster_rows     = cluster_rows,
        cluster_columns  = cluster_columns,
        show_row_dend    = FALSE,
        show_column_dend = FALSE,
        column_title     = "Average across slides"
      )
    }

    if (!is.null(file)) {
      pdf(file, width = width, height = height)
      ComplexHeatmap::draw(hmap)
      dev.off()
      return(invisible(hmap))
    } else {
      return(hmap)
    }
  }

  # ================================================================
  # Difference: per-slide minus average
  # ================================================================
  if (type == "difference") {
    # Compute clustering order from average if requested
    avg.order <- NULL
    if (cluster.by == "average.order") {
      mat.sq <- as.matrix(avg.result$avg.wide)
      hc <- hclust(dist(t(mat.sq)))       # cluster columns (reference)
      avg.order <- colnames(mat.sq)[hc$order]
    }

    hmaps <- list()
    for (sample.name in names(heatmap.list)) {
      diff.mat <- .computeDifference(avg.result$avg.long, sample.name)
      if (is.null(diff.mat)) next

      if (!is.null(avg.order)) {
        # Reorder to match average heatmap's clustering order
        common <- intersect(avg.order, colnames(diff.mat))
        diff.mat <- diff.mat[common, common, drop = FALSE]
        
        hmap <- ComplexHeatmap::Heatmap(
          as.matrix(diff.mat),
          cluster_rows     = FALSE,
          cluster_columns  = FALSE,
          show_row_dend    = FALSE,
          show_column_dend = FALSE,
          column_title     = paste0(sample.name, " vs. average")
        )
      } else if (cluster.by %in% c("reference", "query")) {
        mat.sq <- as.matrix(diff.mat)
        if (cluster.by == "reference") {
          hc            <- hclust(dist(t(mat.sq)))
          ordered.names <- colnames(mat.sq)[hc$order]
          diff.mat <- mat.sq[ordered.names, ordered.names, drop = FALSE]
        } else {
          hc            <- hclust(dist(mat.sq))
          ordered.names <- rownames(mat.sq)[hc$order]
          diff.mat <- mat.sq[ordered.names, ordered.names, drop = FALSE]
        }
        
        hmap <- ComplexHeatmap::Heatmap(
          as.matrix(diff.mat),
          cluster_rows     = FALSE,
          cluster_columns  = FALSE,
          show_row_dend    = FALSE,
          show_column_dend = FALSE,
          column_title     = paste0(sample.name, " vs. average")
        )
      } else {
        hmap <- ComplexHeatmap::Heatmap(
          as.matrix(diff.mat),
          cluster_rows     = cluster_rows,
          cluster_columns  = cluster_columns,
          show_row_dend    = FALSE,
          show_column_dend = FALSE,
          column_title     = paste0(sample.name, " vs. average")
        )
      }
      
      hmaps[[sample.name]] <- hmap

      if (!is.null(file)) {
        f <- sub("\\.pdf$", paste0("_", sample.name, "_diff.pdf"), file)
        pdf(f, width = width, height = height)
        ComplexHeatmap::draw(hmap)
        dev.off()
      }
    }
    
    if (!is.null(file)) {
      return(invisible(hmaps))
    } else {
      return(hmaps)
    }
  }
}



# =====================================================================
# Internal helpers
# =====================================================================

#' Average scaled heatmap matrices across slides
#' @noRd
#' @keywords internal
.averageHeatmaps <- function(heatmap.list) {
  all.celltypes <- unique(unlist(
    lapply(heatmap.list, function(x) colnames(x$heatmap.scaled))
  ))

  df.merge <- expand.grid(ref = all.celltypes, query = all.celltypes,
                           stringsAsFactors = FALSE)
  rownames(df.merge) <- paste0(df.merge$ref, ":", df.merge$query)

  for (sample.name in names(heatmap.list)) {
    df <- heatmap.list[[sample.name]]$heatmap.scaled
    df.long <- df %>%
      tibble::rownames_to_column("query") %>%
      tidyr::pivot_longer(-"query", names_to = "ref", values_to = "value") %>%
      mutate(ref_query = paste0(.data$ref, ":", .data$query))

    df.merge[match(df.long$ref_query, rownames(df.merge)), sample.name] <-
      df.long$value
  }

  df.merge <- df.merge[, !(colnames(df.merge) %in% c("ref", "query")),
                        drop = FALSE]
  df.merge$avg <- rowMeans(df.merge[, names(heatmap.list), drop = FALSE],
                            na.rm = TRUE)

  df.avg.wide <- df.merge %>%
    select("avg") %>%
    tibble::rownames_to_column("tmp") %>%
    tidyr::separate("tmp", c("ref", "query"), sep = ":") %>%
    tidyr::pivot_wider(names_from = "ref", values_from = "avg") %>%
    tibble::column_to_rownames("query")

  list(avg.long = df.merge, avg.wide = df.avg.wide)
}


#' Compute per-slide minus average heatmap
#' @noRd
#' @keywords internal
.computeDifference <- function(avg.long, sample.name) {
  if (!sample.name %in% colnames(avg.long)) {
    warning("Sample '", sample.name, "' not found in averaged data.")
    return(NULL)
  }

  df.diff <- avg.long %>%
    mutate(value = .data[[sample.name]] - .data$avg) %>%
    select("value") %>%
    replace(is.na(.), 0) %>%
    tibble::rownames_to_column("tmp") %>%
    tidyr::separate("tmp", c("ref", "query"), sep = ":") %>%
    tidyr::pivot_wider(names_from = "ref", values_from = "value") %>%
    tibble::column_to_rownames("query")

  df.diff
}
