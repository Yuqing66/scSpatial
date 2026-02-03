


###### functions ######

#' Run FindClusters Across Multiple Seeds
#'
#' Executes Seurat's `FindClusters` repeatedly over a set of random seeds and
#' records the cluster identities for each run.
#'
#' @param srt A Seurat object containing the cells of interest.
#' @param seeds Numeric vector of random seeds to iterate over.
#' @param resolution Resolution parameter passed to `Seurat::FindClusters`.
#' @param ... Additional parameters forwarded to `Seurat::FindClusters`.
#'
#' @return A data frame with cells in rows and one column per seed containing
#'   the `seurat_clusters` assignments produced for that run.
#'
#' @examples
#' \dontrun{
#' seed_clusters <- FindClustersAcrossSeeds(srt, seeds = 0:1, resolution = 0.5)
#' head(seed_clusters)
#' }
#' @export
FindClustersAcrossSeeds <- function(srt, seeds = 0:100, resolution = 0.8, ...){
  df.seeds <- data.frame(row.names = colnames(srt))
  for (i in seq_along(seeds)){
    seed <- seeds[[i]]
    srt <- FindClusters(srt, resolution = resolution, random.seed = seed, ...)
    df.seeds[,i] <- srt$seurat_clusters
  }
  colnames(df.seeds) <- paste0("seed", seeds)
  return(df.seeds)
}

#' Plot Cross-Seed Cluster Membership Proportions
#'
#' Visualises the proportion of cells that transition from the clusters in
#' `group1` to the clusters in `group2` as a tile heatmap, highlighting how
#' stable each cluster remains across two clustering runs.
#'
#' @param group1 A vector of baseline cluster assignments for each cell.
#' @param group2 A vector of comparison cluster assignments for the same cells.
#' @param name.1 Optional axis label describing `group1`.
#' @param name.2 Optional axis label describing `group2`.
#'
#' @return A `ggplot` object that can be further customised or saved.
#'
#' @examples
#' BelongPlot(c(1, 1, 2, 2), c(1, 2, 2, 3))
#' @export
BelongPlot <- function(group1, group2, name.1=NULL, name.2=NULL){
  df <- data.frame(x = factor(group1), y = factor(group2))
  df <- df %>%
    group_by(x,y, .drop=F) %>%
    summarise(z = length(x)) %>%
    group_by(x) %>%
    mutate(tot = sum(z)) %>%
    mutate(prop = z/tot)
  
  # df.wide <- df %>% 
  #   group_by(a,b) %>%
  #   summarise(v = length(a)) %>%
  #   pivot_wider(values_from = v, names_from = b) %>%
  #   replace(is.na(.), 0)
  
  g <- ggplot(df,aes(x=x,y=y))+
    geom_tile(aes(fill=prop),linewidth=0)+
    # color and size are used to remove the white lines between tiles
    theme(axis.ticks = element_blank(),legend.title = element_blank(),
          axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
    xlab(name.1) + ylab(name.2)
  return(g)
}

#' Calculate Cluster-Level Consistency Scores
#'
#' Quantifies the stability of each cluster by comparing how often its cells
#' remain together when reclustered. Scores close to 1 indicate strong
#' agreement between the two cluster assignments.
#'
#' @inheritParams BelongPlot
#'
#' @return A named numeric vector of consistency scores for the clusters in
#'   `group1`.
#'
#' @examples
#' consistency(c(1, 1, 2, 2), c(1, 2, 2, 3))
#' @export
consistency <- function(group1, group2){
  df <- data.frame(x = factor(group1), y = factor(group2))
  df2 <- df %>%
    group_by(x,y, .drop=F) %>%
    summarise(z = length(x)) %>%
    group_by(x) %>%
    mutate(tot = sum(z)) %>%
    mutate(prop = z/tot) %>%
    select(-z,-tot) %>%
    pivot_wider(names_from = x, values_from = prop) %>%
    column_to_rownames(var = "y")
  
  tmp <- apply(df2, 2, function(x) sum(x^2))
  s <- (tmp*nrow(df2) - 1) / (nrow(df2) - 1)
  return(s)
}


#' Summarise Consistency Against a Reference Seed
#'
#' Computes pairwise consistency scores between a reference clustering seed and
#' selected comparison seeds. Optionally returns a long-format summary that is
#' ready for plotting.
#'
#' @param df.seeds Data frame of cluster assignments with cells in rows and
#'   seeds in columns.
#' @param reference_seed Column name used as the baseline clustering.
#' @param compare_seeds Optional vector of column names to compare against; by
#'   default all non-reference seeds are used.
#' @param return_long Logical; if `TRUE`, a tidy data frame of scores is
#'   provided in the `long` element of the returned list.
#'
#' @return A list containing the `scores` matrix (comparison seeds x clusters),
#'   the `reference_seed`, the `compare_seeds` used, and, when requested, a
#'   `long` data frame.
#'
#' @examples
#' df <- data.frame(seed0 = c(1, 1, 2, 2), seed1 = c(1, 2, 2, 3))
#' calcSeedConsistency(df, reference_seed = "seed0", return_long = TRUE)
#' @export
calcSeedConsistency <- function(df.seeds, reference_seed,
                                compare_seeds = NULL, return_long = FALSE){
  if (is.null(compare_seeds)) compare_seeds <- setdiff(colnames(df.seeds), reference_seed)
  ref_clusters <- sort(unique(df.seeds[, reference_seed]))
  scores <- matrix(NA_real_, nrow = length(compare_seeds), ncol = length(ref_clusters),
                   dimnames = list(compare_seeds, ref_clusters))
  for (seed in compare_seeds){
    cons <- consistency(df.seeds[, reference_seed], df.seeds[, seed])
    matched <- intersect(names(cons), ref_clusters)
    scores[seed, matched] <- cons[matched]
  }
  summary <- list(scores = scores, reference_seed = reference_seed,
                  compare_seeds = compare_seeds)
  if (return_long){
    summary$long <- data.frame(
      seed = rep(rownames(scores), each = ncol(scores)),
      cluster = rep(colnames(scores), times = nrow(scores)),
      score = as.vector(scores),
      row.names = NULL
    )
  }
  summary
}


#' Plot Consistency Score Distributions
#'
#' Creates a ridge-density plot that summarises the distribution of consistency
#' scores for each cluster.
#'
#' @param df Data frame containing at least cluster identifiers and scores.
#' @param cluster_col Name of the column with cluster identifiers.
#' @param score_col Name of the column with consistency scores.
#' @param xlim Numeric vector of length two giving x-axis limits (optional).
#' @param bins Number of histogram bins used in the ridge density calculation.
#' @param alpha Transparency applied to the ridge fill.
#' @param fill_palette Optional named vector of colours to use for clusters.
#'
#' @return A `ggplot` object showing per-cluster consistency score
#'   distributions.
#'
#' @examples
#' df <- data.frame(cluster = c("1", "1", "2", "2"), score = c(0.8, 0.9, 0.3, 0.5))
#' plotConsistencyRidges(df)
#' @export
plotConsistencyRidges <- function(df, cluster_col = "cluster",
                                  score_col = "score", xlim = NULL,
                                  bins = 70, alpha = 0.5, fill_palette = NULL){
  plot_df <- df
  clusters <- as.character(plot_df[[cluster_col]])
  cluster_numeric <- suppressWarnings(as.numeric(clusters))
  if (all(!is.na(cluster_numeric))){
    clusters <- factor(clusters, levels = as.character(sort(unique(cluster_numeric))))
  } else {
    clusters <- factor(clusters)
  }
  plot_df$.__cluster__ <- clusters
  plot_df$.__score__ <- as.numeric(plot_df[[score_col]])
  g <- ggplot(plot_df) +
    geom_density_ridges(aes(x = .__score__, y = .__cluster__, fill = .__cluster__),
                        stat = "binline", bins = bins, alpha = alpha) +
    labs(x = "Consistency score", y = "Cluster", fill = "Cluster") +
    NoLegend()
  if (!is.null(fill_palette)){
    g <- g + scale_fill_manual(values = fill_palette)
  }
  if (!is.null(xlim)) g <- g + xlim(xlim)
  g
}


#' Cell-Level Pairwise Retention Proportions
#'
#' Calculates, for each cell, the proportion of its original cluster mates in
#' `group1` that remain grouped together in `group2`. This per-cell metric can
#' be used to measure how cohesive a cluster stays across runs.
#'
#' @inheritParams BelongPlot
#'
#' @return A numeric vector giving, for each cell, the proportion of its
#'   `group1` cluster mates that stay together in `group2`.
#'
#' @examples
#' calcProp(c(1, 1, 2, 2), c(1, 2, 2, 3))
#' @export
calcProp <- function(group1,group2){
  df <- data.frame(x = factor(group1), y = factor(group2))
  df2 <- df %>%
    group_by(x,y, .drop=F) %>%
    summarise(z = length(x)) %>%
    group_by(x) %>%
    mutate(tot = sum(z)) %>%
    mutate(prop = z/tot)
  
  for (i in 1:nrow(df2)){
    ind <- df$x==df2$x[i] & df$y==df2$y[i]
    df$prop[ind] <- df2$prop[i]
  }
  return(df$prop)
}

#' Average Cell Stability Across Seeds
#'
#' Aggregates per-cell retention proportions across multiple clustering runs to
#' produce a stability score that summarises how reliably each cell stays with
#' its original cluster across seeds.
#'
#' @param df.seeds A data frame where each column contains cluster assignments
#'   for a seed and each row corresponds to a cell.
#' @param i.seed1 The column name in `df.seeds` to use as the reference seed.
#' @param i.seed2s Optional vector of column names to compare against. Defaults
#'   to all seeds except `i.seed1`.
#'
#' @return A numeric vector of per-cell stability scores between 0 and 1.
#'
#' @examples
#' df <- data.frame(seed0 = c(1, 1, 2, 2), seed1 = c(1, 2, 2, 3))
#' calcCellScore(df, i.seed1 = "seed0")
#' @export
calcCellScore <- function(df.seeds, i.seed1, i.seed2s=NULL){
  if(is.null(i.seed2s)) i.seed2s <- setdiff(colnames(df.seeds), i.seed1)
  
  props <- data.frame(row.names = rownames(df.seeds))
  for(i in 1:length(i.seed2s)){
    i.seed2 <- i.seed2s[i]
    props[,i] <- calcProp(group1 = df.seeds[,i.seed1], group2 = df.seeds[,i.seed2])
  }
  cellscore <- apply(props, 1, mean)
  return(cellscore)
}

# hist(cellscore, breaks = 100)


# pick one column and generate pair wise comparison against others

###### Clustering robustness tutorial ######

if (FALSE) {
  # Step 0: Define parameters and create output directories
  seed_resolution <- 0.8
  seed_range <- 0:100
  output_root <- "robustness"
  heatmap_dir <- file.path(output_root, "proportionHeatmap")
  dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

  # Optional: recompute neighbours before clustering
  # srt <- FindNeighbors(srt, dims = 1:30, reduction = "integrated.rpca")

  # Step 1: Run clustering across multiple random seeds
  seed_clusters <- FindClustersAcrossSeeds(srt, seeds = seed_range, resolution = seed_resolution)
  write.table(seed_clusters, file = "seuratclusters_pc30_res80_seed0to100.txt")

  # Step 2: Inspect an example pair of seeds
  example_seed1 <- "seed0"
  example_seed2 <- "seed1"
  example_plot <- BelongPlot(seed_clusters[, example_seed1], seed_clusters[, example_seed2],
                             name.1 = example_seed1, name.2 = example_seed2) +
    scale_fill_gradientn(colours = c("grey", "blue", "red", "yellow"))
  print(example_plot)

  # Step 3: Summarise consistency across seeds relative to the reference
  reference_seed <- "seed0"
  consistency_summary <- calcSeedConsistency(seed_clusters, reference_seed = reference_seed,
                                             return_long = TRUE)
  seed_consistency_matrix <- consistency_summary$scores
  consistency_long <- consistency_summary$long

  # Step 4: Visualise per-cluster stability
  consistency_plot <- plotConsistencyRidges(consistency_long, xlim = c(0, 1))
  print(consistency_plot)
  ggsave("robustness/hist_score_seed0_vs_1to100.png", plot = consistency_plot,
         width = 5, height = 8)

  # Step 5: Highlight unstable cells on embeddings
  srt$cellscore <- calcCellScore(seed_clusters, i.seed1 = reference_seed, i.seed2s = NULL)
  cellscore_plot <- FeaturePlot(srt, "cellscore") +
    scale_color_gradientn(colours = c("yellow", "red", "blue", "grey"))
  print(cellscore_plot)
  ggsave("robustness/cellscore_seed0_vs_1to100.png", plot = cellscore_plot,
         width = 7, height = 6)

  cellscore_faceted <- FeaturePlot(srt, "cellscore", split.by = "seurat_clusters",
                                   ncol = 4, cols = c("yellow", "red", "blue", "grey"))
  ggsave("robustness/cellscore_seed0_vs_1to100_split.png", plot = cellscore_faceted,
         width = 73, height = 6, limitsize = FALSE)

  # Optional: Save a reference faceted cluster map
  cluster_reference_plot <- DimPlot.grey(srt, group.by = "seurat_clusters",
                                         split.by = "seurat_clusters")
  ggsave("robustness/seuratclusters_dim50_res80_seed0_split.png",
         plot = cluster_reference_plot, width = 17, height = 12)
}
