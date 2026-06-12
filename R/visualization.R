#' plot the inclusion relationships of the cell type annotations with different granularity
#'
#' @description Cell type annotation can happen at different granularity, e.g. main cell type, sub cell type, and sub sub cell type.
#' To visualize the inclusion relationships of these annotations, a tree plot is generated.
#'
#' @param object Seurat object.
#' @param groups a vector of column names in the metadata.
#' @param width.node Width of Sankey nodes.
#' @param size.text Font size for node labels.
#' @param height Height of the output widget.
#' @param width Width of the output widget.
#' @param save Logical; if TRUE, save the widget as HTML and PDF.
#' @param filename Optional filename for the saved widget.
#' @examples
#' \dontrun{
#' groups <- c("CellType.3","CellType.2","CellType.1","subCellType.2","subCellType.1")
#' p <- SankeyPlot(srt, groups)
#' }
#' @import dplyr
#' @importFrom networkD3 sankeyNetwork saveNetwork
#' @importFrom webshot webshot
#' @export
#'
#'
#'


# webshot::install_phantomjs()



SankeyPlot <- function(object, groups, width.node = 20, size.text = 12, height = NULL, width = NULL,
                       save = TRUE, filename = NULL){
  # Create a data frame with the cell type annotations
  cell_data <- object@meta.data[,groups]
  cell_data <- cell_data %>% arrange(!!!syms(groups))
  
  # Create nodes data frame
  nodes <- data.frame()
  for (i in 1:length(groups)){
    nodes.tmp <- data.frame(
      name = paste0(unique(cell_data[,i]),":col",i))
    nodes <- rbind(nodes, nodes.tmp)
  }
  nodes$node <- 0:(nrow(nodes)-1)
  nodes$label <- gsub(":col[0-9]+$", "", nodes$name)
  
  # Create links data frame
  links <- data.frame()
  for (i in 1:(length(groups)-1)){
    links.tmp <- cell_data %>%
      group_by(.data[[groups[i]]], .data[[groups[i+1]]]) %>%
      summarise(value = n()) %>%
      ungroup() %>%
      as.data.frame()
    links.tmp$source <- match(paste0(links.tmp[,1],":col",i), nodes$name) - 1
    links.tmp$target <- match(paste0(links.tmp[,2],":col",i+1), nodes$name) - 1
    links <- rbind(links, links.tmp[,c("source","target","value")])
  }
  
  links$left <- nodes[links$source + 1, 'name']
  
  if (is.null(height)){
    height <- max(table(gsub(".*:col", "", links$left))) * 35
  }
  
  if (is.null(width)){
    width <- (length(groups) -1) * 250
  }
  
  # Create Sankey diagram
  sankey_diagram <- sankeyNetwork(
    Links = as.data.frame(links),
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "label",
    sinksRight = FALSE,
    nodeWidth = width.node,
    fontSize = size.text,
    LinkGroup = 'left',
    NodeGroup = NULL,
    height = height,
    width = width
  )
  
  # Save the diagram as an HTML file
  if (save){
    if (is.null(filename)){
      filename <- paste0("sankey_",paste0(groups, collapse = "_"),".html")
    }
    saveNetwork(sankey_diagram, filename, selfcontained = TRUE)
    webshot(filename, sub(".html",".pdf",filename), vwidth = width + 50)
  }
  
  return(sankey_diagram)
}


#' @title Scatter plot with margin trend panels
#'
#' @description
#' Create a scatter plot of two variables colored by a third, with loess
#' trend panels on the left and bottom margins showing how the color variable
#' changes along each axis. Useful for visualizing how cell states or
#' transitions relate to spatial positions defined by distance metrics.
#'
#' @param df A data.frame containing the columns to plot.
#' @param x.by Column name for the x-axis.
#' @param y.by Column name for the y-axis.
#' @param color.by Column name for coloring points. Must be numeric
#'   (continuous).
#' @param cols Color palette for the gradient. Defaults to
#'   \code{c("#00E0FF", "yellow", "orange", "red")}.
#' @param point.size Point size for the main scatter plot. Default 0.5.
#' @param point.alpha Point transparency. Default 0.65.
#'
#' @return A \code{patchwork} object combining three panels.
#'
#' @examples
#' \dontrun{
#' plotWithMargin(df, x.by = "minDist_KC_basal",
#'               y.by = "minDist_CD4",
#'               color.by = "transition_rank")
#' }
#' @import ggplot2
#' @importFrom scales squish
#' @import patchwork
#' @export
plotWithMargin <- function(df, x.by, y.by, color.by,
                           cols = c("#00E0FF", "yellow", "orange", "red"),
                           point.size = 0.5, point.alpha = 0.65) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required. Install with:\n",
         "  install.packages('patchwork')")
  }

  global_min <- min(df[[color.by]], na.rm = TRUE)
  global_max <- max(df[[color.by]], na.rm = TRUE)

  x_range <- range(df[[x.by]], na.rm = TRUE)
  x_range <- x_range + c(-0.02, 0.02) * diff(x_range)
  x_breaks <- pretty(x_range, n = 5)
  y_range <- range(df[[y.by]], na.rm = TRUE)
  y_range <- y_range + c(-0.02, 0.02) * diff(y_range)
  y_breaks <- pretty(y_range, n = 5)

  my_color_scale <- scale_color_gradientn(
    colours = cols,
    limits  = c(global_min, global_max),
    oob     = scales::squish,
    name    = color.by
  )

  # Main scatter
  p_main <- ggplot(df, aes(
    x = .data[[x.by]],
    y = .data[[y.by]],
    color = .data[[color.by]]
  )) +
    geom_point(size = point.size, alpha = point.alpha, shape = 16) +
    my_color_scale +
    scale_x_continuous(limits = x_range, breaks = x_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = y_range, breaks = y_breaks, expand = c(0, 0)) +
    theme_bw() +
    theme(
      plot.margin  = margin(0, 0, 0, 0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right"
    ) +
    labs(color = color.by)

  # Bottom margin: loess of color.by ~ x.by
  p_bot <- ggplot(df, aes(x = .data[[x.by]], y = .data[[color.by]])) +
    geom_smooth(
      aes(color = after_stat(y)),
      method  = "loess",
      formula = y ~ x,
      se      = FALSE,
      linewidth = 0.8
    ) +
    scale_x_continuous(limits = x_range, breaks = x_breaks, expand = c(0, 0)) +
    my_color_scale +
    theme_bw() +
    theme(
      plot.margin      = margin(0, 0, 0, 0),
      axis.title.y     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      legend.position  = "none"
    ) +
    xlab(x.by) +
    scale_y_reverse()

  # Left margin: loess of color.by ~ y.by (flipped)
  p_left <- ggplot(df, aes(x = .data[[y.by]], y = .data[[color.by]])) +
    geom_smooth(
      aes(color = after_stat(y)),
      method  = "loess",
      formula = y ~ x,
      se      = FALSE,
      linewidth = 0.8
    ) +
    scale_x_continuous(limits = y_range, breaks = y_breaks, expand = c(0, 0)) +
    coord_flip() +
    my_color_scale +
    theme_bw() +
    theme(
      plot.margin      = margin(0, 0, 0, 0),
      axis.text.y      = element_blank(),
      axis.text.x      = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      axis.ticks.y     = element_blank(),
      axis.title.x     = element_blank(),
      legend.position  = "none"
    ) +
    xlab(y.by) +
    scale_y_reverse()

  # Assemble with patchwork
  layout_design <- "
AB
.C
"
  final <- patchwork::wrap_plots(
    A = p_left,
    B = p_main,
    C = p_bot,
    design = layout_design
  ) +
    patchwork::plot_layout(
      widths  = c(1, 4),
      heights = c(4, 1),
      guides  = "collect"
    )

  return(final)
}
