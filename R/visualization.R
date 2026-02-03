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



