

#' @title Get coordinates of clicks.
#' @description
#' Display the plot and get the coordinates of clicks.
#' @param plot.object An object storing the plot.
#' @return A data frame with x,y coordinates.
#' @examples
#' \dontrun{
#' g <- ImageDimPlot.ssc(srt, fov = fov, group.by = NULL)
#' coords <- getClickCoordinates(g)
#' }
#' @export


getClickCoordinates <- function(plot.object){
  click_coords <- reactiveVal(NULL)

  # Define UI for application
  ui <- fluidPage(
    titlePanel("Click on the plot to get coordinates"),
    plotOutput("plot", click = "plot_click"),
    actionButton("done", "Done"),
    verbatimTextOutput("info")
  )

  # Define server logic
  server <- function(input, output, session) {

    # create plot
    output$plot <- renderPlot({
      plot.object
    })

    # Display click coordinates
    output$info <- renderPrint({
      req(input$plot_click$x, input$plot_click$y)
      paste("x =", input$plot_click$x, "y =", input$plot_click$y)
    })

    # Store coordinates when plot is clicked
    observeEvent(input$plot_click, {
      new_click <- data.frame(
        x = input$plot_click$x,
        y = input$plot_click$y
      )
      # updata stored coordinates
      click_coords(rbind(click_coords(), new_click))
    })

    # Close app when Done button is clicked
    observeEvent(input$done, {
      stopApp(click_coords())
    })
  }

  result <- runApp(shinyApp(ui = ui, server = server))
  return(result)
}


#' Calculate distance of points to a line.
#' @description
#' Calculate the distance of points to a line defined by two points or a slope and an intercept.
#' @param x x coordinates of the points.
#' @param y y coordinates of the points.
#' @param coords.df A data frame with x,y columns for coordinates of the two points defining the line.
#' @param p1 A vector of two elements, the x and y coordinates of the first point defining the line.
#' @param p2 A vector of two elements, the x and y coordinates of the second point defining the line.
#' @param k The slope of the line.
#' @param b The intercept of the line.
#' @return A vector of distances.
#' @examples
#' x <- 1:1
#' y <- 1:5
#' xy1 <- c(0,0)
#' xy2 <- c(5,5)
#' calcDistance.toLine(x, y, p1 = xy1, p2 = xy2)
#' calcDistance.toLine(x, y, k=1, b=0)
#' @export

calcDistance.toLine <- function(x, y, coords.df=NULL, p1=NULL, p2=NULL, k=NULL, b=NULL){
  # coords.df
  if (!is.null(coords.df)){
    p1 <- as.numeric(coords.df[1,])
    p2 <- as.numeric(coords.df[2,])
  }
  # k, b
  if (is.null(p1) && is.null(p2)){
    if (is.null(k) || is.null(b)){
      stop("Please provide xy1 and xy2, or k and b.")
    }
    tmp <- sqrt(k^2 + 1)
    d <- abs(k*x - y + b) / tmp
  }else{
    # p1, p2
    x1 <- p1[1]
    y1 <- p1[2]
    x2 <- p2[1]
    y2 <- p2[2]

    dx <- x2 - x1
    dy <- y2 - y1
    tmp1 <- x2*y1 - y2*x1
    tmp2 <- sqrt(dx^2 + dy^2)
    d <- abs(dy*x - dx*y + tmp1) / tmp2
  }
  return(d)
}



#' Calculate the distance of cells to a shape on an image.
#' @description
#' Calculate the distance of cells on an image to a defined shape such as a line.
#' @param object A Seurat object.
#' @param fov The name of the image fov in the Seurat object.
#' @param shape The shape of the line. Only "line" is supported for now.
#' @param ... Other parameters for calcDistance functions.
#' @return A vector of distances.
#' @examples
#' \dontrun{
#' g <- ImageDimPlot.ssc(srt, fov = fov, group.by = NULL)
#' coords <- getClickCoordinates(g)
#' d <- image_calcDistanceTo(srt, fov = fov, coords.df = coords)
#' }
#' @import Seurat
#' @export


image_calcDistanceTo <- function(object, fov, shape = "line", ...){
  dots <- list(...)
  mc <- match.call(expand.dots = FALSE)[["..."]]

  # get cell coordinates on this image
  cell.coords <- getCoords.cell(object, fov = fov)

  if (shape == "line"){
    d <- calcDistance.toLine(x=cell.coords$x, y=cell.coords$y, ...)
  }else{
    stop("Only line is supported for now.")
  }
  names(d) <- rownames(cell.coords)
  return(d)
}












getCellsInPath <- function(cells.coords, path.coords){
  # label the cells in vs out
  inside <- point.in.polygon(point.x = cells.coords[,1],
                             point.y = cells.coords[,2],
                             pol.x = c(path.coords$x, path.coords$x[1]),
                             pol.y = c(path.coords$y, path.coords$y[1]))
  if (is.null(rownames(cells.coords))){
    res <- inside == 1
  }else{
    res <- rownames(cells.coords)[inside == 1]
  }
  return(res)
}












# need to subset the object to only cells used to calculate the expression / location
# center.cellid should be the cell names. the boolean vector for the specific fov would work but is more error prone.
# if fovs is NULL, calculate for all images
# if include.center, the center.coords should be cells and the expression of it will be included.
calcLocalExpression <- function(srt, center.coords=NULL, center.cellid=NULL, radius, genes,
                                assay=NULL, fovs = NULL, slot = "data", include.center = FALSE){
  if (is.null(assay)){
    assay <- srt@active.assay
  }
  if (is.null(fovs)){
    fovs <- names(srt@images)
  }
  radius_2 <- radius^2
  
  # calculate for each fov separately
  loc.exp <- data.frame()
  for (i in 1:length(fovs)){
    fov <- fovs[i]
    
    # get the coordinates of the cells to calculate for, if not provided
    if (is.null(center.coords)){
      center.coords.fov <- getCoords.cell(srt, fov = fov)
      if (!is.null(center.cellid)){
        center.coords.fov <- center.coords.fov[rownames(center.coords.fov) %in% center.cellid, , drop = FALSE]
      }
    }
    
    # get the expression matrix of specific genes of all cells in srt
    exprs.coords <- getCoords.cell(srt, fov = fov)
    if (!include.center){
      exprs.coords <- exprs.coords[!rownames(exprs.coords) %in% rownames(center.coords.fov), ]
    }
    exprs <- t(srt@assays[[assay]]@layers[[slot]][match(genes, rownames(srt)), match(rownames(exprs.coords), colnames(srt))])
    colnames(exprs) <- genes
    
    # calculate the distance matrix between the center cells and all other cells
    loc.exp.fov <- apply(center.coords.fov, 1, function(r){
      x <- r[1]
      y <- r[2]
      ind <- (exprs.coords[,"x"]-x)^2 + (exprs.coords[,"y"]-y)^2 < radius_2
      if (sum(ind) > 0){
        # calculate the mean expression of the genes in the surrounding cells
        loc.exp.fov <- colSums(exprs[ind, , drop = F], na.rm = TRUE)
      } else {
        loc.exp.fov <- rep(NA, length(genes))
      }
      return(loc.exp.fov)
    })
    loc.exp.fov <- as.data.frame(t(loc.exp.fov))
    rownames(loc.exp.fov) <- rownames(center.coords.fov)
    colnames(loc.exp.fov) <- genes
    loc.exp.fov$fov <- fov
    
    loc.exp <- rbind(loc.exp, loc.exp.fov)
  }
  
  colnames(loc.exp) <- c(paste0("LocalExp_",genes, "_radius", radius, "_", assay, "_", slot), "fov")
  
  return(loc.exp)
}




