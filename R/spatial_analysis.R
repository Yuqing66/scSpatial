

#' @title Get coordinates of clicks.
#' @description
#' Display the plot and get the coordinates of clicks.
#' @param plot.object An object storing the plot.
#' @return A data frame with x,y coordinates.
#' @examples
#' g <- ImageDimPlot.ssc(srt, fov = fov, group.by = NULL)
#' coords <- getClickCoordinates(g)
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
#' calcDistance.toLine(x, y, p1, p2)
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
#' g <- ImageDimPlot.ssc(srt, fov = fov, group.by = NULL)
#' coords <- getClickCoordinates(g)
#' d <- image.calcDistanceTo(srt, fov = fov, coords.df = coords)
#' @import Seurat
#' @export


image.calcDistanceTo <- function(object, fov, shape = "line", ...){
  dots <- list(...)
  mc <- match.call(expand.dots = FALSE)[["..."]]

  # get cell coordinates on this image
  cell.coords <- getCoords.cell(object, fov = fov)

  if (shape == "line"){
    d <- calcDistance.toLine(x=cell.coords$x, y=cell.coords$y, ...)
  }else{
    stop("Only line is supported for now.")
  }
  return(d)
}


