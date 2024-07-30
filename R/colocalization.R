#' Binning 2d points for count in each bin or calculation with weight assigned.
#'
#' @description
#'
#' @param x,y Coordinates.
#' @param weight Assign weight to each point. In the same order as coordinates.
#' @param nbins Number of bins on each axis in the order of c(nbin.x, nbin.y).
#' @param binlim x, y range for bin calculation, e.g. c(min(x),max(x),min(y),max(y)).
#' @param FUN Function for weight.
#' @param num.cutoff Remove bins with less than this number of points.
#' @examples
#' x <- cellCoords[,1]
#' y <- cellCoords[,2]
#' binlim <- c(min(x),max(x),min(y),max(y))
#' bin2d(x=x, y=y, nbins=c(100,100), binlim=binlim)
#' @import
#' @export

bin2d <- function(x, y, weight=NULL, nbins=c(100,100), binlim, FUN=sum, num.cutoff=NULL){
  x.cuts <- seq(from = binlim[1], to = binlim[2], length = nbins[1] + 1)
  y.cuts <- seq(from = binlim[3], to = binlim[4], length = nbins[2] + 1)
  index.x <- cut(x, x.cuts, include.lowest = TRUE)
  index.y <- cut(y, y.cuts, include.lowest = TRUE)
  if (!is.null(num.cutoff)){
    rm.ind <- tapply(1:length(x), list(index.x, index.y), function(ind){
      if (length(ind) < num.cutoff) ind
    })
    rm.ind <- unlist(rm.ind)
    if (length(rm.ind) > 0){
      x <- x[-rm.ind]
      index.x <- index.x[-rm.ind]
      index.y <- index.y[-rm.ind]
      if (!is.null(weight)) weight <- weight[-rm.ind]
    }
  }

  if (is.null(weight)){
    m <- tapply(x, list(index.x, index.y), length) #index.x as rows, index.y as columns
  }else{
    m <- tapply(weight, list(index.x, index.y), FUN = FUN)
  }
  m[is.na(m)] <- 0
  midpoints <- function(x) (x[-1] + x[-length(x)])/2
  bins <- list(x=midpoints(x.cuts), y=midpoints(y.cuts), counts=m)
  return(bins)
}







#' Plot 2d raster density plot
#'
#' @description
#' macOS Preview will do some interpolation of the raster to make the plot in pdf smoother/blurry. Adobe is fine.
#' Raster might have some white line/gap in between.
#'
#' @param x,y x, y are vectors for coordinates of z (optional).
#' @param z A count matrix.
#' @param method "tile" by default, but can change to "raster".
#' @examples
#' plot.mt(wideMatrix)
#' @import ggplot2
#' @export

plot.mt <- function(x=NULL,y=NULL,z,method="tile"){
  wide <- data.frame(z)
  if (is.null(x)) x <- 1:ncol(z)
  if (is.null(y)) y <- 1:nrow(z)
  colnames(wide) <- x
  wide$y <- y
  long <- melt(wide, id.vars = c("y"),
               variable.name = "x",
               value.name = "z")
  long$x <- as.numeric(as.vector(long$x)) # factor to numeric
  if (method == "tile"){
    g <- ggplot(long,aes(x=x,y=y))+
      geom_tile(aes(fill=z,color=z),size=0)
    # color and size are used to remove the white lines between tiles
  }else if (method == "raster"){
    g <- ggplot(long,aes(x=x,y=y))+
      geom_raster(aes(fill=z,color=z))
  }
  g <- g + theme(axis.title = element_blank(),axis.text = element_blank(),
                 axis.ticks = element_blank(),legend.title = element_blank())
  return(g)
}





#' Calculate Wasserstein distance
#'
#' @description
#' Calculate Wasserstein distance between two wide matrices stored in bins1 and bins2 objects.
#'
#' @param bins1,bins2 An object storing multiple tabs of matrices. Raw counts or after smoothing.
#' @param tab The name of a list in bins object.
#' @examples
#' wdist()
#' @import transport
#' @export


#
wdist <- function(bins1, bins2, tab="counts.smooth"){
  value <- wasserstein(getwpp(bins1, tab),getwpp(bins2, tab))
  return(value)
}




#' Constructor for the wpp Class
#'
#' @description
#' Construct an object of class "wpp" from a wide matrix stored in bins object.
#'
#' @param bins An object storing multiple tabs of matrices. Raw counts or after smoothing.
#' @param tab The name of a list in bins object.
#' @examples
#' getwpp(bins, tab="counts.smooth")
#' @import transport
#' @export



# transform the matrix to wpp object.
getwpp <- function(bins, tab){
  mt <- bins[[tab]]
  colnames(mt) <- bins$y
  rownames(mt) <- bins$x
  long <- melt(mt)
  wppobj <- wpp(coordinates = long[,c("Var1","Var2")], mass = long$value)
  return(wppobj)
}













