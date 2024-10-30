#' @name bin2d
#' @title Binning 2d points
#'
#' @description Binning 2d points for count in each bin or calculation with weight assigned to each point.
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




#' @title Calculate Wasserstein distance
#'
#' @description
#' Calculate Wasserstein distance between two wide matrices stored in bins1 and bins2 objects.
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




#' @title Constructor for the wpp Class
#'
#' @description
#' Construct an object of class "wpp" from a wide matrix stored in bins object.
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






#' @title Calculate Euclidean distance
#'
#' @description
#' Extract cell coordinates by name of reduction, or fov name of the slide image.
#'
#' @param coord1 A numeric vector of coordinates for the first point.
#' @param coord2 A numeric vector of coordinates for the other point.
#' @examples
#' calcDistance.euclidean(coord1 = c(1,2,3), coord2 = c(4,5,6))
#' @export

calcDistance.euclidean <- function(coord1=NULL, coord2=NULL){
  dif <- coord1 - coord2
  d <- (sum(dif^2))^(1/length(coord1))
  return(d)
}



#### Gaussian distribution ####


#' @title Gives the density value of a Gaussian distribution
#'
#' @description
#' Calculate the density value of point x at distance d to the center of the Gaussian distribution with parameters sd and h
#' Multiple the y outside of this function if weighing the value by gene expression value of the cell.
#'
#' @param d Distance to the center.
#' @param sd Standard deviation of the Gaussian distribution.
#' @param h Height of the center. If NULL, the AUC is 1.
#' @param scale.factor Scale up the y value to avoid long decimals.
#' @examples
#' x <- seq(-5,5,0.05)
#' plot(x, calcValue.gaussian(d=x, sd=1))
#' @export
#'


calcValue.gaussian <- function(d, sd, h=NULL, scale.factor=1000){
  y <- 1/sd/sqrt(2*pi)*exp((-1)/2*(d/sd)^2)*scale.factor
  return(y)
}



#' @title Find the distance d that the Gaussian distribution has a certain value.
#'
#' @description
#' Find the distance d that the Gaussian distribution has a certain value.
#' For any points further than d, use 0 as approximate to simplify the calculation of value in field.
#'
#' @param sd Standard deviation of the Gaussian distribution.
#' @param value The value of the Gaussian distribution at distance d.
#' @examples
#' findDistance.euclidean(value=0.01, sd=1)
#' @export
#'

findDistance.euclidean <- function(value, sd){
  x <- sqrt((-1)*2*log(value*(sd*sqrt(2*pi))))*sd
  return(x)
}








#### create a field with coordinates and values and parameters ####


#' @title Create a field object
#'
#' @description
#' Simply store the info, for getValueInField function to read.
#'
#' @param coords A data.frame with each column represent one dimension. If there's a "value" column, move that into the weight section.
#' @param method.distr "gaussian" by default.
#' @param method.d "euclidean" by default.
#' @param weight A vector of the same length as coords, scaling the distribution.
#' @param d.cutoff Value for point further than this distance will have value 0.
#' @param ... Other parameters. sd is the standard deviation of the normal distribtion. h is the height of the center point.
#' @examples
#' coords.ifnb1 <- getCoords.transcript(srt, transcript="IFNB1", fov="UV109fov1")
#' findDistance.euclidean(0.00001,sd=1000)
#' field <- createField(coords.ifnb1, sd=1000, h=NULL, d.cutoff=2715, scale.factor=10000)
#' @export


createField <- function(coords, method.distr="gaussian", method.d="euclidean", weight=NULL, d.cutoff, ...){
  dots <- match.call(expand.dots = FALSE)[["..."]]
  if ("value" %in% colnames(coords)){
    ind <- grep("value", colnames(coords))
    weight <- coords[,ind]
    coords <- coords[,-ind]
  }
  field <- list(coords = coords,
                method.distr = method.distr,
                method.d = method.d,
                parameters = c(dots,d.cutoff=d.cutoff),
                weight = weight)
  return(field)
}




#' @title Calculate the value of the field at each query point
#'
#' @description
#' Calculate the value of the field at each query point.
#'
#' @param field A field object created by createField.
#' @param coords A data.frame with each column as a dimension.
#' @examples
#' coords.moDC <- getCoords.cell(srt, ind = srt$subcelltype == "monocyte", fov = "roi1")
#' values <- getValueInField(field, coords.moDC)
#' plot(sort(values))
#' plot(density(values))
#' @export
#'

getValueInField <- function(field, coords){
  if (ncol(coords) != ncol(field$coords)){
    stop("field and coords have different dimentions.")
  }

  if (field$method.distr == "gaussian"){
    f.distr <- calcValue.gaussian
    att <- attributes(f.distr)
    formals(f.distr)[names(field$parameters)] <- field$parameters
    attributes(f.distr) <- att[names(att) != "srcref"]
  }else{
    stop("distribution method not supported")
  }

  if (field$method.d == "euclidean"){
    f.d <- calcDistance.euclidean
  }

  # loop through all query dots
  values.query <- c()
  # leave the is.null(weight) outside of the loop to avoid repetition
  if (!is.null(field$weight)){
    for (i in 1:nrow(coords)){
      coords.p <- as.numeric(coords[i,])
      d.all <- apply(field$coords, 1, function(x){
        d <- f.d(coord1 = as.numeric(x), coord2 = coords.p)
        return(d)
      })
      # set value to 0 if too far
      d.cutoff <- field$parameters$d.cutoff
      ind <- d.all < d.cutoff
      value <- sum(f.distr(d.all[ind]) * field$weight[ind]) # the only difference
      values.query <- c(values.query, value)
    }
  }else{
    for (i in 1:nrow(coords)){
      coords.p <- as.numeric(coords[i,])
      d.all <- apply(field$coords, 1, function(x){
        d <- f.d(coord1 = as.numeric(x), coord2 = coords.p)
        return(d)
      })
      # set value to 0 if too far
      d.cutoff <- field$parameters$d.cutoff
      ind <- d.all < d.cutoff
      value <- sum(f.distr(d.all[ind]))  # the only difference
      values.query <- c(values.query, value)
    }
  }

  return(values.query)
}






#### visualization of the field ####


#' @title Report the size of the image
#'
#' @description
#' Output in the order of: x.min, x.max, y.min, y.max.
#'
#' @param object Seurat object.
#' @examples
#' binlim <- getSize(srt, fov = "roi1")
#' @import Seurat
#' @export
#'

getSize <- function(object, fov){
  tmp.m <- sapply(object@images[[fov]]@molecules$molecules, function(m){
    return(m@bbox)
  })
  tmp.c <- object@images[[fov]]@boundaries$centroids@coords
  x.min <- min(min(tmp.m[1,]),min(tmp.c[,1]))
  x.max <- max(max(tmp.m[3,]),max(tmp.c[,1]))
  y.min <- min(min(tmp.m[2,]),min(tmp.c[,2]))
  y.max <- max(max(tmp.m[4,]),max(tmp.c[,2]))
  res <- c(x.min-1,x.max+1,y.min-1,y.max+1)
  return(res)
}



