#' Plot function and extraction method for GEM-based PLS
#' @name plot.GEMpls
#' @aliases plot.GEM.pls explvar.GEMpls print.GEMpls summary.GEMpls
#' @param x,object \code{GEMpls} object
#' @param y Not used
#' @param ylab \code{character} label for Y axis.
#' @param xlab \code{character} label for X axis.
#' @param main \code{character} main header.
# #' @param panel.first \code{expression} to be evaluated after the plot axes are set up but before any plotting takes place. This can be useful for drawing background grids or scatterplot smooths. Note that this works by lazy evaluation: passing this argument from other plot methods may well not work since it may be evaluated too early.
# #' @param type \code{character} plot type (see ?plot).
#' @param ... additional arguments for \code{plot()}.
#' @return A plot of the PLS model's cross-validated accuracy or the Shaving results if available.
#' @seealso \code{\link{pls}} has examples of how to use this function.
#'
#' @importFrom graphics grid
#' @export
plot.GEMpls <- function(x,y, ylab="error", xlab="nvar", main = "Shaving", ...){
  if(!is.null(x$shave)){
    plot(x$shave, ylab=ylab, xlab=xlab, main=main,...)
  } else {
    if(missing(ylab))
      ylab <- "% correct"
    if(missing(xlab))
      xlab <- "# comp."
    if(missing(main))
      main <- paste0("CV accuracy (",x$effect,")")
    corr <- colMeans(x$classes == as.numeric(x$data$y))*100
    plot(corr, ylab=ylab, xlab=xlab, main=main, ...)
  }
}


#' @rdname plot.GEMpls
#' @export
print.GEMpls <- function(x,...){
  if(!is.null(x$shave))
    print(x$shave, ...)
  else {
    class(x) <- c("mvr", "list")
    print(x, ...)
  }
}

#' @rdname plot.GEMpls
#' @export
summary.GEMpls <- function(object,...){
  if(!is.null(object$shave))
    summary(object$shave, ...)
  else {
    class(object) <- c("mvr", "list")
    summary(object, ...)
  }
}

#' @rdname plot.GEMpls
#' @importFrom HDANOVA scores
#' @export
scores <- function(object, ...){
  if(inherits(object, "asca")){
    HDANOVA::scores(object, ...)
  } else {
    class(object) <- c("mvr", "list")
    pls::scores(object, ...)
  }
}

#' @rdname plot.GEMpls
#' @importFrom HDANOVA scoreplot
#' @export
scoreplot <- function(object, ...){
  if(inherits(object, "asca")){
    HDANOVA::scoreplot(object, ...)
  } else {
    class(object) <- c("mvr", "list")
    pls::scoreplot(object, ...)
  }
}

#' @rdname plot.GEMpls
#' @importFrom HDANOVA loadings
#' @export
loadings <- function(object, ...){
  if(inherits(object, "asca")){
    HDANOVA::loadings(object, ...)
  } else {
    class(object) <- c("mvr", "list")
    pls::loadings(object, ...)
  }
}

#' @rdname plot.GEMpls
#' @importFrom HDANOVA loadingplot
#' @export
loadingplot <- function(object, ...){
  if(inherits(object, "asca")){
    HDANOVA::loadingplot(object, ...)
  } else {
    class(object) <- c("mvr", "list")
    pls::loadingplot(object, ...)
  }
}

#' @rdname plot.GEMpls
#' @export
corrplot <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::corrplot(object, ...)
}

#' @rdname plot.GEMpls
#' @export
R2 <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::R2(object, ...)
}

#' @rdname plot.GEMpls
#' @export
mvrValstats <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::mvrValstats(object, ...)
}

#' @rdname plot.GEMpls
#' @export
explvar <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::explvar(object, ...)
}
