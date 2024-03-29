#' Plot function and extraction method for GEM-based PLS
#' @name plot.GEMpls
#' @aliases plot.GEM.pls explvar.GEMpls print.GEMpls summary.GEMpls
#' @param x,object \code{GEMpls} object
#' @param y Not used
#' @param ylab \code{character} label for Y axis.
#' @param xlab \code{character} label for X axis.
#' @param main  \code{character} main header.
#' @param panel.first \code{expression} to be evaluated after the plot axes are set up but before any plotting takes place. This can be useful for drawing background grids or scatterplot smooths. Note that this works by lazy evaluation: passing this argument from other plot methods may well not work since it may be evaluated too early.
#' @param type \code{character} plot type (see ?plot).
#' @param ... additional arguments for \code{plot()}.
#'
#' @importFrom graphics grid
#' @export
plot.GEMpls <- function(x,y, ...){
  if(!is.null(x$shave)){
    plot(x$shave, ...)
  } else {
    .plot.not.shave(x, ...)
  }
}

.plot.not.shave <- function(x, y, ylab="% correct", xlab="# comp.",
                main=paste0("CV accuracy (",x$effect,")"),
                panel.first=grid(), type="o", ...){
  corr <- colMeans(x$classes == as.numeric(x$data$y))*100
  plot(corr, ylab=ylab, xlab=xlab, main=main, panel.first=panel.first, ...)
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
#' @export
scores <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::scores(object, ...)
}

#' @rdname plot.GEMpls
#' @export
scoreplot <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::scoreplot(object, ...)
}

#' @rdname plot.GEMpls
#' @export
loadings <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::loadings(object, ...)
}

#' @rdname plot.GEMpls
#' @export
loadingplot <- function(object, ...){
  class(object) <- c("mvr", "list")
  pls::loadingplot(object, ...)
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
