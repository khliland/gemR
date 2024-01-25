#' Plot function for GEM-based PLS
#' @param x \code{GEMpls} object
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
plot.GEMpls <- function(x,y, ylab="% correct", xlab="# comp.",
                        main=paste0("CV accuracy (",x$effect,")"),
                        panel.first=grid(), type="o", ...){
  corr <- colMeans(x$classes == as.numeric(x$data$y))*100
  plot(corr, ylab=ylab, xlab=xlab, main=main, panel.first=panel.first, ...)
}
