#' @importFrom utils globalVariables
if(getRversion() >= "2.15.1")  utils::globalVariables(c("A", "Left", "Right"))

#' @aliases confints plot.confints
#' @name confints
#' @title Confidence Intervals of Effect Differences
#'
#' @param X1 \code{data.frame} containing first effect or \code{GEM} object for ER matrix based intervals.
#' @param X2 \code{data.frame} containing second effect.
#' @param confidence Level of confidence, default = 0.95.
#' @param factor (\code{character} or \code{numeric}) indicating which factor to use in ER based intervals (defaul = 1).
#' @param levels \code{vector} (\code{character} or \code{numeric}) indicating which factor levels to. use in ER based intervals (default = c(1,2)).
#' @param x Object of class \code{confint}.
#' @param y Not used.
#' @param xlab X label (\code{character})
#' @param ylab Y label (\code{character})
#' @param sorted Logical indicating if intervals should be sorted according to their mean values, or a vector of indices/labels to sort by.
#' @param labels Logical indicating if sample labels should be used on x axis.
#' @param nonZero Logical indicating if intervals are required not to include zero.
#' @param xlim Limits of the horizontal scale.
#' @param ylim Limits of the vertical scale.
#' @param text.pt Size scaling of text in the plot (default = 16).
#' @param df.used Optional argument indicating how many degrees of freedom have been consumed during deflation. Default = 0.
#' @param ... Further arguments to \code{qplot}.
#'
#' @seealso Analyses using \code{GEM}: \code{\link{elastic}}, \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}, \code{\link{pls}}.
#'
#' @return An object of class \code{confints}, which holds
#' the information needed to perform statistics or plot the
#' confidence intervals is returned from \code{confints}.
#' The plotting routine returns a ggplot structure for plotting.
#'
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom gridExtra grid.arrange
#' @examples
#' data(MS)
#' # Subset to reduce runtime in example
#' MS$proteins <- MS$proteins[,20:70]
#'
#' # Compare MS and non-MS patients within group 1
#' conf <- with(MS, confints(proteins[MS == "yes" & group == 1,],
#'                           proteins[MS == "no"  & group == 1,]))
#' p1 <- plot(conf)
#' p2 <- plot(conf, nonZero = TRUE) # Only intervals without 0.
#' grid.arrange(p1,p2)
#'
#' # Comparison repeated but based on ER matrices
#' gem <- GEM(proteins ~ MS * group, data = MS)
#' print(effs <- colnames(gem$symbolicDesign)) # Inspect factor names
#' confGEM <- confints(gem, factor=effs[3], levels=c("yes.1","no.1"))
#' p1g <- plot(confGEM)
#' p2g <- plot(confGEM, nonZero = TRUE) # Only intervals without 0.
#' grid.arrange(p1g,p2g)
#'
#' # Shorter plot with labels
#' confShort <- conf[1:10,]
#' p1 <- plot(confShort, labels = TRUE)
#' p2 <- plot(confShort, labels = TRUE, nonZero = TRUE)
#' grid.arrange(p1,p2)
#' @export
confints <- function(X1, ...) UseMethod("confints")

#' @rdname confints
#' @export
confints.default <-function(X1, X2, confidence = 0.95, df.used = 0, ...){
  k <- dim(X1)[2]
  A1 <- apply(X1, 2, mean, na.rm=TRUE)
  A2 <- apply(X2, 2, mean, na.rm=TRUE)
  S1 <- apply(X1, 2, sd, na.rm=TRUE)
  S2 <- apply(X2, 2, sd, na.rm=TRUE)
  N1 <- apply(X1, 2, function(i)length(na.omit(i)))
  N2 <- apply(X2, 2, function(i)length(na.omit(i)))
  A  <- A2-A1
  Error <- qt(1-(1-confidence)/2, df = N1+N2-2-df.used) * sqrt(S1^2/N1+S2^2/N2)
  Left  <- A - Error
  Right <- A + Error
  names(A)   <- colnames(X1)
  Out        <- data.frame(A=A, A1=A1, A2=A2, Left=Left, Right=Right)
  class(Out) <- c("confints", "data.frame")
  return(Out)
}

#' @rdname confints
#' @export
confints.GEM <- function(X1, factor = 1, levels = c(1,2), confidence = 0.95, df.used = X1$df.used, ...){
#  dat <- X1$LS[[factor]] + X1$residuals
  dat <- X1$ER.values[[factor]]
  design <- X1$model.frame[-1][[factor]]
  df.used <- df.used + 0
  if(is.numeric(levels)){
    X1 <- dat[design==levels(design)[levels[1]],, drop=FALSE]
    X2 <- dat[design==levels(design)[levels[2]],, drop=FALSE]
  } else {
    X1 <- dat[design==levels[1],, drop=FALSE]
    X2 <- dat[design==levels[2],, drop=FALSE]
  }
  confints(X1, X2, confidence=confidence, df.used=df.used, ...)
}

#' @export
#' @rdname confints
plot.confints <- function(x, y, xlab = '', ylab = 'values',
                          sorted = TRUE, labels = FALSE, nonZero = FALSE,
                          xlim = NULL, ylim = NULL, text.pt = 12, ...){
  if(is.logical(sorted) && sorted){
    x <- x[order(x$A),,drop=FALSE]
  } else {
    x <- x[sorted,,drop=FALSE]
  }
  if(nonZero){
    x <- x[x$Left>=0 | x$Right <= 0,,drop=FALSE]
    if(dim(x)[1] == 0){
      stop('No intervals without zero overlap.')
    }
  }
  x$x <- 1:length(x$A)
  h <- qplot(x, A, data = x, geom = "path", ylab = ylab, xlab = xlab, ...)
  h <- h + geom_ribbon(aes(ymin = Left, ymax = Right), fill = "grey70", color="gray40") + geom_line(aes(y = A))
  if(!is.null(xlim)){
    h <- h + scale_x_continuous(limits = xlim, oob=squish, expand=c(0,0))
  }
  if(!is.null(ylim)){
    h <- h + scale_y_continuous(limits = ylim, oob=squish, expand=c(0,0))
  }
  if(labels){
    h <- h + scale_x_continuous(labels = rownames(x), breaks = 1:nrow(x)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=text.pt))
  }
  h + theme(text = element_text(size=text.pt))
}

#' @export grid.arrange
gridExtra::grid.arrange
