#' @importFrom utils globalVariables
if(getRversion() >= "2.15.1")  utils::globalVariables(c("A", "Left", "Right"))

#' @aliases confints plot.confints
#' @name confints
#' @title Confidence Intervals of Effect Differences
#'
#' @param X1 \code{data.frame} containing first effect.
#' @param X2 \code{data.frame} containing second effect.
#' @param confidence Level of confidence, default = 0.95.
#' @param x Object of class \code{confint}.
#' @param y Not used.
#' @param xlab X label (\code{character})
#' @param ylab Y label (\code{character})
#' @param sorted Logical indicating if intervals should be sorted according to their mean values.
#' @param nonZero Logical indicating if intervals are required not to include zero.
#' @param xlim Limits of the horizontal scale.
#' @param ylim Limits of the vertical scale.
#' @param text.pt Size scaling of text in the plot (default = 16).
#' @param df.used Optional argument indicating how many degrees of freedom have been consumed during deflation. Default = 0.
#' @param ... Further arguments to \code{qplot}.
#'
#' @return An object of class \code{confints}, which holds
#' the information needed to perform statistics or plot the
#' confidence intervals is retunred from \code{confints}.
#' The plotting routine returns a ggplot structure for plotting.
#'
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' data(MS)
#' # Compare MS and non-MS patients within cluster 1
#' conf <- with(MS, confints(proteins[MS == "yes" & cluster == 1,],
#'                           proteins[MS == "no"  & cluster == 1,]))
#' p1 <- plot(conf)
#' p2 <- plot(conf, nonZero = TRUE) # Only intervals without 0.
#' grid.arrange(p1,p2)
confints <-function(X1, X2, confidence = 0.95, df.used = 0){
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

#' @export
#' @rdname confints
plot.confints <- function(x, y, xlab = '', ylab = 'normalised log2', sorted = TRUE, nonZero = FALSE,
                          xlim = NULL, ylim = NULL, text.pt = 16, ...){
  if(sorted)
    x <- x[order(x$A),,drop=FALSE]
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
  h + theme(text = element_text(size=text.pt))
}

#' @export grid.arrange
gridExtra::grid.arrange
