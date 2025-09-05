#' Neural Network by Multilayer Perceptron
#'
#' @param object Object of class \code{GEM}.
#' @param formula A formula specifying the model to be fitted. If not provided, the response variable is taken from the \code{GEM} object.
#' @param factor The factor to be used as response. If \code{formula} is provided, this is ignored.
#' @param hidden Vector with numbers of neurons in the hidden layers.
#' @param linear.output Logical. If \code{TRUE} the output layer is linear, otherwise it is logistic.
#' @param ... Additional arguments passed to \code{neuralnet}.
#'
#' @returns A \code{neuralnet} object that can be inspected and plotted.
#' @importFrom stats formula
#' @seealso Analyses using \code{GEM}: \code{\link{elastic}}, \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}, \code{\link{pls}}.
#' Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
#' @export
#'
#' @examples
#' data(candies, package = "HDANOVA")
#' gemC <- GEM(assessment ~ assessor*candy, data=candies)
#'
#' # Neural network model
#' nn <- neuralnet(gemC, factor = "candy", hidden = c(2))
#' plot(nn, rep="best")
#'
#' # Network weights (input and hidden layers)
#' nn$weights
neuralnet <- function(object, formula, factor = 1, hidden = c(2), linear.output = FALSE, ...){
  if(missing(formula)){
    nam <- names(object$symbolicDesign)
    names(nam) <- nam
    formula <- stats::formula(paste0(nam[factor], "~."))
  } else {
    factor <- as.character(formula[[2]])
  }

  if(is.numeric(factor))
    factor <- names(object$LS)[factor]
  dat <- data.frame(object$ER.values[[factor]])
  dat[[factor]] <- object$symbolicDesign[[factor]]
  if(!requireNamespace("neuralnet", quietly = TRUE))
    stop("Package 'neuralnet' needed for this function to work. Please install it.",
         call. = FALSE)
  nnMod <- neuralnet::neuralnet(formula, data = dat,
                          hidden = hidden, linear.output = linear.output, ...)
  nnMod
}

