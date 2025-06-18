#' Neural Network by Multilayer Perceptron
#'
#' @param object Object of class \code{GEM}.
#' @param factor The factor to be used as response.
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
neuralnet <- function(object, factor = 1, hidden = c(2), linear.output = FALSE, ...){
  if(is.numeric(factor))
    factor <- names(object$LS)[factor]
  dat <- data.frame(object$ER.values[[factor]])
  dat[[factor]] <- object$symbolicDesign[[factor]]
  nnMod <- neuralnet::neuralnet(formula(paste0(factor, "~.")), data = dat,
                          hidden = c(2), linear.output = FALSE)
  nnMod
}

