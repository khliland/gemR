#' Simultaneous Component Analysis
#'
#' @description This function performs Simultaneous Component Analysis (SCA) on a \code{hdanova} object.
#'
#' @param object A \code{hdanova} object.
#'
#' @returns An updated \code{hdanova} object with SCA results.
#' @seealso Analyses using \code{GEM}: \code{\link{elastic}}, \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}, \code{\link{pls}}.
#' Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
#' @examples
#' # Load candies data
#' data(candies, package="HDANOVA")
#'
#' # Basic HDANOVA model with two factors
#' mod <- GEM(assessment ~ candy + assessor, data=candies)
#' mod <- sca(mod)
#' scoreplot(mod)
#' @importFrom HDANOVA sca
#' @export
sca <- function(object){
  obj <- HDANOVA::sca(object)
  obj$effects <- object$effectsOrig
  obj
}
