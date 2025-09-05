#' @aliases pls.GEM
#' @title Partial Least Squares modelling of GEM objects.
#' @param object Object of class \code{GEM}.
#' @param effect The effect to be used as response.
#' @param ncomp Number of PLS components.
#' @param newdata Optional new data matrix for prediction.
#' @param gem2 Second object of class \code{GEM} for comparison.
#' @param validation Optional validation parameters for \code{plsr}.
#' @param jackknife Optional argument specifying if jackknifing should be applied.
#' @param shave Optional argument indicating if variable shaving should be used. \code{shave} should be a list with two elements: the PLS filter method and the proportion to remove. \code{shave = TRUE} uses defaults: \code{list("sMC", 0.2)}.
#' @param df.used Optional argument indicating how many degrees of freedom have been consumed during deflation. Default value from input object.
#' @param ... Additional arguments for \code{plsr}.
#'
#' @description The output of \code{GEM} is used as input to a PLS classification with the selected
#' effect as response. It is possible to compare two models using the \code{gem2} argument. Variable
#' selection is available through Jackknifing (from package \code{pls}) and Shaving (from package \code{plsVarSel}).
#'
#' @details If using the \code{shave} options, the segment type is given as \code{type} instead of \code{segment.type} (see examples).
#'
#' @seealso Analyses using \code{GEM}: \code{\link{elastic}}, \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}.
#' Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
#'
#' @importFrom plsVarSel shaving lda_from_pls lda_from_pls_cv
#' @importFrom pls plsr cvsegments var.jack scores scoreplot loadings loadingplot corrplot R2 mvrValstats explvar
#' @examples
#' data(MS, package = "gemR")
#' # Subset to reduce runtime in example
#' MS$proteins <- MS$proteins[,20:70]
#'
#' gem <- GEM(proteins ~ MS * group, data = MS[-1,])
#'
#' # Simple PLS using interleaved cross-validation
#' plsMod <- pls(gem, 'MS', 6, validation = "CV",
#'               segment.type = "interleaved", length.seg = 25)
#' plot(plsMod)
#' scoreplot(plsMod, labels = "names")
#'
#' # PLS with shaving of variables (mind different variable for cross-validation type)
#' plsModS <- pls(gem, 'MS', 6, validation = "CV",
#'               type = "interleaved", length.seg=25, shave = TRUE)
#' # Error as a function of remaining variables
#' plot(plsModS)
#' # Selected variables for minimum error
#' with(plsModS$shave, colnames(X)[variables[[min.red+1]]])
#'
#' \donttest{ # Time consuming due to leave-one-out cross-validation
#'   plsModJ <- pls(gem, 'MS', 5, validation = "LOO",
#'               jackknife = TRUE)
#'   colSums(plsModJ$classes == as.numeric(MS$MS[-1]))
#'   # Jackknifed coefficient P-values (sorted)
#'   plot(sort(plsModJ$jack[,1,1]), pch = '.', ylab = 'P-value')
#'   abline(h=c(0.01,0.05),col=2:3)
#'
#'   scoreplot(plsModJ)
#'   scoreplot(plsModJ, comps=c(1,3))   # Selected components
#'   # Use MS categories for colouring and clusters for plot characters.
#'   scoreplot(plsModJ, col = gem$symbolicDesign[['MS']],
#'                   pch = 20+as.numeric(gem$symbolicDesign[['group']]))
#'   loadingplot(plsModJ, scatter=TRUE) # scatter=TRUE for scatter plot
#' }
#' @return An object of class \code{GEMpls, mvr, list} containing the fitted PLS model, classifications/predictions, data and optionally Jackknife or Shaving results.
#' @export
pls.GEM <- function(object, effect, ncomp, newdata = NULL, gem2, validation, jackknife = NULL, shave = NULL, df.used = object$df.used, ...){
  classification <- FALSE
  if(!missing(gem2)){
    data <- data.frame(X = I(object$ER.values[[effect]]),
                       y = I(gem2$ER.values[[effect]]))
    lda    <- NULL
  } else {
    if(length(effect) == 1){
      if(is.factor(object$symbolicDesign[[effect]])){
        classification <- TRUE
        data <- data.frame(X = I(object$ER.values[[effect]]),
                           y = I(object$symbolicDesign[[effect]]),
                           Yd = I(model.matrix(~y-1,data.frame(y=object$symbolicDesign[[effect]]))))
      } else {
        data <- data.frame(X = I(object$ER.values[[effect]]),
                           y = I(object$symbolicDesign[[effect]]),
                           Yd = object$symbolicDesign[[effect]])
      }
    } else { # User supplied contrast
      data <- data.frame(X = I(object$ER.values[[effect[[1]]]]),
                         y = effect[[2]],
                         Yd = effect[[2]])
    }
  }
  jack   <- ifelse(is.null(jackknife), FALSE, TRUE)
  shaved <- ifelse(is.null(shave), FALSE, TRUE)
  if(!is.null(newdata)){
    if(missing(gem2)){
      plsMod <- plsr(Yd ~ X, ncomp = ncomp, data = data, ...)
      if(classification)
        lda    <- lda_from_pls(plsMod, data$y, newdata, ncomp)
    } else {
      plsMod <- plsr(y ~ X, ncomp = ncomp, data = data, ...)
    }
  } else {
    if(missing(validation))
      validation <- "LOO"
    if(shaved){
      if(is.logical(shave) && shave)
        shave <- list("sMC", 0.2)
      sh <- shaving(data$y, data$X, ncomp = ncomp, method = shave[[1]], prop = shave[[2]], validation = validation, ...)
      lda <- plsMod <- NULL
    } else {
      if(missing(gem2)){
        plsMod <- plsr(Yd ~ X, ncomp = ncomp, data = data, validation = validation, jackknife = jack, ...)
        if(classification)
          lda    <- lda_from_pls_cv(plsMod, data$X, data$y, ncomp)
      } else {
        plsMod <- plsr(y ~ X, ncomp = ncomp, data = data, validation = validation, jackknife = jack, ...)
      }
      if(jack){
        jack.test <- function(object, ncomp = object$ncomp, use.mean = TRUE, df.used = 0){
          nresp <- dim(object$coefficients)[2]
          sdjack <- sqrt(var.jack(object, ncomp = ncomp, covariance = FALSE,
                                  use.mean = use.mean))
          B <- coef(object, ncomp = ncomp)
          df <- length(object$validation$segments) - 1 - df.used
          tvals <- B/sdjack
          pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
          structure(list(coefficients = B, sd = sdjack, tvalues = tvals,
                         df = df, pvalues = pvals), class = "jacktest")
        }
        if(missing(gem2)){
          jt <- array(0, dim = c(dim(data$X)[2], ifelse(is.null(dim(data$Yd)), 1, dim(data$Yd)[2]), ncomp))
          dimnames(jt) <- list(colnames(data$X), colnames(data$Yd), paste('comp.', 1:ncomp, sep=""))
        } else {
          jt <- array(0, dim = c(dim(data$X)[2], ifelse(is.null(dim(data$y)), 1, dim(data$y)[2]), ncomp))
          dimnames(jt) <- list(colnames(data$X), colnames(data$y), paste('comp.', 1:ncomp, sep=""))
        }
        for(i in 1:ncomp){
          jt[,,i] <- jack.test(plsMod, ncomp = i, df.used = df.used)$pvalues
        }
      }
    }
  }
  object <- plsMod
  if(classification)
    object$classes <- lda
  object$data    <- data
  object$effect  <- effect
  object$call.GEMpls <- match.call()
  object$object     <- object
  # object <- list(classes = lda, data = data, pls = plsMod)
  if(jack)
    object$jack <- jt
  if(shaved){
    object$shave <- sh
    class(object) <- c('GEMpls','list')
  } else {
    class(object) <- c('GEMpls','mvr','list')
  }
  object
}

#' @importFrom HDANOVA pls
#' @export
HDANOVA::pls

