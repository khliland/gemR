#' @aliases pls pls.ER
#' @name pls
#' @title Partial Least Squares modelling of ER objects.
#' @param er Object of class \code{ER}.
#' @param effect The effect to be used as response.
#' @param ncomp Number of PLS components.
#' @param newdata Optional new data matrix for prediction.
#' @param er2 Second object of class \code{ER} for comparison.
#' @param validation Optional validation parameters for \code{plsr}.
#' @param jackknife Optional argument specifying if jackknifing should be applied.
#' @param shave Optional argument indicating if variable shaving should be used. \code{shave} should be a list with two elements: the PLS filter method and the proportion to remove. \code{shave = TRUE} uses defaults: \code{list("sMC", 0.2)}.
#' @param df.used Optional argument indicating how many degrees of freedom have been consumed during deflation. Default value from input object.
#' @param ... Additional arguments for \code{plsr}.
#'
#' @importFrom plsVarSel shaving lda_from_pls lda_from_pls_cv
#' @importFrom pls plsr cvsegments var.jack scores scoreplot loadings loadingplot R2 mvrValstats explvar
#' @examples
#' data(MS, package = "ER")
#' er <- ER(proteins ~ MS * cluster, data = MS[-1,])
#'
#' plsMod <- pls(er, 'MS', 6, validation = "CV",
#'               type = "interleaved", length.seg=25, shave = TRUE)
#' # Error as a function of remaining variables
#' plot(plsMod$shave)
#' # Selected variables for minimum error
#' with(plsMod$shave, colnames(X)[variables[[min.red+1]]])
#'
#' \donttest{
#' plsMod <- pls(er, 'MS', 5, validation = "LOO",
#'               type = "interleaved", length.seg=25, jackknife = TRUE)
#' colSums(plsMod$classes == as.numeric(MS$MS[-1]))
#' # Jackknifed coefficient P-values (sorted)
#' plot(sort(plsMod$jack[,1,1]), pch = '.', ylab = 'P-value')
#' abline(h=c(0.01,0.05),col=2:3)
#'
#' scoreplot(plsMod)
#' }
#'
#' @rdname pls
#' @export
pls <- function(er, ...){
  UseMethod("pls")
}
# setGeneric("pls")

#' @rdname pls
#' @method pls ER
#' @export
pls.ER <- function(er, effect, ncomp, newdata = NULL, er2, validation, jackknife = NULL, shave = NULL, df.used = NULL, ...){
  if(!missing(er2)){
    data <- data.frame(X = I(er$ER.values[[effect]]),
                       y = I(er2$ER.values[[effect]]))
    lda    <- NULL
  } else {
    if(length(effect) == 1){
      data <- data.frame(X = I(er$ER.values[[effect]]),
                         y = I(er$symbolicDesign[[effect]]),
                         Yd = I(model.matrix(~y-1,data.frame(y=er$symbolicDesign[[effect]]))))
    } else { # User supplied contrast
      data <- data.frame(X = I(er$ER.values[[effect[[1]]]]),
                         y = effect[[2]],
                         Yd = effect[[2]])
    }
  }

  if(is.null(df.used)){
    df.used <- er$df.used
  }
  jack   <- ifelse(is.null(jackknife), FALSE, TRUE)
  shaved <- ifelse(is.null(shave), FALSE, TRUE)
  if(!is.null(newdata)){
    if(missing(er2)){
      plsMod <- plsr(Yd ~ X, ncomp = ncomp, data = data, ...)
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
      if(missing(er2)){
        plsMod <- plsr(Yd ~ X, ncomp = ncomp, data = data, validation = validation, jackknife = jack, ...)
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
        if(missing(er2)){
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
  object$classes <- lda
  object$data <- data
  # object <- list(classes = lda, data = data, pls = plsMod)
  if(jack)
    object$jack <- jt
  if(shaved)
    object$shave <- sh
  class(object) <- c('ERpls','mvr','list')
  object
}

#' @export scores
pls::scores

#' @export scoreplot
pls::scoreplot

#' @export loadings
pls::loadings

#' @export loadingplot
pls::loadingplot

#' @export R2
pls::R2

#' @export mvrValstats
pls::mvrValstats

#' @export explvar
pls::explvar
