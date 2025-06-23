#' @aliases elastic elastic.GEM coef.GEMglmnet
#' @name elastic
#' @title Elastic-net modeling of GEM objects.
#' @param gem Object of class \code{GEM}.
#' @param effect The effect to be used as response.
#' @param alpha The elasticnet mixing parameter.
#' @param newdata Optional new data matrix for prediction.
#' @param validation Optional validation parameters.
#' @param segments number of segments or list of segments (optional)
#' @param measure Type of performance summary, default = 'class' (see \code{\link[glmnet]{glmnet}})
#' @param family Type of model response, default = 'multinomial'.
#' @param ... Additional arguments for \code{\link[glmnet]{glmnet}}.
#'
#' @seealso Analyses using \code{GEM}: \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}, \code{\link{pls}}.
#' Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
#'
#' @return An object of class \code{GEMglmnet, cv.glmnet, list} containing the fitted Elastic-net model, classifications/predictions and data.
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom stats rnorm
#' @examples
#' ## Multiple Sclerosis data
#' data(MS, package = "gemR")
#' # Subset to reduce runtime in example
#' MS$proteins <- MS$proteins[,20:70]
#'
#' gem <- GEM(proteins ~ MS * group, data = MS)
#' elasticMod <- elastic(gem, 'MS', validation = "CV")
#' sum(elasticMod$classes == MS$MS)
#' plot(elasticMod)            # Model fit
#' plot(elasticMod$glmnet.fit) # Coefficient trajectories
#'
#' # Select all proteins with non-zeros coefficients
#' coefs     <- coef(elasticMod)
#' (selected <- names(which(coefs[,1] != 0)))
#'
#' \donttest{ # Time consuming due to many variables
#'   ## Diabetes data
#'   data(Diabetes, package = "gemR")
#'   gem.Dia <- GEM(transcriptome ~ surgery * T2D, data = Diabetes)
#'   elasticMod <- elastic(gem.Dia, 'T2D', validation = "LOO")
#' }
#' @export
elastic <- function(gem, ...){
  UseMethod("elastic")
}

#' @rdname elastic
#' @method elastic GEM
#' @export
elastic.GEM <- function(gem, effect, alpha = 0.5, newdata = NULL, validation, segments = NULL,
                          measure = measure, family = family, ...){
  if(length(effect) == 1){
    data <- data.frame(X = I(gem$ER.values[[effect]]),
                       y = gem$symbolicDesign[[effect]],
                       Yd = I(model.matrix(~y-1,data.frame(y=gem$symbolicDesign[[effect]]))))
  } else { # User supplied contrast
    data <- data.frame(X = I(gem$ER.values[[effect[[1]]]]),
                       y = I(effect[[2]]),
                       Yd = I(effect[[2]]))
  }

  if(missing(validation))
    validation <- "LOO"
  if(validation == "LOO"){
    cv <- 1:nrow(data)
  } else {
    if(is.null(segments)){
      k <- 10
      cv <- cvsegments(nrow(data), k, ...)
    } else {
      if(!is.list(segments)){
        if(is.numeric(segments))
          k <- segments
        cv <- cvsegments(nrow(data), k, ...)
      } else {
        cv <- segments
      }
    }
    cv <- unlist(lapply(1:length(cv), function(i) rep(i,length(cv[[i]]))))[order(unlist(cv))]
  }
  # glmnet.data     <- glmnet(data$X,data$y)
  if(is.factor(data$y))
    object     <- cv.glmnet(data$X,data$y, alpha=alpha, foldid = cv, grouped=FALSE, family = "multinomial")
  else
    object   <- cv.glmnet(data$X,data$y, alpha=alpha, foldid = cv, grouped=FALSE, family = "gaussian")
  #co         <- coef(object, s='lambda.min', exact=TRUE)
  #c.vector   <- as.numeric(co[[2]]); names(c.vector) <- rownames(co)
  #temp       <- sort(c.vector)
  #inds       <- which(c.vector!=0)

  if(!is.null(newdata)){
    if(is.factor(data$y))
      classEl <- factor(object$glmnet.fit$classnames)[apply(predict(object, newdata),1,which.max)]
    else
      classEl <- predict(object, newdata)
  } else {
    if(is.factor(data$y))
      object$classes <- factor(object$glmnet.fit$classnames)[apply(predict(object, data$X),1,which.max)]
    else
      object$preds <- predict(object, data$X)
  }

  object$data <- data
  class(object) <- c('GEMglmnet','cv.glmnet','list')
  object
}

#' @export
coef.GEMglmnet <- function(object, s='lambda.min', exact=TRUE, ...){
  do.call(cbind, lapply(glmnet::coef.glmnet(object, s=s, exact=exact, ...), as.matrix))
}
