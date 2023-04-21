#' @aliases elastic elastic.GEM
#' @name elastic
#' @title Elastic-net modeling of GEM objects.
#' @param gem Object of class \code{GEM}.
#' @param effect The effect to be used as response.
#' @param alpha The elasticnet mixing parameter.
#' @param newdata Optional new data matrix for prediction.
#' @param validation Optional validation parameters.
#' @param segments number of segments or list of segments (optional)
#' @param measure Type of performance summary, default = 'class' (see \code{\link{glmnet}})
#' @param family Type of model response, default = 'multinomial'.
#' @param ... Additional arguments for \code{pls::cvsegments}.
#'
#' @seealso \code{\link{GEM}}, \code{\link{pls}} and \code{\link{confints}}.
#'
#' @importFrom glmnet cv.glmnet
#' @examples
#' ## Multiple Sclerosis data
#' data(MS, package = "gemR")
#' gem <- GEM(proteins ~ MS * cluster, data = MS)
#' elasticMod <- elastic(gem, 'MS', validation = "CV")
#' sum(elasticMod$classes == MS$MS)
#' plot(elasticMod)            # Model fit
#' plot(elasticMod$glmnet.fit) # Coefficient trajectories
#'
#' # Select all proteins with non-zeros coefficients
#' coefs     <- coef(elasticMod,s='lambda.min',exact=TRUE)
#' (selected <- rownames(coefs[[1]])[unique(unlist(lapply(coefs,
#'                       function(x)which(as.vector(x) != 0))))][-1])
#'
#' \donttest{
#' ## Diabetes data
#' data(Diabetes, package = "gemR")
#' gem.Dia <- GEM(transcriptome ~ surgery * T2D, data = Diabetes)
#' elasticMod <- elastic(gem.Dia, 'T2D', validation = "LOO")
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
  cv.glm     <- cv.glmnet(data$X,data$y, alpha=alpha, foldid = cv, grouped=FALSE, type.measure='class', family = "multinomial")
  co         <- coef(cv.glm, s='lambda.min', exact=TRUE)
  c.vector   <- as.numeric(co[[2]]); names(c.vector) <- rownames(co)
  temp       <- sort(c.vector)
  inds       <- which(c.vector!=0)

  if(!is.null(newdata)){
    classEl <- factor(cv.glm$glmnet.fit$classnames)[apply(predict(cv.glm, newdata),1,which.max)]
  } else {
    classEl <- factor(cv.glm$glmnet.fit$classnames)[apply(predict(cv.glm, data$X),1,which.max)]
  }

  object <- cv.glm
  object$classes <- classEl
  object$data <- data
  # object <- list(classes = classEl, data = data, glmnet = cv.glm)
  class(object) <- c('GEMglmnet','cv.glmnet','list')
  object
}
