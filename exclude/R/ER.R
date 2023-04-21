#' Effect + Residual Modelling
#'
#' @param formula a model formula specifying features and effects.
#' @param data a \code{data.frame} containing response variables (features) and design factors or other groupings/continuous variables.
#' @param x Object of class \code{ER}.
#' @param y Response name or number.
#' @param what What part of ER to plot; \code{raw} data (default), \code{fits}, \code{residuals} or
#' a named model effect (can be combined with 'mean', see \code{Examples}).
#' @param col Color of points, defaults to grouping. Usually set to a factor name.
#' @param pch Plot character of points, defaults to 1. Usually set to a factor name.
#' @param model.line Include line indicating estimates, default = TRUE. Can be an effect name.
#' @param xlab X label (\code{character})
#' @param ylab Y label (\code{character})
#' @param main Main title, defaults to \code{y} with description from \code{what}.
#' @param ... Additional arguments to \code{plot}
#' @param object ER object.
#' @param variable Numeric for selecting a variable for extraction.
#' @importFrom graphics lines plot
#' @importFrom stats coef lm model.matrix na.omit predict pt qt sd
#'
#' @return An object of class \code{ER} containing effects, ER values,
#' fitted values, residuals, features, coefficients, dummy design, symbolic design, dimensions,
#' highest level interaction and feature names.
#'
#' @references E.F. Mosleth et al. (2020) <doi:10.1016/B978-0-12-409547-2.14882-6>. Comprehensive Chemometrics, 2nd edition; Brown, S., Tauler, R., & Walczak, B. (Eds.). Chapter 4.22. Analysis of Megavariate Data in Functional Omics. Elsevier.
#'
#' @export
#'
#' @examples
#' data(MS, package = "ER")
#' er <- ER(proteins ~ MS * cluster, data = MS)
#' plot(er)                                           # Raw data, first feature
#' plot(er,2)                                         # Raw data, numbered feature
#' plot(er,'Q76L83', col='MS', pch='cluster')         # Selected colour and plot character
#' plot(er,'Q76L83', what='MS', model.line='cluster') # Selected model.line
#' \donttest{
#'   # Example compound plot
#'   old.par <- par(c("mfrow", "mar"))
#'   on.exit(par(old.par))
#'   par(mfrow = c(3,3), mar = c(2,4,4,1))
#'   plot(er,'Q76L83')                                  # Raw data, named feature
#'   plot(er,'Q76L83', what='fits')                     # Fitted values
#'   plot(er,'Q76L83', what='residuals')                # Residuals
#'   plot(er,'Q76L83', what='mean MS')                  # Mean levels
#'   plot(er,'Q76L83', what='mean cluster')             # ----||----
#'   plot(er,'Q76L83', what='mean MS:cluster')          # ----||----
#'   plot(er,'Q76L83', what='MS')                       # ER values
#'   plot(er,'Q76L83', what='cluster')                  # --------||---------
#'   plot(er,'Q76L83', what='MS:cluster')               # --------||---------
#' }
#'
#' # In general there can be more than two, effects, more than two levels, and continuous effects:
#' # MS$three <- factor(c(rep(1:3,33),1:2))
#' # er3    <- ER(proteins ~ MS * cluster + three, data = MS)
ER <- function(formula, data){
  # Handle formulas
  mf <- match.call(expand.dots = FALSE)                         # Bokholderi på input-navn
  mf[[1L]] <- as.name("model.frame")                            # Bytte fra ER til model.frame
  old.opt <- options(contrasts = c("contr.sum","contr.poly"))   # Bytte fra treatment til sum-to-zero
  mf <- eval(mf, parent.frame())                                # Evaluerer input til ER som data.frame utenfor ER-funksjonen
  mm <- model.matrix(mod <- lm(mf))                             # Koder om faktorer til dummy, lager matrise av prediktorer
  mmAssign <- attr(mm,'assign')
  options(old.opt)                                              # Tilbake til tidligere parametrisering
  factorCombinations <- attr(attr(mf, "terms"),"factors")       # Tar ut effektnavn og samspillseffektnavn
  variables <- mf[[1]]                                          # Henter ut responsen(e)
  # Dimensions
  N     <- dim(data)[1]
  p     <- dim(variables)[2]
  N.eff <- dim(factorCombinations)[2]
  dims  <- c(N=N,p=p,N.eff=N.eff)
  # Aggregate mean values and mean effect values
  mr     <- do.call(interaction,mf[-1]) # All subgroups
  u      <- levels(mr)
  design <- data.frame(mr)
  residuals    <- residuals(mod)
  coefficients <- coef(mod)
  effects      <- list()
  corrE        <- list()
  yhat         <- matrix(0.0, N, p)
  for(i in 1:(N.eff+1)){
    effects[[i]] <- mm[,mmAssign == (i-1), drop=FALSE] %*%
      coefficients[mmAssign == (i-1),, drop=FALSE]
    corrE[[i]] <- effects[[i]] + residuals
    yhat <- yhat + effects[[i]]
  }
  names(effects) <- names(corrE) <- c("(Intercept)",colnames(factorCombinations))
  dimnames(yhat) <- dimnames(variables)
  des <- mf[-1]
  interactions <- setdiff(colnames(factorCombinations),colnames(des))
  if(length(interactions) > 0){
    for(i in 1:length(interactions)){
      # if(all(lapply(mf[strsplit(interactions[i],':')[[1]]], class)=="numeric")){
      #   des <- cbind(des, apply(mf[strsplit(interactions[i],':')[[1]]],1,prod))
      # } else {
        des <- cbind(des, do.call(interaction,mf[strsplit(interactions[i],':')[[1]]]))
      # }
      colnames(des)[ncol(des)] <- interactions[i]
    }
  }
  featureNames <- colnames(variables)
  names(featureNames) <- featureNames
  ret <- list(effects = effects,
              ER.values = corrE,
              fitted.values = yhat,
              residuals = residuals,
              features = variables,
              coefficients = coefficients,
              design = mm,
              symbolicDesign = des,
              dimensions = dims,
              highestLevel = design[[1]],
              featureNames = featureNames,
              df.used = dim(design)[1] - mod$df.residual)
  class(ret) <- c("ER", "list")
  ret
}

#' @export
#' @rdname ER
plot.ER <- function(x, y = 1, what = "raw", col = NULL, pch = NULL,
                      model.line = !(what %in% c("fits","residuals")),
                      ylab = '', xlab = '', main = NULL, ...){
  # Set defaults
  if(is.null(col))
    col <- match(x$highestLevel,unique(x$highestLevel))
  if(length(col)==1 && col %in% colnames(x$symbolicDesign))
    col <- match(x$symbolicDesign[,col],unique(x$symbolicDesign[,col]))
  if(is.null(pch))
    pch <- 1
  if(length(pch) == 1 && pch %in% colnames(x$symbolicDesign))
    pch <- match(x$symbolicDesign[,pch],unique(x$symbolicDesign[,pch]))
  has.main <- !is.null(main)

  # Extract for plotting
  if(what %in% c("raw", "fits", "residuals")){
    dataLine   <- x$fitted.values[,y]
    dataPoints <- switch(what,
      raw  = x$features[,y],
      fits = x$fitted.values[,y],
      residuals = x$residuals[,y])
    if(is.null(main))
      main <- c(x$featureNames[y],
                switch(what, raw = "raw data",
                       fits = "fitted values",
                       residuals = "residuals"))
  } else {
    if(what %in% colnames(x$symbolicDesign)){
      dataPoints <- x$ER.values[[what]][,y]
      dataLine   <- x$ER.values[[what]][,y]
      if(is.null(main))
        main <- c(paste(x$featureNames[y],", ",names(x$ER.values[what]), sep=""),
                  "ER values")
    } else {
      if(what %in% paste("mean", colnames(x$symbolicDesign))){
        what <- strsplit( what,"mean ")[[1]][2]
        dataPoints <- x$ER.values[[what]][,y]
        dataLine   <- NULL
        main <- c(paste(x$featureNames[y],", ",names(x$ER.values[what]), sep=""),
                  "effect estimates")
      } else {
        stop('Invalid value for "what" (neither "raw", "fits", "residuals", nor a model effect)')
      }
    }
  }
  # Plot points
  plot(dataPoints, col = col, pch = pch, xlab = xlab, ylab = ylab, main = main, ...)
  # ... and line
  if(is.logical(model.line) && model.line){
    lines(dataLine)
  } else {
    if(!is.logical(model.line)){
      if(model.line %in% colnames(x$symbolicDesign)){
        dataLine   <- x$ER.values[[model.line]][,y]
        lines(dataLine)
      }
      if(model.line == "raw"){
        dataLine   <- x$fitted.values[,y]
        lines(dataLine)
      }
    }
  }
}

#' @export
#' @rdname ER
tableER <- function(object, variable){
  with(object, cbind(symbolicDesign, data.frame(lapply(ER.values,function(i)i[,variable]), check.names = FALSE),
                     fits = fitted.values[,variable], resid = residuals[,variable], resp = features[,variable]))
}
