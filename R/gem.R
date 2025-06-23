#' General Effect Modelling
#'
#' @aliases GEM plot.GEM tableGEM summary.GEM gemR
#' @param formula a model formula specifying features and effects.
#' @param data a \code{data.frame} containing response variables (features) and design factors or other groupings/continuous variables.
#' @param contrasts a \code{character} containing the primary contrasts for use with \code{mixlm::lm} (default = "contr.sum").
#' @param add_residuals Logical indicating if residuals should be added to the ER values (default = TRUE).
#' @param x Object of class \code{GEM}.
#' @param y Response name or number.
#' @param what What part of GEM to plot; \code{raw} data (default), \code{fits}, \code{residuals} or
#' a named model effect (can be combined with 'effect', see \code{Examples}).
#' @param col Color of points, defaults to grouping. Usually set to a factor name or a column name in the input data with custom colours.
#' @param pch Plot character of points, defaults to 1. Usually set to a factor name or a column name in the input data with custom symbols
#' @param model.line Include line indicating estimates, default = TRUE. Can be an effect name.
#' @param ylim Y axis limits (\code{numeric}, but defaults to NULL)
#' @param xlab X label (\code{character})
#' @param ylab Y label (\code{character})
#' @param main Main title, defaults to \code{y} with description from \code{what}.
#' @param extended Extended output in summary (default = TRUE).
#' @param df Show degrees of freedom in summary (default = FALSE).
#' @param digits \code{integer} number of digits for printing.
#' @param object GEM object.
#' @param variable Numeric for selecting a variable for extraction.
#' @param ... Additional arguments to \code{plot}
#' @importFrom graphics lines plot
#' @importFrom stats coef model.matrix na.omit predict pt qt sd aggregate model.frame
#' @importFrom utils data
#' @importFrom mixlm lm
#' @importFrom HDANOVA hdanova
#' @importFrom lme4 getME
#'
#' @return \code{GEM} returns an object of class \code{GEM} containing effects, ER values (effect + residuals),
#' fitted values, residuals, features, coefficients, dummy design, symbolic design, dimensions,
#' highest level interaction and feature names.
#'
#' @references
#' * Mosleth et al. (2021) Cerebrospinal fluid proteome shows disrupted neuronal development in multiple sclerosis. Scientific Report, 11,4087. <doi:10.1038/s41598-021-82388-w>
#'
#' * E.F. Mosleth et al. (2020). Comprehensive Chemometrics, 2nd edition; Brown, S., Tauler, R., & Walczak, B. (Eds.). Chapter 4.22. Analysis of Megavariate Data in Functional Omics. Elsevier. <doi:10.1016/B978-0-12-409547-2.14882-6>
#' @seealso Analyses using \code{GEM}: \code{\link{elastic}}, \code{\link{pca}}, \code{\link{sca}}, \code{\link{neuralnet}}, \code{\link{pls}}.
#' Confidence interval plots: \code{\link{confints}}. Convenience knock-in and knock-out of effects: \code{\link{knock.in}}.
#' @examples
#' ## Multiple Sclerosis
#' data(MS, package = "gemR")
#' # Subset to reduce runtime in example
#' MS$proteins <- MS$proteins[,20:70]
#'
#' gem <- GEM(proteins ~ group * MS, data = MS)
#' print(gem)
#' summary(gem)                                    # Summary of GEM
#' plot(gem)                                       # Raw data, first feature
#' plot(gem,2)                                     # Raw data, numbered feature
#' plot(gem,'Q76L83', col='MS', pch='group')       # Selected colour and plot character
#' plot(gem,'Q76L83', what='effect MS',
#'      model.line='effect group')                 # Comparison of factors (points and lines)
#' print(effs <- colnames(gem$symbolicDesign))     # Inspect factor names
#' eeffs <- paste0("effect ", effs)
#' # Example compound plot
#' old.par <- par(mfrow = c(3,3), mar = c(2,4,4,1))
#' plot(gem,'Q76L83')                         # Raw data, named feature
#' plot(gem,'Q76L83', what='fits')            # Fitted values
#' plot(gem,'Q76L83', what='residuals')       # Residuals
#' plot(gem,'Q76L83', what=eeffs[1])           # Effect levels
#' plot(gem,'Q76L83', what=eeffs[2])           # ----||----
#' plot(gem,'Q76L83', what=eeffs[3])           # ----||----
#' plot(gem,'Q76L83', what=effs[1])            # ER values
#' plot(gem,'Q76L83', what=effs[2])            # --------||---------
#' plot(gem,'Q76L83', what=effs[3])            # --------||---------
#' par(old.par)
#'
#' # Complete overview of GEM
#' tab <- tableGEM(gem, 1)
#'
#' # In general there can be more than two, effects, more than two levels, and continuous effects:
#' MS$three <- factor(c(rep(1:3,33),1:2))
#' gem3    <- GEM(proteins ~ MS * group + three, data = MS)
#'
#'
#' ## Candy assessment
#' data(candies, package = "HDANOVA")
#' gemC <- GEM(assessment ~ assessor*candy, data=candies)
#'
#' # Permutation testing
#' gemC <- permutation(gemC)
#' summary(gemC)
#'
#' # GEM-SCA with ellipsoids in score plots
#' gemSCA <- sca(gemC)
#' scoreplot(gemSCA, factor="candy", ellipsoids="confidence")
#'
#' # GEM-PCA with group colours
#' gemPCA <- pca(gemC)
#' scoreplot(gemPCA, factor="candy",
#'   gr.col=gemPCA$symbolicDesign$candy)
#'
#'
#' ## Lactobacillus
#' data(Lactobacillus, package = "gemR")
#' # Subset to reduce runtime in example
#' Lactobacillus$proteome <- Lactobacillus$proteome[,50:100]
#'
#' gemLac <- GEM(proteome ~ strain * growthrate, data = Lactobacillus)
#' print(gemLac)
#' plot(gemLac)                            # Raw data, first feature
#' plot(gemLac,2)                          # Raw data, numbered feature
#' plot(gemLac,'P.LSA0316', col='strain',
#'     pch='growthrate')                   # Selected colour and plot character
#' plot(gemLac,'P.LSA0316', what='strain',
#'     model.line='growthrate')            # Selected model.line
#'
#'
#' \donttest{ # Don't run this example, it takes too long
#'   ## Diabetes
#'   data(Diabetes, package = "gemR")
#'   gemDia <- GEM(transcriptome ~ surgery * T2D, data = Diabetes)
#'   print(gemDia)
#'   plot(gemDia)                            # Raw data, first feature
#'   plot(gemDia,2)                          # Raw data, numbered feature
#'   plot(gemDia,'ILMN_1720829', col='surgery',
#'       pch='T2D')                          # Selected colour and plot character
#' }
#'
#' @export
GEM <- function(formula, data, contrasts = "contr.sum", add_residuals = TRUE, ...){
  # Use HDANOVA as the work-horse before renaming elements
  HD <- hdanova(formula, data = data, contrasts = contrasts, ...)
  HD$models <- HD$models[1]
  HD$anovas <- HD$anovas[1]

  # Dimensions
  dims <- c(N = ncol(HD$Y), p = nrow(HD$Y), eff = length(HD$LS))

  ER.intercept <- matrix(colMeans(HD$Y), byrow = TRUE, nrow=nrow(HD$Y), ncol=ncol(HD$Y))
  if(add_residuals)
    ER.intercept <- ER.intercept + HD$residuals
  dimnames(ER.intercept) <- dimnames(HD$Y)

  HD$effectsOrig <- HD$effects
  HD$effects <- HD$LS
  if(add_residuals){
    HD$ER.values <- c("(Intercept)"=list(ER.intercept), HD$more$LS_aug)
  } else {
    HD$ER.values <- c("(Intercept)"=list(ER.intercept), HD$LS)
  }
  HD$fitted.values <- HD$Y-HD$residuals
  HD$features <- HD$Y
  HD$design <- HD$X
  HD$symbolicDesign <- HD$model.frame[-1]
  HD$dimensions <- dims
  HD$highestLevel <- data.frame(mf=apply(sapply(HD$model[-1], as.character), 1, paste0, collapse="."))[[1]]
  HD$featureNames <- colnames(HD$Y)
  names(HD$featureNames) <- HD$featureNames
  HD$df.used <- dim(HD$Y)[1] - HD$dfNum["Residuals"]
  HD$data <- data
  HD$call <- match.call()

  class(HD) <- c("GEM", "hdanova", "list")
  HD
}

#' @export
#' @rdname GEM
print.GEM <- function(x, ...){
  cat("Effect + Residual Model\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @export
#' @rdname GEM
plot.GEM <- function(x, y = 1, what = "raw", col = NULL, pch = NULL,
                     model.line = (what %in% c("raw")), ylim = NULL,
                     ylab = '', xlab = '', main = NULL, ...){
  # Set defaults
  if(is.null(col))
    col <- match(x$highestLevel,unique(x$highestLevel))
  if(length(col) == 1 && col %in% colnames(x$symbolicDesign)){
    col <- match(x$symbolicDesign[,col],unique(x$symbolicDesign[,col]))
  } else {
    if(length(col) == 1 && col %in% colnames(x$data)){
      col <- x$data[,col]
    }
  }
  if(is.null(pch))
    pch <- 1
  if(length(pch) == 1 && pch %in% colnames(x$symbolicDesign)){
    pch <- match(x$symbolicDesign[,pch],unique(x$symbolicDesign[,pch]))
  } else {
    if(length(pch) == 1 && pch %in% colnames(x$data)){
      pch <- x$data[,pch]
    }
  }
  has.main <- !is.null(main)
  dataLine   <- NULL

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
      dataLine   <- x$effects[[what]][,y]
      if(is.null(main))
        main <- c(paste(x$featureNames[y],", ",names(x$ER.values[what]), sep=""),
                  "ER values")
    } else {
      if(what %in% paste("effect", colnames(x$symbolicDesign))){
        what <- strsplit( what,"effect ")[[1]][2]
        #        dataPoints <- x$ER.values[[what]][,y]
        #        dataPoints <- tapply(dataPoints, x$symbolicDesign[[what]],mean)[x$symbolicDesign[[what]]]
        dataPoints <- x$effects[[what]][,y]
        dataLine   <- NULL
        main <- c(paste(x$featureNames[y],", ",names(x$ER.values[what]), sep=""),
                  "effect estimates")
      } else {
        stop('Invalid value for "what" (neither "raw", "fits", "residuals", nor a model effect)')
      }
    }
  }
  # Prepare line
  if(is.logical(model.line) && model.line){
    # lines(dataLine)
  } else {
    if(!is.logical(model.line)){
      if(model.line %in% colnames(x$symbolicDesign)){
        dataLine   <- x$ER.values[[model.line]][,y]
        # lines(dataLine)
      }
      if(model.line %in% paste("effect", colnames(x$symbolicDesign))){
        model.line <- strsplit( model.line,"effect ")[[1]][2]
        dataLine <- x$effects[[model.line]][,y]
        # lines(dataLine)
      }
      if(model.line == "raw"){
        dataLine <- x$fitted.values[,y]
        # lines(dataLine)
      }
      if(model.line %in% paste("mean", colnames(x$symbolicDesign))){
        model.line <- strsplit( model.line,"mean ")[[1]][2]
        dataLine <- aggregate(x$data[as.character(x$call[[2]][[2]])],by=list(x$symbolicDesign[[model.line]]),mean)[x$symbolicDesign[[model.line]],-1]
      }
      if(model.line %in% paste("mean", colnames(x$data))){
        model.line <- strsplit( model.line,"mean ")[[1]][2]
        dataLine <- aggregate(unclass(x$data[[as.character(x$call[[2]][[2]])]]),by=list(x$data[[model.line]]),mean)[x$data[[model.line]],-1][,y]
      }
    }
  }
  if(is.null(ylim) && !(is.logical(model.line)&&!model.line)){
    ylim = c(min(min(dataPoints),min(dataLine)),max(max(dataPoints),max(dataLine)))
  }
  # Plot points
  plot(dataPoints, col = col, pch = pch, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  # ... and lines
  if(!(is.logical(model.line)&&!model.line))
    lines(dataLine)
}

#' @export
#' @rdname GEM
tableGEM <- function(object, variable){
  with(object, cbind(symbolicDesign,
                     data.frame(lapply(ER.values,function(i)i[,variable]), check.names = FALSE),
                     fits = fitted.values[,variable],
                     resid = residuals[,variable],
                     resp = features[,variable]))
}

#' @export
#' @rdname GEM
summary.GEM <- function(object, extended=TRUE, df=FALSE, ...){
  dat <- data.frame(SSQ=object$ssq, "Expl.var"=object$explvar*100)
  colnames(dat) <- c("Sum.Sq.", "Expl.var.(%)")
  if(!is.null(object$permute)){
    pvals <- object$permute$pvalues
    pvals[pvals==0] <- 1/object$permute$permutations
    pv <- rep(NA,nrow(dat))
    names(pv) <- rownames(dat)
    pv[names(pvals)] <- pvals
    dat <- cbind(dat, "p-value"=pv)
  }
  mod <- "General Effect Modelling"
  x <- list(dat=dat, model=mod, fit.type=object$fit.type)
  if(extended){
    LS_REML <- "least squares"
    if(!inherits(object$models[[1]],"lm"))
      LS_REML <- ifelse(getME(object$models[[1]],"is_REML"), "REML", "ML")
    ss <- c("I","II","III")
    x$info <- paste0("SS type ", ss[object$SStype], ", ", object$coding, " coding, ",
                     ifelse(object$unrestricted, "unrestricted","restricted"), " model",
                     ", ", LS_REML, " estimation")
    if(!is.null(object$permute))
      x$info <- paste0(x$info, ", ", object$permute$permutations, " permutations")
  }
  if(df){
    x$dat <- cbind(x$dat, "df"=object$dfNum, "df.denom"=object$dfDenom, "err.term"=object$denom)
  }
  class(x) <- c('summary.hdanova')
  x
}

#' @export
#' @rdname GEM
print.summary.GEM <- function(x, digits=2, ...){
  cat(x$mod, "fitted using", x$fit.type, "\n")
  if(!is.null(x$info))
    cat("-", x$info, "\n")
  print(round(x$dat, digits))
  invisible(x$dat)
}
