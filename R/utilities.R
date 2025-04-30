# PCA for internal use
.pca <- function(X, ncomp, scale=FALSE, proj=NULL, ...){
  X <- as.matrix(unclass(X))
  if(!inherits(X,'matrix'))
    stop("'X' must be a matrix")
  if(missing(ncomp))
    ncomp <- min(c(nrow(X)-1, ncol(X)))
  y   <- structure(matrix(rnorm(nrow(X)), ncol=1), rownames = rownames(X))
  dat <- data.frame(y=y, X = I(X))
  mod <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale, ...)
  mod$explvar <- mod$Xvar/mod$Xtotvar*100
  attr(mod$loadings, 'explvar') <- attr(mod$scores, 'explvar') <- mod$explvar
  if(!is.null(proj)){
    mod$projected <- (proj-rep(mod$Xmeans, each=nrow(proj))) %*% mod$loadings
  }
  mod$singulars <- sqrt(mod$Xvar)
  mod$call <- match.call()
  class(mod) <- c('pca')
  return(mod)
}
