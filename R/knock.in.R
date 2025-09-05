#' Knock-in and knock-out of effect(s)
#'
#' @description
#' Convenience functions to apply to GEM objects to isolate/extract (knock-in) one or more effects
#' (and possibly residuals) or to remove (knock-out) one or more effects.
#'
#'
#' @param object GEM object.
#' @param effect Name or number of effect (character or numeric vector).
#' @param residuals Logical indicating if residuals should be added (default = TRUE).
#'
#' @return A \code{data.frame} of ER values where effects have been knocked in
#' (isolated/extracted from data) or knocked out (removed from data).
#' @export
#'
#' @examples
#' data(MS, package = "gemR")
#'
#' # Subset to reduce runtime in example
#' MS$proteins <- MS$proteins[,20:70]
#'
#' gem <- GEM(proteins ~ MS * group, data = MS)
#'
#' # Extract interaction between 'MS' and 'group
#' ER.isolated <- knock.in(gem, 'MS:group')
#'
#' # Remove main effect of 'group'
#' ER.cleaned <- knock.out(gem, 'group')
knock.in <- function(object, effect, residuals = TRUE){
  ER <- object$effects[[effect[1]]]
  if(n.eff<-length(effect) > 1){
    for(j in 2:n.eff)
      ER <- ER + object$effects[[effect[j]]]
  }
  if(residuals)
    ER <- ER + object$residuals
  return(ER)
}

#' @rdname knock.in
#' @export
knock.out <- function(object, effect, residuals = TRUE){
  ER <- object$fitted.values
  for(j in 1:length(effect))
    ER <- ER - object$effects[[effect[j]]]
  if(residuals)
    ER <- ER + object$residuals
  return(ER)
}
