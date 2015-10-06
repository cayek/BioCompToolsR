#' Compute true FDR
#'
#'
#' @export
trueFDR <- function( candidates, outlier ) {
  return(length(which(!(candidates %in% outliers))) / length(candidates))
}


#' Compute true power
#'
#'
#' @export
truePower <- function( candidates, outlier ) {
  return(length(which(candidates %in% outliers)) / length(outliers))
}

