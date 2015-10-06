#' Compute the Fst of haplotype between 2 pops.
#'
#' @param x A vector of allele (0 or 1) value in two pops with same number of individual.
#'
#' @export
Fst2Pop <- function( x ) {
  return((var(x)-var(x[1:(length(x)/2)])/2-var(x[(length(x)/2+1):length(x)])/2)/max(.Machine$double.eps,var(x) ) )
}

#' Compute the Fst of haplotype with OF formula
#'
#' @param freq A matrix of size nbLocus x nbPop
#' @param q A vector or list of size nbPop which contains the normalized size of pop
#'
#' @export
FstOF <- function( freq, q ) {

  sigma2_S = apply( freq * (1 - freq ), MARGIN = 1, FUN = function(x){return(sum(q*x))} )
  sigma2_T = apply( freq , MARGIN = 1, FUN = function(x){return(sum(q*x)*(1-sum(q*x)))} )

  return(1 - sigma2_S / sigma2_T)

}

#' Fdr control with benjamini hochberg
#'
#'
#' @export
fdrControl <- function( p.values, alphaValues ) {
  res = list()
  for (alpha in alphaValues) {
    # expected FDR = alpha
    # print(paste("expected FDR:", alpha))
    # return a list of candidates with expected FDR = alpha
    k = sum( sort(p.values) <  alpha * (1:(length(p.values)))/length(p.values))
    candidates = which( p.values <  alpha * k/length(p.values) )
    res[[paste("alpha",alpha,sep="")]] = candidates

  }
  return( res )

}


#' Fdr control with benjamini hochgerg. This function
#' compute the true TP rate and FDR rate.
#'
#'
#' @export
trueTPandFDR <- function( p.values, alphaValues, outliers ) {
  nbOutlier = length(outliers)
  nbTrue = length(p.values) - length(outliers)
  observedFdr = c()
  observedTP = c()
  for (alpha in alphaValues) {
    # expected FDR = alpha
    # print(paste("expected FDR:", alpha))
    # return a list of candidates with expected FDR = alpha
    k = sum( sort(p.values) <  alpha * (1:(length(p.values)))/length(p.values))
    candidates = which( p.values <  alpha * k/length(p.values) )

    # estimated FDR and True Positif rate
    estimated.FDR = length(which(!(candidates %in% outliers))) / length(candidates)
    observedFdr = c(observedFdr,estimated.FDR)
    estimated.TP = length(which(candidates %in% outliers)) / length(outliers)
    observedTP = c(observedTP,estimated.TP)
    # print(paste("Observed FDR:", estimated.FDR, "Power:", estimated.TP))
  }
  return( data.frame(expectedFdr = alphaValues, observedFdr = observedFdr, observedTP = observedTP) )

}


#' The allele frequency spectrum
#'
#' @param genotype A matrix of genotype of size nbIndiv x nbLocus
#' @param allel 0 or 1 for allele which you want the spectrum
#'
#' @export
SNPspectrum <- function( genotype, allele = 1 ) {
  return(table(apply( data$genotype, MARGIN = 2, FUN = function(x) sum(x==allele))) )
}

#' The folded allele frequency spectrum
#'
#' @param genotype A matrix of genotype of size nbIndiv x nbLocus
#'
#' @export
SNPfoldedspectrum <- function( genotype ) {
  return(table( apply( data$genotype, MARGIN = 2, FUN = function(x) min(sum(x==1),sum(x==0))) ) )
}

#' Power and Fdr
#'
#'
#' @export
power_fdr <- function( p.value, outlier  ) {
  aux = sort(p.value, index.return = TRUE)$ix
  return(data.frame( fdr = sapply(0:length(p.value), function(k) { sum(!(aux[0:k] %in% outlier)) / k } ),
                     power = sapply(0:length(p.value), function(k) { sum(aux[0:k] %in% outlier) / length(outlier) } ) ) )
}




