#' An sNMF wrapper for a good integration to R
#'
#'
#' @export
RunsNMF <- function (genotype,
                     K,
                     ploidy,
                     alpha = 10,
                     max.iteration = 200,
                     tolerance = 1e-5) {

  # create tmp file
  genotype[is.na(genotype)] <- 9
  file.geno <- paste0(tempfile(),".geno")
  LEA::write.geno(genotype, file.geno)

  # run sNMF
  res <- LEA::snmf(file.geno,
                   K,
                   project = "new",
                   repetitions = 1,
                   alpha = alpha,
                   iterations = max.iteration,
                   tolerance = tolerance,
                   ploidy = ploidy,
                   entropy = FALSE)


  return(list(Q = LEA::Q(res,run = 1, K = K), G = LEA::G(res, run = 1, K = K)))

}
