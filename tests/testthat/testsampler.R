library(BioCompToolsR)
context("sampler testing")

test_that("admixture2popSimulationWrapper testing", {
  dataDirectory="/home/cayek/Téléchargements/"
  nbIndiv = 200
  nbLocus = 30000
  prop1 = 0.95
  prop2 = 1 - prop1
  Nm1 = 20
  Nm2 = 0.1

  data2 = admixture2popSimulationWrapper(dataDirectory,nbIndiv, prop2 * nbLocus ,d=1,Nm=Nm2,k=Inf, force = TRUE, admixe = FALSE, write = FALSE)


})
