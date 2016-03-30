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

test_that("sample_correlated_X", {
  c = 0.5
  mu = 0.0
  sd = 1.0
  Y = rnorm(1000,0.5,0.5)
  X = sample_correlated_X(Y,c,mu,sd)

  expect_less_than(abs(cor(X,Y)-c),0.001)
  expect_less_than(abs(mean(X)-mu),0.001)
  expect_less_than(abs(sd(X)-sd),0.001)

  c = 0.02
  mu = 10
  sd = 0.05
  Y = rnorm(1000,0.5,0.5)
  X = sample_correlated_X(Y,c,mu,sd)

  expect_less_than(abs(cor(X,Y)-c),0.001)
  expect_less_than(abs(mean(X)-mu),0.001)
  expect_less_than(abs(sd(X)-sd),0.001)

  qqplot(X,rnorm(1000,mu,sd))
  abline(0,1)

})
