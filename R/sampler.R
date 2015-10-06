#' Sample 2 pop simulated with ms
#'
#'
#' @export
admixture2popSimulationWrapper <- function(dataDirectory,nbIndiv,nbLocus,Nm,d,k = 1.9,force=FALSE, OnlyX=FALSE, admixe = TRUE, write = TRUE, snpsRate = 0.0) {

  # set wd
  wd.old = getwd()

  print("use : ... dataDirectory nbIndiv nbLocus Nm(migration coeficient) d(1 : haploir 2 : diploid) k(cline gradient) force")


  #Haploid or diploid
  if ( d!=1 && d != 2 ) {
    print("d = 1 or d = 2")
    return()
  }

  #test if data already exist
  outputName = paste( "adMixture_K",2,"_Nm",Nm,"_n",nbIndiv,"_L",nbLocus,"_d",d,"_k",k,"_OnlyX", OnlyX,sep="" )
  simulationDirectory=paste( dataDirectory,"/",outputName,sep="" )
  resFile = paste( simulationDirectory,"/simulation.res",sep="" )
  if ( file.exists(file = resFile ) && !force ) {
    print("Data already simulated, we return it!!")
    load(resFile)
    return(res)
  }

  #create directory
  system(paste("rm -rf ", simulationDirectory,sep=""))
  system(paste("mkdir ", simulationDirectory,sep=""))

  setwd(paste(simulationDirectory,sep=""))


  ############################
  #Create 2 island model data#
  ############################

  print("-------Creating islands")

  nbHaplotype = nbIndiv * d

  #launch ms
  #see file:///home/cayek/T%C3%A9l%C3%A9chargements/Thursday_Lab_keys_for_instructors.pdf for more details about ms
  simulate2Islands.py = system.file("python","simulate2Islands.py",package = "BioCompToolsR")
  tmpfile = tempfile()
  system( paste(simulate2Islands.py, as.integer(nbHaplotype), as.integer(nbLocus), Nm, tmpfile, sep=" " ) )
  #Rmk : ms doesnt sample SNIPs

  #load data
  genotype = read.table(file = tmpfile)

  #Generate coordinates
  coord = cbind(sort(c(rnorm(nbHaplotype / 2 , -2, 1), rnorm(nbHaplotype / 2, 2, 1))), runif(nbHaplotype) )# runif because otherwise tess bug...
  if( OnlyX) {

    coord[,2]=0.0

  }
  if( d == 2 ) {
    coord[which((1:nbHaplotype)%%2 == 0 ),] = coord[which((1:nbHaplotype)%%2 == 1 ),]
  }

  dat <- data.frame(coord, genotype)

  #We can remove locus whiwh are not SNIPs
  if(  length(which( apply(genotype,MARGIN=2,sum) < snpsRate * nbHaplotype )) > 0) {
    nbLocus = nbLocus - length(which( apply(genotype,MARGIN=2,sum) < snpsRate * nbHaplotype  ))
    genotype = genotype[ ,- which( apply(genotype,MARGIN=2,sum) < snpsRate * nbHaplotype  ) ]
  }
  dat <- data.frame(coord, genotype)


  #plot pc without nnt SNIPs locus
  pc = prcomp(genotype, scale = FALSE) # Remark when setting T : there are some outliers but if we set F they disappear !!
  plot( pc$x , pch = 19, cex = 2, col = rep(c("blue2","orange"), each = nbHaplotype /2))

  #############################
  #Admixture of the 2 islands#
  #############################
  if (admixe) {
    print("-------Admixing islands")

    # A function for the shape of a cline
    sigmoid <- function(x){ 1/(1 + exp(-x)) }
    p1 <- sigmoid( k * coord[,1] )
    #remark : we only take x


    curve(sigmoid(k * x),-10.,10.)

    # Our admixed genotypes are build from a 2 island population model
    admixed.genotype <- matrix(NA, ncol = nbLocus, nrow = nbHaplotype)
    aux=matrix(NA, ncol = nbLocus, nrow = nbHaplotype)
    for (i in 1:(nbHaplotype/2) ){
      aux[i,]=sample( c(0,1), nbLocus , replace = TRUE, prob = c(p1[i], 1 - p1[i]) )
    }
    for (i in (nbHaplotype/2 + 1):nbHaplotype){
      aux[i,]=sample( c(0,1), nbLocus , replace = TRUE, prob = c(p1[i], 1 - p1[i]) )
    }
    admixed.genotype = rbind(genotype[1:(nbHaplotype/2),],genotype[1:(nbHaplotype/2),]) * (1-aux) +
      rbind(genotype[(nbHaplotype/2+1):nbHaplotype,],genotype[(nbHaplotype/2+1):nbHaplotype,]) * aux
    tess <- data.frame(coord, admixed.genotype)


  } else {
    p1 <- rep(0,nbHaplotype)
    p1[(nbHaplotype/2 + 1):nbHaplotype] <- 1
    admixed.genotype = genotype
    tess <- data.frame(coord, genotype)

  }

  # true Q matrix
  Q=matrix(p1)
  Q=cbind(Q,1-p1)
  barplot(t(Q))

  # true F matrix
  G=apply(genotype[1:(nbHaplotype/2),],2,mean)
  G=matrix(G)
  G=cbind(G,as.vector(apply(genotype[(nbHaplotype/2+1):nbHaplotype,],2,mean)))

  #PCA
  pc = prcomp(admixed.genotype, scale = FALSE) # Remark when setting T : there are some outliers but if we set F they disappear !!
  plot(pc$x, pch = 19, cex = 2, col = rep(c("blue2","orange"), each = nbHaplotype /2))
  points( ( ((t(G)-pc$center)) %*% pc$rotation)[,1],
          ( ((t(G)-pc$center)) %*% pc$rotation)[,2],
          pch = 19, cex = 2, col = c("green") )#Dont forget to center and scale if scale = T

  #############
  #Export Data#
  #############




  #create diploid data
  if( d == 2) {
    admixed.genotype = admixed.genotype[ which((1:nbHaplotype)%%2 == 0 ), ] + admixed.genotype[ which((1:nbHaplotype)%%2 == 1 ), ]
    coord = coord[ which((1:nbHaplotype)%%2 == 0 ), ]
    Q = Q[which((1:nbHaplotype)%%2 == 0 ),]

  }

  res = list()

  res$genotype = admixed.genotype
  res$coord = coord
  res$tessData =  tess
  res$Q =  Q
  res$F =  G
  res$outputName = outputName
  res$nbIndiv =  nbIndiv
  res$nbLocus =  nbLocus
  res$k =  k
  res$Nm =  Nm
  #Compute the Fst on haplotype
  Fst2Pop = function( x ) { return((var(x)-var(x[1:(length(x)/2)])/2-var(x[(length(x)/2+1):length(x)])/2)/max(.Machine$double.eps,var(x) ) ) }
  res$Fst = mean( apply(as.matrix(tess[,c(-1,-2)]),2,Fst2Pop) )
  print(paste("average Fst =",res$Fst))

  if (write) {
    print(paste("-------Exporting data", outputName))

    #export .mat format
    write.table(file = paste(outputName,".mat",sep=""),
                admixed.genotype, row.names = F, col.names = F, quote = F, sep = " ")

    #export .geno format
    write.table(file = paste(outputName,".geno",sep=""),
                t(admixed.genotype), row.names = F, col.names = F, quote = F, sep = "")


    #export .coord file
    write.table(file = paste(outputName,".coord",sep=""),
                coord, row.names = F, col.names = F, quote = F, sep = " ")

    #export .tess format
    write.table(file = paste(outputName,".tess",sep=""),
                tess, row.names = F, col.names = F, quote = F, sep = " ")


    #export Q
    write.table(file = paste(outputName,".trueQ",sep=""),
                Q, row.names = F, col.names = F, quote = F, sep = " ")

    #export G
    write.table(file = paste(outputName,".trueF",sep=""),
                G, row.names = F, col.names = F, quote = F, sep = " ")


    #export param file
    system(paste("echo ", "\"",
                 "nbIndiv = ",nbIndiv,"\n",
                 "nbLocus = ", nbLocus,"\n",
                 "Nm = ", Nm,"\n",
                 "Fst = ", res$Fst,"\n",
                 "k = ",k,'\n',"\"",
                 "> ",outputName, ".param",
                 sep=""))

  }

  save(res,file = resFile)
  # restore old wd
  setwd(wd.old)
  return(res)

}



#' Sample 2 pop simulated with ms. With outlier loci
#'
#'
#' @export
admixture2popSimulationWrapperWithOutlier <- function(nbIndiv,
                                           nbLocus,
                                           prop1,
                                           prop2,
                                           Nm1,
                                           Nm2,
                                           k,
                                           d = 1,
                                           admixe = FALSE) {

  tmpdir = tempdir()
  data.neutre = admixture2popSimulationWrapper(tmpdir,
                                               nbIndiv,
                                               prop1 * nbLocus ,
                                               d = d,
                                               Nm = Nm1,
                                               k = k,
                                               force = TRUE,
                                               admixe = admixe,
                                               write = FALSE)


  data.outlier = admixture2popSimulationWrapper(tmpdir,
                                                nbIndiv,
                                                prop2 * nbLocus ,
                                                d = d,
                                                Nm = Nm2,
                                                k = k,
                                                force = TRUE,
                                                admixe = admixe,
                                                write = FALSE)

  res = list()

  res$genotype = cbind(data.neutre$genotype,data.outlier$genotype)
  res$neutre = 1:(prop1 * nbLocus)
  res$outlier = (prop1 * nbLocus+1):(prop1 * nbLocus + prop2 * nbLocus)
  res$coord = data.neutre$coord

  return(res)
}




#' Sample genotype and environmental gradient.
#'
#' We use lfmm gibs sampler model to generate logit(P_ij) then we generate binary data G_ij
#'
#'
#' @export
lfmm_logit_sample <- function(n,
                               L,
                               K,
                               outlier_prop = 0.0,
                               gamma=1.0,
                               gamma_mu=1.0,
                               gamma_U=1.0,
                               gamma_V=1.0,
                               gamma_B=1.0,
                               mean_B_outlier = 0.5,
                               cor_U1_X = 0.5) {


  G = matrix(0.0,n,L)
  Id_k = diag(K)

  mu = rnorm(L, 0.0, gamma_mu)

  U = MASS::mvrnorm(n, mu = rep(0.0,K), Sigma = gamma_U * Id_k)

  V = MASS::mvrnorm(L, mu = rep(0.0,K), Sigma = gamma_V * Id_k)

  nb_outlier = floor(L * outlier_prop)
  B = c( rnorm( L - nb_outlier, 0.0, gamma_B ), rnorm( nb_outlier, mean_B_outlier, gamma_B ) )
  if (nb_outlier > 0 ) {
    outlier = (L - floor(L * outlier_prop) + 1):L
  } else {
    outlier = c()
  }

  # Generate X
  theta <- acos(cor_U1_X)             # corresponding angle
  aux2    <- rnorm(n, 0.0, 0.5)      # new random data
  aux     <- cbind(U[,1], aux2)         # matrix
  auxctr  <- scale(aux, center=TRUE, scale=FALSE)   # centered columns (mean 0)

  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(auxctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  aux2o  <- (Id-P) %*% auxctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  auxc2  <- cbind(auxctr[ , 1], aux2o)                # bind to matrix
  Y    <- auxc2 %*% diag(1/sqrt(colSums(auxc2^2)))  # scale columns to length 1

  X <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  cor(U[,1], X)                                    # check correlation = rho


  logitP = U %*% t(V) + X %*% t(B)
  logitP = t(apply(logitP, 1, function(row) { return(row+mu) }))
  logitP = logitP +  MASS::mvrnorm(n, mu = rep(0.0,L), Sigma = gamma * diag(L))  # we add epsilon
  P = 1 / (1 + exp(logitP))

  G = apply(P, c(1,2), function(p) { return( rbinom(1,1,p))  })

  return(list( G = G, X = X, outlier = outlier ))

}
