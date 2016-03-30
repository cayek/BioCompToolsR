#' load code on patator
#'
#'
#' @export
runPatator <- function() {

  stop("do not work anymore !! need .Rprofile to be setup...")

  # set wd
  wd.old = getwd()
  setwd( "~/PatatorHomeDir" )


  tmpfile = tempfile(tmpdir = ".")
  if(file.exists(".Rprofile")) {
    file.rename(".Rprofile",tmpfile)
  }

  rprofile <- system.file("bash","Rprofile",package = "BioCompToolsR")
  file.copy(rprofile,".Rprofile")

  # save session
  session::save.session()

  # retrieve wd
  working_dir = unlist(strsplit(getwd(),"/PatatorHomeDir/"))[2]

  # retrieve script file
  scipt <- system.file("bash","runPatator.sh",package = "BioCompToolsR")
  system(paste("xfce4-terminal --hold -e" ,"\"",scipt,working_dir,"\"" ) )
  Sys.sleep(2)

# restore Rprofile
  file.remove(".Rprofile")
  if(file.exists(tmpfile)) {
    file.rename(tmpfile,".Rprofile")
    file.remove(tmpfile)
  }

  # restore old wd
  setwd(wd.old)


}

#' open a R session on patator
#'
#'
#' @export
openPatator <- function() {

  # retrieve wd
  working_dir = unlist(strsplit(getwd(),"/PatatorHomeDir/"))[2]

  # retrieve script file
  scipt <- system.file("bash","runPatator.sh",package = "BioCompToolsR")
  system(paste("xfce4-terminal --hold -e" ,"\"",scipt,working_dir,"\"" ) )
  Sys.sleep(2)


}


#' build the package on patator
#'
#'
#' @export
buildOnPatator <- function(pkgname = NULL) {

  if(is.null(pkgname)) {
    pkgname = c("~/PatatorHomeDir/Projects/AssociationAnalysis/associationr/",
                "~/PatatorHomeDir/Projects/BenchmarkingR/",
                "~/PatatorHomeDir/Projects/BioCompToolsR/")
  }

  for(pkg in pkgname) {
    cmd = paste("R CMD INSTALL --preclean", pkg)
    cat("-> Run on patator :", cmd,"\n")
    res = ssh.utils::run.remote(cmd, remote = "cayek@patator")
    if(res$cmd.error) {
      stop("ERROR")
    }
  }

}

#' Render Rmarkdown on patator
#'
#'
#' @export
renderOnPatator <- function(file.name, dir = "~/Notebook") {

  # retrieve script file
  scipt <- system.file("bash","renderPatator.sh",package = "BioCompToolsR")
  system(paste("xfce4-terminal -T","\"rendering",basename(file.name),"\"","--hold -e" ,"\"",scipt,dir,file.name,"\"" ) )
  Sys.sleep(2)

}
