#' load code on patator
#'
#'
#' @export
runPatator <- function() {

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
  save.session()

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

#' build the package on patator
#'
#'
#' @export
buildOnPatator <- function() {

  cat("TODO")


}
