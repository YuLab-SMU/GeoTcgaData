# .onLoad <- function(libname, pkgname) 
# {
#   library.dynam("GeoTcgaData", pkgname, libname)
# }

GeoTcgaDataStartupMessage <- function()
{

 msg <- "Hello, friend! welcome to use!"
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  msg <- GeoTcgaDataStartupMessage()
  packageStartupMessage(msg)      
  invisible()
}
