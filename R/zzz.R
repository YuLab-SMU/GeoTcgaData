# .onLoad <- function(libname, pkgname) 
# {
#   library.dynam("GeoTcgaData", pkgname, libname)
# }

GeoTcgaDataStartupMessage <- function()
{

  # msg <- "Hello, friend! welcome to use!"
    msg <- paste0(
        " =============================================================\n",
        " Hello, friend! welcome to use GeoTcgaData!                   \n",
        " -------------------------------------------------------------\n",
        " Version:",utils::packageVersion("GeoTcgaData"),"\n",
        " =============================================================\n"
    )
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  msg <- GeoTcgaDataStartupMessage()
  packageStartupMessage(msg)      
  invisible()
}
