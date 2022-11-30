# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("GeoTcgaData", pkgname, libname)
# }

# GeoTcgaDataStartupMessage <- function() {

#   # msg <- "Hello, friend! welcome to use!"
#   msg <- paste0(
#     " =============================================================\n",
#     " Hello, friend! welcome to use GeoTcgaData!                   \n",
#     " -------------------------------------------------------------\n",
#     " Version:", utils::packageVersion("GeoTcgaData"), "\n",
#     " =============================================================\n"
#   )
#   return(msg)
# }

# .onAttach <- function(lib, pkg) {
#   msg <- GeoTcgaDataStartupMessage()
#   packageStartupMessage(msg)
#   invisible()
# }

# This function comes from snpGeneSets package
# .onLoad <- function(libname, pkgname) {
#     extdata <- system.file("extdata", package=pkgname,
#         lib.loc=libname, mustWork=TRUE)
#     if(!file.exists(file.path(extdata, "gene105.db"))) {
#         cat("downloading and installing Gene database of GRCh37/hg19...\n")
#         downloader::download(
#             url="https://www.umc.edu/apps/files/GeneticStudy/gene105.zip",
#                 destfile=file.path(extdata, "gene105.zip"))
#         unzip(file.path(extdata, "gene105.zip"), exdir = extdata)
#         if(file.exists(file.path(extdata, "gene105.zip")))
#             file.remove(file.path(extdata, "gene105.zip"))
#     }
#     if(!file.exists(file.path(extdata, "gene106.db"))) {
#         cat("downloading and installing Gene database of GRCh38/hg38...\n")
#         downloader::download(
#             url="https://www.umc.edu/apps/files/GeneticStudy/gene106.zip",
#                 destfile=file.path(extdata, "gene106.zip"))
#         unzip(file.path(extdata, "gene106.zip"), exdir = extdata)
#         if(file.exists(file.path(extdata, "gene106.zip")))
#             file.remove(file.path(extdata, "gene106.zip"))
#     }
# }
