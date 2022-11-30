
## quiets concerns of R CMD check re: the .'s that appear in pipelines
# if(getRversion() >= "2.15.1") utils::globalVariables("hgnc_file")

#' Gene id conversion types
#'
#' @return a vector
#' @export
#'
#' @examples
#' id_ava()
id_ava <- function() {
    return(colnames(hgnc_file))
}


#' Gene id conversion
#'
#' @param from one of `id_ava()`
#' @param to one of `id_ava()`
#' @param IDs the gene id which needed to convert
#' @param na.rm Whether to remove lines containing NA
#'
#' @return a vector of genes
#' @export
#'
#' @examples
#' id_conversion_vector(
#'     "symbol", "ensembl_gene_id",
#'     c("A2ML1", "A2ML1-AS1", "A4GALT", "A12M1", "AAAS")
#' )
id_conversion_vector <- function(from, to, IDs, na.rm = FALSE) {
    from <- match.arg(from, id_ava())
    to <- match.arg(to, id_ava())
    loc <- which(hgnc_file[, from] %in% IDs)
    if (length(loc) > 0) {
        results_mat <- data.frame(from = hgnc_file[loc, from], 
            to = hgnc_file[loc, to])
        msg <- paste(format(100 * sum(!is.na(results_mat$to)) /
            length(IDs), digits = 4), "%",
        " were successfully converted.",
        sep = ""
        )
        message(msg)
        if (na.rm == TRUE) {
            results_mat <- results_mat[!is.na(results_mat$to), ]
        }
    } else {
        stop("Failed to match")
    }
    return(results_mat)
}
