#' Combine clinical information obtained from TCGA and extract survival data
#'
#' @param Files_dir a dir data
#'
#' @return a matrix, survival time and survival state in TCGA
#' @export
#'
#' @examples
#' tcga_cli_deal(system.file(file.path("extdata", "tcga_cli"), 
#'   package = "GeoTcgaData"))
tcga_cli_deal <- function(Files_dir = "your_clinical_directory") {
  if (Files_dir != "your_clinical_directory") {
    out <- file.path(tempdir(), "clin2.txt")
    Files <- dir(Files_dir)
    file_id <- Files
    # Modify an error
    for (i in seq_len(length(file_id))) {
      file_id[i] <- unlist(strsplit(file_id[i], "\\."))[3]
    }
    time <- "</clin_shared:days_to_last_followup>"
    state <- "</clin_shared:vital_status>"

    days <- "</clin_shared:days_to_death>"
    for (file in Files) {
      # haha means the sample has no survival data
      timee <- "haha"
      statee <- "haha"
      aa <- as.matrix(utils::read.table(file.path(Files_dir, file), 
        sep = "\t", header = TRUE))
      for (i in 1:dim(aa)[1]) {
        if (length(grep(time, aa[i, 1])) > 0) {
          timee <- unlist(strsplit(aa[i, 1], ">"))[2]

          timee <- gsub("</clin_shared:days_to_last_followup", "", timee)
        }
      }

      for (i in 1:dim(aa)[1]) {
        if (length(grep(state, aa[i, 1])) > 0) {
          statee <- unlist(strsplit(aa[i, 1], ">"))[2]

          statee <- gsub("</clin_shared:vital_status", "", statee)
        }
      }
      for (i in 1:dim(aa)[1]) {
        if (length(grep(days, aa[i, 1])) > 0) {
          timee <- unlist(strsplit(aa[i, 1], ">"))[2]

          timee <- gsub("</clin_shared:days_to_death", "", timee)
        }
      }
      file_idd <- unlist(strsplit(file, "\\."))[3]

      utils::write.table(cbind(file_idd, timee, statee), out,
        sep = "\t", row.names = FALSE, col.names = FALSE, 
        append = TRUE, quote = FALSE
      )
    }
    output_cli <- data.table::fread(out, sep = "\t", header = FALSE)
    file.remove(out)
    return(output_cli)
  } else {
    message("please give your directory!")
  }
}
