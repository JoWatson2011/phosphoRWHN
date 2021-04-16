#' Get EnrichR results from API
#'
#'
#' @param genes character vector of HGNC symbols
#' @param enrichrLib EnrichR library to use; see https://maayanlab.cloud/Enrichr/#stats
#' @param pval enrichment pvalue cut off
#'
#' @importFrom httr GET POST http_error
#' @importFrom dplyr filter
#'
#' @return
#'
getEnrichr <- function(genes, enrichrLib, pval) {

  temp <- httr::POST(url = "https://maayanlab.cloud/Enrichr/enrich",
                     body =
                       list(list = paste(genes,
                                         collapse = "\n")))

  if (httr::http_error(temp)) {
    stop(
      sprintf(
        "EnrichR POST API request failed [%s]",
        status_code(STRING_api_post)
      ),
      call. = F
    )
  }

  httr::GET(url = "https://maayanlab.cloud/Enrichr/share")


  r <- httr::GET(
    url = "https://maayanlab.cloud/Enrichr/export",
    query = list(file = "API",
                 backgroundType = enrichrLib)
  )

  if (httr::http_error(r)) {
    stop(
      sprintf(
        "EnrichR POST API request failed [%s]",
        status_code(STRING_api_post)
      ),
      call. = F
    )
  }
  r <- gsub("&#39;", "'", intToUtf8(r$content))
  tc <- textConnection(r)
  enriched <-
    read.table(
      tc,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = ""
    )

  if(nrow(enriched) > 0){
    enriched <- dplyr::filter(enriched, .data$Adjusted.P.value < pval)
  }

  return(enriched)
}
