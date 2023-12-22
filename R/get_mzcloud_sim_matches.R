#' Retrieve mzCloud library similarity matches
#'
#' Queries a Compound Discoverer mass spec alignment database to get the
#' mzCloud spectral library similarity search results for some or all compounds.
#'
#' @param msa An ms_alignment object to query
#' @param ids Compound IDs to get matches for, or NULL to get all matches
#'
#' @return A tibble with the mzCloud similarity results, or NA if this alignment
#'  database does not contain such data. An empty tibble indicates that similarity
#'  search results exist, but none of the requested compounds have such matches.
#'
#' @importFrom rlang .data
#' @export
get_mzcloud_sim_matches <- function(msa, ids = NULL) {
  # check to see if there are any hits at all
  if(!DBI::dbExistsTable(msa$db_connection, "MzCloudSimilaritySearchResultItems")) {
    simpleWarning("This alignment database does not have mzCloud similarity matches!")
    return(NA)
  }
  # if no IDs are specified, return hits for all consolidated compounds
  if(is.null(ids)) {
    ids <- msa$unknown_compound_items$ID
  }

  # get mapping from compounds to hits
  compound_hit_key <- dplyr::tbl(msa$db_connection, "ConsolidatedUnknownCompoundItemsMzCloudSimilaritySearchResultItems")
  # read the mzCloud hit info for all items
  mzcloud_hit_tbl <- dplyr::tbl(msa$db_connection, "MzCloudSimilaritySearchResultItems")

  # get hits and search results associated with compounds of interest
  cloud_matches <- dplyr::filter(compound_hit_key, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  cloud_matches <- dplyr::left_join(x = cloud_matches, y = mzcloud_hit_tbl,
                                    by = c("MzCloudSimilaritySearchResultItemsID" = "ID"))

  # get local copy of the data
  cloud_matches <- dplyr::collect(cloud_matches)

  return(cloud_matches)
}
