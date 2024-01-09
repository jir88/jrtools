#' Retrieve related tables in Compound Discoverer alignment results
#'
#' Queries a Compound Discoverer mass spec alignment database to get a series of
#' related data tables for some or all compounds. If data associated with
#' specific compound IDs are requested, the database traversal path will start
#' at the ConsolidatedUnknownCompoundItems table regardless of the specified
#' starting table, but the additional tables will not be returned.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param from Data table name to start from
#' @param to Target data table name
#' @param via Character vector of data table name(s) that should be included in
#'   database traversal path, in the order they should be traversed
#' @param ids Compound IDs to get matches for, or NULL to get all matches
#'
#' @return A list containing a tibble with ID values relating items within all
#'   returned tables, plus a list containing each table in the traversal path.
#'
#' @importFrom rlang .data
#' @export
get_cdresult_tables <- function(msa, from, to, via = NULL, ids = NULL) {
  # generate database map
  db_map <- jrtools::map_cd_alignment(msa = msa, tables = TRUE)
  # split full path into jumps between desired intermediate tables
  path_segments <- c(from, via, to)
  fixed_from <- path_segments[-length(path_segments)]
  fixed_to <- path_segments[-1]
  # get traversal path for each segment
  db_paths <- mapply(FUN = igraph::shortest_paths,
                     from = fixed_from, to = fixed_to,
                     MoreArgs = list(graph = db_map$data_graph,
                                     output = "both"),
                     SIMPLIFY = FALSE)

  # get tables for all segments
  path_table_names <- purrr::map(db_paths, list("vpath", 1))
  path_table_names <- lapply(X = path_table_names, FUN = igraph::as_ids)
  # combine segments and remove duplicates
  path_table_names <- unique(c(path_table_names, recursive = TRUE))

  # get keys in each segment
  path_key_names <- purrr::map(db_paths, list("epath", 1))
  path_key_names <- lapply(X = path_key_names, FUN = igraph::edge_attr,
                           name = "ConnectedTableName", graph = db_map$data_graph)
  # combine segments and remove duplicates
  path_key_names <- unique(c(path_key_names, recursive = TRUE))

  # add CID table link if needed
  if(!is.null(ids) & path_table_names[1] != "ConsolidatedUnknownCompoundItems") {
    # only need keys for this part
    cid_path <- igraph::shortest_paths(graph = db_map$data_graph,
                                       from = "ConsolidatedUnknownCompoundItems",
                                       to = from, output = "epath")
    cid_key_names <- igraph::edge_attr(graph = db_map$data_graph,
                                       name = "ConnectedTableName",
                                       index = cid_path$epath[[1]])
    # add to existing set of key names
    path_key_names <- c(cid_key_names, path_key_names)
  }

  # collect tables
  path_tbls <- lapply(X = path_table_names, FUN = dplyr::tbl, src = msa$db_connection)

  # extract all ID keys
  key_tbls <- lapply(X = path_key_names, FUN = dplyr::tbl, src = msa$db_connection)
  # include all items in initial table, even if they don't have matches to
  # subsequent tables in the traversal path
  # get table column(s) with ID values
  id_cols <- colnames(key_tbls[[1]])
  id_cols <- id_cols[stringr::str_starts(id_cols, path_table_names[1])]
  id_col_names <- stringr::str_sub(id_cols, start = (stringr::str_length(path_table_names[1]) + 1))
  id_tbl <- dplyr::select(path_tbls[[1]], dplyr::all_of(id_col_names))

  # align ID values into a tibble
  suppressMessages({
    key_tbls <- purrr::reduce(.x = key_tbls, .f = dplyr::left_join)
  })

  # add all items from initial table
  df <- id_cols
  names(df) <- id_col_names
  key_tbls <- dplyr::left_join(x = id_tbl, y = key_tbls, by = df)
  # fix ID column name
  df <- id_col_names
  names(df) <- id_cols
  key_tbls <- dplyr::rename(key_tbls, dplyr::all_of(df))

  # filter tables if CIDs specified
  if(!is.null(ids)) {
    key_tbls <- dplyr::filter(key_tbls, .data$ConsolidatedUnknownCompoundItemsID %in% ids)
  }

  # filter tables to match keys
  df <- mapply(tbl_name = path_table_names, path_tbl = path_tbls, SIMPLIFY = FALSE, FUN = function(tbl_name, path_tbl) {
    ktcn <- colnames(key_tbls)
    id_cols <- ktcn[stringr::str_starts(ktcn, tbl_name)]
    # pull intra-table column names
    id_col_names <- stringr::str_sub(id_cols, start = (stringr::str_length(tbl_name) + 1))
    # extract columns
    id_tbl <- dplyr::select(key_tbls, dplyr::all_of(id_cols))
    # only keep one copy of each ID
    id_tbl <- dplyr::distinct(id_tbl)
    # generate join-by vector
    df <- id_col_names
    names(df) <- id_cols
    # filter the data table
    path_tbl <- dplyr::left_join(x = id_tbl, y = path_tbl, by = df)
    # collect and return data table
    return(dplyr::collect(path_tbl))
  })

  # grab local copy of keys
  key_tbls <- dplyr::collect(key_tbls)

  # construct and return results list
  return(list("keys" = key_tbls, "tables" = df))
}
