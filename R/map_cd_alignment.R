#' Map data table relationships in Compound Discoverer alignment results
#'
#' Gets the relationships between alignment data tables.
#'
#' @details
#' Compound Discoverer alignment results store data across many tables in an
#' SQLite database. The database also contains a map of how these tables relate
#' to each other. This function generates an igraph network from the map.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param tables Also return the data relationship tables?
#'
#' @return An igraph object containing the relationship between data tables. If
#' tables is TRUE, returns a list containing the graph, a table describing all
#' data tables, and a table describing relationships between tables.
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
map_cd_alignment <- function(msa, tables = FALSE) {
  # grab the tables
  data_types <- dplyr::tbl(msa$db_connection, "DataTypes") %>%
    dplyr::collect()
  data_connections <- dplyr::tbl(msa$db_connection, "ConnectedDataTypes") %>%
    dplyr::collect()

  # construct table of relationships
  # rather than include the 'ConnectedData' tables as nodes, we set them as
  # labels on the edge between the two main tables
  df <- data_connections %>%
    dplyr::left_join(y = dplyr::select(data_types, "DataTypeID", "TableName"),
                     by = c("DataTypeID1" = "DataTypeID")) %>%
    dplyr::left_join(y = dplyr::select(data_types, "DataTypeID", "TableName"),
                     by = c("DataTypeID2" = "DataTypeID"),
                     suffix = c(".1", ".2")) %>%
    dplyr::select("TableName.1", "TableName.2", "ConnectedTableName")

  # add unconnected tables
  connected_tables <- c(df$TableName.1, df$TableName.2)
  df2 <- data_types %>%
    dplyr::select("TableName") %>%
    dplyr::filter(!(.data$TableName %in% connected_tables)) %>%
    # make self-loops
    dplyr::mutate(TableName.1 = .data$TableName) %>%
    dplyr::rename(TableName.2 = "TableName")
  df <- dplyr::bind_rows(df, df2)

  # construct the graph
  data_graph <- igraph::graph_from_data_frame(df, directed = FALSE)

  if(tables) {
    return(list("data_graph" = data_graph,
                "data_types" = data_types,
                "connected_data_types" = data_connections))
  } else {
    return(data_graph)
  }
}
