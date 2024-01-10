#' Retrieve factor labels in Compound Discoverer alignment results
#'
#' CompoundDiscoverer alignment databases contain many enum/factor columns
#' encoded as integer values. This function looks up the labels for a given
#' factor in a particular table.
#'
#' @note
#' Compound Discoverer software is produced by Thermo Fisher Scientific. This
#' package is not affiliated with Thermo Fisher Scientific in any way. For an
#' official Python interface to Compound Discoverer alignment files, see
#' \url{https://github.com/thermofisherlsms/pyeds}
#'
#' @param msa An ms_alignment object to query
#' @param table Name of data table where enum is located
#' @param enum Name of the enum/factor column
#'
#' @return A tibble containing information on each value of the enum. If
#'   the column is not an enum, returns NULL with a warning.
#'
#' @importFrom rlang .data
#' @export
get_cdresult_enum <- function(msa, table, enum) {
  # get table listing all tables
  data_types <- dplyr::tbl(src = msa$db_connection, "DataTypes")
  # get table listing all table columns
  data_types_columns <- dplyr::tbl(src = msa$db_connection, "DataTypesColumns")
  # get table listing all enums
  enum_data_type_values <- dplyr::tbl(src = msa$db_connection, "EnumDataTypeValues")

  # get column info
  column_info <- dplyr::filter(data_types, .data$TableName == table)
  column_info <- dplyr::left_join(x = dplyr::select(column_info, "TableName", "DataTypeID"),
                                  y = data_types_columns, by = "DataTypeID")
  column_info <- dplyr::filter(column_info, .data$DBColumnName == enum)
  column_info <- dplyr::collect(column_info)

  # if this column isn't actually an enum
  if(column_info$SpecialValueType != "Enum") {
    warning(paste0("Column `", enum, "` in table `", table, "` is not an enum!"))
    return(NULL)
  }

  # get enum info
  enum_id <- column_info$SpecialValueTypeID[1]
  enum_levels <- dplyr::filter(enum_data_type_values,
                               .data$EnumID == enum_id)

  return(dplyr::collect(enum_levels))
}
