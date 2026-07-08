#' Import the data in a mass spec data alignment
#'
#' Certain mass spec data alignment results are stored as SQLite databases. This
#' function imports relevant result tables for further analysis.
#'
#' @param f Path to file containing an alignment project in SQLite format.
#'
#' @return A list containing the alignment data tables.
#'
#' @export
import_ms_alignment <- function(f) {
  # TODO: rework to only import table data and main unknown compound table?
  # then write other methods to take alignment list/S3 object and desired IDs for extraction
  # probably also pull in WorkflowInputFiles so we have the file IDs
  # TODO: set up query functions to optionally pull subsets of data instead of whole tables
  # use DBI::dbGetQuery for this purpose

  # make a temp directory for unzipping blobs
  # tmp_zip_dir <- tempfile(pattern = "MSA_dir", tmpdir = tempdir())
  # dir.create(tmp_zip_dir)
  con <- NULL

  msa <- tryCatch({
    # connect to the alignment database
    con <- DBI::dbConnect(RSQLite::SQLite(), f)

    # read input files table
    input_files <- DBI::dbReadTable(con, "WorkflowInputFiles")

    # read unknown compound items
    consolidated_unk_comp_items <- DBI::dbReadTable(con, "ConsolidatedUnknownCompoundItems")
    # this table is the main compounds table seen in CD
    # some of the info flags are encoded as binary blobs

    # build S3 object with class ms_alignment
    msa <- list(db_file = f,
                db_connection = con,
                input_files = input_files,
                unknown_compound_items = consolidated_unk_comp_items)
    class(msa) <- "ms_alignment"
    # return MSA
    return(msa)
  }, error = function(e) {
    message("Failed to load CompoundDiscoverer database!")
    # clean up the connection
    if(!is.null(con)) {
      DBI::dbDisconnect(con)
    }
    # re-throw the error
    stop(e)
  })
  # connect to the alignment database
  con <- DBI::dbConnect(RSQLite::SQLite(), f)

  # read input files table
  input_files <- DBI::dbReadTable(con, "WorkflowInputFiles")

  # read unknown compound items
  consolidated_unk_comp_items <- DBI::dbReadTable(con, "ConsolidatedUnknownCompoundItems")
  # this table is the main compounds table seen in CD
  # some of the info flags are encoded as binary blobs

  # build S3 object with class ms_alignment
  msa <- list(db_file = f,
              db_connection = con,
              input_files = input_files,
              unknown_compound_items = consolidated_unk_comp_items)
  class(msa) <- "ms_alignment"

  return(msa)
}

# internal function that parses RT correction data
extract_rt_corrections <- function(rt_corr_tbl) {
  # retention times are 32-bit integers
  # not sure what time format this is, if any
  # first integer is just the length of the array...
  retention_times <- lapply(rt_corr_tbl$RetentionTimes, FUN = function(b) {
    return(readBin(b, "integer", n = length(b)/8, size = 8)[-1])
  })

  corr_values <- lapply(rt_corr_tbl$CorrectionValues, FUN = function(b) {
    return(readBin(b, "integer", n = length(b)/8, size = 8)[-1])
  })

  pred_tols <- lapply(rt_corr_tbl$PredictionTolerances, FUN = function(b) {
    return(readBin(b, "integer", n = length(b)/8, size = 8)[-1])
  })

  df <- lapply(X = 1:length(retention_times), FUN = function(i) {
    return(tibble::tibble(WorkflowId = rt_corr_tbl$WorkflowId[[i]],
                          ID = rt_corr_tbl$ID[[i]],
                          RetentionTimes = retention_times[[i]],
                          CorrectionValues = corr_values[[i]],
                          PredictionTolerances = pred_tols[[i]]))
  })
  return(dplyr::bind_rows(df))
}

#' Fix a 32 bit unsigned integer that has been read as signed
#'
#' This is really just to fix a limitation of readBin/R's 32 bit signed ints
#' @param x Number to be fixed
#' @param adjustment number to be added to convert to uint32 (2^32 by default)
#' @return numeric value of uint32
#' @author jefferis
#' @seealso \code{\link{readBin}}
ConvertIntToUInt<-function(x,adjustment=2^32){
  x=as.numeric(x)
  signs=sign(x)
  x[signs<0]=x[signs<0]+adjustment
  x
}
