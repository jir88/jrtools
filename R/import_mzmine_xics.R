#' Import extracted ion chromatograms exported from MZmine
#'
#' This function loads extracted ion chromatogram (XIC) data exported from MZmine
#' and converts it into a tidy format that works well with the tidyverse.
#'
#' @param file Name of file with XIC data
#'
#' @return A tibble of all the XIC data
#'
#' @importFrom magrittr %>%
#' @export
import_mzmine_xics <- function(file) {
  # XIC headers appear in first row, every other column *facepalm*
  # we initially read in everything as text and pull the stupid headers manually
  raw_xics <- readxl::read_xlsx(file, col_names = FALSE, col_types = c("text"),
                          .name_repair = "unique")
  # pull the headers
  hdrs <- raw_xics %>%
    dplyr::slice(1) %>%
    tidyr::pivot_longer(cols = tidyr::everything())

  # hdrs <- readxl::read_xlsx(file, n_max = 1, col_names = FALSE) %>%
  #   pivot_longer(cols = everything())

  # grab actual header indices
  idx <- which(!is.na(hdrs$value))
  # how many times should actual headers be repeated to properly fill all indices?
  runs <- diff(c(idx, nrow(hdrs) + 1))
  # fill in missing values
  hdrs <- dplyr::mutate(hdrs, value = rep(hdrs$value[idx], times = runs))
  # extract data from headers
  # peak XICs have m/z data embedded in header, but full XICs don't
  hdrs <- tidyr::extract(hdrs, col = value, into = c(NA, "mz", "File"),
                         convert = TRUE, regex = "(m\\/z ([0-9.]+).+: )?(.+)")

  # grab the data column headers and make unique
  data_hdrs <- stringr::str_c(as.character(raw_xics[2, ]), colnames(raw_xics))
  # slice out the actual data and convert to proper types
  raw_xics <- raw_xics %>%
    dplyr::slice(-1:-2) %>%
    readr::type_convert()
  # add column names
  colnames(raw_xics) <- data_hdrs
  # raw_xics <- readxl::read_xlsx("figures/230118 XICs 342mz.xlsx", skip = 1,
  #                               .name_repair = "unique")

  tidy_xics <- raw_xics %>%
    # select(-starts_with("z-axis")) %>%
    dplyr::mutate(Index = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = -Index) %>%
    tidyr::extract(col = "name", into = c("name", "hdr_idx"),
                   regex = "(.+)(\\.{3}\\d+)", convert = TRUE) %>%
    dplyr::left_join(y = hdrs, by = c("hdr_idx" = "name")) %>%
    # left_join(y = sample_key, by = "File") %>%
    dplyr::select(-hdr_idx) %>%
    tidyr::drop_na(value) %>%
    tidyr::pivot_wider(names_from = "name", values_from = "value") %>%
    # z-axis is actually the current m/z, so move that over
    # also anything without a z-axis is a called peak
    dplyr::mutate(mz = dplyr::if_else(is.na(mz), `z-axis`, as.numeric(mz)),
           Called_Peak = is.na(`z-axis`)) %>%
    dplyr::select(-`z-axis`) %>%
    dplyr::rename(RT = `Retention time`, Intensity = `Base peak intensity`)

  return(tidy_xics)
}
