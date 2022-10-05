#' Tidy Compound Discoverer Gap Filling Information
#'
#' Compound Discoverer exports detailed gap filling information in a
#' nearly-useless format. This method converts it into a tidy format for use
#' with the other data exported from Compound Discoverer.
#' @param method_data_file the path to a CSV file containing the raw data
#' @return the data in a tidy tibble
#' @details Compound Discoverer abuses Excel's row grouping feature to export
#'   information on how gaps were filled after peak detection. This method
#'   cleans it up into a format R can use. You MUST export ONLY the Name,
#'   Molecular Weight, and RT columns from the Compounds table and ONLY the Fill
#'   Status and Study File ID columns from the Filled Gaps table. Then open the
#'   XLSX file in Excel and save it as a CSV file.
#' @importFrom magrittr %>%
#' @export
tidy_gap_fill_methods <- function(method_data_file) {
  fill_data <- readr::read_csv(method_data_file);

  # find the first row of each feature
  feature_rows <- which(fill_data[, 3] == "Study File ID") - 1
  # generate feature names
  feature_names <- fill_data %>%
    dplyr::slice(feature_rows) %>%
    dplyr::mutate(Name = dplyr::if_else(is.na(Name),
                                        paste0(`Molecular Weight`, "@", `RT [min]`),
                                        as.character(Name)))
  # fill them in the original data structure
  fill_data <- fill_data %>%
    dplyr::mutate(Name = rep(feature_names$Name,
                      times = diff(append(feature_rows, nrow(fill_data) + 1)))) %>%
    # remove original 'headers'
    dplyr::slice(-feature_rows, -(feature_rows + 1)) %>%
    # rename columns
    dplyr::rename(Fill_Method = `Molecular Weight`, `Study File ID` = `RT [min]`)

  return(fill_data)
}

#' Tidy Compound Discoverer mzCloud Hits
#'
#' Compound Discoverer exports detailed information on mzCloud database matches
#' into a weird nested Excel format that doesn't work well with R. This method
#' converts it into a tidy format for use with the other data exported from
#' Compound Discoverer.
#' @param mzcloud_data a tibble containing the exported mzCloud data
#' @return the data in a tidy tibble with data nesting removed
#' @details Compound Discoverer abuses Excel's row grouping feature to export
#'   information on mzCloud database hits. This method cleans it up into a
#'   format R can use.
#' @importFrom magrittr %>%
#' @export
tidy_mzCloud_hits <- function(mzcloud_data) {
  # watch out for changes to file formatting
  if(sum(mzcloud_data$Formula == "Tags", na.rm = TRUE) !=
     sum(!is.na(mzcloud_data$`Reference Ion`))) {
    warning("Number of 'Tags' entries in Formula column doesn't match number of reference ions! Check exported spreadsheet for formatting changes!")
  }

  # actual feature rows always list reference ion
  feature_idx <- which(!is.na(mzcloud_data$`Reference Ion`))

  mzcloud_features <- mzcloud_data %>%
    dplyr::slice(feature_idx) %>%
    # generate feature IDs
    dplyr::mutate(Feature_ID = paste0(`Calc. MW`, "@", `RT [min]`))

  # pull out column names for the mzCloud hit tables
  hit_cols <- c(as.character(mzcloud_data[2, ]), "Feature_ID")
  # fix NAs that got converted to strings
  hit_cols[which(hit_cols == "NA")] <- NA
  hit_na_cols <- which(is.na(hit_cols))

  df <- mzcloud_data %>%
    # merge feature IDs into the original data structure
    dplyr::mutate(Feature_ID = rep(mzcloud_features$Feature_ID,
                                   times = diff(append(feature_idx, nrow(mzcloud_data) + 1)))) %>%
    # remove original 'headers'
    dplyr::slice(-feature_idx, -(feature_idx + 1)) %>%
    # drop unused columns
    dplyr::select(-dplyr::all_of(hit_na_cols))

  # rename columns
  names(df) <- hit_cols[-hit_na_cols]

  df2 <- dplyr::left_join(mzcloud_features, df, by = "Feature_ID",
                   suffix = c(".feature", ".mzCloud")) %>%
    readr::type_convert(col_types = readr::cols())
  return(df2)
}

#' Export features for import to ID-X method editor
#'
#' The method editor software for the ID-X can import properly-formatted TSV
#' files into the Targeted Mass In/Exclusion nodes. This function accepts rows
#' from the raw Compound Discoverer output format and saves each feature in the
#' TSV format.
#'
#' @param raw_data a data.frame or equivalent formatted like Compound Discoverer
#'   export files
#' @param file Base file name to write to. Ionization polarity and file
#'   extension will be added automatically.
#' @param rt_range Half width of the retention time window around the feature's
#'  observed retention time, in minutes.
#' @param mode Only export features with the given polarity ("positive",
#'  "negative", or "both"). Default is "both".
#'
#' @return Invisibly returns a list with the formatted table(s), separated by
#'  polarity.
#'
#' @section Note:
#'
#' Versions of Compound Discoverer <3.3 don't include raw m/z values or polarity
#' in the exported data, so this method will not work for older datasets.
#'
#' @importFrom magrittr %>%
#' @export
save_target_feature_table <- function(raw_data, file, rt_range = 0.5,
                                      mode = c("both", "positive", "negative")) {
  # check mode option
  mode = match.arg(mode)

  full_table <- dplyr::select(raw_data, Name, Mass = `Calc. MW`, RT = `RT [min]`,
                              `m/z`, Ref_Ion = `Reference Ion`) %>%
    # pull charge and sign
    tidyr::extract(col = Ref_Ion, into = c("Mode", "z"), regex = "([+-])(\\d)$",
            remove = FALSE, convert = TRUE) %>%
    dplyr::mutate(Mode = factor(Mode, levels = c("+", "-"),
                         labels = c("positive", "negative"))) %>%
    # add RT ranges
    dplyr::mutate(`t start (min)` = RT - rt_range, `t stop (min)` = RT + rt_range) %>%
    # change to export format
    dplyr::select(Mode, `m/z`, z, `t start (min)`, `t stop (min)`, Name)

  # create output list
  output <- list()

  # save positive ion table if requested
  if(mode == "both" | mode == "positive") {
    file_name <- paste0(file, "_positive.tsv")
    # pull positive ions
    df <- full_table %>%
      dplyr::filter(Mode == "positive") %>%
      dplyr::select(-Mode)
    readr::write_tsv(df, file = file_name, na = "")
    output[["positive"]] <- df
  }
  if(mode == "both" | mode == "negative") {
    file_name <- paste0(file, "_negative.tsv")
    # pull positive ions
    df <- full_table %>%
      dplyr::filter(Mode == "negative") %>%
      dplyr::select(-Mode)
    readr::write_tsv(df, file = file_name, na = "")
    output[["negative"]] <- df
  }
  return(invisible(output))
}

#' Read MSP file of feature spectra
#'
#' Feature spectra can be exported from Compound Discoverer alignment results by
#' exporting all features to an mzVault library and exporting that library from
#' mzVault as an MSP file. This function reads the spectral data into a useful
#' format.
#'
#' @param path The location of the MSP file to import
#'
#' @return A list containing precursor m/z values, retention times, scan indices,
#' scan modes, raw header data, and spectra for each entry in the MSP file.
#'
#' @importFrom magrittr %>%
#' @export
read_msp <- function(path) {
  # open file
  fd <- file(path, "r", blocking = TRUE)
  # read in the file
  raw_data <- readLines(fd)
  # close the connection
  close(fd)

  # drop comment lines
  df <- raw_data[-which(stringr::str_starts(raw_data, "#"))]
  # find blank lines separating spectra
  blank_lines <- which(stringr::str_length(df) == 0)
  spec_start <- c(1, blank_lines + 1)
  spec_end <- c(blank_lines - 1, length(df))

  # parse each spectrum
  spectra <- lapply(X = seq_along(spec_start), FUN = function(i) {
    spec_lines <- df[spec_start[i]:spec_end[i]]
    # headers aren't consistent, but spectra are...
    frags <- stringr::str_match(spec_lines, "^(\\d+\\.\\d+) (\\d+\\.\\d+)$")
    frags <- frags[, -1]
    colnames(frags) <- c("mz", "intensity")
    # if no frags found, skip this entry
    n_frags <- sum(!is.na(frags[, 1]))
    if(n_frags == 0) {
      return(list(header = NA_character_, fragments = NA_real_))
    }
    # non-matching lines must be headers
    hdr_lines <- spec_lines[which(is.na(frags[, 1]))]
    hdrs <- stringr::str_split(hdr_lines, " = ", n = 2, simplify = TRUE)
    colnames(hdrs) <- c("key", "value")
    # drop non-fragment lines from fragments matrix
    frags <- frags[which(!is.na(frags[, 1])), ]
    # bizarrely, this actually does convert strings to numbers...
    class(frags) <- "numeric"

    return(list(header = hdrs, fragments = frags))
  })

  # extract header data for each spectrum
  df <- lapply(X = purrr::map(spectra, "header"), FUN = function(hdr) {
    if(is.null(dim(hdr))) {
      return(list("mz" = NA_real_, "rt" = NA_real_, "idx" = NA_real_,
                  "mode" = NA_character_))
    }
    mz <- as.numeric(hdr[, "value"][which(hdr[, "key"] == "MS:1000744|Selected Ion m/z")])
    rt <- as.numeric(hdr[, "value"][which(hdr[, "key"] == "MS:1000894|RetentionTime")])
    idx <- as.numeric(hdr[, "value"][which(hdr[, "key"] == "MS:1009001|Spectrum index")])
    mode <- hdr[, "key"][which(stringr::str_starts(hdr[, "key"], "MS:1000130"))]
    # strip off extra mode text
    mode <- stringr::str_sub(mode, start = 12, end = 19)
    return(list("mz" = mz, "rt" = rt, "idx" = idx, "mode" = mode))
  })
  # convert to tibble
  df <- dplyr::bind_rows(df)
  # note any bad headers
  bad_idx <- which(is.na(df$mz))
  df <- df[-bad_idx, ]

  # assemble results
  result <- list(mz = df$mz,
                 rt = df$rt,
                 idx = df$idx,
                 mode = df$mode,
                 headers = purrr::map(spectra, "header")[-bad_idx],
                 spectra = purrr::map(spectra, "fragments")[-bad_idx])
  return(result)
}

#' Search MSP spectra for specific features
#'
#' Feature spectra can be exported from Compound Discoverer alignment results by
#' exporting all features to an mzVault library and exporting that library from
#' mzVault as an MSP file. This function searches an imported MSP library for
#' spectra with particular precursor m/z ratios at specific retention times.
#' This allows one to quickly find all spectra associated with a given feature
#' or set of features.
#'
#' @param msp_data Spectral data generated by [read_msp()]
#' @param mz Vector of m/z values to search for. NOTE this should NOT be neutral
#' masses!
#' @param rt Vector of retention time values to search for
#' @param max_ppm_diff The largest difference allowed between query masses and
#' MSP masses, in PPM
#' @param max_rt_diff The largest difference allowed between query and MSP
#' retention times, in minutes
#'
#' @return A tibble containing all features in the MSP file that match one of
#' the query features.
#'
#' @importFrom magrittr %>%
#' @export
match_msp <- function(msp_data, mz, rt, max_ppm_diff = 25, max_rt_diff = 0.1) {
  # find close masses (within 25 ppm)
  delta_mass <- outer(mz, msp_data$mz,
                      FUN = function(x, y) { return((x - y)/y*1e6) })
  colnames(delta_mass) <- stringr::str_c(msp_data$mz, "@", msp_data$rt)
  rownames(delta_mass) <- stringr::str_c(mz, "@", rt)

  feature_comparison <- tibble::as_tibble(delta_mass, rownames = "Query_Features",
                                  .name_repair = "minimal") %>%
    tidyr::pivot_longer(cols = -Query_Features, names_to = "MSP_Features", values_to = "ppm") %>%
    dplyr::filter(abs(ppm) < max_ppm_diff) %>%
    dplyr::left_join(y = tibble::tibble(MSP_Features = stringr::str_c(msp_data$mz, "@", msp_data$rt),
                                MSP_RT = msp_data$rt, MSP_mz = msp_data$mz,
                                MSP_Index = seq.int(along.with = msp_data$mz)),
                     by = c("MSP_Features")) %>%
    dplyr::mutate(Query_RT = rt[match(Query_Features, rownames(delta_mass))],
                  Query_mz = mz[match(Query_Features, rownames(delta_mass))]) %>%
    dplyr::mutate(dRT = Query_RT - MSP_RT) %>%
    dplyr::filter(abs(dRT) < max_rt_diff)
  return(feature_comparison)
}

#' Search raw MS2 spectra for specific features
#'
#' Thermo RAW files can be converted to mzML format with ProteoWizard and loaded
#' into R with [MSnbase::readMSData()]. This function searches imported MS2 data
#' for spectra with particular precursor m/z ratios at specific retention times.
#' This allows one to quickly find all spectra associated with a given feature
#' or set of features.
#'
#' @param msn_data Spectral data generated by [MSnbase::readMSData()]
#' @param mz Vector of m/z values to search for. NOTE this should NOT be neutral
#' masses!
#' @param rt Vector of retention time values to search for
#' @param mode Vector of +/-1 indicating the polarity of the m/z values
#' @param max_ppm_diff The largest difference allowed between query masses and
#' MS2 precursor ion masses, in PPM
#' @param max_rt_diff The largest difference allowed between query and MS2 scan
#' retention times, in minutes
#'
#' @return A tibble containing data on all scans in the MS2 files that match any
#' of the query features.
#'
#' @importFrom magrittr %>%
#' @export
match_ms2_data <- function(msn_data, mz, rt, mode, max_ppm_diff = 25, max_rt_diff = 0.1) {
  # grab precursor m/z values for every MS2 scan
  prec_mz <- MSnbase::precursorMz(msn_data)
  # grab retention times for every MS2 scan, convert to minutes
  scan_rt <- MSnbase::rtime(msn_data)/60
  # grab scan polarities (1 = ESI+, 0 = ESI-), convert to signs
  scan_polarity <- ifelse(MSnbase::polarity(msn_data) == 1, 1, -1)

  # add signs to m/z values based on scan polarity
  prec_mz <- scan_polarity*prec_mz
  mz <- mode*mz

  # assemble data in a tibble
  df <- tibble::tibble(MS2_Features = stringr::str_c(prec_mz, "@", scan_rt),
                       MS2_RT = scan_rt, MS2_mz = prec_mz,
                       MS2_Index = seq.int(along.with = prec_mz),
                       # don't do this -- it loads all the data to calculate it
                       # MS2_Signal = MSnbase::ionCount(msn_data),
                       MS2_File = MSnbase::fileNames(msn_data)[MSnbase::fromFile(msn_data)])

  # find close masses (within 25 ppm)
  delta_mass <- outer(mz, prec_mz,
                      FUN = function(x, y) { return((x - y)/y*1e6) })
  colnames(delta_mass) <- stringr::str_c(prec_mz, "@", scan_rt)
  rownames(delta_mass) <- stringr::str_c(mz, "@", rt)

  feature_comparison <- tibble::as_tibble(delta_mass, rownames = "Query_Features",
                                          .name_repair = "minimal") %>%
    tidyr::pivot_longer(cols = -Query_Features, names_to = "MS2_Features", values_to = "ppm") %>%
    dplyr::filter(abs(ppm) < max_ppm_diff) %>%
    dplyr::left_join(y = df, by = c("MS2_Features")) %>%
    dplyr::mutate(Query_RT = rt[match(Query_Features, rownames(delta_mass))],
                  Query_mz = mz[match(Query_Features, rownames(delta_mass))]) %>%
    dplyr::mutate(dRT = Query_RT - MS2_RT) %>%
    dplyr::filter(abs(dRT) < max_rt_diff)

  return(feature_comparison)
}


#' Export a Spectrum2 object in FreeStyle CSV format
#'
#' Uses the data in a [MSnbase::Spectrum2] object to write a CSV file in the
#' same format used by FreeStyle to export spectra. These files can be used to
#' search the mzCloud online library.
#'
#' @param spectrum A Spectrum2 object to be exported
#' @param file File to write the spectrum to
#'
#' @return The input spectrum, invisibly
#'
#' @importFrom magrittr %>%
#' @export
export_idx_spectrum <- function(spectrum, file) {
  # assemble stuff to write before writing it, so we don't leave files open

  # construct scan description string
  if(polarity(spectrum) == 1) {
    scan_polarity <- "+"
  } else {
    scan_polarity <- "-"
  }
  # pull m/z values
  mz <- MSnbase::mz(spectrum)
  # pull fragment intensities
  intensities <- MSnbase::intensity(spectrum)
  # calculate range of scan
  mz_range <- range(c(mz, MSnbase::precursorMz(spectrum)))
  mz_range[1] <- floor(mz_range[1])
  mz_range[2] <- ceiling(mz_range[2])
  scan_info <- paste0("FTMS ", scan_polarity, " p ESI d Full ms",
                      MSnbase::msLevel(spectrum), " ",
                      MSnbase::precursorMz(spectrum), "@hcd",
                      MSnbase::collisionEnergy(spectrum),
                      " [", mz_range[1], "-", mz_range[2], "]")
  # grab original scan index
  scan_num <- MSnbase::scanIndex(spectrum)
  # get retention time in minutes
  scan_rt <- MSnbase::rtime(spectrum)/60
  # number of peaks in spectrum
  n_peaks <- MSnbase::peaksCount(spectrum)

  # open file
  fd <- file(file, "w", blocking = TRUE)

  # write fixed header and a fake file path
  writeLines(c("SPECTRUM - MS", "j_random_file.raw"), fd)
  # write scan data
  writeLines(c(scan_info,
               paste0("Scan #: ", scan_num),
               paste0("RT: ", scan_rt),
               "AV: 1",
               paste0("Data points: ", n_peaks)),
             fd)
  # write the spectrum itself
  writeLines("Mass,Intensity", fd)
  writeLines(stringr::str_c(mz, intensities, sep = ","), fd)

  # close the connection
  close(fd)

  return(invisible(spectrum))
}
