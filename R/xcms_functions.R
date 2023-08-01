#' Plot 2D intensity image for MS data
#'
#' Given an MSnExp object, generate a 2D plot with retention time on the x-axis
#' and m/z values on the y-axis. Color indicates signal intensity.
#'
#' @param ms_data a MSnExp object containing one or more samples
#' @param mz_binwidth the m/z bin width (in Da) to use when generating the plot
#' @return list of ggplot objects, one per sample
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @export
plot_spectrum_2d <- function(ms_data, mz_binwidth = 1) {
  n_files <- length(MSnbase::fileNames(ms_data))

  # weirdly, extracting each file individually is much faster
  int_range <- c()
  mz_range <- c()
  int_data <- list()
  for(i in 1:n_files) {
    # pull individual file
    ms_file <- MSnbase::filterFile(ms_data, file = i) %>%
      # drop zero intensities and empty spectra, otherwise bin freaks out
      MSnbase::clean() %>%
      MSnbase::filterEmptySpectra() %>%
      # bin m/z values
      MSnbase::bin(binSize = mz_binwidth)

    # we've binned the mz values, so we just need mz values from one scan
    binned_mz <- MSnbase::mz(ms_file)[[1]]
    # grab scan retention times and calculate scan lengths so they display right
    all_rt <- MSnbase::rtime(ms_file) %>%
      tibble::enframe(name = "Scan", value = "rt") %>%
      dplyr::mutate(rt_width = c(diff(.data$rt), NA))
    # fix last rt bin width
    all_rt$rt_width[length(all_rt$rt_width)] <- all_rt$rt_width[length(all_rt$rt_width) - 1]

    # extract intensities and add m/z and retention times
    all_intensity <- MSnbase::intensity(ms_file) %>% tibble::as_tibble() %>%
      dplyr::mutate(mz = binned_mz) %>%
      tidyr::pivot_longer(cols = -"mz", names_to = "Scan", values_to = "intensity") %>%
      # drop zero intensities
      dplyr::filter(.data$intensity > 0) %>%
      # add retention times
      dplyr::left_join(all_rt, by = "Scan")
    int_data[[i]] <- all_intensity

    # update intensity range
    int_range <- range(int_range, range(all_intensity$intensity))
    # update m/z range
    mz_range <- range(mz_range, range(binned_mz))
  }

  # now generate the actual plots
  plot_list <- list()
  for(i in 1:n_files) {
    # add plot to list
    plot_list[[i]] <- ggplot2::ggplot(int_data[[i]],
                                      ggplot2::aes(x = .data$rt, y = .data$mz,
                                                   width = .data$rt_width,
                                                   fill = log10(.data$intensity + 1))) +
      ggplot2::geom_tile() +
      ggplot2::scale_y_continuous(limits = mz_range) +
      ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, limits = log10(int_range + 1))
  }

return(plot_list)
}
