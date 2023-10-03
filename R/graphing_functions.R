#' Construct fancy-looking ROC plot
#'
#' Generates a ggplot object to display a ROC curve.
#'
#' @param data output from the [calc_roc] function
#' @return a ggplot object containing the plot, so you can fine-tune it as
#'   desired before displaying it
#' @details There isn't really any error-checking here. Put garbage in, you'll
#'   get garbage out.
#'
#' @importFrom rlang .data
#' @export
fancy_roc_plot <- function(data) {
  gp <- ggplot2::ggplot(data = data$curves,
                        mapping = ggplot2::aes(x = (100 - .data$specificity),
                                               y = .data$sensitivity)) +
    ggplot2::geom_ribbon(mapping = ggplot2::aes(x = (100 - .data$specificity),
                                                ymin = .data$`cisens2.5`,
                                                ymax = .data$`cisens97.5`),
                         fill = "grey80") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_path(color = "black", size = 1.0) +
    ggplot2::labs(x = "false-positive rate (%)", y = "true-positive rate (%)") +
    ggplot2::annotate(geom = "text", x = 70, y = 40,
                      label = sprintf(fmt = "%1.1f%%\n(%1.1f%% - %1.1f%%)",
                                      data$auc, data$ci[1], data$ci[3]),
                      hjust = "center", size = 6)

  # if ggthemes is available, add a rangeframe
  if(requireNamespace("ggthemes", quietly = TRUE)) {
    gp <- gp +
      ggthemes::geom_rangeframe(color = "black")
  }

  return(gp)
}

#' Generate plotmath-formatted scientific notation expressions
#'
#' Converts a number into scientific notation using plotmath formatting. The
#' resulting expression can be used in ggplot2 labels.
#'
#' @param l a single numeric value to be formatted
#' @return an expression in plotmath format for displaying this value
#' @details There isn't really any error-checking here. Put garbage in, you'll
#'   get garbage out.
#' @export
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # make zero actually zero
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # return this as an expression
  parse(text=l)
}

#' Pretty scientific notation labeling function for use with ggplot2 scales.
#'
#' Converts supplied breaks into labels using plotmath formatting to generate
#' pretty scientific notation. A prettier version of [scales::label_scientific()].
#'
#' @param digits the maximum number of significant digits to show before the exponent
#' @param prefix string to add before each label
#' @param suffix string to add after each label
#' @param decimal.mark character to be used as the decimal point
#' @param trim whether values should be trimmed or else right-justified to a common width
#' @return A function which takes a numeric vector and converts it to pretty labels
#' @export
label_pretty_scientific <- function(digits = 3, prefix = "", suffix = "",
                                    decimal.mark = ".", trim = TRUE) {

  return(function(x) {
    # handle empty input vector
    if(length(x) == 0) {
      return(character(0))
    }
    # round to the specified maximum number of significant figures
    x <- signif(x, digits = digits)

    # generate formatted output
    output <- format(x, decimal.mark = decimal.mark, trim = trim, #digits = digits,
                     scientific = TRUE)

    # make zero actually zero
    output[x == 0] <- "0"
    # quote the part before the exponent to keep all the digits
    output <- gsub("^(.*)e", "'\\1'e", output)
    # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2)
    output <- gsub("e\\+","e", output)
    # turn the 'e' into plotmath format
    output <- gsub("e", "%*%10^", output)

    # convert 1x10^ or 1.000x10^ -> 10^
    # only do this when ALL values are like this or else 0
    mant_ones <- grepl("^\\'1[\\.0]*\\'", output)
    if(sum(!mant_ones & output != "NA") == 0) {
      output <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", output)
    }

    # add prefix and suffix
    if(nchar(prefix) > 0) {
      output <- paste0('"', prefix, '"~', output)
    }
    if(nchar(suffix) > 0) {
      output <- paste0(output, '~"', suffix, '"')
    }

    # put any NA values back in
    output[is.na(x)] <- NA
    # parse(text = l)
    # return this as an expression
    return(parse(text = output))
  })



}

#' Construct lollipop plots
#'
#' Displays a set of label/value pairs in a lollipop plot.
#' @param labels the text labels to display with each numerical value
#' @param values the numerical values to display
#' @return a ggplot object containing the plot, so you can fine-tune it as
#'   desired before displaying it
#' @details There isn't really any error-checking here. Put garbage in, you'll
#'   get garbage out.
#'
#' @importFrom rlang .data
#' @export
lollipop_plot <- function(labels, values) {
  #TODO: allow switching axes
  #TODO: allow specifying whether or not to sort by values, and the desired order

  df <- tibble::tibble(Label = labels, Value = values)

  gp <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Value, y = .data$Label)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .data$Label), color = "darkgrey", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "darkgrey") +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = .data$Label)) +
    ggplot2::geom_point()

  # if ggthemes is available, add a rangeframe
  if(requireNamespace("ggthemes", quietly = TRUE)) {
    gp <- gp +
      ggthemes::geom_rangeframe(color = "black")
  }

  return(gp)
}
