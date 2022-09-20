#' Construct fancy-looking ROC plot
#'
#' Generates a ggplot object to display a ROC curve.
#'
#' @param data output from the [calc_roc] function
#' @return a ggplot object containing the plot, so you can fine-tune it as
#'   desired before displaying it
#' @details There isn't really any error-checking here. Put garbage in, you'll
#'   get garbage out.
#' @export
fancy_roc_plot <- function(data) {
  return(
    ggplot(data = data$curves, mapping = aes(x = (100 - specificity), y = sensitivity)) +
      geom_ribbon(mapping = aes(x = (100 - specificity),
                                ymin = `cisens2.5`, ymax = `cisens97.5`),
                  fill = "grey80") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_path(color = "black", size = 1.0) +
      geom_rangeframe() +
      labs(x = "false-positive rate (%)", y = "true-positive rate (%)") +
      annotate(geom = "text", x = 70, y = 40,
               label = sprintf(fmt = "%1.1f%%\n(%1.1f%% - %1.1f%%)",
                               data$auc, data$ci[1], data$ci[3]),
               hjust = "center", size = 6)
  )
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

#' Construct lollipop plots
#'
#' Displays a set of label/value pairs in a lollipop plot.
#' @param labels the text labels to display with each numerical value
#' @param values the numerical values to display
#' @return a ggplot object containing the plot, so you can fine-tune it as
#'   desired before displaying it
#' @details There isn't really any error-checking here. Put garbage in, you'll
#'   get garbage out.
#' @export
lollipop_plot <- function(labels, values) {
  #TODO: allow switching axes
  #TODO: allow specifying whether or not to sort by values, and the desired order

  df <- tibble::tibble(Label = labels, Value = values)

  gp <- ggplot2::ggplot(df, ggplot2::aes(x = Value, y = Label)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = Label), color = "darkgrey", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, color = "darkgrey") +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = Label)) +
    ggplot2::geom_point()

  # if ggthemes is available, add a rangeframe
  if(requireNamespace("ggthemes", quietly = TRUE)) {
    gp <- gp +
      ggthemes::geom_rangeframe(color = "black")
  }

  return(gp)
}
