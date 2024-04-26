#' Parse a chemical formula
#'
#' Given a chemical formula string, this function extracts the frequencies of
#' all elements in the formula. Note that formulas with other text besides
#' element abbreviations and numbers will likely not parse correctly. Element
#' abbreviations MUST be properly capitalized.
#'
#' @param cf A character vector with one or more formulas.
#'
#' @return A tibble with one row per formula and one column per element
#'   containing the total number of each element per formula.
#'
#' @importFrom rlang .data
#' @export
parse_chemical_formula <- function(cf) {
  # extract element/number pairs
  el_num <- stringr::str_match_all(cf, pattern = "(?:(?<element>[A-Z][a-z]?)(?<number>\\d*))")
  # construct vectors
  el_num <- lapply(X = el_num, FUN = function(el) {
    # if only NAs, just return "zero Carbon atoms"
    if(sum(!is.na(el[, "number"])) == 0) {
      return(c("C" = 0))
    }
    # convert strings to numbers
    v <- as.integer(el[, "number"])
    # NA values mean there is one of that element
    v[is.na(v)] <- 1
    # add element names
    names(v) <- el[, "element"]
    # collapse any duplicate elements, e.g. COOH or similar
    v <- tapply(X = v, INDEX = names(v), FUN = sum, default = 0)
    return(v)
  })

  # combine all results
  el_num <- dplyr::bind_rows(el_num)
  # NA values should be zero
  el_num <- dplyr::mutate(el_num, dplyr::across(.cols = dplyr::everything(),
                                          .fns = ~ tidyr::replace_na(.x, 0)))
  return(el_num)
}
