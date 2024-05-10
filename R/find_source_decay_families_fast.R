#' Simple source decay family search
#'
#' Some compounds fragment when ionized into the gas phase, forming a source
#' decay family. Mass spec features in these families are essentially duplicates
#' of the parent compound and should be identified and treated specially. This
#' function identifies such families based on very similar retention times and
#' correlated peak areas across all samples. This approach also finds instances
#' where the alignment software has missed adducts, isotopologues, or multiply
#' charged versions of a single compound. It also corrects for occasional
#' alignment software errors.
#'
#' @param feature_rts Named numeric vector with retention times for all features.
#'   Names should be unique identifiers for each feature.
#' @param feature_areas Numeric matrix with one column per feature and one row
#'   per sample. Column names must match the names in \code{feature_rts}. Peak areas
#'   should be log-transformed to reduce the extreme skewness typically associated
#'   with mass spec peak intensities.
#' @param drt_max Maximum retention time difference (in minutes) to consider.
#' @param area_cor_min Minimum peak area correlation to consider.
#'
#' @details
#' This is a heavily simplified version of \link{find_source_decay_families} that
#' only considers retention time and peak area correlation. It is thus much faster,
#' but also potentially prone to make mistakes. Be sure to check the outputs to
#' ensure that SDF calls seem reasonable for the desired downstream application.
#'
#' @return A list containing:\tabular{ll}{
#'    \code{matched_compound_pairs} \tab Tibble of all matching pairs of features \cr
#'    \code{sdf_graph} \tab Network object constructed from feature pairs. \cr
#'    \code{sdf_cluster_members} \tab Tibble of cluster member IDs and associated cluster IDs. \cr
#'    \code{sdf_cluster_centers} \tab Tibble containing the central member of each compound cluster. \cr
#'    \tab \cr
#'}
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#' @importFrom rlang .env
#' @export
find_source_decay_families_fast <- function(feature_rts, feature_areas,
                                            drt_max = 0.01,
                                            area_cor_min = 0.8) {
  # pull feature IDs from RT vector
  feature_ids <- names(feature_rts)
  if(is.null(feature_ids) || length(unique(feature_ids)) < length(feature_ids)) {
    stop(simpleError("feature_rts must have unique names for each feature!"))
  }
  # make sure IDs match in feature_areas
  if(length(setdiff(union(feature_ids, colnames(feature_areas)), feature_ids)) > 0) {
    stop(simpleError("Column names in feature_areas must match the names in feature_rts!"))
  }

  # don't assume that feature_rts and feature_areas have features in same order
  rt2area_idx <- match(feature_ids, colnames(feature_areas))

  # Find similar retention time pairs ----

  delta_rt <- abs(outer(feature_rts, feature_rts, FUN = "-"))
  delta_rt[lower.tri(delta_rt, diag = TRUE)] <- NA
  rt_matches <- which(delta_rt < drt_max, arr.ind = TRUE)

  # Find pairs that also have correlated peak areas ----

  df2 <- mapply(idx1 = rt_matches[, "row"], idx2 = rt_matches[, "col"], FUN = function(idx1, idx2) {
    # convert to feature area index
    idx1 <- rt2area_idx[idx1]
    idx2 <- rt2area_idx[idx2]
    return(stats::cor(x = feature_areas[, idx1], y = feature_areas[, idx2]))
  })

  area_cor_matches <- tibble::tibble(Name.1 = feature_ids[rt_matches[, "row"]],
                                     Name.2 = feature_ids[rt_matches[, "col"]],
                                     dRT = delta_rt[rt_matches],
                                     Correlation = df2) %>%
    dplyr::filter(.data$Correlation > area_cor_min)

  # Generate correlation graph ----

  el <- area_cor_matches %>%
    dplyr::select("Name.1", "Name.2") %>%
    as.matrix()

  sdf_graph <- igraph::graph_from_edgelist(el = el, directed = FALSE)
  # edges are in the order we fed them in, so we can apply the attributes directly
  igraph::E(sdf_graph)$weight <- area_cor_matches$Correlation
  igraph::E(sdf_graph)$rt_diff <- area_cor_matches$dRT
  igraph::V(sdf_graph)$max_area <- apply(feature_areas[, igraph::V(sdf_graph)$name],
                                         MARGIN = 2, FUN = max)
  igraph::V(sdf_graph)$label <- igraph::V(sdf_graph)$name

  # connected components are source decay families
  sdf_comps <- igraph::components(sdf_graph)
  sdf_cluster_members <- tibble::enframe(sdf_comps$membership,
                                         name = "Name", value = "ClusterID")
  sdf_clusters <- lapply(X = seq.int(from = 1, to = sdf_comps$no, by = 1), FUN = function(comp) {
    # get members
    comp_idx <- which(sdf_comps$membership == comp)
    # get component
    comp_graph <- igraph::induced_subgraph(sdf_graph, vids = comp_idx)
    return(comp_graph)
  })

  # calculate eigenvector centrality for all vertices
  for(comp in sdf_clusters) {
    # calculate eigenvector centrality
    comp_eigen_centrality <- igraph::eigen_centrality(comp)$vector
    # apply to the full graph
    igraph::V(sdf_graph)[names(comp_eigen_centrality)]$eigen_centrality <- comp_eigen_centrality
  }

  sdf_cluster_members <- sdf_cluster_members %>%
    dplyr::mutate(EigenCentrality = igraph::V(sdf_graph)$eigen_centrality,
                  MaxArea = igraph::V(sdf_graph)$max_area)

  # calculate centers using eigenvector centrality
  # break ties with max intensity
  sdf_eigen_centers <- sdf_cluster_members %>%
    dplyr::group_by(.data$ClusterID) %>%
    dplyr::slice_max(order_by = tibble::tibble(.env$EigenCentrality, .env$MaxArea), n = 1, with_ties = FALSE)

  return(list(matched_compound_pairs = area_cor_matches,
              sdf_graph = sdf_graph,
              sdf_cluster_members = sdf_cluster_members,
              sdf_cluster_centers = sdf_eigen_centers))
}
