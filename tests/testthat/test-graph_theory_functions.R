test_that("k-plex finding works", {
  min_kc <- 3
  kp <- 2
  test_graph <- igraph::graph_from_edgelist(el = matrix(c("v0", "v1",
                                                          "v0", "v2",
                                                          "v1", "v2",
                                                          "v1", "v4",
                                                          "v2", "v3",
                                                          "v2", "v4",
                                                          "v3", "v4",
                                                          "v4", "v5",
                                                          "v5", "v6",
                                                          "v5", "v7",
                                                          "v5", "v8",
                                                          "v6", "v7",
                                                          "v6", "v8",
                                                          "v7", "v8",
                                                          "v6", "v9",
                                                          "v7", "v9",
                                                          "v8", "v9"
  ),
  ncol = 2,
  byrow = TRUE),
  directed = FALSE)

  test_res <- find_kplexes(test_graph, min_kc = min_kc, kp = kp)

  # check that the found plex has minimum degree == size-kp
  expect_equal(test_res$plex_min_deg[[1]], test_res$plex_size[[1]] - kp)
  # check that it found the right plex
  expect_setequal(as.character(igraph::V(test_res$plex_list[[1]])),
                  c("1", "2", "3", "4", "5"))
})
