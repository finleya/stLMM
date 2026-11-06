test_that("CAR component labels handle empty graphs and reset between components", {
  expect_identical(
    stLMM:::car_connected_components(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), dims = c(0L, 0L)
    )),
    integer(0)
  )
  expect_identical(
    stLMM:::car_connected_components(Matrix::sparseMatrix(
      i = integer(0), j = integer(0), dims = c(4L, 4L)
    )),
    1:4
  )

  # A cycle, a separate pair, and an isolated node, interleaved by node id.
  from <- c(1L, 3L, 5L, 2L)
  to <- c(3L, 5L, 1L, 6L)
  adj <- Matrix::sparseMatrix(
    i = c(from, to), j = c(to, from), x = 1, dims = c(6L, 6L)
  )
  expect_identical(stLMM:::car_connected_components(adj), c(1L, 2L, 1L, 3L, 1L, 2L))
})

test_that("CAR component labels agree with transitive closure for all five-node graphs", {
  n <- 5L
  pairs <- utils::combn(n, 2L)
  for(mask in 0:(2^ncol(pairs) - 1L)){
    keep <- as.logical(intToBits(mask)[seq_len(ncol(pairs))])
    from <- pairs[1L, keep]
    to <- pairs[2L, keep]
    adj <- Matrix::sparseMatrix(
      i = c(from, to), j = c(to, from), x = rep(1, 2L * length(from)),
      dims = c(n, n)
    )

    # Independent reference: connected nodes are mutually reachable.
    reachable <- as.matrix(adj) != 0
    diag(reachable) <- TRUE
    for(k in seq_len(n))
      reachable <- reachable | outer(reachable[, k], reachable[k, ], "&")
    first_node <- apply(reachable, 1L, function(x) which(x)[1L])
    expected <- match(first_node, unique(first_node))

    expect_identical(stLMM:::car_connected_components(adj), expected,
                     info = paste("edge mask", mask))
  }
})

test_that("CAR graph validation still rejects isolated and disconnected areas", {
  pair <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(2L, 1L), x = 1,
                               dims = c(2L, 2L))
  expect_error(car_graph(Matrix::bdiag(pair, pair)), "disconnected components")
  expect_error(car_graph(Matrix::bdiag(pair, Matrix::Diagonal(1L, 0))), "isolated area")

  # A high-degree hub exercises a broad queue; explicit zero edges are ignored.
  n <- 100L
  adj <- Matrix::sparseMatrix(
    i = c(rep(1L, n - 1L), 2:n, 2L, 3L),
    j = c(2:n, rep(1L, n - 1L), 3L, 2L),
    x = c(rep(1, 2L * (n - 1L)), 0, 0), dims = c(n, n),
    dimnames = list(as.character(seq_len(n)), as.character(seq_len(n)))
  )
  expect_identical(stLMM:::car_connected_components(adj), rep(1L, n))
  graph <- car_graph(adj)
  expect_identical(graph$degree, c(n - 1, rep(1, n - 1L)))
  expect_equal(graph$adjacency, Matrix::drop0(adj))
})
