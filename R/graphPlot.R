utils::globalVariables(c("degree", "id", "x", "xend", "y", "yend"))

graph_plot_requires <- function(){
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("error: package 'ggplot2' is required for graph plotting")
}

graph_plot_layout <- function(ids, geometry = NULL){
  n <- length(ids)

  if(!is.null(geometry)){
    if(!requireNamespace("sf", quietly = TRUE))
      stop("error: package 'sf' is required to plot graph geometry")
    pts <- suppressWarnings(sf::st_point_on_surface(geometry))
    coords <- sf::st_coordinates(pts)
    return(data.frame(
      id = ids,
      x = coords[, 1],
      y = coords[, 2],
      stringsAsFactors = FALSE
    ))
  }

  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  data.frame(
    id = ids,
    x = cos(theta),
    y = sin(theta),
    stringsAsFactors = FALSE
  )
}

graph_plot_edges <- function(layout, from, to, directed = FALSE,
                             bridge = rep(FALSE, length(from))){
  from_match <- match(as.character(from), layout$id)
  to_match <- match(as.character(to), layout$id)
  keep <- !is.na(from_match) & !is.na(to_match)

  data.frame(
    x = layout$x[from_match[keep]],
    y = layout$y[from_match[keep]],
    xend = layout$x[to_match[keep]],
    yend = layout$y[to_match[keep]],
    directed = directed,
    bridge = bridge[keep],
    stringsAsFactors = FALSE
  )
}

graph_plot_arrow_backoff <- function(edges, layout){
  if(is.null(edges) || !nrow(edges))
    return(edges)

  x_range <- range(layout$x, finite = TRUE)
  y_range <- range(layout$y, finite = TRUE)
  span <- sqrt(diff(x_range)^2 + diff(y_range)^2)
  if(!is.finite(span) || span <= 0)
    return(edges)

  dx <- edges$xend - edges$x
  dy <- edges$yend - edges$y
  len <- sqrt(dx^2 + dy^2)
  keep <- is.finite(len) & len > 0
  if(!any(keep))
    return(edges)

  gap <- 0.018 * span
  frac <- pmin(gap / len[keep], 0.35)
  edges$xend[keep] <- edges$xend[keep] - frac * dx[keep]
  edges$yend[keep] <- edges$yend[keep] - frac * dy[keep]

  edges
}

plot_graph_base <- function(ids,
                            edges,
                            geometry = NULL,
                            degree = NULL,
                            title = NULL,
                            show_ids = FALSE,
                            show_nodes = TRUE,
                            polygon_fill = "grey95",
                            polygon_color = "white",
                            edge_color = "grey25",
                            bridge_color = "#c44e52",
                            node_color = "#1f2937",
                            arrow = FALSE,
                            edge_alpha = 0.65,
                            node_size = 1.8,
                            ...){
  graph_plot_requires()

  layout <- graph_plot_layout(ids, geometry)
  if(!is.null(degree))
    layout$degree <- as.numeric(degree)

  p <- ggplot2::ggplot()

  if(!is.null(geometry)){
    sf_dat <- sf::st_sf(
      id = ids,
      degree = if(is.null(degree)) NA_real_ else as.numeric(degree),
      geometry = geometry
    )
    if(is.null(degree)){
      p <- p + ggplot2::geom_sf(
        data = sf_dat,
        fill = polygon_fill,
        color = polygon_color,
        linewidth = 0.25
      )
    } else {
      p <- p + ggplot2::geom_sf(
        data = sf_dat,
        ggplot2::aes(fill = degree),
        color = polygon_color,
        linewidth = 0.25
      ) +
        ggplot2::scale_fill_viridis_c(name = "degree", option = "C")
    }
  }

  if(!is.null(edges) && nrow(edges)){
    if(isTRUE(arrow))
      edges <- graph_plot_arrow_backoff(edges, layout)
    normal_edges <- edges[!edges$bridge, , drop = FALSE]
    bridge_edges <- edges[edges$bridge, , drop = FALSE]
    edge_arrow <- if(isTRUE(arrow)) {
      grid::arrow(length = grid::unit(0.08, "inches"), type = "closed")
    } else {
      NULL
    }

    if(nrow(normal_edges)){
      p <- p + ggplot2::geom_segment(
        data = normal_edges,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = edge_color,
        alpha = edge_alpha,
        linewidth = 0.35,
        arrow = edge_arrow
      )
    }
    if(nrow(bridge_edges)){
      p <- p + ggplot2::geom_segment(
        data = bridge_edges,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = bridge_color,
        alpha = 0.9,
        linewidth = 0.55,
        linetype = "22",
        arrow = edge_arrow
      )
    }
  }

  if(is.null(geometry) && !is.null(degree)){
    p <- p + ggplot2::geom_point(
      data = layout,
      ggplot2::aes(x = x, y = y, color = degree),
      size = node_size
    ) +
      ggplot2::scale_color_viridis_c(name = "degree", option = "C")
  } else if(isTRUE(show_nodes)) {
    p <- p + ggplot2::geom_point(
      data = layout,
      ggplot2::aes(x = x, y = y),
      color = node_color,
      size = node_size
    )
  }

  if(isTRUE(show_ids)){
    p <- p + ggplot2::geom_text(
      data = layout,
      ggplot2::aes(x = x, y = y, label = id),
      size = 2.8,
      vjust = -0.7
    )
  }

  if(!is.null(geometry)){
    p <- p + ggplot2::coord_sf()
  } else {
    p <- p + ggplot2::coord_equal()
  }

  p <- p +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0, face = "bold")
    )

  p
}

plot.stLMM_car_graph <- function(x,
                                 show_ids = FALSE,
                                 show_nodes = TRUE,
                                 show_bridges = TRUE,
                                 color_by_degree = FALSE,
                                 ...){
  adj <- Matrix::drop0(x$adjacency)
  trip <- Matrix::summary(adj)
  trip <- trip[trip$i < trip$j, , drop = FALSE]

  bridge <- rep(FALSE, nrow(trip))
  if(isTRUE(show_bridges) && nrow(x$island_added_edges)){
    pair_key <- function(a, b){
      a <- as.character(a)
      b <- as.character(b)
      paste(pmin(a, b), pmax(a, b), sep = "\r")
    }
    edge_key <- pair_key(x$ids[trip$i], x$ids[trip$j])
    bridge_key <- pair_key(x$island_added_edges$from, x$island_added_edges$to)
    bridge <- edge_key %in% bridge_key
  }

  layout <- graph_plot_layout(x$ids, x$geometry)
  edges <- graph_plot_edges(
    layout,
    from = x$ids[trip$i],
    to = x$ids[trip$j],
    bridge = bridge
  )

  plot_graph_base(
    ids = x$ids,
    edges = edges,
    geometry = x$geometry,
    degree = if(isTRUE(color_by_degree)) x$degree else NULL,
    title = "Areal adjacency graph",
    show_ids = show_ids,
    show_nodes = show_nodes,
    arrow = FALSE,
    ...
  )
}

plot.stLMM_graph <- function(x,
                             show_ids = FALSE,
                             show_nodes = TRUE,
                             color_by_degree = FALSE,
                             ...){
  graph_type <- x$graph_type %||% "unknown"

  if(graph_type %in% c("dagar", "dagar_time")){
    ids <- as.character(x$ids)
    geometry <- x$geometry %||% NULL
    layout <- graph_plot_layout(ids, geometry)

    from <- character(0)
    to <- character(0)
    for(i in seq_along(x$parent_count)){
      m <- x$parent_count[i]
      if(m > 0L){
        idx <- seq.int(x$parent_start[i] + 1L, length.out = m)
        parents <- x$parent_index[idx] + 1L
        from <- c(from, ids[parents])
        to <- c(to, rep(ids[i], m))
      }
    }
    edges <- graph_plot_edges(layout, from = from, to = to, directed = TRUE)

    return(plot_graph_base(
      ids = ids,
      edges = edges,
      geometry = geometry,
      degree = if(isTRUE(color_by_degree)) x$degree else NULL,
      title = paste0(
        if(graph_type == "dagar_time") "Ordered DAGAR-time graph" else "Ordered DAGAR graph",
        " (",
        x$ordering %||% "unknown",
        ")"
      ),
      show_ids = show_ids,
      show_nodes = show_nodes,
      arrow = TRUE,
      ...
    ))
  }

  if(graph_type %in% c("car", "car_time")){
    ids <- as.character(x$ids)
    geometry <- x$geometry %||% NULL
    layout <- graph_plot_layout(ids, geometry)
    edges <- graph_plot_edges(
      layout,
      from = ids[x$edge_i],
      to = ids[x$edge_j]
    )

    return(plot_graph_base(
      ids = ids,
      edges = edges,
      geometry = geometry,
      degree = if(isTRUE(color_by_degree)) x$degree else NULL,
      title = paste0(toupper(graph_type), " graph"),
      show_ids = show_ids,
      show_nodes = show_nodes,
      arrow = FALSE,
      ...
    ))
  }

  stop("error: graph plotting is implemented for CAR, CAR-time, DAGAR, and DAGAR-time graphs")
}
