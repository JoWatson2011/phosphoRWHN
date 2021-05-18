#' Calculate RWHN
#'
#' @param hetNet List; output of constructHetNet()
#' @param transMat List; output of calculateTransitionMatrix()
#' @param vertices Data frame; vertices in heterogeneous network, as containined in output of constructHetNet()
#' @param seeds Character vector of seed nodes
#' @param transitionProb Integer; transition probability
#' @param restart Integer; restart probability
#' @param eps Integer; steady state defined by distance between p0 and p1
#' @param random Logical; randomly permute edges? For control cases
#' @param eta_xy Integer; weighting on protein layer
#' @param eta_yz Integer; weighting on function layer
#' **For visualisation:**
#' @param rwhn_output output of calculateRWHN() (or getRWHN())
#' @param database character; name of annotation database
#' @param colours character vector of length 2 with colours for scale
#' @param removeCommon logical; if rwhn_output is a list, remove functions ranked at the same position in all cases. This is useful for removing noise.
#'
#' @importFrom stats reorder
#' @importFrom igraph bipartite.mapping graph_from_data_frame simplify graph_from_adjacency_matrix get.adjacency graph_from_incidence_matrix get.incidence
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter mutate group_by arrange desc ungroup group_split bind_rows n
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient scale_x_discrete guides guide_colourbar xlab ylab theme  margin element_blank element_rect element_text unit
#'
#' @return A data frame of ranked nodes from the functional layer
#'
#' @export calculateRWHN
calculateRWHN <- function(hetNet,
                          seeds,
                          transitionProb,
                          restart,
                          eta_xy,
                          eta_yz,
                          eps = 1 / 10 ^ 6,
                          random = F,
                          filterFunctions = T) {

  if(class(hetNet) != "list" & all(names(hetNet) %in% c("v", "edgelists"))){
    stop("hetNet must be output from constructHetNet()")
  }

  edgelists <- hetNet$edgelists
  vertices <- hetNet$v

  inter_el <- lapply(edgelists[c("x", "y", "z")], function(i) {
    nw <- i %>%
      igraph::graph_from_data_frame(directed = F) %>%
      igraph::simplify()

    if (random) {
      names <- unique(c(i$from, i$to))
      adj <- matrix(
        data = sample.int(2,
                          length(names) * length(names),
                          TRUE),
        nrow = length(names),
        ncol = length(names),
        dimnames = list(names, names)
      )
      randomnw <- igraph::graph_from_adjacency_matrix(adj)
      mat <- igraph::get.adjacency(randomnw, sparse = T)
    } else{
      mat <- igraph::get.adjacency(nw, sparse = T)
    }
    return(as.matrix(mat))
  })

  intra_el <-
    lapply(edgelists[c("xy", "yx", "yz", "zy")], function(i) {
      nw <- i %>%
        igraph::graph_from_data_frame(directed = F) %>%
        igraph::simplify()
      bimap <- igraph::bipartite.mapping(nw)

      if (bimap[[1]] == T) {
        if (random) {
          incd <- matrix(
            data = sample.int(2,
                              length(unique(i$from)) * length(unique(i$to)),
                              TRUE),
            nrow = length(unique(i$from)),
            ncol = length(unique(i$to)),
            dimnames = list(unique(i$from), unique(i$to))
          )
          randomnw <- igraph::graph_from_incidence_matrix(incd)
          bimap_random <- igraph::bipartite.mapping(randomnw)
          mat <- igraph::get.incidence(
            randomnw,
            types = bimap_random$type,
            attr = NULL,
            names = TRUE,
            sparse = T
          )
        } else{
          mat <- igraph::get.incidence(
            nw,
            types = bimap$type,
            attr = NULL,
            names = TRUE,
            sparse = T
          )
        }
      }
      return(as.matrix(mat))
    })

  heterogenousNetwork <- c(inter_el, intra_el)

  vertices$transition <- ifelse(vertices$layer == "phos",
                                "x",
                                ifelse(
                                  vertices$layer == "prot",
                                  "y",
                                  ifelse(vertices$layer == "func", "z", NA)
                                ))

  transitionMat <-
    calculateTransitionMatrix(transitionProb = transitionProb,
                              vertices = vertices,
                              hetNet = heterogenousNetwork)

  rwhn <- getRWHN(
    transMat = transitionMat,
    restart = restart,
    eta_xy = eta_xy,
    eta_yz = eta_yz,
    seeds = seeds,
    eps = eps
  )

  if(filterFunctions){
    rwhn <- rwhn %>%
    dplyr::filter(.data$name %in% vertices[vertices$layer == "func", ]$v) %>%
    dplyr::mutate(rank = 1:n())
  }

  return(rwhn)
}

#' @export
#' @rdname calculateRWHN
calculateTransitionMatrix <-
  function(transitionProb, hetNet, vertices) {
    xy_tp <- intra_tm(
      intra = hetNet[["xy"]],
      vCols =  vertices[vertices$transition == "y", ]$v,
      vRows = vertices[vertices$transition == "x", ]$v,
      transitionProb = transitionProb
    )
    yx_tp <- t(xy_tp)
    yz_tp <- intra_tm(
      intra = hetNet[["yz"]],
      vCols =  vertices[vertices$transition == "z", ]$v,
      vRows = vertices[vertices$transition == "y", ]$v,
      transitionProb = transitionProb
    )
    zy_tp <- t(yz_tp)

    x_tp <- inter_tm(
      inter = hetNet[["x"]],
      intra_names = colnames(hetNet[["yx"]]),
      transitionProb = transitionProb
    )

    y_tp <- inter_tm(
      inter = hetNet[["y"]],
      intra_names = unique(c(colnames(hetNet[["xy"]]), colnames(hetNet[["zy"]]))),
      transitionProb = transitionProb
    )

    z_tp <- inter_tm(
      inter = hetNet[["z"]],
      intra_names = colnames(hetNet[["yz"]]),
      transitionProb = transitionProb
    )


    # # create the full transition matrix
    # # |siteMatrix (x)   Site2Prot (yx)         0         |
    # # |Prot2Site (xy)   protMatrix (y)   prot2func (zy)  |
    # # |     0           func2prot (yz)   funcMatrix(z)   |

    top0 <- matrix(
      0,
      nrow = nrow(x_tp),
      ncol = ncol(z_tp),
      dimnames = list(rownames(x_tp), colnames(z_tp))
    )
    btm0 <- matrix(
      0,
      nrow = nrow(z_tp),
      ncol = ncol(x_tp),
      dimnames = list(rownames(z_tp), colnames(x_tp))
    )

    yz_tp <-
      yz_tp[, colnames(yz_tp)[colnames(yz_tp) %in% colnames(z_tp)]]
    zy_tp <-
      zy_tp[rownames(zy_tp)[rownames(zy_tp) %in% rownames(z_tp)], ]

    # tmp1 <- cbind(x_tp, xy_tp[rownames(x_tp), colnames(y_tp), drop = F], top0)
    # tmp2 <- cbind(yx_tp[rownames(y_tp), colnames(x_tp), drop = F], y_tp, yz_tp[rownames(y_tp), colnames(z_tp), drop = F])
    # tmp3 <- cbind(btm0, zy_tp[rownames(z_tp), colnames(y_tp), drop = F], z_tp)
    tmp1 <- cbind(x_tp, xy_tp, top0)
    tmp2 <- cbind(yx_tp, y_tp, yz_tp)
    tmp3 <- cbind(btm0, zy_tp, z_tp)

    M1 <- rbind(tmp1, tmp2, tmp3)

    return(list(
      M1 = M1,
      x_tp = x_tp,
      y_tp = y_tp,
      z_tp = z_tp
    ))
  }

#' @export
#' @rdname calculateRWHN
getRWHN <- function(transMat,
                    restart,
                    eta_xy,
                    eta_yz,
                    seeds,
                    eps) {
  M1 <- transMat[["M1"]]
  x <- transMat[["x_tp"]]
  y <- transMat[["y_tp"]]
  z <- transMat[["z_tp"]]

  seedScores <- as.data.frame(rownames(x))
  seedScores$Scores <- ifelse(seedScores[, 1] %in% seeds,
                              1 / sum(seedScores[, 1] %in% seeds),
                              0)

  probabilityVector <- c(seedScores$Score,
                         rep.int((1 / nrow(M1) * eta_xy),
                                 nrow(y)), rep.int((1 / nrow(M1) * eta_yz),
                                                   nrow(z)))

  iter <- 0
  p0 <- probabilityVector
  p1 <- c(rep(0, length(probabilityVector)))

  while (sum(abs(p0 - p1)) > eps) {
    p0 <- p1
    p1 <-
      ((1 - restart) * t(M1)) %*% p0 + restart * probabilityVector
    p1 <- p1 / sum(p1)

    if (!any(is.na(p1))) {
      iter <- iter + 1
    } else{
      p1 <- p0
      break
    }
  }

  pi <- data.frame(
    V1 = p1[, 1],
    name = c(rownames(x), rownames(y), rownames(z)),
    stringsAsFactors = F
  )
  pi <- pi[order(pi[, 1], decreasing = T), ]

  pi$rank <- 1:nrow(pi)

  return(pi)
}

#' @export
#' @rdname calculateRWHN

heatmap_RWHN <- function(rwhn_output,
                         database = "",
                         colours = c(low = "#e6e4f8",
                                     high = "#46009e"),
                         removeCommon = T) {
  if (class(rwhn_output) == "list") {
    if (is.null(names(rwhn_output))) {
      loop <- 1:length(rwhn_output)
    } else{
      loop <- names(rwhn_output)
    }
    rwhn_flt <-  lapply(loop, function(i) {
      rwhn_output[[i]] %>%
        dplyr::mutate(seed = i,
                      rank = 1:n())
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(.data$name) %>%
      dplyr::mutate(
        rank_dif = abs(.data$rank - mean(.data$rank)),
        color = ifelse(.data$rank_dif == 0, T, NA),
        V1 = signif(.data$V1, digits = 2)
      ) %>%
      dplyr::arrange(desc(.data$V1))# %>%
      #filter(rank_dif > (0.05 / max(rank_dif))) %>%

      # {
      #   if (removeCommon)
      #
      #   else
      #     .data
      # } %>%

    if(removeCommon){# Filter terms that appear in the same position in all conditions
      rwhn_flt <- dplyr::filter(rwhn_flt, .data$rank_dif > 1)
    }

    rwhn_flt <- rwhn_flt %>%
      dplyr::ungroup() %>%
      dplyr::group_split(.data$seed) %>%
      sapply(function(i) {
        if (nrow(i > 0)) {
          pct <-
            i[1, ]$V1                             # Filter top 5% of terms ( with a messy loop!!)
          df <- i[1, ]
          if (nrow(i) > 1) {
            for (x in 2:nrow(i)) {
              if (pct < 0.05) {
                df <- rbind(df, i[x, ])
                pct <- pct + i[x, ]$V1
              } else{
                break
              }
            }
          }
          df$rank <- 1:nrow(df)
        } else{
          df <- data.frame()
        }

        return(df)
      }, simplify = F, USE.NAMES = T) %>%
      dplyr::bind_rows() %>%
      dplyr::filter(!grepl("0\\.", .data$name))
  } else{
    rwhn_flt <- rwhn_output %>%
      dplyr::mutate(rank = 1:n(),
                    seed = "") %>%
      dplyr::filter(!grepl("0\\.", .data$name))

    if (nrow(rwhn_flt > 0)) {
      pct <-
        rwhn_flt[1, ]$V1                             # Filter top 5% of terms ( with a messy loop!!)
      df <- rwhn_flt[1, ]
      if (nrow(rwhn_flt) > 1) {
        for (x in 2:nrow(rwhn_flt)) {
          if (pct < 0.05) {
            df <- rbind(df, rwhn_flt[x, ])
            pct <- pct + rwhn_flt[x, ]$V1
          } else{
            break
          }
        }
      }
      df$rank <- 1:nrow(df)
    } else{
      df <- data.frame()
    }
    rwhn_flt <- dplyr::arrange(df, rank)

  }

  sighm <-
    ggplot2::ggplot(rwhn_flt, ggplot2::aes(x = as.factor(.data$seed),
                                           y = reorder(.data$name, .data$rank))) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$rank), colour = "white") +
    ggplot2::scale_fill_gradient(
      breaks = seq(
        from = 1,
        to = max(rwhn_flt$rank),
        by = 5
      ),
      low = colours[1],
      high = colours[2]
    ) +
    ggplot2::guides(color = FALSE,
                    fill =  ggplot2::guide_colourbar(title = "Rank", reverse = T)) +
    ggplot2::xlab("Seed nodes") +
    ggplot2::ylab(database) +
    ggplot2::theme(
      legend.key.size = ggplot2::unit(.25, "cm"),
      legend.title = ggplot2::element_text(size = 5),
      legend.text = ggplot2::element_text(size = 5),
      axis.text.y = ggplot2::element_text(size = 5),
      axis.title = ggplot2::element_text(size = 6),
      panel.background = ggplot2::element_rect(fill = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(0, 0, 0, 0, "cm")
    ) +
    scale_x_discrete(position = "top")

  return(sighm)
}




### -- Helper functions

inter_tm <- function(inter, intra_names, transitionProb) {
  ### Which nodes in inter have degree 0 in intra ?

  notHetNames <-
    colnames(inter)[!(colnames(inter) %in% intra_names)]
  notHet <-
    inter[rownames(inter) %in% notHetNames, colnames(inter) %in% notHetNames]

  # eq2 when (B)i,j = 0
  if (class(notHet)[1] %in% c("data.frame", "matrix")) {
    notHet_tp <- (notHet) / rowSums(notHet)
  } else{
    notHet_tp <- data.frame()
  }

  # Which nodes are in inter and intra?
  hetNames <- colnames(inter)[colnames(inter) %in% intra_names]
  het <-
    inter[rownames(inter) %in% hetNames, colnames(inter) %in% hetNames]

  # eq2 when (B)i,j > 0
  if (class(het)[1] %in% c("data.frame", "matrix")) {
    het_tp <- ((1 - transitionProb) * het) / rowSums(het)
  } else if (length(hetNames) != 0) {
    het_tp <- data.frame()
    for (i in 1:length(hetNames)) {
      het_tp[1, i] <- 1 - transitionProb
      colnames(het_tp)[i] <- rownames(het_tp)[i] <- hetNames[i]
    }
  } else {
    het_tp <- data.frame()
  }

  # which nodes aren't in inter but are in intra?
  onlyHetNames <- intra_names[!(intra_names %in% colnames(inter))]
  onlyHet <- matrix(
    rep.int(0, (length(onlyHetNames) ^ 2)),
    nrow = length(onlyHetNames),
    ncol = (length(onlyHetNames)),
    dimnames = list(onlyHetNames, onlyHetNames)
  )

  # combine
  tp <- data.frame()

  if (ncol(notHet_tp) != 0 & ncol(het_tp) != 0) {
    row <- nrow(notHet_tp)
    tp[1:row, 1:row] <- notHet_tp
    row <- row + 1

    tp[row:(row + nrow(het_tp) - 1), row:(row + nrow(het_tp) - 1)] <-
      het_tp
    row <- row + nrow(het_tp)

    if (ncol(onlyHet) != 0) {
      tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <-
        onlyHet

      colnames(tp) <-
        c(colnames(notHet_tp),
          colnames(het_tp),
          colnames(onlyHet))
      rownames(tp) <-
        c(rownames(notHet_tp),
          rownames(het_tp),
          rownames(onlyHet))

    } else{
      colnames(tp) <- c(colnames(notHet_tp), colnames(het_tp))
      rownames(tp) <- c(rownames(notHet_tp), rownames(het_tp))
    }
  } else if (ncol(notHet_tp) == 0 & ncol(het_tp) != 0) {
    row <- nrow(het_tp)

    tp[1:row, 1:row] <- het_tp
    row <- row + 1

    if (ncol(onlyHet) != 0) {
      tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <-
        onlyHet
      colnames(tp) <- c(colnames(het_tp), colnames(onlyHet))
      rownames(tp) <- c(rownames(het_tp), rownames(onlyHet))
    } else {
      colnames(tp) <- colnames(het_tp)
      rownames(tp) <- rownames(het_tp)
    }
  } else if (ncol(notHet_tp) != 0 & ncol(het_tp) == 0) {
    row <- nrow(notHet_tp)

    tp[1:row, 1:row] <- notHet_tp
    row <- row + 1

    if (ncol(onlyHet) != 0) {
      tp[row:(row + nrow(onlyHet) - 1), row:(row + nrow(onlyHet) - 1)] <-
        onlyHet
      colnames(tp) <- c(colnames(notHet_tp), colnames(onlyHet))
      rownames(tp) <- c(rownames(notHet_tp), rownames(onlyHet))
    } else {
      colnames(tp) <- colnames(notHet_tp)

      rownames(tp) <- rownames(notHet_tp)
    }
  } else{
    print("Situation not accounted for")
  }

  tp <- as.matrix(tp)
  tp[is.na(tp)] <- 0
  tp <- tp[order(match(rownames(tp),
                       intra_names)),
           order(match(colnames(tp),
                       intra_names))]

  return(tp)
}

intra_tm <- function(intra, vCols, vRows, transitionProb) {
  # Which nodes have degree > 0 in bipartite network?
  present <- intra[rowSums(intra) > 0, ]

  # eq1
  present_tp <- (transitionProb * present) / rowSums(present)

  # rbind with nodes of degree = 0 in bipartite network
  tp <- rbind(present, intra[rowSums(intra) == 0, ])

  # return(tp)

  # add rows for nodes in inter but not intra networks
  rows <- vRows[!(vRows %in% rownames(intra))]
  cols <- vCols[!(vCols %in% colnames(intra))]

  fin_tp <- matrix(0,
                   nrow = (nrow(tp) + length(rows)),
                   ncol = (ncol(tp) + length(cols)))
  fin_tp[1:nrow(tp), 1:ncol(tp)] <- tp

  colnames(fin_tp) <- c(colnames(tp), cols)
  rownames(fin_tp) <- c(rownames(tp), rows)

  return(fin_tp)
}
