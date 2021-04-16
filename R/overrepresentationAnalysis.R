#' Overrepresentation analysis
#'
#' @param clustering A named vector; names correspond to sites in format *Gene name*_*AminoAcid+Position* (e.g. MAPK1_T185), values correspond to cluster. This is
#' @param database Annotation database to calculate enrichment from. See [enrichR website](https://maayanlab.cloud/Enrichr/#stats)
#' @param simplify If database = "GO_Biological_Process_2018", remove redundant terms
#' @param vis logical; visualise output?
#' @param colours if vis = T, a character vector of length 2 for colour scale
#' @param RWHN_sig optional; results of RWHN if wanting to perform comparison,
#'
#' @importFrom dplyr arrange desc mutate filter select n bind_rows
#' @importFrom tidyr separate
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient guides guide_colourbar theme unit element_text element_rect element_blank margin scale_x_discrete xlab ylab geom_point
#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
#'
#' @return If vis = T, a ggplot object. If vis = F, a data frame containing enriched terms.
#' @export

overrepresentationAnalysis <- function(clustering,
                                       database = "GO_Biological_Process_2018",
                                       simplify = T,
                                       vis = F,
                                       colours = c(low = "#439981", high = "#b2d4cd"),
                                       RWHN_sig = NULL
                                       ){


  enrichedTerms <- lapply(unique(clustering)[order(unique(clustering))], function(i){
    ids <- names(clustering[clustering == i ] )
    cl_prots <- unique(gsub("_.*", "", ids))
    cl_prots <- cl_prots[cl_prots != ""]

    enriched <- getEnrichr(cl_prots, database, 0.05) %>%
      dplyr::mutate(cluster = i)

    return(enriched)
  }) %>%
    dplyr::bind_rows()

  if(grepl("GO_", database)){
    enrichedTerms <-  enrichedTerms %>%
      tidyr::separate(.data$Term,
               into = c("Term", "GOID"),
               sep = " \\(",
               extra = "drop") %>%
      dplyr::mutate(GOID = sub("\\)",
                        "",
                        .data$GOID))
  }


  if(grepl("GO_Biological_Process", database) & simplify){

    enrichedTerms_flt <- lapply(unique(clustering)[order(unique(clustering))],
                                function(i){
      df <- enrichedTerms[enrichedTerms$cluster == i,]

      keepID <- simplifyGO(GOID = df$GOID, simplifyData = simpleGO)

      df_flt <- df %>%
        dplyr::filter(.data$GOID %in% keepID) %>%
        dplyr::arrange(dplyr::desc(.data$Adjusted.P.value))
      if(nrow(df_flt) > 0){
        df_flt <- dplyr::mutate(df_flt, rank = 1:n())
      }
      return(df_flt)
    }) %>%
      bind_rows() %>%
      dplyr::mutate(V1 = signif(.data$Adjusted.P.value, digits = 2),
             name = factor(.data$Term, unique(.data$Term)))

  } else {

    enrichedTerms_flt <- lapply(unique(clustering)[order(unique(clustering))], function(i){
      df <- enrichedTerms[enrichedTerms$cluster ==i, c("Term", "Adjusted.P.value", "cluster")] %>%
        dplyr::arrange(desc(.data$Adjusted.P.value))

      if(nrow(df) > 0){
        df <- dplyr::mutate(df, rank = 1:n())
      }

      return(df)
    }) %>%
      do.call(rbind, .data) %>%
      mutate(V1 = signif(.data$Adjusted.P.value, digits = 2),
             name = factor(.data$Term, unique(.data$Term)))

  }

  if(!is.null(RWHN_sig)){
    enrichedTerms_flt$Term <- tolower(enrichedTerms_flt$Term)
    enrichedTerms_flt <- merge(enrichedTerms_flt,RWHN_sig$data[,c("name", "seed")],
                               by.x = "Term", by.y = "name", all.x = T)
    enrichedTerms_flt$rwhn <- ifelse(enrichedTerms_flt$cluster == enrichedTerms_flt$seed, T, NA)
    enrichedTerms_flt <- dplyr::select(enrichedTerms_flt, -.data$seed) %>%  unique()
  }

  if(vis == T){
    gg <- ggplot2::ggplot(enrichedTerms_flt,
                          ggplot2::aes(y = .data$name,
                                       x = as.factor(.data$cluster)
                          )
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$Adjusted.P.value), color = "white") +
      ggplot2::scale_fill_gradient(breaks = seq(0, 0.05, 0.01), limits = c(0, 0.05),
                                   low = colours[1], high = colours[2]) +
      ggplot2::guides(color = FALSE,
                      fill = ggplot2::guide_colourbar(title="FDR", reverse = T)) +
      ggplot2::theme(legend.key.size = ggplot2::unit(.25, "cm"),
                     legend.title = ggplot2::element_text(size = 5),
                     legend.text = ggplot2::element_text(size = 5),
                     axis.text.y = ggplot2::element_text(size = 5),
                     axis.title = ggplot2::element_text(size = 6),
                     panel.background = ggplot2::element_rect(fill = "black"),
                     panel.grid = ggplot2::element_blank(),
                     legend.margin = ggplot2::margin(0,0,0,0, "cm")
      ) +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::xlab("Cluster") +
      ggplot2::ylab(database)

    return(gg)

    if(!is.null(RWHN_sig)){
      gg <- gg +
        ggplot2::geom_point(ggplot2::aes(shape = .data$rwhn), show.legend = F)
    }
  }else{
    return(enrichedTerms_flt)
  }

}
