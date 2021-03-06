#'  Construct Heterogeneous network
#'
#' @param phosphoData Optional; data frame with IDs in format
#' *Gene name*_*AminoAcid+Position* (e.g. MAPK1_T185) in column 1,
#'  quantitative values in subsequent columns
#' @param clustering A named vector; names correspond to sites
#' in format *Gene name*_*AminoAcid+Position* (e.g. MAPK1_T185), values correspond to cluster
#' @param enrichrLib Library to extract functional annotations on.
#' See https://maayanlab.cloud/Enrichr/#stats
#' @param stringConf numeric; value between 0 and 1, definining confidence of STRING interactions for the formation of the protein layer
#' @param modules Calculate enrichment based on modules ? otherwise will be calculated on clusters
#' @param pval Integer; enrichment cut-off. Default 0.05
#' @param simplify Logical; if enrichR library is GO_Biological_Process_2018 or KEGG pathway, slim list of functional annotations?
#'
#' @importFrom tibble column_to_rownames
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter select mutate bind_rows
#' @importFrom stats cor na.omit
#' @importFrom tidyr pivot_longer separate separate_rows
#' @importFrom igraph graph_from_data_frame simplify as_data_frame
#' cluster_louvain membership
#' @importFrom httr POST content http_error status_code
#' @importFrom GOSemSim godata goSim
#' @importFrom purrr pluck
#' @importFrom AnnotationDbi select
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom GO.db GO.db GOBPCHILDREN
#' @importFrom rlang .data
#'
#' @examples
#' cl <- c("MAPK1_T185" = 1, "MAPK3_T202" = 1,
#'         "RAF1_S338" = 2, "RAF1_T491" = 2)
#' constructHetNet(cl)
#'
#' @return A list containing edgelists and vertices of heterogeneous network
#' @export constructHetNet

constructHetNet <- function(clustering, phosphoData = NULL,
                            enrichrLib = "GO_Biological_Process_2018",
                            stringConf = c(0.4,0.7),
                            modules = T, pval = 0.05, simplify = T){

  if(is.null(names(clustering))){
    stop("Vector must be named with gene names + modified residue
         (see documentation)")
  }

  ## Phospho
  phos <- lapply(unique(clustering)[order(unique(clustering))],
                 function(x){
                   cl <- names(clustering[clustering == x])

                   if(is.null(phosphoData)){
                     xpnd <- expand.grid(cl, cl)
                   }else{
                     if(colnames(phosphoData)[1] != "id" |
                        class(phosphoData$id) != "character"){
                       stop("phosphoData must have id
             (gene names + modified residue,
             as described in documentation)
             as first column")
                     }

                     xpnd <- phosphoData %>%
                       dplyr::filter(.data$id %in% cl) %>%
                       tibble::column_to_rownames(var = "id") %>%
                       t() %>%
                       stats::cor() %>%
                       as.data.frame() %>%
                       tibble::rownames_to_column(var = "phos1") %>%
                       tidyr::pivot_longer(-.data$phos1,
                                           values_to = "value",
                                           names_to = "phos2")

                     if(anyNA(xpnd)){
                       xpnd <- dplyr::select(xpnd, -.data$value)
                     } else {
                       xpnd <- xpnd %>%
                         dplyr::filter(.data$value < as.double(1.0),
                                       .data$value > 0.99) %>%
                         unique()
                     }
                   }
                   return(xpnd)
                 }) %>%
    dplyr::bind_rows() %>%
    igraph::graph_from_data_frame(directed = F) %>%
    igraph::simplify() %>%
    igraph::as_data_frame(what = "edges") %>%
    dplyr::select(phos1 = .data$from, phos2 = .data$to)

  ## Phospho >> Prot
  phos_prot <- igraph::graph_from_data_frame(phos, directed = F) %>%
    igraph::as_data_frame("vertices") %>%
    tidyr::separate(.data$name,
                    into = c("prot", "site"),
                    sep = "_",
                    remove = F) %>%
    dplyr::select(phos = .data$name , .data$prot)


  ## Prot
  prots <- unique(phos_prot$prot)
  prot <- if(stringConf[1] == 0.4){
    STRING_lowconf
  }else if(stringConf[1] == 0.7){
    STRING_highconf
  }else{
    #ERRORMSGHERE
    #Store error messages as intrnal data frame to be called??
  }

  prot <- prot %>%
    dplyr::filter(.data$protein1 %in% prots |
                    .data$protein2 %in% prots) %>%
    dplyr::filter(.data$experimental >= stringConf) %>%
    dplyr::select(.data$protein1, .data$protein2) %>%
    na.omit() %>%
    igraph::graph_from_data_frame(directed = F) %>%
    igraph::simplify(remove.multiple = F,remove.loops = T) %>%
    igraph::as_data_frame() %>%
    dplyr::select(prot1 = .data$from, prot2 = .data$to)

  ## Protein >> Func

  enrichedTerms <- #if(length(prots) > 25){
    lapply(unique(clustering)[order(unique(clustering))], function(i){
      ids <- names(clustering[clustering == i ] )
      cl_prots <- gsub("_.*", "", ids)

      enrichedTerms <- if(modules){
        el <- prot %>%
          dplyr::filter(.data$prot1 %in% cl_prots &
                          .data$prot2 %in% cl_prots) %>%
          na.omit()
        nw <- igraph::graph_from_data_frame(el,
                                            directed = F) %>%
          igraph::simplify(remove.multiple = F,
                           remove.loops = T) %>%
          igraph::cluster_louvain()
        cluster_nw <- data.frame(gene.symbol = names(igraph::membership(nw)),
                                 module = as.integer(igraph::membership(nw)),
                                 stringsAsFactors = F)

        cluster_nw <- lapply(unique(cluster_nw$module), function(x){
          sample <- cluster_nw %>%
            dplyr::filter(.data$module == x)

          enrichedTerms <- getEnrichr(sample$gene.symbol,
                                      enrichrLib = enrichrLib,
                                      pval = pval)

          if(nrow(enrichedTerms) == 0 ){
            enrichedTerms <- data.frame()
          }

          return(enrichedTerms)
        }) %>%
          dplyr::bind_rows()
      }else{
        enrichedTerms <- getEnrichr(cl_prots,
                                    enrichrLib = enrichrLib,
                                    pval = pval)
      }
    }) %>%
    dplyr::bind_rows()




  if(grepl("GO", enrichrLib) & !is.null(enrichedTerms)){
    enrichedTerms <- enrichedTerms %>%
      tidyr::separate(.data$Term,
                      into = c("Term", "GOID"),
                      sep = " \\(",
                      extra = "drop") %>%
      dplyr::mutate(GOID = sub("\\)",
                               "",
                               .data$GOID))
    if(simplify == T){

      ont <-  ifelse(grepl("Biological", enrichrLib),"BP",
                     ifelse(grepl("Molecular", enrichrLib),"MF",
                            ifelse(grepl("Cellular", enrichrLib), "CC", NA
                            )
                     )
      )

      hsGO <- suppressMessages(
        GOSemSim::godata('org.Hs.eg.db',
                         ont = ont,
                         computeIC=FALSE
        )
      )

      keepID <- simplifyGO(GOID = enrichedTerms$GOID,
                           semData = hsGO)
    }else{
      keepID <- enrichedTerms$GOID
    }

    prot_func <- enrichedTerms %>%
      dplyr::filter(.data$GOID %in% keepID) %>%
      tidyr::separate_rows(.data$Genes,
                           sep = ";") %>%
      dplyr::select(func = .data$Term,
                    prot = .data$Genes,
                    ID = .data$GOID)
  }else{
    enrichedTerms <- mutate(enrichedTerms,
                            Term = tolower(.data$Term))

    prot_func <- enrichedTerms %>%
      separate_rows(.data$Genes, sep = ";") %>%
      dplyr::select(func = .data$Term, prot = .data$Genes)
  }

  ## Func
  if(grepl("GO", enrichrLib)){

    ont <-  ifelse(grepl("Biological", enrichrLib),"BP",
                   ifelse(grepl("Molecular", enrichrLib),"MF",
                          ifelse(grepl("Cellular", enrichrLib), "CC", NA
                          )
                   )
    )

    if(!exists("hsGO")){
      hsGO <- suppressMessages(
        GOSemSim::godata('org.Hs.eg.db',
                         ont = ont,
                         computeIC=FALSE
        )
      )
    }
    semSim <- expand.grid(prot_func$ID,
                          prot_func$ID,
                          stringsAsFactors = F)

    semSim$sim <- apply(semSim, 1,
                        function(i) GOSemSim::goSim(i[1],
                                                    i[2],
                                                    semData = hsGO,
                                                    measure = "Wang")
    )

    func <- semSim %>%
      dplyr::filter(.data$sim < 1 & .data$sim > 0.7)
    func$func1 <- (AnnotationDbi::select(GO.db,
                                         keys = func$Var1, columns = "TERM",
                                         keytype = "GOID"))$TERM
    func$func2 <- (AnnotationDbi::select(GO.db,
                                         keys = func$Var2, columns = "TERM",
                                         keytype = "GOID"))$TERM
    func <- dplyr::select(func, .data$func1, .data$func2) %>%
      unique() %>%
      na.omit()
  }else{

    func <- pathwaySimilarities %>%
      dplyr::filter(.data$func1 %in% enrichedTerms$Term &
                      .data$func2 %in% enrichedTerms$Term) %>%
      dplyr::select(.data$func1, .data$func2)
  }


  v <- rbind(
    data.frame(v = unique(phos_prot$phos),
               layer = "phos",
               stringsAsFactors = F),
    data.frame(v = unique(phos$phos1),
               layer = "phos",
               stringsAsFactors = F),
    data.frame(v = unique(phos$phos2),
               layer = "phos",
               stringsAsFactors = F),
    data.frame(v = unique(prot_func$prot),
               layer = "prot",
               stringsAsFactors = F),
    data.frame(v = unique(phos_prot$prot),
               layer = "prot",
               stringsAsFactors = F),
    data.frame(v = unique(prot$prot1),
               layer = "prot",
               stringsAsFactors = F),
    data.frame(v = unique(prot$prot2),
               layer = "prot",
               stringsAsFactors = F),
    data.frame(v = unique(prot_func$func),
               layer = "func",
               stringsAsFactors = F)
  )

  if(length(unique(func$func1)) != 0){
    v <- rbind(
      v,
      data.frame(v = unique(func$func1),
                 layer = "func",
                 stringsAsFactors = F)
    )
  }
  if(length(unique(func$func2)) != 0){
    v <- rbind(
      v,
      data.frame(v = unique(func$func2),
                 layer = "func",
                 stringsAsFactors = F)
    )
  }

  v <- v %>%
    na.omit() %>%
    unique()

  edgelists <- list(x = phos, y = prot, z = func,
                    xy = phos_prot[,c("phos", "prot")],
                    yx = phos_prot[,c("prot", "phos")],
                    yz = prot_func[,c("prot", "func")],
                    zy = prot_func[,c("func", "prot")]
  ) %>%
    lapply(dplyr::rename, "to" = 2, "from" = 1)

  return(list(v = v,
              edgelists = edgelists)
  )
}
