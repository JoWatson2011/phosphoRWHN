# Internal Data generation

####
# STRING networks
####
#Read in STRING edges
STRING_flt <- function(path, conf){
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(biomaRt))

  STRING.expmt <- readr::read_table2(path) %>%
    dplyr::select(1:2,7) %>%
    filter(experimental >= conf) %>%
    tibble::rowid_to_column("ID")
  STRING.expmt$protein1 <- sub("9606.", "", STRING.expmt$protein1, fixed=T)
  STRING.expmt$protein2 <- sub("9606.", "", STRING.expmt$protein2, fixed=T)


  #ENSEMBL --> GENE name
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl" )
  genes1 <- unique(c(STRING.expmt$protein1, STRING.expmt$protein2))
  G_list <- getBM(filters = "ensembl_peptide_id",
                  attributes= c("ensembl_peptide_id", "hgnc_symbol"),
                  values=genes1,
                  mart= ensembl)

  #How many duplicates? How many missing gene names? How many not matched
  G_list.duplicates <- G_list[duplicated(G_list$ensembl_peptide_id),]
  sum(G_list$ensembl_peptide_id %in% G_list.duplicates[,1]) # each dupl. occurs twice (most likely)
  G_list.unfound <- G_list[G_list$hgnc_symbol == "",1]
  G_list.unincluded <- genes1[!(genes1 %in% G_list$ensembl_peptide_id)]
  table(G_list.unfound %in% G_list.unincluded)

  #match ensembl ids to gene symbols
  STRING.prot1 <- merge(STRING.expmt[,c(1:2,4)],G_list,by.x="protein1",by.y ="ensembl_peptide_id",all = T, sort = T)
  STRING.prot2 <- merge(STRING.expmt[,c(1,3:4)],G_list,by.x="protein2",by.y ="ensembl_peptide_id",all = T, sort = F)
  STRING.expmt.gene <- merge(STRING.prot1, STRING.prot2, by = "ID") %>%
    filter(!duplicated(ID))

  #verify that it worked...
  nrow(STRING.expmt) == nrow(STRING.expmt.gene)
  identical(STRING.expmt$experimental, STRING.expmt.gene$experimental.x)
  identical(STRING.expmt$experimental, STRING.expmt.gene$experimental.y)

  #remove ensembl ids
  STRING.expmt.gene <- STRING.expmt.gene[,c(1,4,7,6)]
  colnames(STRING.expmt.gene) <- c("ID", "protein1","protein2","experimental")

  return(STRING.expmt.gene)
}

path <- "" #PATH TO LOCAL DOWNLOAD OF STRING FILE
# Available from: https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens
# Version used in original paper:
# 9606.protein.links.detailed.v11.0.txt


STRING_lowconf <- STRING_flt("data/9606.protein.links.detailed.v11.0.txt", 400)
saveRDS(STRING_lowconf, "data/STRINGexpmtgene_lowconf.rds")
STRING_highconf <- STRING_flt("data/9606.protein.links.detailed.v11.0.txt", 700)
saveRDS(STRING_highconf, "data/STRINGexpmtgene_highconf.rds")


####
# GO Similarities
####
library("org.Hs.eg.db")
library("GO.db")
library("GOSemSim")

gene2GOi<-AnnotationDbi::select(org.Hs.eg.db,keys(org.Hs.egGO2EG),c("ENTREZID", "SYMBOL"),  "GOALL") %>% filter(ONTOLOGYALL == "BP")
gene2GO<-unstack(gene2GOi[,c(1,5)])

GOterms <- AnnotationDbi::select(GO.db, keys=gene2GOi$GOALL, columns=c("TERM"), keytype = "GOID")

GO2Gene<-unstack(stack(gene2GO)[2:1])

freq <- sapply(GO2Gene,length)
freqCutOff <- length(gene2GO)*0.05
highFreqTerms <- names(freq[freq>freqCutOff])


# ####
# # Pathway Similarities
# ####
# # Using the `all edges bma wang.txt` file
# # (found in the supplementary materials of:
# #
# # Stoney RA, Ames RM, Nenadic G, Robertson DL, Schwartz JM.
# # Disentangling the multigenic and pleiotropic nature of molecular function.
# # BMC Syst Biol. 2015;9 Suppl 6(Suppl 6):S3. doi:10.1186/1752-0509-9-S6-S3)
# #
# # and the following code:
#
# pathwaySimilarities <- readr::read_tsv("all edges bma wang.txt",
#                          col_names = c("sim", "func1", "func2")) %>%
#   mutate(func1 = tolower(func1),
#          func2 = tolower(func2)) %>%
#   mutate(func1 = gsub("_", " ", func1),
#          func2 = gsub("_", " ", func2)
#   )
#
usethis::use_data(STRING_highconf,
                  STRING_lowconf,
                  highFreqTerms,
                 # pathwaySimilarities,
                  internal = T,
                  overwrite = T)

