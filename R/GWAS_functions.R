# A collection of functions for working with GWAS tables
# Niklas Schandry, Patrick Hüther 2019-2020


#' The araGenome object.
araGenome <- biomaRt::getBM(c("ensembl_gene_id",
                             "chromosome_name",
                             "start_position",
                             "end_position",
                             "strand",
                             "description",
                             "transcript_biotype"),
                           mart = biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org'))

#' araGenome as a nicely formatted GRanges object
araGenes <- araGenome %>%
  dplyr::rename(GeneId=ensembl_gene_id) %>%
  dplyr::mutate(strand=dplyr::case_when(strand == 1 ~ "+",
                                        strand == -1 ~ "-"),
                start = start_position,
                end = end_position) %>%
  dplyr::filter(stringr::str_detect(chromosome_name,"[1-9]")) %>%
  dplyr::mutate(seqnames=paste0("Chr",chromosome_name)) %>%
  plyranges::as_granges()

### Function Zone

#' Read an arbitrary GWAS table.
#' @details This function is a flexible version of \code{\link{read_gwas}}
#' @param gwas_path A gwas result table
#' @param p_col Name of the column containing p-values
#' @param chrom_col Name of the column containing chromosomes
#' @param pos_col Name of the column containing positions

format_gwas <- function(gwas_path, p_col = NULL, chrom_col = NULL, pos_col = NULL){
  gwas_object <- vroom::vroom(gwas_path)
  gwas_object <- dplyr::mutate(
                 pv = eval(p_col),
                 chrom = eval(chrom_col),
                 pos = eval(pos_col),
                 .data = gwas_object)
  fdr_corr <- qvalue::qvalue(p = gwas_object$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05 )]) %>% max() %>% log10() %>% abs()
  bf_corr <- (0.05/nrow(gwas_object)) %>% log10() %>% abs()
  gwas_object %<>%
    dplyr::mutate(chrom = as.double(stringr::str_extract(chrom, "[0-9]")),
                  Significant = dplyr::case_when(-log10(pv) > bf_corr ~ "Bonferroni",
                                                 -log10(pv) > fdr_thresh ~ "FDR",
                                                 TRUE ~ "Not"),
                  log10_p = -log10(pv),
                  fdr_thresh = fdr_thresh,
                  bf_corr = bf_corr,
                  model = dplyr::case_when(grepl(pattern = "_specific_", gwas_path) ~ "specific",
                                           grepl(pattern = "_common_", gwas_path) ~ "common",
                                           grepl(pattern = "_any_", gwas_path) ~ "any",
                                           TRUE ~ "Unknown"
                  )) %>%
    dplyr::arrange(dplyr::desc(log10_p))
  return(gwas_object)
}

#' Read a GWAS object produced by the sommer package.
#' @details This function is a flexible version of \code{\link{read_gwas}}
#' @param sommer_output output of the sommer::gwas() function.

format_sommer_gwas <- function(sommer_output) {
  sommer_scores <- sommer_output %$%
    scores %>%
    t %>%
    as.data.frame()
  score_col <- sommer_scores %>% dplyr::select(tidyselect::ends_with("score")) %>% colnames()
  gwas_table <- sommer_scores %>%
    tibble::rownames_to_column("Locus") %>%
    tidyr::separate("Locus", c("chrom","pos"), sep = "_" ) %>%
    dplyr::mutate(log10_p = .data[[paste(score_col)]],
                  pv = 10^(-(log10_p)),
                  chrom = as.numeric(chrom),
                  pos = as.numeric(pos),
                  mac = 5# pseudomac
    )
  fdr_corr <- qvalue::qvalue(p = gwas_table$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05)]) %>% max() %>% log10() %>% abs()
  bf_corr <- (0.05/nrow(gwas_table)) %>% log10() %>% abs()
  gwas_table %<>% dplyr::mutate(
    Significant = dplyr::case_when(-log10(pv) > bf_corr ~ "Bonferroni",
                                   -log10(pv) > fdr_thresh ~ "FDR",
                                   TRUE ~ "Not")) %>%
  dplyr::arrange(dplyr::desc(log10_p))
  return(gwas_table)
}

#' Read the GWAS tables produced by gwas-flow
#' @details A special case of the \code{\link{read_gwas}} function
#' @param gwasflow_out Results table (csv)
#' @param permutation_results Results from permutation test (if done)
#' @seealso \code{\link{read_gwas}}

format_gwasflow <- function(gwasflow_out, permutation_results = NULL) {
  gwas_table <- vroom::vroom(gwasflow_out) %>%
    dplyr::mutate(log10_p = -log10(pval),
           pv = pval,
           chrom = as.numeric(chr),
           pos = as.numeric(pos),
    )
  fdr_corr <- qvalue::qvalue(p = gwas_table$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05)]) %>% max() %>% log10() %>% abs()
  bf_corr <- (0.05/nrow(gwas_table)) %>% log10() %>% abs()
  if(!is.null(permutation_results)) {
    sig_perm <- vroom::vroom(permutation_results) %$% -log10(min_p) %>% min
    gwas_table %<>% dplyr::mutate(
      Significant = dplyr::case_when(-log10(pv) > sig_perm ~ "Permutation",
                                     -log10(pv) > bf_corr ~ "Bonferroni",
                                     -log10(pv) > fdr_thresh ~ "FDR",
                                     TRUE ~ "Not"))
    cat(paste("Permutation:", sig_perm ,"\nBF:", bf_corr,"\nFDR:", fdr_thresh))
  } else{
    gwas_table %<>% dplyr::mutate(
      Significant = dplyr::case_when(-log10(pv) > bf_corr ~ "Bonferroni",
                                     -log10(pv) > fdr_thresh ~ "FDR",
                                     TRUE ~ "Not"))
    cat(paste("\nBF:", bf_corr,"\nFDR:", fdr_thresh))
  }

  return(gwas_table)
}




#' Read a GWAS result file from limix (csv) and put it into a specific (somewhat arbitrary) format.
#' @details `read_gwas` returns a tibble containing SNP positions, -log10 transformed p-values and information whether a particular SNP passes a significance threshold after multiple testing correction (Bonferroni or Benjamini-Hochberg)
#' @param gwas_path Path to a gwas result file containing SNP positions and p-values
#' @seealso \code{\link{format_gwas}}

read_gwas <- function(gwas_path){
  gwas_object <- vroom::vroom(gwas_path)
  fdr_corr <- qvalue::qvalue(p = gwas_object$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05 )]) %>% max() %>% log10() %>% abs()
  bf_corr <- (0.05/nrow(gwas_object)) %>% log10() %>% abs()
  gwas_object %<>%
    dplyr::mutate(chrom = as.double(stringr::str_extract(chrom, "[0-9]")),
                          Significant = dplyr::case_when(-log10(pv) > bf_corr ~ "Bonferroni",
                                                         -log10(pv) > fdr_thresh ~ "FDR",
                                                         TRUE ~ "Not"),
                          log10_p = -log10(pv),
                          fdr_thresh = fdr_thresh,
                          bf_corr = bf_corr,
                          model = dplyr::case_when(grepl(pattern = "_specific_", gwas_path) ~ "specific",
                                                   grepl(pattern = "_common_", gwas_path) ~ "common",
                                                   grepl(pattern = "_any_", gwas_path) ~ "any",
                                                   TRUE ~ "Unknown"
                          )) %>%
    dplyr::arrange(dplyr::desc(log10_p))
  return(gwas_object)
}

#' Query the 1001genomes API to obtain accession ids that carry a SNP of interest
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}

get_polymorph_acc <- function(gwas_table, SNPrank){
  return(gwas_table %>%
           dplyr::arrange(dplyr::desc(log10_p)) %>%
           dplyr::slice(SNPrank) %>%
           dplyr::select(chrom, pos) %>%
           {httr::content(httr::GET(paste0("http://tools.1001genomes.org/api/v1.1/variants.csv?type=snps;accs=all;chr=",
                                           .$chrom,
                                           ";pos=",
                                           .$pos)), col_types = "iicccccc") }
  )
}

#' Plot an interactive map of accessions that carry a SNP of interest, needs a list of accessions with lng, lat fields
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_path The path to a table of accessions.
#' @seealso \code{\link{read_gwas}}

plot_acc_map <- function(gwas_table, SNPrank, acc_path = "~/labshare/lab/accessions/1001genomes-accessions.csv"){
  allAccessions <- readr::read_csv(acc_path)
  allAccessions %>%
    dplyr::filter(id %in% get_polymorph_acc(gwas_table, SNPrank)$strain) %>%
    leaflet::leaflet(data=.) %>%
      leaflet::addTiles() %>%
      leaflet::addCircleMarkers(lng=~longitude,
                                lat=~latitude,
                                label=~as.character(id),
                                popup=~as.character(name),
                                stroke=FALSE, radius=5, fillOpacity=0.8, color="#007243")
}



# get_expression <- function(GeneID){
#   if(is.null(GeneID)) {
#     stop("No GeneID supplied")
#   }
#   if(!exists("Kawakatsu_dat")){
#     Kawakatsu_dat <- readr::read_csv("https://arapheno.1001genomes.org/static/rnaseq/Epigenomic_Diversity_in_A._thaliana_(Kawakatsu_et_al._2016).csv") %>%
#       tidyr::as_tibble() %>%
#       magrittr::set_colnames(c("ACC_ID", colnames(.)[2:ncol(.)]))
#   }
#   if(!paste(GeneID) %in% colnames(Kawakatsu_dat)){
#     stop(paste0("No expression data in Kawakatsu2016 for GeneID: ",GeneID))
#   }
#     Kawakatsu_dat %>%
#       dplyr::select(ACC_ID, eval(GeneID)) %>%
#       tidyr::pivot_longer(dplyr::starts_with("AT"), names_to = "GeneID", values_to = "phenotype_value")
#
# }
# httr::content(httr::GET("https://arapheno.1001genomes.org/rest/rnaseq/52/ATMG01410/values/")) %>% {data.frame(matrix(unlist(.), nrow=length(.), byrow=T, stringsAsFactors=FALSE, ))}
#
# httr::content(httr::GET("https://arapheno.1001genomes.org/rest/rnaseq/52/ATMG01410/values/")) %>%
#   unlist %>%
#   plyr::ldply() %>%
#   pivot_wider(names_from =  .id, values_from = V1)

#' Get expression data for a gene of interest
#' @param GeneID An Arabidopsis thaliana gene identifier
#' @param study (default 52) the study of interest at arapheno, see list here
#' @param list_studies If TRUE, will return a list of available studies from arapheno.

get_expression <- function(GeneID = NULL, study = 52, list_studies = FALSE){
  if(list_studies){
    data.table::rbindlist(
      httr::content(httr::GET("https://arapheno.1001genomes.org/rest/study/list/"))
    ) %>% as.data.frame()
  }
  if(is.null(GeneID)) {
    stop("No GeneID supplied")
  }
  data.table::rbindlist(
    httr::content(httr::GET(paste0("https://arapheno.1001genomes.org/rest/rnaseq/",
                                   study,
                                   "/",
                                   GeneID,
                                   "/values")))
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(ACC_ID = accession_id)
 }


#' Based on a GeneID and a GWAS table and a rank, returns a table of expression values for that gene,
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param GeneID An Arabidopsis thaliana gene identifier
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{read_gwas}}

intersect_expression_snp <- function(GeneID, gwas_table, SNPrank){
  get_expression(GeneID = GeneID) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table,SNPrank)$strain ~ TRUE,
                                            TRUE ~ FALSE))
}

#' Based on a GWAS table and a rank, returns a table of expression values for the gene that is closest to that SNP
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{get_nearest_genes}}

retrieve_counts <- function(gwas_table, SNPrank){
  genes <- get_nearest_genes(gwas_table, SNPrank) %>%
             dplyr::slice(SNPrank) %>%
             .$GeneId
  get_expression(GeneID = genes) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                                             TRUE ~ FALSE))
}

#' Intersect SNP presence with Gene expression.s
#' @details Based on a GWAS table and a rank, returns a table of expression values for the gene that is closest to that SNP (or an arbitrary gene, see GeneID)where accessions that contain the SNP have TRUE in hasSNP
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param GeneID (optional) if not NULL, the counts for this gene will be plotted.
#' @param nobees Set to true to disable beeswarm geom
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{get_nearest_genes}}
#' @seealso \code{\link{retrieve_counts}}

plot_intersect_expression_snp <- function(gwas_table, SNPrank, GeneID = NULL, nobees = FALSE){
  if(nobees){
    overplot_geom <- geom_point(alpha = 0.3)
  } else{
    overplot_geom <- ggbeeswarm::geom_beeswarm(alpha = 0.3)
  }
  if(is.null(GeneID)){
  message("No GeneID supplied, plotting values for Gene closest to SNP")
  p <-  retrieve_counts(gwas_table, SNPrank) %>%
    ggplot(aes(x = hasSNP, y = phenotype_value)) +
    geom_boxplot(aes(fill = hasSNP)) +
    overplot_geom +
    labs(title = paste0("Expression of nearest gene by SNP presence"),
         caption = "Expression data from araPheno",
         x = "SNP present",
         y = "Value") +
    theme_bw() +
    facet_wrap(~phenotype_name)
  return(p)
  } else {
    get_expression(GeneID = GeneID) %>%
      dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table = gwas_table, SNPrank = SNPrank)$strain ~ TRUE,
                                              TRUE ~ FALSE)) %>%
      ggplot(aes(x = hasSNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = hasSNP)) +
      overplot_geom +
      labs(title = paste0("Expression of nearest gene by SNP presence"),
           caption = "Expression data from araPheno",
           x = "SNP present",
           y = "Value") +
      theme_bw() +
      facet_wrap(~phenotype_name)
  }
}



#' Subset a larger wide-format phenotype table to a defined phenotype.
#' @details This function takes a phenotype table in wide format, and returns that phenotype. Optionally includes a "specific" variable.
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.

get_phenotype <- function(phenotype_table, phenotype, acc_col = "ACC_ID", specific = NULL){
  if(!paste(phenotype) %in% colnames(phenotype_table)){
    stop(paste0("No data in the phenotype table for: ", phenotype))
  }
  if(acc_col != "ACC_ID"){
    message("Adding ACC_ID column")
    dplyr::mutate(ACC_ID = eval(acc_col), .data = phenotype_table)
  }
  if(is.null(specific)){
  phenotype_table %>%
    dplyr::select(ACC_ID, eval(phenotype)) %>%
    tidyr::pivot_longer(tidyselect::matches(eval(phenotype)), names_to = "Phenotype", values_to = "phenotype_value")
  } else {
  phenotype_table %>%
      dplyr::select(ACC_ID, eval(phenotype), eval(specific)) %>%
      tidyr::pivot_longer(tidyselect::matches(eval(phenotype)), names_to = "Phenotype", values_to = "phenotype_value")
  }

}

#' Get intersection between phenotype, and variant table.
#' @details  Based on a table of phenotypes, a phenotype name, a GWAS table and a rank, returns a table of phenotype values for that gene, where accessions that contain the SNP have TRUE in hasSNP
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_expression_snp}}


intersect_phenotype_snp <- function(phenotype_table, phenotype, gwas_table, SNPrank, acc_col = "ACC_ID", specific = NULL) {
  if(!is.null(specific)){
  get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col, specific = specific) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                              TRUE ~ FALSE))
  } else {
    get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col) %>%
      dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                                              TRUE ~ FALSE))
  }
}

#' Split phenotype table by SNP presence and plot
#' @details Based on a Phenotype table, the name of the phenotype and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @param nobees Set to true to disable beeswarm geom
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_phenotype_snp}}

plot_intersect_phenotype_snp <- function(phenotype_table, phenotype, gwas_table, SNPrank, acc_col = "ACC_ID", specific = NULL, nobees = FALSE){
  if(nobees){
    overplot_geom <- geom_point(alpha = 0.3)
  } else{
    overplot_geom <- ggbeeswarm::geom_beeswarm(alpha = 0.3)
  }
  if(is.null(specific)){
    p <-  intersect_phenotype_snp(phenotype_table = phenotype_table,
                                  phenotype = phenotype,
                                  gwas_table = gwas_table,
                                  SNPrank = SNPrank,
                                  specific = specific) %>%
      ggplot(aes(x = hasSNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = hasSNP)) +
      overplot_geom +
      labs(title = paste0("Phenotype values by SNP presence"),
           x = "SNP present",
           y = "Value") +
      theme_bw()
  } else{
    p <-  intersect_phenotype_snp(phenotype_table = phenotype_table,
                                  phenotype = phenotype,
                                  gwas_table = gwas_table,
                                  SNPrank = SNPrank,
                                  specific = specific) %>%
      ggplot(aes(x = hasSNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = hasSNP)) +
      overplot_geom +
      labs(title = paste0("Phenotype values by SNP presence"),
           x = "SNP present",
           y = "Value") +
      theme_bw() +
      facet_grid(reformulate( specific, "Phenotype")) #????
  }

  return(p)
}

#' Find genes that are closest to a SNP.
#' @details Based on a GWAS table, returns the gene annotation that is closest to each SNP for the number of SNP specified. This function always returns the closest annotation. To limit the lookup range use \code{\link{get_overlapping_genes}}. Lookup is done via ensembl plants; requires internet connection.
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param n_hit The number of SNPs that should be looked up.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_overlapping_genes}}



get_nearest_genes <- function(gwas_table = NULL, n_hit = 1){
  if(is.null(gwas_table)){
    stop("GWAS output missing")
  }
  # get annotation info from plants_mart
  if(!exists("araGenes")){

        # Get genome, rename columns, make into granges
    araGenes <- araGenome %>%
      dplyr::rename(GeneId=ensembl_gene_id) %>%
      dplyr::mutate(strand=dplyr::case_when(strand == 1 ~ "+",
                                            strand == -1 ~ "-"),
                    start = start_position,
                    end = end_position) %>%
      dplyr::filter(stringr::str_detect(chromosome_name,"[1-9]")) %>%
      dplyr::mutate(seqnames=paste0("Chr",chromosome_name)) %>%
      plyranges::as_granges()
  }


  snp <- gwas_table %>%
    dplyr::arrange(dplyr::desc(-log10(pv))) %>%
    dplyr::slice(1:n_hit) %>%
    dplyr::mutate(seqnames = paste0("Chr", chrom),
                  start = pos,
                  end = pos) %>%
    tibble::rownames_to_column("SNP_rank") %>%
    plyranges::as_granges()

  plyranges::join_nearest(snp, araGenes) %>%
    tidyr::as_tibble() %>%
    dplyr::select(SNP_rank,
                  chrom,
                  pos,
                  Significant,
                  log10_p,
                  mac,
                  GeneId,
                  description,
                  transcript_biotype,
                  start_position,
                  end_position#, .drop_ranges = TRUE
                  ) %>%
    dplyr::mutate(description = dplyr::case_when(description == "" ~ "N/A",
                                                 TRUE ~ description))

  #araGenes %>% filter_by_overlaps(snps, maxgap = distance)
}

#' Find genes that overlap with a SNP.
#' @details Based on a GWAS table, returns the gene annotation that is overlapping with each SNP for the number of SNP specified. If there is no overlap, the SNP is returned, but the GeneID column is empty. Lookup is done via ensembl plants; requires internet connection.
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param n_hit The number of SNPs that should be looked up.
#' @param distance The maximum distance from SNP to annotation, can be varied to lookup genes within a specific distance
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_nearest_genes}}


get_overlapping_genes <- function(gwas_table = NULL, n_hit = 1, distance = -1){
  if(is.null(gwas_table)){
    stop("GWAS output file missing")
  }
# get annotation info from plants_mart
  if(!exists("araGenes")){

    # Get genome, rename columns, make into granges
    araGenes <- araGenome %>%
      dplyr::rename(GeneId=ensembl_gene_id) %>%
      dplyr::mutate(strand=dplyr::case_when(strand == 1 ~ "+",
                                            strand == -1 ~ "-"),
                    start = start_position,
                    end = end_position) %>%
      dplyr::filter(stringr::str_detect(chromosome_name,"[1-9]")) %>%
      dplyr::mutate(seqnames=paste0("Chr",chromosome_name)) %>%
      plyranges::as_granges()
  }


  snp <- gwas_table %>%
    dplyr::arrange(dplyr::desc(-log10(pv))) %>%
    dplyr::slice(1:n_hit) %>%
    dplyr::mutate(seqnames = paste0("Chr", chrom),
                  start = pos,
                  end = pos,
                  max_distance = distance) %>%
    tibble::rownames_to_column("SNP_rank") %>%
    plyranges::as_granges()

  plyranges::join_overlap_inner(snp,araGenes, maxgap = distance) %>%
    plyranges::select(SNP_rank, chrom, pos, Significant, log10_p, mac, GeneId, max_distance, description, transcript_biotype, start_position, end_position, .drop_ranges = TRUE) %>%
    tidyr::as_tibble()
}

#' Based on a GWAS table, generates a manhatten plot.
#' @details For performance reasons, everything with a log10(p) smaller than p_filter is filtered out.
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param title Specify plot title
#' @param subtitle Specify plot subtitle
#' @param p_filter everything with a log10(p) below this value will not be included in the plot
#' @param mac_filter everything with a mac (minor allele count) below this will not be plotted. This requires a mac column.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{plot_annotated_gwas}}



plot_gwas <- function(gwas_table, title = "No Title", subtitle = NULL, p_filter = 2, mac_filter = 0 ){
  ggplot(aes(x=pos, y=log10_p), data = gwas_table %>% dplyr::filter(log10_p > p_filter) %>% dplyr::filter(mac > mac_filter)) +
    # geom_hline(linetype = "dotted", yintercept = bf_corr) +
    facet_grid(~chrom, scales = "free_x", switch = "x") +
    #  geom_label_repel(aes(x=pos, y= -log10(pv),label = gene), data = ft_specific_limix_mac5 %>% filter(gene != "NA")) +
    ggthemes::theme_few() +
    #  geom_mark_rect(aes(filter = gene != 'NA', label = gene)) +
    geom_point(aes(color = Significant)) +
    scale_color_viridis_d() +
    #  ylab("-log10(p)") +
    scale_y_continuous("-log10(pvalue)", breaks =  c(2,4,6,8), labels =  c("2","4","6","8")) +
    ggtitle(title, sub = subtitle) +
    #  guides(color = guide_legend("Adjusted Signifcance")) +
    theme(text = element_text(size=18), # Base text size
          plot.title = element_text(hjust = 0.5), #plot title, centered
          plot.subtitle = element_text(hjust = 0.5), # plot subtitle, centered
          axis.title.x = element_blank(),# axis tick labels.
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          panel.spacing.x = unit(0, "lines"), #remove spacing between facets
          legend.position = "bottom" #remove legend
    )
}

#' Plot annotated manhattan plot
#' @details Based on a GWAS table, generates a manhatten plot, with gene annotations. Number of annotations can be toggled by changing nlabels. By default only plots annotations for SNPs that are within a gene annotation. For performance reasons, everything with a log10(p) smaller than p_filter is filtered out.
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param title Specify plot title
#' @param subtitle Specify plot subtitle
#' @param nlabels How many labels should be plotted (default: 5)
#' @param labeltype Which entry from the lookup should be displayed (character, defaults to "GeneId"). Options: SNP_rank, chrom, pos, Significant, log10_p, mac, GeneId, max_distance, description, transcript_biotype, start_position, end_position.
#' @param match_nearest Defaults to FALSE; if TRUE will use find_nearest_genes instead of find_overlapping_genes for annotations
#' @param p_filter everything with a log10(p) below this value will not be included in the plot (default: 2)
#' @param mac_filter everything with a mac (minor allele count) below this will not be plotted.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{plot_gwas}}
#' @seealso \code{\link{get_nearest_genes}}
#' @seealso \code{\link{get_overlapping_genes}}

plot_annotated_gwas <- function(gwas_table,
                                title ="No Title",
                                subtitle = NULL,
                                nlabels = 5,
                                labeltype = "GeneId",
                                match_nearest = FALSE,
                                p_filter = 2,
                                mac_filter = 0) {
  if(!match_nearest){
    annotations <- gwas_table %>% get_overlapping_genes(nlabels) %>% tidyr::unite(labs, SNP_rank, labeltype, sep = " :")
  } else {
    annotations <- gwas_table %>% get_nearest_genes(nlabels) %>% tidyr::unite(labs, SNP_rank, labeltype, sep = " :")
  }
  #Step1: Bind those tables together
  gwas_table %>%
    dplyr::filter(abs(log10_p) > p_filter) %>%
    dplyr::filter(mac > mac_filter) %>%
    #Step 2: Plot
    ggplot(aes(x=pos, y=log10_p)) +
    geom_point(aes(color = Significant)) +
    ggrepel::geom_text_repel(aes(label = labs), data = annotations) +
    facet_grid(~chrom, scales = "free_x", switch = "x") +
    ggthemes::theme_few() +
    scale_color_viridis_d() +
    scale_y_continuous("-log10(pvalue)", breaks =  c(-2,-4,-6,-8, 2 ,4, 6, 8), labels =  c("2","4","6","8", "2","4","6","8")) +
    ggtitle(title, sub = subtitle) +
    #  guides(color = guide_legend("Adjusted Signifcance")) +
    theme(text = element_text(size=18), # Base text size
          plot.title = element_text(hjust = 0.5), #plot title, centered
          plot.subtitle = element_text(hjust = 0.5), # plot subtitle, centered
          axis.title.x = element_blank(),# axis tick labels.
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          panel.spacing.x = unit(0, "lines"), #remove spacing between facets
          legend.position = "none" #remove legend
    )
}


## Linkage Related Functions.
# Niklas Schandry

#' Calculate linkage around SNP of interest.
#' @details Details are fuzzy right now
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param rank Rank of the SNP of interest
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_depth Maximum SNP distance to calculate LD for (only relevant when anchored = FALSE)
#' @param ld_stats The LD statistics, see SNPStats::ld. Default is D.prime
#' @param ld_symmetric Should a symmetric matrix be returned? see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param use_all_acc If true, will override use_phenotype_table and use all 1135 sequenced accessions.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param anchored Perform anchored analysis. If TRUE Linkage will only be estimated for the SNP of interest, and the surrounding ones, but not between surrounding SNPs.

snp_linkage <- function(gwas_table,
                        rank,
                        nuc_range = 50000,
                        ld_depth = "1000",
                        ld_stats = "D.prime",
                        ld_symmetric = FALSE,
                        use_phenotype_table = NULL,
                        acc_col = "ACC_ID",
                        anchored = FALSE,
                        use_all_acc = FALSE) {

  gwas_table <- gwas_table %>%
    dplyr::arrange(dplyr::desc(log10_p))

  # Step 1: Download variant table
  ## Define region for API call

  region_lower <- gwas_table %>%
    dplyr::slice(rank) %>%
    {.$pos - (nuc_range / 2)}

  region_upper <- gwas_table %>%
    dplyr::slice(rank) %>%
    {.$pos + (nuc_range / 2)}

  if(region_lower < 1) {
    region_lower <- 1
    message("Nucleotide range out of bounds (negative start). Set start to 1.")
  }

  region <- gwas_table %>%
    dplyr::slice(rank) %>%
    {paste0("Chr", .$chrom, ":", region_lower, ".." , region_upper)}

  ## Define genotypes for API call
  if(is.null(use_phenotype_table)){
    message("Calculating LD based on strains that carry SNP")
    genotypes <- stringr::str_flatten(get_polymorph_acc(gwas_table = gwas_table, SNPrank = rank)$strain, collapse = ",")
  }
  if(!is.null(use_phenotype_table)){
    message("Calculating LD based on all strains in phenotype table")
    genotypes <- stringr::str_flatten(levels(as.factor(unlist(use_phenotype_table[, eval(acc_col)]))), collapse = ",")
    ## construct call
  }
  if(use_all_acc){
    if(!is.null(use_phenotype_table)){
      message("NOT using supplied genotype table because use_all_acc is TRUE. \n
            Calculating LD based on all sequenced strains (1135 genomes)")}
    if(is.null(use_phenotype_table)){
      message("Calculating LD based on all sequenced strains (1135 genomes)")}
    genotypes <- c("88,108,139,159,265,350,351,403,410,424,428,430,
    470,476,484,504,506,531,544,546,628,630,680,681,685,687,728,
    742,763,765,766,768,772,801,853,854,867,868,870,915,932,991,
    992,997,1002,1006,1061,1062,1063,1066,1070,1158,1166,1254,
    1257,1313,1317,1552,1612,1622,1651,1652,1676,1684,1739,1741,
    1756,1757,1793,1797,1819,1820,1829,1834,1835,1851,1852,1853,
    1872,1890,1925,1942,1943,1954,2016,2017,2031,2053,2057,2081,
    2091,2106,2108,2141,2159,2166,2171,2191,2202,2212,2239,2240,
    2276,2278,2285,2286,2317,2370,2412,4779,4807,4826,4840,4857,
    4884,4900,4939,4958,5023,5104,5151,5165,5210,5236,5249,5253,
    5276,5279,5349,5353,5395,5486,5577,5644,5651,5717,5718,5720,
    5726,5741,5748,5757,5768,5772,5776,5779,5784,5798,5800,5811,
    5822,5830,5831,5832,5836,5837,5856,5860,5865,5867,5874,5890,
    5893,5907,5921,5950,5984,5993,6008,6009,6010,6011,6012,6013,
    6016,6017,6019,6020,6021,6022,6023,6024,6025,6030,6034,6035,
    6036,6038,6039,6040,6041,6042,6043,6046,6064,6069,6070,6071,
    6073,6074,6076,6077,6085,6086,6087,6088,6090,6091,6092,6094,
    6095,6096,6097,6098,6099,6100,6101,6102,6104,6105,6106,6107,
    6108,6109,6111,6112,6113,6114,6115,6118,6119,6122,6123,6124,
    6125,6126,6128,6131,6132,6133,6134,6136,6137,6138,6140,6141,
    6142,6145,6148,6149,6150,6151,6153,6154,6163,6166,6169,6172,
    6173,6174,6177,6180,6184,6188,6189,6191,6192,6193,6194,6195,
    6198,6201,6202,6203,6209,6210,6214,6216,6217,6218,6220,6221,
    6231,6235,6237,6238,6240,6241,6242,6243,6244,6252,6255,6258,
    6268,6276,6284,6296,6390,6396,6413,6424,6434,6445,6680,6739,
    6740,6744,6749,6750,6805,6806,6814,6830,6897,6898,6900,6901,
    6903,6904,6907,6908,6909,6911,6913,6915,6917,6918,6919,6920,
    6922,6923,6924,6926,6927,6929,6931,6932,6933,6938,6940,6943,
    6944,6945,6951,6956,6957,6958,6959,6960,6961,6963,6966,6967,
    6968,6969,6970,6971,6973,6974,6975,6976,6979,6981,6982,6984,
    6986,6987,6989,6990,6992,6997,7000,7002,7003,7008,7013,7014,
    7025,7026,7028,7031,7033,7036,7058,7061,7062,7063,7064,7067,
    7068,7071,7072,7077,7081,7092,7094,7096,7102,7103,7106,7107,
    7109,7111,7117,7119,7120,7125,7126,7127,7130,7133,7138,7143,
    7147,7158,7160,7161,7162,7163,7164,7165,7169,7177,7181,7183,
    7186,7192,7199,7202,7203,7207,7208,7209,7213,7217,7218,7223,
    7231,7236,7244,7248,7250,7255,7258,7268,7273,7276,7280,7282,
    7287,7288,7296,7298,7305,7306,7307,7314,7316,7319,7320,7322,
    7323,7327,7328,7332,7333,7337,7342,7343,7344,7346,7347,7349,
    7350,7353,7354,7356,7358,7359,7372,7373,7377,7378,7382,7383,
    7384,7387,7394,7396,7404,7411,7413,7415,7416,7417,7418,7419,
    7424,7427,7430,7460,7461,7471,7475,7477,7514,7515,7516,7517,
    7520,7521,7523,7525,7529,7530,7566,7568,7717,7757,7767,7917,
    7947,8037,8057,8077,8132,8171,8214,8222,8227,8230,8231,8233,
    8234,8235,8236,8237,8238,8239,8240,8241,8242,8243,8244,8246,
    8247,8249,8256,8258,8259,8264,8283,8284,8285,8290,8297,8306,
    8307,8311,8312,8326,8334,8335,8337,8343,8351,8354,8357,8365,
    8366,8369,8376,8386,8387,8419,8420,8422,8424,8426,8427,8464,
    8483,8699,8723,9027,9057,9058,9067,9069,9070,9075,9078,9079,
    9081,9084,9085,9089,9091,9095,9099,9100,9102,9103,9104,9105,
    9106,9111,9113,9114,9115,9121,9125,9128,9130,9131,9133,9134,
    9298,9312,9314,9321,9323,9332,9336,9339,9343,9352,9353,9363,
    9369,9370,9371,9380,9381,9382,9383,9386,9388,9390,9391,9392,
    9394,9395,9399,9402,9404,9405,9407,9408,9409,9412,9413,9416,
    9421,9427,9433,9436,9437,9442,9450,9451,9452,9453,9454,9455,
    9470,9471,9476,9481,9503,9506,9507,9508,9509,9510,9511,9512,
    9513,9514,9515,9517,9518,9519,9520,9521,9522,9523,9524,9525,
    9526,9527,9528,9529,9530,9531,9532,9533,9534,9535,9536,9537,
    9539,9540,9541,9542,9543,9544,9545,9546,9547,9548,9549,9550,
    9551,9552,9553,9554,9555,9556,9557,9558,9559,9560,9561,9562,
    9564,9565,9567,9568,9569,9571,9573,9574,9576,9577,9578,9579,
    9581,9582,9583,9584,9585,9586,9587,9588,9589,9590,9591,9592,
    9593,9594,9595,9596,9597,9598,9599,9600,9601,9602,9606,9607,
    9608,9609,9610,9611,9612,9613,9615,9616,9617,9619,9620,9621,
    9622,9624,9625,9626,9627,9628,9629,9630,9631,9632,9633,9634,
    9635,9636,9637,9638,9639,9640,9641,9642,9643,9644,9645,9646,
    9647,9648,9649,9651,9653,9655,9656,9657,9658,9659,9660,9661,
    9663,9664,9665,9666,9667,9668,9669,9670,9671,9672,9673,9676,
    9677,9678,9679,9680,9681,9682,9683,9684,9685,9686,9687,9689,
    9690,9691,9692,9693,9694,9695,9696,9697,9698,9699,9700,9701,
    9703,9704,9705,9706,9707,9708,9709,9710,9711,9712,9713,9714,
    9716,9717,9718,9719,9720,9721,9722,9723,9725,9726,9727,9728,
    9729,9730,9731,9732,9733,9735,9736,9737,9738,9739,9741,9743,
    9744,9745,9747,9748,9749,9754,9755,9756,9757,9758,9759,9761,
    9762,9764,9766,9768,9769,9770,9771,9772,9774,9775,9776,9777,
    9778,9779,9780,9781,9782,9783,9784,9785,9786,9787,9788,9789,
    9790,9791,9792,9793,9794,9795,9796,9797,9798,9799,9800,9801,
    9802,9803,9804,9805,9806,9807,9808,9809,9810,9811,9812,9813,
    9814,9815,9816,9817,9819,9820,9821,9822,9823,9824,9825,9826,
    9827,9828,9830,9831,9832,9833,9834,9835,9836,9837,9838,9839,
    9840,9841,9843,9844,9845,9846,9847,9848,9849,9850,9851,9852,
    9853,9854,9855,9856,9857,9858,9859,9860,9861,9862,9864,9866,
    9867,9868,9869,9870,9871,9873,9874,9875,9876,9877,9878,9879,
    9880,9881,9882,9883,9885,9886,9887,9888,9890,9891,9892,9894,
    9895,9897,9898,9899,9900,9901,9902,9903,9904,9905,9906,9908,
    9909,9910,9911,9912,9914,9915,9917,9918,9920,9921,9924,9925,
    9926,9927,9928,9929,9930,9932,9933,9935,9937,9938,9939,9941,
    9942,9943,9944,9945,9946,9947,9948,9949,9950,9951,9952,9953,
    9955,9956,9957,9958,9959,9960,9962,9963,9964,9965,9966,9968,
    9969,9970,9971,9972,9973,9974,9975,9976,9978,9979,9980,9981,
    9982,9983,9984,9985,9986,9987,9988,9990,9991,9993,9995,9996,
    9997,9998,9999,10001,10002,10004,10005,10006,10008,10009,
    10010,10011,10012,10013,10014,10015,10017,10018,10020,10022,
    10023,10027,14312,14313,14314,14315,14318,14319,15560,15591,
    15592,15593,18694,18696,19949,19950,19951") # Very elegant.
  }
  subset_url <- paste0("http://tools.1001genomes.org/api/v1/vcfsubset/strains/",
                       stringr::str_remove_all(genotypes, "[\n| ]"),
                       "/regions/",
                       region,
                       "/type/fullgenome/format/vcf")

  ## Make tempfile
  tmp <- tempfile(fileext = ".vcf")

  ## Write vcf table to tempfile

  writeLines( httr::content( httr::GET(subset_url ) ), tmp)

  message(paste0("Downloaded genotype vcf from ",subset_url," to: ", tmp,"\n"))
  # Step 2 Read vcf table from tempfile

  tmp_vcf <- VariantAnnotation::readVcf(tmp)
  tmp_vcf <- tmp_vcf[VariantAnnotation::isSNV(tmp_vcf)]
  tmp_SM <- VariantAnnotation::genotypeToSnpMatrix(tmp_vcf)

  # Step 3 convert to SNPMatrix

  SM_for_linkage <- tmp_SM$genotypes

  # Step 4 calculate ld on genotypes table of SM; return
  if(!anchored){
    ld_tab <- snpStats::ld(SM_for_linkage, depth = ld_depth, stats = ld_stats, symmetric = ld_symmetric )
  }

  # The anchored approach: Calculate LD not for all vs all in a region, but for a region against one specific SNP.

  if(anchored){
      region <- gwas_table %>%
      dplyr::slice(rank) %>%
      {paste0("Chr", .$chrom, ":", .$pos , ".." , .$pos)}
    message(paste0("Anchored analysis, with ", region,  " as anchor.\n"))
    # Start with getting the anchor information, similar to above, but the region is only one position.
    anchor_url <-  paste0("http://tools.1001genomes.org/api/v1/vcfsubset/strains/",
                          genotypes,
                          "/regions/",
                          region,
                          "/type/fullgenome/format/vcf")

    anchor_tmp <- tempfile(fileext = ".vcf")

    writeLines( httr::content( httr::GET(anchor_url ) ), anchor_tmp)

    message(paste0("Downloaded anchor vcf from", anchor_url ," to ", anchor_tmp,"\n"))

    anchor_vcf <- VariantAnnotation::readVcf(anchor_tmp)

    #anchor_vcf <- anchor_vcf[VariantAnnotation::isSNV(anchor_vcf)]

    anchor_SM <- VariantAnnotation::genotypeToSnpMatrix(anchor_vcf)

    anchor_SM <- anchor_SM$genotypes

    ld_tab <- snpStats::ld(x = anchor_SM , y = SM_for_linkage, stats = ld_stats)
  }
  return(ld_tab)
}

#' Calculate linkage around SNP of interest.
#' @details Details are fuzzy right now
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param rank Rank of the SNP of interest
#' @param snpmatrix_path the path to the SNPmatrix in fst format.
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_depth Maximum SNP distance to calculate LD for (only relevant when anchored = FALSE)
#' @param ld_stats The LD statistics, see SNPStats::ld. Default is D.prime
#' @param ld_symmetric Should a symmetric matrix be returned? see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions in the SNPmatrix are used. Accessions are assumed to be named as numbers.
#' @param acc_col the column name of the accession column in the phenotype table
#' @param anchored Perform anchored analysis. If TRUE Linkage will only be estimated for the SNP of interest, and the surrounding ones, but not between surrounding SNPs.

snp_linkage_snpmatrix <- function(gwas_table,
                        rank,
                        snpmatrix_path = NULL,
                        nuc_range = 50000,
                        ld_depth = "1000",
                        ld_stats = "D.prime",
                        ld_symmetric = FALSE,
                        use_phenotype_table = NULL,
                        acc_col = "ACC_ID",
                        anchored = FALSE) {

  if(!file.exists(snpmatrix_path) | !(tools::file_ext(snpmatrix_path) == "fst" )) {
    stop("Please point to snpmatrix in fst format")
  }

  gwas_table <- gwas_table %>%
    dplyr::arrange(dplyr::desc(log10_p))

  # Step 1: Download variant table

  ## Define region for API call

  region_lower <- gwas_table %>%
    dplyr::slice(rank) %>%
    {.$pos - (nuc_range / 2)}

  region_upper <- gwas_table %>%
    dplyr::slice(rank) %>%
    {.$pos + (nuc_range / 2)}

  chrom <- gwas_table %>%
    dplyr::slice(rank) %>%
    dplyr::select(chrom) %>%
    unique() %>%
    paste()

  if(region_lower < 1) {
    region_lower <- 1
    message("Nucleotide range out of bounds (negative start). Set start to 1.")
  }

  ## Define genotypes for subsetting
  if(is.null(use_phenotype_table)){
    message("Calculating LD based on strains in SNPmatrix. Strains are assumed to be numeric.")
    genotypes <- colnames(fst::fst(snpmatrix_path)) %>% {suppressWarnings(as.numeric(.))} %>% na.omit()
  }

  if(!is.null(use_phenotype_table)){
    message("Calculating LD based on all strains in phenotype table,")
    genotypes <- levels(as.factor(unlist(use_phenotype_table[, eval(acc_col)])))
    ## construct call
  }

  pos = c(region_lower:region_upper)

  # Very base R subsetting
  tmpmatrix <- fst::fst(snpmatrix_path)
  tmpmatrix <- tmpmatrix[(tmpmatrix$chrom == chrom) & (tmpmatrix$pos %in% pos), paste(c(genotypes, "chrom", "pos"))]
  rownames(tmpmatrix) <- paste(tmpmatrix$chrom, tmpmatrix$pos, sep = "_")
  tmpmatrix <- tmpmatrix[, paste(genotypes)]
  tmpmatrix <- tmpmatrix %>%
    as.matrix() %>%
    t()
  storage.mode(tmpmatrix) <- "raw"

  SM_for_linkage <- methods::new("SnpMatrix", tmpmatrix)

  # Step 4 calculate ld on genotypes table of SM; return
  if(!anchored){
    ld_tab <- snpStats::ld(SM_for_linkage, depth = ld_depth, stats = ld_stats, symmetric = ld_symmetric )
  }

  # The anchored approach: Calculate LD not for all vs all in a region, but for a region against one specific SNP.

  if(anchored){
    region <- gwas_table %>%
      dplyr::slice(rank) %>%
      {paste0("Chr", .$chrom, ":", .$pos , ".." , .$pos)}
    message(paste0("Anchored analysis, with ", region,  " as anchor.\n"))
    # Start with getting the anchor information, similar to above, but the region is only one position.
    anchor_pos <- gwas_table %>%
      dplyr::slice(rank) %>%
      .$pos
    anchor_tmpmatrix <- fst::fst(snpmatrix_path)
    anchor_tmpmatrix <- anchor_tmpmatrix[(anchor_tmpmatrix$chrom == chrom) & (anchor_tmpmatrix$pos == pos), paste(c(genotypes, "chrom", "pos"))]
    rownames(anchor_tmpmatrix) <- paste(anchor_tmpmatrix$chrom, anchor_tmpmatrix$pos, sep = "_")
    anchor_tmpmatrix <- anchor_tmpmatrix[, paste(genotypes)]
    anchor_tmpmatrix <- anchor_tmpmatrix %>%
      as.matrix() %>%
      t()
    storage.mode(anchor_tmpmatrix) <- "raw"

    anchor_SM <-  methods::new("SnpMatrix", anchor_tmpmatrix)

    ld_tab <- snpStats::ld(x = anchor_SM , y = SM_for_linkage, stats = ld_stats)
  }
  return(ld_tab)
}

#' Produce a plot of linked SNPs
#' @details This function takes a gwas table and a rank and produces a plot that illustrates, SNPs that are within nuc_range around the SNP, the degree of linkage, the impact of the other SNPs and a track of gene annotations to easily identify linked SNPs that have a high impact.
#' @details Impact according to 1001genomes.org
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param rank Rank of the SNP of interest
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_stats The LD statistics, see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param use_all_acc If true, will override use_phenotype_table and use all 1135 sequenced accessions.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param linkage_cutoff Only SNPs that have an LD values >= this will be plotted (default 0.8)
#' @param LD_legend Toggle color legend for LD. Off by default.
#' @param data_only If true does not plot, only returns the data used for plotting.

plot_anchored_ld <-   function(gwas_table,
                               rank,
                               nuc_range = 50000,
                               ld_stats = c("D.prime") ,
                               use_phenotype_table = NULL,
                               use_all_acc = FALSE,
                               acc_col = "ACC_ID",
                               linkage_cutoff = 0.8,
                               LD_legend = FALSE,
                               data_only = FALSE) {
  if(length(ld_stats) != 1) {
    stop("The number of LD stats must be exactly one.")
  }

  # Calculate anchored LD

  anc_ld <- snp_linkage(gwas_table = gwas_table,
                        rank = rank,
                        nuc_range = nuc_range,
                        anchored = TRUE,
                        use_phenotype_table = use_phenotype_table,
                        ld_stats = ld_stats)

  # Check if this actually produced results

  if(anc_ld %>% na.omit() %>% nrow() == 0) {
    stop("LD calculation returned only NAs.")
  }

  # Define SNP details.

  snp_pos <- gwas_table %>%
    dplyr::slice(rank) %>% .$pos

  chrom = gwas_table %>%
    dplyr::slice(rank) %>%
    .$chrom

  start_pos = gwas_table %>%
    dplyr::slice(rank) %>%
    {.$pos - as.numeric(nuc_range) / 2}

  if(start_pos < 1) {
    start_pos <- 1 # No message, because snp_linkage will issue a message.
  }

  end_pos = gwas_table %>%
    dplyr::slice(rank) %>% {.$pos + as.numeric(nuc_range) / 2}

  ## Define genotypes for API call
  if(is.null(use_phenotype_table)){
    message("Calculating LD based on strains that carry SNP. This could be a bad idea!")
    genotypes <- stringr::str_flatten(get_polymorph_acc(gwas_table = gwas_table, SNPrank = rank)$strain, collapse = ",")
  }
  if(!is.null(use_phenotype_table)){
    message("Calculating LD based on all strains in phenotype table")
    genotypes <- stringr::str_flatten(levels(as.factor(unlist(use_phenotype_table[, eval(acc_col)]))), collapse = ",")

  }
  ## construct call (this directly reads the csv from 1001genomes)
  snp_impacts <- httr::content(
    httr::GET(
      paste0(
        "http://tools.1001genomes.org/api/v1.1/effects.csv?accs=", genotypes,
        ";chr=", chrom,
        ";start=", start_pos,
        ";end=", end_pos,
        ";type=snps")
    ), col_types = "iifccccccccccc")
  ## Labels for the gene plots

  gene_labels <- araGenome %>%
    dplyr::filter(chromosome_name == chrom, start_position %in% c(start_pos:end_pos), end_position %in% c(start_pos:end_pos)) %>%
    dplyr::mutate(gene = ensembl_gene_id,
                  molecule = 0,
                  start = start_position,
                  end = end_position,
                  direction = strand,
                  strand = dplyr::case_when(strand == -1 ~ "reverse", TRUE ~ "forward"),
                  layer_var = "ORF")

  ## Transforming LD table into more plotable table.

  plot_data <- anc_ld %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(tidyselect::matches(":")) %>%
    na.omit %>%
    dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(name, "[:|_]",3)[,2])) %>%
    dplyr::left_join(., snp_impacts, by = c("pos")) %>%
    dplyr::mutate(layer_var = effect_impact)

  ## Build Plot

  p <- ggplot(data = plot_data) +
    # Tile Geom for LD values
    geom_tile(aes(x = pos , color = value, y = 0, height = 2),
              data = plot_data %>% dplyr::mutate(layer_var = "LD") %>% dplyr::filter(layer_var == "LD")) +
    # Line denoting SNP of interest
    geom_vline(aes(xintercept = snp_pos), color = "darkred") +
    # Modifier SNPs
    geom_point(aes(x= pos, y = 0, shape = effect_impact), size = 2.5, alpha = 0.8,
               data = plot_data %>% dplyr::filter(value >= linkage_cutoff)) +

    # color scale
    scale_color_viridis_c(option = "plasma", direction = -1) +
    # Gene Arrows
    gggenes::geom_gene_arrow(aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = direction), data = gene_labels) +
    # Gene labels
    gggenes::geom_gene_label(aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene), data = gene_labels) + # Gene labels
    # Minimal Theme
    theme_minimal() +
    # ylim(-0.1,0.1) +
    facet_grid(rows = vars(layer_var %>% forcats::fct_relevel("HIGH", "MODERATE", "LOW", "MODIFIER", "LD", "ORF")), switch = "y") +
    # Theme adjusments, mainly removing the y-axis
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          plot.background = element_rect(fill = "white"),
          legend.position = "bottom",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  if(isFALSE(LD_legend)){
    p <- p +
      guides(color = FALSE,
             shape = FALSE,
             fill = guide_legend("Gene ID")) # Color guide
  } else
    p <- p +
    guides(color = guide_colorbar(title = "Linkage"),
           shape = FALSE,
           fill = guide_legend("Gene ID")) # Color guide

  if(data_only){
    return(plot_data)
  } else{
    return(p)
  }

}

#' Produce a plot of linked SNPs based on a SNPMatrix
#' @details This function takes a gwas table and a rank and produces a plot that illustrates, SNPs that are within nuc_range around the SNP, the degree of linkage, the impact of the other SNPs and a track of gene annotations to easily identify linked SNPs that have a high impact.
#' @details Impact according to 1001genomes.org
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param rank Rank of the SNP of interest
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param snpmatrix_path filepath to SNPmatrix in fst format.
#' @param ld_stats The LD statistics, see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param linkage_cutoff Only SNPs that have an LD values >= this will be plotted (default 0.8)
#' @param LD_legend Toggle color legend for LD. Off by default.
#' @param data_only If true does not plot, only returns the data used for plotting.

plot_anchored_ld_snpmatrix <-   function(gwas_table,
                               rank,
                               nuc_range = 50000,
                               snpmatrix_path = NULL,
                               ld_stats = c("D.prime") ,
                               use_phenotype_table = NULL,
                               acc_col = "ACC_ID",
                               linkage_cutoff = 0.8,
                               LD_legend = FALSE,
                               data_only = FALSE) {

  if(length(ld_stats) != 1) {
    stop("The number of LD stats must be exactly one.")
  }

  # Calculate anchored LD

  anc_ld <- snp_linkage_snpmatrix(gwas_table = gwas_table,
                        rank = rank,
                        nuc_range = nuc_range,
                        snpmatrix_path = snpmatrix_path,
                        anchored = TRUE,
                        use_phenotype_table = use_phenotype_table,
                        acc_col = acc_col,
                        ld_stats = ld_stats
                        )

  # Check if this actually produced results

  if(anc_ld %>% na.omit() %>% nrow() == 0) {
    stop("LD calculation returned only NAs.")
  }

  # Define SNP details.

  snp_pos <- gwas_table %>%
    dplyr::slice(rank) %>% .$pos

  chrom = gwas_table %>%
    dplyr::slice(rank) %>%
    .$chrom

  start_pos = gwas_table %>%
    dplyr::slice(rank) %>% {.$pos - as.numeric(nuc_range) / 2}
  if(start_pos < 1) {
    start_pos <- 1 # No message, because snp_linkage will issue a message.
  }
  end_pos = gwas_table %>%
    dplyr::slice(rank) %>% {.$pos + as.numeric(nuc_range) / 2}

  ## Define genotypes

  if(is.null(use_phenotype_table)){
    message("Calculating LD based on strains in SNPmatrix. Strains are assumed to be numeric.")
    genotypes <- colnames(fst::fst(snpmatrix_path)) %>% {suppressWarnings(as.numeric(.))} %>% na.omit() %>% stringr::str_flatten(.,",")
  }

  if(!is.null(use_phenotype_table)){
    message("Calculating LD based on all strains in phenotype table,")
    genotypes <- stringr::str_flatten(levels(as.factor(unlist(use_phenotype_table[, eval(acc_col)]))), collapse = ",")
    ## construct call
  }
  ## construct call (this directly reads the csv from 1001genomes)
  snp_impacts <- httr::content(
    httr::GET(
      paste0(
        "http://tools.1001genomes.org/api/v1.1/effects.csv?accs=", genotypes,
        ";chr=", chrom,
        ";start=", start_pos,
        ";end=", end_pos,
        ";type=snps")
    ), col_types = "iifccccccccccc")
  ## Labels for the gene plots

  gene_labels <- araGenome %>%
    dplyr::filter(chromosome_name == chrom, start_position %in% c(start_pos:end_pos), end_position %in% c(start_pos:end_pos)) %>%
    dplyr::mutate(gene = ensembl_gene_id,
                  molecule = 0,
                  start = start_position,
                  end = end_position,
                  direction = strand,
                  strand = dplyr::case_when(strand == -1 ~ "reverse", TRUE ~ "forward"),
                  layer_var = "ORF")

  ## Transforming LD table into more plotable table.

  plot_data <- anc_ld %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(tidyselect::matches(":")) %>%
    na.omit %>%
    dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(name, "[:|_]",3)[,2])) %>%
    dplyr::left_join(., snp_impacts, by = c("pos")) %>%
    dplyr::mutate(layer_var = effect_impact)

  ## Build Plot

  p <- ggplot(data = plot_data) +
    # Tile Geom for LD values
    geom_tile(aes(x = pos , color = value, y = 0, height = 2),
              data = plot_data %>% dplyr::mutate(layer_var = "LD") %>% dplyr::filter(layer_var == "LD")) +
    # Line denoting SNP of interest
    geom_vline(aes(xintercept = snp_pos), color = "darkred") +
    # Modifier SNPs
    geom_point(aes(x= pos, y = 0, shape = effect_impact), size = 2.5, alpha = 0.8,
               data = plot_data %>% dplyr::filter(value >= linkage_cutoff)) +

    # color scale
    scale_color_viridis_c(option = "plasma", direction = -1) +
    # Gene Arrows
    gggenes::geom_gene_arrow(aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = direction), data = gene_labels) +
    # Gene labels
    gggenes::geom_gene_label(aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene), data = gene_labels) + # Gene labels
    # Minimal Theme
    theme_minimal() +
    # ylim(-0.1,0.1) +
    facet_grid(rows = vars(layer_var %>% forcats::fct_relevel("HIGH", "MODERATE", "LOW", "MODIFIER", "LD", "ORF")), switch = "y") +
    # Theme adjusments, mainly removing the y-axis
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing.y = unit(0, "lines"),
          plot.background = element_rect(fill = "white"),
          legend.position = "bottom",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())
  if(isFALSE(LD_legend)){
    p <- p +
      guides(color = FALSE,
             shape = FALSE,
             fill = guide_legend("Gene ID")) # Color guide
  } else
    p <- p +
    guides(color = guide_colorbar(title = "Linkage"),
           shape = FALSE,
           fill = guide_legend("Gene ID")) # Color guide

  if(data_only){
    return(plot_data)
  } else{
    return(p)
  }

}

#' Obtain accession ids that carry a SNP of interest from a local SNPmatrix. This is a sister of get_polymorph_acc
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param SNPmatrix The SNPmatrix to use for identifying accessions that carry the relevant SNP.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_polymorph_acc}}

get_SNPmatrix_acc <- function(gwas_table, SNPrank, SNPmatrix){
  if(is.null(SNPmatrix)) {
    stop("No SNPmatrix specified")
  }
  dats <- gwas_table %>%
    dplyr::arrange(dplyr::desc(log10_p)) %>%
    dplyr::slice(SNPrank) %>%
    dplyr::select(chrom, pos)

  return(SNPmatrix %>%
           dplyr::filter(chrom ==  dats$chrom, pos == dats$pos) %>%
           dplyr::select(-chrom, -pos) %>%
           t %>%
           magrittr::set_colnames("SNP") %>%
           as.data.frame() %>%
           tibble::rownames_to_column("ACC_ID") %>%
           dplyr::filter(SNP == 1) )
}

#' Get intersection between phenotype, and local SNPmatrix SNP calls.
#' @details  Based on a table of phenotypes, a phenotype name, a GWAS table and a rank, and a SNPmatrix returns a table of phenotype values for that gene, where accessions that contain the SNP have TRUE in hasSNP
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @param SNPmatrix The SNPmatrix to use for identifying accessions that carry the relevant SNP.
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_expression_snp}}
#' @seealso \code{\link{get_SNPmatrix_acc}}

intersect_phenotype_snpmatrix <- function(phenotype_table, phenotype, gwas_table, SNPrank, acc_col = "ACC_ID", specific = NULL, SNPmatrix) {
  if(!is.null(specific)){
    get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col, specific = specific) %>%
      dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_SNPmatrix_acc(gwas_table, SNPrank)$ACC_ID ~ TRUE,
                                              TRUE ~ FALSE))
  } else {
    get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col) %>%
      dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_SNPmatrix_acc(gwas_table, SNPrank)$ACC_ID ~ TRUE,
                                              TRUE ~ FALSE))
  }
}

#' Split phenotype table by SNP presence and plot
#' @details Based on a Phenotype table, the name of the phenotype and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{read_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param SNPmatrix The SNPmatrix to use for identifying accessions that carry the relevant SNP.
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @param nobees Set to true to disable beeswarm geom
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_phenotype_snp}}

plot_intersect_phenotype_snpmatrix <- function(phenotype_table, phenotype, gwas_table, SNPrank, SNPmatrix, acc_col = "ACC_ID", specific = NULL, nobees = FALSE){
  if(nobees){
    overplot_geom <- geom_point(alpha = 0.3)
  } else{
    overplot_geom <- ggbeeswarm::geom_beeswarm(alpha = 0.3)
  }
  if(is.null(specific)){
    p <-  intersect_phenotype_snpmatrix(phenotype_table = phenotype_table,
                                        phenotype = phenotype,
                                        gwas_table = gwas_table,
                                        SNPrank = SNPrank,
                                        specific = specific) %>%
      ggplot(aes(x = hasSNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = hasSNP)) +
      stat_summary(color = "purple") +
      overplot_geom +
      labs(title = paste0("Phenotype values by SNP presence"),
           x = "SNP present",
           y = "Value") +
      theme_bw()
  } else{
    p <-  intersect_phenotype_snpmatrix(phenotype_table = phenotype_table,
                                        phenotype = phenotype,
                                        gwas_table = gwas_table,
                                        SNPrank = SNPrank,
                                        specific = specific) %>%
      ggplot(aes(x = hasSNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = hasSNP)) +
      stat_summary(color = "purple") +
      overplot_geom +
      labs(title = paste0("Phenotype values by SNP presence"),
           x = "SNP present",
           y = "Value") +
      theme_bw() +
      facet_grid(reformulate( specific, "Phenotype")) #????
  }

  return(p)
}
