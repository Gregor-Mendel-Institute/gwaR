# A collection of functions for working with GWAS tables
# Niklas Schandry, Patrick Hüther 2019


#' The araGenes object.
araGenome <- biomaRt::getBM(c("ensembl_gene_id",
                             "chromosome_name",
                             "start_position",
                             "end_position",
                             "strand",
                             "description",
                             "transcript_biotype"),
                           mart = biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org'))

#' A nicely formatted GRanges object
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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

#' Plot an interactive map of accessions that carry a SNP of interest
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}

plot_acc_map <- function(gwas_table, SNPrank){
  allAccessions <- readr::read_csv("~/labshare/lab/accessions/1001genomes-accessions.csv")
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
    ) %>% as_data_frame()
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
    dplyr::as_data_frame() %>%
    dplyr::mutate(ACC_ID = accession_id)
 }


#' Based on a GeneID and a GWAS table and a rank, returns a table of expression values for that gene,
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param GeneID An Arabidopsis thaliana gene identifier
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
  print(p)
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
  }
}

#' Split phenotype table by SNP presence and plot
#' @details Based on a Phenotype table, the name of the phenotype and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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

  print(p)
}

#' Find genes that are closest to a SNP.
#' @details Based on a GWAS table, returns the gene annotation that is closest to each SNP for the number of SNP specified. This function always returns the closest annotation. To limit the lookup range use {/link get_overlapping_genes()}  Lookup is done via ensembl plants; requires internet connection.
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
#' @param title Specify plot title
#' @param subtitle Specify plot subtitle
#' @param p_filter everything with a log10(p) below this value will not be included in the plot
#' @param mac_filter everything with a mac (minor allele count) below this will not be plotted.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{plot_annotated_gwas}}



plot_gwas <- function(gwas_table, title = "No Title", subtitle = NULL, p_filter = 2, mac_filter = 0 ){
  color_gmi_light <- ("#abd976")
  color_gmi_dark <- ("#007243")
  GWAS_colors <- c(color_gmi_dark, color_gmi_light, "grey50")
  names(GWAS_colors) <- c("Bonferroni", "FDR", "Not")
  ggplot(aes(x=pos, y=log10_p), data = gwas_table %>% dplyr::filter(log10_p > p_filter) %>% dplyr::filter(mac > mac_filter)) +
    # geom_hline(linetype = "dotted", yintercept = bf_corr) +
    facet_grid(~chrom, scales = "free_x", switch = "x") +
    #  geom_label_repel(aes(x=pos, y= -log10(pv),label = gene), data = ft_specific_limix_mac5 %>% filter(gene != "NA")) +
    ggthemes::theme_few() +
    #  geom_mark_rect(aes(filter = gene != 'NA', label = gene)) +
    geom_point(aes(color = Significant)) +
    scale_color_manual(values = GWAS_colors) +
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
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
  color_gmi_light <- ("#abd976")
  color_gmi_dark <- ("#007243")
  GWAS_colors <- c(color_gmi_dark, color_gmi_light, "grey50")
  names(GWAS_colors) <- c("Bonferroni", "FDR", "Not")
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
    scale_color_manual(values = GWAS_colors) +
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
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
#' @param rank Rank of the SNP of interest
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_depth Maximum SNP distance to calculate LD for (only relevant when anchored = FALSE)
#' @param ld_stats The LD statistics, see SNPStats::ld
#' @param ld_symmetric Should a symmetric matrix be returned? see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param anchored Perform anchored analysis. If TRUE Linkage will only be estimated for the SNP of interest, and the surrounding ones, but not between surrounding SNPs.

snp_linkage <- function(gwas_table,
                        rank,
                        nuc_range = 50000,
                        ld_depth = "1000",
                        ld_stats = c("D.prime", "R.squared"),
                        ld_symmetric = FALSE,
                        use_phenotype_table = NULL,
                        acc_col = "ACC_ID",
                        anchored = FALSE) {

  # Step 1: Download variant table
  ## Define region for API call

  region <- gwas_table %>%
    dplyr::slice(rank) %>%
    {paste0("Chr", .$chrom, ":", .$pos - as.numeric(nuc_range) / 2, ".." , .$pos + as.numeric(nuc_range) / 2)}

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
  subset_url <- paste0("http://tools.1001genomes.org/api/v1/vcfsubset/strains/",
                       genotypes,
                       "/regions/",
                       region,
                       "/type/fullgenome/format/vcf")

  ## Make tempfile

  tmp <- tempfile(fileext = ".vcf")

  ## Write vcf table to tempfile

  writeLines( httr::content( httr::GET(subset_url ) ), tmp)

  message(paste0("Downloaded genotype vcf to ", tmp))
  # Step 2 Read vcf table from tempfile
  # vcf_info <- pegas::VCFloci(tmp)


  #variants_for_linkage <- pegas::read.vcf(tmp, to = nuc_range + 1)
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
    message(paste0("Anchored analysis, with ", region,  " as anchor"))
    # Start with getting the anchor information, similar to above, but the region is only one position.
    anchor_url <-  paste0("http://tools.1001genomes.org/api/v1/vcfsubset/strains/",
                          genotypes,
                          "/regions/",
                          region,
                          "/type/fullgenome/format/vcf")

    anchor_tmp <- tempfile(fileext = ".vcf")

    writeLines( httr::content( httr::GET(anchor_url ) ), anchor_tmp)

    message(paste0("Downloaded anchor vcf to ", anchor_tmp))

    anchor_vcf <- VariantAnnotation::readVcf(anchor_tmp)

    anchor_vcf <- anchor_vcf[VariantAnnotation::isSNV(anchor_vcf)]

    anchor_SM <- VariantAnnotation::genotypeToSnpMatrix(anchor_vcf)

    anchor_SM <- anchor_SM$genotypes

    ld_tab <- snpStats::ld(x = anchor_SM , y = SM_for_linkage, stats = ld_stats)
  }
  return(ld_tab)
}

#' Produce a plot of linked SNPs
#' @details This function takes a gwas table and a rank and produces a plot that illustrates, SNPs that are within nuc_range around the SNP, the degree of linkage, the impact of the other SNPs and a track of gene annotations to easily identify linked SNPs that have a high impact.
#' @details Impact according to 1001genomes.org
#' @param gwas_table Object returned from \code{\link{read_gwas()}} function
#' @param rank Rank of the SNP of interest
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_stats The LD statistics, see SNPStats::ld
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param linkage_cutoff Only SNPs that have an LD values >= this will be plotted (default 0.8)
#' @param LD_legend Toggle color legend for LD. Off by default.
#' @param data_only If true does not plot, only returns the data used for plotting.

plot_anchored_ld <-  function(gwas_table,
                              rank,
                              nuc_range = 50000,
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

  anc_ld <- snp_linkage(gwas_table = gwas_table, rank = rank, nuc_range = nuc_range, anchored = TRUE, use_phenotype_table = use_phenotype_table, ld_stats = ld_stats)

  # Define SNP details.

  snp_pos <- gwas_table %>%
    dplyr::slice(rank) %>% .$pos

  chrom = gwas_table %>%
    dplyr::slice(rank) %>%
    .$chrom

  start_pos = gwas_table %>%
    dplyr::slice(rank) %>% {.$pos - as.numeric(nuc_range) / 2}

  end_pos = gwas_table %>%
    dplyr::slice(rank) %>% {.$pos + as.numeric(nuc_range) / 2}

  ## Define genotypes for API call
  if(is.null(use_phenotype_table)){
    message("Calculating LD based on strains that carry SNP")
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
    ), col_types = "iifcccccccccccc")

  ## Labels for the gene plots

  gene_labels <- araGenome %>%
    dplyr::filter(chromosome_name == chrom, start_position %in% c(start_pos:end_pos), end_position %in% c(start_pos:end_pos)) %>%
    dplyr::mutate(gene = ensembl_gene_id,
           molecule = -0.25,
           start = start_position,
           end = end_position,
           direction = strand,
           strand = dplyr::case_when(strand == -1 ~ "reverse", TRUE ~ "forward"))

  ## Transforming LD table into more plotable table.

  plot_data <- anc_ld %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(tidyselect::matches(":")) %>%
    na.omit %>%
    dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(name, "[:|_]",3)[,2])) %>%
    dplyr::left_join(., snp_impacts)

  ## Build Plot

  p <- ggplot(data = plot_data) +
    # Tile Geom for LD values
    geom_tile(aes(x = pos , color = value, y = 0, height = 0.04)) +
    # Line denoting SNP of interest
    geom_vline(aes(xintercept = snp_pos), color = "darkred") +
    # Modifier SNPs
    geom_point(aes(x= pos, y = -0.2, shape = effect_impact), data = {. %>% dplyr::filter(value >= linkage_cutoff, effect_impact == "MODIFIER")}) +
    # LOW SNPs
    geom_point(aes(x= pos, y = -0.15, shape = effect_impact), data = {. %>% dplyr::filter(value >= linkage_cutoff, effect_impact == "LOW")}) +
    # MODERATE SNPs
    geom_point(aes(x= pos, y = -0.1, shape = effect_impact), data = {. %>% dplyr::filter(value >= linkage_cutoff, effect_impact == "MODERATE")}) +
    # HIGH SNPs
    geom_point(aes(x= pos, y = -0.05, shape = effect_impact), data = {. %>% dplyr::filter(value >= linkage_cutoff, effect_impact == "HIGH")}) +
    # color scale
    scale_color_viridis_c(option = "plasma", direction = -1) +
    # Gene Arrows
    gggenes::geom_gene_arrow(aes(xmin = start, xmax = end, y = molecule, fill = gene, forward = direction), data = gene_labels) +
    # Gene labels
    gggenes::geom_gene_label(aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene), data = gene_labels) + # Gene labels
    # Minimal Theme
    theme_minimal() +
    # Theme adjusments, mainly removing the y-axis
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(color = LD_legend) # Color guide

    if(data_only){
      return(plot_data)
    } else{
      return(p)
    }
}




