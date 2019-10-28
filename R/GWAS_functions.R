# A collection of functions for working with GWAS tables
# Niklas Schandry, Patrick Hüther 2019



#' This function has no real purpose.
#' @param gwas A gwas result table

format_gwas <- function(gwas){
  gwas %<>% dplyr::mutate(chrom = as.double(stringr::str_extract(chrom, "[0-9]")),
                   log10_p = score,
  )
  return(gwas_object)
}

#' Read a GWAS result file from limix (csv) and put it into a specific (somewhat arbitrary) format.
#'
#' `read_gwas` returns a tibble containing SNP positions, -log10 transformed p-values and information whether
#' a particular SNP passes a significance threshold after multiple testing correction (Bonferroni or Benjamini-Hochberg)
#' @param gwas_path Path to a gwas result file containing SNP positions and p-values

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
#' @param gwas_table Object returned from read_gwas() function
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
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}

plot_acc_map <- function(gwas_table, SNPrank){
  allAccessions <- readr::read_csv("~/labshare/lab/accessions/1001genomes-accessions.csv")
  allAccessions %>%
    filter(id %in% get_polymorph_acc(gwas_table, SNPrank)$strain) %>%
    leaflet::leaflet(data=.) %>%
      leaflet::addTiles() %>%
      leaflet::addCircleMarkers(lng=~longitude,
                                lat=~latitude,
                                label=~as.character(id),
                                popup=~as.character(name),
                                stroke=FALSE, radius=5, fillOpacity=0.8, color="#007243")
}

#' Get expression data for a gene of interest
#' @param GeneID An Arabidopsis thaliana gene identifier

get_expression <- function(GeneID){
  if(!exists("Kawakatsu_dat")){
    Kawakatsu_dat <- readr::read_csv("https://arapheno.1001genomes.org/static/rnaseq/Epigenomic_Diversity_in_A._thaliana_(Kawakatsu_et_al._2016).csv") %>%
      tidyr::as_tibble() %>%
      magrittr::set_colnames(c("ACC_ID", colnames(.)[2:ncol(.)]))
  }
  if(!paste(GeneID) %in% colnames(Kawakatsu_dat)){
    stop(paste0("No expression data in Kawakatsu2016 for GeneID: ",GeneID))
  }
    Kawakatsu_dat %>%
      dplyr::select(ACC_ID, eval(GeneID)) %>%
      tidyr::pivot_longer(dplyr::starts_with("AT"), names_to = "GeneID", values_to = "phenotype_value")
}

#' Based on a GeneID and a GWAS table and a rank, returns a table of expression values for that gene,
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param GeneID An Arabidopsis thaliana gene identifier
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso  \code{\link{get_expression}}
#' @seealso  \code{\link{read_gwas}}

intersect_expression_snp <- function(GeneID, gwas_table, SNPrank){
  get_expression(GeneID) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table,SNPrank)$strain ~ TRUE,
                                            TRUE ~ FALSE))
}

#' Based on a GWAS table and a rank, returns a table of expression values for the gene that is closest to that SNP
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{get_nearest_genes}}

retrieve_counts <- function(gwas_table, SNPrank){
  genes <- get_nearest_genes(gwas_table, SNPrank) %>%
             dplyr::slice(SNPrank) %>%
             .$GeneId
  get_expression(genes) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                                             TRUE ~ FALSE))
}

#' Based on a GWAS table and a rank, returns a table of expression values for the gene that is closest to that SNP
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{get_nearest_genes}}
#' @seealso \code{\link{retrieve_counts}}

plot_intersect_expression_snp <- function(gwas_table, SNPrank){

  p <-  retrieve_counts(gwas_table, SNPrank) %>%
    ggplot(aes(x = hasSNP, y = phenotype_value)) +
    geom_boxplot(aes(fill = hasSNP)) +
    ggbeeswarm::geom_beeswarm(alpha = 0.3) +
    labs(title = paste0("Expression of nearest gene by SNP presence"),
         caption = "Expression data from araPheno",
         x = "SNP present",
         y = "Value") +
    theme_bw() +
    facet_wrap(~GeneID)
  print(p)
}


#################
# Below are functions that are generalizations of the expression specific functions.
# They use a phenotype table, which comes in wide format and extract phenotype values.
# These phenotypes are matched to a specific SNP in a GWAS table with intersect_phenotype_snp, giving a presence | abscence factor
# Plotting is available via plot_intersect_phenotype_snp.
#################

#get_phenotype is a more general version of get_expression, that works with custom phenotype tables, in wide format.

#' Subset a larger wide-format phenotype table to a defined phenotype.
#' Optionally includes a "specific" variable.
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

#' Get intersection between phenotype, and variant tale.
#' Based on a table of phenotypes, a phenotype name, a GWAS table and a rank, returns a table of phenotype values for that gene,
#' where accessions that contain the SNP have TRUE in hasSNP
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_expression_snp}}


intersect_phenotype_snp <- function(phenotype_table, phenotype, gwas_table, SNPrank, acc_col = "ACC_ID", specific = NULL) {
  get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col, specific = specific) %>%
    dplyr::mutate(hasSNP = dplyr::case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                              TRUE ~ FALSE))
}

#' Based on a GeneID and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from read_gwas() function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{intersect_phenotype_snp}}

plot_intersect_phenotype_snp <- function(phenotype_table, phenotype, gwas_table, SNPrank, acc_col = "ACC_ID", specific = NULL){
  p <-  intersect_phenotype_snp(phenotype_table, phenotype , gwas_table, SNPrank, specific = specific)  %T>% print %>%
    ggplot(aes(x = hasSNP, y = phenotype_value)) +
    geom_boxplot(aes(fill = hasSNP)) +
    ggbeeswarm::geom_beeswarm(alpha = 0.3) +
    labs(title = paste0("Phenotype values by SNP presence"),
         x = "SNP present",
         y = "Value") +
    theme_bw() +
    facet_grid(reformulate( specific, "Phenotype")) #????
  print(p)
}


#' Based on a GWAS table, returns the gene annotation that is closest to each SNP for the number of SNP specified.
#' This function always returns the closest annotation. To limit the lookup range use get_overlapping_genes()
#' Lookup is done via ensembl plants; requires internet connection.
#' !This function will assign araGenes in the Global Environment!
#' @param gwas_table Object returned from read_gwas() function
#' @param n_hit The number of SNPs that should be looked up.
#' @seealso \code{\link{read_gwas}}
#' @seealso \code{\link{get_overlapping_genes}}



get_nearest_genes <- function(gwas_table = NULL, n_hit = 1){
  if(is.null(gwas_table)){
    stop("GWAS output missing")
  }
  # get annotation info from plants_mart
  if(!exists("araGenes")){
    ara <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
    # Get genome, rename columns, make into granges
    araGenes <<- biomaRt::getBM(c("ensembl_gene_id",
                                  "chromosome_name",
                                  "start_position",
                                  "end_position",
                                  "strand",
                                  "description",
                                  "transcript_biotype"),
                                mart=ara) %>%
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

#' Based on a GWAS table, returns the gene annotation that is overlapping with each SNP for the number of SNP specified.
#' If there is no overlap, the SNP is returned, but the GeneID column is empty.
#' Lookup is done via ensembl plants; requires internet connection.
#' !This function will assign araGenes in the Global Environment!
#' @param gwas_table Object returned from read_gwas() function
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
    ara <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
    # Get genome, rename columns, make into granges
    araGenes <<- biomaRt::getBM(c("ensembl_gene_id",
                                  "chromosome_name",
                                  "start_position",
                                  "end_position",
                                  "strand",
                                  "description",
                                  "transcript_biotype"),
                                mart=ara) %>%
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

  #araGenes %>% filter_by_overlaps(snps, maxgap = distance)
}

#' Based on a GWAS table, generates a manhatten plot.
#' For performance reasons, everything with a log10(p) smaller than p_filter is filtered out.
#' @param gwas_table Object returned from read_gwas() function
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
  ggplot(aes(x=pos, y=log10_p), data = gwas_table %>% dplyr::filter(log10_p > p_filter, mac > mac_filter)) +
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

#' Based on a GWAS table, generates a manhatten plot, with gene annotations. Number of annotations can be toggled
#' by changing nlabels. By default only plots annotations for SNPs that are within a gene annotation.
#' For performance reasons, everything with a log10(p) smaller than p_filter is filtered out.
#' @param gwas_table Object returned from read_gwas() function
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
    filter(abs(log10_p) > p_filter) %>%
    filter(mac > mac_filter) %>%
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
