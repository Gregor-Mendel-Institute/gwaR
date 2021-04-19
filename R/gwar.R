# A collection of functions for working with GWAS tables
# Niklas Schandry, Patrick Hüther 2019-2021

#' A list of sequenced accessions as of 2020
#' @format A single column of accesion ids
#' \describe{
#' \item{ACC_ID}{Accession IDs}
#' }
sequenced_accessions <- read.csv("data/accessions.csv")

#' limix example output; 500 random rows
#' @format 500 rows
#'  \describe{
#'   \item{chrom}{Chromosome, numeric}
#'   \item{pos}{Position on the chromosome}
#'   \item{pv}{p-value}
#'   \item{maf}{minor allele frequency}
#'   \item{mac}{minor allele count}
#' }
limix500rows <- read.csv("data/example_limix.csv")

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
#' @details This function is a flexible version of \code{\link{format_gwas}}
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
#' @details `format_gwas` returns a tibble containing SNP positions, -log10 transformed p-values and information whether a particular SNP passes a significance threshold after multiple testing correction (Bonferroni or Benjamini-Hochberg)
#' @param gwas_path Path to a gwas result file containing SNP positions and p-values
#' @seealso \code{\link{format_gwas}}

read_limix <- function(gwas_path){
  gwas_object <- vroom::vroom(gwas_path)
  if(nrow(gwas_object) < 10000){
  message("Table has less than 10000 rows, no FDR correction")
  fdr_thresh <- 10^6
  } else{
  fdr_corr <- qvalue::qvalue(p = gwas_object$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05 )]) %>% max() %>% log10() %>% abs()
  }
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

#' Find genes that are closest to a SNP.
#' @details Based on a GWAS table, returns the gene annotation that is closest to each SNP for the number of SNP specified. This function always returns the closest annotation. To limit the lookup range use \code{\link{get_overlapping_genes}}. Lookup is done via ensembl plants; requires internet connection.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param n_hit The number of SNPs that should be looked up.
#' @seealso \code{\link{format_gwas}}
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
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param n_hit The number of SNPs that should be looked up.
#' @param distance The maximum distance from SNP to annotation, can be varied to lookup genes within a specific distance
#' @seealso \code{\link{format_gwas}}
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

#' Plot annotated manhattan plot
#' @details Based on a GWAS table, generates a manhatten plot, with gene annotations. Number of annotations can be toggled by changing nlabels. By default only plots annotations for SNPs that are within a gene annotation. For performance reasons, everything with a log10(p) smaller than p_filter is filtered out.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param title Specify plot title
#' @param subtitle Specify plot subtitle
#' @param p_filter everything with a log10(p) below this value will not be included in the plot (default: 2)
#' @param mac_filter everything with a mac (minor allele count) below this will not be plotted. (default: 0)
#' @param annotated Should the plot contain annotations for some peaks (controls see below)? Default is TRUE
#' @param nlabels How many labels should be plotted (default: 5)
#' @param labeltype Which entry from the lookup should be displayed (character, defaults to "GeneId"). Options: SNP_rank, chrom, pos, Significant, log10_p, mac, GeneId, max_distance, description, transcript_biotype, start_position, end_position.
#' @param match_nearest Defaults to FALSE; if TRUE will use find_nearest_genes instead of find_overlapping_genes for annotations
#' @seealso \code{\link{format_gwas}}
#' @seealso \code{\link{plot_gwas}}
#' @seealso \code{\link{get_nearest_genes}}
#' @seealso \code{\link{get_overlapping_genes}}

plot_gwas <- function(gwas_table,
                      title ="No Title",
                      subtitle = NULL,
                      p_filter = 2,
                      mac_filter = 0,
                      nlabels = 5,
                      annotated = TRUE,
                      labeltype = "GeneId",
                      match_nearest = FALSE) {
  if(!match_nearest){
    annotations <- gwas_table %>% get_overlapping_genes(nlabels) %>% tidyr::unite(labs, SNP_rank, labeltype, sep = " :")
  } else {
    annotations <- gwas_table %>% get_nearest_genes(nlabels) %>% tidyr::unite(labs, SNP_rank, labeltype, sep = " :")
  }
  if(annotated){
  plot_annots <- ggrepel::geom_text_repel(aes(label = labs), data = annotations)
  } else{
      plot_annots <- NULL
  }
  gwas_table %>%
    dplyr::filter(abs(log10_p) > p_filter) %>%
    dplyr::filter(mac > mac_filter) %>%
    #Step 2: Plot
    ggplot(aes(x=pos, y=log10_p)) +
    geom_point(aes(color = Significant)) +
    plot_annots +
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



#' Obtain accession ids that carry a SNP.
#' This function takes a table of gwas pvalues per chromosome and numeric position and retrieves the accessions that carry that SNP.
#' @details If no SNPmatrix is provided, the information will be obtained from 1001genomes.org
#' @details This function returns a single column named ACC_ID, that contains the accessions that contain the snp.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param SNPmatrix A path pointing to a SNPmatrix in .fst format.

get_accessions <- function(gwas_table, SNPrank, SNPmatrix = NULL){

  ## Define chromosome and position

  dats <- gwas_table %>%
    dplyr::arrange(dplyr::desc(log10_p)) %>%
    dplyr::slice(SNPrank) %>%
    dplyr::select(chrom, pos)

  ## If there is no SNPmatrx path provided, use polymorph

  if(is.null(SNPmatrix)){
    message("No SNPmatrix provided, retrieving information from 1001genomes.org")
    acc_table <- httr::content(httr::GET(paste0("http://tools.1001genomes.org/api/v1.1/variants.csv?type=snps;accs=all;chr=",
                                                dats$chrom,
                                                ";pos=",
                                                dats$pos)), col_types = "iicccccc") %>%
      dplyr::mutate(ACC_ID = strain) %>%
      dplyr::select(ACC_ID)
  } else{

    # Check SNPmatrix

    if(!file.exists(eval(SNPmatrix))) {
      stop("Please point to snpmatrix in fst format")
    }

    #Subsetting

    SNPmat <- fst::fst(SNPmatrix)
    SNPmat <- SNPmat[(SNPmat$chrom == dats$chrom) & (SNPmat$pos == dats$pos),]
    acc_table <- SNPmat %>%
      dplyr::as_tibble() %>%
      dplyr::select(-chrom,-pos) %>%
      t
    acc_table <- acc_table %>%
      magrittr::set_colnames("SNP") %>%
      as.data.frame() %>%
      tibble::rownames_to_column("ACC_ID") %>%
      dplyr::filter(SNP == 1) %>%
      dplyr::select(-"SNP")
  }

  return(acc_table)

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
      dplyr::select("ACC_ID", eval(phenotype)) %>%
      tidyr::pivot_longer(tidyselect::matches(eval(phenotype)), names_to = "Phenotype", values_to = "phenotype_value")
  } else {
    phenotype_table %>%
      dplyr::select("ACC_ID", eval(phenotype), eval(specific)) %>%
      tidyr::pivot_longer(tidyselect::matches(eval(phenotype)), names_to = "Phenotype", values_to = "phenotype_value")
  }

}

#' Get intersection between phenotype, and SNP presence.
#' @details Based on a table of phenotypes, a phenotype name, a GWAS table and a rank, and a SNPmatrix returns a table of phenotype values for that gene, where accessions that contain the SNP have TRUE in SNP.
#' @details Also provides a simple default plot, toggled via the plot argument
#' @param phenotype_table a table containing phenotyping measurements, and accession ids (see below)
#' @param phenotype a specific phenotype from the phenotype table. Must match to a column name of the phenotype table
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param acc_col the column that contains accession identifiers.
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @param plot Default: FALSE. Toggle if a plot should be returned
#' @param nobees Set to true to disable beeswarm geom (only relevant if plot = TRUE)
#' @param SNPmatrix, Provide SNPmatrix to control if accessions are based on 1001genomes.org or local matrix (default: NULL -> use 1001genomes).
#' @param ... Arguments passed to other functions. Most relevant is
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{format_gwas}}
#' @seealso \code{\link{expression_by_snp}}
#' @seealso \code{\link{get_accessions}}

phenotype_by_snp <- function(phenotype_table,
                             phenotype,
                             gwas_table,
                             SNPrank,
                             acc_col = "ACC_ID",
                             specific = NULL,
                             plot = FALSE ,
                             nobees = FALSE,
                             SNPmatrix = NULL,
                             ...) {

  #Join phenotypes and snp information

  if(!is.null(specific)){
   result <- get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col, specific = specific) %>%
      dplyr::mutate(SNP = dplyr::case_when(ACC_ID %in% get_accessions(gwas_table, SNPrank, SNPmatrix = SNPmatrix)$ACC_ID ~ TRUE,
                                              TRUE ~ FALSE))
  } else {
  result <-  get_phenotype(phenotype_table = phenotype_table, phenotype = phenotype, acc_col = acc_col) %>%
      dplyr::mutate(SNP = dplyr::case_when(ACC_ID %in% get_accessions(gwas_table, SNPrank, SNPmatrix = SNPmatrix)$ACC_ID ~ TRUE,
                                              TRUE ~ FALSE))
  }

  #Return table if plot = F

  if(!plot){
    return(result)
  } else {

    # If plot is true
    # Define overplot geom (bees or nobees)
    if(nobees){
    overplot_geom <- geom_point(alpha = 0.3)
    } else {
    overplot_geom <- ggbeeswarm::geom_beeswarm(alpha = 0.3)
      }
    # Define facetting, this depends on specific
    if(is.null(specific)){
      plot_facets <- NULL
    } else {
      plot_facets <- facet_grid(reformulate(specific, "Phenotype"))
    }
    # Plot
    p <- result %>%
      ggplot(aes(x = SNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = SNP)) +
      overplot_geom +
      labs(title = paste0("Phenotype values by SNP presence"),
           x = "SNP present",
           y = "Value") +
      theme_bw() +
      plot_facets
    return(p)
  }
}

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

#' Based on a GWAS table and a rank, returns a table of expression values for a gene where accessions that contain the SNP have TRUE in SNP.
#' If no GeneID is specified, the gene that is closest to the SNP will be looked up.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param GeneID If NULL (default), retrieves counts for the gene that is nearest. If not NULL it should be a TAIR GeneID ("ATXGNNNNN"). If a GeneID is provided expression of that gene will be intersected with the SNP presence.
#' @param plot If TRUE, will return a plot instead of a table.
#' @param nobees If TRUE, observations will be plotted with geom_point.
#' @param SNPmatrix SNPmatrix, to control if accessions are based on 1001genomes.org or local matrix.
#' @seealso \code{\link{format_gwas}}
#' @seealso \code{\link{get_expression}}
#' @seealso \code{\link{get_nearest_genes}}

expression_by_snp <- function(gwas_table,
                              SNPrank,
                              GeneID = NULL,
                              plot = FALSE,
                              nobees = FALSE,
                              SNPmatrix = NULL,
                              ...){
 if(is.null(GeneID)){
  genes <- get_nearest_genes(gwas_table, SNPrank) %>%
    dplyr::slice(SNPrank) %>%
    .$GeneId
    } else{
    genes <- GeneID
    }
  results <- get_expression(GeneID = genes) %>%
    dplyr::mutate(SNP = dplyr::case_when(ACC_ID %in% get_accessions(gwas_table, SNPrank, SNPmatrix = SNPmatrix)$ACC_ID ~ TRUE,
                                            TRUE ~ FALSE))
  if(!plot){
    return(results)
  } else {
    if(nobees){
      overplot_geom <- geom_point(alpha = 0.3)
    } else {
      overplot_geom <- ggbeeswarm::geom_beeswarm(alpha = 0.3)
    }
    p <- results %>%
      ggplot(aes(x = SNP, y = phenotype_value)) +
      geom_boxplot(aes(fill = SNP)) +
      overplot_geom +
      labs(title = paste0("Expression of ", genes," by SNP presence"),
           caption = "Expression data from araPheno",
           x = "SNP present",
           y = "Value") +
      theme_bw() +
      facet_wrap(~phenotype_name)
    return(p)
  }
}


#' Plot an interactive map of accessions that carry a SNP of interest, needs a list of accessions with lng, lat fields
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param SNPrank The (-log10(p)) rank of the SNP of interest
#' @param SNPmatrix provide path to control if accessions are based on 1001genomes.org or local matrix.
#' @seealso \code{\link{format_gwas}}

plot_acc_map <- function(gwas_table, SNPrank, SNPmatrix){

  sequenced_accessions %>%
    dplyr::filter(id %in% get_accessions(gwas_table, SNPrank, SNPmatrix = SNPmatrix)$ACC_ID) %>%
    leaflet::leaflet(data=.) %>%
    leaflet::addTiles() %>%
    leaflet::addCircleMarkers(lng=~longitude,
                              lat=~latitude,
                              label=~as.character(id),
                              popup=~as.character(name),
                              stroke=FALSE, radius=5, fillOpacity=0.8, color="#007243")
}

#' Retrieve impacts of SNPs within a certain range around a GWAS peak. Uses snpmatrix for genotype lookup
#' @details Based on a Phenotype table, the name of the phenotype and a GWAS table and a rank, produces a boxplot of that phenotype, grouped by presence of that SNP.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param rank The (-log10(p)) rank of the SNP of interest
#' @param nuc_range Window around the SNP of interest to retrieve information for.
#' @param SNPmatrix The path to the snpmatrix to use for identifying accessions that carry the relevant SNP.
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param acc_col the column that contains accession identifiers.
#' @seealso \code{\link{get_phenotype}}
#' @seealso \code{\link{format_gwas}}
#' @seealso \code{\link{phenotype_by_snp}}

polymorph_impact <- function(gwas_table, rank, nuc_range, SNPmatrix, use_phenotype_table = NULL, acc_col){

  gwas_table <- gwas_table %>%
    dplyr::arrange(dplyr::desc(log10_p))

  snp_pos <- gwas_table %>%
    dplyr::slice(rank) %>% .$pos

  chrom = gwas_table %>%
    dplyr::slice(rank) %>%
    .$chrom

  start_pos = gwas_table %>%
    dplyr::slice(rank) %>% {
      .$pos - as.numeric(nuc_range) / 2
    }

  if (start_pos < 1) {
    start_pos <-
      1 # No message, because snp_linkage will issue a message.
  }
  end_pos = gwas_table %>%
    dplyr::slice(rank) %>% {
      .$pos + as.numeric(nuc_range) / 2
    }

  ## Define genotypes

  if (is.null(use_phenotype_table)) {
    message("Retrieving SNP impacts for strains in SNPmatrix. Strains are assumed to be numeric.")
    genotypes <-
      colnames(fst::fst(SNPmatrix)) %>% {
        suppressWarnings(as.numeric(.))
      } %>% na.omit() %>% stringr::str_flatten(., ",")
  }

  if (!is.null(use_phenotype_table)) {
    message("Retrieving SNP impacts for strains in Phenotype table.")
    genotypes <-
      stringr::str_flatten(levels(as.factor(unlist(
        use_phenotype_table[, eval(acc_col)]
      ))), collapse = ",")
  }

  ## Retrieve SNPimpacts from polymorph DB
  snp_impacts <- httr::content(httr::GET(
    paste0(
      "http://tools.1001genomes.org/api/v1.1/effects.csv?accs=",
      genotypes,
      ";chr=",
      chrom,
      ";start=",
      start_pos,
      ";end=",
      end_pos,
      ";type=snps"
    )
  ), col_types = "iifccccccccccc") %>%
    dplyr::select(chr, pos,strain,type, effect_impact)
  return(snp_impacts)
}


## Linkage Related Functions.
# Niklas Schandry

#' Calculate linkage around SNP of interest.
#' @details This function calculates LD in a specific region.
#' @details The SNPs per accession in the region are retrieved from polymorph or from a provided SNPmatrix.
#' @details Sometimes, the ld function produces no result, especially when used with anchored = T.
#' @details anchored analysis will use one SNP can compute LD of that SNP against a region.
#' @details Default settings return the LD table.
#' @param gwas_table Object returned from \code{\link{format_gwas}} function
#' @param rank Rank of the SNP of interest
#' @param anchored Perform anchored analysis. If TRUE Linkage will only be estimated for the SNP of interest, and the surrounding ones, but not between surrounding SNPs.
#' @param plot Default FALSE. Should a plot be returned, instead of the table? Only supported for anchored = T
#' @param SNPmatrix default NULL. Needs to provided as a path to an fst file. Will be used instead of polymorph.
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param ld_depth Maximum SNP distance to calculate LD for (only relevant when anchored = FALSE)
#' @param ld_stats The LD statistics, see SNPStats::ld. Default is LLR
#' @param ld_symmetric Should a symmetric matrix be returned? see SNPStats::ld
#' @param ld_cutoff Only SNPs that have an LD values >= this will be plotted (default 10)
#' @param add_ld_band add a layer showing the LD as a lot of lines. (default FALSE)
#' @param ld_legend Toggle color legend for LD. Off by default.
#' @param use_phenotype_table If supplied: Genotypes listed here will be used for linkage analysis. Otherwise, all accessions that carry this SNP will be included.
#' @param use_all_acc If true, will override use_phenotype_table and use all 1135 sequenced accessions.
#' @param acc_col Name of the column containing accession identifiers in the phenotype table
#' @param data_only Return the data used for plotting, instead of the plot, defaults to FALSE.

snp_linkage <- function(gwas_table,
                        rank,
                        anchored = FALSE,
                        plot = FALSE,
                        SNPmatrix = NULL,
                        nuc_range = 50000,
                        ld_depth = 1000,
                        ld_stats = "LLR",
                        ld_symmetric = FALSE,
                        ld_cutoff = 10,
                        ld_legend = FALSE,
                        add_ld_band = FALSE,
                        use_phenotype_table = NULL,
                        acc_col = "ACC_ID",
                        use_all_acc = FALSE,
                        data_only = FALSE
                        #debug.return.sm = F,
                        #debug.return.sm.anc = F,
                        #debug.return.gt = F
                        ) {

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

  chrom <- gwas_table %>%
    dplyr::slice(rank) %>%
    dplyr::select(chrom) %>%
    unique() %>%
    paste()

  pos = c(region_lower:region_upper)

  if(is.null(SNPmatrix)){

    ## Define genotypes for API call
    if(is.null(use_phenotype_table)){
      message("Calculating LD based on strains that carry SNP")
      genotypes <- stringr::str_flatten(get_accessions(gwas_table = gwas_table, SNPrank = rank)$ACC_ID, collapse = ",")
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
      genotypes <-sequenced_accessions %>% unlist %>%  stringr::str_flatten(collapse = ",")
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

      anchor_vcf <- anchor_vcf[VariantAnnotation::isSNV(anchor_vcf)]

      anchor_SM <- VariantAnnotation::genotypeToSnpMatrix(anchor_vcf)

      anchor_SM <- anchor_SM$genotypes
    }
  }

  # Using a SNPmatrix

  if(!is.null(SNPmatrix)){
    if(!file.exists(SNPmatrix) | !(tools::file_ext(SNPmatrix) == "fst" )) {
      stop("Please point to snpmatrix in fst format")
    }
    ## Define genotypes for subsetting
    if(is.null(use_phenotype_table)){
      message("Calculating LD based on strains in SNPmatrix. Strains are assumed to be numeric.")
      genotypes <- colnames(fst::fst(SNPmatrix)) %>% {suppressWarnings(as.numeric(.))} %>% na.omit()
    }

    if(!is.null(use_phenotype_table)){
      message("Calculating LD based on all strains in phenotype table,")
      genotypes <- levels(as.factor(unlist(use_phenotype_table[, eval(acc_col)])))
      ## construct call
    }
    # Very base R subsetting
    tmpmatrix <- fst::fst(SNPmatrix)

    tmpmatrix <-
      tmpmatrix[(tmpmatrix$chrom == chrom) & (tmpmatrix$pos %in% pos),
                paste(c(genotypes, "chrom", "pos"))]

    rownames(tmpmatrix) <- paste(tmpmatrix$chrom, tmpmatrix$pos, sep = "_")

    tmpmatrix <- tmpmatrix[, paste(genotypes)]

    tmpmatrix <- tmpmatrix %>%
      as.matrix() %>%
      t()

    ## Recode matrix. In snpStats, 0 is coding for missing. In the SNPmatrix 0 is coding for "allele 1".
    ## Contrary to what the documentation seems to say, in a numeric matrix 0 is missing, 1 is allele 1, 2 is hetero, 3 is allele 2.
    storage.mode(tmpmatrix) <- "double"
    tmpmatrix[tmpmatrix == 1] <- 3
    tmpmatrix[tmpmatrix == 0] <- 1
    storage.mode(tmpmatrix) <- "raw"

    ## Change storage mode so coercion works

    SM_for_linkage <- as(tmpmatrix, "SnpMatrix")

    if(anchored){
      region <- gwas_table %>%
        dplyr::slice(rank) %>%
        {paste0("Chr", .$chrom, ":", .$pos , ".." , .$pos)}
      message(paste0("Anchored analysis, with ", region,  " as anchor.\n"))
      # Start with getting the anchor information, similar to above, but the region is only one position.
      anchor_pos <- gwas_table %>%
        dplyr::slice(rank) %>%
        .$pos

      anchor_tmpmatrix <- fst::fst(SNPmatrix)

      anchor_tmpmatrix <- anchor_tmpmatrix[(anchor_tmpmatrix$chrom == chrom) & (anchor_tmpmatrix$pos == anchor_pos),
                                           paste(c(genotypes, "chrom", "pos"))]

      rownames(anchor_tmpmatrix) <- paste(anchor_tmpmatrix$chrom, anchor_tmpmatrix$pos, sep = "_")

      anchor_tmpmatrix <- anchor_tmpmatrix[, paste(genotypes)]

      anchor_tmpmatrix <- anchor_tmpmatrix %>%
        as.matrix() %>%
        t()

      ## Recode matrix. see above
      storage.mode(anchor_tmpmatrix) <- "double"
      anchor_tmpmatrix[anchor_tmpmatrix == 1] <- 3
      anchor_tmpmatrix[anchor_tmpmatrix == 0] <- 1
      storage.mode(anchor_tmpmatrix) <- "raw"

      anchor_SM <-  as(anchor_tmpmatrix, "SnpMatrix")
    }

  }
  # if(debug.return.gt){
  #   return(genotypes)
  # }
  # if(debug.return.sm) {
  #   return(SM_for_linkage)
  # }
  # if(debug.return.sm.anc){
  #   return(anchor_SM)
  # }
  if(anchored){
    ld_tab <- snpStats::ld(SM_for_linkage, anchor_SM, stats = ld_stats)
  }
  if(!anchored){
    ld_tab <- snpStats::ld(SM_for_linkage, depth = ld_depth, stats = ld_stats, symmetric = ld_symmetric)
  }


  if(!plot){
    return(ld_tab)
  }
  if(!anchored & plot){
    message("Non-anchored plots are not supported, returning table")
    return(ld_tab)
  }
  if(anchored & plot){
    if(nrow(na.omit(ld_tab)) == 0) {
      stop("LD calculation returned only NAs.")
    }
    snp_pos <- gwas_table %>%
      dplyr::slice(rank) %>%
      .$pos
    gene_labels <- araGenome %>%
      dplyr::filter(chromosome_name == chrom,
                    start_position %in% c(region_lower:region_upper),
                    end_position %in% c(region_lower:region_upper)) %>%
      dplyr::mutate(gene = ensembl_gene_id,
                    molecule = 0,
                    start = start_position,
                    end = end_position,
                    direction = strand,
                    strand = dplyr::case_when(strand == -1 ~ "reverse", TRUE ~ "forward"),
                    layer_var = "ORF")
    impacts <- polymorph_impact(gwas_table = gwas_table,
                                rank = rank,
                                nuc_range = nuc_range,
                                SNPmatrix = SNPmatrix)
    plot_data <- ld_tab %>%
      as.data.frame() %>%
      magrittr::set_colnames("value") %>%
      tibble::rownames_to_column("name")  %>%
      dplyr::mutate(pos = as.numeric(stringr::str_split_fixed(name, "[:|_]",3)[,2])) %>%
      dplyr::left_join(., impacts, by = c("pos")) %>%
      dplyr::mutate(layer_var = effect_impact) %>%
      dplyr::filter(!is.na(layer_var)) %>%
      dplyr::filter(value >= ld_cutoff)

    if(data_only){
      return(plot_data)
    }
    ## Build Plot

    p <- ggplot(data = plot_data) +

      # Line denoting SNP of interest
      geom_vline(aes(xintercept = snp_pos), color = "darkred") +
      # SNPs
      geom_point(aes(x= pos, y = 0, shape = effect_impact, color = value), size = 2.5, alpha = 0.8,
                 data = plot_data) +

      # color scale
      scale_color_viridis_c(option = "plasma", direction = -1) +
      scale_fill_viridis_c(option = "plasma", direction = -1) +
      # Gene Arrows
      gggenes::geom_gene_arrow(aes(xmin = start,
                                   xmax = end,
                                   y = molecule,
                                   #fill = gene,
                                   forward = direction), data = gene_labels) +
      # Gene labels
      gggenes::geom_gene_label(aes(xmin = start,
                                   xmax = end,
                                   y = molecule,
                                   #fill = gene,
                                   label = gene), data = gene_labels) + # Gene labels
      # Minimal Theme
      theme_minimal() +
      # ylim(-0.1,0.1) +
      facet_grid(rows = vars(
        #forcats::fct_relevel("HIGH", "MODERATE", "LOW", "MODIFIER", "LD", "ORF"))
        ordered(layer_var, levels = c("HIGH", "MODERATE", "LOW","MODIFIER","ORF"))),
        switch = "y") +
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
    if(add_ld_band){
      # Tile Geom for LD values
      p <- p +
        geom_tile(aes(x = pos , fill = value, color = value , y = 0, height = 2),
                data = plot_data %>%
                  dplyr::select(pos,value) %>%
                  dplyr::group_by(pos,value) %>%
                  dplyr::distinct() %>%
                  dplyr::mutate(layer_var = "LD")) +
        facet_grid(rows = vars(
          #forcats::fct_relevel("HIGH", "MODERATE", "LOW", "MODIFIER", "LD", "ORF"))
          ordered(layer_var, levels = c("HIGH", "MODERATE", "LOW","MODIFIER", "LD","ORF"))),
          switch = "y")
    }
    if(!ld_legend){
      p <- p +
        guides(color = FALSE,
               shape = FALSE,
               fill = FALSE) # Color guide
    } else{
      p <- p +
        guides(color = FALSE,
               shape = FALSE,
               fill = guide_colorbar(title = "Linkage")) # Color guide
    }
    return(p)
  }

}

#' Calculate linkage around SNP of interest and obtain phenotypes from phenotype table.
#' @details This function calculates LD in a specific region.
#' @details The SNPs per accession in the region are retrieved from polymorph or from a provided SNPmatrix.
#' @details For all SNPs that pass some criteria (see params), the phenotypes are then retrieved and split by presence of each SNP (SNP = TRUE or FALSE).
#' @param phenotype_table A table containing phenotypes that should be retrieved. Column with accession IDs has to be name 'ACC_ID'
#' @param phenotype Name of the phenotype of interest
#' @param chrom Chromosome where the SNP of interest is
#' @param pos Position of the SNP of interest on the chromosome
#' @param nuc_range Range of nucleotides that will be analyzed (total, split evenly up and downstream of the SNP)
#' @param SNPmatrix default "~/SNPmatrix.fst". If Needs to provided as a path to an fst file. Will be used for linkage calculation instead of polymorph. To use polymorph set to NULL
#' @param ld_stats The LD statistics, see SNPStats::ld. Default is LLR
#' @param ld_cutoff Only SNPs that have an LD values >= this will be plotted (default 10)
#' @param specific (optional) treatment column that was used to split samples for specific GWAS.
#' @param impacts SNP impacts to subset for can ("MODIFIER","LOW","MODERATE","HIGH"; default: c("MODERATE","HIGH"))

linkage_phenotypes <- function(phenotype_table,
                               phenotype,
                               chrom,
                               pos,
                               nuc_range,
                               ld_stats = "LLR",
                               ld_cutoff = 10,
                               SNPmatrix = "~/SNPmatrix.fst",
                               specific = c("bx","dag"),
                               impacts =  c("MODERATE","HIGH")){

  snps_in_range <- data.frame("chrom" = chrom, "pos" = pos, log10_p = 10) %>%
    gwaR::snp_linkage(anchored = T,
                plot = T,
                rank = 1,
                SNPmatrix = "~/SNPmatrix.fst",
                ld_stats =  "LLR",
                nuc_range = nuc_range,
                ld_legend = T,
                ld_cutoff = ld_cutoff,
                data_only = T) %>%
    dplyr::filter(effect_impact %in% impacts) %>%
    dplyr::select(chr, pos, effect_impact, type) %>%
    unique()

  message(paste0("Retrieving ", phenotype, " for ", nrow(snps_in_range), " linked SNPs (", ld_stats, " > ", ld_cutoff,")"))

  ### Retreive phenotypes

  pheno_by_SNP <- data.frame()
  for(i in 1:nrow(snps_in_range)){
    currsnp <- snps_in_range$pos[i]
    currtype <- snps_in_range$type[i]
    #cat(paste(i,"\n"))
    if(tryCatch(
      gwaR::phenotype_by_snp(phenotype_table,
                             phenotype = phenotype,
                             gwas_table = data.frame(chrom = chrom, pos = currsnp , log10_p = 10),
                             specific = specific,
                             SNPrank = 1,
                             plot = F,
                             SNPmatrix = SNPmatrix) ,
      error = function(e)
        return(TRUE) ) == T ){
      cat(paste("Skipping SNP",i, "at", currsnp ,"\n"))
    } else{
      tmp_df <- gwaR::phenotype_by_snp(phenotype_table,
                                         phenotype = phenotype,
                                         gwas_table = data.frame(chrom = chrom, pos = currsnp , log10_p = 10),
                                         specific = specific,
                                         SNPrank = 1,
                                         plot = F,
                                         SNPmatrix = SNPmatrix)  %>%
        dplyr::mutate(pos = currsnp,
               type = currtype)
      pheno_by_SNP <- rbind(pheno_by_SNP, tmp_df)
    }
  }
  ### Add original SNP back in
  pheno_by_SNP <- rbind(gwaR::phenotype_by_snp(phenotype_table,
                                               phenotype = phenotype,
                                               gwas_table = data.frame(chrom = chrom, pos = pos , log10_p = 10),
                                               specific = specific,
                                               SNPrank = 1,
                                               plot = F,
                                               SNPmatrix = SNPmatrix) %>%
                          dplyr::mutate(pos = pos,
                                 type = "Hit"),
                        pheno_by_SNP)
  pheno_by_SNP
}

#' Get information from thalemine.
#' @details This retrieves information from the InterMineR API. Takes a geneID, returns publications, functional annotations and GO Term
#' @param GeneID The GeneID of interest

mine_gene <- function(GeneID) {
  # Define Mine
  thalemine <- InterMineR::initInterMine(mine=InterMineR::listMines()["ThaleMine"])

  # First Query
  pubs <- InterMineR::getTemplateQuery(thalemine, "gene_publications")
  pubs$where[[1]][["value"]] <- GeneID
  ## Run
  pubs_info  <- InterMineR::runQuery(thalemine, pubs) %>% tibble::as_tibble()
  if(nrow(pubs_info > 0)){
    pubs_info <- pubs_info %>%
      dplyr::select("Gene.primaryIdentifier",
                    "Gene.symbol",
                    "Gene.publications.title", "Gene.publications.pubMedId") %>%
      magrittr::set_colnames(c("Identifier", "Symbol", "Title", "PMID"))
  }
  # Second Query
  functional <- InterMineR::getTemplateQuery(thalemine, "gene_generif")
  functional$where[[1]][["value"]] <- GeneID
  ## Run
  functionals <- InterMineR::runQuery(thalemine, functional) %>%
    tibble::as_tibble()
  if(nrow(functionals > 0)) {
    functionals <- functionals %>% dplyr::select("Gene.primaryIdentifier",
                                                 "Gene.symbol",
                                                 "Gene.geneRifs.annotation", "Gene.geneRifs.publication.title", "Gene.geneRifs.publication.pubMedId") %>%
      magrittr::set_colnames(c("Identifier", "Symbol", "Annotation", "Title", "PMID"))
  }

  # Third
  go_key <- InterMineR::getTemplateQuery(thalemine, "Keyword_GO_genes")
  go_key$where[[1]][["value"]] <- GeneID
  go_info  <- InterMineR::runQuery(thalemine, go_key) %>%
    tibble::as_tibble()

  # Return
  return(list("functions" = functionals,
              "publications" = pubs_info,
              "GO-Term Keywords" = go_info))

}
