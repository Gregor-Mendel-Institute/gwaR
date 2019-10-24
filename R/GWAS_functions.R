# A collection of functions for working with GWAS tables
# Niklas Schandry, Patrick Hüther 2019


require(tidyverse)
require(ggthemes)
require(ggrepel)
require(feather)
require(biomaRt)
require(plyranges)
require(leaflet)

#This function has no real purpose.

format_gwas <- function(gwas){
  gwas %<>% mutate(chrom = as.double( str_extract(chrom, "[0-9]")),
                   log10_p = score,
  )
  return(gwas_object)
}

#This function reads GWAS output from limix.

read_gwas <- function(gwas_path){
  gwas_object <- vroom::vroom(gwas_path)
  fdr_corr <- qvalue::qvalue(p = gwas_object$pv)
  fdr_thresh <- cbind(fdr_corr$pvalues[which(fdr_corr$qvalues < 0.05 )]) %>% max() %>% log10() %>% abs()
  bf_corr <- (0.05/nrow(gwas_object)) %>% log10() %>% abs()
  gwas_object %<>% mutate(chrom = as.double( str_extract(chrom, "[0-9]")),
                          Significant = case_when(-log10(pv) > bf_corr ~ "Bonferroni",
                                                  -log10(pv) > fdr_thresh ~ "FDR",
                                                  TRUE ~ "Not"),
                          log10_p = -log10(pv),
                          fdr_thresh = fdr_thresh,
                          bf_corr = bf_corr,
                          model = case_when(grepl(pattern = "_specific_", gwas_path) ~ "specific",
                                            grepl(pattern = "_common_", gwas_path) ~ "common",
                                            grepl(pattern = "_any_", gwas_path) ~ "any",
                                            TRUE ~ "Unknown"
                          )
  )
  return(gwas_object)
}

# This function takes a GWAS table (from read_gwas) and queries 1001genomes to obtain the
# IDs of all accessions that carry the SNP ranked at SNPrank.

get_polymorph_acc <- function(gwas_table, SNPrank){
  return(gwas_table %>%
           arrange(desc(log10_p)) %>%
           dplyr::slice(SNPrank) %>%
           dplyr::select(chrom, pos) %>%
           {httr::content(httr::GET(paste0("http://tools.1001genomes.org/api/v1.1/variants.csv?type=snps;accs=all;chr=",
                                           .$chrom,
                                           ";pos=",
                                           .$pos)), col_types = "iicccccc") }
  )
}

#This function plots an interative map of accessions that carry the SNP of interest;
#needs GWAS table and the rank of the SNP of interest

plot_acc_map <- function(gwas_table, SNPrank){
  allAccessions <- readr::read_csv("~/labshare/lab/accessions/1001genomes-accessions.csv")
  allAccessions %>%
    filter(id %in% get_polymorph_acc(gwas_table, SNPrank)$strain) %>%
    leaflet(data=.) %>%
    addTiles() %>%
    addCircleMarkers(lng=~longitude,
                     lat=~latitude,
                     label=~as.character(id),
                     popup=~as.character(name),
                     stroke=FALSE, radius=5, fillOpacity=0.8, color="#007243")
}

# Get expression data for a specific GeneID (format: ATXGNNNNN)

get_expression <- function(GeneID){
  if(!exists("Kawakatsu_dat")){
    Kawakatsu_dat <- read_csv("https://arapheno.1001genomes.org/static/rnaseq/Epigenomic_Diversity_in_A._thaliana_(Kawakatsu_et_al._2016).csv") %>%
      as_tibble() %>% set_colnames(c("ACC_ID", colnames(.)[2:ncol(.)]))
  }
  if(!paste(GeneID) %in% colnames(Kawakatsu_dat)){
    stop(paste0("No expression data in Kawakatsu2016 for GeneID: ",GeneID))
  }
  Kawakatsu_dat %>% dplyr::select(ACC_ID, eval(GeneID)) %>%
    pivot_longer(starts_with("AT"), names_to = "GeneID", values_to = "phenotype_value")
}

# Based on a GeneID and a SNPtable, returns a table of expression values for that gene,
# where accessions that contain SNP have TRUE in hasSNP

intersect_expression_snp <- function(GeneID, gwas_table, SNPrank){
  get_expression(GeneID) %>%
    mutate(hasSNP = case_when(ACC_ID %in% get_polymorph_acc(gwas_table,SNPrank)$strain ~ TRUE,
                              TRUE ~ FALSE))
}

# Better than the above, simply works with a GWAS table and the rank, makes use of get_nearest_genes to find closest gene.

retrieve_counts <- function(gwas_table, SNPrank){
  genes <-  get_nearest_genes(gwas_table, SNPrank) %>%
    slice(SNPrank) %>%
    .$GeneId
  get_expression(genes) %>%
    mutate(hasSNP = case_when(ACC_ID %in% get_polymorph_acc(gwas_table, SNPrank)$strain ~ TRUE,
                              TRUE ~ FALSE))
}

# Convenience quick plot, using data from retrieve_counts

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


# Finding nearest genes for GWAS hits

# This function takes a table as it comes out of read_gwas and extracts the gene annotations that match the top n_hit SNPs
# from the arabidopsis genome.
# This function will assign araGenes in the Global Environment!

get_nearest_genes <- function(GWAS = NULL, n_hit = 1){
  if(is.null(GWAS)){
    stop("GWAS output missing")
  }
  if(!require(plyranges)) {
    stop("Please install plyranges")
  }
  if(!require(biomaRt)) {
    stop("Please install biomaRt")
  }
  if(!require(tidyverse)) {
    stop("Please install tidyverse")
  }
  require(plyranges)
  require(biomaRt)
  require(tidyverse)
  # get annotation info from plants_mart
  if(!exists("araGenes")){
    ara <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
    # Get genome, rename columns, make into granges
    araGenes <<- getBM(c("ensembl_gene_id",
                         "chromosome_name",
                         "start_position",
                         "end_position",
                         "strand",
                         "description",
                         "transcript_biotype"),
                       mart=ara) %>%
      dplyr::rename(GeneId=ensembl_gene_id) %>%
      dplyr::mutate(strand=case_when(strand == 1 ~ "+",
                                     strand == -1 ~ "-"),
                    start = start_position,
                    end = end_position) %>%
      dplyr::filter(str_detect(chromosome_name,"[1-9]")) %>%
      dplyr::mutate(seqnames=paste0("Chr",chromosome_name)) %>%
      as_granges()
  }


  snp <- GWAS %>%
    dplyr::arrange(desc(-log10(pv))) %>%
    dplyr::slice(1:n_hit) %>%
    dplyr::mutate(seqnames = paste0("Chr", chrom),
                  start = pos,
                  end = pos) %>%
    tibble::rownames_to_column("SNP_rank") %>%
    as_granges()

  plyranges::join_nearest(snp, araGenes) %>%
    as_tibble() %>%
    dplyr::select(SNP_rank, chrom, pos, Significant, log10_p, mac, GeneId, description, transcript_biotype, start_position, end_position#, .drop_ranges = TRUE
    ) %>%
    mutate(description = case_when(
      description == "" ~ "N/A",
      TRUE ~ description
    ))


  #araGenes %>% filter_by_overlaps(snps, maxgap = distance)
}

# Find overlapping genes

# This function takes a table as it comes out of read_gwas and extracts the gene annotations that overlap with
# the top n_hit SNPs from the arabidopsis genome. Can be used for nearest matching by changing the distance parameter.
# This function will assign araGenes in the Global Environment!

get_overlapping_genes <- function(GWAS = NULL, n_hit = 1, distance = -1){
  if(is.null(GWAS)){
    stop("GWAS output file missing")
  }
  if(!require(plyranges)) {
    stop("Please install plyranges")
  }
  if(!require(biomaRt)) {
    stop("Please install bioMart")
  }
  if(!require(tidyverse)) {
    stop("Please install tidyverse")
  }
  require(plyranges)
  require(biomaRt)
  require(tidyverse)
  # get annotation info from plants_mart
  if(!exists("araGenes")){
    ara <- biomaRt::useMart(biomart = "plants_mart", dataset = "athaliana_eg_gene", host = 'plants.ensembl.org')
    # Get genome, rename columns, make into granges
    araGenes <<- getBM(c("ensembl_gene_id",
                         "chromosome_name",
                         "start_position",
                         "end_position",
                         "strand",
                         "description",
                         "transcript_biotype"),
                       mart=ara) %>%
      dplyr::rename(GeneId=ensembl_gene_id) %>%
      dplyr::mutate(strand=case_when(strand == 1 ~ "+",
                                     strand == -1 ~ "-"),
                    start = start_position,
                    end = end_position) %>%
      dplyr::filter(str_detect(chromosome_name,"[1-9]")) %>%
      dplyr::mutate(seqnames=paste0("Chr",chromosome_name)) %>%
      as_granges()
  }


  snp <- GWAS %>%
    dplyr::arrange(desc(-log10(pv))) %>%
    dplyr::slice(1:n_hit) %>%
    dplyr::mutate(seqnames = paste0("Chr", chrom),
                  start = pos,
                  end = pos,
                  max_distance = distance) %>%
    tibble::rownames_to_column("SNP_rank") %>%
    as_granges()

  join_overlap_inner(snp,araGenes, maxgap = distance) %>%
    plyranges::select(SNP_rank, chrom, pos, Significant, log10_p, mac, GeneId, max_distance, description, transcript_biotype, start_position, end_position, .drop_ranges = TRUE) %>%
    as_tibble()

  #araGenes %>% filter_by_overlaps(snps, maxgap = distance)
}

# Plot GWAS
## Convenience function to plot manhattan plots, via ggplot2, using GMI branding

plot_gwas <- function(x, title = "No Title", subtitle = NULL){
  color_gmi_light <- ("#abd976")
  color_gmi_dark <- ("#007243")
  GWAS_colors <- c(color_gmi_dark, color_gmi_light, "grey50")
  names(GWAS_colors) <- c("Bonferroni", "FDR", "Not")
  ggplot(aes(x=pos, y=log10_p), data = x) +
    # geom_hline(linetype = "dotted", yintercept = bf_corr) +
    facet_grid(~chrom, scales = "free_x", switch = "x") +
    #  geom_label_repel(aes(x=pos, y= -log10(pv),label = gene), data = ft_specific_limix_mac5 %>% filter(gene != "NA")) +
    theme_few() +
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

# Plot manhattan plot with annotations
## Defaults to only plotting overlapping annotations, can be toggled by "match_nearest = TRUE"


plot_annotated_gwas <- function(gwas, title ="No Title", subtitle = NULL, nlabels = 5, labeltype = "GeneId", match_nearest = FALSE) {
  if(!match_nearest){
    annotations <- gwas %>% get_overlapping_genes(nlabels) %>%  unite(labs, SNP_rank, labeltype, sep = " :")
  } else {
    annotations <- gwas %>% get_nearest_genes(nlabels) %>%  unite(labs, SNP_rank, labeltype, sep = " :")
  }

  color_gmi_light <- ("#abd976")
  color_gmi_dark <- ("#007243")
  GWAS_colors <- c(color_gmi_dark, color_gmi_light, "grey50")
  names(GWAS_colors) <- c("Bonferroni", "FDR", "Not")
  #Step1: Bind those tables together
  gwas %>%
    filter(abs(log10_p) > 0.5) %>%
    #Step 2: Plot
    ggplot(aes(x=pos, y=log10_p)) +
    geom_point(aes(color = Significant)) +
    geom_text_repel(aes(label = labs), data = annotations) +
    facet_grid(~chrom, scales = "free_x", switch = "x") +
    theme_few() +
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
