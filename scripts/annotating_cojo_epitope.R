
# libraries
library(tidyverse)
library(data.table)
library(future) # for parallelization
library(furrr)  # for parallel 

#----------#
# inputs (locus breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_vep_extract <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/"
path_cojo <- "16-Dec-24_collected_independent_snps.csv"

# outputs
path_cojo_epitop <- paste0("/scratch/dariush.ghasemi/projects/pqtl_annotation/cojo_epitope_symbol_matching.tsv")
out_lb_epitop_cojo <- paste0(path_freez, "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_epitope_symbol_matching.tsv")

#----------#
# first extract the files inside the zip, 
# then take the list of them, and match seq_loci 
# to read corresponding file

#unzip(path_vep, exdir = path_vep_extract)
files_annot <- list.files(paste0(path_vep_extract, "snps_ld_in_meta_annot"), full.names = TRUE)


files_split <- files_annot %>%
  as_tibble() %>%
  # extract lead SNP, locus and seqid from annotation file name
  dplyr::mutate(
    txtpath = value,
    txtname = basename(txtpath) %>% str_remove("_snps_.*+"),
    snp     = txtname %>% str_remove_all("_seq.*+") %>% str_replace_all("_", ":"),
    seqid   = txtname %>% str_remove_all("_chr.*+") %>% str_remove("^\\d+_\\d+_[ATCG]+_[ATCG]+_"),
    locus   = txtname %>% str_remove_all(".+_chr")
  ) %>%
  dplyr::select(seqid, locus, snp, txtpath)


#----------#
# read LB file
lb_cistrans <- data.table::fread(paste0(path_freez, path_lb_cistrans))
cojo <- data.table::fread(paste0(path_freez, path_cojo))


cojo_annot <- cojo %>%
  dplyr::mutate(
    locus = str_remove_all(locus, "chr"),
  ) %>%
  left_join(
    files_split,
    join_by(study_id == seqid, locus, SNP == snp)
  ) %>%
  left_join(
    lb_cistrans %>% mutate(
      locus = str_c(chr, start, end, sep = "_"),
      Ensemble_noisoform = str_remove_all(Ensembl, "\\.\\d+") # remove isoform number from Ensembl_id
      ) %>%
      dplyr::select(phenotype_id, locus, cis_or_trans, symbol, Ensemble_noisoform),
    join_by(study_id == phenotype_id, locus)
    )


#----------------------------------------#
#----    Epitope for seqid and all   ----
#----------------------------------------#

# Steps to Implement:
# Append path to annotation file (from Michele F) to COJO results.
# Next, loop through each row, read the annotation file.
# Depending on the conditions below, define these columns:
#    1. epitope_inclusive: if any COJO or proxy variants has moderate- or 
#       high-impact consequences on any protein coding genes.
#    2. genes_with_epitope: All unique gene symbols where epitope was observed.
#    3. epitope: for cis only, if any COJO or proxy variants has moderate- or 
#       high-impact consequences on the gene encoding protein seqid (with gene symbol matching). 
#    4. epitope_snp: causal variants causing epitope effect.
#    5. epitope_high: for trans only, if any variant in annotated file has 
#       high-impact consequences on any coding genes.


mid_impact <- c(
  "inframe_insertion",
  "inframe_deletion",
  "missense_variant", 
  "missense_variant,splice_region_variant",
  "protein_altering_variant"
)

high_impact <- c(
  "transcript_ablation",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "stop_gained",
  "frameshift_variant",
  "stop_lost",
  "start_lost",
  "transcript_amplification",
  "feature_elongation",
  "feature_truncation"
)

epitope_consequences <- c(mid_impact, high_impact)

# Find epitope effect for COJO variants
find_epitope <- function(symbol, cis_or_trans, txtpath) {
  
  # If path is missing or empty or even not exists, 
  # return NA for all output columns
  # in case, gene id is empty, return NA
  if (
    is.na(txtpath) || txtpath == "" || 
    !file.exists(txtpath) || 
    (cis_or_trans == "cis" && symbol == "")
  ) {
    return(
      data.frame(
        txtpath = txtpath,
        epitope = NA,
        epitope_snp = NA,
        epitope_inclusive = NA,
        epitope_high = NA,
        genes_with_epitope = NA,
        genes_matched = NA
      ))
  }
  
  # Read the annotated file
  annot_df <- data.table::fread(txtpath)
  
  # filter for epitope consequences
  epitope_rows <- annot_df[Consequence %in% epitope_consequences]
  
  # filter for high-impact consequences on any protein coding genes
  high_impact_rows <- epitope_rows %>% 
    dplyr::filter(
      Consequence %in% high_impact,
      BIOTYPE == "protein_coding"
    )
  
  # take only protein cloding genes
  epitope_inclusive_rows <- epitope_rows[BIOTYPE == "protein_coding"]
  
  # Default output columns
  genes_matched <- NA
  epitope_high  <- NA
  
  # Only match gene symbol when locus/cojo is cis
  if (cis_or_trans == "cis") {
    
    # separte gene symbols if multiple symbols exist for protein target
    gene_names <- strsplit(symbol, "\\|\\s*")[[1]]
    # filter for genes matching with protein target
    gene_rows  <- epitope_rows %>% dplyr::filter(SYMBOL %in% gene_names)
    
    # define epitope artifact
    epitope <- ifelse(nrow(gene_rows) > 0, "Yes", "No")
    
    # report snps with middle/high impact on protein target gene
    epitope_snp <- ifelse(
      nrow(gene_rows) > 0,
      paste(unique(gene_rows$SNPID), collapse = ";"),
      NA
    )
    
    # flag variable to indicate rows with mismatched genes
    matched_rows <- annot_df[symbol %in% gene_names]
    genes_matched  <- ifelse(nrow(matched_rows) > 0, "Yes", "No")
    
  } else {
    
    # now define epitope for trans loci
    epitope <- NA
    epitope_snp <- NA
    
    # defining epitope only using high-impact annotations
    epitope_high <- ifelse(nrow(high_impact_rows) > 0, "Yes", "No")
    
  }
  
  
  # Find epitope effect for any protein coding genes at locus
  epitope_inclusive <- ifelse(nrow(epitope_inclusive_rows) > 0, "Yes", "No")
  
  # Report multiple affected genes in a single row
  genes_with_epitope <- ifelse(
    nrow(epitope_rows) > 0,
    paste(unique(epitope_rows$SYMBOL), collapse = ", "),
    NA
  )
  
  # shape final output
  return(
    data.frame(
      txtpath,
      epitope,
      epitope_snp,
      epitope_inclusive,
      epitope_high,
      genes_with_epitope,
      genes_matched
    )
  )
  
}


# find number of available cores
parallel::detectCores()
# set cores number for parallel analysis
future::plan(multicore, workers = 32) # default is parallelly::availableCores()



# Apply function to each row using pmap in purrr
# which allows named arguments and avoids atomic vector issues
results_epitope <- pmap_dfr(
  cojo_annot %>% dplyr::select(symbol, cis_or_trans, txtpath),
  find_epitope
  )

# Combine results with COJO
cojo_annot_epitop <- cojo_annot %>%
  inner_join(results_epitope, join_by(txtpath)) %>%
  dplyr::select(- c(txtpath, symbol, Ensemble_noisoform))

# save COJO file with epitope effect
data.table::fwrite(
  cojo_annot_epitop, 
  file = path_cojo_epitop, 
  quote = F, 
  row.names = F, 
  sep = "\t"
  )


#----------------------------------------#
#----    Epitope for COJO variants   ----
#----------------------------------------#

# prepare combined results for join with LB
cojo_annot_epitop_4join <- cojo_annot_epitop %>%
  dplyr::select(
    study_id,
    locus,
    SNP,
    epitope,
    epitope_snp,
    epitope_inclusive,
    genes_with_epitope,
    epitope_high
    ) %>%
  group_by(study_id, locus) %>%
  summarise(
    cis_epitope_effect = all(epitope == "Yes"),
    epitope_cojo_any = any(epitope == "Yes"),
    epitope_cojo_high = any(epitope_high == "Yes"),
    impact_variant_vep = all(epitope_inclusive == "Yes"),
    total_cojo_snps = n(),
    epitope_yes = sum(epitope == "Yes"),
    epitope_status = paste0(epitope_yes, "of", total_cojo_snps),
    prop_epitope_yes = round(epitope_yes / total_cojo_snps, 2),
    epitope_causing_cojo = paste(unique(epitope_snp), collapse = "; "),
    genes_with_epitope = paste(unique(genes_with_epitope), collapse = ";")
  ) %>%
  ungroup()

#----------#
# Combine results with LB
lb_epitope_cojo <- lb_cistrans %>%
  dplyr::mutate(locus = str_c(chr, start, end, sep = "_")) %>%
  left_join(
    cojo_annot_epitop_4join,
    join_by(phenotype_id == study_id, locus)
  )

#----------#
# save LB file with epitope effect for COJO SNPs
data.table::fwrite(
  lb_epitope_cojo,
  file = out_lb_epitop_cojo,
  quote = F,
  row.names = F,
  sep = "\t"
)

