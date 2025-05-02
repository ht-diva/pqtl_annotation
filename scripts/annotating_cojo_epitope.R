
# libraries
library(tidyverse)
library(data.table)

#----------#
# inputs (locus breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_vep.tsv"
path_vep_extract <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/"
path_cojo <- "16-Dec-24_collected_independent_snps.csv"

# outputs
path_lb_epitop <- paste0(path_freez, "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_vep_epitope_high_moderate.tsv")


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
  dplyr::select(seqid, locus, snp, txtpath) #%>% 
  #distinct(seqid, locus, .keep_all = T)


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
      dplyr::select(phenotype_id, locus, Ensemble_noisoform),
    join_by(study_id == phenotype_id, locus)
    )


#----------------------------------------#
#----    Epitope for seqid and all   ----
#----------------------------------------#

# Steps to Implement:
# Read the LB file which has path to each annotated TXT file.
# Loop through each row, read the corresponding annotated TXT file, and filter relevant rows.
# Check if the Consequence column contains any epitope effect's consequences.
# If there’s a match, store TRUE and the corresponding SNPID; otherwise, store FALSE and NA.
# In the end, we repot 4 new variables:
#    1. epitope_effect_all: Indicates if any variant in annotated file has an epitope effect, regardless of gene match.
#    2. genes_with_epitope_effects: Lists all unique gene symbols where an epitope effect was observed.
#    3. epitope_effect: indicating epitope effect with seqid gene-matching only for rows where cis_or_trans == "cis". 
#    4. epitope_causing_variant: causal variants in a single row (concatenated with ;) causing epitope effect.


epitope_consequences <- c(
  "missense_variant", 
  "missense_variant,splice_region_variant",
  "start_lost",
  "stop_lost",
  "stop_gained",
  "protein_altering_variant",
  "inframe_insertion",
  "inframe_deletion",
  "feature_truncation",
  "feature_elongation",
  "transcript_amplification",
  "frameshift_variant",
  "splice_donor_variant",
  "splice_acceptor_variant",
  "transcript_ablation"
)


# Function to find epitope effect of variants in locus breaker results
find_epitope <- function(Ensemble_noisoform, cis_or_trans, txtpath) {
  
  # Check if the file path is empty or not exists
  # Handles missing files gracefully
  # If path is missing, return NA for all output columns
  if (is.na(txtpath) || txtpath == "" || !file.exists(txtpath)) {
    return(
      data.frame(
        missing_ld = "Yes",
        epitope_effect = NA,
        epitope_causing_variant = NA,
        epitope_effect_all = NA,
        genes_with_epitope_effects = NA
      ))
    }
  
  # Read the annotated file
  annot_df <- data.table::fread(txtpath)
  
  # Check for any epitope effects in the whole file (epitope_effect_all)
  epitope_rows_all <- annot_df %>% dplyr::filter(Consequence %in% epitope_consequences)
  
  # Default values for epitope_effect (gene-specific check)
  epitope_effect <- "No"
  epitope_causing_variant <- NA
  
  # Only check Ensembl_id match if "cis"
  if (cis_or_trans == "cis") {
    epitope_rows <- annot_df %>%
      dplyr::filter(Gene == Ensemble_noisoform) %>% # take annotations where Gene ID matches with Ensembl_id of seqid
      dplyr::filter(Consequence %in% epitope_consequences) # Check if any variant has epitope consequences

    if (nrow(epitope_rows) > 0) {
      # report multiple causal variants in a single row
      epitope_effect <- "Yes"
      epitope_causing_variant <- paste(unique(epitope_rows$SNPID), collapse = ";")
    }
  }
  
  # Find if there is epitope effect for any genes at locus
  epitope_effect_all <- ifelse(nrow(epitope_rows_all) > 0, "Yes", "No")
  
  # Report multiple affected genes in a single row
  genes_with_epitope_effects <- ifelse(
    nrow(epitope_rows_all) > 0,
    paste(unique(epitope_rows_all$SYMBOL), collapse = ";"),
    NA
    )
  
  # flag here that LD path is not missing
  missing_ld <- "No"
  
  # shape final output
  res <- data.frame(
    missing_ld,
    epitope_effect_all,
    genes_with_epitope_effects,
    epitope_effect,
    epitope_causing_variant
  )
  
  return(res)
}


# Apply function to each row using pmap in purrr
# which allows named arguments and avoids atomic vector issues
results_epitope <- pmap_dfr(
  lb_cistrans_annot %>% dplyr::select(Ensemble_noisoform, cis_or_trans, txtpath),
  find_epitope
  )



# Combine results with LB
lb_cistrans_annot_epitop <- lb_cistrans_annot %>%
  dplyr::select(- c(locus, Ensemble_noisoform, txtpath)) %>%
  cbind(results_epitope)



# save Lb file with epitope effect
data.table::fwrite(lb_cistrans_annot_epitop, file = path_lb_epitop, quote = F, row.names = F, sep = "\t")
