
# libraries
library(dplyr)
library(data.table)
library(purrr)

#----------#

# inputs (locus_breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "Archived/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_vep_extract <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/"
path_cojo <- "16-Dec-24_collected_independent_snps.csv"

# outputs
path_st3 <- paste0(path_freez, "05-Aug-25_suppl_table_3_cojo_variants_vep_annotated.tsv")

#----------#

# To find genes impacted by middle/high impact COJO or proxy variants:
#    1. List annotated files extracted before from Michele F's zip file
#    2. Extract locus to align ST3 with loci list ST2
#    3. Append path to annotation files to COJO
#    4. Define a function to find impacted genes and consequences
#    5. Add annotations by iterating the function row-wise

#----------#
# read data files
lb_cistrans <- fread(paste0(path_freez, path_lb_cistrans))
cojo <- read.csv(paste0(path_freez, path_cojo))


# 1. List the annotated files in zip
files_annot <- list.files(paste0(path_vep_extract, "snps_ld_in_meta_annot"), full.names = TRUE)

# 2. Extract locus to align ST3 with loci list ST2
files_split <- files_annot %>%
  as_tibble() %>%
  # extract COJO SNP, locus, and seqid from filename
  mutate(
    txtpath = value,  # easy rename; as_tibble used 'value' as the column name
    txtname = basename(txtpath) %>% str_remove("_snps_.*+"), # get rid of full path, then remove target variant from file name
    snp     = txtname %>% str_remove_all("_seq.*+") %>% str_replace_all("_", ":"),
    seqid   = txtname %>% str_remove_all("^\\d+_\\d+_[ATCG]+_[ATCG]+_") %>% str_remove_all("_(chr\\d+_\\d+_\\d+)|_(\\d+_\\d+_\\d+)"),
    locus   = txtname %>% str_remove_all("^\\d+_\\d+_[ATCG]+_[ATCG]+_(seq.\\d+.\\d+)_") %>% str_remove_all("chr")
  ) %>%
  dplyr::select(seqid, locus, snp, txtpath) # take required column for joining COJO file


# 3. Append path to annotation files to COJO
cojo_annot <- cojo %>%
  dplyr::mutate(
    locus = str_remove_all(locus, "chr")
  ) %>%
  # take annotation file corresponding to the locus and COJO variant
  left_join(
    files_split,
    join_by(study_id == seqid, locus, SNP == snp)
  ) %>%
  left_join(
    lb_cistrans %>% mutate(
      locus = str_c(chr, start, end, sep = "_") # define locus for join
    ) %>%
      dplyr::select(phenotype_id, locus, cis_or_trans),
    join_by(study_id == phenotype_id, locus)
  )


#----------#

# 4. Function to find find impacted genes and consequences
find_annotation <- function(txtpath) {
  
  # If path is missing or empty or even not exists, returns NA
  # in case, gene id is empty, returns NA
  if (is.na(txtpath) || txtpath == "" || !file.exists(txtpath)) {
    return(data.table(txtpath = txtpath, genes_with_consequence_types = NA_character_))
  }
  
  # Read the annotated file
  annot_dt <- data.table::fread(txtpath)
  
  # Subset to variants with mid/high impact on any protein-coding genes
  cojo_rows <- annot_dt %>%
    dplyr::filter(
      Consequence %in% epitope_consequences,
      BIOTYPE == "protein_coding"
    )
  
  # If any relevant rows found, collapse consequences per gene
  if (nrow(cojo_rows) > 0) {
    
    combined_impacts <- cojo_rows %>%
      group_by(SNP, SYMBOL) %>%
      summarise(
        vep_annot = paste(unique(Consequence), collapse = ","),
        .groups = "drop"
      )
    
    # column will be used in Suppl. Table 3
    genes_with_consequence_types <- combined_impacts %>%
      mutate(label = paste0(SYMBOL, ": ", vep_annot)) %>%
      pull(label) %>%
      unique() %>%
      paste(collapse = "; ")
    
  } else {
    # otherwise return NA
    genes_with_consequence_types <- NA_character_
  }
  
  return(
    data.frame(
      txtpath,
      genes_with_consequence_types,
      stringsAsFactors = FALSE
    )
  )
}

#----------#
# Supplementary Table 3

# 5. iterate function to read annotated file and find annotation row-wise
results_annot <- pmap_dfr(
  cojo_annot %>% dplyr::select(txtpath),
  find_annotation
)


# Combine results with COJO
st3 <- cojo_annot %>%
  right_join(results_annot, join_by(txtpath)) %>%
  dplyr::select(- txtpath)

#----------#
# save output
data.table::fwrite(
  st3,
  file = path_st3,
  sep = "\t",
  quote = F,
  row.names = F
  )


