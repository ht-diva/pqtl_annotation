
# libraries
library(tidyverse)
library(data.table)

#----------#

# inputs (locus_breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "Archived/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_vep_extract <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/"
path_cojo <- "16-Dec-24_collected_independent_snps.csv"

# outputs
path_st3 <- paste0(path_freez, "17-Jul-25_suppl_table_3_cojo_variants_vep_annotated.tsv")

#----------#

# To add annotated gene name mapping to lead SNP of the locus in LB results:
#    1. List annotated files extracted before from zip file sent by Michele F
#    2. 
#    3. Extract columns needed for join with LB
#    4. Append annotated file name to LB file
#    5. Write a function to find closest gene name encoding a protein
#    6. Write a function reading annotated file, filter for lead SNP and find gene name
#    7. Add annotated gene name by iterating the functions in row-wise manner

#----------#
# read data files
lb_cistrans <- fread(paste0(path_freez, path_lb_cistrans))
cojo <- read.csv(paste0(path_freez, path_cojo))


# 1. List the annotated files in zip
files_annot <- list.files(paste0(path_vep_extract, "snps_ld_in_meta_annot"), full.names = TRUE)

# 3. Extract columns needed for join with LB
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
  dplyr::select(seqid, locus, snp, txtpath) # take required column for joining LB file


# 4. Append annotated file name to COJO
cojo_annot <- cojo %>%
  dplyr::mutate(
    locus = str_remove_all(locus, "chr")
  ) %>%
  # since there are more than one annotated file for some loci, join by 'SNPID == snp'
  # ensures that we only take annotation file corresponding to the locus lead variant
  # among the list of annotation files belonged to all independent snps at a locus
  left_join(
    files_split,
    join_by(study_id == seqid, locus, SNP == snp)
  ) %>%
  left_join(
    lb_cistrans %>% mutate(
      locus = str_c(chr, start, end, sep = "_"), # rename or reshape ID columns for join
      Ensemble_noisoform = str_remove_all(Ensembl, "\\.\\d+") # remove isoform number from Ensembl_id
    ) %>%
      dplyr::select(phenotype_id, locus, cis_or_trans, Ensemble_noisoform),
    join_by(study_id == phenotype_id, locus)
  )


#----------#

# 5. Function to find annotation of COJO variants
find_annotation <- function(Ensemble_noisoform, cis_or_trans, cojo_snp, txtpath) {
  
  # If path is missing or empty or even not exists, returns NA
  # in case, gene id is empty, returns NA
  if (is.na(txtpath) || txtpath == "" || !file.exists(txtpath) || 
      (cis_or_trans == "cis" && Ensemble_noisoform == "")) {
    return(data.table(consequence = NA_character_, affected_gene = NA_character_))
  }
  
  # Read the annotated file
  annot_dt <- data.table::fread(txtpath)
  
  # Subset to COJO variant
  cojo_rows <- annot_dt[SNP == cojo_snp]
  
  # Default values for consequence (gene-specific check)
  #consequence <- "No"
  
  get_result <- function(rows) {
    data.frame(
      consequence   = paste(unique(rows$Consequence), collapse = "; "),
      affected_gene = paste(unique(rows$SYMBOL), collapse = "; ")
    )
  }
  
  # Only check Ensembl_id match if "cis"
  if (cis_or_trans == "cis") {
    # see if Gene ID matches with one of the ids in Ensembl_id for the seqid
    gene_ids  <- strsplit(Ensemble_noisoform, ";\\s*")[[1]]
    gene_rows <- cojo_rows %>% dplyr::filter(Gene %in% gene_ids)
    
    return(if (nrow(gene_rows) > 0) get_result(gene_rows) else get_result(cojo_rows))

    # Now COJO variants in trans loci
    } else {
      coding_rows    <- cojo_rows[BIOTYPE   == "protein_coding"]
      canonical_rows <- cojo_rows[CANONICAL == "YES"]
      
      if (nrow(coding_rows)    > 0) return(get_result(coding_rows))
      if (nrow(canonical_rows) > 0) return(get_result(canonical_rows))
    
    return(get_result(cojo_rows))
  }
}

#----------#
# Supplementary Table 3

# 6. iterate function to read annotated file, subset data to COJO SNP,
# and find annotation in row-wise manner
results_annot <- pmap(
  cojo_annot %>% dplyr::select(Ensemble_noisoform, cis_or_trans, cojo_snp, txtpath),
  find_annotation
)

# Combine results with COJO
st3 <- cojo_annot %>%
  dplyr::select(- c(txtpath)) %>%
  cbind(results_annot)

#----------#
# save output
data.table::fwrite(
  st3,
  file = path_st3,
  sep = "\t",
  quote = F,
  row.names = F
  )


lb_annot_missing <- lb_annot %>% 
  filter(!complete.cases(.)) %>%
  dplyr::rename(Chr = chr, bp = POS, SNP = SNPID, study_id = phenotype_id) %>%
  dplyr::select(Chr, bp, SNP, locus, study_id)


write.csv(lb_annot_missing, "/scratch/dariush.ghasemi/projects/SNPS_IN_LD/missing_loci.csv", quote = F, row.names = F)

