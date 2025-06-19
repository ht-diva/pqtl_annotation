
# libraries
library(tidyverse)
library(data.table)

#----------#

# inputs (locus_breaker or LB results)
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_cistrans <- "Archived/mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann.csv"
path_vep <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/snps_ld_in_meta_annot.zip"
path_vep_extract <- "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/unzipped/"
path_vep_missing <- "/scratch/dariush.ghasemi/projects/snps_ld_in_meta_annot"

# outputs
path_lb_gene <- paste0(path_freez, "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_vep.tsv")

#----------#

# To add annotated gene name mapping to lead SNP of the locus in LB results:
#    1. Extract annotated files inside the zip, 
#    2. List the annotated files
#    3. Extract columns needed for join with LB
#    4. Append annotated file name to LB file
#    5. Write a function to find closest gene name encoding a protein
#    6. Write a function reading annotated file, filter for lead SNP and find gene name
#    7. Add annotated gene name by iterating the functions in row-wise manner

#----------#
# read data files
lb_cistrans <- data.table::fread(paste0(path_freez, path_lb_cistrans))


# 1. Extract zip (run it once)
#unzip(path_vep, exdir = path_vep_extract)


# 2. List the annotated files in zip
files_annot <- list.files(paste0(path_vep_extract, "snps_ld_in_meta_annot"), full.names = TRUE)

# annotations for the 55 missing loci
files_annot_missing <- list.files(path_vep_missing, full.names = TRUE) #%>% gsub("_(\\d+_\\d+)_", "_chr\\1_", .)


# 3. Extract columns needed for join with LB
files_split <- files_annot %>%
  c(files_annot_missing) %>%
  as_tibble() %>%
  # extract lead SNP, locus and seqid from annotation file name
  mutate(
    txtpath = value,  # easy rename; as_tibble used 'value' as the column name
    txtname = basename(txtpath) %>% str_remove("_snps_.*+"), # get rid of full path, then remove target variant from file name
    snp     = txtname %>% str_remove_all("_seq.*+") %>% str_replace_all("_", ":"),
    seqid   = txtname %>% str_remove_all("^\\d+_\\d+_[ATCG]+_[ATCG]+_") %>% str_remove_all("_(chr\\d+_\\d+_\\d+)|_(\\d+_\\d+_\\d+)"),
    locus   = txtname %>% str_remove_all("^\\d+_\\d+_[ATCG]+_[ATCG]+_(seq.\\d+.\\d+)_") %>% str_remove_all("chr")
  ) %>%
  dplyr::select(seqid, locus, snp, txtpath) # take required column for joining LB file


# 4. Append annotated file name to LB file
lb_annot <- lb_cistrans %>%
  # rename or reshape ID columns for join
  dplyr::mutate(
    locus = str_c(chr, start, end, sep = "_")
  ) %>%
  # since there are more than one annotated file for some loci, join by 'SNPID == snp'
  # ensures that we only take annotation file corresponding to the locus lead variant
  # among the list of annotation files belonged to all independent snps at a locus
  left_join(
    files_split,
    join_by(phenotype_id == seqid, locus, SNPID == snp)
  )


#----------#

# 5. Function to find name of the closest gene 
# encoding a protein for each row in LB results
grep_gene_lb <- function(df) {
  
  # Check if any protein coding gene are present
  if ("protein_coding" %in% df$BIOTYPE) {
    
    # filter out any lincRNA, ncRNAs, etc.
    df_prot <- df %>% dplyr::filter(BIOTYPE == "protein_coding")
    
    # Consider in case all distances were missing
    if(all(is.na(df_prot$DISTANCE))){
      # Report unique genes
      return(
        df_prot %>%
          summarise(annotated_genes = paste(unique(SYMBOL), collapse = "; ")) %>%
          pull(annotated_genes)
      )
    } else{
      # Report closest gene by sorting the distance
      return(
        df_prot %>%
          top_n(DISTANCE, n = -1) %>% # take the closest gene in vicinity of the lead variant
          # Collapses target genes into a string
          summarise(annotated_genes = paste(unique(SYMBOL), collapse = "; ")) %>%
          pull(annotated_genes)
        )
      }
    } else {
    # Check if any canonical transcripts are present
    if ("YES" %in% df$CANONICAL) {
      # Report canonical transcripts
      return(
        df %>% 
          dplyr::filter(CANONICAL == "YES") %>%
          # Collapses target genes into a string
          summarise(annotated_genes = paste(unique(SYMBOL), collapse = "; ")) %>% 
          pull(annotated_genes)
      ) 
    }
    else {
    # Return empty if no relevant genes found
    return(data.frame())  
  }
  }
}


# 6. Function to read annotated file, subset data 
# to lead SNP of the locus and find gene name
take_annot <- function(filepath, variant){
  data.table::fread(file = filepath) %>% 
    dplyr::filter(SNPID == variant) %>%
    grep_gene_lb(.)
}

#----------#

# 7. iterating the functions
lb_with_genes <- lb_annot %>%
  dplyr::filter(!is.na(txtpath)) %>% # in case there is a missing annotation for any locus
  dplyr::mutate(
    annot_gene_vep = map2_chr(txtpath, SNPID, take_annot) # iterate function to retrieve gene names
    ) %>%
  dplyr::select(- c(txtpath, locus)) # remove columns Solene suggested


#----------#
# save output
data.table::fwrite(
  lb_with_genes,
  file = path_lb_gene,
  sep = "\t",
  quote = F,
  row.names = F
  )


lb_annot_missing <- lb_annot %>% 
  filter(!complete.cases(.)) %>%
  dplyr::rename(Chr = chr, bp = POS, SNP = SNPID, study_id = phenotype_id) %>%
  dplyr::select(Chr, bp, SNP, locus, study_id)


write.csv(lb_annot_missing, "/scratch/dariush.ghasemi/projects/SNPS_IN_LD/missing_loci.csv", quote = F, row.names = F)

