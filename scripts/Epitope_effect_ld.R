locus_breker_file = "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_frozen_version_23_08_24/locus_breaker_cis_trans_uniprot_gene_rsid_hotspotted_closest_gene_annotated.csv"

x=read.delim(locus_breker_file, header=T, sep=",")
# get the single-SNP regions or the 2-SNP regions:
library(dplyr)
library(tidyr)
library(data.table)

# add VEP annotation:
vep_file = "/exchange/healthds/pQTL/pQTL_workplace/annotations/VEP/data/snps_ld_in_meta_annot.tar.gz"
library(archive)
library(readr)
read_csv(archive_read(vep_file, file = 1), col_types = cols(), sep="\t")
vep = fread(vep_file, data.table=F)
vep <- read.delim(file = untar(vep_file,compressed="gzip"),sep="\t")
vep$locus_id = paste(vep$chr, vep$start, vep$end, sep="_")
# have 2 columns called SNPID, rename one?
colnames(vep)[20] = "SNPID_rep"


# select only the case you want to explore
# Step 1: Filter rows with matching phenotype_id and locus_id
# Perform the join, keeping only SeqId and Ensembl_Gene_ID from mapping
# In cases of multiple matches, left_join will automatically duplicate the rows in vep for each match.
vep <- vep %>%
  inner_join(x, by = c("phenotype_id", "locus_id"))


# Count unique combinations
unique_combinations <- vep %>%
  distinct(phenotype_id, locus_id) %>%
  count()



# add the TARGET_GENE_NAME fir the seqid:
mapping_file = "/exchange/healthds/pQTL/pQTL_workplace/annotations/mapping_file/20241108_mapped_gene_file_GRCh37.txt"
mapping = fread(mapping_file, data.table=F)
mapping$SeqId = paste("seq.", mapping$SeqId, sep="")
mapping$SeqId = gsub("-", ".", mapping$SeqId)

# add the TARGET_GENE_NAME fir the seqid:
mapping_file = "/exchange/healthds/pQTL/pQTL_workplace/annotations/mapping_file/20241108_mapped_gene_file_GRCh37.txt"
mapping = fread(mapping_file, data.table=F)
mapping$SeqId = paste("seq.", mapping$SeqId, sep="")
mapping$SeqId = gsub("-", ".", mapping$SeqId)


# Perform the join, keeping only SeqId and Ensembl_Gene_ID from mapping
# In cases of multiple matches, left_join will automatically duplicate the rows in vep for each match.


vep_with_genes <- vep %>%
  left_join(mapping %>% select(SeqId, Ensembl_Gene_ID),
            by = c("phenotype_id" = "SeqId"))


# Rename the Ensembl_Gene_ID column to "target_gene"
vep_with_genes <- vep_with_genes %>%
  rename(target_gene = Ensembl_Gene_ID)

# Add the two new columns
vep_with_genes <- vep_with_genes %>%
  mutate(
    epitope_effect_target = Gene == target_gene & epitope_effect, # TRUE if match and epitope_effect is TRUE
    epitope_effect_all = ave(epitope_effect, locus_id, FUN = any) # TRUE if any epitope_effect is TRUE in the locus_id
  )


vep_with_genes <- vep_with_genes %>%
  mutate(
    rare_SNP_MINF = if_else(MINF < 0.01 | (1 - MINF) < 0.01, TRUE, FALSE),
    rare_SNP_MAXF = if_else(MAXF < 0.01 | (1 - MAXF) < 0.01, TRUE, FALSE)
  )


# Collapse by locus_id, summarizing columns as needed
collapsed_data <- vep_with_genes %>%
  group_by(phenotype_id, locus_id) %>%
  summarize(
    target_genes = paste(unique(target_gene), collapse = ", "), # Collapses target genes into a string
    epitope_effect_target = any(epitope_effect_target), # TRUE if any target row matches
    epitope_effect_all = any(epitope_effect),          # TRUE if any row has epitope_effect
    .groups = "drop"                                   # Removes grouping in the output
  ) %>% data.frame()


table(collapsed_data$epitope_effect_target)
table(collapsed_data$epitope_effect_all)