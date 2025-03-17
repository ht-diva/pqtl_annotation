
library(tidyverse)
library(data.table)

# Handling Missing Annotations
# - For the **55 loci** where no annotation file was found:
#   - Each missing **SNPID** was manually checked in **Locus Breaker**.
# - If an annotation file was found for the **same SNPID but associated with a different locus or seqid**, it was used to fill in missing gene names.
# - This step recovered **16 additional SNPIDs**, reducing the missing count to **39**.
# - The remaining **39 SNPIDs** still lacked annotation and required additional data from Michele F.

# Resolving Multiple Gene Names for a Single SNPID
# - For SNPIDs associated with multiple gene names:
#   - **GeneCards.org** was used to verify **gene aliases**.
# - When multiple names were listed, priority was given to the **older HGNC-approved gene name**.
# - This ensured consistency in gene annotation across all variants.


lb_filled <- lb_annot %>%
  # add annotated genes for the loci with available annotation file from Michele F
  left_join(
    lb_with_genes %>% dplyr::select(SNPID, phenotype_id, annot_gene_vep),
    join_by(SNPID, phenotype_id)
    ) %>%
  dplyr::mutate(
    annot_gene_vep = case_when(
      SNPID == "19:51628529:A:G" ~ "SIGLEC9", # SNP's annotation file is NOT present for a seqid but present for other seqids
      SNPID == "19:54759361:C:T" ~ "LILRA6",  # the same for below lead SNPs
      SNPID == "19:43793256:A:C" ~ "PSG9",
      SNPID == "19:51628529:A:G" ~ "SIGLEC9",
      SNPID == "19:54759361:C:T" ~ "LILRA6",
      SNPID == "16:1297810:A:G" ~ "TPSAB1",
      SNPID == "14:21385991:C:T" ~ "RNASE3",
      SNPID == "11:12073116:A:T" ~ "DKK3",
      SNPID == "9:136141870:C:T" ~ "SURF6",
      SNPID == "7:99971313:C:T" ~ "PILRB",
      SNPID == "5:96252432:A:G" ~ "ERAP2",
      SNPID == "3:49721532:A:G" ~ "APEH",
      SNPID == "3:186378731:A:C" ~ "HRG",
      SNPID == "1:154425456:A:G" ~ "IL6R",    # till here, SNPs with missing annotation file
      SNPID == "1:155201064:C:T" ~ "GBA1",    # GBA; GBA1       --> 1st of 32 SNPs with multiple gene names
      SNPID == "1:159182745:A:G" ~ "ACKR1",   # DARC; ACKR1
      SNPID == "1:213054712:A:G" ~ "SPATA45", # C1orf227; SPATA45
      SNPID == "2:277003:A:G" ~ "ALKAL2",     # FAM150B; ALKAL2
      SNPID == "2:160635077:A:G" ~ "MARCH7",  # MARCH7; MARCHF7
      SNPID == "2:207038421:G:T" ~ "GPR1",    # GPR1; CMKLR2
      SNPID == "2:233286654:C:T" ~ "ALPG",    # ALPPL2; ALPG
      SNPID == "2:241813496:C:T" ~ "MAB21L4", # C2orf54; MAB21L4
      SNPID == "5:102896752:A:T" ~ "MACIR",   # C5orf30; MACIR
      SNPID == "6:42940673:A:G" ~ "GNMT",     # GNMT; CNPY3-GNMT
      SNPID == "7:99332948:G:T" ~ "CYP3A7",   #CYP3A7; CYP3A7-CYP3A51P
      SNPID == "8:134198301:C:T" ~ "WISP1",   #WISP1; CCN4
      SNPID == "9:74265096:A:G" ~ "TMEM2",    #TMEM2; CEMIP2
      SNPID == "10:5922670:A:G" ~ "FBH1",     #FBXO18; FBH1
      SNPID == "10:115848372:A:C" ~ "CCDC186", #C10orf118; CCDC186
      SNPID == "11:63318067:C:T" ~ "PLAAT2",    #HRASLS2; PLAAT2
      SNPID == "11:119207446:C:T" ~ "MFRP",   #C1QTNF5; MFRP
      SNPID == "12:89921860:A:T" ~ "GALNT4",  #GALNT4; POC1B-GALNT4
      SNPID == "12:89927081:C:T" ~ "GALNT4",  #GALNT4; POC1B-GALNT4
      SNPID == "12:117176552:A:G" ~ "SPRING1", #C12orf49; SPRING1
      SNPID == "12:122356436:C:T" ~ "WDR66",  #WDR66; CFAP251
      SNPID == "12:122420459:G:T" ~ "WDR66",  #WDR66; CFAP251
      SNPID == "13:103701690:A:G" ~ "ERCC5",  #ERCC5; BIVM-ERCC5
      SNPID == "14:21155108:A:G" ~ "RNASE4",  #ANG; RNASE4
      SNPID == "14:21155108:A:G" ~ "RNASE4",  #ANG; RNASE4
      SNPID == "14:95027722:A:G" ~ "SERPINA4", #SERPINA5; SERPINA4
      SNPID == "17:17850095:A:G" ~ "DRC3",    #LRRC48; DRC3
      SNPID == "17:60759347:A:G" ~ "MARCH10", #MARCH10; MARCHF10
      SNPID == "17:73034193:C:G" ~ "ATP5H",   #ATP5H; ATP5PD
      SNPID == "19:19271498:C:G" ~ "BORCS8",  #MEF2BNB; BORCS8
      SNPID == "19:28063549:A:G" ~ "SLC25A1P5", #LINC00662; CTC-459F4.1; SLC25A1P5; ...
      TRUE ~ annot_gene_vep
    )
  ) %>%
  dplyr::select(- c(txtpath, locus)) # remove columns Solene suggested


# save filled results
data.table::fwrite(lb_filled, file = path_lb_gene, quote = F, row.names = F, sep = "\t")

