# Annotation of Loci
This is side repository to assign a gene name to each locus lead variant and to assess presence of epitope effect within the pQTL meta-analysis using annotations from Michele F.

### Gene Name Annotation Process

This process was designed to determine the gene names for **7,870 variants** mapped to different genetic regions associated with various protein sequence IDs and was done in several steps:
- A filename was generated for each variant using **seqid, locus, and SNPID** columns.
- These filenames were compared against **12,095 annotated files** provided by Michele F.
- The `grep_gene_lb()` algorithm was applied to extract the correct **gene names** where a match was found.
- Using this approach, **7,815 out of 7,870 variants** were successfully annotated.



### Identification of Epitope Effect
Multiple steps were implemented to define epitpe-related information:
- Read the LB file which has path to each annotated TXT file.
- Loop through each row, read the corresponding annotated TXT file, and filter relevant rows.
- Check if the consequence column contains mid- or high-impact consequences.
- If there’s a match, store TRUE and the corresponding SNPID; otherwise, store FALSE and NA.

Finally, we reported below variables along with annotated LB results:
- **epitope_effect_cojo**: Yes/No/NA, Yes (only for cis) if **all COJO variants** (or proxies by r² > 0.8) have VEP-derived moderate or high impact consequence (as defined by Gene symbol) on the gene encoding the measured protein. No if this was not the case. NA for trans-pQTLs.
- **cis_epitope_effect**: Yes/No/NA, Yes (only for cis) if **any COJO variants** (or proxies by r² > 0.8) have VEP-derived moderate or high impact consequence (as defined by Gene symbol) on the gene encoding the measured protein. No if this was not the case. NA for trans-pQTLs.
- **epitope_causing_snp**: List of variants causing epitope effect for a locus.
- **affected_genes_from_vep**: List of all unique gene symbols for which an epitope effect was observed.
- **impact_variant_vep**: Yes/No, Yes if all COJO variants (or proxies by r² > 0.8) have VEP-derived moderate or high IMPACT consequence on any protein-coding gene. No otherwise.
