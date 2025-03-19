# Annotation of Loci
This is side repository to use annotation files from Michele F to assign a gene name to each lead variants or loci within the pQTL meta-analysis pipeline project.

## Gene Name Annotation Process

This process was designed to determine the gene names for **7,870 variants** mapped to different genetic regions associated with various protein sequence IDs.

### Step 1: Identifying Gene Names from Annotated Files
- A filename was generated for each variant using **seqid, locus, and SNPID** columns.
- These filenames were compared against **12,095 annotated files** provided by Michele F.
- The `grep_gene_lb()` algorithm was applied to extract the correct **gene names** where a match was found.
- Using this approach, **7,815 out of 7,870 variants** were successfully annotated.

### Step 2: Handling Missing Annotations
- For the **55 loci** where no annotation file was found:
  - Each missing **SNPID** was manually checked in **Locus Breaker**.
  - If an annotation file was found for the **same SNPID but associated with a different locus or seqid**, it was used to fill in missing gene names.
  - This step recovered **16 additional SNPIDs**, reducing the missing count to **39**.
- The remaining **39 SNPIDs** still lacked annotation and required additional data from Michele F.

### Step 3: Resolving Multiple Gene Names for a Single SNPID
- For SNPIDs associated with multiple gene names:
  - **GeneCards.org** was used to verify **gene aliases**.
  - When multiple names were listed, priority was given to the **older HGNC-approved gene name**.
- This ensured consistency in gene annotation across all variants.

This structured approach maximized gene annotation accuracy while addressing missing data through manual verification and external resources.
