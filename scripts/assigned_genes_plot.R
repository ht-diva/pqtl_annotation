

# Here we find the number of annotated genes to each 
# lead variant at the identified locus in meta-analysis.
# we report number of genes:
#   1. in general for all 7,815 annotated loci
#   2. separated for cis/trans loci
#   3. separated for cis/trans loci split by chromosome

#------------#

out_plt_hist  <- "19-Mar-25_histogram_n_assigned_genes.png"
out_plt_panel <- "19-Mar-25_histogram_n_assigned_genes_cistrans_2panels.png"
out_plt_multi <- "19-Mar-25_histogram_n_assigned_genes_cistrans_multipanels.png"

#------------#

# function to find number of assigned genes to each lead SNP 
find_n_gene <- function(filepath, variant){
  data.table::fread(file = filepath) %>% 
    dplyr::filter(SNPID == variant) %>%
    dplyr::select(SYMBOL) %>% 
    n_distinct()
}


# LB with annotated file name n=7,815
lb_non_na <- lb_annot %>% dplyr::filter(!is.na(txtpath))


# iterate the function
n_genes <- map2_dbl(
  lb_non_na$txtpath, 
  lb_non_na$SNPID, 
  find_n_gene
  )


#----------------------------------------#
#-----      No. Genes Histogram     -----
#----------------------------------------#

# histogram of no. genes
lb_non_na %>%
  cbind(n_genes) %>% # join no. of genes to LB results
  ggplot() +
  geom_histogram(aes(n_genes),color = "steelblue2", fill = "steelblue3") +
  scale_x_continuous(breaks = seq(0, 250, 10))+
  scale_y_continuous(breaks = seq(0, 1200, 100))+
  annotate(geom="text", size = 5, x=200, y=1100, label="Median (min-max)") +
  annotate(geom="text", size = 5, x=204, y=1000, label="56 (4-233)") +
  labs(
    x = "\nNumber of genes to which a locus lead SNP is assigned",
    y = "Number of loci\n"
  )+
  theme_classic()

ggsave(out_plt_hist, plot=last_plot(), height=5.5, width=10, dpi=150)


#----------------------------------------#
#-----      No. Genes Cis/Trans     -----
#----------------------------------------#

lb_non_na %>%
  cbind(n_genes) %>%
  ggplot() +
  geom_histogram(aes(n_genes, fill = cis_or_trans), show.legend = F, color = "grey40") +
  scale_fill_manual(values = c("#405d27", "#f4a688"))+
  facet_grid(~cis_or_trans)+
  scale_x_continuous(breaks = seq(0, 250, 25))+
  scale_y_continuous(breaks = seq(0, 1200, 100))+
  labs(
    x = "\nNumber of genes to which a locus lead SNP is assigned",
    y = "Number of loci\n"
  )+
  theme_bw()

ggsave(out_plt_panel, plot=last_plot(), height=5.5, width=10, dpi=150)


#----------------------------------------#
#-----      No. Genes per CHRs      -----
#----------------------------------------#

lb_non_na %>%
  cbind(n_genes) %>%
  ggplot() +
  geom_density(aes(n_genes, fill = cis_or_trans), show.legend = F, color = NA, alpha = .6) +
  scale_fill_manual(values = c("#405d27", "#f4a688"))+
  facet_wrap(~chr, nrow = 4, scales = "free") +
  labs(
    x = "\nNumber of genes to which a locus lead SNP is assigned",
    y = "Number of loci\n"
  )+
  theme_bw()+
  theme(
    strip.background = element_blank()
  )

ggsave(out_plt_multi, plot=last_plot(), height=5.5, width=12, dpi=150)

