
# number of loci with epitope

library(tidyverse)

#----------#
# input
path_freez <- "/exchange/healthds/pQTL/results/META_CHRIS_INTERVAL/Locus_breaker_cojo_frozen_version_1812024/"
path_lb_epitop <- paste0(path_freez, "mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann_vep_epitope.tsv")

# output
out_plt_epitope     <- "23-Apr-25_epitope_loci_per_chromosome.png"
out_plt_epitope_all <- "23-Apr-25_epitope_loci_per_chromosome_ld_proxies.png"

#----------#
# number and percentage of cis loci with epitope effect
data.table::fread(path_lb_epitop) %>%
  dplyr::filter(cis_or_trans == "cis" & epitope_effect != "") %>%
  count(epitope_effect, name = "Count") %>%
  dplyr::mutate(Percent = Count / 1734)

#----------#
# read data and draw bar plot for cis loci
data.table::fread(path_lb_epitop) %>%
  dplyr::filter(cis_or_trans == "cis" & epitope_effect != "") %>%
  #Ensure chr is treated as an ordered factor
  dplyr::mutate(chr = factor(chr, levels = as.character(1:22))) %>%
  count(name = "Count", chr, epitope_effect) %>%
  ggplot(aes(x = chr, y = Count, fill = epitope_effect)) +
  #geom_col(aes(x = Chromosome, y = Count, fill = epitope_effect)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("No" = "#405d27", "Yes" = "#f4a688"))+
  scale_y_continuous(breaks = seq(0, 250, 10)) +
  labs(
    x = "\nChromosome",
    y = "Number of cis Loci with Epitope\n", 
    fill = "Epitope Effect"
    ) +
  theme_classic() + 
  theme(legend.position = c(0.92, 0.85))


# save the plot in "/home/dariush.ghasemi"
ggsave(out_plt_epitope, plot = last_plot(), height=5.5, width=10, dpi=150)

#----------#
# number of loci with epitope all
data.table::fread(path_lb_epitop) %>%
  dplyr::filter(epitope_effect_all != "") %>%
  count(epitope_effect_all, name = "Count") %>%
  dplyr::mutate(Percent = Count / sum(Count))


# bar plot for all loci
data.table::fread(path_lb_epitop) %>%
  dplyr::filter(epitope_effect_all != "") %>%
  dplyr::mutate(chr = factor(chr, levels = as.character(1:22))) %>%
  count(name = "Count", chr, epitope_effect_all) %>%
  ggplot(aes(x = chr, y = Count, fill = epitope_effect_all)) +
  geom_col() +
  #scale_fill_manual(values = c("No" = "#1f77b4", "Yes" = "#ff7f0e")) +
  scale_fill_manual(values = c("No" = "#00758F", "Yes" = "#f2ae72")) +
  scale_y_continuous(breaks = seq(0, 1200, 100)) +
  labs(
    x = "\nChromosome",
    y = "#Loci with Epitope (LD proxies)\n", 
    fill = "Epitope Effect"
  ) +
  theme_classic() + 
  theme(legend.position = c(0.92, 0.85))


# save the plot in "/home/dariush.ghasemi"
ggsave(out_plt_epitope_all, plot = last_plot(), height=5.5, width=10, dpi=150)

#----------#
