


library(showtext)
font_add_google('Libre Franklin', 'franklin')
font_add_google('Domine', 'domine')
showtext_opts(dpi = 150)
showtext_auto()

lb_epitope_cojo <- data.table::fread(out_lb_epitop_cojo)

lb_epitope_cojo %>% 
  count(cis_or_trans, epitope_effect_cojo) %>% 
  #filter(cis_or_trans == "cis") %>%
  group_by(cis_or_trans) %>%
  summarise(n = sum(epitope_effect_cojo == "Yes", na.rm = T))

# summary table of three created columns
lb_epitope_cojo %>%
  pivot_longer(cols = starts_with("epitope_effect_cojo"),
               names_to = "source", values_to = "class") %>%
  count(cis_or_trans, class, source) %>%
  pivot_wider(names_from = source, values_from = n) %>%
  DT::datatable()

# trans loci with epitope reported in literature
lb_epitope_cojo %>%
  count(epitope_effect_cojo_high, uniprot_match) %>%
  DT::datatable()

# 
lb_epitope_cojo %>%
  filter(
    cis_or_trans == "trans",
    uniprot_match == "YES",
    epitope_effect_cojo_high
    ) %>%
  dplyr::select(
    phenotype_id, locus, SNPID, cis_or_trans, uniprot_match, unip_matching_study,
    epitope_effect_cojo_high, epitope_status, epitope_causing_cojo
    ) %>% View()
  #filter(str_detect(SNPID, "^19:")) %>% distinct(locus)
  #fwrite("trans_loci_with_epitope.tsv", sep = "\t", row.names = F)

#----------#
# contingency tables
lb_epitope_cojo %>%
  filter(cis_or_trans == "cis") %>%
  count(epitope_effect_cojo) %>%
  ggplot(aes(x = epitope_effect_cojo, y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n, vjust = -0.5), size = 6) +
  scale_y_continuous(breaks = seq(0, 1600, 100)) +
  labs(
    x = "\nEpitope Effect", 
    y = "#Loci", 
    title = "<span >*cis*</span> Loci with Epitope-Related COJO SNPs",
    subtitle = "For these loci, we found epitope effect for all COJO variants"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "franklin"),
    plot.title = ggtext::element_markdown(family = "domine", size = 17, face = 'bold', lineheight = 2),
    axis.ticks.length = unit(1.75, 'mm'),
    axis.ticks = element_line(size = .9),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    plot.margin = margin(t = 10, l = 10, b = 10, r = 5)
  )

ggsave(filename = "06-May-25_loci_with_epitope_cojo_variants.png", height = 5.5, width = 7.5)



#----------#
lb_epitope_cojo %>% #head(50) %>%
  ggplot(aes(x = reorder(cis_or_trans, prop_epitope_yes), y = prop_epitope_yes)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Locus", y = "Proportion of Epitope SNPs", title = "Epitope-Related SNPs by Locus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



lb_epitope_cojo %>% head(50) %>%
  ggplot(aes(x = chr, y = phenotype_id, fill = prop_epitope_yes)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "Proportion") +
  labs(title = "Epitope SNP Proportion by Study and Locus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



