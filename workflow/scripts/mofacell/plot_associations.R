library(magrittr)
library(ggplot2)

if (exists("snakemake")){
  associations <- read.csv(snakemake@input$associations)
  plt1 <- snakemake@output$plt1
  plt1png <- snakemake@output$plt1png
  plt2 <- snakemake@output$plt2
  plt2png <- snakemake@output$plt2png
} else {
  associations <- read.csv('results/mofacell/associations~Julio_OCEAN_Nereid~LIANA_all~all.csv', row.names = 1)
  plt1 <- 'plots/mofacell/associations_clr~Julio_OCEAN_Nereid~LIANA_all~all.pdf'
  plt1png <- 'plots/mofacell/associations_clr~Julio_OCEAN_Nereid~LIANA_all~all.png'
  plt2 <- 'plots/mofacell/associations_cov~Julio_OCEAN_Nereid~LIANA_all~all.pdf'
  plt2png <- 'plots/mofacell/associations_cov~Julio_OCEAN_Nereid~LIANA_all~all.png'
}

(
  p1 <- associations %>%
    dplyr::filter(stringr::str_detect(variable, '_clr')) %>%
    dplyr::mutate(variable = stringr::str_remove(variable, '_clr'),
                  order = as.numeric(stringr::str_remove(factor, "Factor")),
                  star = ifelse(pval < 1e-5, "*", "")) %>%
    dplyr::arrange(order) %>%
    dplyr::mutate(factor = factor(factor, unique(factor))) %>%
    dplyr::mutate(variable = factor(variable, levels = rev(unique(variable)))) %>%
    ggplot(aes(y = variable, x = factor, fill = -log10(pval))) +
    scale_fill_viridis_c(option = "plasma")  + theme_minimal() +
    geom_tile() +
    geom_text(aes(label=star), color = "white", size = 2.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = 'Factor', y = 'Cell type',
         caption = "* indicates pval < 1e-5",
         subtitle = "ANOVA test for association between Factor score and CLR-transformed cell type proportions")
)

ggsave(plot = p1, filename = plt1, device = "pdf", width = 20, height = 12, units = "cm")
ggsave(plot = p1, filename = plt1png, device = "png", width = 20, height = 12, units = "cm")

exclude <- c("WhatIsComplimentaryEdgarID", "ID", "ProjectID", "Prep", "Method",
             "ContactCensusGeoID", "ContactCensusState", "ContactCensusCounty",
             "ContactCensusTract", "WGS_flag", "RNAseq_Glom", "RNAseq_TI",
             "snRNAseq_10x_flag", "snRNAseq_Parse_flag", "Somascan_plasma",
             "Somascan_urine", "DaysBxtoBL", "NEPTUNE_SAFID")

toplot2 <- associations %>%
  dplyr::filter(!stringr::str_detect(variable, '_clr')) %>%
  dplyr::filter(!variable %in% exclude) %>%
  dplyr::mutate(order = as.numeric(stringr::str_remove(factor, "Factor")),
                star = ifelse(pval < 0.001, "*", "")) %>%
  dplyr::arrange(order) %>%
  dplyr::mutate(factor = factor(factor, unique(factor))) %>%
  dplyr::mutate(variable = factor(variable, levels = rev(unique(variable))))

variables <- toplot2 %>%
  dplyr::filter(star == "*") %>%
  dplyr::select(variable) %>% unlist()

outliers <- c("Factor7", "Factor8",
              "Factor11", "Factor12",
              "Factor13", "Factor14",
              "Factor15", "Factor18",
              "Factor19", "Factor20")
(
  p2 <- toplot2 %>%
    dplyr::filter(variable %in% variables) %>%
    dplyr::filter(!factor %in% outliers) %>%
    ggplot(aes(y = variable, x = factor, fill = -log10(pval))) +
    scale_fill_viridis_c(option = "plasma")  + theme_minimal() +
    geom_tile() +
    geom_text(aes(label=star), color = "white", size = 2.5) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = rel(0.7))) +
    labs(x = 'Factor', y = 'Covariate',
         caption = "* indicates pval < 0.001",
         subtitle = "ANOVA test for association with metadata covariates")
)

ggsave(plot = p2, filename = plt2, device = "pdf", width = 20, height = 12, units = "cm")
ggsave(plot = p2, filename = plt2png, device = "png", width = 20, height = 12, units = "cm")

