library(magrittr)
library(ggplot2)

if (exists("snakemake")){
  scores <- read.csv(snakemake@input$scores, row.names = 1)
  r2 <- read.csv(snakemake@input$r2)
  plt <- snakemake@output$plt
  ggs <- snakemake@output$ggs
} else {
  scores <- read.csv('results/mofacell/scores~Julio_OCEAN_Nereid~sn10x.csv', row.names = 1)
  r2 <- read.csv("results/mofacell/r2_by_sample~Julio_OCEAN_Nereid~sn10x.csv")
  plt <- 'plots/mofacell/scores_r2_boxplot~Julio_OCEAN_Nereid~sn10x.pdf'
  ggs <- 'results/mofacell/scores_r2_boxplot~Julio_OCEAN_Nereid~sn10x.rds'
}

scores <- scores %>% dplyr::rename(ID = Sample)
r2 <- r2 %>% dplyr::rename(ID = Group) %>%
  dplyr::filter(!Factor == 'All')
to_plot <- r2 %>%
  dplyr::left_join(scores, by = c('ID', 'Factor'))

plot_factor_r2 <- function(f){
  p <- to_plot %>%
    dplyr::filter(Factor == f) %>%
    ggplot(aes(x = Cohort, y = Score)) + 
    facet_wrap(vars(View), nrow = 3) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = R2), width = 0.3, height = 0, size = 1) +
    scale_color_viridis_c() +
    coord_flip() +
    theme_minimal() +
    labs(y = "Factor score", 
         x = "Disease group",
         subtitle = f)
  return(p)
}

my_plots <- lapply(unique(to_plot$Factor), plot_factor_r2)
pdf(plt, width = 12)
my_plots
dev.off()

saveRDS(my_plots, file = ggs)
