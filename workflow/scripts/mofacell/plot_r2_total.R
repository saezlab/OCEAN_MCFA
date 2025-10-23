library(ggplot2)
library(magrittr)
library(dplyr)

if (exists("snakemake")){
  r2 <- read.csv(snakemake@input$r2)
  data <- snakemake@output$data
  plt <- snakemake@output$plt
  gg <- snakemake@output$gg
} else {
  r2 <- read.csv('results/mofacell/r2_total~Julio_OCEAN_Nereid~sn10x.csv')
  data <- "results/mofacell/r2_total_heatmap~Julio_OCEAN_Nereid~sn10x.csv"
  plt <- "plots/mofacell/r2_total_heatmap~Julio_OCEAN_Nereid~sn10x.pdf"
  gg <- "results/mofacell/r2_total_heatmap~Julio_OCEAN_Nereid~sn10x.Rds"
}

to_plot <- r2 %>%
  dplyr::filter(Factor != "All") %>%
  dplyr::mutate(order = as.numeric(stringr::str_remove(Factor, "Factor"))) %>%
  dplyr::arrange(order) %>%
  dplyr::mutate(Factor = factor(Factor, unique(Factor)))

to_plot_classic <- to_plot %>%
  dplyr::filter(!stringr::str_detect(View, "&"))
to_plot_liana <- to_plot %>%
  dplyr::filter(stringr::str_detect(View, "&"))
(
  p1 <- to_plot_classic %>%
  dplyr::mutate(View = factor(View, levels = rev(unique(View)))) %>%
  ggplot(aes(x = View, y = Factor, fill = R2))
  + geom_tile()
  + geom_text(aes(label = sprintf("%.1f", R2)), color = "white", size = 3) 
  + scale_fill_viridis_c(option = "plasma", limits = c(0, max(to_plot$R2) + 10))  + theme_minimal()
  + labs(fill = expr(R^2), x = 'View')
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  + coord_flip()
)

if (nrow(to_plot_liana) > 1){
  (
  p2 <- to_plot_liana %>%
        dplyr::mutate(View = factor(View, levels = rev(unique(View)))) %>%
        ggplot(aes(x = View, y = Factor, fill = R2))
        + geom_tile()
        + geom_text(aes(label = sprintf("%.1f", R2)), color = "white", size = 3) 
        + scale_fill_viridis_c(option = "plasma", limits = c(0, max(to_plot$R2) + 10))  + theme_minimal()
        + labs(fill = expr(R^2), x = 'View')
        + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        + coord_flip()
  )
  my_plots <- list(p1, p2)
} else {
  my_plots <- list(p1)
}

write.csv(to_plot, data)
pdf(plt, width = 30, height = 20)
my_plots
dev.off()

saveRDS(my_plots, file = gg)

#saveRDS(p1, file = gg)
#ggsave(plot = p1, filename = plt, device = "pdf", width = 30, height = 20, units = "cm")

