library(ggplot2)
library(magrittr)
library(dplyr)

if (exists("snakemake")){
  r2 <- read.csv(snakemake@input$r2)[,-1] # skip rownames (contains duplicates)
  data <- snakemake@output$data
  plt <- snakemake@output$plt1
  gg <- snakemake@output$gg1
} else {
  r2 <- read.csv('results/mofacell/r2_by_sample~Julio_OCEAN_Nereid~all.csv')
  data <- "results/mofacell/r2_by_sample_heatmap~Julio_OCEAN_Nereid~all.csv"
  plt <- "plots/mofacell/r2_by_sample_heatmap~Julio_OCEAN_Nereid~all.pdf"
  gg <- "results/mofacell/r2_by_sample_heatmap~Julio_OCEAN_Nereid~all.Rds"
}

to_plot <- r2 %>%
  dplyr::filter(Factor != "All") %>%
  dplyr::filter(Factor %in% c("Factor1", "Factor2", "Factor3"))

(
  p1 <- ggplot(to_plot, aes(x = View, y = Group, fill = R2))
  + geom_tile()
  + facet_wrap(.~Factor, ncol = 3)
  #+ geom_text(aes(label = sprintf("%.1f", R2)), color = "white", size = 3) 
  + scale_fill_viridis_c(option = "plasma", limits = c(-3, 100))  + theme_minimal()
  + labs(fill = expr(R^2), subtitle = expr(paste(R^2, " by sample")), x = element_blank(), y = element_blank())
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
)

write.csv(to_plot, data)
saveRDS(p1, file = gg)
ggsave(plot = p1, filename = plt, device = "pdf", width = 35, height = 22, units = "cm")

