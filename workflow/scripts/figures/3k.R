library(magrittr)
library(ggplot2)
library(patchwork)

loadings <- read.csv("results/mofacell/loadings~combined_sn10xLRCluster2.csv",
                     row.names = 1)

loadings_f3 <- loadings %>% dplyr::filter(Factor == "Factor3")
to_plot_classic <- loadings_f3 %>%
  dplyr::filter(!stringr::str_detect(view, "&"))
to_plot_liana <- loadings_f3 %>%
  dplyr::filter(stringr::str_detect(view, "&"))


n = 10
loadings_to_plot <- to_plot_classic %>%
  dplyr::filter(view %in% c("DCT", "PT", "DTL_aPT")) %>%
  dplyr::group_by(view) %>%
  dplyr::slice_max(abs(Loading), n = 2) %>%
  dplyr::ungroup() %>%
  dplyr::select(variable) %>% unlist()

selected_variables <- c("HAVCR1", "PROM1", "VCAM1", "VIM", "TMSB4X", "TMSB10", "WFDC2", "UMOD", "TRPM6")

(
  p1b <- loadings_f3 %>%
    dplyr::filter(variable %in% selected_variables) %>%
    dplyr::filter(view %in% c("DCT", "PT", "DTL_aPT")) %>%
    dplyr::group_by(view) %>%
    dplyr::arrange(desc(Loading), .by_group = TRUE) %>%
    dplyr::mutate(variable = forcats::fct_inorder(variable)) %>%
    dplyr::mutate(variable = forcats::fct_rev(variable)) %>%
    ggplot(aes(y = variable, x = view, fill = Loading))
  + geom_tile()
  + scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    limits = c(-2, 2)
  )  + theme_classic()
  + labs(x = "View", y = element_blank(), subtitle = 'Transcriptome views')
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)), axis.text.y = element_text(size = rel(1)))
)

to_plot_liana <- loadings %>%
  dplyr::filter(Factor == "Factor3") %>%
  dplyr::filter(view %in% c("FIB&DCT", "FIB&DTL_aPT", "FIB&PT"))
n = 10

loadings_to_plot <- to_plot_liana %>%
  dplyr::slice_max(abs(Loading), n = n) %>%
  dplyr::select(variable) %>% unlist()

top_loadings <- to_plot_liana %>%
  dplyr::filter(variable %in% loadings_to_plot) %>%
  dplyr::group_by(view) %>%
  dplyr::arrange(desc(Loading), .by_group = TRUE) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable)) %>%
  dplyr::mutate(variable = forcats::fct_rev(variable))


(
  p2b <- ggplot(top_loadings, aes(x = view, y = variable, fill = Loading))
  + geom_tile()
  + scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    limits = c(-2, 2)
  )  + theme_classic()
  + labs(x = "View", y = element_blank(), subtitle = 'L-R views')
  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)), axis.text.y = element_text(size = rel(1)))
)


design <- "
111222"
p1b + p2b + plot_layout(design=design, guides = "collect") + plot_annotation(title = "Factor 3: Selected loadings")
