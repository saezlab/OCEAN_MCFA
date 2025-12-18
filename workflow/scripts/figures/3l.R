library(magrittr)
library(ggplot2)
library(patchwork)

factor_scores <- read.csv("results/mofacell/widescores~combined_sn10xLRCluster2.csv")
clr <- read.csv("results/preprocessing/clr.csv")
df <- dplyr::full_join(factor_scores, clr)

# R and padj from statistical testing
factor_columns_to_test <- c("Factor3")
clr_metadata <- setdiff(colnames(clr), c("OCEAN_MS_Exp_ID"))
# Pearson correlation for continuous metadata
corr_results <- data.frame(Factor=character(), Metadata=character(), Correlation=numeric(), P_value=numeric())
for (factor in factor_columns_to_test) {
  for (meta in clr_metadata) {
    if (sum(!is.na(df[[factor]])) > 1 & sum(!is.na(df[[meta]])) > 1) {
      test <- cor.test(df[[factor]], df[[meta]], use="complete.obs", method="pearson")
      corr_results <- rbind(corr_results, data.frame(Factor=factor, Metadata=meta, Correlation=test$estimate, P_value=test$p.value))
    }
  }
}
# Apply FDR correction for Pearson correlation p-values
corr_results$Adjusted_P_value <- p.adjust(corr_results$P_value, method="fdr")
#corr_results$annotation <- ifelse(corr_results$Adjusted_P_value < 0.05, "*", "")

r_immune <- corr_results %>%
  dplyr::filter(Metadata == "Immune") %>%
  dplyr::pull(Correlation)

padj_immune <- corr_results %>%
  dplyr::filter(Metadata == "Immune") %>%
  dplyr::pull(Adjusted_P_value)


p1 <- df %>%
  ggplot(aes(x = Factor3, y = Immune)) +
  ylim(c(-3,5)) +
  geom_point(color = "darkgray") + theme_classic() +
  scale_colour_brewer(palette = "Set2") + labs(subtitle = "Immune", x = "Factor 3", y = "CLR-transformed %") +
  geom_smooth(method='lm', se = FALSE, color = "black") +
  annotate("text", x = 0.2, y = 4.5,   label = paste0("italic(R) == ", signif(r_immune, 2),
                                                      " * \", padj = \" * ", signif(padj_immune, 2)),
           parse = TRUE
  )

r_fib <- corr_results %>%
  dplyr::filter(Metadata == "FIB") %>%
  dplyr::pull(Correlation)

padj_fib <- corr_results %>%
  dplyr::filter(Metadata == "FIB") %>%
  dplyr::pull(Adjusted_P_value)

p2 <- df %>%
  ggplot(aes(x = Factor3, y = FIB)) +
  ylim(c(-3,5)) +
  geom_point(color = "darkgray") + theme_classic() +
  scale_colour_brewer(palette = "Set2") + labs(subtitle = "FIB", x = "Factor 3", y = "CLR-transformed %") +
  geom_smooth(method='lm', se = FALSE, color = "black") +
  annotate("text", x = 0.2, y = 4.5,   label = paste0("italic(R) == ", signif(r_fib, 2),
                                                    " * \", padj = \" * ", signif(padj_fib, 2)),
           parse = TRUE
  )

p1/p2
