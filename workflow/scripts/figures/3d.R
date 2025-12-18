metadata <- read.csv("data/metadata.csv")
scores <- read.csv("results/mofacell/widescores~combined_sn10xLRCluster2.csv") %>%
  dplyr::mutate(Factor5=-Factor5)
clr <- read.csv("results/preprocessing/clr.csv")

df <- dplyr::left_join(scores, metadata, by = c("OCEAN_MS_Exp_ID")) %>%
  dplyr::left_join(clr, by = c("OCEAN_MS_Exp_ID")) %>%
  dplyr::mutate("Project (Batch)" = Project1,
                Age,
                Pediatric,
                Sex,
                "Cohort" = Cohort_Manuscript,
                "Self-ID Race" = Race,
                "eGFR" = NEPTUNE_eGFR_at_Bx,
                "Steroid at V3" = SteroidAtV3_NEPTUNE,
                "Interstitial Fibrosis" = InterstitialFibrosis_NEPTUNE,
                "Tubular Atrophy" = TubularAtrophy_NEPTUNE,
                "log(uPCR)" = log(NEPTUNE_UPCR_at_Bx))

df %>%
  tidyr::pivot_longer(cols = starts_with("Factor"), names_to = "Factor", values_to = "Score") %>%
  dplyr::mutate(factor_num = as.integer(stringr::str_remove(Factor, "Factor"))) %>%
  dplyr::mutate(Factor = forcats::fct_reorder(Factor, factor_num)) %>%
  ggplot(aes(x = Score, y = Cohort_Manuscript)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.3) + facet_wrap(vars(Factor), scales = "free_x") + theme_classic()

# Assign Factors 8-20 as outlier-driven factors
factor_columns_to_plot <- c("Factor3", "Factor5")
additional_factors_tested <- c("Factor1", "Factor2", "Factor4", "Factor6", "Factor7")
metadata_columns_to_plot <- c("Project (Batch)", "Age", "Cohort", "Self-ID Race", "eGFR",
                      "Steroid at V3", "Interstitial Fibrosis",
                      "Tubular Atrophy", "log(uPCR)", "Sex", "Pediatric")
additional_tested_metadata_columns <- c("IgA_Present", "Allele_Markup_APOL1", "antiNephrin_Ever",
                                        "aNephrin_Bx", "Ethnicity", "N264K_APOL1", "number_of_APOL1_risk_alleles",
                                        "RAASBlockPreBx", "TNF_Grouping_JASN", "nGloms_NEPTUNE", "SMHpct_NEPTUNE",
                                        "GMHpct_NEPTUNE", "SEHpct_NEPTUNE", "GEHpct_NEPTUNE", "MESHYPERpct_NEPTUNE",
                                        "SEGEPILESpct_NEPTUNE", "GLOEPILESpct_NEPTUNE", "CELLCRSpct_NEPTUNE",
                                        "FIBROCELLCRSpct_NEPTUNE", "MEST.C_Score_M", "MEST.C_Score_E",
                                        "MEST.C_Score_S", "MEST.C_Score_T", "MEST.C_SC0re_C", "MEST.C_TotalScore")
metadata_columns <- c(metadata_columns_to_plot, additional_tested_metadata_columns)
factor_columns_to_test <- c(factor_columns_to_plot, additional_factors_tested)
# Separate continuous and categorical metadata
continuous_metadata <- metadata_columns[sapply(df[metadata_columns], is.numeric)]
categorical_metadata <- metadata_columns[sapply(df[metadata_columns], is.factor) | sapply(df[metadata_columns], is.character) | sapply(df[metadata_columns], is.logical)]
clr_metadata <- setdiff(colnames(clr), c("OCEAN_MS_Exp_ID"))

# Pearson correlation for continuous metadata
corr_results <- data.frame(Factor=character(), Metadata=character(), Correlation=numeric(), P_value=numeric())

for (factor in factor_columns_to_test) {
  for (meta in continuous_metadata) {
    if (sum(!is.na(df[[factor]])) > 1 & sum(!is.na(df[[meta]])) > 1) {
      test <- cor.test(df[[factor]], df[[meta]], use="complete.obs", method="pearson")
      corr_results <- rbind(corr_results, data.frame(Factor=factor, Metadata=meta, Correlation=test$estimate, P_value=test$p.value))
    }
  }
}

# ANOVA test for categorical metadata
anova_results <- data.frame(Factor=character(), Metadata=character(), F_statistic=numeric(), P_value=numeric())

for (factor in factor_columns_to_test) {
  for (meta in categorical_metadata) {
    print(meta)
    print(factor)
    if (length(unique(df[[meta]])) > 1) { # Ensure multiple groups exist
      model <- aov(df[[factor]] ~ as.factor(df[[meta]]), data=df)
      test <- summary(model)
      p_value <- test[[1]][["Pr(>F)"]][1]
      f_statistic <- test[[1]][["F value"]][1]
      dof <- test[[1]][["Df"]][1]
      anova_results <- rbind(anova_results,
                             data.frame(Factor=factor, Metadata=meta, F_statistic=f_statistic,
                                        P_value=p_value, DoF=dof))
    }
  }
}
corr_results$Test <- "Pearson Correlation"
anova_results$Test <- "ANOVA"
results <- dplyr::full_join(corr_results, anova_results, by = c("Metadata", "Factor", "P_value", "Test"))
# Apply FDR correction 
results$Adjusted_P_value <- p.adjust(results$P_value, method="fdr")
results$annotation <- ifelse(results$Adjusted_P_value < 0.01, "*", "")

corr_results <- results %>%
  dplyr::filter(Test == "Pearson Correlation", Metadata %in% metadata_columns_to_plot, Factor %in% factor_columns_to_plot)
# Plot Pearson Correlation Heatmap
p1 <- ggplot(corr_results, aes(y=Metadata, x=Factor, fill=Correlation)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="white", size = 5) +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-1, 1)) +
  labs(x="Metadata", y=element_blank(), fill="Pearson Correlation") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Plot anova Test Results Heatmap
anova_results <- results %>%
  dplyr::filter(Test == "ANOVA", Metadata %in% metadata_columns_to_plot, Factor %in% factor_columns_to_plot)
p2 <- ggplot(anova_results, aes(y=Metadata, x=Factor, fill=F_statistic)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="black", size = 5) +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  labs(x="Factor", y=element_blank(), fill="F statistic", caption="* Adjusted p-value < 0.01") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


design <- "
1
1
1
2
2
2
2
"

patchwork::wrap_plots(p1, p2, design=design) + patchwork::plot_layout(guides = "collect")

# Without legend
p1 <- ggplot(corr_results, aes(y=Metadata, x=Factor, fill=Correlation)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="white", size = 5) +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, limits=c(-1, 1)) +
  labs(x=element_blank(), y=element_blank(), fill="Pearson Correlation") +
  theme_minimal() + 
  theme(axis.text.y = element_text(size = 12), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = "none")

# Plot anova Test Results Heatmap
p2 <- ggplot(anova_results, aes(y=Metadata, x=Factor, fill=F_statistic)) +
  geom_tile() +
  geom_text(aes(label=annotation), color="black", size = 5) +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0) +
  labs(x=element_blank(), y=element_blank(), fill="F statistic") +
  theme_minimal() + theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12), legend.position = "none")


design <- "
1
1
1
2
2
2
2
"

patchwork::wrap_plots(p1, p2, design=design) + patchwork::plot_layout(guides = "collect")

write.csv(results, "results/downstream/factor_metadata_associations.csv")
