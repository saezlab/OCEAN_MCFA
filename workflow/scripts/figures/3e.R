metadata <- read.csv("data/metadata.csv")
scores <- read.csv("results/mofacell/widescores~combined_sn10xLRCluster2.csv") %>%
  dplyr::mutate(Factor5=-Factor5)

df <- dplyr::left_join(scores, metadata, by = c("OCEAN_MS_Exp_ID")) %>%
  dplyr::mutate("Project (Batch)" = Project1,
                Age,
                "Cohort" = Cohort_Manuscript,
                "Self-ID Race" = Race,
                "eGFR" = NEPTUNE_eGFR_at_Bx,
                "Steroid at V3" = SteroidAtV3_NEPTUNE,
                "Interstitial Fibrosis" = InterstitialFibrosis_NEPTUNE,
                "Tubular Atrophy" = TubularAtrophy_NEPTUNE,
                "log(uPCR)" = log(NEPTUNE_UPCR_at_Bx))

pE1 <- ggplot(df, aes(x = Factor3, y = Factor5)) +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), alpha = 0.2) +
  geom_point(aes(color = NEPTUNE_eGFR_at_Bx)) +
  theme_classic() + labs(colour = "eGFR at biopsy") +
  scale_color_viridis_c(option = "magma", limits = c(0,190))

pE2 <- ggplot(df, aes(x = Factor3, y = Factor5)) +
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  geom_vline(aes(xintercept = 0), alpha = 0.2) +
  geom_point(aes(color = Age)) +
  theme_classic() + labs(colour = "Age") +
  scale_color_viridis_c(direction = -1) 


pE <- pE1 + pE2 & theme(legend.position = "top")

pE
