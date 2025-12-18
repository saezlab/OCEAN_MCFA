library(magrittr)

print("INFO: Job running...")
if (exists('snakemake')){
  metadata <- read.csv(snakemake@input$metadata, row.names = 1)
  metadata_clean_csv <- snakemake@output$metadata_clean
} else {
  metadata <- readxl::read_xlsx('data/OCEAN_Metadata_v14_MiKTMC_Charlotte.xlsx')
  metadata_clean_csv <- 'data/metadata.csv'
}

print(head(metadata))
simplify_cohort <- function(cohort){
  ifelse(stringr::str_detect(cohort, "IgAN"), "IgAN",
         ifelse(cohort %in% c("AIN/T2D", "C1QN", "C3_Glomerulopathy",
                              "COVAN", "Immune_Cmpx_GN",
                              "MN_Vax", "MPGN"), "Other",
                ifelse(stringr::str_detect(cohort, "FSGS"), "FSGS", cohort)))
}

metadata_clean <- metadata %>%
  dplyr::mutate(across(where(is.character), ~dplyr::na_if(., "N/A"))) %>%
  dplyr::mutate(across(c(Age, NEPTUNE_eGFR_at_Bx, NEPTUNE_UPCR_at_Bx,
                         number_of_APOL1_risk_alleles, nGloms_NEPTUNE, SMHpct_NEPTUNE,
                         GMHpct_NEPTUNE, SEHpct_NEPTUNE, GEHpct_NEPTUNE, MESHYPERpct_NEPTUNE, SEGEPILESpct_NEPTUNE,
                         GLOEPILESpct_NEPTUNE, CELLCRSpct_NEPTUNE, FIBROCELLCRSpct_NEPTUNE,
                         InterstitialFibrosis_NEPTUNE, TubularAtrophy_NEPTUNE), ~as.numeric(.))) %>%
  dplyr::mutate(Cohort_Charlotte = simplify_cohort(Cohort))

write.csv(metadata_clean, metadata_clean_csv)
