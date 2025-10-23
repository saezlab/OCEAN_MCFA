library(SeuratDisk)

if (exists("snakemake")){
  data <- snakemake@input$data
} else {
  data <- "data/Julio18.rds"
}
print(data)
stopifnot(file.exists(data))
seurat <- readRDS(data)
if (exists("snakemake")){
  SaveH5Seurat(seurat, snakemake@output$seurat, 
               overwrite = TRUE, verbose = TRUE)
}

