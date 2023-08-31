filename <- ""

sce <- readRDS(paste0("./test_data/rds/", filename, ".rds"))

sceasy::convertFormat(sce, from = "sce", to = "anndata", outFile = paste0("./test_data/", filename, ".h5ad"))
