filename <- "BATCH"

sce <- readRDS(paste0("./test/test_data/rds/", filename, ".rds"))

sceasy::convertFormat(sce, from = "sce", to = "anndata", main_layer = "counts", outFile = paste0("./test/test_data/", filename, ".h5ad"))
