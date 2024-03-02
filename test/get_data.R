filename <- "MOBSC"

# sce  <- readRDS(paste0("./test_data/rds/", filename, ".rds"))

# sce <- readRDS((url("https://figshare.com/ndownloader/files/40581983")))

sceasy::convertFormat(sce, from = "sce", to = "anndata", main_layer = "counts", outFile = paste0("./test_data/", filename, ".h5ad"))

