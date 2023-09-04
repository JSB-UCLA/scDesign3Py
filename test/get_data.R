filename <- "Zhengmix4eq"

sce <- readRDS(paste0("./test_data/rds/", filename, ".rds"))

sceasy::convertFormat(Zhengmix4eq_sce, from = "sce", to = "anndata", main_layer = "counts", transfer_layers = c("logcounts","normcounts"), outFile = paste0("./test/test_data/", filename, ".h5ad"))
