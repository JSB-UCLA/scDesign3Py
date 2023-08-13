library(scDesign3)

sce <- readRDS("./test/test_data/VISIUM.rds")

bpparam <- BiocParallel::SnowParam(stop.on.error = FALSE)

const_data <- construct_data(
    sce[1:10, 1:10],
    corr_by = "1",
    celltype = "cell_type",
    spatial = c("spatial1", "spatial2"),
    pseudotime = NULL,
    other_covariates = NULL,
)
const_data

marginal <- fit_marginal(
    data = const_data,
    mu_formula = "s(spatial1,spatial2,bs='gp',k=2)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE,
    # n_cores = 3,
    # parallelization = "bpmapply",
    # BPPARAM = bpparam,
)
marginal

copula <- fit_copula(
    sce = sce,
    assay_use = "counts",
    input_data = const_data$dat,
    family_use = "nb",
    n_cores = 1,
    marginal_list = marginal,
    # n_cores = 3,
    # parallelization = "bpmapply",
    # BPPARAM = bpparam,
)
copula
