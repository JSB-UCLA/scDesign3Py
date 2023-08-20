library(scDesign3)

sce <- readRDS("F:/00工作/000清华/000上课/科研/暑研/JingYi Li/scDesign_Rpy/test/test_data/VISIUM.rds")

sce <- sce[1:10, 1:10]

bpparam <- BiocParallel::SnowParam(stop.on.error = FALSE)

const_data <- construct_data(
    sce,
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

paras <- extract_para(
    sce = sce,
    assay_use = "counts",
    marginal_list = marginal,
    new_covariate = const_data$newCovariate,
    data = const_data$dat,
    family_use = "nb",
    n_cores = 1
  )
paras

# a <- paras$mean_mat
# colnames(a) <- NULL
# row.names(a) <- NULL
# 
# b <- paras$sigma_mat
# colnames(b) <- NULL
# row.names(b) <- NULL
# 
# c <- paras$zero_mat
# colnames(c) <- NULL
# row.names(c) <- NULL


newcount <- simu_new(
  sce = sce,
  mean_mat = paras$mean_mat,
  sigma_mat = paras$sigma_mat,
  zero_mat = paras$zero_mat,
  # mean_mat = a,
  # sigma_mat = b,
  # zero_mat = c,
  quantile_mat = NULL,
  copula_list = copula$copula_list,
  n_cores = 1,
  family_use = "nb",
  input_data = const_data$dat,
  new_covariate = const_data$new_covariate,
  important_feature = copula$important_feature
)
newcount

whole <- scdesign3(sce = sce, corr_formula = "1",
                   celltype = "cell_type",
                   spatial = c("spatial1", "spatial2"),
                   pseudotime = NULL,
                   other_covariates = NULL,mu_formula = "s(spatial1,spatial2,bs='gp',k=2)",
                   sigma_formula = "1",
                   family_use = "nb",
                   n_cores = 1,
                   usebam = FALSE,
                   return_model = TRUE)
whole
