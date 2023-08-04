import itertools
from typing import Literal, Optional, Union

import anndata as ad
from rpy2.robjects import NULL, r
from rpy2.robjects.packages import importr

from ._utils._format import InputError, _anndata_to_sce, _other_to_list, _parapara, _strvec_none2ri, _typecheck


class scDesign3:
    def __init__(self) -> None:
        """

        Arguments:
        ----------
        sce: SingleCellExperiment R object changed from the given Anndata object

        assay_use: The name of the assay used for modeling

        all_covar: All covariates (explainary variables) used to model gene expression pattern

        bpparam: If @parallelization is bpmapply, this R object specifies the parameters in the parallelization function.

        construct_data_res: Result of calling @construct_data
        """
        # import r package
        self._rscdesign = importr("scDesign3")

    # get_bpparam
    @staticmethod
    def get_bpparam(**kwargs):
        """Get your parallelization parameters robject.

        Check R function `BiocParallel::MulticoreParam` for more informarion on how to set the parameters. Current available parameters including "workers","tasks","stop_on_error","progressbar","RNGseed","timeout","exportglobals","log","threshold","logdir","resultdir","jobname","force_GC","fallback","manager_hostname" and "manager_port".
        """
        para_list = [
            "workers",
            "tasks",
            "stop_on_error",
            "progressbar",
            "RNGseed",
            "timeout",
            "exportglobals",
            "log",
            "threshold",
            "logdir",
            "resultdir",
            "jobname",
            "force_GC",
            "fallback",
            "manager_hostname",
            "manager_port",
        ]
        para_dict = {key: value for key, value in kwargs.items() if key in para_list}
        if not para_dict:
            raise InputError("Please check R function BiocParallel::MulticoreParam and check parameter names.")
        else:
            return r("BiocParallel::MulticoreParam")(**para_dict)

    # construct data
    @_typecheck(
        anndata=ad.AnnData,
        corr_by=Union[str, list],
        assay_use=Optional[str],
        default_assay_name=Optional[str],
        celltype=Optional[str],
        pseudotime=Optional[Union[str, list]],
        spatial=Optional[list],
        other_covariates=Optional[Union[str, list]],
        ncell=Optional[int],
        parallelization=str,
        show=bool,
    )
    def construct_data(
        self,
        anndata: ad.AnnData,
        corr_by: Union[str, list[str]],
        assay_use: Optional[str] = None,
        default_assay_name: Optional[str] = None,
        celltype: Optional[str] = None,
        pseudotime: Optional[Union[str, list[str]]] = None,
        spatial: Optional[list[str]] = None,
        other_covariates: Optional[Union[str, list[str]]] = None,
        ncell: Optional[int] = None,
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "mcmapply",
        bpparam: Optional[r] = None,
        show: bool = False,
    ) -> None:
        """Construct the input data (covaraite matrix and expression matrix)

        This function constructs the input data for fit_marginal.

        This function takes a `anndata.AnnData` object as the input. Based on users' choice, it constructs the matrix of covaraites (explainary variables) and the expression matrix (e.g., count matrix for scRNA-seq).

        Parameters:
        ----------
        anndata: `anndata.AnnData`
            `anndata.AnnData` object to store the single cell experiment information.
        corr_by: `str` or `list[str]`
            Indicates the groups for correlation structure. If '1', all cells have one estimated corr. If 'ind', no corr (features are independent). If others, this variable decides the corr structures.
        assay_use: `str` (default: None)
            Indicates the assay you will use. If None, please specify a name for the assay stored in `anndata.AnnData.X` in @default_assay_name.
        default_assay_name: `str` (default: None)
            Specified only when @assay_use is None. Asign a name to your default single cell experiment.get_bpparamone)
            The name of cell type variable in the `anndata.AnnData.obs`.
        pseudotime: `str` or `list[str]` (default: None)
            The name of pseudotime and (if exist) multiple lineages in the `anndata.AnnData.obs`.
        spatial: `list[str]` (default: None)
            The names of spatial coordinates in the `anndata.AnnData.obs`.
        other_covariates: `str` or `list[str]` (default: None)
            The other covaraites you want to include in the data.
        ncell: `int` (default: None)
            The number of cell you want to simulate. Default is None, which means only the provided cells in the `anndata.AnnData` object will be used. If an arbitrary number is provided, the fucntion will use Vine Copula to simulate a new covaraite matrix.
        parallelization: `str` (default: 'mcmapply')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam.
        bpparam: `robject` (default: None)
            If @parallelization is 'bpmapply', first call method @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', remain default None.
        """

        # modify input
        corr_by, celltype, pseudotime, other_covariates = _other_to_list(
            corr_by, celltype, pseudotime, other_covariates
        )

        if ncell is None:
            ncell = anndata.n_obs

        self.all_covar = list(itertools.chain(*filter(None, [celltype, pseudotime, spatial, other_covariates])))
        if not self.all_covar:
            raise InputError("At least one covariate should be specified.")

        # Construct R SingleCellExperiment object
        self.sce, self.assay_use = _anndata_to_sce(
            data=anndata,
            assay_use=assay_use,
            default_assay_name=default_assay_name,
            covar=self.all_covar,
        )

        # change other parameters to R form
        celltype, pseudotime, spatial, other_covariates, corr_by = _strvec_none2ri(
            celltype, pseudotime, spatial, other_covariates, corr_by
        )

        # check parallelization parameter
        bpparam = _parapara(parallelization, bpparam)

        # Call scDesign3::construct_data
        self.construct_data_res = self._rscdesign.construct_data(
            sce=self.sce,
            assay_use=self.assay_use,
            celltype=celltype,
            pseudotime=pseudotime,
            spatial=spatial,
            other_covariates=other_covariates,
            ncell=ncell,
            corr_by=corr_by,
            parallelization=parallelization,
            BPPARAM=bpparam,
        )
