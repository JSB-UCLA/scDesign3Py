import itertools
from typing import Literal, Optional, Union

import anndata as ad
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.packages import importr

from ._utils._errors import InputError, SequentialError
from ._utils._format import _anndata2sce, _bpparamcheck, _other2list, _strvec_none2ri, _typecheck


class scDesign3:
    def __init__(self) -> None:
        """
        :sce:
            SingleCellExperiment R object changed from the given Anndata object

        :assay_use:
            The name of the assay used for modeling

        :all_covar:
            All covariates (explainary variables) used to model gene expression pattern

        :construct_data_res:
            Result of calling @construct_data

        :fit_marginal_res:
            Result of calling @fit_marginal

        :fit_copula_res:
            Result of calling @fit_copula

        :model_paras:
            Result of calling @extract_para

        :simu_res:
            Result of calling @simu_new
        """
        # import r package
        self._rscdesign = importr("scDesign3")

    # get_bpparam
    @staticmethod
    def get_bpparam(name=Literal["MulticoreParam", "SnowParam"], show=True, **kwargs):
        """Get your parallelization parameters robject.

        Check R function `BiocParallel::MulticoreParam` and `BiocParallel::SnowParam` for more informarion on how to set the parameters.

        Arguments:
        ----------
        name: `str`
            The type of your selected parallel parameter. If windows, choose SnowParam. If linux or mac, choose MulticoreParam or SnowParam.

        show: `bool` (default: True)
            Whether to print the constructed onject information on the screen.
        """

        if name == "MulticoreParam":
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
        elif name == "SnowParam":
            para_list = [
                "workers",
                "type",
                "tasks",
                "stop_on_error",
                "progressbar",
                "RNGseed",
                "timeout",
                "exportglobals",
                "exportvariables",
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
        else:
            raise InputError("Currently only support MulticoreParam for linux/mac and SnowParam for windows.")

        para_dict = {key: value for key, value in kwargs.items() if key in para_list}
        if not para_dict:
            raise InputError(
                "Please check R document BiocParallel for argument details. https://www.bioconductor.org/packages/devel/bioc/manuals/BiocParallel/man/BiocParallel.pdf"
            )
        elif name == "MulticoreParam":
            para_object = r("BiocParallel::MulticoreParam")(**para_dict)
        elif name == "SnowParam":
            para_object = r("BiocParallel::SnowParam")(**para_dict)

        if show:
            print(para_object)
        return para_object

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
        bpparam=Optional[ro.methods.RS4],
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
        bpparam: Optional[ro.methods.RS4] = None,
    ) -> ro.vectors.ListVector:
        """Construct the input data (covaraite matrix and expression matrix)

        This function constructs the input data for @fit_marginal.

        This function takes a `anndata.AnnData` object as the input. Based on users' choice, it constructs the matrix of covaraites (explainary variables) and the expression matrix (e.g., count matrix for scRNA-seq).

        Arguments:
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

        bpparam: `rpy2.robject.methods.RS4` (default: None)
            If @parallelization is 'bpmapply', first call method @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', remain default None.
        """

        # modify input
        corr_by, celltype, pseudotime, other_covariates = _other2list(
            corr_by, celltype, pseudotime, other_covariates
        )

        if ncell is None:
            ncell = anndata.n_obs

        self.all_covar = list(itertools.chain(*filter(None, [celltype, pseudotime, spatial, other_covariates])))
        if not self.all_covar:
            raise InputError("At least one covariate should be specified.")

        # Construct R SingleCellExperiment object
        self.sce, self.assay_use = _anndata2sce(
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
        bpparam = _bpparamcheck(parallelization, bpparam)

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

        return self.construct_data_res

    # fit_marginal
    @_typecheck(
        mu_formula=str,
        sigma_formula=str,
        family_use=Union[str, list],
        n_cores=int,
        usebam=bool,
        data=Optional[ro.vectors.ListVector],
        predictor=str,
        parallelization=str,
        bpparam=Optional[ro.methods.RS4],
        trace=bool,
    )
    def fit_marginal(
        self,
        mu_formula: str,
        sigma_formula: str,
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]],
        n_cores: int,
        usebam: bool,
        data: Optional[ro.vectors.ListVector] = None,
        predictor: str = "gene",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "mcmapply",
        bpparam: Optional[ro.methods.RS4] = None,
        trace: bool = False,
    ) -> ro.vectors.ListVector:
        """Fit the marginal models

        @fit_marginal fits the per-feature regression models.

        The function takes the result from @construct_data as the input, and fit the regression models for each feature based on users' specification.

        Arguments:
        ----------
        mu_formula: `str`
            A string of the mu parameter formula

        sigma_formula: `str`
            A string of the sigma parameter formula

        family_use: `str`
            A string or a vector of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution', 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution' and 'gaussian distribution' respectively.

        n_cores: `int`
            The number of cores to use.

        usebam: `bool`
            If True, call R function mgcv::bam for calculation acceleration.

        data: `rpy2.robject.vectors.ListVector` (default: None)
            The result of @construct_data. Default is None, using the class property @construct_data_res.

        predictor: `str` (default: 'gene')
            A string of the predictor for the gam/gamlss model. This is essentially just a name.

        parallelization: `str` (default: 'mcmapply')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam.

        bpparam: `rpy2.robject.methods.RS4` (default: None)
            If @parallelization is 'bpmapply', first call method @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', remain default None.

        trace: `bool` (default: False)
            If True, the warning/error log and runtime for gam/gamlss will be returned.
        """

        # use constructed data
        try:
            if data is None:
                data = self.construct_data_res
        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        # change other parameters to R form
        (family_use,) = _other2list(family_use)
        (family_use,) = _strvec_none2ri(family_use)

        # check parallelization parameter
        bpparam = _bpparamcheck(parallelization, bpparam)

        self.fit_marginal_res = self._rscdesign.fit_marginal(
            data=data,
            predictor=predictor,
            mu_formula=mu_formula,
            sigma_formula=sigma_formula,
            family_use=family_use,
            n_cores=n_cores,
            usebam=usebam,
            parallelization=parallelization,
            BPPARAM=bpparam,
            trace=trace,
        )

        return self.fit_marginal_res

    # fit copula
    @_typecheck()
    def fit_copula():
        pass

    @_typecheck()
    def extract_para():
        pass

    @_typecheck()
    def simu_new():
        pass
