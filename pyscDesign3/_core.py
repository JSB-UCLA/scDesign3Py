import itertools
import warnings
from typing import Literal, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.rlike.container import OrdDict
from rpy2.robjects import NULL
from rpy2.robjects.packages import importr

from ._utils._errors import InputError, SequentialError
from ._utils._format import _addname, _anndata2sce, _bpparamcheck, _other2list, _strvec_none2ri, _typecheck


class scDesign3:
    """Python interface for scDesign3

    All functions are arranged in scDesign3 class for easy reuse generated results.

    Attributes:
    -----------
    :sce:
        SingleCellExperiment R object changed from the given Anndata object

    :assay_use:
        The name of the assay used for modeling

    :all_covar:
        All covariates (explainary variables) used to model gene expression pattern

    :family_use:
        The set marginal distribution of each gene.

    :copula:
        The copula type you set.

    :n_cores:
        The number of cores used for model fitting.

    :parallelization:
        The parallelization method.

    :bpparam:
        If @parallelization is 'bpmapply', the corresponding R object to set the parallelization parameters.

    :return_py:
        Whether the functions will return a result easy for manipulation.

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

    :whole_pipeline_res:
        Result of calling @scdesign3
    """

    @_typecheck(
        n_cores=int,
        parallelization=str,
        bpparam=Optional[ro.methods.RS4],
        return_py=bool,
    )
    def __init__(
        self,
        n_cores: int = 1,
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "mcmapply",
        bpparam: Optional[ro.methods.RS4] = None,
        return_py: bool = True,
    ) -> None:
        """Decide basic settings when running the class methods.

        Details:
        ----------
        Decide basic settings including the parallel computation method, number cores to use and whether to return the more pythonic result.

        Arguments:
        ----------
        n_cores: `int` (default: 1)
            The number of cores to use.

        parallelization: `str` (default: 'mcmapply')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam.

        bpparam: `rpy2.robject.methods.RS4` (default: None)
            If @parallelization is 'bpmapply', first call method @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', remain default None.

        return_py: `bool` (default: True)
            If @return_py is True, a `rpy2.rlike.container.OrdDict` object is returned, which can be manipulated basically like the `dict` object and has `__getitem__()`, `get()`, `items()` method. What's more, values hosted in the object are transferred to `pandas`, `numpy` object. Else, a `rpy2.robjects.ListVector` object is returned, which also has the `items()` method, but if you want to get the value of a key, use `rx2()` method instead of `__getitem__()` and `get()` method. All values hosted basically remained as `rpy2.robjects.vectors` objects.

            It is recommended that if there's no need to manipulate the result, or if you are familiar with `rpy2` package, set the @return_py as False, else setting @return_py as True may help you better manipulate the results. Default is True.
        """
        # import r package
        self._rscdesign = importr("scDesign3")
        self.parallelization = parallelization
        self.n_cores = n_cores
        self.bpparam = _bpparamcheck(parallelization, bpparam)
        self.return_py = return_py
        self.__convert = ro.numpy2ri.converter + ro.default_converter + ro.pandas2ri.converter

    # construct data
    @_typecheck(
        anndata=ad.AnnData,
        corr_formula=Union[str, list],
        assay_use=Optional[str],
        default_assay_name=Optional[str],
        celltype=Optional[str],
        pseudotime=Optional[Union[str, list]],
        spatial=Optional[list],
        other_covariates=Optional[Union[str, list]],
        ncell=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        return_py=Union[bool, str],
    )
    def construct_data(
        self,
        anndata: ad.AnnData,
        corr_formula: Union[str, list[str]],
        assay_use: Optional[str] = None,
        default_assay_name: Optional[str] = None,
        celltype: Optional[str] = None,
        pseudotime: Optional[Union[str, list[str]]] = None,
        spatial: Optional[list[str]] = None,
        other_covariates: Optional[Union[str, list[str]]] = None,
        ncell: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        return_py: bool = "default",
    ) -> ro.vectors.ListVector:
        """Construct the input data

        This function constructs the input data  (covaraite matrix and expression matrix) for @fit_marginal.

        Details:
        ----------
        This function takes a `anndata.AnnData` object as the input. Based on users' choice, it constructs the matrix of covaraites (explainary variables) and the expression matrix (e.g., count matrix for scRNA-seq).

        Arguments:
        ----------
        anndata: `anndata.AnnData`
            `anndata.AnnData` object to store the single cell experiment information.

        corr_formula: `str` or `list[str]`
            Indicates the groups for correlation structure. If '1', all cells have one estimated corr. If 'ind', no corr (features are independent). If others, this variable decides the corr structures.

        assay_use: `str` (default: None)
            Indicates the assay you will use. If None, please specify a name for the assay stored in `anndata.AnnData.X` in @default_assay_name.

        default_assay_name: `str` (default: None)
            Specified only when @assay_use is None. Asign a name to your default single cell experiment.

        celltype: `str` (default: None)
            The name of cell type variable in the `anndata.AnnData.obs`.

        pseudotime: `str` or `list[str]` (default: None)
            The name of pseudotime and (if exist) multiple lineages in the `anndata.AnnData.obs`.

        spatial: `list[str]` (default: None)
            The names of spatial coordinates in the `anndata.AnnData.obs`.

        other_covariates: `str` or `list[str]` (default: None)
            The other covaraites you want to include in the data.

        ncell: `int` (default: 'default')
            The number of cell you want to simulate. Default is 'default', which means only the provided cells in the `anndata.AnnData` object will be used. If an arbitrary number is provided, the fucntion will use Vine Copula to simulate a new covaraite matrix.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        A `dict` like object.

        count_mat: `pandas.DataFrame`
            The expression matrix.

            The row corresponds to the observations and the column corresponds to the genes.

        dat: `pandas.DataFrame`
            The original covariate matrix.

        newCovariate: `pandas.DataFrame`
            The simulated new covariate matrix. If @ncell is default, the @newCovariate is basically the same as the @dat only without the `corr_group` column.

        filtered_gene: `None` or `list[str]`
            The genes that are excluded in the marginal and copula fitting steps because these genes only express in less than two cells.
        """

        if return_py == "default":
            return_py = self.return_py

        # modify input
        corr_formula, celltype, pseudotime, other_covariates = _other2list(
            corr_formula, celltype, pseudotime, other_covariates
        )

        if ncell == "default":
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
        celltype, pseudotime, spatial, other_covariates, corr_formula = _strvec_none2ri(
            celltype, pseudotime, spatial, other_covariates, corr_formula
        )

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if bpparam == "default":
            bpparam = self.bpparam
        else:
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
            corr_by=corr_formula,
            parallelization=parallelization,
            BPPARAM=bpparam,
        )

        if return_py:
            with self.__convert.context():
                res = ro.conversion.get_conversion().rpy2py(self.construct_data_res)
            res["count_mat"] = _addname(
                array=res["count_mat"],
                row_name=self.construct_data_res.rx2("count_mat").rownames,
                col_name=self.construct_data_res.rx2("count_mat").colnames,
            )
            if res["filtered_gene"] is NULL:
                res["filtered_gene"] = None
            else:
                res["filtered_gene"] = [i for i in res["filtered_gene"]]
            return res
        else:
            return self.construct_data_res

    # fit_marginal
    @_typecheck(
        mu_formula=str,
        sigma_formula=str,
        family_use=Union[str, list],
        usebam=bool,
        data=Union[ro.vectors.ListVector, OrdDict, str, dict],
        predictor=str,
        n_cores=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        trace=bool,
        return_py=Union[bool, str],
    )
    def fit_marginal(
        self,
        mu_formula: str,
        sigma_formula: str,
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]],
        usebam: bool,
        data: Union[ro.vectors.ListVector, OrdDict, dict] = "default",
        predictor: str = "gene",
        n_cores: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        trace: bool = False,
        return_py: bool = "default",
    ) -> ro.vectors.ListVector:
        """Fit the marginal models

        @fit_marginal fits the per-feature regression models.

        Details:
        ----------
        The function takes the result from @construct_data as the input, and fit the regression models for each feature based on users' specification.

        Arguments:
        ----------
        mu_formula: `str`
            A string of the mu parameter formula

        sigma_formula: `str`
            A string of the sigma parameter formula

        family_use: `str` or `list[str]`
            A string or a list of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution', 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution' and 'gaussian distribution' respectively.

        usebam: `bool`
            If True, call R function mgcv::bam for calculation acceleration.

        data: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict` (default: 'default')
            The result of @construct_data. Default is 'default', using the class property @construct_data_res.

        predictor: `str` (default: 'gene')
            A string of the predictor for the gam/gamlss model. This is essentially just a name.

        n_cores: `int` (default: 'default')
            The number of cores to use. Default is 'default', use the setting when initializing.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        trace: `bool` (default: False)
            If True, the warning/error log and runtime for gam/gamlss will be returned.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        A `dict` like object. Each key corresponds to one gene name and the length is equal to the total gene number.

        Every gene has the following keys.

        fit: `rpy2.rlike.container.OrdDict`
            The fitted regression models.

        removed_cell: `numpy.ndarray`
            The removed cell (observation) when fitting the marginal regression model.

        warning: `rpy2.rlike.container.OrdDict`
            If @trace = True, this key is returned. Hosting the warning/error log when fitting the marginal models.

        time: `numpy.ndarray`
            If @trace = True, this key is returned. The runtime for each marginal model. The first value corresponds to runtime for gam model and the second one corresponds to runtime for gamlss model.
        """

        if return_py == "default":
            return_py = self.return_py

        # use constructed data
        try:
            if data == "default":
                data = self.construct_data_res
            elif isinstance(data, OrdDict) or isinstance(data, dict):
                data = data.copy()
                try:
                    (tmp,) = _other2list(data["filtered_gene"])
                    (tmp,) = _strvec_none2ri(tmp)
                    data["filtered_gene"] = tmp
                except KeyError:
                    raise InputError("The result given should be the (modified) output of @construct_data.")

                with self.__convert.context():
                    data = ro.conversion.get_conversion().py2rpy(data)
        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        # change other parameters to R form
        self.family_use = family_use
        (family_use,) = _other2list(family_use)
        (family_use,) = _strvec_none2ri(family_use)

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if n_cores == "default":
            n_cores = self.n_cores
        if bpparam == "default":
            bpparam = self.bpparam
        else:
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

        if return_py:
            # with self.__convert.context():
            #     res = ro.conversion.get_conversion().rpy2py(self.fit_marginal_res)
            #     return res
            warnings.warn(
                "There's an unfixed problem in changing the @fit_marginal_res back to a more pythonic version. The return type is rpy2.robjects.vectors.ListVector, which is assigned method @items which is similar to dict. If you want to get the value using the key name, instead of using res[key_name], please use res.rx2(key_name)."
            )
            return self.fit_marginal_res
        else:
            return self.fit_marginal_res

    # fit copula
    @_typecheck(
        input_data=Union[str, pd.DataFrame],
        copula=str,
        empirical_quantile=bool,
        marginal_dict=Union[ro.vectors.ListVector, str, OrdDict, dict],
        family_use=Union[str, list],
        dt=bool,
        pseudo_obs=bool,
        epsilon=float,
        family_set=Union[str, list],
        important_feature=Union[str, list],
        n_cores=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        return_py=Union[bool, str],
    )
    def fit_copula(
        self,
        input_data: Union[Literal["count_mat", "dat", "newCovariate"], pd.DataFrame] = "dat",
        copula: Literal["gaussian", "vine"] = "gaussian",
        empirical_quantile: bool = False,
        marginal_dict: Union[ro.vectors.ListVector, OrdDict, dict] = "default",
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]] = "default",
        dt: bool = True,
        pseudo_obs: bool = False,
        epsilon: float = 1e-6,
        family_set: Union[str, list[str]] = ["gaussian", "indep"],
        important_feature: Union[Literal["all", "auto"], list[bool]] = "all",
        n_cores: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        return_py: bool = "default",
    ) -> ro.vectors.ListVector:
        """Fit the copula model

        @fit_copula fits the copula model.

        Details:
        ----------
        This function takes the result from @fit_marginal as the input and fit the copula model on the residuals.

        Arguments:
        ----------
        input_data: `str` or `pandas.DataFrame` (default: 'dat')
            One of the output of @construct_data, only need to specify the name if directly use the @construct_data_res. An alternative is to directly provided the corresponding DataFrame.

        copula: `str` (default: 'gaussian')
            A string of the copula choice. Must be one of 'gaussian' or 'vine'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.

        empirical_quantile: `bool` (default: False)
            Please only use it if you clearly know what will happen! If True, DO NOT fit the copula and use the EMPIRICAL CDF values of the original data; it will make the simulated data fixed (no randomness). Only works if ncell is the same as your original data.

        marginal_dict: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict` (default: 'default')
            The result of @fit_marginal. Default is 'default', using the class property @fit_marginal_res.

        family_use: `str` or `list[str]` (default: 'default')
            A string or a list of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'. Default is 'default', use the class property @family_use.

        dt: `bool` (default: True)
            If True, perform the distributional transformation to make the discrete data continuous. This is useful for discrete distributions (e.g., Poisson, NB). Note that for continuous data (e.g., Gaussian), DT does not make sense and should be set as False.

        pseudo_obs: `bool` (default: False)
            If True, use the empirical quantiles instead of theoretical quantiles for fitting copula.

        epsilon: `float` (default: 1e-6)
            A numeric variable for preventing the transformed quantiles to collapse to 0 or 1.

        family_set: `str` or `list[str]` (default: ['gaussian', 'indep'])
            A string or a string list of the bivariate copula families.

        important_feature: `str` or `list[bool]` (default: 'all')
            A string or list which indicates whether a gene will be used in correlation estimation or not. If this is a string, then this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a list, then this should be a logical vector with length equal to the number of genes in @sce. True in the logical vector means the corresponding gene will be included in gene-gene correlation estimation and False in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.

        n_cores: `int` (default: 'default')
            The number of cores to use. Default is 'default', use the setting when initializing.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        A `dict` like object.

        model_aic: `pandas.Series`
            An array with three values. In order, they are the marginal AIC, the copula AIC, the total AIC.

        model_bic: `pandas.Series`
            An array with three values. In order, they are the marginal BIC, the copula BIC, the total BIC.

        important_feature: `list[bool]`
            A vector showing the genes regarded as the inportant genes.

        copula_list: `rpy2.rlike.container.OrdDict`
            A dict of the fitted copula model. If using Gaussian copula, a dict of correlation matrices; if vine, a dict of vine objects.

            Caution that though it's name is `list`, it is actually a `dict` like object.
        """

        if return_py == "default":
            return_py = self.return_py

        if copula in ["gaussian", "vine"]:
            self.copula = copula
        else:
            raise InputError("The copula type could only be either 'gaussian' or 'vine'.")

        # check sce and assay use
        try:
            sce = self.sce
            assay_use = self.assay_use
        except AttributeError:
            raise SequentialError("Please run @construct_data first to get the R sce object.")

        # change other parameters to R form
        if family_use == "default":
            family_use = self.family_use
        family_use, family_set = _other2list(family_use, family_set)
        family_use, family_set = _strvec_none2ri(family_use, family_set)

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if n_cores == "default":
            n_cores = self.n_cores
        if bpparam == "default":
            bpparam = self.bpparam
        else:
            bpparam = _bpparamcheck(parallelization, bpparam)

        # check important feature
        if isinstance(important_feature, list):
            important_feature = ro.vectors.BoolVector(important_feature)
        elif not (important_feature == "all" or important_feature == "auto"):
            raise InputError("important_feature should be 'all', 'auto', or a list.")

        # use construct data res
        try:
            if isinstance(input_data, pd.DataFrame):
                with self.__convert.context():
                    input_data = ro.conversion.get_conversion().py2rpy(input_data)
            elif input_data in ["count_mat", "dat", "newCovariate"]:
                input_data = self.construct_data_res.rx2(input_data)

        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        # use fit marginal res
        try:
            if marginal_dict == "default":
                marginal_dict = self.fit_marginal_res
            elif isinstance(marginal_dict, OrdDict) or isinstance(marginal_dict, dict):
                with self.__convert.context():
                    marginal_dict = ro.conversion.get_conversion().py2rpy(marginal_dict)
        except AttributeError:
            raise SequentialError("Please first run @fit_marginal.")

        self.fit_copula_res = self._rscdesign.fit_copula(
            sce=sce,
            assay_use=assay_use,
            input_data=input_data,
            empirical_quantile=empirical_quantile,
            marginal_list=marginal_dict,
            family_use=family_use,
            copula=copula,
            DT=dt,
            pseudo_obs=pseudo_obs,
            epsilon=epsilon,
            family_set=family_set,
            important_feature=important_feature,
            n_cores=n_cores,
            parallelization=parallelization,
            BPPARAM=bpparam,
        )

        if return_py:
            with self.__convert.context():
                res = ro.conversion.get_conversion().rpy2py(self.fit_copula_res)
            res["model_aic"] = _addname(
                array=res["model_aic"],
                row_name=self.fit_copula_res.rx2("model_aic").names,
            )
            res["model_bic"] = _addname(
                array=res["model_bic"],
                row_name=self.fit_copula_res.rx2("model_bic").names,
            )
            res["important_feature"] = [i for i in res["important_feature"]]

            if copula == "gaussian":
                res["copula_list"] = {
                    key: _addname(
                        array=value,
                        row_name=self.fit_copula_res.rx2("copula_list").rx2(key).rownames,
                        col_name=self.fit_copula_res.rx2("copula_list").rx2(key).colnames,
                    )
                    for key, value in res["copula_list"].items()
                }

            return res
        else:
            return self.fit_copula_res

    @_typecheck(
        data=Union[str, pd.DataFrame],
        new_covariate=Union[str, pd.DataFrame, None],
        marginal_dict=Union[ro.vectors.ListVector, str, OrdDict, dict],
        family_use=Union[str, list],
        n_cores=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        return_py=Union[bool, str],
    )
    def extract_para(
        self,
        data: pd.DataFrame = "dat",
        new_covariate: Optional[pd.DataFrame] = "newCovariate",
        marginal_dict: Union[ro.vectors.ListVector, OrdDict, dict] = "default",
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]] = "default",
        n_cores: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        return_py: bool = "default",
    ) -> ro.vectors.ListVector:
        """Extract the parameters of each cell's distribution

        This function generates parameter matricies which determine each cell's distribution

        Details:
        ----------
        The function takes the new covariate (if use) from @construct_data and marginal models from @fit_marginal.

        Arguments:
        ----------
        data: `pandas.DataFrame` (default: 'dat')
            The data used for fitting the gene marginal models. Default is 'dat', use the @construct_data_res, 'dat' output.

        new_covariate: `pandas.DataFrame` or None (default: 'newCovariate')
            The new covariates to simulate new gene expression data using the gene marginal models. Default is 'newCovariate', use the @construct_data_res, 'newCovariate' output.

        marginal_dict: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict` (default: 'default')
            The result of @fit_marginal. Default is 'default', using the class property @fit_marginal_res.

        family_use: `str` or `list[str]` (default: 'default')
            A string or a list of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'. Default is 'default', use the class property @family_use.

        n_cores: `int` (default: 'default')
            The number of cores to use. Default is 'default', use the setting when initializing.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        A `dict` like object.

        mean_mat: `pandas.DataFrame`
            A matrix of the mean parameter.

            The row corresponds to the observations and the column corresponds to the genes.

        sigma_mat: `pandas.DataFrame`
            A matrix of the sigma parameter (for Gaussian, the variance; for NB, the dispersion.).

            The row corresponds to the observations and the column corresponds to the genes.

        zero_mat: `pandas.DataFrame`
            A matrix of the zero-inflation parameter (only non-zero for ZIP and ZINB).

            The row corresponds to the observations and the column corresponds to the genes.
        """

        if return_py == "default":
            return_py = self.return_py

        # check sce and assay use
        try:
            sce = self.sce
            assay_use = self.assay_use
        except AttributeError:
            raise SequentialError("Please run @construct_data first to get the R sce object.")

        # check what is data and new_covariate
        try:
            if isinstance(data, pd.DataFrame):
                with self.__convert.context():
                    data = ro.conversion.get_conversion().py2rpy(data)
            elif data == "dat":
                data = self.construct_data_res.rx2(data)

        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        try:
            if isinstance(new_covariate, pd.DataFrame):
                with self.__convert.context():
                    new_covariate = ro.conversion.get_conversion().py2rpy(new_covariate)
            elif new_covariate == "newCovariate":
                new_covariate = self.construct_data_res.rx2(new_covariate)
            elif new_covariate is None:
                new_covariate = NULL
        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        # use fit marginal res
        try:
            if marginal_dict == "default":
                marginal_dict = self.fit_marginal_res
            elif isinstance(marginal_dict, OrdDict) or isinstance(marginal_dict, dict):
                with self.__convert.context():
                    marginal_dict = ro.conversion.get_conversion().py2rpy(marginal_dict)
        except AttributeError:
            raise SequentialError("Please first run @fit_marginal.")

        # change other parameters to R form
        if family_use == "default":
            family_use = self.family_use
        (family_use,) = _other2list(family_use)
        (family_use,) = _strvec_none2ri(family_use)

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if n_cores == "default":
            n_cores = self.n_cores
        if bpparam == "default":
            bpparam = self.bpparam
        else:
            bpparam = _bpparamcheck(parallelization, bpparam)

        self.model_paras = self._rscdesign.extract_para(
            sce=sce,
            assay_use=assay_use,
            marginal_list=marginal_dict,
            family_use=family_use,
            new_covariate=new_covariate,
            data=data,
            n_cores=n_cores,
            parallelization=parallelization,
            BPPARAM=bpparam,
        )

        if return_py:
            with self.__convert.context():
                res = ro.conversion.get_conversion().rpy2py(self.model_paras)
            res["mean_mat"] = _addname(
                array=res["mean_mat"],
                row_name=self.model_paras.rx2("mean_mat").rownames,
                col_name=self.model_paras.rx2("mean_mat").colnames,
            )
            res["sigma_mat"] = _addname(
                array=res["sigma_mat"],
                row_name=self.model_paras.rx2("sigma_mat").rownames,
                col_name=self.model_paras.rx2("sigma_mat").colnames,
            )

            res["zero_mat"] = ro.r["as.matrix"](self.model_paras.rx2("zero_mat"))
            with self.__convert.context():
                tmp = ro.conversion.get_conversion().rpy2py(res["zero_mat"])
            res["zero_mat"] = _addname(
                array=tmp,
                row_name=res["zero_mat"].rownames,
                col_name=res["zero_mat"].colnames,
            )

            return res
        else:
            return self.model_paras

    @_typecheck(
        mean_mat=Union[str, np.ndarray, pd.DataFrame],
        sigma_mat=Union[str, np.ndarray, pd.DataFrame],
        zero_mat=Union[str, np.ndarray, pd.DataFrame],
        quantile_mat=Optional[Union[np.ndarray, pd.DataFrame]],
        copula_dict=Union[ro.vectors.ListVector, str, None, OrdDict, dict],
        input_data=Union[str, pd.DataFrame],
        new_covariate=Union[str, pd.DataFrame],
        family_use=Union[str, list],
        important_feature=Union[str, list],
        fastmvn=bool,
        nonnegative=bool,
        nonzerovar=bool,
        filtered_gene=Union[None, list, str],
        n_cores=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        return_py=Union[bool, str],
    )
    def simu_new(
        self,
        mean_mat: Union[np.ndarray, pd.DataFrame] = "mean_mat",
        sigma_mat: Union[np.ndarray, pd.DataFrame] = "sigma_mat",
        zero_mat: Union[np.ndarray, pd.DataFrame] = "zero_mat",
        quantile_mat: Optional[Union[np.ndarray, pd.DataFrame]] = None,
        copula_dict: Optional[Union[ro.vectors.ListVector, OrdDict, dict]] = "default",
        input_data: pd.DataFrame = "dat",
        new_covariate: pd.DataFrame = "newCovariate",
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]] = "default",
        important_feature: Union[Literal["all", "auto"], list[bool]] = "all",
        fastmvn: bool = False,
        nonnegative: bool = True,
        nonzerovar: bool = False,
        filtered_gene: Optional[list[str]] = "default",
        n_cores: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        return_py: bool = "default",
    ) -> ro.vectors.FloatMatrix:
        """Simulate new data

        Generate new simulated data based on fitted marginal and copula models.

        Details:
        ----------
        The function takes the new covariate (if use) from @construct_data, parameter matricies from @extract_para and multivariate Unifs from @fit_copula.

        Arguments:
        ----------
        mean_mat: `numpy.ndarray` or `pandas.DataFrame` (default: 'mean_mat')
            A matrix of the mean parameter. Default is 'mean_mat', use the @model_paras, mean_mat output.

        sigma_mat: `numpy.ndarray` or `pandas.DataFrame` (default: 'sigma_mat')
            A matrix of the sigma parameter. Default is 'sigma_mat', use the @model_paras, sigma_mat output.

        zero_mat: `numpy.ndarray` or `pandas.DataFrame` (default: 'zero_mat')
            A matrix of the zero-inflation parameter. Default is 'zero_mat', use the @model_paras, zero_mat output.

        quantile_mat: `numpy.ndarray` or `pandas.DataFrame` (default: None)
            A matrix of the multivariate quantile. Default is None, if parameter @copula_dict is provided.

        copula_dict: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict` (default: 'default')
            Copulas for generating the multivariate quantile matrix. Default is 'default', use the @fit_copula_res, copula_list output.

        data: `str` or `pandas.DataFrame` (default: 'dat')
            An input count matrix. Default is 'dat', use the @construct_data_res, 'dat' output.

        new_covariate: `str` or `pandas.DataFrame` (default: 'newCovariate')
            A dataframe which contains covariates of targeted simulated data from @construct_data. Default is 'newCovariate', use the @construct_data_res, 'newCovariate' output.

        family_use: `str` or `list[str]` (default: 'default')
            A string or a list of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'. Default is 'default', use the class property @family_use.

        important_feature: `str` or `list[bool]` (default: 'all')
            A string or list which indicates whether a gene will be used in correlation estimation or not. If this is a string, then this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a list, then this should be a logical vector with length equal to the number of genes in @sce. True in the logical vector means the corresponding gene will be included in gene-gene correlation estimation and False in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.

        fastmvn: `bool` (default: False)
            If True, the sampling of multivariate Gaussian is done by R function mvnfast, otherwise by R function mvtnorm.

        nonnegative: `bool` (default: True)
            If True, values < 0 in the synthetic data will be converted to 0. Default is True, since the expression matrix is nonnegative.

        nonzerovar: `bool` (default: False)
            If True, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.

        filtered_gene: `None` or `list[str]` (default: 'default')
            None or a list which contains genes that are excluded in the marginal and copula fitting steps because these genes only express in less than two cells. Default is 'default', use the @construct_data_res, 'filtered_gene' output.

        n_cores: `int` (default: 'default')
            The number of cores to use. Default is 'default', use the setting when initializing.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        `pandas.DataFrame`
            The new simulated count (expression) matrix.

            The row corresponds to the observations and the column corresponds to the genes.
        """

        if return_py == "default":
            return_py = self.return_py

        # check what is data and new_covariate
        try:
            if isinstance(input_data, pd.DataFrame):
                with self.__convert.context():
                    input_data = ro.conversion.get_conversion().py2rpy(input_data)
            elif input_data == "dat":
                input_data = self.construct_data_res.rx2(input_data)

        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        try:
            if isinstance(new_covariate, pd.DataFrame):
                with self.__convert.context():
                    new_covariate = ro.conversion.get_conversion().py2rpy(new_covariate)
            elif new_covariate == "newCovariate":
                new_covariate = self.construct_data_res.rx2(new_covariate)

        except AttributeError:
            raise SequentialError("Please first run @construct_data.")

        # check the matrix
        try:
            if isinstance(mean_mat, np.ndarray) or isinstance(mean_mat, pd.DataFrame):
                with self.__convert.context():
                    mean_mat = ro.conversion.get_conversion().py2rpy(mean_mat)
            elif mean_mat == "mean_mat":
                mean_mat = self.model_paras.rx2(mean_mat)
        except AttributeError:
            raise SequentialError("mean_mat is not provided. Please first run @extract_para.")

        try:
            if isinstance(sigma_mat, np.ndarray) or isinstance(sigma_mat, pd.DataFrame):
                with self.__convert.context():
                    sigma_mat = ro.conversion.get_conversion().py2rpy(sigma_mat)
            elif sigma_mat == "sigma_mat":
                sigma_mat = self.model_paras.rx2(sigma_mat)
        except AttributeError:
            raise SequentialError("sigma_mat is not provided. Please first run @extract_para.")

        try:
            if isinstance(zero_mat, np.ndarray) or isinstance(zero_mat, pd.DataFrame):
                with self.__convert.context():
                    zero_mat = ro.conversion.get_conversion().py2rpy(zero_mat)
            elif zero_mat == "zero_mat":
                zero_mat = self.model_paras.rx2(zero_mat)
        except AttributeError:
            raise SequentialError("zero_mat is not provided. Please first run @extract_para.")

        # check the copula lists
        if quantile_mat is None:
            quantile_mat = NULL
        else:
            with self.__convert.context():
                quantile_mat = ro.conversion.get_conversion().py2rpy(quantile_mat)

        try:
            if copula_dict == "default":
                copula_dict = self.fit_copula_res.rx2("copula_list")

            elif isinstance(copula_dict, OrdDict) or isinstance(copula_dict, dict):
                if self.copula == "gaussian":
                    with self.__convert.context():
                        # first change each value to R data.frame
                        tmp = {
                            key: ro.conversion.get_conversion().py2rpy(value) for key, value in copula_dict.items()
                        }
                    with ro.default_converter.context():
                        # then change R data.frame to R named matrix
                        copula_dict = ro.conversion.get_conversion().py2rpy(
                            {key: ro.r["as.matrix"](value) for key, value in tmp.items()}
                        )

            elif copula_dict is None:
                copula_dict = NULL
        except AttributeError:
            raise SequentialError("Please first run @fit_copula.")

        # check sce and assay use
        try:
            sce = self.sce
            assay_use = self.assay_use
        except AttributeError:
            raise SequentialError("Please run @construct_data first to get the R sce object.")

        # change other parameters to R form
        if family_use == "default":
            family_use = self.family_use
        (family_use,) = _other2list(family_use)
        (family_use,) = _strvec_none2ri(family_use)

        if filtered_gene == "default":
            filtered_gene = self.construct_data_res.rx2("filtered_gene")
        else:
            (filtered_gene,) = _other2list(filtered_gene)
            (filtered_gene,) = _strvec_none2ri(filtered_gene)

        # check important feature
        if isinstance(important_feature, list):
            important_feature = ro.vectors.BoolVector(important_feature)
        elif not (important_feature == "all" or important_feature == "auto"):
            raise InputError("important_feature should be 'all', 'auto', or a list.")

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if n_cores == "default":
            n_cores = self.n_cores
        if bpparam == "default":
            bpparam = self.bpparam
        else:
            bpparam = _bpparamcheck(parallelization, bpparam)

        self.simu_res = self._rscdesign.simu_new(
            sce=sce,
            assay_use=assay_use,
            mean_mat=mean_mat,
            sigma_mat=sigma_mat,
            zero_mat=zero_mat,
            quantile_mat=quantile_mat,
            copula_list=copula_dict,
            family_use=family_use,
            fastmvn=fastmvn,
            nonnegative=nonnegative,
            nonzerovar=nonzerovar,
            filtered_gene=filtered_gene,
            input_data=input_data,
            new_covariate=new_covariate,
            important_feature=important_feature,
            n_cores=n_cores,
            parallelization=parallelization,
            BPPARAM=bpparam,
        )

        if return_py:
            with self.__convert.context():
                res = ro.conversion.get_conversion().rpy2py(self.simu_res)
            res = _addname(
                array=res.T,
                row_name=self.simu_res.colnames,
                col_name=self.simu_res.rownames,
            )
            return res
        else:
            return self.simu_res

    @_typecheck(
        anndata=ad.AnnData,
        corr_formula=Union[str, list],
        mu_formula=str,
        sigma_formula=str,
        family_use=Union[str, list],
        usebam=bool,
        assay_use=Optional[str],
        default_assay_name=Optional[str],
        celltype=Optional[str],
        pseudotime=Optional[Union[str, list]],
        spatial=Optional[list],
        other_covariates=Optional[Union[str, list]],
        ncell=Union[int, str],
        copula=str,
        empirical_quantile=bool,
        family_set=Union[str, list],
        important_feature=Union[str, list],
        dt=bool,
        pseudo_obs=bool,
        fastmvn=bool,
        nonnegative=bool,
        nonzerovar=bool,
        return_model=bool,
        n_cores=Union[int, str],
        parallelization=str,
        bpparam=Optional[Union[ro.methods.RS4, str]],
        trace=bool,
        return_py=Union[bool, str],
    )
    def scdesign3(
        self,
        anndata: ad.AnnData,
        corr_formula: Union[str, list[str]],
        mu_formula: str,
        sigma_formula: str,
        family_use: Union[Literal["binomial", "poisson", "nb", "zip", "zinb", "gaussian"], list[str]],
        usebam: bool,
        assay_use: Optional[str] = None,
        default_assay_name: Optional[str] = None,
        celltype: Optional[str] = None,
        pseudotime: Optional[Union[str, list[str]]] = None,
        spatial: Optional[list[str]] = None,
        other_covariates: Optional[Union[str, list[str]]] = None,
        ncell: int = "default",
        copula: Literal["gaussian", "vine"] = "gaussian",
        empirical_quantile: bool = False,
        family_set: Union[str, list[str]] = ["gaussian", "indep"],
        important_feature: Union[Literal["all", "auto"], list[bool]] = "all",
        dt: bool = True,
        pseudo_obs: bool = False,
        fastmvn: bool = False,
        nonnegative: bool = True,
        nonzerovar: bool = False,
        return_model: bool = False,
        n_cores: int = "default",
        parallelization: Literal["mcmapply", "bpmapply", "pbmcmapply"] = "default",
        bpparam: Optional[ro.methods.RS4] = "default",
        trace: bool = False,
        return_py: bool = "default",
    ) -> ro.vectors.ListVector:
        """The wrapper for the whole scDesign3 pipeline

        Arguments:
        ----------
        anndata: `anndata.AnnData`
            `anndata.AnnData` object to store the single cell experiment information.

        corr_formula: `str` or `list[str]`
            Indicates the groups for correlation structure. If '1', all cells have one estimated corr. If 'ind', no corr (features are independent). If others, this variable decides the corr structures.

        assay_use: `str` (default: None)
            Indicates the assay you will use. If None, please specify a name for the assay stored in `anndata.AnnData.X` in @default_assay_name.

        default_assay_name: `str` (default: None)
            Specified only when @assay_use is None. Asign a name to your default single cell experiment.get_bpparamone)
            The name of cell type variable in the `anndata.AnnData.obs`.

        celltype: `str` (default: None)
            The name of cell type variable in the `anndata.AnnData.obs`.

        pseudotime: `str` or `list[str]` (default: None)
            The name of pseudotime and (if exist) multiple lineages in the `anndata.AnnData.obs`.

        spatial: `list[str]` (default: None)
            The names of spatial coordinates in the `anndata.AnnData.obs`.

        other_covariates: `str` or `list[str]` (default: None)
            The other covaraites you want to include in the data.

        ncell: `int` (default: 'default')
            The number of cell you want to simulate. Default is 'default', which means only the provided cells in the `anndata.AnnData` object will be used. If an arbitrary number is provided, the fucntion will use Vine Copula to simulate a new covaraite matrix.

        mu_formula: `str`
            A string of the mu parameter formula

        sigma_formula: `str`
            A string of the sigma parameter formula

        family_use: `str` or `list[str]`
            A string or a list of strings of the marginal distribution. Must be one of 'binomial', 'poisson', 'nb', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution', 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution' and 'gaussian distribution' respectively.

        usebam: `bool`
            If True, call R function mgcv::bam for calculation acceleration.

        copula: `str` (default: 'gaussian')
            A string of the copula choice. Must be one of 'gaussian' or 'vine'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.

        empirical_quantile: `bool` (default: False)
            Please only use it if you clearly know what will happen! If True, DO NOT fit the copula and use the EMPIRICAL CDF values of the original data; it will make the simulated data fixed (no randomness). Only works if ncell is the same as your original data.

        family_set: `str` or `list[str]` (default: ['gaussian', 'indep'])
            A string or a string list of the bivariate copula families.

        important_feature: `str` or `list[bool]` (default: 'all')
            A string or list which indicates whether a gene will be used in correlation estimation or not. If this is a string, then this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a list, then this should be a logical vector with length equal to the number of genes in @sce. True in the logical vector means the corresponding gene will be included in gene-gene correlation estimation and False in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.

        dt: `bool` (default: True)
            If True, perform the distributional transformation to make the discrete data continuous. This is useful for discrete distributions (e.g., Poisson, NB). Note that for continuous data (e.g., Gaussian), DT does not make sense and should be set as False.

        pseudo_obs: `bool` (default: False)
            If True, use the empirical quantiles instead of theoretical quantiles for fitting copula.

        fastmvn: `bool` (default: False)
            If True, the sampling of multivariate Gaussian is done by R function mvnfast, otherwise by R function mvtnorm.

        nonnegative: `bool` (default: True)
            If True, values < 0 in the synthetic data will be converted to 0. Default is True, since the expression matrix is nonnegative.

        nonzerovar: `bool` (default: False)
            If True, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.

        return_model: `bool` (default: False)
            If True, the marginal models and copula models will be returned.

        n_cores: `int` (default: 'default')
            The number of cores to use. Default is 'default', use the setting when initializing.

        parallelization: `str` (default: 'default')
            The specific parallelization function to use. If 'bpmapply', first call method @get_bpparam. Default is 'default', use the setting when initializing.

        bpparam: `rpy2.robject.methods.RS4` (default: 'default')
            If @parallelization is 'bpmapply', first call function @get_bpparam to get the robject. If @parallelization is 'mcmapply' or 'pbmcmapply', it should be None. Default is 'default', use the setting when initializing.

        trace: `bool` (default: False)
            If True, the warning/error log and runtime for gam/gamlss will be returned.

        return_py: `bool` (default: 'default')
            If True, functions will return a result easy for manipulation in python. Default is 'default', use the setting when initializing.

        Output:
        --------
        A `dict` like object.

        new_count: `pandas.DataFrame`
            A matrix of the new simulated count (expression) matrix.

            The row corresponds to the observations and the column corresponds to the genes.

        new_covariate: `pandas.DataFrame`
            A dataframe of the new covariate matrix.

        model_aic: `pandas.Series`
            An array with three values. In order, they are the marginal AIC, the copula AIC, the total AIC.

        model_bic: `pandas.Series`
            An array with three values. In order, they are the marginal BIC, the copula BIC, the total BIC.

        marginal_list: `rpy2.rlike.container.OrdDict`
            A dict of marginal regression models if return_model = True.

            Caution that though it's name is `list`, it is actually a `dict` like object.

        corr_list: `rpy2.rlike.container.OrdDict`
            A dict of correlation models (conditional copulas) if return_model = True.

            Caution that though it's name is `list`, it is actually a `dict` like object.
        """

        if return_py == "default":
            return_py = self.return_py

        # modify input
        self.family_use = family_use
        corr_formula, celltype, pseudotime, other_covariates, family_use, family_set = _other2list(
            corr_formula, celltype, pseudotime, other_covariates, family_use, family_set
        )

        if ncell == "default":
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
        celltype, pseudotime, spatial, other_covariates, corr_formula, family_use, family_set = _strvec_none2ri(
            celltype, pseudotime, spatial, other_covariates, corr_formula, family_use, family_set
        )

        # check important feature
        if isinstance(important_feature, list):
            important_feature = ro.vectors.BoolVector(important_feature)
        elif not (important_feature == "all" or important_feature == "auto"):
            raise InputError("important_feature should be 'all', 'auto', or a list.")

        # check parallelization parameter
        if parallelization == "default":
            parallelization = self.parallelization
        if n_cores == "default":
            n_cores = self.n_cores
        if bpparam == "default":
            bpparam = self.bpparam
        else:
            bpparam = _bpparamcheck(parallelization, bpparam)

        self.whole_pipeline_res = self._rscdesign.scdesign3(
            sce=self.sce,
            assay_use=self.assay_use,
            corr_formula=corr_formula,
            celltype=celltype,
            pseudotime=pseudotime,
            spatial=spatial,
            other_covariates=other_covariates,
            ncell=ncell,
            mu_formula=mu_formula,
            sigma_formula=sigma_formula,
            family_use=family_use,
            usebam=usebam,
            empirical_quantile=empirical_quantile,
            copula=copula,
            fastmvn=fastmvn,
            DT=dt,
            pseudo_obs=pseudo_obs,
            family_set=family_set,
            important_feature=important_feature,
            nonnegative=nonnegative,
            nonzerovar=nonzerovar,
            return_model=return_model,
            n_cores=n_cores,
            parallelization=parallelization,
            BPPARAM=bpparam,
            trace=trace,
        )

        if return_py:
            with self.__convert.context():
                res = ro.conversion.get_conversion().rpy2py(self.whole_pipeline_res)
            res["model_aic"] = _addname(
                array=res["model_aic"],
                row_name=self.whole_pipeline_res.rx2("model_aic").names,
            )
            res["model_bic"] = _addname(
                array=res["model_bic"],
                row_name=self.whole_pipeline_res.rx2("model_bic").names,
            )
            res["new_count"] = _addname(
                array=res["new_count"].T,
                row_name=self.whole_pipeline_res.rx2("new_count").colnames,
                col_name=self.whole_pipeline_res.rx2("new_count").rownames,
            )
            if res["marginal_list"] is NULL:
                res["marginal_list"] = None
            if res["corr_list"] is NULL:
                res["corr_list"] = None
            return res
        else:
            return self.whole_pipeline_res
