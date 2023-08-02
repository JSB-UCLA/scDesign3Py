from typing import Literal, Optional, Union

import anndata as ad

from ._utils._format import InputError, _anndata_to_sce, _other_to_list, _typecheck


class scDesign3:
    def __init__(self) -> None:
        pass

    # @_typecheck(
    #     anndata=ad.AnnData,
    #     corr_by=Union[str, list],
    #     assay_use=Optional[str],
    #     default_assay_name=Optional[str],
    #     celltype=Optional[str],
    #     pseudotime=Union[str, list],
    #     spatial=list,
    #     other_covariates=Union[str, list],
    #     ncell=Optional[int],
    # )
    def construct_data(
        self,
        anndata: ad.AnnData,
        corr_by: Union[str, list[str]],
        assay_use: Optional[str] = None,
        default_assay_name: Optional[str] = None,
        celltype: Optional[str] = None,
        pseudotime: Union[str, list[str]] = [],
        spatial: list[str] = [],
        other_covariates: Union[str, list[str]] = [],
        ncell: Optional[int] = None,
    ):
        # modify input
        corr_by, pseudotime, other_covariates = _other_to_list(corr_by, pseudotime, other_covariates)

        self._all_covar = pseudotime + spatial + other_covariates
        if not (celltype is None):
            self._all_covar.append(celltype)
        if not self._all_covar:
            raise InputError("At least one covariate should be specified.")

        # Construct R SingleCellExperiment object
        self._sce, self._assay_use = _anndata_to_sce(
            data=anndata,
            assay_use=assay_use,
            default_assay_name=default_assay_name,
            covar=self._all_covar,
        )
