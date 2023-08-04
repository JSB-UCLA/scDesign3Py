import warnings
from collections import OrderedDict
from functools import wraps
from inspect import signature
from typing import Optional, Union, get_args, get_origin

import anndata as ad
import numpy as np
import pandas as pd
from rpy2.robjects import NULL, default_converter, numpy2ri, pandas2ri, r
from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector, ListVector, Matrix, StrVector


class ConvertionError(RuntimeError):
    pass


class InputError(RuntimeError):
    pass


def _recurse_r_tree(data):
    """
    step through an R object recursively and convert the types to python types as appropriate.
    Leaves will be converted to e.g. numpy arrays or lists as appropriate and the whole tree to a dictionary.
    """
    r_dict_types = [DataFrame, ListVector]
    r_array_types = [FloatVector, IntVector, Matrix]
    r_list_types = [StrVector]
    if type(data) in r_dict_types:
        return OrderedDict(zip(data.names, [_recurse_r_tree(elt) for elt in data]))
    elif type(data) in r_list_types:
        return [_recurse_r_tree(elt) for elt in data]
    elif type(data) in r_array_types:
        return np.array(data)
    else:
        if hasattr(data, "rclass"):  # An unsupported r class
            raise KeyError(
                "Could not proceed, type {} is not defined"
                "to add support for this type, just add it to the imports "
                "and to the appropriate type list above".format(type(data))
            )
        else:
            return data  # We reached the end of recursion


def _anndata_to_sce(data: ad.AnnData, assay_use=None, default_assay_name=None, covar=None):
    """Extract anndata info used for scDesign3 and change into R SingleCellExperiment"""

    ## check unique cell names and gene names
    if not data.obs.index.is_unique:
        warnings.warn(
            "Warning: The input AnnData object may have duplicate cell names. AnnData.obs_names_make_unique is called.",
            UserWarning,
        )
        data.obs_names_make_unique()
    if not data.var.index.is_unique:
        warnings.warn(
            "Warning: The input AnnData object may have duplicate gene names. AnnData.var_names_make_unique is called.",
            UserWarning,
        )
        data.var_names_make_unique()

    ## check gene info
    try:
        gene_info = data.var.iloc[:, [0]]
    except:
        gene_names = data.var_names
        gene_info = pd.DataFrame([None for i in range(len(gene_names))], index=gene_names, columns=["default"])

    ## check cell info dtype
    if covar is None:
        cell_info = data.obs
    else:
        try:
            cell_info = data.obs[covar]
        except KeyError:
            raise InputError(
                "At least one of the input covariate: {} doesn't exist in your Anndata.obs dataframe.".format(
                    ", ".join(covar)
                )
            )

    for _, dtype in cell_info.dtypes.items():
        if (dtype != "category") and (not np.issubdtype(dtype, np.number)):
            raise ConvertionError(
                "Please make sure the given Anndata.obs dataframe has only category and numeric dtypes."
            )

    ## check input assay use
    if assay_use is None:
        count_matrix = data.X.toarray().T
        if default_assay_name is None:
            raise ConvertionError(
                "If use Anndata.X as the default count matrix, please specify the default_assay_name argument used in R SingleCellExperiment as the default assay name."
            )
        assay_use = default_assay_name

    else:
        try:
            count_matrix = data.layers[assay_use].toarray().T
        except:
            raise ConvertionError("The specified assay name doesn't exist in anndata object. Please check.")

    with (default_converter + pandas2ri.converter + numpy2ri.converter).context():
        sce = r("SingleCellExperiment::SingleCellExperiment")(
            assay=ListVector({assay_use: count_matrix}),
            colData=cell_info,
            rowData=gene_info,
        )

    return sce, assay_use


def _other_to_list(*args):
    """
    If None, return None

    If list, return it

    If other type, return a list with that as a element
    """
    return [x if isinstance(x, list) or (x is None) else [x] for x in args]


def _strvec_none2ri(*args):
    """
    Para list[str] and None should call this func to change format

    If str or list[str], return StrVector

    If None, return NULL
    """
    return [NULL if x is None else StrVector(x) for x in args]

def _parapara(parallelization, bpparam):
    if (bpparam is None) and ((parallelization == 'pbmcmapply') or (parallelization == 'mcmapply')):
        return NULL
    elif (not (bpparam is None)) and (parallelization == 'bpmapply'):
        return bpparam
    else:
        raise InputError('If parallelization is pbmcmapply or mcmapply, bpparam should be None. If parallelization is bpmapply, please run get_bpparam to get bpparam.')


def _typecheck(*type_args, **type_kwargs):
    """
    This function supports check types including type in python built in library and third party library. It also supports several types packaged in typing.Union
    """

    def decorator(func):
        sig = signature(func)
        # mapping relationship between function parameters and decorator convention parameter
        bound_types = sig.bind_partial(*type_args, **type_kwargs).arguments

        @wraps(func)
        def wrapper(*args, **kwargs):
            # Obtain the actual value passed in during function execution
            bound_values = sig.bind(*args, **kwargs)
            # type check
            for name, value in bound_values.arguments.items():
                if name in bound_types:
                    from_typing = get_origin(bound_types[name])
                    if (from_typing is Union) or (from_typing is Optional):
                        supported = get_args(bound_types[name])
                        if not any(map(lambda x: isinstance(value, x), supported)):
                            right_types = " or ".join(
                                map(
                                    lambda x: x.__name__,
                                    supported,
                                )
                            )
                            raise TypeError(f"Argument {name} must be {right_types}")
                    else:
                        if not isinstance(value, bound_types[name]):
                            right_types = bound_types[name].__name__
                            raise TypeError(f"Argument {name} must be {right_types}")
            return func(*args, **kwargs)

        return wrapper

    return decorator
