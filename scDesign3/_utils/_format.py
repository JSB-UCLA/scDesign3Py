import warnings
from collections import OrderedDict
from functools import wraps
from inspect import signature
from typing import Optional, Union, get_args, get_origin

import anndata as ad
import numpy as np
import pandas as pd
from rpy2.robjects import numpy2ri, pandas2ri, r
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

    pandas2ri.activate()
    numpy2ri.activate()
    sce = r("SingleCellExperiment::SingleCellExperiment")(
        assay=ListVector({assay_use: count_matrix}),
        colData=cell_info,
        rowData=gene_info,
    )

    return sce, assay_use


def _other_to_list(*args):
    result = []
    for x in args:
        if isinstance(x, list) or (x is None):
            result.append(x)
        else:
            result.append([x])

    return result


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
                    if not isinstance(value, bound_types[name]):
                        from_typing = get_origin(bound_types[name])
                        if (from_typing is Union) or (from_typing is Optional):
                            right_types = " or ".join(
                                map(
                                    lambda x: x.__name__,
                                    get_args(bound_types[name]),
                                )
                            )
                            raise TypeError(f"Argument {name} must be {right_types}")

                        else:
                            right_types = bound_types[name].__name__
                            raise TypeError(f"Argument {name} must be {right_types}")
            return func(*args, **kwargs)

        return wrapper

    return decorator
