from typing import Literal

from rpy2.robjects import r

from .._utils._errors import InputError


def get_bpparam(name=Literal["MulticoreParam", "SnowParam"], show=True, **kwargs):
    """Get your parallelization parameters robject

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
