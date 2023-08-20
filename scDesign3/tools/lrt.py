import rpy2.robjects as ro

from .._utils._errors import NotInplementedError


def perform_lrt(
    alter_marginal: ro.vectors.ListVector,
    null_marginal: ro.vectors.ListVector,
):
    try:
        alter_marginal = alter_marginal.rx2("fit")
        null_marginal = null_marginal.rx2("fit")

    except:
        raise NotInplementedError("Currently the input of this function can only be the result of fit_marginal_res")
