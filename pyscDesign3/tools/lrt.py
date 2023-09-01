from typing import Union

import rpy2.robjects as ro
from rpy2.rlike.container import OrdDict

from .._utils._format import convert


def perform_lrt(
    alter_marginal: Union[ro.vectors.ListVector, OrdDict, dict],
    null_marginal: Union[ro.vectors.ListVector, OrdDict, dict],
):
    """Perform the likelihood ratio test

    Perform the likelihood ratio test to compare two list of marginal models.

    Details:
    ----------
    The function takes two lists of marginal models (by default, the first list is the alternative and the second is the null) from @fit_margnial. Note that LRT only makes sense for NESTED models. This can be quite tricky if you use penalized-splines (e.g., for trajectory data).

    Arguments:
    ----------
    alter_marginal: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict`
        A dict of marginal models from the alternative hypothesis.

    null_marginal: `rpy2.robject.vectors.ListVector` or `rpy2.rlike.container.OrdDict` or `dict`
        A dict of marginal models from the null hypothesis. It must be strictly nested in the alternative model.

    Output:
    ----------
    `pandas.DataFrame`
        A dataframe of the LRT result.
    """

    with convert.context():
        lrt_func = ro.r("scDesign3::perform_lrt")
        lrt_res = lrt_func(alter_marginal, null_marginal)

        return lrt_res
