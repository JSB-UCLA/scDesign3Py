from typing import Union

import rpy2.robjects as ro
from rpy2.rlike.container import OrdDict

from .._utils._errors import NotInplementedError
from .._utils._format import convert


def perform_lrt(
    alter_marginal: Union[ro.vectors.ListVector, OrdDict],
    null_marginal: Union[ro.vectors.ListVector, OrdDict],
):
    with convert.context():
        lrt_func = ro.r("scDesign3::perform_lrt")
        lrt_res = lrt_func(alter_marginal, null_marginal)

        return lrt_res
