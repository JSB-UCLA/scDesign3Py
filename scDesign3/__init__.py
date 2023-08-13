import os

R_HOME = os.environ.get("R_HOME")

if R_HOME is None:
    raise ImportError(
        "R must be installed and R_HOME should be added to the system environment to use this interface."
    )
else:
    print(f"The R project used is located at {R_HOME}")


import rpy2.robjects.packages as rpackages
from rpy2.robjects import r
from rpy2.robjects.vectors import StrVector

# set encoding method in case decoding warning occur
r("invisible(suppressWarnings(Sys.setlocale('LC_ALL', 'en_US.UTF-8')))")

# check if all dependencies installed
r_requirements = {
    "cran": ["devtools"],
    "dev": ["scDesign3"],
    "bioc": [],
}

cran_install = [x for x in r_requirements["cran"] if not rpackages.isinstalled(x)]
bioc_install = [x for x in r_requirements["bioc"] if not rpackages.isinstalled(x)]
dev_install = [x for x in r_requirements["dev"] if not rpackages.isinstalled(x)]
all_install = cran_install + bioc_install + dev_install

try:
    if cran_install:
        utils = rpackages.importr("utils")
        # select a mirror for R packages
        utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
        utils.install_packages(StrVector(cran_install))

    if bioc_install:
        biocmanager_install = r("BiocManager::install")
        biocmanager_install(StrVector(bioc_install), update=False)

    if dev_install:
        if "scDesign3" in dev_install:
            r('devtools::install_github("SONGDONGYUAN1994/scDesign3")')

except:
    raise ImportError(
        "Dependencies auto installation failed. Please check in R manually if packages {} are all correctly installed.".format(
            ", ".join(all_install)
        )
    )

if all_install:
    print("Dependencies {} are installed in R.".format(", ".join(all_install)))


from scDesign3.tools.get_bpparam import get_bpparam
from scDesign3.tools.lrt import perform_lrt
from scDesign3.tools.plot import plot_reduceddim

from ._core import scDesign3
