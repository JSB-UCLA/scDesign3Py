# Installation

```{eval-rst}
.. Attention::
    The **pyscDesign3** package is developed based on `rpy2` and utilize the R interpreter on user's computer.
```

## Step 1: Install R

Instruction for installing **R** can be found at [r-project](https://www.r-project.org/).

If **R** has been installed, this step can be skipped.

## Step 2: Install scDesign3 R package (optional)

Run the following code in **R** console.

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
```

```{eval-rst}
.. Note::
    Even if you forget to install the **scDesign3** **R** package, the interface will try to install **scDesign3** automatically, which may cost extra time when importing the package for the first time. However, if there's unexpected error occurred during automatical installation, please install the **scDesign3** **R** package manually in **R**.
```

## Step 3: Install pyscDesign3

To install the latest version from [Github](https://github.com/DILIU899/pyscDesign3),

run:

```powershell
git clone https://github.com/DILIU899/pyscDesign3.git
cd pyscDesign3
pip install pyscDesign3
```

Or quick install from PyPI:

```powershell
pip install pyscDesign3
```

```{eval-rst}
.. Note::
    If there's any problem in installing `rpy2`, please refer to `rpy2 installation tutorial <https://pypi.org/project/rpy2/>`_.
```
