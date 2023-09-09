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

## Step 3: Specify `R_HOME`

When installing pyscDesign3, `rpy2` package will be installed using pip. To make sure `rpy2` can be correctly installed, you need to first specify the intended **R** interpreter by add the corresponding `R_HOME` to SYSTEM variable.

To get `R_HOME`, run the following code in **R** console.

```r
R.home()
```

- For linux users

**Please specify the `R_HOME` before you install the `rpy2` !**

```bash
export R_HOME=your_R_HOME_path
```

- For windows users

Please specify the `R_HOME` before you use the package and it's fine to set the SYSTEM variable after installing the `rpy2`.

```powershell
setx R_HOME "your_R_HOME_path"
```

Or from GUI

1. Press **Windows + R** to open the Windows Run prompt.
2. Type in **sysdm.cpl** and click **OK**.
3. Open the **Advanced** tab and click on the **Environment Variables** button in the System Properties window.
4. Click the **New…** button under the user-specific section.
5. Set variable name as **R_HOME** and value as **your_R_HOME_path**.

## Step 4: Install pyscDesign3

To install the latest version from [Github](https://github.com/DILIU899/pyscDesign3),

run:

```powershell
git clone https://github.com/DILIU899/pyscDesign3.git
cd pyscDesign3
pip install .
```

Or quick install from PyPI:

```powershell
pip install pyscDesign3
```

```{eval-rst}
.. Note::
    If there's any problem in installing `rpy2`, please refer to `rpy2 installation tutorial <https://pypi.org/project/rpy2/>`_.
```

## Check for installation

Run the following code in python:

```python
import pyscDesign3
```

If successfully installed and imported, the R loacation will be printed.
