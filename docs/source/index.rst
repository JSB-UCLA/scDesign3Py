.. scDesign3Py documentation master file, created by
   sphinx-quickstart on Tue Aug 29 21:15:19 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**scDesign3Py**: python interface for **scDesign3** **R** package
====================================================================

The **R** package **scDesign3** is a unified probabilistic framework that generates realistic in silico high-dimensional single-cell omics data of various cell states, including discrete cell types, continuous trajectories, and spatial locations by learning from real datasets. 

**scDesign3Py** is the python interface of the **R** package **scDesign3**, enabling python users to easily use the **scDesign3** functions and share almost the same usage as the **R** version.

Summary of usage of **scDesign3**
-----------------------------------

.. figure:: _static/scDesign3_illustration.png
   :align: center
   :width: 600  
   :alt: scDesign3 overview

|

To find out more details about **scDesign3**, you can check out manuscript on Nature Biotechnology:

`Song, D., Wang, Q., Yan, G. et al. scDesign3 generates realistic in silico data for multimodal single-cell and spatial omics. Nat Biotechnol (2023). <https://www.nature.com/articles/s41587-023-01772-1>`_

For tutorial of **R** package **scDesign3**, you can check the `tutorial website <https://songdongyuan1994.github.io/scDesign3/docs/index.html>`_.

.. toctree::
   :maxdepth: 1
   :caption: Get Started
   
   Installation <set_up/Installation>
   API <set_up/API>
   Quick Start <tutorial/quick_start> 
   About scDesign3Py <set_up/info>

.. toctree::
   :maxdepth: 1
   :caption: Tutorial
   
   All in one simulation <tutorial/all_in_one>
   Step by step simulation <tutorial/stepwise>
   Get BPPARAM <tutorial/bpparam>
   Perform likelihood ratio test <tutorial/lrt>
   Plot results <tutorial/plot>

.. toctree::
   :maxdepth: 1
   :caption: More examples
   
   Model parameter <examples/model_para/index> 
   Data simulation <examples/data_simu/index> 
   Model selection <examples/model_select/index> 
   Model alteration <examples/model_alter/index> 
