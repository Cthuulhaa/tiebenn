Requirements
============

* **Python 3.9** or later (successfully tested with Python 3.10 and 3.12).
* **SeisBench**, the ML model toolbox used for phase picking and denoising.
* **NonLinLoc**, a suite of C programs for probabilistic hypocenter estimation.
* **GMT** and **PyGMT**, for map generation.
* **PyOcto**, phase associator after `Münchmeyer (2024) <https://seismica.library.mcgill.ca/article/view/1130>`_.
* **Pyrocko**, open-source seismology toolbox and library.
* **GaMMA**, phase associator after `Zhu et al. (2022) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021JB023249>`_.

Optional
~~~~~~~~

* **NLLGrid**, a Python class for handling NonLinLoc grid files. Hosted `here <https://github.com/claudiodsf/nllgrid>`_. Useful if you would eventually like to try 3D NonLinLoc grids in Python.
* **Sphinx**, in case you want to generate TieBeNN's documentation.
