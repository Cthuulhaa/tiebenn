# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TieBeNN'
copyright = '2025, C. Ramos'
author = 'C. Ramos'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]

autosummary_generate = True

autodoc_mock_imports = [
    "gamma",         # GaMMA
    "nllgrid",       # NonLinLoc grid tools
    "joblib",        # parallelization
    "matplotlib",    # visualization
    "cartopy",       # optional in some visual modules
    "obspy",         # waveform handling
    "pyrocko",       # optional dependency
    "pyocto",        # for PyOcto-based associator
    "torch",         # needed for seisbench-based modules
    "seisbench",
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
    'imported-members': True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/tiebenn_logo.svg'
html_css_files = ['css/custom.css']
html_favicon = '_static/tiebenn_favicon.png'
