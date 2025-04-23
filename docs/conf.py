# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime
sys.path.insert(0, os.path.abspath('..'))

import logging # Temporal!
logging.basicConfig(level=logging.DEBUG) # Temporal!

import tiebenn.tools.utm # Temporal!!!
print("utm imported:", dir(tiebenn.tools.utm)) # Temporal!
import tiebenn.tools.sb_tools # Ttemporal
print("sb_tools imported:", dir(tiebenn.tools.sb_tools)) # Temporal!

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TieBeNN'
copyright = f'{datetime.now().year}, C. Ramos'
author = 'C. Ramos'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
    "sphinx_autodoc_typehints",
]

autodoc_mock_imports = ["pygmt"]

autodoc_typehints = "description"
autodoc_typehints_description_target = "documented"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/tiebenn_logo.svg'
html_css_files = ['css/custom.css']
html_favicon = '_static/tiebenn_favicon.png'
