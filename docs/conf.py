# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import pathlib
import sys

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
# see https://pypi.org/project/setuptools-scm/ for details
from pkg_resources import get_distribution


root = pathlib.Path(__file__).parent.parent.absolute()
os.environ["PYTHONPATH"] = str(root)
sys.path.insert(0, str(root))
print("python exec:", sys.executable)
print("sys.path:", sys.path)

import seaflux  # isort:skip

# -- Project information -----------------------------------------------------

project = "seaflux"
copyright = "2021, Luke Gregor"
author = "Luke Gregor"

release = get_distribution("seaflux").version
# for example take major/minor
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "nbsphinx",
    "recommonmark",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "numpydoc",
    "nbsphinx",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "sphinxcontrib.srclinks",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

html_theme_options = {
    "logo_only": True,
    "display_version": False,
    "style_nav_header_background": "#1d3550",
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
# html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
# html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_favicon = "favicon.ico"
html_logo = "img/seaflux_logo.png"
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

html_css_files = [
    "css/custom.css",
]

# -- nbsphinx specific options ----------------------------------------------
# this allows notebooks to be run even if they produce errors.
nbsphinx_allow_errors = True
