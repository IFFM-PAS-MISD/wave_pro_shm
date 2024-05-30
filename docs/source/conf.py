# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------

project = "WaveProSHM"
copyright = "Pawel Kudela"
author = "Pawel Kudela"

# The full version, including alpha/beta/rc tags
release = "1.2.1"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinxcontrib.matlab",
    "sphinx.ext.autodoc",
    "sphinx_copybutton",
    "myst_parser",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.inkscapeconverter",
    "sphinx_design",
]
myst_enable_extensions = ["colon_fence"]
# Bibliography
bibtex_bibfiles = ['refs.bib']
bibtex_default_style = 'plain'
# latex
# no chapter titles
latex_elements = {
    'fncychap': '\\usepackage[Conny]{fncychap}',
    'preamble': r' \makeatletter \renewcommand{\@chapapp}{} \makeatother',
#    'tableofcontents': ' ', # no table of content
}
# If false, no module index is generated.
latex_domain_indices = False
# Matlab setup
matlab_src_dir = os.path.dirname(os.path.abspath("../../src/"))
primary_domain = "mat"
matlab_short_links = True

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['*.mlapp','build','source/_static/MathJax/README.md']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# The master toctree document.
master_doc = "index"

source_suffix = ['.rst', '.md']
# source_suffix = ".rst"


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# html_logo = '_static/logo.png'

# html_theme_options = {
#     'collapse_navigation': False,
#     'display_version': False,
#     'navigation_depth': 4,
# }
# html configuration
html_show_sourcelink = False

#html_sidebars = {
#    "**": [
#        "about.html",
#        "navigation.html",
#        "relations.html",  # needs 'show_related': True theme option to display
#        "searchbox.html",
#        "donate.html",
#    ]
#}
html_js_files = [
    'js/custom.js',
     'MathJax/es5/tex-mml-chtml.js',
]
