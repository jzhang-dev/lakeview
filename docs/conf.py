# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Lakeview"
copyright = "2022, Jia-Yuan Zhang"
author = "Jia-Yuan Zhang"
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",
    "myst_nb",
    "IPython.sphinxext.ipython_console_highlighting",
]
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
}

# templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/jzhang-dev/lakeview",
    "use_repository_button": True,
    "use_fullscreen_button": False,
}
# html_static_path = ['_static']

# -- Options for autodoc -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#module-sphinx.ext.autodoc

autodoc_typehints = "description"
autodoc_class_signature = "separated"

autodoc_default_options = {
    "member-order": "bysource",
    "show-inheritance": True,
}
# autodoc_type_aliases = {'GroupIdentifier': "lakeview._type_alias.GroupIdentifier"}

# -- Options for Intersphinx -------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pysam": ("https://pysam.readthedocs.io/en/latest", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None)
}


# -- Options for mysb-nb -------------------------------------------------
# https://myst-nb.readthedocs.io/en/latest/computation/execute.html

nb_execution_mode = "cache"
nb_execution_timeout = 600
