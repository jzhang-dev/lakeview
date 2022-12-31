# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Lakeview'
copyright = '2022, Jia-Yuan Zhang'
author = 'Jia-Yuan Zhang'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc', "myst_nb", "IPython.sphinxext.ipython_console_highlighting"]
source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
}

#templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
#html_static_path = ['_static']

# -- Options for autodoc -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#module-sphinx.ext.autodoc

autodoc_typehints = 'description'

# -- Options for mysb-nb -------------------------------------------------
# https://myst-nb.readthedocs.io/en/latest/computation/execute.html

nb_execution_mode = "cache"
nb_execution_timeout = 600
