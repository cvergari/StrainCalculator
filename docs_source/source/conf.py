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
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Strain Calculator'
copyright = '2024, Claudio Vergari'
author = 'Claudio Vergari'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinxcontrib.matlab',               
              'sphinx.ext.napoleon',
              'sphinx.ext.mathjax',
              'sphinx.ext.autosectionlabel',
              'sphinx.ext.autosummary',
              'sphinx.ext.todo']

primary_domain = 'mat'
matlab_src_dir = os.path.abspath('./../../..')
print('matlab_src_dir: {}'.format(matlab_src_dir))
assert(os.path.exists(matlab_src_dir))

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, matlab_src_dir)

autodoc_member_order = 'alphabetical'
autosummary_imported_members = True
autosummary_generate = True  # Turn on sphinx.ext.autosummary


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    }
    
# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'
html_theme = 'pyramid'

html_show_sourcelink = False  # Hides the "show source" hyperlink

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

todo_include_todos  = True
