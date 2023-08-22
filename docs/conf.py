# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
from pathlib import Path
project = 'sdof'
copyright = '2023, Claudio M. Perez'
author = 'Claudio M. Perez'
description = "Lightning-fast integration for single-degree-of-freedom systems."
abstract = """
This package solves scalar differential equations of the form

   $$m \ddot{u} + c \dot{u} + k u = f(t)$$

Integration is carried out using a Generalized - :math:`\alpha`
integrator that is implemented under the hood in highly optimized
multi-threaded C code.
"""

version = '0.0.0'
release = '0.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'nbsphinx',
    'sphinx_design'
]

templates_path   = ['_templates']
exclude_patterns = []

source_suffix = '.rst'
root_doc = 'index'
language = 'en'

# -- Options for HTML output -------------------------------------------------

html_title = project
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_favicon = '_static/favicon.ico'
html_css_files = [
    "css/peer.css",
] + [
    'css/css/'+str(file.name) for file in (Path(__file__).parents[0]/"_static/css/css/").glob("*.css")
]
html_additional_pages = {'index': 'home.html'}
html_context = {
    'description': description,
    'highlights': {
        "": "Because <code>sdof</code> is implemented in standard C, linking should never be a problem."
    },
    **globals()
}
html_show_sourcelink = False
html_theme_options = {
    "github_url": f"https://github.com/BRACE2/{project}",
}

autodoc_member_order = 'bysource'
