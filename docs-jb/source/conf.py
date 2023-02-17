# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'openCADD'
author = 'Volkamer Lab'
copyright = '2023 Volkamer Lab'
release = '2.0.0'


# -- Internationalization ------------------------------------------------
# specifying the natural language populates some key tags
language = "en"

# -- Templates
# Add any paths that contain templates here, relative to this directory (i.e. source).
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_baseurl = ''
html_favicon = ''
html_sourcelink_suffix = ''
html_title = ''
html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_theme = 'pydata_sphinx_theme'
html_logo = 'logo.svg'


# This is added due to this issue:
#  https://github.com/pydata/pydata-sphinx-theme/issues/1094#issuecomment-1368264928
html_theme_options = {
    "logo": {
        # "text": "an open-source Python library for computer-aided drug design",
        "image_light": "logo.svg",
        "image_dark": "logo.svg",
    },
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_align": "content",  # alt: "left" or "right"
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_persistent": ["search-button"],
    "icon_links": [
            {
                # Label for this link
                "name": "GitHub",
                # URL where the link will redirect
                "url": "https://github.com/volkamerlab/opencadd/",  # required
                # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
                "icon": "fa-brands fa-square-github",
                # The type of image to be used (see below for details)
                "type": "fontawesome",
            }
       ],
    "search_bar_text": "search openCADD ...",
    "primary_sidebar_end": ["indices", "sidebar-ethical-ads"],
    "secondary_sidebar_items": ["page-toc"],  # "edit-this-page", "sourcelink",
    "show_prev_next": True,
    "footer_items": ["navbar-logo", "copyright"],  # "sphinx-version", "theme-version"
}

html_sidebars = {
    "**": ["search-field", "sidebar-nav-bs"]
}

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# -- Extensions
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'numpydoc',
    'myst_nb',
    'jupyter_sphinx',
    'sphinx_design',
    'sphinxcontrib.bibtex',
    # To show a copy button next to code blocks; see: https://sphinx-copybutton.readthedocs.io/en/latest/
    'sphinx_copybutton',
    'sphinx_togglebutton',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    # 'sphinx.ext.napoleon',
]




autosummary_generate = True
autosummary_ignore_module_all = False
autosummary_imported_members = False

numpydoc_attributes_as_param_list = False
numpydoc_citation_re = 'r"[\\w-]+"'
numpydoc_class_members_toctree = False
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False
numpydoc_use_plots = True
# numpydoc_xref_ignore = "all"
numpydoc_xref_param_type = True

# # Napoleon settings (ref: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html)
# napoleon_google_docstring = False
# napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = False
# napoleon_include_private_with_doc = False
# napoleon_include_special_with_doc = True
# napoleon_use_admonition_for_examples = False
# napoleon_use_admonition_for_notes = False
# napoleon_use_admonition_for_references = False
# napoleon_use_ivar = False
# napoleon_use_param = True
# napoleon_use_rtype = True
# napoleon_preprocess_types = False
# napoleon_type_aliases = None
# napoleon_attr_annotations = True

myst_enable_extensions = [
    'colon_fence',
    'dollarmath',
    # 'linkify',
    'substitution',
    'tasklist'
]
myst_url_schemes = ['mailto', 'http', 'https']


# EPUB options
epub_show_urls = 'footnote'


bibtex_bibfiles = ['references.bib']




# ------ From jupyter-book

# add_module_names = False
#
# comments_config = {'hypothesis': False, 'utterances': False}
# execution_allow_errors = False
# execution_excludepatterns = []
# execution_in_temp = False
# execution_timeout = 30
# extensions = [
#     'sphinx_thebe',
#     'sphinx_comments',
#     'sphinx.ext.duration',
#     'sphinx.ext.doctest',
#     'sphinx_jupyterbook_latex'
# ]
#
# jupyter_cache = ''
# jupyter_execute_notebooks = 'force'
# latex_engine = 'pdflatex'
#
# nb_output_stderr = 'show'
# numfig = True
#
# pygments_style = 'sphinx'
# suppress_warnings = ['myst.domains']
# use_jupyterbook_latex = True
# use_multitoc_numbering = True
