# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'openCADD'
copyright = '2023, Volkamer Lab'
author = 'Volkamer Lab'
release = '2.0.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# -- Extensions
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    # To show a copy button next to code blocks
    #  see: https://sphinx-copybutton.readthedocs.io/en/latest/
    'sphinx_copybutton',
    'sphinx_design',
]

autosummary_generate = True


# Napoleon settings (ref: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html)
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# EPUB options
epub_show_urls = 'footnote'



# -- Templates
# Add any paths that contain templates here, relative to this directory (i.e. source).
templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_css_files = ['css/custom.css']
html_theme = 'pydata_sphinx_theme'

# This is added due to this issue:
#  https://github.com/pydata/pydata-sphinx-theme/issues/1094#issuecomment-1368264928
html_theme_options = {
    "logo": {
        # "text": "an open-source Python library for computer-aided drug design",
        "image_light": "logo.png",
        "image_dark": "logo.png",
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
