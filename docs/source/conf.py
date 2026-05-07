"""Sphinx configuration for polartools.

Built with pydata-sphinx-theme + MyST + sphinx-autoapi. Sources are Markdown
(MyST) under ``docs/source/``; the API reference is auto-generated from the
package docstrings (numpydoc style) by sphinx-autoapi.
"""

import pathlib

_pkg_dir = pathlib.Path(__file__).parents[2] / "polartools"

import polartools  # noqa: E402

# -- Project information -----------------------------------------------------

project = "polartools"
copyright = "2020-2026, UChicago Argonne, LLC"
author = "Argonne National Laboratory"
version = polartools.__version__
release = polartools.__version__

# -- General configuration ---------------------------------------------------

extensions = [
    "autoapi.extension",
    "myst_parser",
    "nbsphinx",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_tabs.tabs",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "matplotlib.sphinxext.plot_directive",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "numpydoc",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

templates_path = ["_templates"]

# Suppress noise from autoapi's static cross-module import resolution: the
# hoisted _larch / _pyrixs modules are deliberately excluded from autoapi
# (they're upstream re-exports), so cross-references to them from absorption /
# process_images can't be resolved without changing this list of warnings.
suppress_warnings = ["autoapi.python_import_resolution"]

# numpydoc plays nicely with napoleon when napoleon is disabled for numpy
# style — keep numpydoc for the existing docstring style across the package.
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
napoleon_numpy_docstring = False  # numpydoc handles numpy-style; avoid double-parse
napoleon_google_docstring = True

# -- sphinx-autoapi ----------------------------------------------------------

autoapi_dirs = [str(_pkg_dir)]
autoapi_root = "api"
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_ignore = [
    "*/tests/*",
    "*/conftest*",
    "*/_version*",
    "*/_larch*",
    "*/_pyrixs*",
]
autoapi_keep_files = True
autoapi_add_toctree_entry = False
autoapi_python_class_content = "both"

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "lmfit": ("https://lmfit.github.io/lmfit-py/", None),
}

# -- nbsphinx ----------------------------------------------------------------

nbsphinx_execute = "never"

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = "polartools"
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/APS-4ID-POLAR/polartools",
    "use_edit_page_button": True,
    "show_toc_level": 2,
    "navbar_align": "left",
    "navigation_with_keys": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version"],
    "logo": {"text": "polartools"},
    "icon_links": [
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/polartools/",
            "icon": "fa-brands fa-python",
        },
    ],
}

html_context = {
    "github_user": "APS-4ID-POLAR",
    "github_repo": "polartools",
    "github_version": "main",
    "doc_path": "docs/source",
}

# Configuration options for plot_directive.
plot_html_show_source_link = False
plot_html_show_formats = False
