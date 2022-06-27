# This code is part of Qiskit.
#
# (C) Copyright IBM 2018.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# pylint: disable=invalid-name

"""
Sphinx documentation builder
"""

# General options:

project = "Qiskit QEC"
copyright = "2022, Qiskit Development Team"  # pylint: disable=redefined-builtin
author = "Qiskit Development Team"

# The short X.Y version
version = ""
# The full version, including alpha/beta/rc tags
release = "0.0.0"

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.extlinks",
    "jupyter_sphinx",
    "sphinx_autodoc_typehints",
    "reno.sphinxext",
    "nbsphinx",
]
templates_path = ["_templates"]
numfig = True
numfig_format = {"table": "Table %s"}
language = None
pygments_style = "colorful"
add_module_names = False
modindex_common_prefix = ["qiskit_qec."]

# html theme options
html_theme = "qiskit_sphinx_theme"
html_last_updated_fmt = "%Y/%m/%d"
html_theme_options = {
    "logo_only": True,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
}
html_css_files = ["gallery.css"]
html_static_path = ["_static"]

# autodoc/autosummary options
autosummary_generate = True
autosummary_generate_overwrite = False
autoclass_content = "both"

# nbsphinx options (for tutorials)
nbsphinx_timeout = 180
# TODO: swap this with always if tutorial execution is too slow for ci and needs
# a separate job
# nbsphinx_execute = os.getenv('QISKIT_DOCS_BUILD_TUTORIALS', 'never')
nbsphinx_execute = "never"
nbsphinx_widgets_path = ""
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
