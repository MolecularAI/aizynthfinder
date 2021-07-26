import os
import sys

sys.path.insert(0, os.path.abspath("."))

project = "aizynthfinder"
copyright = "2020, Molecular AI group"
author = "Molecular AI group"
release = "3.0.0"

# This make sure that the cli_help.txt file is properly formated
with open("cli_help.txt", "r") as fileobj:
    lines = fileobj.read().splitlines()
with open("cli_help.txt", "w") as fileobj:
    fileobj.write(".. code-block::\n\n")
    fileobj.write("  " + "\n  ".join(lines))

extensions = [
    "sphinx.ext.autodoc",
]
autodoc_member_order = "bysource"
autodoc_typehints = "description"

html_theme = "alabaster"
html_theme_options = {
    "description": "A fast, robust and flexible software for retrosynthetic planning",
    "fixed_sidebar": True,
}
