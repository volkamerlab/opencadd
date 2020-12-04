"""
Utilities helpful in the context of Jupyter Notebooks
"""

from IPython.display import display, HTML, Markdown


def print_html(*args, **kwargs):
    """
    Displays rich HTML output
    """
    return display(HTML(*args, **kwargs))


hprint = print_html


def print_markdown(*args, **kwargs):
    """
    Displays rich Markdown output
    """
    return display(Markdown(*args, **kwargs))


mprint = print_markdown