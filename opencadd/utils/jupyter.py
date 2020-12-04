"""
Utilities helpful in the context of Jupyter Notebooks
"""

from IPython.display import display, HTML, Markdown


def print_html(*args, sep=" ", **kwargs):
    """
    Displays rich HTML output

    Parameters
    ----------
    args : str
        HTML text to be rendered
    sep : str, optional=" "
        Default separator between items in `args`
    """
    contents = sep.join([str(a) for a in args])
    return display(HTML(contents, **kwargs))


hprint = print_html


def print_markdown(*args, sep=" ", alert=None, hidden=None, **kwargs):
    """
    Displays rich Markdown output

    Parameters
    ----------
    args : str
        Markdown-formatted text to be rendered
    sep : str, optional=" "
        Default separator between items in `args`
    alert : str, optional=None
        Colored box style. Must be one of info, warning, success or danger.
    hidden : str, optional=None
        If set, use value as the title of a collapsed box that can be revealed
        on click.
    """
    contents = sep.join([str(a) for a in args])
    if alert is not None:
        contents = "\n\n".join(
            [f'<div class="alert alert-block alert-{alert}">', contents, "</div>"]
        )
    if hidden is not None:
        contents = "\n\n".join(
            [f"<details>", f"<summary>{hidden}</summary>", contents, "</details>"]
        )
    return display(Markdown(contents, **kwargs))


mprint = print_markdown