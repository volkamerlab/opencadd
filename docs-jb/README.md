# openCADD Documentation

## Compiling opencadd's Documentation

The docs for this project are built with [Sphinx](http://www.sphinx-doc.org/en/master/).
To compile the docs, first ensure that Sphinx and the ReadTheDocs theme are installed.


```bash
conda install sphinx sphinx_rtd_theme
```


Once installed, you can use the `Makefile` in this directory to compile static HTML pages by
```bash
make html
```

The compiled docs will be in the `_build` directory and can be viewed by opening `index.html` (which may itself
be inside a directory called `html/` depending on what version of Sphinx is installed).


## `/source`

### `/source/_static`

**Static Doc Directory**

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
templates_path = ['_static']
```
**Examples of file to add to this directory:**
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images


### `/source/_templates`

**Templates Doc Directory**

Add any paths that contain templates here, relative to  
the `conf.py` file's directory.
They are copied after the builtin template files,
so a file named "page.html" will overwrite the builtin "page.html".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
html_static_path = ['_templates']
```

**Examples of file to add to this directory:**
* HTML extensions of stock pages like `page.html` or `layout.html`
