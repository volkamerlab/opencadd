Our object model
================

Explain how Atomium can be used for all kinds of IO and structure manipulation.

Detail some naming conventions used in Atomium and how to do some basic tasks

What's the difference between a Model and a Molecule?

Some useful methods in each of them?

Questions:
- How to access all alpha carbons in the protein?
- How to remove all waters
- How to remove all ligands
- How to select a subset of the protein (all isoleucines)

Alignment API
=============

Agree on the required arguments (and their names) for every alignment engine,
and the returned objects so they are consistent across modules.

For this we will define the following API:

- Input parameters:
    - atomium Models [link to docs]
    - optional kwargs
- Returns:
    - superposed structures [using atomium Models/Molecules]
    - scores

How the package is structured
==============================

Explain the directory tree and what each subpackage should contain and why


What each PR should contain
============================

* Docs
* Tests
* Examples
* Little benchmark
* Explain added 3rd party dependencies

Maybe create a GitHub Pull Request template

How to add a new test
=====================

Little example