For developers
==============


How the package will be structured
----------------------------------

Contributing to opencadd
-----------------------------------

First of all, thank you for considering making this project better.

You can help to make opencadd better by raising an issue to report a bug or suggest a
feature which is not implemented yet.
You can also create a new feature and add it via a Pull Request.

How Pull Requests work
----------------------

1. Fork the ``opencadd`` repository into your profile. Now you have your own copy of the repository.
2. Clone the repository. Run ``clone https://github.com/<your_username>/opencadd`` in your terminal.
3. Create a new branch to work on: ``git checkout -b meaningful-branch-name``.
4. Make your changes.
5. Add, commit and push your change to your forked repository.
6. Click on ``Create pull request`` Button.
7. Follow the template given there. Be aware of the requirements_.

The ``opencadd`` team will review the pull request and if there are no flaws it will be merged
to the ``master`` branch. If there are some flaws, the team will request you to fix them.

.. _requirements:

**************************************
What should each Pull Request contain?
**************************************

* Documentation (should be in ``rst`` format)
* Tests_ with Pytest
* Examples for the use of the new functions
* Little benchmark
* If 3rd party dependencies have been added, rationale behind that decision
* It should be formatted by black with ``-l 99``
* Short summary of the changes that were made
* If you implemented a new method or bigger feature, please porvide a tutorial for this method. For this you can follow this template.

Some people have trouble with NGLview. If that is the case for you, follow this `troubleshooting guide
<https://github.com/SBRG/ssbio/wiki/Troubleshooting#tips-for-nglview>`_.

LINK IS MISSING, WILL BE ADDED AS SOON AS TEMPLATE IS IN MASTER.

Formatting with black
---------------------

**1. Option:**

* Apply "black -l 99" on your code before committing::

        $> black -l 99 opencadd/


**2. Option:**

* Configuring your IDE to that effect
* Example in VS Code:

    * go to the settings
    * search for "python formatting"
    * set "Python › Formatting: Provider" to black
    * add "-l 99" to "Python › Formatting: Black Args"
    * activate "Editor: Format On Save"


.. _Tests:

How to add a new test
---------------------

- write a unit test for the new (or changed) function with `pytest
  <https://docs.pytest.org/en/latest/>`_.
- add new dependencies to ``test_env.yaml``


Steps made by the Github Actions Workflow
-----------------------------------------

The actions are executed automatically for every Pull Request submitted,
and for every commit pushed to ``master``. Steps run in Ubuntu and MacOS are:

* Report additional information about the test-build.
* Fixing conda in MacOS (to get the project).
* Creating the environment and getting all necessary dependencies.
* Installing the package in this environment.
* Running the tests.

The formating check is done in ubuntu.

* Checkout the code.
* Installing the linter (pylint) and the formatter (black).
* Running pylint (using  configuration at ``.pylintrc``).
* Running ``black -l 99 --check`` (check mode).
