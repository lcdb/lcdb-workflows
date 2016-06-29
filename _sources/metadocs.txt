Building documentation
======================

From the top-level dir of the repo:

.. code:: bash

  bash docs/build-docs.sh

Then check http://lcdb.github.io/lcdb-workflows to see the updates.

Specifically, this script:

  - build the sphinx documentation in its current state in the current
    directory using ``(cd docs && make html)``
  - clone the repo from github into a temp dir
  - checkout the `gh-pages
    <https://help.github.com/articles/creating-project-pages-manually/>`_
    branch
  - copy the built documentation over to the temp repo's ``gh-pages`` branch
  - add everything that was copied
  - make a git commit (referring to the current commit of the entire repo)
  - push the changes to github

There are some custom additions to the standard sphinx Makefile as dependencies
of the ``html`` rule, s

    - the testing config file is copied over to the docs dir if it doens't
      exist or has been updated
    - if any of the workflows' Snakefiles have been updated, a rulegraph and
      full dag is built for that workflow and the corresponding images are
      saved in the docs/source/images dir
